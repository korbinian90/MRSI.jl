function lipid_decontamination_L1(csi, lipid_basis, brain_mask, n_loops=5; args...)
    FT_mask = elliptical_mask(size(csi)[1:2], [1, 1], size(csi, 1) / 2 + 1) # why unsymmetrical?
    FT0 = FT(csi, FT_mask)
    for i in 1:n_loops
        println("Decontamination Loop $i of $n_loops. Norm Results: ")
        csi = lipid_suppression(csi, lipid_basis, brain_mask, FT_mask, FT0; args...)
    end
    return csi
end

FT(x, mask) = 1 / eltype(x).(sqrt(length(x))) * fftshift(fft(ifftshift(x))) .* mask
FT_adjoint(x) = eltype(x).(sqrt(length(x))) * fftshift(ifft(ifftshift(x)))
lipid_xfm(spectrum, lipid_basis) = lipid_basis * conj.(spectrum)

function lipid_suppression(x0::AbstractArray{Complex{T}}, lipid_basis, brain_mask, FT_mask, FT0; maximum_iterations=10, maximum_line_iterations=150, gradient_tolerance=1f-30, alpha=0.01f0, beta=0.6f0, t0=1, p_norm=1, xfm_weight=1f-3, L1_smooth=1f-15) where T
    x = copy(x0)
    g0 = w_gradient(x0, FT0, lipid_basis, FT_mask, brain_mask, xfm_weight)
    dx = -g0
    f1, consis, lip = zero(T), zero(T), zero(T)

    for iteration in 0:maximum_iterations
        FT_x = FT(x, FT_mask)
        FT_dx = FT(dx, FT_mask)

        f0 = first(objective(FT_x, FT_dx, x, dx, 0, FT0, brain_mask, lipid_basis, p_norm, xfm_weight, L1_smooth))
        t = t0

        line_iter = 0
        for i in 1:maximum_line_iterations
            f1, consis, lip = objective(FT_x, FT_dx, x, dx, t, FT0, brain_mask, lipid_basis, p_norm, xfm_weight, L1_smooth)
            if f1 ≤ f0 - alpha * t * abs(g0 ⋅ dx)
                break
            end
            line_iter = i
            t *= beta
        end

        if line_iter > 2
            t0 *= beta
        elseif line_iter < 1
            t0 /= beta
        elseif line_iter == maximum_line_iterations
            println("Reached max line search,.... not so good... might have a bug in operators. exiting... ")
            return
        end

        @show iteration f1 consis lip

        x += t * dx
        g1 = w_gradient(x, FT0, lipid_basis, FT_mask, brain_mask, xfm_weight)
        bk = (g1 ⋅ g1) / ((g0 ⋅ g0) + eps(Float32))
        g0 = g1
        dx = -g1 + bk * dx

        if (norm(dx) < gradient_tolerance)
            break
        end
    end

    println("L2 of consistency: $(sqrt(consis))")
    println("L1 of lipid: $(lip/xfm_weight)")
    return x
end

function objective(FT_x::AbstractArray{Complex{T}}, FT_dx, x, dx, t, FT0, mask, lipid_basis, p_norm, xfm_weight, L1_smooth) where T
    obj = norm(FT_x + t * FT_dx - FT0)^2
    G = zero(T)
    
    if xfm_weight > 0
        for J in CartesianIndices(mask)
            if mask[J]
                x_update = x[J, :] + t * dx[J, :]
                Dx = lipid_xfm(x_update, lipid_basis)
                D2 = real.(Dx .* conj.(Dx))
                G += sum((D2 .+ L1_smooth) .^ T(p_norm / 2))
            end
        end
    end

    XFM = G * xfm_weight
    res = obj + XFM
    return res, obj, XFM
end

function w_gradient(x, FT0, lipid_basis, FT_mask, brain_mask, xfm_weight)
    grad_obj = 2 * FT_adjoint(FT(x, FT_mask) - FT0)

    grad_xfm = if xfm_weight != 0
        lipid_gradient(x, lipid_basis, brain_mask)
    else
        0
    end

    grad = grad_obj + xfm_weight * grad_xfm
    return grad
end

function lipid_gradient(x, lipid_basis, mask)
    grad = similar(x)
    for J in CartesianIndices(mask)
        if mask[J]
            Dx = lipid_xfm(x[J, :], lipid_basis)
            grad[J, :] = sign.(Dx)' * lipid_basis
        end
    end
    return grad
end

function elliptical_mask(size, coefficients, radius, epsilon=1e-6)
    kspace_center = (size .÷ 2) .+ 1
    is_inside(I) = sum(((Tuple(I) .- kspace_center) ./ coefficients) .^ 2) ≤ radius^2 + epsilon
    mask = is_inside.(CartesianIndices(size))
    return mask
end
