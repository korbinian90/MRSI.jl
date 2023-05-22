function lipid_suppression!(csi, mask, lipid_mask=trues(size(csi)[1:3]); type=:L2, L2_beta=0.1, L1_n_loops=5, args...)
    for slice in axes(csi, 3)
        csi[:, :, slice, :] = lipid_suppression_slice(csi[:, :, slice, :], mask[:, :, slice], lipid_mask[:, :, slice], type, L2_beta, L1_n_loops, args...)
    end
end

function lipid_suppression_slice(csi, mask, lipid_mask, type, L2_beta, L1_n_loops, args...)
    fft_dim = 3
    csi = fftshift(fft(csi, fft_dim), fft_dim)

    lipid_basis = get_lipid_basis(csi, lipid_mask)
    if type == :L2
        csi = lipid_decontamination_L2(csi, lipid_basis, mask, L2_beta)
    elseif type == :L1
        csi = lipid_decontamination_L1(csi, lipid_basis, mask, L1_n_loops; args...)
    end

    csi = ifft(ifftshift(csi, fft_dim), fft_dim)
    return csi
end

function get_lipid_basis(csi, lipid_mask)
    count = 0
    lipid_basis = zeros(eltype(csi), sum(lipid_mask), size(csi, 3))

    for J in CartesianIndices(lipid_mask)
        if lipid_mask[J]
            count += 1
            lipid_basis[count,:] = csi[J, :]
        end
    end

    return lipid_basis
end

function lipid_decontamination_L2(csi, lipid_basis, mask, beta)
    beta *= 1e-15
    lipid_basis = transpose(lipid_basis)
    lipid_inv = inv(I + beta * (lipid_basis * lipid_basis'))

    for J in CartesianIndices(mask)
        if mask[J]
            csi[J, :] = lipid_inv * csi[J, :]
        end
    end

    csi .*= mask
    return csi
end
