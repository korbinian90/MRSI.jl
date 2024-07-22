function coil_combine(data, refscan; ref_point_for_combine=5)
    pt = max(1, min(ref_point_for_combine, size(refscan, 4))) # between 1 and length of time dimension
    coil_weight = conj(refscan[:, :, :, [pt], :])
    combined = combine(data, coil_weight)
    return combined .* scaling(coil_weight)
end

function scaling(coil_weight)
    scale = 1 ./ sum(abs.(coil_weight) .^ 2; dims=5)
    scale[.!isfinite.(scale)] .= 0
    return scale
end

function combine(data, weight)
    combined = zeros(eltype(data), size(data)[1:4])
    for cha in axes(data, 5)
        combined .+= data[:, :, :, :, cha] .* weight[:, :, :, :, cha]
    end
    return combined
end
