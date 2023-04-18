function coil_combine(data, refscan)
    coil_weight = conj(refscan[:, :, :, [1], :]) # refscan - only use one timepoint and one dim6
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
