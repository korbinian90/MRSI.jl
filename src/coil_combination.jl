function coil_combine(data, refscan)
    coil_weight = conj(refscan[:, :, :, [1], :]) # refscan - only use one timepoint and one dim6
    # @show sum(.!isfinite.(coil_weight))
    # @show sum(.!isfinite.(data))
    combined = sum(data .* coil_weight; dims=5) # might need loop for RAM
    # @show sum(.!isfinite.(scaling(coil_weight)))
    combined .*= scaling(coil_weight)

    combined = dropdims(combined; dims=5)
    return combined
end

function scaling(coil_weight)
    scale = 1 ./ sum(abs.(coil_weight) .^ 2; dims=5)
    scale[.!isfinite.(scale)] .= 0
    return scale
end
