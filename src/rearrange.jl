function rearrange(circle_data, headers, n_channels)
    s = calculate_additional_data(headers, n_channels)
    # @show size(circle_data)
    permuted1 = permute_adc_line_after_adc_points(circle_data, s)
    # @show size(permuted1)
    # @show s
    cut = cut_unused_adc_points(permuted1, s)
    # @show size(cut)
    circle = reshape_into_circles(cut, s)
    # @show size(circle)
    permuted2 = permutedims(circle, (3, 1, 2, 4, 5))
    # @show size(permuted2)
    return permuted2
end

function permute_adc_line_after_adc_points(data, s)
    data = permutedims(data, (1, 2, 4, 3, 5))
    data = reshape(data, s[:temporal_interleaves], s[:part], s[:adc_points] * s[:adcs], s[:n_channels])
    return data
end
cut_unused_adc_points(circle, s) = circle[:, :, 1:s[:useful_adc_points], :]
reshape_into_circles(circle, s) = reshape(circle, s[:temporal_interleaves], s[:part], s[:points_on_circle], s[:n_fid], s[:n_channels])
