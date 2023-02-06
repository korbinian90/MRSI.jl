function rearrange(slice_data, scaninfo)
    return vcat((rearrange_circle(d, i) for (d, i) in zip(slice_data, scaninfo))...)
end

function rearrange_circle(circle_data, scaninfo)
    s = calculate_additional_info(scaninfo)
    # input dims: (adc_points, adc_line, temporal_interleaves, n_channels)
    adc_merged = merge_adc_lines(circle_data, s)
    cut = cut_unused_adc_points(adc_merged, s)
    circle = reshape_into_circles(cut, s)
    ordered = merge_temporal_interleaves(circle, s)
    # output dims: (points_on_circle, n_fid, n_channels)
    return ordered
end

merge_adc_lines(data, s) = reshape(data, s[:adc_points] * s[:adcs], s[:temporal_interleaves], s[:n_channels])
cut_unused_adc_points(data, s) = data[1:s[:useful_adc_points], :, :, :]
reshape_into_circles(data, s) = reshape(data, s[:points_on_circle], s[:n_fid] ÷ s[:temporal_interleaves], s[:temporal_interleaves], s[:n_channels])
function merge_temporal_interleaves(data, s)
    data = permutedims(data, (1, 3, 2, 4,))
    data = reshape(data, s[:points_on_circle], s[:n_fid], s[:n_channels])
    return data
end
