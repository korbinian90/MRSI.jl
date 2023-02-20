function rearrange(slice_data, info)
    return [rearrange_circle(circle, info) for circle in slice_data]
end

function rearrange_circle(circle_data, info)
    # input dims: (n_adc_points, adc_line, n_temporal_interleaves, n_channels)
    adc_merged = merge_adc_lines(circle_data, info)
    cut = cut_unused_adc_points(adc_merged, info)
    circle = reshape_into_circles(cut, info)
    ordered = merge_temporal_interleaves(circle, info)
    # output dims: (n_points_on_circle, n_fid, n_channels)
    return ordered
end

merge_adc_lines(data, info) = reshape(data, info[:n_adc_points] * info[:n_adcs], info[:n_temporal_interleaves], info[:n_channels])
cut_unused_adc_points(data, info) = data[1:info[:n_useful_adc_points], :, :]
reshape_into_circles(data, info) = reshape(data, info[:n_points_on_circle], info[:n_fid] ÷ info[:n_temporal_interleaves], info[:n_temporal_interleaves], info[:n_channels])
function merge_temporal_interleaves(data, info)
    data = permutedims(data, (1, 3, 2, 4,))
    data = reshape(data, info[:n_points_on_circle], info[:n_fid], info[:n_channels])
    return data
end
