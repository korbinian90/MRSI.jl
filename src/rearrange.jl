function rearrange(slice_data, scaninfo)
    return vcat((rearrange_circle(slice_data[i], scaninfo[i]) for i in eachindex(slice_data))...)
end

function rearrange_circle(circle_data, scaninfo)
    n_channels = read_n_channels(scaninfo)
    # dims: (adc_points, adc_line, part, TI, channels)
    s = calculate_additional_data(scaninfo)
    s[:n_channels] = n_channels
    # @show size(circle_data)
    adc_merged = merge_adc_lines(circle_data, s)
    # @show size(adc_merged)
    # @show s
    cut = cut_unused_adc_points(adc_merged, s)
    # @show size(cut)
    circle = reshape_into_circles(cut, s)
    # @show size(circle)
    ordered = merge_temporal_interleaves(circle, s)
    # @show size(ordered)
    return ordered
end

merge_adc_lines(data, s) = reshape(data, s[:adc_points] * s[:adcs], s[:temporal_interleaves], s[:part], s[:n_channels])
cut_unused_adc_points(data, s) = data[1:s[:useful_adc_points], :, :, :]
reshape_into_circles(data, s) = reshape(data, s[:points_on_circle], s[:n_fid] ÷ s[:temporal_interleaves], s[:temporal_interleaves], s[:part], s[:n_channels])
function merge_temporal_interleaves(data, s)
    data = permutedims(data, (1, 3, 2, 4, 5))
    data = reshape(data, s[:points_on_circle], s[:n_fid], s[:part], s[:n_channels])
    return data
end
