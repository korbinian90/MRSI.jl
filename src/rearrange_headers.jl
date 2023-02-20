"""Returns headers as nested array [slice, part][circle][adc,TI]"""
function rearrange_headers(headers::AbstractArray{ScanHeaderVD}, info)
    slices = unique(h.dims[SLI] for h in headers)
    rearranged = Array{Any}(undef, maximum(slices), info[:n_part])
    for slice in slices
        slice_headers = [h for h in headers if h.dims[SLI] == slice]
        rearranged[slice, :] .= rearrange_headers_slice(slice_headers, info[:n_part], info[:max_n_circles])
    end
    return rearranged
end

function rearrange_headers_slice(slice_headers, n_part, max_n_circles)
    circles = unique(h.dims[LIN] for h in slice_headers)
    (part_order, circle_order, circles_per_part) = calculate_part_order(n_part, max_n_circles)

    rearranged = [Array{Any}(undef, n) for n in circles_per_part]
    for num in circles # oldADC 3D
        circle_of_part = circle_order[num]
        part = part_order[num]

        circle_headers = [h for h in slice_headers if h.dims[LIN] == num]
        rearranged[part][circle_of_part] = order_adcs_and_temporal_interleaves(circle_headers)
    end
    return rearranged
end

function order_adcs_and_temporal_interleaves(circle_headers)
    n_temporal_interleaves = maximum(h.dims[IDB] for h in circle_headers)
    n_adcs = maximum(h.dims[IDA] for h in circle_headers)

    ordered = Array{ScanHeaderVD}(undef, n_adcs, n_temporal_interleaves)
    for head in circle_headers
        TI = head.dims[IDB]
        adc = head.dims[IDA]
        ordered[adc, TI] = head
    end
    return ordered
end

function calculate_n_circles_per_part(part_max, max_n_circles, i)
    if part_max == 0 # 2D
        return max_n_circles
    end
    return max(2, Int(ceil(sqrt(max_n_circles^2 - (i * max_n_circles / part_max)^2))))
end

function calculate_part_order(n_part, max_n_circles)
    part_max = n_part ÷ 2
    circles_per_part = [calculate_n_circles_per_part(part_max, max_n_circles, i) for i in -part_max:part_max]
    part_order = vcat((repeat([i], n) for (i, n) in enumerate(circles_per_part))...)
    circle_order = vcat((1:c for c in circles_per_part)...)
    return (part_order, circle_order, circles_per_part)
end
