function initialize_header_storage(info)
    return [[Circle(info) for _ in 1:info[:circles_per_part][i]] for i in 1:info[:n_part]]
end

function store!(header_array, header, info)
    circle = get_circle(header_array, header, info)
    store!(circle, header)
    return circle
end

function get_circle(header_array, header, info)
    num = header.dims[LIN]
    circle = info[:circle_order][num]
    part = info[:part_order][num]
    return header_array[part][circle]
end

function store!(c::Circle, h::ScanHeaderVD)
    if isnothing(c.headers)
        c.headers = [ScanHeaderVD[]]
    end
    while length(c.headers) < h[:TI]
        push!(c.headers, ScanHeaderVD[])
    end
    push!(c.headers[h[:TI]], h)
end

function is_complete(c::Circle)
    required_length = c[:n_fid] * c[:n_points_on_circle] ÷ c[:n_TI]
    is_full(heads) = length(heads) * c[:n_adc_points] ≥ required_length
    return all(is_full.(c.headers))
end
