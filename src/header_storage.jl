"""
    headers = initialize_header_storage(info)

Returs empty Circle array of size [nPart][nCircles]
Each Circle can store all headers corresponding to that part and circle.
"""
function initialize_header_storage(info)
    [[Circle(info) for _ in 1:info[:circles_per_part][i]] for i in 1:info[:n_part]]
end

function store!(header_array, header, info)
    circle = get_circle(header_array, header, info)
    store!(circle, header)
    return circle
end

function get_circle(header_array, header, info)
    header_array[part_from_one(header, info)][header[:circle]]
end

# store headers as [TI][adc][average]
function store!(c::Circle, h::ScanHeaderVD)
    if isnothing(c.headers)
        c.headers = [[ScanHeaderVD[]]]
    end
    while length(c.headers) < h[:TI]
        push!(c.headers, [ScanHeaderVD[]])
    end
    while length(c.headers[h[:TI]]) < h[:adc]
        push!(c.headers[h[:TI]], ScanHeaderVD[])
    end
    push!(c.headers[h[:TI]][h[:adc]], h)
end

function sort_into_circles(headers, info)
    headers_sorted_into_circles = initialize_header_storage(info)
    for head in headers
        store!(headers_sorted_into_circles, head, info)
    end
    return headers_sorted_into_circles
end
