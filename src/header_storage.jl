"""
    headers = initialize_header_storage(info)

Returs empty Circle array of size [nPart][nCircles]
Each Circle can store all headers corresponding to that part and circle.
"""
function initialize_header_storage(info)
    return [[Circle(info) for _ in 1:info[:circles_per_part][i]] for i in 1:info[:n_part]]
end

function store!(header_array, header, info)
    circle = get_circle(header_array, header, info)
    store!(circle, header)
    return circle
end

function get_circle(header_array, header, info)
    return header_array[part_from_one(header, info)][header[:circle]]
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
