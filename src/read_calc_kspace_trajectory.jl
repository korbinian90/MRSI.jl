function read_kspace_coordinates(slice_headers)
    return vcat((construct_circle_coordinates(first(c), c) for c in slice_headers)...)
end

construct_circle_coordinates(c::Circle) = construct_circle_coordinates(first(first(c.headers)), c)
function construct_circle_coordinates(header, info)
    xy = read_first_kspace_coordinate_normalized(header, info)
    coords = construct_circle_coordinates(xy, info[:n_points_on_circle])
    return coords
end

function construct_circle_coordinates(xy::Number, points_on_circle::Number)
    angle_increment = 2pi / points_on_circle
    r = abs(xy)
    angle_first_point = angle(xy) - pi / 2
    phi = angle_first_point .- (0:points_on_circle-1) * angle_increment
    coords = r .* exp.(im .* phi)
    return coords
end

function read_first_kspace_coordinate_normalized(head::ScanHeaderVD, info)
    x = parse_to_float(head.ice_param[1], head.ice_param[2])
    y = parse_to_float(head.ice_param[3], head.ice_param[4])
    return Complex(x, y) ./ 2info[:max_r]
end
parse_to_float(a::Int16, b::Int16) = a + b / 10_000

radius_normalized(head::ScanHeaderVD, info) = abs(read_first_kspace_coordinate_normalized(head, info))
