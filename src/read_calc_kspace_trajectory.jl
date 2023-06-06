construct_circle_coordinates(c::Circle) = construct_circle_coordinates(first(first(first(c.headers))), c)
function construct_circle_coordinates(header::ScanHeaderVD, info)
    xy = read_first_kspace_coordinate_normalized(header, info)
    return construct_circle_coordinates(xy, info)
end

function construct_circle_coordinates(xy::Complex, info)
    n = info[:n_points_on_circle]
    angle_increment = 2pi / n
    r = abs(xy)
    angle_first_point = -angle(xy) - pi / 2 # identical to matlab
    phi = angle_first_point .- (0:n-1) * angle_increment
    coords = r .* exp.(im .* phi)
    return coords
end

function read_first_kspace_coordinate_normalized(head::ScanHeaderVD, info)
    x = parse_to_float(head.ice_param[1], head.ice_param[2])
    y = parse_to_float(head.ice_param[3], head.ice_param[4])
    return Complex(x, y) ./ 2info[:max_r]
end
parse_to_float(a::Int16, b::Int16) = a + b / 10_000

radius_normalized(c::Circle) = radius_normalized(first(first(c.headers)), c)
radius_normalized(head::ScanHeaderVD, info) = abs(read_first_kspace_coordinate_normalized(head, info))
