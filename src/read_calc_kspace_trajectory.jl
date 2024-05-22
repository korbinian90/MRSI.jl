function construct_circle_coordinates(c::Circle, gradient_delay_us)
    xy = read_first_kspace_coordinate_normalized(c)
    gradient_delay_rad = gradient_delay_to_rad(gradient_delay_us, c)
    return construct_circle_coordinates(xy, c[:n_points_on_circle], gradient_delay_rad)
end

function construct_circle_coordinates(xy::Complex, n, gradient_delay_rad)
    angle_increment = 2pi / n
    r = abs(xy)
    angle_first_point = -angle(xy) - pi / 2 # identical to matlab
    phi = angle_first_point .- (0:n-1) * angle_increment
    coords(phi) = r * Complex(real(exp(im * (phi + real(gradient_delay_rad)))), imag(exp(im * (phi + imag(gradient_delay_rad)))))
    return coords.(phi)
end

function read_first_kspace_coordinate_normalized(c::Circle)
    return read_first_kspace_coordinate_normalized(c[:first_head], c)
end

function read_first_kspace_coordinate_normalized(head::ScanHeaderVD, info)
    x = parse_to_float(head.ice_param[1], head.ice_param[2])
    y = parse_to_float(head.ice_param[3], head.ice_param[4])
    return Complex(x, y) ./ 2info[:max_r]
end
parse_to_float(a::Int16, b::Int16) = a + b / 10_000

radius_normalized(c::Circle) = abs(read_first_kspace_coordinate_normalized(c))
radius_normalized(head::ScanHeaderVD, info) = abs(read_first_kspace_coordinate_normalized(head, info))

function gradient_delay_to_rad(gradient_delay_us, c::Circle)
    gradient_delay_us
    return -gradient_delay_us[min(c[:n_TI], 3)] * 10^-6 * 2pi * (1E9 / c[:dwelltime]) / c[:n_TI]
end
