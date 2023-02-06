function kspace_coordinates(sliceinfo)
    n_grid = sliceinfo.twix[:n_frequency]
    fov_read = sliceinfo.twix[:fov_readout]
    coords = read_kspace_coordinates(sliceinfo)
    return normalize_kspace(coords, n_grid, fov_read)
end

function read_kspace_coordinates(slice_headers)
    coordinates = Complex[]
    for circle_headers in slice_headers
        s = calculate_additional_info(circle_headers)
        xy = read_kspace_coordinate(circle_headers[1])
        coords = construct_coordinates_circle(xy, s[:points_on_circle])
        append!(coordinates, coords)
    end
    return coordinates
end

function construct_coordinates_circle(xy, points_on_circle)
    angle_increment = 2pi / points_on_circle
    r = abs(xy)
    angle_first_point = angle(xy) - pi / 2
    phi = angle_first_point .- (0:points_on_circle-1) * angle_increment
    coords = r .* exp.(im .* phi)
    return coords
end

function normalize_kspace(coordinates, n_grid, fov_read)
    return coordinates ./ 2max_r(n_grid, fov_read)
end

function read_kspace_coordinate(head::ScanHeaderVD)
    x = parse_to_float(head.ice_param[1], head.ice_param[2])
    y = parse_to_float(head.ice_param[3], head.ice_param[4])
    return Complex(x, y)
end
parse_to_float(a::Int16, b::Int16) = a + b / 10_000

function radius_normalized(xy, n_grid, fov_read)
    r = abs(xy)
    delta_gm = 1e6 / (fov_read * GYRO_MAGNETIC_RATIO_OVER_TWO_PI)
    maxR = delta_gm * n_grid / 2
    return r / (2maxR)
end

function max_r(n_grid, fov_read)
    delta_gm = 1e6 / (fov_read * GYRO_MAGNETIC_RATIO_OVER_TWO_PI)
    max_r = delta_gm * n_grid / 2
    return max_r
end
