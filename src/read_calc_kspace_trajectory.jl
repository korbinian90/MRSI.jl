function kspace_coordinates(sliceinfo)
    n_grid = sliceinfo[:n_frequency]
    fov_read = sliceinfo[:fov_readout]
    coords = read_kspace_coordinates(sliceinfo)
    return normalize_kspace(coords, n_grid, fov_read)
end

function read_kspace_coordinates(slice_headers)
    coordinates = Complex[]
    for circle_headers in slice_headers
        s = calculate_additional_info(circle_headers)
        xy = read_first_kspace_coordinate(first(circle_headers))
        coords = construct_coordinates_circle(xy, s[:points_on_circle])
        append!(coordinates, coords)
    end
    return coordinates
end

function construct_coordinates_circle(c::CircleTI)
    xy = read_first_kspace_coordinate(first(c.headers))
    coords = construct_coordinates_circle(xy, c[:n_points_on_circle])
    return normalize_kspace(coords, c[:n_frequency], c[:fov_readout])
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

function read_first_kspace_coordinate(head::ScanHeaderVD)
    x = parse_to_float(head.ice_param[1], head.ice_param[2])
    y = parse_to_float(head.ice_param[3], head.ice_param[4])
    return Complex(x, y)
end
parse_to_float(a::Int16, b::Int16) = a + b / 10_000

function radius_normalized(xy, info)
    r = abs(xy)
    delta_gm = 1e6 / (info[:fov_readout] * GYRO_MAGNETIC_RATIO_OVER_TWO_PI)
    r_norm = delta_gm * info[:n_frequency]
    return r / r_norm
end
function radius_normalized(head::ScanHeaderVD, info)
    return radius_normalized(read_first_kspace_coordinate(head), info)
end

function max_r(n_grid, fov_read)
    delta_gm = 1e6 / (fov_read * GYRO_MAGNETIC_RATIO_OVER_TWO_PI)
    max_r = delta_gm * n_grid / 2
    return max_r
end
