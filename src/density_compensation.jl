function density_compensation!(data, info; do_hamming_filter=true, ice=false)
    areas = calculate_effective_area_per_circle(get_radii(info; ice))

    dcf = areas[info[:circle]] * number_of_points_correction(info; ice)

    if do_hamming_filter
        dcf *= hamming_filter(radius_normalized(info))
    end

    data .*= dcf
end

function calculate_effective_area_per_circle(radii)
    n = length(radii)
    radii = radii ./ maximum(radii)
    push!(radii, 2radii[end] - radii[end-1]) # extrapolate one radius

    areas_mid = [areas_to_middle_of_circles(radii, i) for i in 1:n]
    effective_areas = areas_mid .- [0, areas_mid[1:end-1]...]

    return effective_areas
end

function areas_to_middle_of_circles(r, i)
    r_mid = (r[i] + r[i+1]) / 2
    return pi * r_mid^2
end

function points_on_circle_ratio(c::Circle)
    return c[:max_n_points_on_circle] ./ c[:n_points_on_circle]
end

# requires normalized radius
hamming_filter(r) = 0.54 + 0.46 * cos(2π * r)

function get_radii(info; ice=false)
    if ice
        collect(range(1; step=2, length=info[:max_n_circles]))
    else
        info[:radii]
    end
end

function number_of_points_correction(info; ice=false)
    if ice
        1 / info[:n_points_on_circle]
    else
        points_on_circle_ratio(info)
    end
end
