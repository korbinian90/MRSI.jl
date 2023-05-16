function density_compensation!(data, info)
    areas = calculate_effective_area_per_circle(info[:radii])
    number_of_points_correction = points_on_circle_ratio(info)

    dcf = areas[info[:circle]] * number_of_points_correction

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

# Uses calculated radii instead of values read from the headers
function density_compensation_ice!(data, info)
    radii = collect(range(1; step=2, length=info[:max_n_circles]))
    areas = calculate_effective_area_per_circle(radii)
    number_of_points_correction = 1 / info[:n_points_on_circle]

    dcf = areas[info[:circle]] * number_of_points_correction
    
    data .*= dcf
end
