function density_compensation!(data, info)
    areas = calculate_effective_area_per_circle(info)
    number_of_points_correction = points_on_circle_ratio(info)

    dcf = (areas .* number_of_points_correction)

    data .*= dcf[info[:circle]]
end

function calculate_effective_area_per_circle(info)
    radii = copy(info[:radii]) ./ maximum(info[:radii])
    push!(radii, 2radii[end] - radii[end-1]) # extrapolate one radius

    areas_mid = [areas_to_middle_of_circles(radii, i) for i in 1:length(info[:radii])]
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
