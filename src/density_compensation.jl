function density_compensation!(circles, sliceinfo, max_r, max_points)
    areas = calculate_effective_area_per_circle(sliceinfo, max_r)
    number_of_points_correction = points_on_circle_ratio(sliceinfo, max_points)

    @show dcf = (areas .* number_of_points_correction)

    circles .*= dcf

end

function calculate_effective_area_per_circle(sliceinfo, max_r)
    n_circles = length(sliceinfo)

    radii = [radius_normalized(first(circle_headers), sliceinfo) for circle_headers in sliceinfo]
    radii ./= max_r
    push!(radii, 2radii[end] - radii[end-1]) # extrapolate one radius

    areas_mid = [areas_to_middle_of_circles(radii, i) for i in 1:n_circles]

    effective_areas = areas_mid .- [0, areas_mid[1:end-1]...]

    return effective_areas
end

function areas_to_middle_of_circles(r, i)
    r_mid = (r[i] + r[i+1]) / 2
    return pi * r_mid^2
end

function points_on_circle_ratio(sliceinfo, max_points)
    points = [calculate_additional_info(circleinfo)[:points_on_circle] for circleinfo in sliceinfo]
    return max_points ./ points
end
