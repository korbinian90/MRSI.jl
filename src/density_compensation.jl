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

# function SamplingDensityVoronoi(K,MaxRad)
#     K=K - mean(K,1);
#     if nargin==1
#         MaxRad=max(sqrt(K(:,1).^2+ K(:,2).^2));
#     end
#     [uK,ia,ic] = unique(K,'rows');
#     %Retrieve K by using command (K =) uK(ic),:
#     [v,c] = voronoin(double(uK),{'Qbb','QJ'}) ;
    
#     VArea = zeros(length(c),1) ;
    
#     for i = 1:length(c)
#         v1 = v(c{i},1)  ;
#         v2 = v(c{i},2) ;
        
#         VertOut=find((v1.^2+v2.^2)>(MaxRad^2));% finding vertices put outside the encoding disk (dealing with Voronoi domains of outter k-space points)
#         if ~isempty(VertOut)
#             v1(VertOut)=inf;v2(VertOut)=inf; % assign them to inf (later their area will be assigned to 1 )
#             %VertIn = find((v1.^2+v2.^2)<=(MaxRad^2));% finding vertices put outside the encoding disk (dealing with Voronoi domains of outter k-space points)
#             %v1 = v1(VertIn);
#             %v2 = v2(VertIn);
#         end
#         VRadius(i) = sqrt(mean(v1).^2+mean(v2).^2);
#         VArea(i) = polyarea(v1,v2) ; % if one of the vertex coordinate is inf, the area becomes NaN
#     end
#     %VArea = VArea /max(VArea(:)); This should not be done here for 3D
    
#     if sum(~isnan(VArea(:)))==0
#         maxVArea = 1; %all VArea are NaN, there is no max, set values to 1
#     else
#         maxVArea = max(VArea(~isnan(VArea(:))));
#     end
    
#     VArea(isnan(VArea(:))) = maxVArea; % basically assign value equivalent to the k-space center to all outer points for which domain computation was not possible
    
#     %  figure
#     %  hold on
#     %  voronoi(uK(:,1), uK(:,2))
#     %  figure
#     %  scatter(VRadius(ic),VArea(ic))
    
    
#     W = VArea(ic);
#     return W
# end    