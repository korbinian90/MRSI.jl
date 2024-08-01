using Statistics, ImageFiltering

function brain_mask_from_csi(csi, water, info)
    lipid_mask = lipid_mask_from_mrsi(csi, info)
    water_mask = intensity_mask(water)
    combined_mask = water_mask .& .!lipid_mask
    smoothed_mask = imfilter(combined_mask, Kernel.gaussian((3, 3, 3))) .> 0.5
    return lipid_mask, water_mask, combined_mask, smoothed_mask
end

function intensity_mask(dat_file)
    info = read_scan_info(dat_file)
    image = reconstruct(info[:PATREFSCAN])

    return intensity_mask(image)
end

function intensity_mask(image::AbstractArray)
    summed = sqrt.(sum(abs.(image[:, :, :, 2:end, :]) .^ 2; dims=4:5))[:, :, :]
    cut_to_ellipse!(summed)
    mask_mrsi(summed; factor=1.5)
end

function mask_mrsi(weight::AbstractArray; factor=1, threshold=nothing)
    if threshold isa Nothing
        threshold = determine_threshold(weight)
    end
    mask = weight .> (threshold * factor)
    # remove small holes
    mask = imfilter(mask, Kernel.gaussian((1.5, 1.5, 1.5))) .> 0.5

    return mask
end

function determine_threshold(weight)
    w = weight[isfinite.(weight) .& weight .> 0]
    q05, q15, q80, q99 = quantile.(Ref(w), (0.05, 0.15, 0.8, 0.99))
    high_intensity = mean(w[q80.<=w.<=q99])
    noise = mean(w[w.<=q15])
    if noise > high_intensity / 10
        noise = mean(w[w.<=q05])
        if noise > high_intensity / 10
            noise = 0 # no noise detected
        end
    end
    threshold = max(5noise, high_intensity / 5)

    return threshold
end

function cut_to_ellipse!(image)
    len = size(image)[1:2]
    center = (len .+ 1) ./ 2
    (a, b) = len ./ 2
    for I in CartesianIndices(image)
        x = I[1] .- center[1]
        y = I[2] .- center[2]
        if x^2 / a^2 + y^2 / b^2 > 1
            image[I] = 0
        end
    end
end