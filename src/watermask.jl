using Statistics, ImageFiltering

function watermask(dat_file)
    info = read_scan_info(dat_file)
    image = reconstruct(info[:PATREFSCAN]; combine=true)
    cut_to_ellipse!(image)
    summed = sqrt.(sum(abs.(image[:, :, :, 2:end, :]) .^ 2; dims=4:5))[:, :, :, 1, 1]

    mask_mrsi(summed; factor=1.5)
end

function mask_mrsi(weight::AbstractArray; factor=1, threshold=nothing)
    if threshold isa Nothing
        w = weight[:]
        q05, q15, q8, q99 = quantile.(Ref(w), (0.05, 0.15, 0.8, 0.99))
        high_intensity = mean(filter(isfinite, w[q8.<=w.<=q99]))
        noise = mean(filter(isfinite, w[w.<=q15]))
        if noise > high_intensity / 10
            noise = mean(filter(isfinite, w[w.<=q05]))
            if noise > high_intensity / 10
                noise = 0 # no noise detected
            end
        end
        threshold = max(5noise, high_intensity / 5)
    end
    mask = weight .> (threshold * factor)
    # remove small holes
    mask = imfilter(mask, Kernel.gaussian((1,1,1))) .> 0.5

    return mask
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