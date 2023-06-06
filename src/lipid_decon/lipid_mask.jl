function lipid_mask_from_mrsi(data, info; fat_range_ppm=[1.8, 0.5], threshold=0.1)
    dims = 4
    spectrum = reverse(fftshift(fft(data, dims), dims));
    fat_range = ppm_to_vecsize_point.(Ref(info), fat_range_ppm)
    fat = sum(abs.(selectdim(spectrum, dims, fat_range[1]:fat_range[2])); dims)
    all = sum(abs.(spectrum); dims)
    return (fat ./ all) .> threshold
end

function ppm_to_vecsize_point(info, ppm, center_around_ppm=4.65)
    bandwidth_frequency = 1e9 / info[:dwelltime]
    step_frequency = bandwidth_frequency / info[:n_fid]
    offset = ceil(info[:n_fid] / 2) # k-space definition, where for even matrix sizes, k-space is also not centered around 0, but is (-N/2):1:(N/2-1)

    vecsize_point = -1e-6 * (ppm - center_around_ppm) * info[:larmor_frequency] / step_frequency + offset
    return round(Int, vecsize_point)
end
