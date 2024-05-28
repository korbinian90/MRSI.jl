# Idea to improve: take maximum instead of sum for fat
function lipid_mask_from_mrsi(data, info; fat_range_ppm=[1.8, 0.5], threshold=0.1)
    dims = 4
    spectrum = reverse(fftshift(fft(data, dims), dims))
    fat_range = ppm_to_vecsize_point.(Ref(info), fat_range_ppm)
    fat = sum(abs.(selectdim(spectrum, dims, fat_range[1]:fat_range[2])); dims)
    all = sum(abs.(spectrum); dims)
    return (fat ./ all) .> threshold
end

# Problem: requires brain mask
function lipid_mask_from_mrsi_pipeline(csi, brain_mask; fat_range_ppm=[1.8, 0.5], threshold=0.1)
    dims = 4
    spectrum = fftshift(fft(csi, dims), dims)

    f_range = size(csi, 4)
    f_range_start = ceil(f_range * 0.82)
    f_range_end = ceil(f_range * 0.96)

    for I in eachindex(spectrum)
        ftspectramax[I] = max(spectrum[I, f_range_start:f_range_end])
    end

    lipid_mask = ftspectramax .>= 0.3 * maximum(ftspectramax)

    tumor_mask = imerode(brain_mask, strel_disk(1))
    lipid_mask = lipid_mask .- tumor_mask .= 1

    corner_mask = .!imdilate(brain_mask, strel_disk(6))
    lipid_mask = lipid_mask .- corner_mask .= 1

    lipid_mask_total = lipid_mask_total + lipid_mask
    return lipid_mask_total
end

function strel_disk(r)
    struct_elm = ones(Int, 2 * r + 1, 2 * r + 1)
    for i in 1:(2*r+1)
        for j in 1:(2*r+1)
            if ((i - r - 1)^2 + (j - r - 1)^2) > r^2
                struct_elm[i, j] = 0
            end
        end
    end
    return struct_elm
end

function ppm_to_vecsize_point(info, ppm, center_around_ppm=4.65)
    bandwidth_frequency = 1e9 / info[:dwelltime]
    step_frequency = bandwidth_frequency / info[:n_fid]
    offset = ceil(info[:n_fid] / 2) # k-space definition, where for even matrix sizes, k-space is also not centered around 0, but is (-N/2):1:(N/2-1)

    vecsize_point = -1e-6 * (ppm - center_around_ppm) * info[:larmor_frequency] / step_frequency + offset
    return round(Int, vecsize_point)
end

function get_masks(csi, info; brain_mask=nothing, lipid_mask=nothing, kw...)
    sz = size(csi)[1:3]
    brain = read_raw(brain_mask, sz) .!= 0
    lipid = if lipid_mask == :from_spectrum
        lipid_mask_from_mrsi(csi, info)
    else
        read_raw(lipid_mask, sz) .!= 0
    end
    return brain, lipid
end

function read_raw(fn, sz; datatype=Float32)
    raw = ones(datatype, sz)
    if !isnothing(fn)
        read!(fn, raw)
    end
    return raw
end
