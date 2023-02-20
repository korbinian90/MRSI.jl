function reconstruct(filename, type=:ONLINE)
    info = extract_twix(read_twix_protocol(filename))
    info[:filename] = filename
    headers = read_scan_headers(info)[type]
    calculate_complete_info!(info, headers)
    headers = rearrange_headers(headers, info)

    image = mmaped_image(info)

    for (i, slice_headers) in enumerate(headers)
        image[:, :, i, :, :] .= MRSI.reconstruct_slice(slice_headers, info)
    end

    fft_slice_dim!(image)

    return image
end

function reconstruct_slice(slice_headers, info)
    kspace_data = read_rearrange_correct(slice_headers, info)
    kspace_points = read_kspace_coordinates(info)
    n_grid = info[:n_frequency]

    fov_shift!(kspace_data, kspace_points, info)

    return fourier_transform(kspace_data, kspace_points, n_grid)
end

function read_rearrange_correct(slice_headers, info)
    slice = open(info[:filename]) do io
        read_slice(io, slice_headers, info)
    end
    circles = rearrange(slice, info)

    # @show max_r = maximum_radius(info)
    # max_points = 120
    density_compensation!(circles, info)

    return vcat(circles...)
end
