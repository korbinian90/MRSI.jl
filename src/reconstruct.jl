function reconstruct(filename, type=:ONLINE)
    scaninfo = ScanInfo(filename, type)
    scaninfo.info[:n_fid] = MRSI.calculate_additional_info(scaninfo[1][1])[:n_fid]

    image = mmaped_image(scaninfo)

    for (i, sliceinfo) in enumerate(scaninfo)
        image[:, :, i, :, :] .= MRSI.reconstruct_slice(filename, sliceinfo)
    end

    fft_slice_dim!(image)

    return image
end

function reconstruct_slice(filename, sliceinfo)
    kspace_data = read_rearrange_correct(filename, sliceinfo)
    kspace_points = read_kspace_coordinates(sliceinfo)
    n_grid = sliceinfo[:n_frequency]

    fov_shift!(kspace_data, kspace_points, sliceinfo)

    return fourier_transform(kspace_data, kspace_points, n_grid)
end

function read_rearrange_correct(filename, sliceinfo)
    slice = open(filename) do io
        read_slice(io, sliceinfo)
    end
    circles = rearrange(slice, sliceinfo)

    # @show max_r = maximum_radius(sliceinfo)
    # max_points = 120
    # density_compensation!(circles, sliceinfo, max_r, max_points)

    return vcat(circles...)
end
