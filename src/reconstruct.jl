function reconstruct(filename, type=:ONLINE)
    scaninfo = ScanInfo(filename, type)
    image = mmaped_image(scaninfo)

    for (i, sliceinfo) in enumerate(scaninfo)
        image[:, :, i, :, :] .= MRSI.reconstruct_slice(filename, sliceinfo)
    end

    # for i in axes(image, 4), j in axes(image, 5)
    #     image[:, :, :, i, j] .= fft3D(image[:, :, :, i, j])
    # end

    return image
end

function reconstruct_slice(filename, sliceinfo)
    kspace_data = read_rearrange_correct(filename, sliceinfo)
    kspace_points = MRSI.kspace_coordinates(sliceinfo)
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
