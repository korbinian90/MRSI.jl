function reconstruct(filename, type=:ONLINE)
    scaninfo = ScanInfo(filename, type)
    image = mmaped_image(scaninfo)

    for (i, sliceinfo) in enumerate(scaninfo)
        image[:, :, i, :, :] .= MRSI.reconstruct_slice(filename, sliceinfo)
    end

    for i in axes(image, 4), j in axes(image, 5)
        image[:, :, :, i, j] .= fft3D(image[:, :, :, i, j])
    end

    return image
end

function mmaped_image(scaninfo, name=tempname())
    n_fid = MRSI.calculate_additional_info(first(first(scaninfo)))[:n_fid]
    n_slices = prod(size(scaninfo))
    sz = (scaninfo.twix[:n_frequency], scaninfo.twix[:n_frequency], n_slices, n_fid, scaninfo.twix[:n_channels])
    return mmap(name, Array{ComplexF64,5}, sz)
end

function reconstruct_slice(filename, scaninfo)
    kspace_data = read_rearrange_slice(filename, scaninfo)
    kspace_points = MRSI.kspace_coordinates(scaninfo)
    n_grid = scaninfo.twix[:n_frequency]

    fov_shift!(kspace_data, kspace_points, scaninfo)

    return fourier_transform(kspace_data, kspace_points, n_grid)
end

function read_rearrange_slice(filename, scaninfo)
    slice = open(filename) do io
        read_slice(io, scaninfo)
    end
    return rearrange(slice, scaninfo)
end
