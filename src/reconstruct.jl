function reconstruct(filename, type=:ONLINE)
    scaninfo = ScanInfo(filename, type)

    info = MRSI.calculate_additional_info(first(first(scaninfo)))

    sz = (scaninfo.twix[:n_frequency], scaninfo.twix[:n_frequency], prod(size(scaninfo)), info[:n_fid], scaninfo.twix[:n_channels])
    image = write_emptynii(sz, tempname(); datatype=ComplexF32)

    for (i, si) in enumerate(scaninfo)
        image.raw[:, :, i, :, :] .= MRSI.reconstruct_slice(filename, si)
    end

    return image
end

function reconstruct_slice(filename, scaninfo)
    slice = open(filename) do io
        read_slice(io, scaninfo)
    end
    kspace_data = rearrange(slice, scaninfo)

    n_grid = scaninfo.twix[:n_frequency]
    kspace_points = MRSI.kspace_coordinates(scaninfo)

    return fourier_transform(kspace_data, kspace_points, n_grid)
end
