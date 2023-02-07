function reconstruct(filename, type=:ONLINE)
    scaninfo = ScanInfo(filename, type)

    info = MRSI.calculate_additional_info(first(first(scaninfo)))

    sz = (scaninfo.twix[:n_frequency], scaninfo.twix[:n_frequency], prod(size(scaninfo)), info[:n_fid], scaninfo.twix[:n_channels])
    image = mmap(tempname(), Array{ComplexF32,5}, sz)

    for (i, si) in enumerate(scaninfo)
        image[:, :, i, :, :] .= MRSI.reconstruct_slice(filename, si)
    end

    # image = fft_partitions!(image)

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
