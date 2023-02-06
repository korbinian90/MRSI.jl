function reconstruct(filename, type=:ONLINE)
    scaninfo = ScanInfo(filename, type)
    image = cat((reconstruct_slice(filename, si) for si in scaninfo)...; dims=3)
    return image
end

function reconstruct_slice(filename, scaninfo)
    n_grid = scaninfo.twix[:n_frequency]

    slice = open(filename) do io
        read_slice(io, scaninfo)
    end

    kspace_data = rearrange(slice, scaninfo)
    kspace_points = MRSI.kspace_coordinates(scaninfo)
    return fourier_transform(kspace_data, kspace_points, n_grid)
end
