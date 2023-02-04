function reconstruct(filename; n_channels=1, fov_read)
    headers = read_data_headers(filename)
    n_grid = 2 * length(headers[:ONLINE][1])

    image = cat((reconstruct_slice(filename, h, n_channels, n_grid, fov_read) for h in headers[:ONLINE])...; dims=3)
    return image
end

function reconstruct_slice(filename, headers, n_channels, n_grid, fov_read)
    slice = open(filename) do io
        read_slice(io, headers, n_channels)
    end

    kspace_data = rearrange(slice, headers, n_channels)
    kspace_points = MRSI.kspace_coordinates(headers, n_grid, fov_read)
    return fourier_transform(kspace_data, kspace_points, n_grid)
end
