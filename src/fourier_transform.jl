function fourier_transform(slice, kspace_points, n_grid)
    @assert length(kspace_points) == size(slice, 1)
    slice_matrix = reshape(slice, length(kspace_points), :)
    dft_matrix = ones(n_grid^2, length(kspace_points))
    image = dft_matrix * slice_matrix
    return reshape(image, n_grid, n_grid, size(slice)[2:end]...)
end