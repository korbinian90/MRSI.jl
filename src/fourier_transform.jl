function fourier_transform(slice, kspace_coordinates, n_grid)
    @assert length(kspace_coordinates) == size(slice, 1)
    slice_matrix = reshape(slice, length(kspace_coordinates), :)
    dft_matrix = calculate_dft_matrix(kspace_coordinates, n_grid)
    image_matrix = dft_matrix * slice_matrix
    image = reshape(image_matrix, n_grid, n_grid, size(slice)[2:end]...) # one slice
    return image
end

function calculate_dft_matrix(kspace_coordinates, n_grid)
    r_range = n_grid / 2 - 0.5
    rx = ry = -r_range:r_range
    ry = reshape(ry, 1, :)
    ks = reshape(kspace_coordinates, 1, 1, :)
    f(rx, ry, k) = exp(im * -2pi * (rx * real(k) + ry * imag(k)))
    dft_matrix = f.(rx, ry, ks) # size: (n_grid, n_grid, n_kspace)
    return reshape(dft_matrix, n_grid * n_grid, :)
end

fft3D(image; dims=3) = fftshift(fft(ifftshift(image, dims), dims), dims)
