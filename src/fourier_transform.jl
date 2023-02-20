"""
    fourier_transform(kdata, kspace_coordinates, n_grid)
Fourier transform from kspace to a slice. Can be one circle or whole slice

"""
function fourier_transform(kdata, kspace_coordinates, n_grid)
    @assert length(kspace_coordinates) == size(kdata, 1)
    data_matrix = reshape(kdata, length(kspace_coordinates), :)
    dft_matrix = calculate_dft_matrix(kspace_coordinates, n_grid)
    image_matrix = dft_matrix * data_matrix
    image = reshape(image_matrix, n_grid, n_grid, size(kdata)[2:end]...) # one slice
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

function fft_slice_dim!(image::AbstractArray{T,5}) where {T}
    for i in axes(image, 4), j in axes(image, 5)
        image[:, :, :, i, j] .= fft_slice_dim(image[:, :, :, i, j])
    end
end
fft_slice_dim(image; dims=3) = fftshift(fft(ifftshift(image, dims), dims), dims)
