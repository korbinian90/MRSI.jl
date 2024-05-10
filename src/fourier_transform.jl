"""
    fourier_transform(kdata, kspace_coordinates, n_grid)
Fourier transform from kspace to a slice. Can be one circle or whole slice
"""
function fourier_transform(kdata::AbstractArray{T}, c::Circle) where T
    n_grid = c[:n_frequency]
    dft_matrix = get_dft_matrix(c, n_grid, T)
    data_matrix = reshape(kdata, size(dft_matrix, 2), :)
    image_matrix = dft_matrix * data_matrix
    image = reshape(image_matrix, n_grid, n_grid, size(kdata)[2:end]...) # one slice
    return image
end

"""
    get_dft_matrix(circle::Circle, n_grid)
Retrieves the DFT matrix for a given circle. Memoizes the matrix for the same circle
"""
function get_dft_matrix(circle::Circle, n_grid, datatype)
    xy = datatype(read_first_kspace_coordinate_normalized(circle))
    return calculate_dft_matrix(xy, circle[:n_points_on_circle], n_grid)    
end

@memoize function calculate_dft_matrix(xy::Complex, n_points_on_circle, n_grid)
    kspace_coordinates = construct_circle_coordinates(xy, n_points_on_circle)
    return calculate_dft_matrix(kspace_coordinates, n_grid)
end

function calculate_dft_matrix(kspace_coordinates::AbstractArray{T}, n_grid) where T
    r_range = n_grid / 2 - 0.5
    rx = ry = -r_range:r_range
    ry = reshape(ry, 1, :)
    ks = reshape(kspace_coordinates, 1, 1, :)
    f(rx, ry, k) = exp(T(im * -2pi * (rx * real(k) + ry * imag(k))))
    dft_matrix = f.(rx, ry, ks) # size: (n_grid, n_grid, n_kspace)
    return reshape(dft_matrix, n_grid * n_grid, :)
end

function fft_slice_dim!(image::AbstractArray{T,5}) where {T}
    for i in axes(image, 4), j in axes(image, 5)
        image[:, :, :, i, j] .= fft_slice_dim(image[:, :, :, i, j])
    end
end
fft_slice_dim(image; dims=3) = fftshift(fft(ifftshift(image, dims), dims), dims)
