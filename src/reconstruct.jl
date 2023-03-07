"""
    reconstruct(filename, type=:ONLINE)
    
Reconstructs a SIEMENS dat file
"""
function reconstruct(filename, type=:ONLINE)
    data_headers, info = read_scan_info(filename, type)
    image = mmaped_image(info)

    read_and_reconstruct_image_per_circle!(image, data_headers, info)
    fft_slice_dim!(image)

    return image
end

function read_and_reconstruct_image_per_circle!(image, headers, info)
    headers_sorted_into_circles = initialize_header_storage(info)
    for head in headers
        circle = store!(headers_sorted_into_circles, head, info)
        if is_complete(circle)
            image[:, :, circle[:part], :, :] .+= reconstruct(circle)
        end
    end
end

# Returns [n_freq, n_phase, n_points, n_channels]
function reconstruct(c::Circle)
    kspace_coordinates = construct_circle_coordinates(c)
    kdata = read_data(c)

    fov_shift!(kdata, kspace_coordinates, c)
    frequency_offset_correction!(kdata, c)
    density_compensation!(kdata, c)
    image = fourier_transform(kdata, kspace_coordinates, c[:n_frequency])

    return image
end
