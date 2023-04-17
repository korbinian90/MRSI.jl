"""
    reconstruct(filename, type=:ONLINE)
    
Reconstructs a SIEMENS dat file
"""
function reconstruct(filename, type=:ONLINE; datatype=ComplexF32, kw...)
    data_headers, info = read_scan_info(filename, type)
    image = mmaped_image(info, datatype)

    read_and_reconstruct_image_per_circle!(image, data_headers, info; kw...)
    fft_slice_dim!(image)

    return image
end

function read_and_reconstruct_image_per_circle!(image, headers, info; kw...)
    headers_sorted_into_circles = initialize_header_storage(info)
    for head in headers
        circle = store!(headers_sorted_into_circles, head, info)
        if is_complete(circle)
            image[:, :, circle[:part], :, :] .+= reconstruct(circle; kw...)
        end
    end
end

# Returns [n_freq, n_phase, n_points, n_channels]
function reconstruct(c::Circle; ice=false)
    kspace_coordinates = construct_circle_coordinates(c)
    kdata = read_data(c)

    fov_shift!(kdata, kspace_coordinates, c)
    kdata = conj.(kdata)
    frequency_offset_correction!(kdata, c)
    
    if ice
        density_compensation_ice!(kdata, c)
    else
        density_compensation!(kdata, c)
    end
    image = fourier_transform(kdata, kspace_coordinates, c[:n_frequency])

    return image
end

function full_reconstruct(file; kw...)
    image = reconstruct(file; kw...)
    refscan = reconstruct(file, :PATREFSCAN; kw...)
    combined = coil_combine(image, refscan)

    combined = reverse(combined; dims=1) # LR flip

    return combined
end
