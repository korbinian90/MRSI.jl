function reconstruct(filename, type=:ONLINE)
    data_headers, info = read_scan_info(filename, type)

    image = mmaped_image(info)
    read_and_reconstruct_image_per_circle!(image, data_headers, info)
    fft_slice_dim!(image)

    return image
end

function read_scan_info(filename, type)
    checkVD(filename) || error("only implemented for VD files")
    info = extract_twix(read_twix_protocol(filename))
    info[:filename] = filename
    headers = read_scan_headers(info)[type]

    ## Info cheating (obtained from all headers instead of streamed like in ICE)
    info[:n_fid] = calculate_fid(info, headers) # fine, is available in ICE
    info[:radii] = [radius_normalized(headers[findfirst(h -> info[:circle_order][h.dims[LIN]] == i, headers)], info) for i in 1:info[:max_n_circles]] # problem, required for dcf
    info[:max_n_points_on_circle] = maximum([get_n_points_on_circle(h, info[:oversampling_factor]) for h in headers]) # can be corrected afterwards
    ##
    return headers, info
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

## Functions to store headers into circles
initialize_header_storage(info) = [[Circle(info) for _ in 1:info[:circles_per_part][i]] for i in 1:info[:n_part]]

function store!(header_array, head, info)
    circle = get_circle(header_array, head, info)
    store!(circle, head)
    return circle
end

function get_circle(header_array, head, info)
    num = head.dims[LIN]
    circle = info[:circle_order][num]
    part = info[:part_order][num]
    return header_array[part][circle]
end

function store!(c::Circle, h::ScanHeaderVD)
    if isnothing(c.headers)
        c.headers = [ScanHeaderVD[]]
    end
    while length(c.headers) < h[:TI]
        push!(c.headers, ScanHeaderVD[])
    end
    push!(c.headers[h[:TI]], h)
end

function is_complete(c::Circle)
    required_length = c[:n_fid] * c[:n_points_on_circle] ÷ c[:n_TI]
    is_full(heads) = length(heads) * c[:n_adc_points] ≥ required_length
    return all(is_full.(c.headers))
end
