function reconstruct_ice(filename, type=:ONLINE)
    info = extract_twix(read_twix_protocol(filename))
    info[:filename] = filename
    headers = read_scan_headers(filename, info[:n_channels])[type]

    ## Info cheating
    info[:n_fid] = MRSI.calculate_additional_info(ScanInfo(filename, :ONLINE)[1][2])[:n_fid]
    ##
    calculate_additional_info!(info)

    image = mmaped_image(info)

    read_and_reconstruct_image_per_circle_ti!(image, headers, info)

    fft_slice_dim!(image)

    return image
end

function read_and_reconstruct_image_per_circle_ti!(image, headers, info)
    header_array = [[CircleTI[] for _ in 1:info[:circles_per_part][i]] for i in 1:info[:n_part]]
    for head in headers
        num = head.dims[LIN]
        circle = info[:circle_order][num]
        part = info[:part_order][num]
        TI = head[:TI]
        circle_ti_arr = header_array[part][circle]
        while length(circle_ti_arr) < TI
            push!(circle_ti_arr, CircleTI(info))
        end
        c = circle_ti_arr[TI]
        if add_header_and_test_if_last(c, head)
            image[:, :, part, c[:fid_points_for_TI], :] .+= reconstruct_circle(c)
        end
    end
end

function add_header_and_test_if_last(c::CircleTI, h::ScanHeaderVD)
    push!(c.headers, h)
    required_length = c[:n_fid] * c[:n_points_on_circle] ÷ c[:n_TI]
    return length(c.headers) * h[:n_adc_points] ≥ required_length
end

# Returns (n_freq, n_phase, n_points, n_channels)
function reconstruct_circle(c::CircleTI)
    data = read_data(c)
    data = data[begin:c[:n_useful_adc_points], :]
    data = reshape(data, c[:n_points_on_circle], :, c[:n_channels])
    coords = construct_circle_coordinates(c)

    fov_shift!(data, coords, c)

    return fourier_transform(data, coords, c[:n_frequency])
end

function read_data(c::CircleTI)
    open(c[:filename]) do io
        vcat((read_adc(io, head.data_position, c[:n_channels]) for head in c.headers)...)
    end
end
