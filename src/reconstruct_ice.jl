function reconstruct_ice(filename, type=:ONLINE)
    info = extract_twix(read_twix_protocol(filename))
    info[:filename] = filename
    calculate_additional_info!(info)
    headers = read_scan_headers(filename, info[:n_channels])[type]

    ## Info cheating
    info[:n_fid] = MRSI.calculate_additional_info(ScanInfo(filename, :ONLINE)[1][2])[:n_fid] # fine, is available in ICE
    info[:radii] = [radius_normalized(headers[findfirst(h -> info[:circle_order][h.dims[LIN]] == i, headers)], info) for i in 1:info[:max_n_circles]] # problem
    info[:max_n_points_on_circle] = maximum([get_n_points_on_circle(h, info[:oversampling_factor]) for h in headers]) # TODO
    ##


    image = mmaped_image(info)

    read_and_reconstruct_image_per_circle!(image, headers, info)

    fft_slice_dim!(image)

    return image
end

function read_and_reconstruct_image_per_circle!(image, headers, info)
    header_array = [[Circle(info) for _ in 1:info[:circles_per_part][i]] for i in 1:info[:n_part]] # [part][circle]
    for head in headers
        c = store!(header_array, head, info)
        if is_complete(c)
            image[:, :, c[:part], :, :] .+= reconstruct_circle(c)
        end
    end
    return header_array
end

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
        c.headers = [ScanHeaderVD[] for _ in 1:h[:n_TI]]
    end
    push!(c.headers[h[:TI]], h)
end

function is_complete(c::Circle)
    required_length = c[:n_fid] * c[:n_points_on_circle] ÷ c[:n_TI]
    is_full(heads) = length(heads) * c[:n_adc_points] ≥ required_length
    return all(is_full.(c.headers))
end

# Returns (n_freq, n_phase, n_points, n_channels)
function reconstruct_circle(c::Circle)
    data = read_data(c)
    coords = construct_circle_coordinates(c)

    fov_shift!(data, coords, c)
    frequency_offset_correction!(data, c)

    density_compensation!(data, c)
    return fourier_transform(data, coords, c[:n_frequency])
end

function read_data(c::Circle)
    data = zeros(ComplexF64, c[:n_points_on_circle], c[:n_fid], c[:n_channels])
    open(c[:filename]) do io
        for (TI, heads) in enumerate(c.headers)
            data[:, TI:c[:n_TI]:end, :] .= read_reshape_one_TI(io, heads, c)
        end
    end
    return data
end

function read_reshape_one_TI(io, heads, info)
    data = read_one_TI(io, heads, info[:n_channels])
    data = @view data[1:info[:n_useful_adc_points], :]
    return reshape(data, info[:n_points_on_circle], :, info[:n_channels])
end

function read_one_TI(io, headers, n_channels)
    vcat((read_adc(io, head.data_position, n_channels) for head in headers)...)
end
