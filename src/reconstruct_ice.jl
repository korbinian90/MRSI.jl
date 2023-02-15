function reconstruct_ice(filename, type=:ONLINE)
    info = extract_twix(read_twix_protocol(filename))
    info[:filename] = filename

    headers = read_scan_headers(filename, info[:n_channels])[type]
    # set_info!(info, headers)
    # headers = rearrange_headers(headers, info[:n_part], info[:max_n_circles])
    # ScanInfo(headers, info)

    ## Info cheating
    info[:n_fid] = MRSI.calculate_additional_info(ScanInfo(filename, :ONLINE)[1][2])[:n_fid]
    # check max_r


    max_n_circles = info[:n_frequency] ÷ 2
    # image = mmaped_image_ice(info)
    image = zeros(ComplexF32, info[:n_frequency], info[:n_frequency], info[:n_part], info[:n_fid], info[:n_channels])

    (part_order, circle_order, circles_per_part) = calculate_part_order(info[:n_part], max_n_circles)
    # n_adcs = maximum(h.dims[IDA] for h in headers)

    header_array = [[CircleTI[] for j in 1:circles_per_part[i]] for i in 1:info[:n_part]]

    # n_fid = round(Int, s[:adc_points] * s[:n_adcs] / s[:n_points_on_circle] * s[:temporal_interleaves] - 0.5)


    for head in headers
        num = head.dims[LIN]
        circle = circle_order[num]
        part = part_order[num]
        TI = head[:TI]
        circle_ti_arr = header_array[part][circle]
        while length(circle_ti_arr) < TI
            push!(circle_ti_arr, CircleTI(info))
        end
        c = circle_ti_arr[TI]
        if add_header_and_test_if_last(c, head)
            slice = reconstruct_circle(c) # (n_freq, n_phase, n_points, n_channels)
            @show fid_points = c[:fid_points_for_TI]
            image[:, :, part, c[:fid_points_for_TI], :] .+= slice
        end
    end

    # for i in axes(image, 4), j in axes(image, 5)
    #     image[:, :, :, i, j] .= fft3D(image[:, :, :, i, j])
    # end

    return image
end

function Base.getindex(h::ScanHeaderVD, s::Symbol)
    if s == :n_adc_points
        h.dims[COL]
    elseif s == :n_TI
        @show h.dims[IDD] - 1
    elseif s == :TI
        h.dims[IDB]
    elseif s == :adc
        h.dims[IDA]
    end
end

function Base.getindex(c::CircleTI, s::Symbol)
    if haskey(c.info, s) # look in info
        c.info[s]
    elseif s == :n_adcs # calculate here
        maximum(h.dims[IDA] for h in c.headers)
    elseif s == :n_points_on_circle
        get_n_points_on_circle(first(c.headers), c.info[:oversampling_factor])
    elseif s == :n_fid
        get_n_fid(c[:n_adc_points], c[:n_adcs], c[:n_TI], c[:n_points_on_circle] - 0.5) # why 0.5?
    elseif s == :fid_points_for_TI
        c[:TI]:c[:n_TI]:c[:n_fid]
    elseif s == :n_useful_adc_points
        (c[:n_fid] * c[:n_points_on_circle]) ÷ c[:n_TI]
    else # look in ScanHeaderVD
        first(c.headers)[s]
    end
end


# get_n_adc_points(h::ScanHeaderVD) = h.dims[COL]
# get_n_TI(h::ScanHeaderVD) = h.dims[IDD]
# get_TI(h::ScanHeaderVD) = h.dims[IDB]
# function get_n_fid(c::CircleTI)
#     h = first(c.headers)
#     n_adcs = maximum(h.dims[IDA] for h in c.headers)
#     return get_n_fid(get_n_adc_points(h), n_adcs, get_n_TI(h), get_n_points_on_circle(h))
# end
get_n_fid(n_adc_points, n_adcs, n_TI, n_points_on_circle) = round(Int, n_adc_points * n_adcs * n_TI / n_points_on_circle - 0.5) # why 0.5?
function get_n_points_on_circle(h::ScanHeaderVD, oversampling_factor)
    return @show round(Int, max(h.dims[IDC] - 1, h.ice_param[6]) * oversampling_factor)
end

# get_fid_points_for_TI(c::CircleTI) = get_fid_points_for_TI(get_TI(c), get_n_TI(c), get_n_fid(c))
# get_fid_points_for_TI(TI, n_TI, n_fid) = TI:n_TI:n_fid
function add_header_and_test_if_last(c::CircleTI, h::ScanHeaderVD)
    push!(c.headers, h)
    # push!(c.order, h[:adc])
    @show total_length = c[:n_fid] * c[:n_points_on_circle] ÷ c[:n_TI]
    return length(c.headers) * h[:n_adc_points] ≥ total_length
end

function mmaped_image_ice(info, name=tempname())
    n_slices = info[:n_part]
    sz = (info[:n_frequency], info[:n_frequency], n_slices, info[:n_fid], info[:n_channels])
    return mmap(name, Array{ComplexF64,5}, sz)
end



function reconstruct_circle(c::CircleTI)
    info = c.info

    @show info[:n_channels]
    data = open(info[:filename]) do io
        vcat((read_adc(io, head.data_position, info[:n_channels]) for head in c.headers)...)
    end
    data = reshape(data[begin:c[:n_useful_adc_points]], c[:n_points_on_circle], :)
    coords = construct_coordinates_circle(c)
    # return data
    return fourier_transform(data, coords, info[:n_frequency])
end
