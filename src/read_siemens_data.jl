const GYRO_MAGNETIC_RATIO_OVER_TWO_PI = 42.57747892; # value taken from MATLAB script read_ascconv.m

function read_data_headers(filename)
    if !checkVD(filename)
        error("only VD implemented")
    end
    headers = read_scan_headers(filename)
    return rearrange_headers(headers)
end

function read_scan_headers(filename; scan=1)  # only with one scan implemented for now
    scan_headers = Dict()
    open(filename) do io
        file_header = read(io, HeaderInfo)
        offset = file_header.meas_offset[scan]
        seek(io, offset)

        # skip twix header
        header_length = read(io, UInt32)
        seek(io, offset + header_length)

        while true
            scan = read(io, ScanHeaderVD)
            if scan.mask[1] # signals ACQEND
                break
            end

            start, n_bytes = scan.data_position
            if haskey(scan_headers, scan.type)
                push!(scan_headers[scan.type], scan)
            else
                scan_headers[scan.type] = [scan]
            end
            seek(io, start + n_bytes)
        end
    end
    return scan_headers
end

function rearrange_headers(scan_headers)
    headers_rearranged = Dict()
    for (type, headers) in scan_headers
        slices = unique(h.dims[SLI] for h in headers)
        headers_rearranged[type] = Array{Any}(undef, maximum(slices))
        for slice in slices
            slice_headers = [h for h in headers if h.dims[SLI] == slice]
            circles = unique(h.dims[LIN] for h in slice_headers)

            headers_rearranged[type][slice] = Array{Any}(undef, maximum(circles))
            for circle in circles
                circle_headers = [h for h in headers if h.dims[LIN] == circle]

                temporal_interleaves = maximum(h.dims[IDB] for h in circle_headers)
                part_encodings = length(unique(h.dims[SEG] for h in circle_headers))
                adcs = maximum(h.dims[IDA] for h in circle_headers)

                scan_headers_reshaped = Array{ScanHeaderVD}(undef, temporal_interleaves, part_encodings, adcs)
                for head in circle_headers
                    TI = head.dims[IDB]
                    part = head.dims[SEG]
                    adc = head.dims[IDA]
                    scan_headers_reshaped[TI, part, adc] = head
                end

                headers_rearranged[type][slice][circle] = scan_headers_reshaped
            end
        end
    end
    return headers_rearranged
end

function calculate_additional_data(circle_headers)
    s = Dict{Symbol,Int}()
    s[:adc_points] = maximum(h.dims[COL] for h in circle_headers)
    s[:adcs] = maximum(h.dims[IDA] for h in circle_headers)
    s[:part] = length(unique(h.dims[SEG] for h in circle_headers))

    s[:points_on_circle] = maximum(h.ice_param[6] for h in circle_headers)
    if s[:points_on_circle] == 0
        s[:points_on_circle] = maximum(h.dims[IDC] - 1 for h in circle_headers)
    end
    s[:temporal_interleaves] = maximum(h.dims[IDB] for h in circle_headers)

    s[:n_fid] = round(Int, s[:adc_points] * s[:adcs] / s[:points_on_circle] * s[:temporal_interleaves] - 0.5)
    s[:useful_adc_points] = (s[:n_fid] * s[:points_on_circle]) ÷ s[:temporal_interleaves]
    return s
end

function checkVD(filename)
    startints = zeros(UInt32, 2)
    read!(filename, startints)
    return startints[1] < 10_000 && startints[2] <= 64
end

get_col(scan) = 1 * scan.mdh_info[1]

function get_scan_mask(bytes)
    str = join(reverse.(bitstring.(bytes)))
    return [c == '1' for c in str]
end

function scan_info(mask)
    if mask[26]
        return :NOISADJSCAN
    elseif mask[24]
        return :PATREFANDIMASCAN
    elseif mask[23]
        return :PATREFSCAN
    elseif mask[4]
        return :ONLINE
    elseif mask[1]
        return :ACQEND
    end
    return :OTHER
end

function read_adc(io::IO, (start, bytes), channels)
    data_length = bytes ÷ 8
    adc = zeros(ComplexF32, data_length, channels)
    seek(io, start)
    read!(io, adc)
    return view(adc, 5:data_length, :)
end

function read_slice(io, headers, n_channels)
    return [read_circle(io, h, n_channels) for h in headers]
end

function read_circle(io, headers, n_channels)
    points_per_adc = headers[1].dims[COL]
    output = zeros(ComplexF32, size(headers)..., points_per_adc, n_channels)
    for I in CartesianIndices(headers)
        seek(io, 5)
        output[I, :, :] .= read_adc(io, headers[I].data_position, n_channels)
    end
    return output
end

function kspace_coordinates(slice_headers, n_grid, fov_read)
    coordinates = Complex[]
    for (i_circle, headers) in enumerate(slice_headers)
        s = calculate_additional_data(headers)
        points_on_circle = s[:points_on_circle]
        angle_increment = 2pi / points_on_circle
        xy = read_kspace_coordinate(headers[1])
        r = radius_normalized(xy, n_grid, fov_read)
        angle_first_point = angle(xy) - pi / 2
        for i_phi in 0:(points_on_circle-1)
            phi = angle_first_point - i_phi * angle_increment
            push!(coordinates, r * exp(im * phi))
        end

    end
    return coordinates
end

function read_kspace_coordinate(head::ScanHeaderVD)
    x = parse_to_float(head.ice_param[1], head.ice_param[2])
    y = parse_to_float(head.ice_param[3], head.ice_param[4])
    return ComplexF32(x, y)
end
parse_to_float(a::Int16, b::Int16) = a + b / 10_000

function radius_normalized(xy, n_grid, fov_read)
    r = abs(xy)
    delta_gm = 1e6 / (fov_read * GYRO_MAGNETIC_RATIO_OVER_TWO_PI)
    maxR = delta_gm * n_grid / 2
    return r / (2maxR)
end

#





## Adapted from Philipp Ehses read_twix_hdr.m
function read_twix_protocol(filename, start)
    protocol = Dict()
    open(filename) do io
        seek(io, start)
        n_entries = read(io, UInt32)
        for i in 1:n_entries
            chars = [read(io, Char) for _ in 1:10]
            name = parse_entry_name(chars)
            seek(io, position(io) + length(name) - 9)
            len = read(io, UInt32)
            buffer = [read(io, Char) for _ in 1:len]
            # TODO
        end
    end
end

function parse_entry_name(chars)
    return match(r"^\w*", join(chars)).match
end

# n = read_twix_protocol(f, p)
