const GYRO_MAGNETIC_RATIO_OVER_TWO_PI = 42.57747892; # value taken from MATLAB script read_ascconv.m

function read_data_headers(filename, n_channels)
    if !checkVD(filename)
        error("only VD implemented")
    end
    headers = read_scan_headers(filename, n_channels)
    return rearrange_headers(headers)
end

function read_scan_headers(filename, n_channels; scan=1)  # only with one scan implemented for now
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

            start, n_adc = scan.data_position
            if haskey(scan_headers, scan.type)
                push!(scan_headers[scan.type], scan)
            else
                scan_headers[scan.type] = [scan]
            end
            n_bytes = n_adc * n_channels * 8
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

                scan_headers_reshaped = Array{ScanHeaderVD}(undef, adcs, temporal_interleaves, part_encodings)
                for head in circle_headers
                    TI = head.dims[IDB]
                    part = head.dims[SEG]
                    adc = head.dims[IDA]
                    scan_headers_reshaped[adc, TI, part] = head
                end

                headers_rearranged[type][slice][circle] = scan_headers_reshaped
            end
        end
    end
    return headers_rearranged
end

function calculate_additional_data(scaninfo)
    s = Dict{Symbol,Int}()
    s[:adc_points] = maximum(h.dims[COL] for h in scaninfo)
    s[:adcs] = maximum(h.dims[IDA] for h in scaninfo)
    s[:part] = length(unique(h.dims[SEG] for h in scaninfo))

    s[:points_on_circle] = maximum(h.ice_param[6] for h in scaninfo)
    if s[:points_on_circle] == 0
        s[:points_on_circle] = maximum(h.dims[IDC] - 1 for h in scaninfo)
    end
    s[:points_on_circle] *= read_oversampling_factor(scaninfo)

    s[:temporal_interleaves] = maximum(h.dims[IDB] for h in scaninfo)

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

function read_adc(io::IO, (start, adc_length), channels)
    adc = zeros(ComplexF32, adc_length, channels)
    seek(io, start)
    read!(io, adc)
    return view(adc, 5:adc_length, :)
end

function read_slice(io, scaninfo)
    return [read_circle(io, si) for si in scaninfo]
end

function read_circle(io, scaninfo)
    headers = scaninfo.scan_headers
    n_channels = read_n_channels(scaninfo)
    points_per_adc = headers[1].dims[COL]

    output = zeros(ComplexF32, points_per_adc, size(headers)..., n_channels)
    for I in CartesianIndices(headers)
        seek(io, 5) # skip the first 5 elements
        output[:, I, :] .= read_adc(io, headers[I].data_position, n_channels) # dims: (adc_points, adc_line, TI, part, channels)
    end
    return output
end

function kspace_coordinates(sliceinfo)
    n_grid = read_n_grid(sliceinfo)
    fov_read = read_fov_readout(sliceinfo)
    coords = read_kspace_coordinates(sliceinfo)
    return normalize_kspace(coords, n_grid, fov_read)
end

function read_kspace_coordinates(slice_headers)
    coordinates = Complex[]
    for circle_headers in slice_headers
        s = calculate_additional_data(circle_headers)
        xy = read_kspace_coordinate(circle_headers[1])
        coords = construct_coordinates_circle(xy, s[:points_on_circle])
        append!(coordinates, coords)
    end
    return coordinates
end

function construct_coordinates_circle(xy, points_on_circle)
    angle_increment = 2pi / points_on_circle
    r = abs(xy)
    angle_first_point = angle(xy) - pi / 2
    phi = angle_first_point .- (0:points_on_circle-1) * angle_increment
    coords = r .* exp.(im .* phi)
    return coords
end

function normalize_kspace(coordinates, n_grid, fov_read)
    return coordinates ./ 2max_r(n_grid, fov_read)
end

function read_kspace_coordinate(head::ScanHeaderVD)
    x = parse_to_float(head.ice_param[1], head.ice_param[2])
    y = parse_to_float(head.ice_param[3], head.ice_param[4])
    return Complex(x, y)
end
parse_to_float(a::Int16, b::Int16) = a + b / 10_000

function radius_normalized(xy, n_grid, fov_read)
    r = abs(xy)
    delta_gm = 1e6 / (fov_read * GYRO_MAGNETIC_RATIO_OVER_TWO_PI)
    maxR = delta_gm * n_grid / 2
    return r / (2maxR)
end

function max_r(n_grid, fov_read)
    delta_gm = 1e6 / (fov_read * GYRO_MAGNETIC_RATIO_OVER_TWO_PI)
    max_r = delta_gm * n_grid / 2
    return max_r
end
