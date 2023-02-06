function read_rearrange_data_headers(filename, type, n_channels, n_part, max_n_circles)
    if !checkVD(filename)
        error("only VD implemented")
    end
    headers = read_scan_headers(filename, n_channels)[type]
    return rearrange_headers(headers, n_part, max_n_circles)
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

function checkVD(filename)
    startints = zeros(UInt32, 2)
    read!(filename, startints)
    return startints[1] < 10_000 && startints[2] <= 64
end

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
