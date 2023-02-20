function read_scan_headers(info)
    scan_headers = Dict()
    open(info[:filename]) do io
        seek_to_first_scan_header!(io)
        while true
            scan = read(io, ScanHeaderVD)
            if scan.mask[1] # signals ACQEND
                break
            end
            store_header!(scan_headers, scan)
            seek_to_next_header!(io, scan.data_position, info[:n_channels])
        end
    end
    return scan_headers
end

function seek_to_first_scan_header!(io; scan=1)
    file_header = read(io, HeaderInfo)
    offset = file_header.meas_offset[scan]
    seek(io, offset)
    # skip twix header
    header_length = read(io, UInt32)
    seek(io, offset + header_length)
end

function store_header!(headers, head)
    if haskey(headers, head.type)
        push!(headers[head.type], head)
    else
        headers[head.type] = [head]
    end
end

function seek_to_next_header!(io, (start, n_adc), n_channels)
    n_bytes = n_adc * n_channels * 8
    seek(io, start + n_bytes)
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
