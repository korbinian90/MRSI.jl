@enum DIM COL = 1 LIN = 3 AVE = 4 SLI = 5 PAR = 6 ECO = 7 PHS = 8 REP = 9 SET = 10 SEG = 11 IDA = 12 IDB = 13 IDC = 14 IDD = 15 IDE = 16
Base.to_index(d::DIM) = Int(d)

# Define ScanInfo type that iterates over headers [slice][circle]
struct ScanInfo
    headers::AbstractArray
    twix::Dict
end
function ScanInfo(f::AbstractString, type::Symbol)
    twix = extract_twix(read_twix_protocol(f))
    n_channels = twix[:n_channels]
    n_part = twix[:n_part]
    n_grid = twix[:n_frequency]
    max_n_circles = n_grid ÷ 2
    ScanInfo(read_rearrange_data_headers(f, type, n_channels, n_part, max_n_circles), twix)
end
function Base.getindex(s::ScanInfo, i)
    if s.headers[i] isa ScanHeaderVD
        return s.headers[i]
    end
    return ScanInfo(s.headers[i], s.twix)
end
function Base.iterate(s::ScanInfo, state=1)
    if state > length(s.headers)
        return nothing
    end
    return (s[state], state + 1)
end
Base.length(s::ScanInfo) = length(s.headers)
Base.size(s::ScanInfo, dim...) = size(s.headers, dim...)

struct HeaderInfo
    n_scans
    meas_id
    file_id
    meas_offset
    meas_length
end
function Base.read(io::IO, ::Type{HeaderInfo})
    first = read(io, UInt32)
    n_scans = read(io, UInt32)
    meas_id = read(io, UInt32)
    file_id = read(io, UInt32)
    n_scans = n_scans
    meas_offset = zeros(Int, n_scans)
    meas_length = zeros(Int, n_scans)
    for i in 1:n_scans
        meas_offset[i] = read(io, UInt64)
        meas_length[i] = read(io, UInt64)
        seek(io, position(io) + 152 - 16)
    end
    meas_offset = meas_offset
    meas_length = meas_length
    HeaderInfo(n_scans, meas_id, file_id, meas_offset, meas_length)
end

struct ScanHeaderVD
    mask
    type
    dims
    ice_param
    data_position
end
function Base.read(io::IO, ::Type{ScanHeaderVD})
    header_start = position(io)
    mdh_byte_length = 192
    mask_offset = 40
    channel_mdh_offset_64 = 4
    ice_param_offset = 48 + (48 - 8) * 2

    seek(io, position(io) + mask_offset)
    mask_bytes = zeros(UInt8, 8)
    read!(io, mask_bytes)
    mask = get_scan_mask(mask_bytes)

    dims = zeros(Int16, 16)
    read!(io, dims)
    dims = Int.(dims)
    dims[3:end] .+= 1

    ice_param = zeros(Int16, 24)
    seek(io, header_start + ice_param_offset)
    read!(io, ice_param)

    data_start = header_start + mdh_byte_length
    adc_length = dims[COL] + channel_mdh_offset_64

    ScanHeaderVD(mask, scan_info(mask), dims, ice_param, (data_start, adc_length))
end
