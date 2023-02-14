function extract_twix(twix)
    position = try
        twix["MeasYaps"]["sSpecPara"]["sVoI"]["sPosition"] |> sPos -> (sPos["dCor"], sPos["dSag"], sPos["dTra"])
    catch
        (0, 0, 0) # Parameter doesn't exist
    end
    return Dict(
        :oversampling_factor => twix["Dicom"]["flReadoutOSFactor"],
        :fov_readout => twix["MeasYaps"]["sSpecPara"]["sVoI"]["dReadoutFOV"],
        :fov_phase => twix["MeasYaps"]["sSpecPara"]["sVoI"]["dPhaseFOV"],
        :n_channels => length(map(coil -> coil["lRxChannelConnected"], twix["MeasYaps"]["sCoilSelectMeas"]["aRxCoilSelectData"][1]["asList"])),
        :n_part => twix["MeasYaps"]["sKSpace"]["lPartitions"],
        :n_frequency => twix["MeasYaps"]["sKSpace"]["lBaseResolution"],
        :n_phase_encoding => twix["MeasYaps"]["sKSpace"]["lPhaseEncodingLines"],
        :position => position,
    )
end

function set_info!(info, headers)
    info[:max_n_circles] = info[:n_frequency] ÷ 2
    info[:max_r] = maximum(radius_normalized(h, info) for h in headers)
    info[:max_points_on_circle] = maximum(points_on_circle(h, info) for h in headers)
    info[:n_adc_points] = maximum(h.dims[COL] for h in headers)
    info[:n_part] = length(unique(h.dims[SEG] for h in headers))
    # info[:n_fid] = get_fid(first(first(headers)), info) # calculated on one circle
end

function get_fid(circle_headers, info)
    n_adcs = maximum(h.dims[IDA] for h in circle_headers)
    n_points_on_circle = points_on_circle(circle_headers, info)
    n_temporal_interleaves = maximum(h.dims[IDB] for h in circle_headers)
    return round(Int, info[:n_adc_points] * n_adcs * n_temporal_interleaves / n_points_on_circle - 0.5) # why 0.5?
end

function points_on_circle(header, info)
    n = header.ice_param[6]
    if n == 0
        n = header.dims[IDC] - 1
    end
    return n * info[:oversampling_factor]
end

# Define ScanInfo type that iterates over headers [slice][circle]
struct ScanInfo
    headers::AbstractArray
    info::Dict
end
function ScanInfo(f::AbstractString, type::Symbol)
    if !checkVD(f)
        error("only VD implemented")
    end
    info = extract_twix(read_twix_protocol(f))
    headers = read_scan_headers(f, info[:n_channels])[type]
    set_info!(info, headers)
    headers = rearrange_headers(headers, (@show info[:n_part]), (@show info[:max_n_circles]))
    ScanInfo(headers, info)
end
function Base.getindex(s::ScanInfo, i::Integer)
    if s.headers[i] isa ScanHeaderVD
        return s.headers[i]
    end
    return ScanInfo(s.headers[i], s.info)
end
function Base.getindex(s::ScanInfo, i::Symbol)
    return s.info[i]
end
function Base.iterate(s::ScanInfo, state=1)
    if state > length(s.headers)
        return nothing
    end
    return (s[state], state + 1)
end
Base.length(s::ScanInfo) = length(s.headers)
Base.size(s::ScanInfo, dim...) = size(s.headers, dim...)
Base.lastindex(s::ScanInfo) = length(s)
