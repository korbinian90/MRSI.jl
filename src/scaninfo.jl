const GYRO_MAGNETIC_RATIO_OVER_TWO_PI = 42.57747892; # value taken from MATLAB script read_ascconv.m

function extract_twix(twix)
    position = try
        twix["MeasYaps"]["sSpecPara"]["sVoI"]["sPosition"] |> sPos -> (sPos["dCor"], sPos["dSag"], sPos["dTra"])
    catch
        (0, 0, 0) # Parameter doesn't exist
    end
    info = Dict(
        :oversampling_factor => twix["Dicom"]["flReadoutOSFactor"],
        :fov_readout => twix["MeasYaps"]["sSpecPara"]["sVoI"]["dReadoutFOV"],
        :fov_phase => twix["MeasYaps"]["sSpecPara"]["sVoI"]["dPhaseFOV"],
        :n_channels => length(map(coil -> coil["lRxChannelConnected"], twix["MeasYaps"]["sCoilSelectMeas"]["aRxCoilSelectData"][1]["asList"])),
        :n_part => twix["MeasYaps"]["sKSpace"]["lPartitions"],
        :n_frequency => twix["MeasYaps"]["sKSpace"]["lBaseResolution"],
        :n_phase_encoding => twix["MeasYaps"]["sKSpace"]["lPhaseEncodingLines"],
        :position => position,
    )
    calculate_additional_info!(info)
    return info
end

function calculate_additional_info!(info)
    info[:max_n_circles] = info[:n_frequency] ÷ 2
    info[:max_r] = max_r(info[:n_frequency], info[:fov_readout])
    info[:part_order], info[:circle_order], info[:circles_per_part] = calculate_part_order(info[:n_part], info[:max_n_circles])
    # info[:max_points_on_circle] = maximum(points_on_circle(h, info) for h in headers)
    # info[:n_adc_points] = maximum(h.dims[COL] for h in headers)
    # info[:n_seg] = length(unique(h.dims[SEG] for h in headers))
    # info[:n_fid] = get_fid(first(first(headers)), info) # calculated on one circle
end

function max_r(n_grid, fov_read)
    delta_gm = 1e6 / (fov_read * GYRO_MAGNETIC_RATIO_OVER_TWO_PI)
    max_r = delta_gm * n_grid / 2
    return max_r
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
    headers = rearrange_headers(headers, info[:n_part], info[:max_n_circles])
    ScanInfo(headers, info)
end
function Base.getindex(s::ScanInfo, i::Integer)
    if s.headers[i] isa ScanHeaderVD
        return s.headers[i]
    end
    return ScanInfo(s.headers[i], s.info)
end
function Base.getindex(s::ScanInfo, i::Symbol)
    if i == :n_adc_points
        first(s.headers)[:n_adc_points]
    elseif i == :n_points_on_circle
        get_n_points_on_circle(first(s.headers), s.info[:oversampling_factor])
    else
        s.info[i]
    end
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

# Headers for one circle (adcs accumulated)
mutable struct Circle
    headers # [TI][adc]
    info
end
Circle(info) = Circle(nothing, info)
function Base.getindex(c::Circle, s::Symbol)
    if haskey(c.info, s) # look in info
        c.info[s]
    elseif s == :part
        c[:part_order][c[:LIN]]
    elseif s == :circle
        c[:circle_order][c[:LIN]]
    elseif s == :n_adcs # calculate here
        maximum(h.dims[IDA] for h in c.headers)
    elseif s == :n_points_on_circle
        get_n_points_on_circle(first(first(c.headers)), c.info[:oversampling_factor])
    elseif s == :n_fid
        get_n_fid(c[:n_adc_points], c[:n_adcs], c[:n_TI], c[:n_points_on_circle] - 0.5) # why 0.5?
    elseif s == :fid_points_for_TI
        c[:TI]:c[:n_TI]:c[:n_fid]
    elseif s == :n_useful_adc_points
        (c[:n_fid] * c[:n_points_on_circle]) ÷ c[:n_TI]
    else # look in ScanHeaderVD
        first(first(c.headers))[s]
    end
end

get_n_fid(n_adc_points, n_adcs, n_TI, n_points_on_circle) = round(Int, n_adc_points * n_adcs * n_TI / n_points_on_circle - 0.5) # why 0.5?
get_n_points_on_circle(h::ScanHeaderVD, oversampling_factor) = round(Int, max(h.dims[IDC] - 1, h.ice_param[6]) * oversampling_factor)
