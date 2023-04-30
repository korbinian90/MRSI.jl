const GYRO_MAGNETIC_RATIO_OVER_TWO_PI = 42.57747892; # value taken from MATLAB script read_ascconv.m

function read_scan_info(filename, type)
    checkVD(filename) || error("only implemented for VD files")
    info = extract_twix(filename, type)
    info[:filename] = filename
    headers = read_scan_headers(info)[type]

    ## Info cheating (obtained from all headers instead of streamed like in ICE)
    info[:radii] =  read_radii(headers, info) # problem, required for dcf
    info[:max_n_points_on_circle] = maximum([get_n_points_on_circle(h, info[:oversampling_factor]) for h in headers]) # can be corrected afterwards
    ##
    return headers, info
end

function create_dict_from_twix(twix, type)
    return Dict(
        :oversampling_factor => twix["Dicom"]["flReadoutOSFactor"],
        :fov_readout => twix["MeasYaps"]["sSpecPara"]["sVoI"]["dReadoutFOV"],
        :fov_phase => twix["MeasYaps"]["sSpecPara"]["sVoI"]["dPhaseFOV"],
        :n_channels => length(twix["MeasYaps"]["sCoilSelectMeas"]["aRxCoilSelectData"][1]["asList"]),
        :n_part => twix["MeasYaps"]["sKSpace"]["lPartitions"],
        :n_frequency => twix["MeasYaps"]["sKSpace"]["lBaseResolution"],
        :n_phase_encoding => twix["MeasYaps"]["sKSpace"]["lPhaseEncodingLines"],
        :position => pos_as_vector(twix["MeasYaps"]["sSliceArray"]["asSlice"][1], "sPosition"),
        :slice_normal => pos_as_vector(twix["MeasYaps"]["sSliceArray"]["asSlice"][1], "sNormal"),
        :larmor_frequency => twix["MeasYaps"]["sTXSPEC"]["asNucleusInfo"][1]["lFrequency"],
        :dwelltime => twix["MeasYaps"]["sRXSPEC"]["alDwellTime"][1],
        :n_fid => get_fid(twix, type),
    )
end

function get_fid(twix, type)
    n_fid = if type == :PATREFSCAN
        twix["Meas"]["alICEProgramPara"][8]
    else
        twix["Meas"]["alICEProgramPara"][7]
    end
    if n_fid == 0 # fix for certain datasets
        n_fid = twix["MeasYaps"]["sSpecPara"]["lVectorSize"]
    end
    return n_fid
end

function pos_as_vector(entry, type)
    if haskey(entry, type)
        [get(entry[type], "dCor", 0), get(entry[type], "dSag", 0), get(entry[type], "dTra", 0)]
    else
        [0, 0, 0]
    end
end

function calculate_twix_info!(info)
    info[:max_n_circles] = info[:n_frequency] ÷ 2
    info[:max_r] = max_r(info[:n_frequency], info[:fov_readout])
    info[:circles_per_part] = calculate_n_circles_per_part(info[:n_part], info[:max_n_circles])
end

read_radii(headers, info) = [radius_normalized(headers[findfirst(h -> h[:circle] == i, headers)], info) for i in 1:info[:max_n_circles]]
calculate_radii(headers, info) = collect(1:2:info[:n_frequency]) * radius_normalized(headers[findfirst(h -> h[:circle] == 1, headers)], info)

function extract_twix(filename, type)
    twix = read_twix_protocol(filename)
    info = create_dict_from_twix(twix, type)
    calculate_twix_info!(info)
    return info
end

function calculate_n_circles_per_part(n_part, max_n_circles)
    part_max = n_part ÷ 2
    if part_max == 0 # 2D
        return max_n_circles
    end
    n_circles(part) = max(2, Int(ceil(sqrt(max_n_circles^2 - (part * max_n_circles / part_max)^2))))
    return [n_circles(part) for part in -part_max:part_max]
end

function max_r(n_grid, fov_read)
    delta_gm = 1e6 / (fov_read * GYRO_MAGNETIC_RATIO_OVER_TWO_PI)
    max_r = delta_gm * n_grid / 2
    return max_r
end

# Headers for one circle (adcs accumulated)
mutable struct Circle
    headers # [TI][adc]
    info
end
Circle(info) = Circle(nothing, info)
function Base.getindex(c::Circle, s::Symbol)
    h = first(first(c.headers))
    if haskey(c.info, s) # look in info
        c.info[s]
    elseif s == :n_TI
        length(c.headers)
    elseif s == :part
         part_from_one(h, c)
    elseif s == :n_points_on_circle
        get_n_points_on_circle(h, c.info[:oversampling_factor])
    elseif s == :fid_points_for_TI
        c[:TI]:c[:n_TI]:c[:n_fid]
    elseif s == :n_useful_adc_points
        (c[:n_fid] * c[:n_points_on_circle]) ÷ c[:n_TI]
    else # look in ScanHeaderVD
        h[s]
    end
end

function Base.getindex(h::ScanHeaderVD, s::Symbol)
    if s == :n_adc_points
        h.dims[COL]
    elseif s == :n_TI
        h.dims[IDD] - 1
    elseif s == :TI
        h.dims[IDB]
    elseif s == :adc
        h.dims[IDA]
    elseif s == :circle
        h.dims[LIN]
    elseif s == :SEG
        h.dims[SEG] - 1 # -n/2 : n/2
    end
end
get_n_points_on_circle(h::ScanHeaderVD, oversampling_factor) = round(Int, max(h.dims[IDC] - 1, h.ice_param[6]) * oversampling_factor)
part_from_one(h::ScanHeaderVD, info) = h[:SEG] + info[:n_part] ÷ 2 + 1
