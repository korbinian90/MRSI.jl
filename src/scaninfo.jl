const GYRO_MAGNETIC_RATIO_OVER_TWO_PI = 42.57747892; # value taken from MATLAB script read_ascconv.m

function extract_twix(twix)
    info = Dict(
        :oversampling_factor => twix["Dicom"]["flReadoutOSFactor"],
        :fov_readout => twix["MeasYaps"]["sSpecPara"]["sVoI"]["dReadoutFOV"],
        :fov_phase => twix["MeasYaps"]["sSpecPara"]["sVoI"]["dPhaseFOV"],
        :n_channels => length(twix["MeasYaps"]["sCoilSelectMeas"]["aRxCoilSelectData"][1]["asList"]),
        :n_part => twix["MeasYaps"]["sKSpace"]["lPartitions"],
        :n_frequency => twix["MeasYaps"]["sKSpace"]["lBaseResolution"],
        :n_phase_encoding => twix["MeasYaps"]["sKSpace"]["lPhaseEncodingLines"],
        :position => pos_as_vector(twix["MeasYaps"]["sSliceArray"]["asSlice"][1], "sPosition"),
        :slice_normal => pos_as_vector(twix["MeasYaps"]["sSliceArray"]["asSlice"][1], "sNormal"),
    )
    calculate_twix_info!(info)
    return info
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
    info[:part_order], info[:circle_order], info[:circles_per_part] = calculate_part_order(info[:n_part], info[:max_n_circles])
end

# only for oldADC
function calculate_part_order(n_part, max_n_circles)
    part_max = n_part ÷ 2
    circles_per_part = [calculate_n_circles_per_part(part_max, max_n_circles, i) for i in -part_max:part_max]
    part_order = vcat((repeat([i], n) for (i, n) in enumerate(circles_per_part))...)
    circle_order = vcat((1:c for c in circles_per_part)...)
    return (part_order, circle_order, circles_per_part)
end

function calculate_n_circles_per_part(part_max, max_n_circles, i)
    if part_max == 0 # 2D
        return max_n_circles
    end
    return max(2, Int(ceil(sqrt(max_n_circles^2 - (i * max_n_circles / part_max)^2))))
end


function calculate_fid(info, headers)
    max_adc_and_ti_header = argmax(h -> h[:adc] * h[:TI], headers)

    n_adcs = max_adc_and_ti_header[:adc]
    n_adc_points = max_adc_and_ti_header[:n_adc_points]
    n_points_on_circle = get_n_points_on_circle(max_adc_and_ti_header, info[:oversampling_factor])
    n_TI = max_adc_and_ti_header[:TI]

    return get_n_fid(n_adc_points, n_adcs, n_TI, n_points_on_circle)
end
get_n_fid(n_adc_points, n_adcs, n_TI, n_points_on_circle) = round(Int, n_adc_points * n_adcs * n_TI / n_points_on_circle - 0.5) # why 0.5?

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
    if haskey(c.info, s) # look in info
        c.info[s]
    elseif s == :n_TI
        length(c.headers)
    elseif s == :part
        c[:part_order][c[:LIN]] # only oldADC
    elseif s == :circle
        c[:circle_order][c[:LIN]] # only oldADc
    elseif s == :n_points_on_circle
        get_n_points_on_circle(first(first(c.headers)), c.info[:oversampling_factor])
    elseif s == :fid_points_for_TI
        c[:TI]:c[:n_TI]:c[:n_fid]
    elseif s == :n_useful_adc_points
        (c[:n_fid] * c[:n_points_on_circle]) ÷ c[:n_TI]
    else # look in ScanHeaderVD
        first(first(c.headers))[s]
    end
end

function Base.getindex(h::ScanHeaderVD, s::Symbol)
    if s == :n_adc_points
        h.dims[COL]
    elseif s == :n_TI # only for oldADC
        h.dims[IDD] - 1
    elseif s == :TI
        h.dims[IDB]
    elseif s == :adc
        h.dims[IDA]
    elseif s == :LIN
        h.dims[LIN]
    end
end
get_n_points_on_circle(h::ScanHeaderVD, oversampling_factor) = round(Int, max(h.dims[IDC] - 1, h.ice_param[6]) * oversampling_factor)
