const GYRO_MAGNETIC_RATIO_OVER_TWO_PI = 42.57747892; # value taken from MATLAB script read_ascconv.m

function read_slice(io, scaninfo)
    return [read_circle(io, si) for si in scaninfo]
end

function read_circle(io, scaninfo)
    headers = scaninfo.headers
    n_channels = scaninfo.twix[:n_channels]
    points_per_adc = headers[1].dims[COL]

    output = zeros(ComplexF32, points_per_adc, size(headers)..., n_channels)
    for I in CartesianIndices(headers)
        seek(io, 5) # skip the first 5 elements
        output[:, I, :] .= read_adc(io, headers[I].data_position, n_channels) # dims: (adc_points, adc_line, TI, part, channels)
    end
    return output
end

function read_adc(io::IO, (start, adc_length), channels)
    adc = zeros(ComplexF32, adc_length, channels)
    seek(io, start)
    read!(io, adc)
    return view(adc, 5:adc_length, :)
end

function calculate_additional_info(scaninfo)
    s = Dict{Symbol,Int}()
    s[:adc_points] = maximum(h.dims[COL] for h in scaninfo)
    s[:adcs] = maximum(h.dims[IDA] for h in scaninfo)
    s[:part] = length(unique(h.dims[SEG] for h in scaninfo))

    s[:points_on_circle] = maximum(h.ice_param[6] for h in scaninfo)
    if s[:points_on_circle] == 0
        s[:points_on_circle] = maximum(h.dims[IDC] - 1 for h in scaninfo)
    end
    s[:points_on_circle] *= scaninfo.twix[:oversampling_factor]

    s[:temporal_interleaves] = maximum(h.dims[IDB] for h in scaninfo)

    s[:n_fid] = round(Int, s[:adc_points] * s[:adcs] / s[:points_on_circle] * s[:temporal_interleaves] - 0.5)
    s[:useful_adc_points] = (s[:n_fid] * s[:points_on_circle]) ÷ s[:temporal_interleaves]
    s[:n_channels] = scaninfo.twix[:n_channels]
    return s
end
