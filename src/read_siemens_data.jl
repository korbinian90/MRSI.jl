function read_slice(io, slice_headers, info)
    return [read_circle(io, s, info) for s in slice_headers]
end

function read_circle(io, circle_headers, info; type=ComplexF64)
    info = calculate_circle_info(info, circle_headers)
    output = zeros(type, info[:n_adc_points], size(circle_headers)..., info[:n_channels])
    for I in CartesianIndices(circle_headers) # loops over [adc_line, TI, part]
        seek(io, 5) # skip the first 5 elements
        output[:, I, :] .= read_adc(io, circle_headers[I].data_position, info[:n_channels]) # dims: [adc_points, adc_line, TI, part, channels]
    end
    return output
end

function read_adc(io::IO, (start, adc_length), channels)
    adc = zeros(ComplexF32, adc_length, channels)
    seek(io, start)
    read!(io, adc)
    return view(adc, 5:adc_length, :)
end
