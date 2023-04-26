function read_data(c::Circle, datatype)
    data = zeros(datatype, c[:n_points_on_circle], c[:n_fid], c[:n_channels])
    open(c[:filename]) do io
        for (TI, heads) in enumerate(c.headers)
            data[:, TI:c[:n_TI]:end, :] .= read_reshape_one_TI(io, heads, c)
        end
    end
    return data
end

function read_reshape_one_TI(io, heads, info)
    data = read_one_TI(io, heads, info[:n_channels])
    data = @view data[1:info[:n_useful_adc_points], :]
    return reshape(data, info[:n_points_on_circle], :, info[:n_channels])
end

function read_one_TI(io, headers, n_channels)
    vcat((read_adc(io, head.data_position, n_channels) for head in headers)...)
end

function read_adc(io::IO, (start, adc_length), channels)
    adc = zeros(ComplexF32, adc_length, channels)
    seek(io, start)
    read!(io, adc)
    return view(adc, 5:adc_length, :)
end
