function read_data(c::Circle, datatype, noise_matrix_cholesky=nothing)
    data = zeros(datatype, c[:n_points_on_circle], c[:n_fid], c[:n_channels])
    open(c[:filename]) do io
        for (TI, heads) in enumerate(c.headers) # c.headers [TI][adc][average]
            data[:, TI:c[:n_TI]:end, :] .= read_reshape_one_TI(io, heads, c, noise_matrix_cholesky)
        end
    end
    return data
end

function read_reshape_one_TI(io, heads, info, noise_matrix_cholesky=nothing)
    data = read_one_TI(io, heads, info[:n_channels])
    data = @view data[1:info[:n_useful_adc_points], :]
    if !isnothing(noise_matrix_cholesky)
        data /= noise_matrix_cholesky
    end
    data = reshape(data, info[:n_points_on_circle], :, info[:n_channels])
    return data
end

# headers [adc][average]
read_one_TI(io, headers, n_channels) = vcat(read_adc_average.(io, headers, n_channels)...)

# headers [average]
read_adc_average(io, headers, n_channels) = sum(read_adc(io, head.data_position, n_channels) for head in headers) ./ length(headers)

function read_adc(io::IO, (start, adc_length), channels)
    adc = zeros(ComplexF32, adc_length, channels)
    seek(io, start)
    read!(io, adc)
    return @view adc[5:adc_length, :]
end

function read_noise(info)
    open(info[:filename], "r") do io
        vcat((MRSI.read_adc(io, head.data_position, info[:n_channels]) for head in info[:headers])...)
    end    
end
