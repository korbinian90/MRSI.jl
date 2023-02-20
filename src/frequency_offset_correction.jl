function frequency_offset_correction!(data, info; dim=2)
    freq = calculate_frequency(info)
    freq = to_dim(freq, dim)

    timeoffset = 0:info[:n_points_on_circle]-1

    correction = exp.(2pi * im .* freq .* timeoffset)

    kdata = ifftshift(ifft(data, dim), dim)
    kdata .*= correction
    data .= fft(fftshift(kdata, dim), dim)
end

calculate_frequency(info) = ((0:info[:n_fid]-1) ./ info[:n_fid] .+ 0.5 * (1 / info[:n_fid] - 1)) ./ info[:n_points_on_circle] .* info[:n_TI]

to_dim(arr, dim) = reshape(arr, ones(Int, dim - 1)..., :)
