function noise_decorrelation(noiseinfo)

    noise = read_noise(noiseinfo)

    noise_correlation = calculate_noise_correlation(noise)

    # noise_scaling_factor = noise_dwelltime / (reco_dwelltime * nTempIntsPerAngInt / datasize{1}(1))
    noise_scaling_factor = 1.25

    noise_matrix_cholesky = cholesky(noise_correlation * noise_scaling_factor / 2).U
    
    return noise_matrix_cholesky
end

function calculate_noise_correlation(data, cut_first_n_points=19)
    cut_first_n_points = min(max(cut_first_n_points, 1), size(data, 1)) # between 1 and size(data, 2)
    data = @view data[cut_first_n_points:end, :]

    correlation = 1 / size(data, 1) * (data' * data)

    return correlation
end
