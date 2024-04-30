function noise_decorrelation(scan_info)
    noise = read_noise(scan_info[:NOISADJSCAN])

    noise_correlation = calculate_noise_correlation(noise)

    noise_matrix_cholesky = cholesky(noise_correlation * noise_scaling_factor(scan_info) / 2).U
    
    return noise_matrix_cholesky
end

function calculate_noise_correlation(data, start_at_point=19)
    start_at_point = min(max(start_at_point, 1), size(data, 1)) # between 1 and size(data, 1)
    data = @view data[start_at_point:end, :]

    correlation = 1 / size(data, 1) * (data' * data)

    return correlation
end

function noise_scaling_factor(scan_info)
    header = header_first_circle(scan_info)

    noise_dwelltime = 10000 # currently hard-coded since it is not written in sequence
    n_TI = max(header[:n_TI], 1) # workaround: n_TI is not set in header for old sequences
    
    online_dwelltime = scan_info[:ONLINE][:dwelltime]
    points = get_n_points_on_circle(header, scan_info[:ONLINE][:oversampling_factor])    

    return  noise_dwelltime / (online_dwelltime * n_TI / points)
end

function header_first_circle(scan_info)
    for head in scan_info[:ONLINE][:headers]
        if head[:circle] == 1
            return head
        end
    end
end
