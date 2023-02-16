function fov_shift!(slice, kspace_coordinates, info)
    position = info[:position] # (coronal, sagital, transversal)
    # TODO rotation
    x_coordinates = real.(kspace_coordinates)
    y_coordinates = imag.(kspace_coordinates)

    fov_shift1 = calc_shift.(-x_coordinates, position[2], info[:n_phase_encoding], info[:fov_readout])
    fov_shift2 = calc_shift.(y_coordinates, position[1], info[:n_frequency], info[:fov_phase])
    slice .*= (fov_shift1 .* fov_shift2)
end

calc_shift(coord, position, len, fov) = exp(2pi * im * coord * position * len / fov)
