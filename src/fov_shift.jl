function fov_shift!(slice, kspace_coordinates, scaninfo)
    n_phase = scaninfo[:n_phase_encoding]
    n_read = scaninfo[:n_frequency]
    fov_read = scaninfo[:fov_readout]
    fov_phase = scaninfo[:fov_phase]
    position = scaninfo[:position] # (coronal, sagital, transversal)
    # TODO rotation
    x_coordinates = real.(kspace_coordinates)
    y_coordinates = imag.(kspace_coordinates)

    fov_shift1 = calc_shift.(-x_coordinates, position[2], n_phase, fov_read)
    fov_shift2 = calc_shift.(y_coordinates, position[1], n_read, fov_phase)
    slice .*= (fov_shift1 .* fov_shift2)
end

function calc_shift(coord, position, len, fov)
    exp(2pi * im * coord * position * len / fov)
end
