function fov_shift!(slice, kspace_coordinates, scaninfo)
    n_phase = scaninfo.twix[:n_phase_encoding]
    n_read = scaninfo.twix[:n_frequency]
    fov_read = scaninfo.twix[:fov_readout]
    fov_phase = scaninfo.twix[:fov_phase]
    position = scaninfo.twix[:position] # (coronal, sagital, transversal)

    x_coordinates = real.(kspace_coordinates)
    y_coordinates = imag.(kspace_coordinates)

    fov_shift1 = calc_shift.(-x_coordinates, position[2], n_phase, fov_read)
    fov_shift2 = calc_shift.(y_coordinates, position[1], n_read, fov_phase)

    slice .*= (fov_shift1 .* fov_shift2)
end

function calc_shift(coord, position, len, fov)
    exp(2pi * im * coord * position * len / fov)
end
