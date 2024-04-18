function fov_shift!(circle, kspace_coordinates, info)
    R1 = rotation_between(info[:slice_normal], [0, 0, 1])
    R2 = RotZ(-info[:in_plane_rotation])
    R = R2 * R1

    rotated_position = R * info[:position] # (coronal, sagital, transversal)

    x_coordinates = real.(kspace_coordinates)
    y_coordinates = imag.(kspace_coordinates)

    fov_shift1 = calc_shift.(-x_coordinates, rotated_position[2], info[:n_phase_encoding], info[:fov_readout])
    fov_shift2 = calc_shift.(y_coordinates, rotated_position[1], info[:n_frequency], info[:fov_phase])
    circle .*= (fov_shift1 .* fov_shift2)
end

calc_shift(coord, position, len, fov) = exp(2pi * im * coord * position * len / fov)
