function hamming_radii(n_loops, info)
    max_gm_distance = 1e3 / info[:fov_readout] / GYRO_MAGNETIC_RATIO_OVER_TWO_PI
    kmax = max_gm_distance / 2 * 2 * floor(info[:n_phase_encoding] / 2 - 1)

    f(k, kmax) = (-0.46 * kmax^2 / pi^2 + 0.27 * k^2 + 0.46 * kmax^2 * cos(k * pi / kmax) / pi^2 + 0.46 * k * kmax * sin(k * pi / kmax) / pi) / (kmax * (0.27 - 2 * 0.46 / pi^2))
    nyquist_radial_offset = f(max_gm_distance / 2, kmax)
    # enlarge kmax because we want to keep it compareable to shifted natural
    kmax = kmax + max_gm_distance / 2

    # Coefficients for the Inverse Series of the Series Expansion of the function f(k, kmax)
    poly_coeffs = [0.5946167018324481, 0.11931056466536513, 0.060659083690310554, 0.03950646580294113, 0.028860886242479126, 0.022514570108001085, 0.018317955376041764, 0.015340571927137355, 0.013117598365113984, 0.011392266426157073, 0.010011731169853292, 0.008879593655298277, 0.007932186533645164, 0.007125847703410522, 0.006429682485159981, 0.005821244858204183, 0.005283854264506816, 0.004804870085466081, 0.004374548421712105, 0.003985264704965461, 0.00363097, 0.00330682, 0.00300887, 0.00273391, 0.00247927, 0.00224273, 0.00202243, 0.0018168, 0.00162451, 0.00144441]
    function radius(i)
        if i == 1
            return max_gm_distance / 2
        end
        # shift now every index such that we dont need to omeasure the k-space center
        scaled = (nyquist_radial_offset + (kmax - nyquist_radial_offset) / n_loops * (i - 1))
        return sum(c * scaled^((2i - 1) / 2) / kmax^((2i - 3) / 2) for (i, c) in enumerate(poly_coeffs))
    end

    return radius.(1:n_loops)
end
