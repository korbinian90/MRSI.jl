function lipid_suppression!(csi::AbstractArray{T,4}, mask; type=:L2, L2_beta=0.1f0, L1_n_loops=5, args...)
    lipid_mask = lipid_mask_from_mrsi_pipeline(csi, info; fat_range_ppm=[1.8, 0.5], threshold=0.1)
    scale = 1.1047e+09/norm(csi(:))
    for slice in axes(csi, 3)
        csi_slice = csi[:,:,slice,:] * scale
        csi[:, :, slice, :] = lipid_suppression_slice(csi_slice, mask[:, :, slice], lipid_mask[:, :, slice], type, L2_beta, L1_n_loops, args...) / scale
    end
end

function lipid_suppression!(csi::AbstractArray{T,5}, mask; args...) where T
    for channel in axes(csi, 5)
        csi[:, :, :, :, channel] = lipid_suppression!(csi[:, :, :, :, channel], mask; args...)
    end
end

function lipid_suppression_slice(csi, mask, lipid_mask, type, L2_beta, L1_n_loops, args...)
    fft_dim = 3 # spectrum
    csi = fftshift(fft(csi, fft_dim), fft_dim)

    lipid_basis = generate_lipid_basis(csi, lipid_mask) # Q: kann man hier auch gleich lipid basis definieren, anstatt die lipid_basis aus der lipid_mask zu generieren? Die lipid_mask wird ja dann auch aus dem Spektrum generiert
    if type == :L2
        csi = lipid_decontamination_L2(csi, lipid_basis, mask, L2_beta) # brain_mask looks unimportant (only for efficiency?)
    elseif type == :L1
        csi = lipid_decontamination_L1(csi, lipid_basis, mask, L1_n_loops; args...)
    end

    csi = ifft(ifftshift(csi, fft_dim), fft_dim)
    return csi
end

function generate_lipid_basis(csi, lipid_mask)
    count = 0
    lipid_basis = zeros(eltype(csi), sum(lipid_mask), size(csi, 3))

    for J in CartesianIndices(lipid_mask)
        if lipid_mask[J]
            count += 1
            lipid_basis[count, :] = csi[J, :]
        end
    end

    return lipid_basis
end
