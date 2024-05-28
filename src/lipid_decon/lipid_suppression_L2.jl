function lipid_decontamination_L2(csi, lipid_basis, mask, beta)
    beta *= 1.0f-15
    lipid_basis = Matrix(transpose(lipid_basis))
    lipid_inv = inv(I + beta * (lipid_basis * lipid_basis'))

    Threads.@threads for J in CartesianIndices(mask)
        if mask[J]
            csi[J, :] = lipid_inv * csi[J, :]
        end
    end

    csi .*= mask
    return csi
end
