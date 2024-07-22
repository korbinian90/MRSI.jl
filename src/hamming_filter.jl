# requires normalized radius r in (-0.5, 0.5)
hamming_filter(r) = 0.54 + 0.46 * cos(2π * r)

function hamming_filter_z(csi)
    csi_k = ifft_slice_dim(csi)

    n_slices = size(csi, 3)
    slice_radii = range(-0.5, 0.5, n_slices)
    filter = hamming_filter.(slice_radii)
    filter = reshape(filter, 1, 1, :) # to apply the filter in the 3rd dimension
    
    csi_k .*= filter

    return fft_slice_dim(csi_k)
end
