function mmaped_image(scaninfo, name=tempname())
    n_fid = MRSI.calculate_additional_info(first(first(scaninfo)))[:n_fid]
    n_slices = prod(size(scaninfo))
    sz = (scaninfo[:n_frequency], scaninfo[:n_frequency], n_slices, n_fid, scaninfo[:n_channels])
    return mmap(name, Array{ComplexF64,5}, sz)
end
