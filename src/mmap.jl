function mmaped_image(info, name=tempname())
    n_slices = info[:n_part]
    sz = (info[:n_frequency], info[:n_frequency], n_slices, info[:n_fid], info[:n_channels])
    return mmap(name, Array{ComplexF64,5}, sz)
end
