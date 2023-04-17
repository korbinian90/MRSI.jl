function mmaped_image(info, datatype, name=tempname())
    n_slices = info[:n_part]
    sz = (info[:n_frequency], info[:n_frequency], n_slices, info[:n_fid], info[:n_channels])
    return mmap(name, Array{datatype,5}, sz)
    # return zeros(ComplexF64, sz) # to avoid mmapping
end
