# A memory mapped object behaves as normal objects stored in RAM, however, 
# it is dynamically stored to disk when the RAM limit is reached.
# This allows working with extremely large objects as if they would fit in memory
function mmaped_image(info, datatype, name=tempname())
    n_slices = info[:n_part]
    sz = (info[:n_frequency], info[:n_frequency], n_slices, info[:n_fid], info[:n_channels])
    return mmap(name, Array{datatype,5}, sz)
    # return zeros(datatype, sz) # to avoid mmapping
end
