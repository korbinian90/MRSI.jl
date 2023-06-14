# A memory mapped object behaves as normal objects stored in RAM, however, 
# it is dynamically stored to disk when the RAM limit is reached.
# This allows working with extremely large objects as if they would fit in memory
function mmaped_image(info, datatype, arg)
    sz = (info[:n_frequency], info[:n_frequency], info[:n_part], info[:n_fid], info[:n_channels])

    if arg == false # no mmap
        return zeros(datatype, sz)
    end

    name = get_name(arg)
    return mmap(name, Array{datatype,5}, sz; shared=false)
end

get_name(::Bool) = tempname() # arg == true
function get_name(arg::AbstractString)
    parts = splitext(arg)
    isfolder = isempty(parts[2])
    if isfolder
        mkpath(arg)
        return tempname(arg)
    end
    mkpath(splitdir(arg)[1])
    return arg
end
