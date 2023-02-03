module MRSI

include("headers.jl")
include("read_siemens_data.jl")
include("rearrange.jl")

export read_data_headers,
    read_slice,
    rearrange

end
