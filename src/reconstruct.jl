read_oversampling_factor(s::ScanInfo) = s.twix_header["Dicom"]["flReadoutOSFactor"]
read_fov_readout(s::ScanInfo) = s.twix_header["MeasYaps"]["sSpecPara"]["sVoI"]["dReadoutFOV"]
read_n_channels(s::ScanInfo) = s.twix_header["MeasYaps"]["sCoilSelectMeas"]["aRxCoilSelectData"][1]["asList"][1]["lRxChannelConnected"]
function read_n_grid(s::ScanInfo)
    n_frequency = s.twix_header["MeasYaps"]["sKSpace"]["lBaseResolution"]
    n_phase_encoding = s.twix_header["MeasYaps"]["sKSpace"]["lPhaseEncodingLines"]
    @assert n_frequency == n_phase_encoding
    return n_frequency
end

function reconstruct(filename)
    scaninfo = ScanInfo(filename, :ONLINE)
    image = cat((reconstruct_slice(filename, si) for si in scaninfo)...; dims=3)
    return image
end

function reconstruct_slice(filename, scaninfo)
    n_grid = read_n_grid(scaninfo)

    slice = open(filename) do io
        read_slice(io, scaninfo)
    end

    kspace_data = rearrange(slice, scaninfo)
    kspace_points = MRSI.kspace_coordinates(scaninfo)
    return fourier_transform(kspace_data, kspace_points, n_grid)
end
