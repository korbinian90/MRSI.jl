"""
    image = reconstruct(filename; combine=:auto, datatype=ComplexF32, do_fov_shift=true, do_freq_cor=true, do_dens_comp=true, conj_in_beginning=true, ice=false, old_headers=false, mmap=true)

Reconstructs a SIEMENS dat file

Coil combine is performed only for AC datasets. 
Set combine to `false`|`true` to force coil combination (including normalization by the patrefscan). 
Set ice to `true` to use calculated radii for density compensation instead of reading them from the headers. 
Set `old_headers` to `true` for older dat files with part and circle stored in LIN
Set `mmap` to `false` to compute in RAM or a filename or folder for storing the temporary result array.

    reconstruct(filename, type; options...)

Reconstruction without coil combination.
`type` can be :ONLINE or :PATREFSCAN
"""
function reconstruct(file; combine=:auto, kw...)
    image = reconstruct(file, :ONLINE; kw...)

    if combine == false || combine == :auto && size(image, 5) == 1
        return image
    end

    refscan = reconstruct(file, :PATREFSCAN; kw...)
    combined = coil_combine(image, refscan)
    return combined
end

function reconstruct(filename, type; datatype=ComplexF32, old_headers=false, mmap=true, kw...)
    data_headers, info = read_scan_info(filename, type, old_headers)
    image = mmaped_image(info, datatype, mmap)
    circle_array = sort_headers(data_headers, info)

    for circle in vcat(circle_array...)
        selectdim(image, 3, circle[:part]) .+= reconstruct(circle; datatype, kw...)
    end
    fft_slice_dim!(image)

    return image
end

# Returns [n_freq, n_phase, n_points, n_channels]
function reconstruct(c::Circle; datatype=ComplexF32, ice=false, do_fov_shift=true, do_freq_cor=true, do_dens_comp=true, conj_in_beginning=true)
    kspace_coordinates = datatype.(construct_circle_coordinates(c))
    kdata = read_data(c, datatype)

    if do_fov_shift; fov_shift!(kdata, kspace_coordinates, c) end
    if conj_in_beginning; kdata = conj.(kdata) end
    if do_freq_cor; frequency_offset_correction!(kdata, c) end

    if do_dens_comp
        if ice
            density_compensation_ice!(kdata, c)
        else
            density_compensation!(kdata, c)
        end
    end

    image = fourier_transform(kdata, kspace_coordinates, c[:n_frequency])

    image = reverse(image; dims=1) # LR flip

    return image
end
