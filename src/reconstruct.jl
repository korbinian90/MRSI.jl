"""
    csi = reconstruct(filename; combine=:auto, datatype=ComplexF32, ice=false, old_headers=false, mmap=true, lipid_decon=false, lipid_mask, brain_mask, do_fov_shift=true, do_freq_cor=true, do_dens_comp=true, conj_in_beginning=true)

Reconstructs a SIEMENS dat file

Coil combine is performed automatically for AC datasets. 
Set `combine` to `false`|`true` to force coil combination (including normalization by the patrefscan). 
Set `ice` to `true` to use calculated radii for density compensation instead of reading them from the headers. 
Set `old_headers` to `true` for older dat files with part and circle stored in LIN
Set `mmap` to `false` to compute in RAM or a filename or folder for storing the temporary result array.
Set `lipid_decon` to `:L1` or `:L2` to perform lipid decontamination. 
For lipid decontamination, pass the filenames of Float32 raw files via `lipid_mask` and `brain_mask`.
Alternatively, use `lipid_mask=:from_spectrum`. If nothing is provided, a full mask is used.

    csi, info = reconstruct(filename, type; options...)

Reconstruction without coil combination.
`type` can be `:ONLINE` or `:PATREFSCAN`
`info` is a `Dict` containing scan information.
"""
function reconstruct(file; combine=:auto, zero_fill=false, lipid_decon=nothing, kw...)
    csi, info = reconstruct(file, :ONLINE; kw...)

    if combine == true || combine == :auto && size(csi, 5) > 1
        refscan, _ = reconstruct(file, :PATREFSCAN; kw...)
        csi = coil_combine(csi, refscan)
    end

    if !isnothing(lipid_decon)
        mask, lipid_mask = get_masks(csi, info; kw...)
        lipid_suppression!(csi, mask, lipid_mask; type=lipid_decon)
    end

    if zero_fill
        csi = PaddedView(0, csi, (size(csi)[1:3]..., info[:vec_size], size(csi, 5)))
    end

    return csi
end

function reconstruct(filename, type; datatype=ComplexF32, old_headers=false, mmap=true, kw...)
    data_headers, info = read_scan_info(filename, type, old_headers)
    csi = mmaped_image(info, datatype, mmap)
    circle_array = sort_headers(data_headers, info)

    p = Progress(sum(length.(circle_array)))
    lk = ReentrantLock()
    for part in circle_array
        Threads.@threads for circle in part
            rec = reconstruct(circle; datatype, kw...)
            lock(lk) do
                # This needs to be guarded with a lock because of threaded addition
                selectdim(csi, 3, circle[:part]) .+= rec
            end
            next!(p)
        end
    end
    fft_slice_dim!(csi)

    return csi, info
end

# Returns [n_freq, n_phase, n_points, n_channels]
function reconstruct(c::Circle; datatype=ComplexF32, ice=false, do_fov_shift=true, do_freq_cor=true, do_dens_comp=true, conj_in_beginning=true, kw...)
    kspace_coordinates = datatype.(construct_circle_coordinates(c))
    kdata = read_data(c, datatype)

    if do_fov_shift
        fov_shift!(kdata, kspace_coordinates, c)
    end
    if conj_in_beginning
        kdata = conj.(kdata)
    end
    if do_freq_cor
        frequency_offset_correction!(kdata, c)
    end

    if do_dens_comp
        if ice
            density_compensation_ice!(kdata, c)
        else
            density_compensation!(kdata, c)
        end
    end

    csi = fourier_transform(kdata, kspace_coordinates, c[:n_frequency])

    csi = reverse(csi; dims=1) # LR flip

    return csi
end
