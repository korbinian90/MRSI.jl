"""
    csi = reconstruct(filename; combine=:auto, datatype=ComplexF32, ice=false, old_headers=false, mmap=true, lipid_decon=false, lipid_mask, brain_mask, do_fov_shift=true, do_freq_cor=true, do_dens_comp=true, conj_in_beginning=true, do_hamming_filter=true, do_hamming_filter_z=true, gradient_delay_us=[0,0,0])

Reconstructs a SIEMENS dat file

Coil combine is performed automatically for AC datasets. 
Set `combine` to `false`|`true` to force coil combination (including normalization by the patrefscan). 
Set `ice` to `true` to use calculated radii for density compensation instead of reading them from the headers. 
Set `old_headers` to `true` for older dat files with part and circle stored in LIN
Set `mmap` to `false` to compute in RAM or a filename or folder for storing the temporary result array.
Set `lipid_decon` to `:L1` or `:L2` to perform lipid decontamination.
Set `gradient_delay` to [[x1,y1], [x2,y2], [x3,y3]] to apply a gradient delay correction
For lipid decontamination, pass the filenames of Float32 raw files via `lipid_mask` and `brain_mask`.
Alternatively, use `lipid_mask=:from_spectrum`. If nothing is provided, a full mask is used.
Set settings for lipid decontamination: `L2_beta=0.1f0`, `L1_n_loops=5`

    csi, info = reconstruct(filename, type; options...)

Reconstruction without coil combination.
`type` can be `:ONLINE` or `:PATREFSCAN`
`info` is a `Dict` containing scan information.
"""
function reconstruct(file::AbstractString; time_point=5, combine=:auto, do_noise_decorrelation=false, zero_fill=false, lipid_decon=nothing, old_headers=false, do_hamming_filter_z=true, kw...)
    scan_info = read_scan_info(file, old_headers)

    if do_noise_decorrelation
        kw = Dict{Symbol,Any}(kw...) # required to modify kw
        kw[:noise_matrix_cholesky] = noise_decorrelation(scan_info)
    end

    csi = reconstruct(scan_info[:ONLINE]; kw...)

    if combine == true || combine == :auto && size(csi, 5) > 1
        refscan = reconstruct(scan_info[:PATREFSCAN]; kw...)
        csi = coil_combine(csi, refscan; time_point)
    end

    if do_hamming_filter_z && size(csi, 3) > 1
        csi = hamming_filter_z(csi)
    end

    if !isnothing(lipid_decon)
        mask, lipid_mask = get_masks(csi, scan_info[:ONLINE]; kw...)
        lipid_suppression!(csi, mask, lipid_mask; type=lipid_decon, kw...)
    end

    if zero_fill
        vec_size = scan_info[:ONLINE][:vec_size]
        csi = PaddedView(0, csi, (size(csi)[1:3]..., vec_size, size(csi, 5)))
    end

    return csi
end

function reconstruct(info::Dict; datatype=ComplexF32, mmap=true, kw...)
    csi = mmaped_image(info, datatype, mmap)
    circle_array = sort_into_circles(info[:headers], info)

    p = Progress(sum(length.(circle_array))) # Progress bar
    for part in circle_array
        for circle in part
            rec = reconstruct(circle; datatype, kw...)
            selectdim(csi, 3, circle[:part]) .+= rec
            next!(p) # Progress bar
        end
    end

    fft_slice_dim!(csi)

    return csi
end

# Returns [n_freq, n_phase, n_points, n_channels]
function reconstruct(c::Circle; noise_matrix_cholesky=nothing, datatype=ComplexF32, ice=false, do_fov_shift=true, do_freq_cor=true, do_dens_comp=true, do_hamming_filter=true, conj_in_beginning=true, gradient_delay_us=[0, 0, 0], kw...)
    kdata = read_data(c, datatype, noise_matrix_cholesky)

    if do_fov_shift
        fov_shift!(kdata, c, gradient_delay_us)
    end
    if conj_in_beginning
        kdata = conj.(kdata)
    end
    if do_freq_cor
        frequency_offset_correction!(kdata, c)
    end
    if do_dens_comp
        density_compensation!(kdata, c; do_hamming_filter, ice)
    end

    csi = fourier_transform(kdata, c, gradient_delay_us)

    csi = reverse(csi; dims=1) # LR flip

    return csi
end
