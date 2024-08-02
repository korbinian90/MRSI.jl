import Pkg
Pkg.activate(@__DIR__)

using Revise, MRSI, LinearAlgebra, MriResearchTools, Mmap, FFTW, Statistics

invivo_datfile = "/home/korbi/data/MRSI/fireICE_invivo1_1/dat/meas_MID00143_FID25434_csi_fidesi_crt_commit0dce.dat"

# reconstruct(invivo_datfile; lipid_decon=:L2, save=savenii)

outfolder = "/home/korbi/data/MRSI/fireICE_invivo1_1/lipid_regularisation_test/"

scan_info = read_scan_info(invivo_datfile)
info = scan_info[:ONLINE]

csi = reconstruct(scan_info[:ONLINE])
# io = open("/home/korbi/julia_tmp_mmap/csi.nii", "r")
# csi = mmap(io, Array{ComplexF32,5}, (info[:n_frequency], info[:n_frequency], info[:n_part], info[:n_fid], info[:n_channels]))
refscan = reconstruct(scan_info[:PATREFSCAN])

brain_mask = MRSI.brain_mask_from_csi(csi, refscan, scan_info[:ONLINE])
combined = coil_combine(csi, refscan; ref_point_for_combine)
savenii(comb)


MRSI.lipid_suppression!(csi, brain_mask, scan_info[:ONLINE]; type=lipid_decon, kw...)

combined = coil_combine(csi, refscan; ref_point_for_combine)

##

## single channel
outfolder = "/home/korbi/data/MRSI/fireICE_invivo1_1/masking_channelwise/"
mkpath(outfolder)
data = csi[:,:,:,:,1]

dims = 4
spectrum = abs.(fftshift(fft(data, dims), dims))

fat_range_ppm = [1.8, 0.5]
fat_range = MRSI.ppm_to_vecsize_point.(Ref(info), fat_range_ppm)
fat_max = maximum(selectdim(spectrum, dims, fat_range[1]:fat_range[2]); dims)
savenii(fat_max, "fat_max_1.nii.gz", outfolder)

water_range_ppm = [5.0, 4.3]
water_range = MRSI.ppm_to_vecsize_point.(Ref(info), water_range_ppm)
water_max = maximum(selectdim(spectrum, dims, water_range[1]:water_range[2]); dims)
savenii(water_max, "water_max_1.nii.gz", outfolder)


all = maximum(spectrum; dims)
savenii(all, "all_1.nii.gz", outfolder)

savenii(water_max ./ all, "water_ratio_1.nii.gz", outfolder)
savenii(fat_max ./ all, "fat_ratio_1.nii.gz", outfolder)
savenii(water_max ./ fat_max, "fat_water_ratio_1.nii.gz", outfolder)

fat_ratio_mask = fat_max ./ all .< 0.1
MRSI.cut_to_ellipse!(fat_ratio_mask)
savenii(fat_ratio_mask, "fat_ratio_mask_1.nii.gz", outfolder)

savenii(sum(spectrum; dims), "all_sum_1.nii.gz", outfolder)

savenii(spectrum, "spec.nii.gz", outfolder)
savenii(abs.(data), "fids.nii.gz", outfolder)


## multi-channels
function easy_combine(csi; time_point=5)
    combined = zeros(eltype(csi), size(csi)[1:4])
    for channel in axes(csi, 5)
        csi_channel = selectdim(csi, 5, channel)
        combined .+= csi_channel ./ (csi_channel[:,:,:,time_point] ./ abs.(csi_channel[:,:,:,time_point]))
    end
    return combined
end

data = easy_combine(csi)

outfolder = "/home/korbi/data/MRSI/fireICE_invivo1_1/masking_easy_comb/"
mkpath(outfolder)

dims = 4
spectrum = abs.(fftshift(fft(data, dims), dims))

fat_range_ppm = [1.8, 0.5]
fat_range = MRSI.ppm_to_vecsize_point.(Ref(info), fat_range_ppm)
fat_max = maximum(selectdim(spectrum, dims, fat_range[1]:fat_range[2]); dims)
savenii(fat_max, "fat_max_1.nii.gz", outfolder)

water_range_ppm = [5.0, 4.3]
water_range = MRSI.ppm_to_vecsize_point.(Ref(info), water_range_ppm)
water_max = maximum(selectdim(spectrum, dims, water_range[1]:water_range[2]); dims)
savenii(water_max, "water_max_1.nii.gz", outfolder)


all = maximum(spectrum; dims)
savenii(all, "all_1.nii.gz", outfolder)

savenii(water_max ./ all, "water_ratio_1.nii.gz", outfolder)
savenii(fat_max ./ all, "fat_ratio_1.nii.gz", outfolder)
savenii(water_max ./ fat_max, "fat_water_ratio_1.nii.gz", outfolder)

fat_ratio_mask = fat_max ./ all .< 0.1
MRSI.cut_to_ellipse!(fat_ratio_mask)
savenii(fat_ratio_mask, "fat_ratio_mask_1.nii.gz", outfolder)

savenii(sum(spectrum; dims), "all_sum_1.nii.gz", outfolder)

savenii(spectrum, "spec.nii.gz", outfolder)
savenii(abs.(data), "fids.nii.gz", outfolder)


## Load CombinedCSI.mat
using MAT

mat = matopen("/home/korbi/data/MRSI/fireICE_invivo1_1/CombinedCSI.mat")
mat_csi = read(mat, "csi")
close(mat)


data = mat_csi["Data"]
outfolder = "/home/korbi/data/MRSI/fireICE_invivo1_1/masking_mat/"
mkpath(outfolder)

dims = 4
spectrum = abs.(fftshift(fft(data, dims), dims))

fat_range_ppm = [1.8, 0.5]
fat_range = MRSI.ppm_to_vecsize_point.(Ref(info), fat_range_ppm)
fat_max = maximum(selectdim(spectrum, dims, fat_range[1]:fat_range[2]); dims)
savenii(fat_max, "fat_max_1.nii", outfolder)

water_range_ppm = [5.0, 4.3]
water_range = MRSI.ppm_to_vecsize_point.(Ref(info), water_range_ppm)
water_max = maximum(selectdim(spectrum, dims, water_range[1]:water_range[2]); dims)
savenii(water_max, "water_max_1.nii", outfolder)


all = maximum(spectrum; dims)
savenii(all, "all_1.nii", outfolder)

savenii(water_max ./ all, "water_ratio_1.nii", outfolder)
savenii(fat_max ./ all, "fat_ratio_1.nii", outfolder)
savenii(water_max ./ fat_max, "fat_water_ratio_1.nii", outfolder)

fat_ratio_mask = fat_max ./ all .< 0.1
MRSI.cut_to_ellipse!(fat_ratio_mask)
savenii(fat_ratio_mask, "fat_ratio_mask_1.nii", outfolder)

savenii(sum(spectrum; dims), "all_sum_1.nii", outfolder)

savenii(spectrum, "spec.nii", outfolder)
savenii(abs.(data), "fids.nii", outfolder)

##
header = niread("/home/korbi/data/MRSI/fireICE_invivo1_1/maps/Orig/Glu_amp_map.nii").header
for channel in 1:4:32
    csi_channel = selectdim(csi, 5, channel)
    # savenii(abs.(fftshift(fft(csi_channel, 4), 4)), "spectrum_$channel.nii", "/home/korbi/data/MRSI/fireICE_invivo1_1/lipid_mask_channelwise/")
    lipid, all = MRSI.lipid_mask_from_mrsi(csi_channel, info)
    lipid_mask1 = lipid ./ all .> 0.9
    
    lipid_mask2 = lipid .> @show 0.8 * mean(lipid[lipid_mask1])

    lipid_mask_easy = lipid .> mean(lipid)
    lipid_mask_easy2 = lipid .> mean(lipid_mask_easy)

    combined_mask = lipid_mask1 .& lipid_mask2
    MRSI.cut_to_ellipse!(combined_mask)
    MRSI.cut_to_ellipse!(lipid_mask1)
    MRSI.cut_to_ellipse!(lipid_mask2)
    MRSI.cut_to_ellipse!(lipid_mask_easy)
    MRSI.cut_to_ellipse!(lipid_mask_easy2)
    

    savenii(lipid, "max_lipid_$channel.nii.gz", "/home/korbi/data/MRSI/fireICE_invivo1_1/lipid_mask_channelwise/", header)
    savenii(lipid_mask1, "lipid_mask1_$channel.nii.gz", "/home/korbi/data/MRSI/fireICE_invivo1_1/lipid_mask_channelwise/", header)
    savenii(lipid_mask2, "lipid_mask2_$channel.nii.gz", "/home/korbi/data/MRSI/fireICE_invivo1_1/lipid_mask_channelwise/", header)
    savenii(lipid_mask_easy, "lipid_mask_easy_$channel.nii.gz", "/home/korbi/data/MRSI/fireICE_invivo1_1/lipid_mask_channelwise/", header)
    savenii(lipid_mask_easy2, "lipid_mask_easy2_$channel.nii.gz", "/home/korbi/data/MRSI/fireICE_invivo1_1/lipid_mask_channelwise/", header)
    savenii(combined_mask, "lipid_mask_combined_$channel.nii.gz", "/home/korbi/data/MRSI/fireICE_invivo1_1/lipid_mask_channelwise/", header)
end

https://niivue.github.io/niivue-vscode/?images=https://korbinian90.github.io/Shared_Data/lipid_mask_easy2_1.nii.gz,https://korbinian90.github.io/Shared_Data/lipid_mask_easy2_5.nii.gz,https://korbinian90.github.io/Shared_Data/lipid_mask_easy2_9.nii.gz,https://korbinian90.github.io/Shared_Data/lipid_mask_easy2_13.nii.gz
https://niivue.github.io/niivue-vscode/?images=https://korbinian90.github.io/Shared_Data/lipid_mask_combined_1.nii.gz,https://korbinian90.github.io/Shared_Data/lipid_mask_combined_5.nii.gz,https://korbinian90.github.io/Shared_Data/lipid_mask_combined_9.nii.gz,https://korbinian90.github.io/Shared_Data/lipid_mask_combined_13.nii.gz

https://niivue.github.io/niivue-vscode/?images=https://korbinian90.github.io/Shared_Data/lipid_mask_easy2_1.nii.gz,https://korbinian90.github.io/Shared_Data/lipid_mask_easy2_5.nii.gz,https://korbinian90.github.io/Shared_Data/lipid_mask_easy2_9.nii.gz,https://korbinian90.github.io/Shared_Data/lipid_mask_easy2_13.nii.gz,https://korbinian90.github.io/Shared_Data/lipid_mask_combined_1.nii.gz,https://korbinian90.github.io/Shared_Data/lipid_mask_combined_5.nii.gz,https://korbinian90.github.io/Shared_Data/lipid_mask_combined_9.nii.gz,https://korbinian90.github.io/Shared_Data/lipid_mask_combined_13.nii.gz
# datfile = "/data/ICE/virtualShare/mrsi/dats/noise_decor/meas_MID00156_FID107701_csi_fidesi_crt_OldADC_Test8_8x8_AC_2D.dat"
# datfile = "/home/korbi/data/ICE/virtualShare/mrsi/dats/merged/AC2D/meas_MID00127_FID125462_lh_csi_fid_crt_meldBerni_5mai_AC_2D.dat"
datfile = "/home/korbi/data/ICE/virtualShare/mrsi/dats/dats_rotated/1/meas_MID00153_FID127085_lh_csi_fid_crt_meldBerni_4mai_vc3d_rot.dat"

scan_info = read_scan_info(datfile)[:ONLINE];
scan_info[:max_r]
scan_info[:radii]

reconstruct(datfile; do_hamming_filter=true);


## search protocol
prot = MRSI.read_twix_protocol(datfile)

# search all entries of the prot dictionary if it contains the string "MemBlock". It is a nested dictionary
for (k, v) in prot
    for (kk, vv) in v
        if occursin("MemBlock", string(kk))
            println(k, " ", kk, " ", vv)
        end 
    end
end

prot["MeasYaps"]["sWipMemBlock"][""]
##
scan_info = read_scan_info(datfile)

scan_info[:NOISADJSCAN][:headers]

[h[:n_TI] for h in MRSI.sort_into_circles(scan_info[:ONLINE][:headers], scan_info[:ONLINE])[1]]

scan_info[:ONLINE][:headers][1][:n_points_on_circle]

circles = MRSI.sort_into_circles(scan_info[:ONLINE][:headers], scan_info[:ONLINE])
circles[1][1][:n_points_on_circle]

res = open(datfile, "r") do io
    MRSI.read_adc(io, scan_info[:NOISADJSCAN][:headers][1].data_position, 32)
    vcat((MRSI.read_adc(io, head.data_position, 32) for head in scan_info[:NOISADJSCAN][:headers])...)
end




noise = MRSI.read_noise(scan_info[:NOISADJSCAN])

@time noise_correlation = MRSI.calculate_noise_correlation(noise)

# noise_scaling_factor = noise_dwelltime / (reco_dwelltime * nTempIntsPerAngInt/datasize{1}(1))
noise_scaling_factor = 1.25

scan_info[:NOISADJSCAN][:dwelltime]
scan_info[:ONLINE][:dwelltime]

##

noise_matrix_cholesky = cholesky(noise_correlation * noise_scaling_factor / 2).U

data = rand(ComplexF64, 840, 32)

circle_array = MRSI.sort_into_circles(scan_info[:ONLINE][:headers], scan_info[:ONLINE])

MRSI.read_data(circle_array[1][1], ComplexF64, noise_matrix_cholesky)

##
reconstruct(datfile; do_noise_decorrelation=true)
