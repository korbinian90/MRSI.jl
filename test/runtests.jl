using MRSI
using Test
using TestItemRunner

@run_package_tests

include("read_headers.jl")

@testitem "MRSI.jl" begin
    f_old16 = "C:/ICE/virtualShare/mrsi/dats/dats1/meas_MID00175_FID80315_csi_fidesi_crt_OldADC_Test8_16x16.dat"
    data_headers, info = read_scan_info(f_old16, :ONLINE)
    image = MRSI.mmaped_image(info)

    MRSI.read_and_reconstruct_image_per_circle!(image, data_headers, info)
    MRSI.fft_slice_dim!(image)

end

@testitem "Fourier Transform" begin
    f_old16 = "C:/ICE/virtualShare/mrsi/dats/dats1/meas_MID00175_FID80315_csi_fidesi_crt_OldADC_Test8_16x16.dat"
    headers, info = read_scan_info(f_old16, :ONLINE)

    all_circles = MRSI.initialize_header_storage(info)
    for h in headers
        MRSI.store!(all_circles, h, info)
    end

    slice = 1
    # reconstruct all circles of one slice in one fourier transform (MATLAB order)
    kdata = vcat((MRSI.read_data(c) for c in all_circles[slice])...)
    kspace_coordinates = vcat((MRSI.construct_circle_coordinates(c) for c in all_circles[slice])...)
    image = MRSI.fourier_transform(kdata, kspace_coordinates, info[:n_frequency])

    # reconstruct circle per circle (ICE order)
    image2 = similar(image)
    for c in all_circles[slice]
        kdata = MRSI.read_data(c)
        kspace_coordinates = MRSI.construct_circle_coordinates(c)
        image2 .+= MRSI.fourier_transform(kdata, kspace_coordinates, info[:n_frequency])
    end

    @test size(image) == (16, 16, 840, 1) # (n_grid, n_grid, slice, n_fid, n_channels)
    @test all(image .≈ image2)
end
