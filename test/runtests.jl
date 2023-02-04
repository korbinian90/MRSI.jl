using MRSI
using Test
using TestItemRunner

@run_package_tests

include("read_headers.jl")

@testitem "MRSI.jl" begin
    f_old16 = "C:/ICE/virtualShare/mrsi/dats/dats1/meas_MID00175_FID80315_csi_fidesi_crt_OldADC_Test8_16x16.dat"
    channels = 1
    headers = read_data_headers(f_old16)

    h = headers[:ONLINE][1]
    slice = open(f_old16) do io
        read_slice(io, h, 1)
    end
    rearrange(slice, h, 1)


end

@testitem "Fourier Transform" begin
    f_old16 = "C:/ICE/virtualShare/mrsi/dats/dats1/meas_MID00175_FID80315_csi_fidesi_crt_OldADC_Test8_16x16.dat"
    channels = 1
    n_grid = 16
    fov_read = 220
    headers = read_data_headers(f_old16)

    h = headers[:ONLINE][1]
    slice = open(f_old16) do io
        read_slice(io, h, 1)
    end

    ordered_kspace = rearrange(slice, h, 1)

    kspace_points = MRSI.kspace_coordinates(h, n_grid, fov_read)
    ft = fourier_transform(ordered_kspace, kspace_points, n_grid)
    @test size(ft) == (16, 16, 1, 840, 1, 1) # (n_grid, n_grid, slice, n_fid, part, n_channels)
end
