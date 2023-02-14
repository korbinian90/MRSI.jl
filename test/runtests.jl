using MRSI
using Test
using TestItemRunner

@run_package_tests

include("read_headers.jl")

@testitem "MRSI.jl" begin
    f_old16 = "C:/ICE/virtualShare/mrsi/dats/dats1/meas_MID00175_FID80315_csi_fidesi_crt_OldADC_Test8_16x16.dat"
    scaninfo = ScanInfo(f_old16, :ONLINE)

    sliceinfo = first(scaninfo)
    slice = open(f_old16) do io
        read_slice(io, sliceinfo)
    end
    rearrange(slice, sliceinfo)


end

@testitem "Fourier Transform" begin
    f_old16 = "C:/ICE/virtualShare/mrsi/dats/dats1/meas_MID00175_FID80315_csi_fidesi_crt_OldADC_Test8_16x16.dat"
    scaninfo = ScanInfo(f_old16, :ONLINE)
    n_grid = scaninfo[:n_frequency]

    sliceinfo = first(scaninfo)
    slice = open(f_old16) do io
        read_slice(io, sliceinfo)
    end

    ordered_kspace = vcat(rearrange(slice, sliceinfo)...)

    kspace_points = MRSI.kspace_coordinates(sliceinfo)
    ft = fourier_transform(ordered_kspace, kspace_points, n_grid)
    @test size(ft) == (16, 16, 840, 1) # (n_grid, n_grid, slice, n_fid, n_channels)
end
