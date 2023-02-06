using TestItemRunner
@testitem "read_data" begin
    f_old16 = "C:/ICE/virtualShare/mrsi/dats/dats1/meas_MID00175_FID80315_csi_fidesi_crt_OldADC_Test8_16x16.dat"
    scaninfo = ScanInfo(f_old16, :ONLINE)

    sliceinfo = first(scaninfo)
    slice = open(f_old16) do io
        read_slice(io, sliceinfo)
    end

    re = rearrange(slice, sliceinfo)
    @test size(re) == (960, 840, 1) # (points_per_slice, n_fid, n_channels)
end