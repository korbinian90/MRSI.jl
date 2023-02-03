using TestItemRunner
@testitem "read_data" begin
    f_old16 = "C:/ICE/virtualShare/mrsi/dats/dats1/meas_MID00175_FID80315_csi_fidesi_crt_OldADC_Test8_16x16.dat"
    channels = 1
    headers = read_data_headers(f_old16)

    h = headers[:ONLINE][1]
    slice = open(f_old16) do io
        read_slice(io, h, 1)
    end
    re = MRSI.rearrange_circle(slice[6], h[6], 1)
    @show size(re)

    re = rearrange(slice, h, 1)
    @test size(re) == (960, 1, 840, 1) # (points_per_slice, par, n_fid, n_channels)
end