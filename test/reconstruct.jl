using TestItemRunner
@testitem "reconstruct" begin
    f_old16 = "C:/ICE/virtualShare/mrsi/dats/dats1/meas_MID00175_FID80315_csi_fidesi_crt_OldADC_Test8_16x16.dat"
    n_channels = 1
    fov_read = 220
    image = reconstruct(f_old16; n_channels, fov_read)
    @test size(image) == (16, 16, 1, 840, 1, 1)
end
