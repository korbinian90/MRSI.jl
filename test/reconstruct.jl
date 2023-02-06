using TestItemRunner
@testitem "reconstruct" begin
    f_old16 = "C:/ICE/virtualShare/mrsi/dats/dats1/meas_MID00175_FID80315_csi_fidesi_crt_OldADC_Test8_16x16.dat"
    image = reconstruct(f_old16)
    @test size(image) == (16, 16, 1, 840, 1)
end
