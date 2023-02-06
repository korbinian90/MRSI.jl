using TestItemRunner
@testitem "read_headers" begin
    f_old16 = "C:/ICE/virtualShare/mrsi/dats/dats1/meas_MID00175_FID80315_csi_fidesi_crt_OldADC_Test8_16x16.dat"
    scaninfo = ScanInfo(f_old16, :ONLINE)
    @test length(scaninfo[1]) == 8 # n_circles
    @test size(scaninfo[1][1]) == (26, 1) # ()
    @test scaninfo[1][1][1] isa MRSI.ScanHeaderVD
end
