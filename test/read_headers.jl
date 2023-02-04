using TestItemRunner
@testitem "read_headers" begin
    f_old16 = "C:/ICE/virtualShare/mrsi/dats/dats1/meas_MID00175_FID80315_csi_fidesi_crt_OldADC_Test8_16x16.dat"
    headers = read_data_headers(f_old16)
    @test length(headers[:ONLINE][1]) == 8 # n_circles
    @test size(headers[:ONLINE][1][1]) == (26, 1, 1) # ()
    @test headers[:ONLINE][1][1][1] isa MRSI.ScanHeaderVD
end
