using TestItemRunner
@testitem "read_headers" begin
    f_old16 = "C:/ICE/virtualShare/mrsi/dats/dats1/meas_MID00175_FID80315_csi_fidesi_crt_OldADC_Test8_16x16.dat"
    headers, info = read_scan_info(f_old16, :ONLINE)

    @test length(headers) == 208 # n_adcs in scan
    @test info[:n_fid] == 840
    @test headers[1] isa MRSI.ScanHeaderVD
end
