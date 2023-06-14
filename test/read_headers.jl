using TestItemRunner
@testitem "read_headers" begin
    file = "$(@__DIR__)/data/test.dat"
    info = read_scan_info(file)[:ONLINE]
    headers = info[:headers]


    @test length(headers) == 96 # n_adcs in scan
    @test info[:n_fid] == 144
    @test headers[1] isa MRSI.ScanHeaderVD

    @test info[:max_n_circles] == 8
    @test info[:position] == [4.0, 10.0, 17.0]
end
