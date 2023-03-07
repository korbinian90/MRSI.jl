using TestItemRunner
@testitem "read_data" begin
    f_old16 = "C:/ICE/virtualShare/mrsi/dats/dats1/meas_MID00175_FID80315_csi_fidesi_crt_OldADC_Test8_16x16.dat"
    headers, info = read_scan_info(f_old16, :ONLINE)

    all_circles = MRSI.initialize_header_storage(info)
    for h in headers
        MRSI.store!(all_circles, h, info)
    end

    kdata = [vcat((MRSI.read_data(c) for c in part_circles)...) for part_circles in all_circles]

    @test size(kdata) == (1,) # one slice
    @test size(kdata[1]) == (960, 840, 1) # (points_per_slice, n_fid, n_channels)
end
