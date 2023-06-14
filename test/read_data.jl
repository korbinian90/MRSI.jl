using TestItemRunner
@testitem "read_data" begin
    file = "$(@__DIR__)/data/test.dat"
    info = read_scan_info(file)[:ONLINE]

    all_circles = MRSI.initialize_header_storage(info)
    for h in info[:headers]
        MRSI.store!(all_circles, h, info)
    end
    kdata = [vcat((MRSI.read_data(c, ComplexF32) for c in part_circles)...) for part_circles in all_circles]

    @test size(kdata) == (1,) # one slice
    @test size(kdata[1]) == (960, 144, 1) # (points_per_slice, n_fid, n_channels)
end
