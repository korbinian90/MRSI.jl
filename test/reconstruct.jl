using TestItemRunner
@testitem "reconstruct" begin
    file = "$(@__DIR__)/data/test.dat"
    image = reconstruct(file)
    @test size(image) == (16, 16, 1, 144, 1)
end
