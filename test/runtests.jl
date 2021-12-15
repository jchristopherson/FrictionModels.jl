using FrictionModels
using Test

@testset "CoulombModel" begin
    v = -0.5
    mdl = CoulombModel(0.25)
    N = 100.0
    F = friction(mdl, N, v)
    @test F.f == sign(v) * N * mdl.coefficient
end

