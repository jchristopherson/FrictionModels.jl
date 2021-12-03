using FrictionModels
using Test

@testset "CoulombModel" begin
    v = -0.5
    mdl = CoulombModel(0.25)
    N = 100.0
    F = friction(mdl, N, v)
    @test F.force == -sign(v) * N * mdl.coefficient
end

@testset "LuGreModel" begin
   v = -0.5
   mdl = LuGreModel(
       0.25,
       0.15,
       0.01,
       1.0e6,
       1.0e3,
       0.0
   )
   N = 100.0
   z = 0.0
   F = friction(mdl, N, v, z)
   @test F.force == z * mdl.bristle_stiffness
end
