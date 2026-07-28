using Test
using Random
using DTVerify

Random.seed!(20260521)

@testset "DTVerify" begin
    @testset "extract" begin include("test_extract.jl") end
    @testset "tests"   begin include("test_tests.jl") end
    @testset "decisions" begin include("test_decisions.jl") end
end
