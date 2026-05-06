using Test, NCTSSoS

@testset "Moment linear form cache primitives" begin
    keys = [Int[2], Int[], Int[1, 2], Int[1]]
    @test sort(keys; lt=NCTSSoS.key_lt) == [Int[], Int[1], Int[1, 2], Int[2]]
    @test NCTSSoS.key_lt(Int[1], Int[1, 0])
    @test !NCTSSoS.key_lt(Int[1, 0], Int[1])
    @test NCTSSoS.key_isequal(Int[3, 4], Int[3, 4])

    form = NCTSSoS.LinearMomentForm{Vector{Int},Float64}([
        Int[2] => 1,
        Int[1] => 2,
        Int[2] => -0.25,
        Int[1] => -2,
        Int[3] => 0,
    ])

    @test collect(form) == [Int[2] => 0.75]
    @test length(form) == 1
    @test !isempty(form)
    @test NCTSSoS._assert_linear_moment_form_invariants(form) === nothing

    complex_form = NCTSSoS.LinearMomentForm{Vector{Int},ComplexF64}([
        Int[2] => 1,
        Int[1] => im,
        Int[1] => -im,
    ])
    @test collect(complex_form) == [Int[2] => 1.0 + 0.0im]
end
