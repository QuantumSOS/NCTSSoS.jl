function mixed_contracts()
    @testset "Mixed products retain phases and contractions" begin
        _,parts = Next.mixed_variables(spin=(Pauli(),1),matter=(Fermionic(),1),osc=(Bosonic(),1))
        x,y,z = only.(parts.spin)
        a,ad = only.(parts.matter)
        b,bd = only.(parts.osc)
        # CAR and CCR each produce two terms, while XY contributes a phase.
        p = x*a*b*y*ad*bd
        @test p == im*z + im*z*bd*b - im*z*ad*a - im*z*ad*a*bd*b
        @test length(Next.terms(p)) == 4
        @test x*a*b*x*a*bd == zero(x)
        @test Next.degree(x*ad*a*bd*b) == 5
    end
end
