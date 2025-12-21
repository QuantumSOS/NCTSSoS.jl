# StateWord and NCStateWord Basis Generation Tests
#
# These tests verify that basis generation for trace polynomial optimization
# produces the correct number and uniqueness of basis elements.

using Test
using NCTSSoS.FastPolynomials
using NCTSSoS.FastPolynomials:
    StateWord,
    NCStateWord,
    MaxEntangled,
    Arbitrary,
    get_state_basis,
    get_ncbasis,
    cyclic_symmetric_canon,
    _state_canon,
    degree,
    tr

@testset "StateWord hash consistency" begin
    reg, (x,) = create_unipotent_variables([("x", 1:2)])
    
    @testset "Basic equality" begin
        sw1 = tr(x[1])
        sw2 = tr(x[1])
        @test sw1 == sw2
        @test hash(sw1) == hash(sw2)
    end
    
    @testset "After multiplication" begin
        sw1 = tr(x[1])
        sw3 = sw1 * sw1  # tr(x1) * tr(x1) = tr(x1)^2
        sw4 = tr(x[1]) * tr(x[1])
        @test sw3 == sw4
        @test hash(sw3) == hash(sw4)
    end
    
    @testset "Cyclic equivalence for MaxEntangled" begin
        # For MaxEntangled (trace) states, cyclic permutations should be equivalent
        # tr(x1*x2) = tr(x2*x1)
        m12 = x[1] * x[2]  # This is Monomial multiplication
        m21 = x[2] * x[1]
        
        sw12 = tr(m12)
        sw21 = tr(m21)
        
        @test sw12 == sw21
        @test hash(sw12) == hash(sw21)
    end
end

@testset "NCStateWord hash consistency" begin
    reg, (x,) = create_unipotent_variables([("x", 1:2)])
    
    @testset "Basic equality" begin
        sw = tr(x[1])
        nc = x[2]
        ncsw1 = NCStateWord(sw, nc)
        ncsw2 = NCStateWord(sw, nc)
        @test ncsw1 == ncsw2
        @test hash(ncsw1) == hash(ncsw2)
    end
    
    @testset "After multiplication" begin
        sw = tr(x[1])
        nc = x[2]
        ncsw1 = NCStateWord(sw, nc)
        ncsw3 = ncsw1 * ncsw1
        ncsw4 = NCStateWord(sw, nc) * NCStateWord(sw, nc)
        @test ncsw3 == ncsw4
        @test hash(ncsw3) == hash(ncsw4)
    end
    
    @testset "Different NCStateWords" begin
        sw = tr(x[1])
        ncsw1 = NCStateWord(sw, x[1])
        ncsw2 = NCStateWord(sw, x[2])
        @test ncsw1 != ncsw2
    end
end

@testset "get_state_basis uniqueness" begin
    reg, (vars,) = create_unipotent_variables([("v", 1:4)])
    
    basis = get_state_basis(reg, 2; state_type=MaxEntangled)
    
    @testset "No duplicates in basis" begin
        @test length(basis) == length(unique(basis))
    end
    
    @testset "Hash consistency for all pairs" begin
        found_equal = false
        for i in 1:length(basis), j in i+1:length(basis)
            if basis[i] == basis[j]
                @test hash(basis[i]) == hash(basis[j])
                found_equal = true
            end
        end
        # We shouldn't find any equal pairs since unique! was applied
        @test !found_equal
    end
end

@testset "get_state_basis size matches NCTSSOS" begin
    # This test verifies that our basis generation matches NCTSSOS
    # for the trace polynomial optimization case
    
    @testset "4 variables, order 2" begin
        # Using all variables on same site to match NCTSSOS behavior
        reg, (vars,) = create_unipotent_variables([("v", 1:4)])
        
        basis = get_state_basis(reg, 2; state_type=MaxEntangled)
        
        # NCTSSOS produces wbasis of size 53 for n=4, d=2 with binary=true
        # This is calculated by:
        # - tbasis: 21 elements (trace basis up to degree 2)
        # - basis: 17 elements (nc-word basis up to degree 2)
        # - Combined with degree constraint: 53 total
        @test length(basis) == 53
    end
    
    @testset "2 variables, order 2" begin
        reg, (vars,) = create_unipotent_variables([("v", 1:2)])
        
        basis = get_state_basis(reg, 2; state_type=MaxEntangled)
        
        # For n=2, d=2: should have specific count
        # tbasis: 1 (identity) + 2 (single) + 1 (two-length) + 3 (pairs) = 7
        # basis: 1 + 2 + 2 = 5 (with U^2=I simplification)
        # Combined: 15
        @test length(basis) == 15
    end
end

@testset "Unipotent simplification with site-based commutation" begin
    # Test that operators on different sites commute (sorted by site)
    
    reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
    
    @testset "Cross-site monomials are distinct before simplify" begin
        # y[1]*x[1] and x[1]*y[1] should be different operators before simplification
        m_yx = y[1] * x[1]
        m_xy = x[1] * y[1]
        
        @test m_yx != m_xy  # These should be different monomials before simplify
    end
    
    @testset "Cross-site monomials become equal after simplify (site commutation)" begin
        using NCTSSoS.FastPolynomials: simplify
        
        # After simplification, operators on different sites commute
        # so y[1]*x[1] and x[1]*y[1] should become the same (sorted by site)
        m_yx = y[1] * x[1]
        m_xy = x[1] * y[1]
        
        m_yx_simp = simplify(m_yx)
        m_xy_simp = simplify(m_xy)
        
        @test m_yx_simp == m_xy_simp  # Should be equal after site-based sorting
    end
    
    @testset "Basis size with site commutation" begin
        # With site-based commutation, basis sizes depend on site structure:
        # - Single site (4 vars): all operators within-site, no cross-site commutation
        # - Multi site (2+2 vars): cross-site operators commute, reducing some basis elements
        
        reg_single, (vars,) = create_unipotent_variables([("v", 1:4)])
        basis_single = get_state_basis(reg_single, 2; state_type=MaxEntangled)
        
        reg_multi, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
        basis_multi = get_state_basis(reg_multi, 2; state_type=MaxEntangled)
        
        # Multi-site has smaller basis due to cross-site commutation
        @test length(basis_single) == 53
        @test length(basis_multi) == 49
        @test length(basis_multi) < length(basis_single)
    end
end
