# Lift interface: witnesses are expressed in the original Pauli algebra, not
# in realified/quotiented moment or solver coordinates. Verification does not
# trust a solver objective, PSD tolerance, or any claimed lift map.

function _wf_gram_psd!(Q)
    size(Q,1) == size(Q,2) && LinearAlgebra.ishermitian(Q) ||
        throw(ArgumentError("A rational Gram witness must be square and Hermitian."))
    for i in axes(Q,1)
        pivot = real(Q[i,i])
        pivot >= 0 || throw(ArgumentError("Gram witness is not positive semidefinite."))
        if iszero(pivot)
            all(iszero, view(Q,i+1:size(Q,1),i)) ||
                throw(ArgumentError("A zero Gram pivot has a nonzero column."))
        else
            for j in i+1:size(Q,1), k in i+1:size(Q,1)
                Q[j,k] -= Q[j,i]*conj(Q[k,i])/pivot
            end
        end
    end
    return nothing
end

function _wf_witness_polynomial(reference, p)
    p isa NC.Polynomial && p.registry === reference.registry ||
        throw(ArgumentError("A Pauli witness polynomial uses the wrong registry."))
    return _wf_exact(p)
end

function _wf_certificate_bound(p, level, sectors; squares=(), grams=(), sector_multipliers=(),
                               inequalities=(), window_multipliers=())
    λ = _PauliQ(level)
    residual = p - λ
    for f in squares
        a = _wf_witness_polynomial(p,f)
        residual -= a'*a
    end
    for witness in grams
        rows = [_wf_witness_polynomial(p,v) for v in witness.rows]
        Q = Matrix{_PauliCQ}(witness.gram)
        size(Q) == (length(rows),length(rows)) || throw(ArgumentError("Gram dimensions do not match its rows."))
        _wf_gram_psd!(copy(Q))
        for j in eachindex(rows), i in eachindex(rows)
            iszero(Q[i,j]) || (residual -= Q[i,j]*(rows[i]'*rows[j]))
        end
    end
    if !isempty(sector_multipliers)
        length(sector_multipliers) == length(sectors) || throw(ArgumentError("Provide one multiplier per sector."))
        for (s,u) in zip(sectors, sector_multipliers)
            a = s*_wf_witness_polynomial(p,u)
            residual -= (a+a')*(1//2)
        end
    end
    if !isempty(window_multipliers)
        length(window_multipliers) == length(inequalities) ||
            throw(ArgumentError("Provide one multiplier per energy-window inequality."))
        for (g, μ) in zip(inequalities, window_multipliers)
            q = _PauliQ(μ)
            q >= 0 || throw(ArgumentError("Energy-window multipliers must be nonnegative."))
            residual -= q*g
        end
    end
    NC._checked_hermitian(residual, "Certificate residual")
    # Every Hermitian Pauli is unitary; all residual coefficients are REAL
    # rationals. The norm enclosure and the returned bound are exact.
    enclosure = sum((abs(real(c)) for (c,_) in NC.terms(residual)); init=zero(_PauliQ))
    return (; bound=λ-enclosure, level=λ, residual_l1=enclosure,
            residual=Tuple(NC.terms(residual)), status=:verified_rational)
end

_wf_signature(p) = (p.registry, Tuple(NC.terms(p)))

struct PauliSOSCertificate{S,T,R}
    hamiltonian::S
    target::T
    bound::_PauliQ
    level::_PauliQ
    residual_l1::_PauliQ
    residual::R
    status::Symbol
    function PauliSOSCertificate(H, level; target_mz_sector=nothing, target_spin_j=nothing,
                                squares=(), grams=(), sector_multipliers=())
        HQ = _wf_exact(NC._checked_hermitian(H,"Hamiltonian"))
        target = _wf_target(HQ; mz=target_mz_sector, spin_j=target_spin_j)
        sig = _wf_signature(HQ)
        label = (; mz=target.mz, spin_j=target.spin_j)
        result = _wf_certificate_bound(HQ,level,target.sectors; squares,grams,sector_multipliers)
        return new{typeof(sig),typeof(label),typeof(result.residual)}(
            sig,label,result.bound,result.level,result.residual_l1,result.residual,result.status)
    end
end

"""
    certify_rational_sos(H, level; squares=(), grams=(), target_mz_sector=nothing,
                         target_spin_j=nothing, sector_multipliers=())

Verify a lifted SOS witness in the original Pauli algebra. `squares` contains
polynomials f for f†f; each `grams` entry is `(;rows,gram)` with a Hermitian
PSD rational matrix tested by exact Schur elimination. Binary floats are
interpreted as their exact stored values, not as nearby rational numbers. Use
rational factors to avoid introducing indefinite rounded Gram matrices. The
optional one-sided sector terms are Hermitian parts of (Q-q)u. The corrected
bound subtracts the exact ℓ₁ residual; it does not certify the uncorrected
level. An empty witness gives a norm bound, not a lifted solver dual.
"""
certify_rational_sos(H,level;kwargs...) = PauliSOSCertificate(H,level;kwargs...)

# Apply a Pauli polynomial exactly, without constructing a 2^N × 2^N matrix.
function _wf_apply(p, ket)
    out = zeros(_PauliCQ,length(ket))
    region = 1:p.registry.nsites
    for (c,w) in NC.terms(p)
        x,z = _wf_local_masks(w,region)
        for j in 0:length(ket)-1
            out[(j ⊻ x)+1] += c*_wf_pauli_phase(x,z,j)*ket[j+1]
        end
    end
    return out
end

struct PauliEnergyWindow{S,T}
    hamiltonian::S
    target::T
    lower::_PauliQ
    upper::_PauliQ
    function PauliEnergyWindow(H, cert::PauliSOSCertificate, ket::AbstractVector;
                               target_mz_sector=cert.target.mz, target_spin_j=cert.target.spin_j)
        HQ = _wf_exact(NC._checked_hermitian(H,"Hamiltonian"))
        signature = _wf_signature(HQ)
        signature == cert.hamiltonian || throw(ArgumentError("Energy certificate belongs to another Hamiltonian."))
        target = _wf_target(HQ; mz=target_mz_sector, spin_j=target_spin_j)
        label = (; mz=target.mz,spin_j=target.spin_j)
        cert.target in ((;mz=nothing,spin_j=nothing),label) ||
            throw(ArgumentError("Energy lower certificate does not cover this target."))
        big(2)^HQ.registry.nsites == length(ket) ||
            throw(ArgumentError("A variational ket needs one entry per computational basis vector."))
        a = _PauliCQ[x for x in ket]
        normsq = sum(abs2,a)
        normsq > 0 || throw(ArgumentError("The variational ket must be nonzero."))
        all(s -> all(iszero,_wf_apply(s,a)), target.sectors) ||
            throw(ArgumentError("The variational ket is outside the requested sector."))
        ub = real(LinearAlgebra.dot(a,_wf_apply(HQ,a)))/normsq
        return new{typeof(signature),typeof(label)}(signature,label,cert.bound,ub)
    end
end

"Construct a ground-energy enclosure from an exact lower certificate and an exact target-compatible Rayleigh quotient."
certify_energy_window(H, cert, ket; kwargs...) = PauliEnergyWindow(H,cert,ket;kwargs...)

_wf_window_bounds(window, H, target) = throw(ArgumentError("Pass a PauliEnergyWindow from certify_energy_window."))
function _wf_window_bounds(window::PauliEnergyWindow, H, target)
    _wf_signature(H) == window.hamiltonian || throw(ArgumentError("Energy window belongs to another Hamiltonian."))
    (;mz=target.mz,spin_j=target.spin_j) == window.target ||
        throw(ArgumentError("Energy window belongs to another sector."))
    return window.lower, window.upper
end

"""
    certify_observable_interval(H, O, window; lower, upper)

Verify endpoint witnesses for O and -O, respectively. Supply NamedTuples
with `level` and optional `squares`, `grams`, `sector_multipliers` and
`window_multipliers` (in order H-lo, hi-H). Returns exact rationals valid for
all sector-compatible states in the energy window, including non-invariant
states. The certificate never uses symmetry equalities, stationarity, or KMS.
No solver estimate is promoted to a verified endpoint without this check.
"""
function certify_observable_interval(H, O, window::PauliEnergyWindow; lower, upper)
    HQ = _wf_exact(NC._checked_hermitian(H,"Hamiltonian"))
    OQ = _wf_witness_polynomial(HQ, NC._checked_hermitian(O,"Observable"))
    target = _wf_target(HQ; window.target...)
    lo, hi = _wf_window_bounds(window,HQ,target)
    inequalities = (HQ-lo,hi-HQ)
    function endpoint(p; level, squares=(), grams=(), sector_multipliers=(), window_multipliers=())
        return _wf_certificate_bound(p,level,target.sectors;
            inequalities,squares,grams,sector_multipliers,window_multipliers)
    end
    lb, negub = endpoint(OQ;lower...),endpoint(-OQ;upper...)
    return (; lower=lb.bound, upper=-negub.bound, status=:verified_rational,
            energy_target=:energy_window, lower_certificate=lb, negated_upper_certificate=negub)
end
