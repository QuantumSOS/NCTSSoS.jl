"""
    real_program(program)

Convert the reduced Hermitian blocks to real PSD cones. Already real blocks
keep their side. A genuinely complex block uses [Re(A) -Im(A); Im(A) Re(A)].
This fallback doubles that reduced block; it does not claim an antiunitary
fixed-basis reduction. No unreduced site-space matrix is constructed.
"""
function real_program(p::MomentProgram)
    out = MomentProgram{Float64}(p.nvars)
    append!(out.objective_ids, p.objective_ids)
    append!(out.objective_coefs, p.objective_coefs)
    for b in p.blocks
        b.kind in (PSD, HPSD) || throw(ArgumentError("The M1 real conversion expects only PSD blocks."))
        is_real = all(isreal(p.forms.coefs[t]) for f in form_range(b) for t in term_range(p.forms,f))
        n, side = b.side, is_real ? b.side : 2b.side
        for j in 1:side, i in j:side
            row, col = mod1(i,n), mod1(j,n)
            imaginary = i > n && j <= n
            conjugated = row < col
            f = b.first_form + tri_offset(n, max(row,col), min(row,col))
            for t in term_range(p.forms,f)
                c = p.forms.coefs[t]
                coefficient = imaginary ? (conjugated ? -imag(c) : imag(c)) : real(c)
                iszero(coefficient) || push_term!(out.forms, coefficient, p.forms.ids[t])
            end
            close_form!(out.forms)
        end
        push_block!(out, PSD, side, b.label)
    end
    validate(out)
end

_anti_real_program(p, generators) = any(is_antiunitary, generators) ? real_program(p) : p
