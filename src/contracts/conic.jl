# =============================================================================
# Contract 2: blocks of linear forms on moment ids
# =============================================================================
#
# Provider: the pipeline layer.
# User:     the conic layer (moment form, SOS form, solver writers).
#
# This is a data format. The conic layer has no knowledge of algebras, words or
# symmetries. A change to this file needs the agreement of the architect thread.

"""
    LinearForms{C}

A list of affine forms on the real moment variables `y[1], ..., y[nvars]`, in
compressed row storage. Form `k` is

    sum(coefs[j] * y[ids[j]] for j in offsets[k]:(offsets[k + 1] - 1))

`y[1]` is the constant 1 (the moment of the identity), so a term with id 1 is
the constant part of a form. A coefficient type `C <: Complex` gives complex
forms on real variables.

All forms of all blocks of a program are in one `LinearForms`, so a build makes
three arrays and no object for each entry.
"""
struct LinearForms{C<:Number}
    offsets::Vector{Int}
    ids::Vector{Int32}
    coefs::Vector{C}
end

LinearForms{C}() where {C<:Number} = LinearForms{C}(Int[1], Int32[], C[])

Base.length(f::LinearForms) = length(f.offsets) - 1
term_range(f::LinearForms, k::Integer) = f.offsets[k]:(f.offsets[k + 1] - 1)

"""
    push_term!(forms::LinearForms, coefficient, id::Integer)

Add one term to the form that is open at the end of `forms`. `close_form!` ends
the open form. An empty form is the zero form.
"""
function push_term!(f::LinearForms{C}, coefficient::Number, id::Integer) where {C}
    push!(f.ids, Int32(id))
    push!(f.coefs, convert(C, coefficient))
    return f
end

"""
    close_form!(forms::LinearForms) -> Int

End the open form and return its index.
"""
function close_form!(f::LinearForms)
    push!(f.offsets, length(f.ids) + 1)
    return length(f)
end

"""
    ConeKind

- `ZERO`: each form is zero (real part and imaginary part).
- `NONNEGATIVE`: each form is real and at least zero.
- `SOC`: the forms `(t, x...)` are real and satisfy `norm(x) <= t`.
- `PSD`: the real symmetric matrix of the forms is positive semidefinite.
- `HPSD`: the complex Hermitian matrix of the forms is positive semidefinite.
"""
@enum ConeKind::UInt8 ZERO NONNEGATIVE SOC PSD HPSD

"""
    ConeBlock(kind, side, first_form, label)

One cone constraint on consecutive forms of a program.

- `ZERO`, `NONNEGATIVE`, `SOC`: `side` forms.
- `PSD`, `HPSD`: the lower triangle of a `side × side` matrix in column-major
  order, `side * (side + 1) ÷ 2` forms: `(1,1), (2,1), ..., (side,1), (2,2), ...`.
  Entries are stored without scaling. For `HPSD`, entry `(i, j)` with `i > j` is
  the stored form and entry `(j, i)` is its complex conjugate.

`label` is an index into the label list that the pipeline returns with the
program (block origin, symmetry sector, strengthening kind). The conic layer
copies it to its reports and does not read it.
"""
struct ConeBlock
    kind::ConeKind
    side::Int
    first_form::Int
    label::Int32
end

"""
    nforms(kind::ConeKind, side::Integer) -> Int

The number of forms that a block of this kind and side uses.
"""
nforms(kind::ConeKind, side::Integer) = kind in (PSD, HPSD) ? (side * (side + 1)) >> 1 : Int(side)
nforms(b::ConeBlock) = nforms(b.kind, b.side)
form_range(b::ConeBlock) = b.first_form:(b.first_form + nforms(b) - 1)

"""
    tri_offset(n, i, j) -> Int

The 0-based position of entry `(i, j)`, `i >= j`, in the stored lower triangle
of a matrix of side `n`. The form index is `block.first_form + tri_offset(n, i, j)`.
"""
tri_offset(n::Integer, i::Integer, j::Integer) = (j - 1) * n - ((j - 1) * (j - 2)) >> 1 + (i - j)

"""
    MomentProgram{C}

A moment relaxation as a conic program:

    minimise    sum(objective_coefs[j] * y[objective_ids[j]])
    subject to  y[1] == 1,
                the forms of each block are in the cone of the block.

The variables `y[2:nvars]` are free real numbers. A complex moment uses two
variables (real part and imaginary part); `MomentLabel` records which.

A program can grow: `push_block!` adds a block after a solve, which the cone
separation loop needs.
"""
mutable struct MomentProgram{C<:Number}
    nvars::Int
    const objective_ids::Vector{Int32}
    const objective_coefs::Vector{Float64}
    const forms::LinearForms{C}
    const blocks::Vector{ConeBlock}
end

MomentProgram{C}(nvars::Integer=1) where {C<:Number} =
    MomentProgram{C}(Int(nvars), Int32[], Float64[], LinearForms{C}(), ConeBlock[])

nvars(p::MomentProgram) = p.nvars

"""
    push_block!(program, kind, side, label = 0) -> ConeBlock

Declare that the last `nforms(kind, side)` closed forms of `program.forms` are
one block. The forms of a block must be pushed and closed first, in the order
that `ConeBlock` gives.
"""
function push_block!(p::MomentProgram, kind::ConeKind, side::Integer, label::Integer=0)
    n = nforms(kind, side)
    first_form = length(p.forms) - n + 1
    last_used = isempty(p.blocks) ? 0 : last(form_range(last(p.blocks)))
    first_form == last_used + 1 || throw(ArgumentError(
        "push_block!: a $kind block of side $side needs $n new forms; $(length(p.forms) - last_used) forms were closed after the last block."))
    block = ConeBlock(kind, Int(side), first_form, Int32(label))
    push!(p.blocks, block)
    return block
end

"""
    MomentLabel{W}(word, part)

What one moment variable means: `part` is `0x00` for the constant `y[1]`, `0x01`
for the real part of the moment of `word`, and `0x02` for its imaginary part.
With a symmetry, `word` is the representative of an orbit. The pipeline returns
a `Vector{MomentLabel{W}}` of length `nvars` with each program; the conic layer
does not read it.
"""
struct MomentLabel{W}
    word::W
    part::UInt8
end

"""
    validate(program::MomentProgram) -> program

Check the invariants of the format and throw an `ArgumentError` that names the
first broken one. The conic layer can assume a validated program.
"""
function validate(p::MomentProgram{C}) where {C}
    n = nvars(p)
    n >= 1 || throw(ArgumentError("A program needs nvars >= 1: y[1] is the constant 1."))
    f = p.forms
    length(f.ids) == length(f.coefs) || throw(ArgumentError("forms.ids and forms.coefs differ in length."))
    first(f.offsets) == 1 && issorted(f.offsets) && last(f.offsets) == length(f.ids) + 1 ||
        throw(ArgumentError("forms.offsets is not a valid row pointer; a form is open or the offsets are not sorted."))
    all(id -> 1 <= id <= n, f.ids) || throw(ArgumentError("A form uses a moment id outside 1:$n."))
    length(p.objective_ids) == length(p.objective_coefs) ||
        throw(ArgumentError("objective_ids and objective_coefs differ in length."))
    all(id -> 1 <= id <= n, p.objective_ids) || throw(ArgumentError("The objective uses a moment id outside 1:$n."))
    next = 1
    for (k, b) in enumerate(p.blocks)
        b.side >= 1 || throw(ArgumentError("Block $k has side $(b.side)."))
        b.first_form == next || throw(ArgumentError("Block $k starts at form $(b.first_form); expected $next. Blocks must use consecutive forms."))
        next += nforms(b)
        b.kind == HPSD && !(C <: Complex) && throw(ArgumentError("Block $k is HPSD but the coefficient type $C is real; use PSD."))
        if C <: Complex && b.kind in (NONNEGATIVE, SOC, PSD)
            all(j -> iszero(imag(f.coefs[j])), f.offsets[b.first_form]:(f.offsets[last(form_range(b)) + 1] - 1)) ||
                throw(ArgumentError("Block $k is $(b.kind) but has a coefficient with an imaginary part."))
        end
        if b.kind == HPSD
            for i in 1:b.side
                r = term_range(f, b.first_form + tri_offset(b.side, i, i))
                all(j -> iszero(imag(f.coefs[j])), r) ||
                    throw(ArgumentError("Block $k is HPSD but diagonal entry $i has a coefficient with an imaginary part."))
            end
        end
    end
    next == length(f) + 1 || throw(ArgumentError("$(length(f) - next + 1) forms at the end belong to no block."))
    return p
end
