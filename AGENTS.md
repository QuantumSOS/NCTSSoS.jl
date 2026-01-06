# AGENTS.md

## Notes
- Dropping `Term` usually recreates it as `PauliMonomial` (scalar × basis word); expect extra allocation unless you keep a 1-term fast path.
- CCR/CAR (boson/fermion): product of two basis words can expand to a sum; keep moment *labels* as `Monomial{A,T}` basis, and represent matrix entries as linear combinations.
- If you keep separate `Monomial`/`PauliMonomial`/`PhysicsMonomial`, you still need mixed ops (`Monomial×PhysicsMonomial`, distributive `*`, `+`, `adjoint`, sort/dedup/zero-drop) to stay closed.
- `Dict` keys: `Vector` uses content `isequal`/`hash` (identity only for `IdDict`), but mutating a vector after insertion breaks the dict; tuple-words give type-instability for variable-length words → keep `Vector{T}` storage and use copy-on-insert / an internal `WordKey{T}` wrapper when you need “frozen” keys.

## Conversation Transcript

`````text
User:
# AGENTS.md instructions for /Users/exaclior/projects/NCTSSoS.jl-review

<INSTRUCTIONS>
# AGENTS.md

At start, say hi + 1 inspire line.
Owner: yushengzhao.
Style: *telegraph; noun-phrases ok; drop filler; min tokens; ship small diffs; verify before handoff.*
Rule: *ask for clarification if unclear*

## Scope
Codex CLI config tracked in git. Path: `home/configs/codex/`.

## Key Files
- `home/configs/codex/config.toml`: main config (model, MCP, trust).
- `home/codex.nix`: links `~/.codex/config.toml` → `home/configs/codex/config.toml` (`mkOutOfStoreSymlink`).

## Agent Protocol
- Safe default: `git status`, `git diff`, `git log`.
- Deletes: use `trash …` (no `rm` unless explicit).
- No destructive git unless explicit (`reset --hard`, `clean`, `restore`, …).
- Prefer root-cause fix; add small regression check when it fits.
- External knowledge: prefer `gm -p "Search {question} online"` (if `gm` exists); else Exa (`mcp__exa__web_search_exa`, `mcp__exa__get_code_context_exa`). Quote exact errors; prefer 2024–2026 sources. Network may require approval.
- “Make a note” => edit this file.

--- project-doc ---

## Skills
These skills are discovered at startup from multiple local sources. Each entry includes a name, description, and file path so you can open the source for full instructions.
- skill-creator: Guide for creating effective skills. This skill should be used when users want to create a new skill (or update an existing skill) that extends Codex's capabilities with specialized knowledge, workflows, or tool integrations. (file: /Users/exaclior/.codex/skills/.system/skill-creator/SKILL.md)
- skill-installer: Install Codex skills into $CODEX_HOME/skills from a curated list or a GitHub repo path. Use when a user asks to list installable skills, install a curated skill, or install a skill from another repo (including private repos). (file: /Users/exaclior/.codex/skills/.system/skill-installer/SKILL.md)
- Discovery: Available skills are listed in project docs and may also appear in a runtime "## Skills" section (name + description + file path). These are the sources of truth; skill bodies live on disk at the listed paths.
- Trigger rules: If the user names a skill (with `$SkillName` or plain text) OR the task clearly matches a skill's description, you must use that skill for that turn. Multiple mentions mean use them all. Do not carry skills across turns unless re-mentioned.
- Missing/blocked: If a named skill isn't in the list or the path can't be read, say so briefly and continue with the best fallback.
- How to use a skill (progressive disclosure):
  1) After deciding to use a skill, open its `SKILL.md`. Read only enough to follow the workflow.
  2) If `SKILL.md` points to extra folders such as `references/`, load only the specific files needed for the request; don't bulk-load everything.
  3) If `scripts/` exist, prefer running or patching them instead of retyping large code blocks.
  4) If `assets/` or templates exist, reuse them instead of recreating from scratch.
- Description as trigger: The YAML `description` in `SKILL.md` is the primary trigger signal; rely on it to decide applicability. If unsure, ask a brief clarification before proceeding.
- Coordination and sequencing:
  - If multiple skills apply, choose the minimal set that covers the request and state the order you'll use them.
  - Announce which skill(s) you're using and why (one short line). If you skip an obvious skill, say why.
- Context hygiene:
  - Keep context small: summarize long sections instead of pasting them; only load extra files when needed.
  - Avoid deeply nested references; prefer one-hop files explicitly linked from `SKILL.md`.
  - When variants exist (frameworks, providers, domains), pick only the relevant reference file(s) and note that choice.
- Safety and fallback: If a skill can't be applied cleanly (missing files, unclear instructions), state the issue, pick the next-best approach, and continue.
</INSTRUCTIONS>

User:
<environment_context>
  <cwd>/Users/exaclior/projects/NCTSSoS.jl-review</cwd>
  <approval_policy>on-request</approval_policy>
  <sandbox_mode>workspace-write</sandbox_mode>
  <network_access>restricted</network_access>
  <shell>zsh</shell>
</environment_context>

User:
A lot of trouble are caused due to Monomials not in their normal form. What if I force `Monomial` to be in their normal form by construction? This would guard against a huge range of problems, like being able to compare two Monomials not in their normal form. But the problem is I need to support Fermionic and Bosonic Monomial. A monomial's normal form could be a polynomial. What will be the most elegant way to put this?

Assistant:
Hi.
Small invariants, big leverage.

- Key obstruction: boson/fermion “normal ordering” is `Monomial(word) -> Polynomial(...)` (δ / contraction terms). Julia constructors can’t change the return type, so “`Monomial` always normal by construction” can’t hold globally.
- Clean split: `Monomial` = *basis word only* (guaranteed already normal-ordered); `NormalForm` = `Polynomial` (sparse linear combination) = the canonical algebra element.
- Unify API: `nf(m::Monomial{A})::Polynomial{A}` that wraps `simplify(m)` (if `simplify` returns `Monomial`/`Term`, lift to 1-term `Polynomial`; if it returns `Polynomial`, use it). Then make user-facing arithmetic live on `Polynomial`/`NormalForm` so you never compare raw words.
- Extra guard: for `BosonicAlgebra`/`FermionicAlgebra`, make `Monomial{A}(word)` *reject* non-normal-ordered words (“use `nf`/`simplify`”), so `==/hash/isless` stay meaningful.

Clarify: do you want `*(m1::Monomial, m2::Monomial)` to start returning a `Polynomial` (breaking but safest), or keep `*` as raw concatenation and introduce a “safe” multiply/constructor (`nf`, `mul_nf`, etc.)?

User:
What I meant to say is: given a vector of integers each representing a non-commuting variable hencing representing a monomial in its entierity, the constructor `Monomial` should find the normal form of the monomial (i.e after physical site based commutation and simplification via unipotency etc) and construct such a Monomial. The only concern I have is: normal form of bosonic operator string is a weighted sum of bosnonic operator string (the normal order ones) is there an elegant, julian way to represent Monomial of Bosonic Operator String? If unclear, keep asking me questions.

Assistant:
Hi.
Invariant-first; bugs nowhere to hide.

- Bosonic case: “normal form” of a word is **not a monomial**; it’s a weighted sum of normal‑ordered words ⇒ represent it as `Polynomial{BosonicAlgebra}` (your current `Polynomial` is already the right container: coeff + `Monomial{BosonicAlgebra}` terms).
- Make `Monomial{BosonicAlgebra}` mean **basis word only** (already normal‑ordered). Enforce in constructor: if the input word is not normal‑ordered, throw and tell the user to call `simplify`/`normal_form`.
- Add a single entrypoint like `normal_form(::Type{A}, word)::Polynomial{A}` that does “commute-by-site + algebra simplify”; for algebras where `simplify(::Monomial)` returns `Monomial`/`Term`, just lift to a 1‑term `Polynomial`.

Clarify: do you want `Monomial{PauliAlgebra}(word)` to also absorb the ±i phase (so it can’t return a `Monomial` either), or is “normal form” for you only the commuting/site ordering part?

User:
How about this "struct Monomial{A,T} <: AbstractMonomial
    normal_word::Vector{T}
end

struct PauliMonomial{T} <: AbstractMonomial
    phase::ComplexF64
    normal_word::Vector{T}
end

# Suggest a better name, more mathematical
struct PhysicsMonomial{A,T} <: Abstractmonomial
    weight::Vector{Float64}
    normal_words::Vector{Vector{T}}
end
" where A is algebras and T are index types for variable

Assistant:
Hi.
Name things by algebra, not by hope.

- Your `PauliMonomial{T}(phase, normal_word)` is exactly `Term{Monomial{PauliAlgebra,T},ComplexF64}` (scalar × basis word). I’d reuse `Term`, not add a new “monomial” type.
- Your `PhysicsMonomial{A,T}(weights, normal_words)` is exactly a `Polynomial{A,T,C}` (linear combination of basis words). More mathematical name: `NormalForm` or `LinearCombination`; but really it *is* a polynomial.
- Keep `Monomial{A,T}` = **canonical basis word only** (normal-ordered word, no coefficients, no sums). Enforce in constructor/check.
- Provide one “from raw word” API that always works: `normal_form(::Type{A}, word)::Polynomial{A}` (even for Pauli; returns 1-term poly if no expansion).

Clarify: do you want to *forbid* constructing non-normal-ordered `Monomial{BosonicAlgebra}` (throw), and force users through `normal_form/simplify` for boson/fermion?

User:
but lacking PauliMonomial and PhysicsMonomial lacks the ability to clearly mark when is an object a member in the basis and when is it really a polynomial

Assistant:
Hi.
Clear types; clear algebra.

- You already have the “basis vs polynomial” markers:
  - basis element (single word, coef = 1): `Monomial{A,T}`
  - “PauliMonomial” (phase × one word): `Term{Monomial{PauliAlgebra,T},ComplexF64}`
  - “PhysicsMonomial” (Σ weightᵢ * wordᵢ): `Polynomial{A,T,C}` (this is the mathematical object; don’t name it “monomial”)

- Most Julian way to get your naming without duplicating structs: type aliases + API names:
  ```julia
  const BasisWord{A,T} = Monomial{A,T}
  const ScaledWord{A,T,C} = Term{Monomial{A,T},C}
  const NormalForm{A,T,C} = Polynomial{A,T,C}
  ```
  Then expose `normal_form(::Monomial)::NormalForm` (always a polynomial; 1-term for most algebras).

- If you also want “raw word” vs “already-normal” tracked by type (not just by convention), you need a *wrapper* (`RawWord`/`NormalWord`) or a separate `RawMonomial` type; adding `PauliMonomial/PhysicsMonomial` just reimplements `Term/Polynomial`.

Question: do you need the distinction for dispatch/safety (compiler-enforced), or just for readability in the public API/docs?

User:
I need it for for dispatch/safety (compiler-enforced). If I can guarantee each monomial is always in normal form, then comparision defined on monomials will be less error prone

Assistant:
Hi.
Strong invariants; weak bugs.

- You already have the 3 compiler-enforced “buckets” you want:
  - **basis element**: `Monomial{A,T}` (one normal-ordered word)
  - **scaled basis element** (Pauli phase etc): `Term{Monomial{A,T},C}`
  - **sum** (boson/fermion normal form): `Polynomial{A,T,C}`
  - So `PauliMonomial` ≅ `Term{Monomial{PauliAlgebra,T},ComplexF64}`, and “PhysicsMonomial” ≅ `Polynomial`.

- The missing piece is **preventing non-normal words from ever becoming `Monomial`**. For that you need a separate “raw” carrier type, e.g. `RawWord{A,T}` (or a `NormalizationTag` type parameter).
  - `Monomial{A,T}`: only constructible if `word` already satisfies the basis invariants (site-sorted / normal-ordered / idempotency applied / etc). Otherwise throw.
  - `normal_form(::Type{A}, word)::Polynomial{A}`: the *smart constructor* that rewrites a raw word into the canonical algebra element:
    - returns `Monomial` for algebras with 1-word normal form (projector/unipotent/…)
    - returns `Term` for Pauli
    - returns `Polynomial` for boson/fermion

- Then comparisons stay safe by construction: `==/hash/isless` defined on `Monomial` assumes “already normal”; raw words don’t get those methods (or they error).

Two clarifiers:
1) Do you want `*(::Monomial, ::Monomial)` to **return canonical** (i.e. `Term/Polynomial` depending on algebra), or keep `*` fast and force an explicit `normal_form` step?
2) For fermions: should “nilpotent words” be **unconstructible as `Monomial`** (since they’re zero), or allowed but treated specially?

User:
Hold on, what about this design: "struct Monomial{A,T} <: AbstractMonomial
    # always guaranteed to be normal
    # 1. For non-commutative meaning site sorted
    # 2. For unipotent, projective algebra, in addition to site sorted, meaning simplified to no-consecutive
    # 3. For Pauli algebra, in addition to site sorted, meaning at most 1 pauli operator for each site
    # 4. For fermionic and bosonic algebra: meaning all creation operators on the left of anhilation operators and creation/anhilation operators are sorted according to physical site
    normal_word::Vector{T}
end

struct ComposedMonomial{AS} <:AbstractMonomial
    normal_monos::Tuple{AS} # not sure if this is correct
end

# Type needs to be more careful, perhaps only 2 MT <: AbstractMonomial , could be either Monomial{A,T} or ComposedMonomial{AS}
struct Polynomial{A,T,F} <: AbstractPolynomial
    # We lose the notion of Term here
    coeffs::Vector{F}
    monos::Vector{Monomial{A,T}} # Vector{ComposedMonomial{AS}} in general
end

struct PauliMonomial{T}
    phase::ComplexF64
    normal_word::Vector{Monomial{PauliAlgebra,T}}
end

# Suggest a better name, more mathematical
struct PhysicsMonomial{A,T} # here A can only be either FermionicAlgebra or BosonicAlgebra
    weight::Vector{Float64}
    normal_words::Vector{Monomial{A,T}}
end
"

Assistant:
Hi.
Type-invariants; fearless refactors.

- Rejecting non-basis words for `Pauli/Fermion/Boson` is consistent: their “normal form” lives in coefficients/sums, not in `Monomial`.
- Consequence: `*(::Monomial{Pauli/Boson/Fermion}, ::Monomial{...})` must return a **element** type (your `Polynomial`), not a `Monomial`, otherwise you can’t even build `σx₁*σy₁` or `c₁*c₁†`.
- If you remove `Term`, then `PauliMonomial` is just “a Term for Pauli”; and `PhysicsMonomial` is just `Polynomial` restricted to `{Bosonic,Fermionic}`. You’ll either (a) allocate 1-term `Polynomial`s everywhere, or (b) reinvent `Term` anyway.
- Constructor policy:
  - `Monomial{NonCommutative/Projector/Unipotent}(word)` can canonicalize (site-sort + local rule) and store.
  - `Monomial{Pauli}(word)` should *check* “already reduced basis” (site-sorted + ≤1 op per site) else throw.
  - `Monomial{Boson/Fermion}(word)` should *check* `is_normal_ordered(word)` else throw.

Two questions:
1) OK with algebra-dependent `*` return type (monomial for some, polynomial for others), or do you want `*` to always return `Polynomial` for uniformity?
2) “Make a note”: where should I write it (`CLAUDE.md`, or should I add a repo-root `AGENTS.md`)?

User:
Always normal could be enforced by simplification process in the constructor for algebra that simplies from monomial to monomial. I am willing to reject inputs that are not normal from pauli, fermi, bosonic algebra. I will remove Term and simplify Polynomial accordingly. I will remove Term, I like my chances of keeping them invariant (but do make a note about this). There's no need for raw monomial, its just a thin wrapper around vector of integers. There is no need for normal form, constructor forces normal form.

Assistant:
Hi.
Type-invariants; fearless refactors.

- If `m1*m2` must always return a “monomial-ish” object, then for boson/fermion that object **must be a linear combination** (your `PhysicsMonomial`), because δ/commutator terms are real.
- Then you’re forced into full distributive algebra anyway (`+(…), *(…), adjoint, simplify, combine-like-terms, sort/dedup/zero-drop`) on `PhysicsMonomial`, plus mixed ops (`Monomial * PhysicsMonomial`, etc.). Easiest: make **one canonical sum type** (call it `NormalForm`/`NCElement`/`LinearCombination`) and treat `Monomial` and `PauliMonomial` as 0/1-term special cases (or wrappers) that **promote** to it for `+`/general `*`.
- Warning: removing `Term` means you will reimplement it as `PauliMonomial` (1-term fast path). If you skip the fast path, you’ll allocate 1-term vectors everywhere.

Clarify (core): when you say “entries of moment matrices are monomials”, do you mean the **moment labels** `u†v` (single basis word) or the **symbolic matrix entries** after applying CCR/CAR (often sums)?

User:
1. yes 2. Why would I need to support Monomial * PhysicsMonomial? Elaborate more on sum type 3. What do you mean by fast path? I don't expect to store PauliMonomial in Polynomial. In Polynomial, we will only store Monomial  4. I meant the moment labels u†v

Assistant:
Hi.
Clean invariants; scalable algebra.

- Moment labels `u†v`: for CCR/CAR, `u†v` in the **quotient algebra** can equal a *linear combination* of normal‑ordered words (e.g. `c c† = c† c + 1`). So the *moment variables* should still be keyed by **basis** `Monomial{A,T}`; the matrix entry is then a linear combo of those basis moments.
- `Monomial * PhysicsMonomial`: closure pressure. Localizing entries look like `u† * (Σ cₖ wₖ) * v`. Even if you start with `Monomial*Monomial`, you quickly get a sum, then you still need `sum * Monomial` (or you must special-case “3-arg multiply then simplify” everywhere).
- “Sum type”: one `AbstractElement` with concrete variants `{basis word, scaled word, linear-combo}` + promotion rules for mixed ops. That’s exactly the role of `Term` (scaled word) + `Polynomial` (linear combo) in your current code.
- “Fast path”: avoid allocating length‑1 `coeffs/monos` vectors. Your `PauliMonomial` is that fast path (i.e. it *is* `Term` under a new name).

Clarifier: do you want the JuMP monomap/dictionary keys to stay `Monomial{A,T}` only (recommended), or do you intend to create independent moment variables for `PhysicsMonomial` sums too?

User:
in summary, do you advice for or against such rewrite?

Assistant:
Hi.
Small diff; big invariant win.

- Against full rewrite (drop `Term`, add `PauliMonomial`/`PhysicsMonomial`): you’ll reimplement `Term`/`Polynomial` anyway, plus N×M mixed-method surface (`Monomial×PhysicsMonomial`, `+`, `adjoint`, hashing, sorting, combine-like-terms, zero-drop).
- For incremental change: keep `Monomial{A,T}` = **basis word, always canonical**, and enforce that in constructors for algebras where “normalize” stays 1 word (NC/Projector/Unipotent).
- For Pauli/boson/fermion: `Monomial{A}(word)` should be **checked-only** (reject non-basis words); algebraic multiplication should return an **element** type (your current `Term` for Pauli, `Polynomial` for boson/fermion).
- Moment labels `u†v`: keep as `Monomial{A,T}` keys; matrix entries become linear combos of those labels (already what your `simplify` pipeline models).
- If you really want “monomial-ish closure”, reuse existing types but rename via aliases (basis/scaled/sum). Don’t delete `Term`; it’s the needed 1-term fast path.

Where should I write the “note” you asked for earlier (`CLAUDE.md` vs a new `AGENTS.md` at repo root)?

User:
to a new `AGENTS.md` at repo root

Assistant:
Hi.
Notes now; refactors later.

- Added repo-root note in `AGENTS.md:1` re: dropping `Term` (expect reintroducing 1-term fast path + mixed-op closure pressure).
- FYI: untracked `imagined.jl:1` appeared (your sketch); tell me if you want it kept/ignored/trashed.

User:
Export our whole conversation history into `AGENTS.md`
````

Assistant:
Combining like terms in fermionic/bosonic normal-ordering now relies on dictionaries keyed by mutable vectors, so identical monomials are not merged and coefficients are wrong in multiply/adjoint/simplify paths. This correctness issue affects core algebra operations.

Full review comments:

- [P1] Duplicate physics monomials not merged after multiplication — src/types/physics_monomial.jl:292-299
  In `PhysicsMonomial` multiplication the like-term aggregation uses `Dict{Vector{T},Int}` keyed by raw word vectors, but arrays are hashed/compared by object identity, not contents, so identical words produced from different contraction paths end up as separate entries instead of being summed. This yields wrong coefficients whenever simplification generates duplicate words (e.g., distributing `(a₁+a₂)*(a₁+a₂†)` can leave multiple identical normal-ordered words unmerged). Use an immutable, value-compared key such as the `Monomial` itself or `Tuple(word)` to combine terms; the same pattern appears in the adjoint aggregation block nearby.

- [P1] Fermionic/bosonic simplify fails to combine identical terms — src/types/physics_monomial.jl:569-577
  `_combine_physics_terms` groups results with a `Dict{Vector{T},Int}` keyed by mutable vectors, which are compared by identity rather than contents. If normal-ordering produces the same canonical word more than once (common for fermionic Wick expansions), the duplicates are kept as separate terms instead of summing coefficients, giving incorrect PhysicsMonomial outputs. Switching the key to a content-equal type (e.g., `Monomial{A,T}` or `Tuple{T,...}`) is needed to correctly merge like terms across simplify calls.

User:
are you sure hashing are different? I had : "ulia> v1 = [1,2,3]
3-element Vector{Int64}:
 1
 2
 3

julia> v2 = [1,2,3]
3-element Vector{Int64}:
 1
 2
 3

julia> v1h = hash(v1)
0x902ea5284b019f89

julia> v2h = hash(v2)
0x902ea5284b019f89"

Assistant:
Hi.
Measure twice; patch once.

`Dict` in Julia uses `isequal` + `hash`; for `Vector`, both are content-based, so `Dict{Vector,...}` will merge `v1=[1,2,3]` and `v2=[1,2,3]` (your `hash` check is consistent). Identity-based behavior is `IdDict`, not `Dict`.

Real footgun remains: vector keys are mutable—if any word vector is mutated after it’s used as a key (or reused/aliased scratch vectors get mutated later), the dict can “lose” entries and you can get wrong combining.

Do you ever mutate `m.word` / returned `w` after inserting, or reuse the same `Vector` object across terms?

User:
I don't mutate `Monomial`'s word and I have thought about making them a tuple, but that would make `Monomial` type unstable right?

Assistant:
Hi.
Keep values fixed; proofs easy.

- `Dict` is content-equality (`isequal` + `hash`); your `hash(v1)==hash(v2)` is expected for equal vectors. Identity-equality is `IdDict`.
- Tuple storage: variable-length tuple ⇒ length in the type. Either (a) param by length → many `Monomial` concrete types / type-unstable containers, or (b) store as `Tuple`/`Tuple{Vararg{T}}` → abstract field. So for variable-length words, `Vector{T}` is the right representation.
- If you only wanted “immutable key” safety: keep `Vector{T}` words, but key combining dicts by `Monomial{A,T}` (you already define `==/hash`) or `copy(w)` on insert (costly but bulletproof).

`````
