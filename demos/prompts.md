# Demo Prompts for polyopt-guide

Copy-paste one of these into Claude Code to start a guided session.

---

## Prompt 1: Polyopt expert, NCTSSoS novice

```
I want to compute the maximum quantum violation of the I₃₃₂₂ Bell inequality using NCTSSoS. I know how to set this up mathematically — it's nc polynomial optimization over dichotomic observables A₁,A₂,A₃ (Alice) and B₁,B₂ (Bob) with Aᵢ²=Bⱼ²=I and [Aᵢ,Bⱼ]=0. The Bell functional is I = A₁(B₁+B₂) + A₂(B₁+B₂) + A₃(B₁-B₂) - A₁ - 2B₁ - B₂. I've done this before with ncpol2sdpa at NPA level 2 — how do I translate that to NCTSSoS? I also want to compare bounds with and without sparsity exploitation.
```

---

## Prompt 2: Physicist, no optimization background

```
I'm studying a 4-site Fermi-Hubbard chain with periodic boundaries, t=1, U=4 at half-filling. The Hamiltonian is H = -t Σᵢ(c†ᵢσ cᵢ₊₁σ + h.c.) + U Σᵢ nᵢ↑ nᵢ↓. I know from DMRG the half-filled ground state energy is around -2.83, but I also know the bare Hamiltonian by itself is grand-canonical unless I add particle-number constraints. How do I set this up in NCTSSoS so the bound is really for half-filling, and how would I extract the nearest-neighbor spin correlation ⟨Sᶻ₁Sᶻ₂⟩ from the solved moments?
```
