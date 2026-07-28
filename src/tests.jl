## ============================================================
## tests.jl — Statistical axis of the DTV framework
## ============================================================
##
## Implements the per-measurement fidelity gap and the one-sided
## paired-t equivalence test of `3_framework.tex` §III
## ("Verification as a Likelihood Equivalence Test"), plus the
## three-region joint verdict.
##
## All functions take only J* arrays (or scalars) plus the
## measurement-set size. They do not depend on PMDSE / KPIUtils.
##
## Sign convention. For any likelihood-admissible M-estimator the
## profile log-likelihood satisfies
##
##     ln L_M = -J*(M) + C_ρ
##
## (§III, "State Estimation as the Likelihood Generator"). The
## estimator-specific constant C_ρ is identical for the reference and
## the candidate and cancels in the paired difference, so the fidelity
## gap below is well-defined for WLS, WLAV, rWLAV, etc.
## ============================================================

# Declared verdict knobs (deployment-dependent; the values below are the
# defaults used in the paper's results, 4_results.tex §IV).
const DELTA_EQ_DEFAULT = 0.05   # equivalence margin δ_eq (loss units per measurement)
const ALPHA_EQ_DEFAULT = 0.05   # one-sided significance α_eq of the equivalence test
const ALPHA_OP_DEFAULT = 0.10   # operational tolerance α on the worst-alarm disagreement η_∞

# ── Helpers ────────────────────────────────────────────────────────────────

_as_vector(x::Real) = [Float64(x)]
_as_vector(x::AbstractVector) = Float64.(collect(x))

# ── Fidelity gap d̄ ──────────────────────────────────────────────────────────

"""
    fidelity_gap(J_r, J_c, M) → (d_bar, s_d, d_j)

Per-measurement excess negative log-likelihood of the candidate relative to
the reference (`3_framework.tex` eq. for `d_j` and `d̄`):

    d_j = (J*_c(j) − J*_r(j)) / |M_j|,     d̄ = mean_j d_j,     s_d = std_j d_j.

`J_r`, `J_c` are paired per-realisation objectives (scalars for a single
realisation, or equal-length vectors for `n` Monte-Carlo realisations). `M` is
the measurement-set size |M|: a scalar applied to every realisation, or a
per-realisation vector. `s_d` is the sample standard deviation (undefined and
returned as `NaN` for `n < 2`). Non-positive `M` entries yield `NaN` gaps.
"""
function fidelity_gap(J_r, J_c, M)
    Jr = _as_vector(J_r)
    Jc = _as_vector(J_c)
    length(Jr) == length(Jc) || error(
        "DTVerify.fidelity_gap: J_r and J_c must have equal length " *
        "(got $(length(Jr)) vs $(length(Jc))).")
    n = length(Jr)
    Mv = M isa Real ? fill(Float64(M), n) : Float64.(collect(M))
    length(Mv) == n || error(
        "DTVerify.fidelity_gap: M must be a scalar or a length-$n vector " *
        "(got length $(length(Mv))).")
    d = similar(Jr)
    @inbounds for j in 1:n
        d[j] = Mv[j] > 0 ? (Jc[j] - Jr[j]) / Mv[j] : NaN
    end
    d_bar = mean(d)
    s_d = n < 2 ? NaN : std(d)   # sample std (denominator n−1)
    return (d_bar = d_bar, s_d = s_d, d_j = d)
end

# ── Equivalence test (one-sided paired t) ───────────────────────────────────

"""
    equivalence_test(d_bar, s_d, n; δ_eq=$(DELTA_EQ_DEFAULT), α_eq=$(ALPHA_EQ_DEFAULT))
        → (holds, bound, t_crit, δ_eq, α_eq)

One-sided, paired equivalence test on the fidelity gap (`3_framework.tex`
eq. for the equivalence confidence bound). Equivalence is accepted when the
upper one-sided `(1−α_eq)` confidence bound lies below the margin:

    d̄ + t_{1−α_eq, n−1} · s_d/√n  <  δ_eq.

`holds` is `true` when the candidate is statistically equivalent to the
reference at confidence `1−α_eq`. Requires `n ≥ 2` with finite `s_d`;
otherwise the variance is undefined, the null cannot be rejected, and
`(holds=false, bound=NaN, …)` is returned with a `@warn`.
"""
function equivalence_test(d_bar, s_d, n::Integer;
                          δ_eq::Real = DELTA_EQ_DEFAULT,
                          α_eq::Real = ALPHA_EQ_DEFAULT)
    if n < 2 || !isfinite(s_d) || !isfinite(d_bar)
        @warn "equivalence_test: need n ≥ 2 with finite d̄, s_d (got n=$n, " *
              "d̄=$d_bar, s_d=$s_d); cannot reject H₀."
        return (holds = false, bound = NaN, t_crit = NaN, δ_eq = δ_eq, α_eq = α_eq)
    end
    t_crit = quantile(TDist(n - 1), 1 - α_eq)
    bound = d_bar + t_crit * s_d / sqrt(n)
    return (holds = bound < δ_eq, bound = bound, t_crit = t_crit,
            δ_eq = δ_eq, α_eq = α_eq)
end

# ── Three-region joint verdict ──────────────────────────────────────────────

"""
    joint_verdict(equivalence_holds, η_inf; α=$(ALPHA_OP_DEFAULT)) → Symbol

Two-axis verdict of `3_framework.tex` §III ("Joint Verdict"):

  • `:consequential`  when `η_inf > α`                       (decision-changing),
  • `:verified`       when `η_inf ≤ α` and equivalence holds (statistically equivalent),
  • `:tolerable`      when `η_inf ≤ α` and equivalence fails (detectable, decision-preserving).

`equivalence_holds` is the boolean from [`equivalence_test`](@ref); `η_inf` is
the worst-alarm operational disagreement `max_D (1 − ACC_D)`; `α` is the
operational tolerance.
"""
function joint_verdict(equivalence_holds::Bool, η_inf::Real;
                       α::Real = ALPHA_OP_DEFAULT)
    (isfinite(η_inf) && η_inf > α) && return :consequential
    return equivalence_holds ? :verified : :tolerable
end
