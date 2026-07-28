## ============================================================
## api.jl — Main user-facing entry points
## ============================================================
##
## Two functions, one NamedTuple schema:
##   verify(se_r, math_r, se_c, math_c; …)              — single realisation
##   verify_mc(se_r_vec, math_r, se_c_vec, math_c; …)   — n realisations
##
## Both return the framework's two-axis result: the statistical fidelity
## gap d̄ with its equivalence test, the operational worst-alarm
## disagreement η_∞, the diagnostic MAPE_U, and the joint verdict
## (:verified / :tolerable / :consequential) of `3_framework.tex` §III.
## ============================================================

# Inline solution-dict extractor: PMDSE wraps the bus dict either at the
# top level (`sol["bus"]`) or under `sol["solution"]["bus"]`.
function _sol(se)
    haskey(se, "bus") && return se
    haskey(se, "solution") && return se["solution"]
    return se
end

# Disagreement counts (FP + FN) per decision from a pooled `raw` NamedTuple.
function _counts(raw)
    m(c) = c.fp + c.fn
    vv, nev, vuf, pvur = m(raw.vv), m(raw.nev), m(raw.vuf), m(raw.pvur)
    return (vv = vv, nev = nev, vuf = vuf, pvur = pvur,
            total = vv + nev + vuf + pvur)
end

# Assemble the shared result NamedTuple from the two axes.
function _result(; verdict, fg, eq, op, raw, mape_u, maxape_u,
                   n, n_meas, kr, kc, δ_eq, α_eq, α, Jr, Jc)
    return (
        verdict           = verdict,
        d_bar             = fg.d_bar,
        s_d               = fg.s_d,
        equivalence_holds = eq.holds,
        eq_bound          = eq.bound,
        eta_inf           = op.eta_inf,
        eta               = op.eta,
        acc               = op.acc,
        rec               = op.rec,
        mape_u            = mape_u,
        maxape_u          = maxape_u,
        counts            = _counts(raw),
        raw               = raw,
        n_realizations    = n,
        n_meas            = n_meas,
        k_r               = kr,
        k_c               = kc,
        δ_eq              = δ_eq,
        α_eq              = α_eq,
        α                 = α,
        J_star_r          = Jr,
        J_star_c          = Jc,
    )
end

"""
    verify(se_r, math_r, se_c, math_c;
           δ_eq=$(DELTA_EQ_DEFAULT), α_eq=$(ALPHA_EQ_DEFAULT), α=$(ALPHA_OP_DEFAULT),
           k_r=nothing, k_c=nothing, M=nothing, U_nom=NaN,
           vv_lo=0.9, vv_hi=1.1, nev_thresh=10.0, vuf_thresh=2.0, pvur_thresh=2.0)
        → NamedTuple

Single-realisation digital-twin verification of candidate `M_c` against
reference `M_r`. Computes the per-measurement fidelity gap `d̄`, the
operational worst-alarm disagreement `η_∞`, `MAPE_U`, and the joint verdict.

A single realisation carries no variance, so the equivalence test cannot be
run (`equivalence_holds = false`) and the verdict is therefore never
`:verified` — use [`verify_mc`](@ref) with `n ≥ 2` for a full statistical
verdict. The measurement count `M = |M|` defaults to `n_measurements(math_r)`;
pass `M` explicitly if the math dict carries no `meas` block.
"""
function verify(
    se_r::AbstractDict, math_r::AbstractDict,
    se_c::AbstractDict, math_c::AbstractDict;
    δ_eq::Real = DELTA_EQ_DEFAULT,
    α_eq::Real = ALPHA_EQ_DEFAULT,
    α::Real = ALPHA_OP_DEFAULT,
    k_r::Union{Integer, Nothing} = nothing,
    k_c::Union{Integer, Nothing} = nothing,
    M::Union{Real, Nothing} = nothing,
    U_nom::Float64 = NaN,
    vv_lo::Float64 = 0.9,
    vv_hi::Float64 = 1.1,
    nev_thresh::Float64 = 10.0,
    vuf_thresh::Float64 = 2.0,
    pvur_thresh::Float64 = 2.0,
)
    kr = isnothing(k_r) ? n_state_variables(math_r) : Int(k_r)
    kc = isnothing(k_c) ? n_state_variables(math_c) : Int(k_c)

    Jr = J_star(se_r)
    Jc = J_star(se_c)

    Mval = isnothing(M) ? n_measurements(math_r) : M
    (Mval isa Real && Mval > 0) ||
        @warn "verify: measurement count |M| = $Mval; the fidelity gap will be " *
              "NaN. Pass `M = <n_measurements>` explicitly."

    fg = fidelity_gap(Jr, Jc, Mval)
    # n = 1 ⇒ variance undefined; equivalence is not testable (no @warn spam).
    eq = (holds = false, bound = NaN, t_crit = NaN, δ_eq = δ_eq, α_eq = α_eq)

    d = decision_disagreement(
        _sol(se_r), math_r, _sol(se_c), math_c;
        U_nom=U_nom, vv_lo=vv_lo, vv_hi=vv_hi,
        nev_thresh=nev_thresh, vuf_thresh=vuf_thresh, pvur_thresh=pvur_thresh,
    )
    op = operational_axis(d.raw)
    verdict = joint_verdict(eq.holds, op.eta_inf; α=α)

    return _result(; verdict=verdict, fg=fg, eq=eq, op=op, raw=d.raw,
                     mape_u=d.mape_u, maxape_u=d.maxape_u,
                     n=1, n_meas=Mval, kr=kr, kc=kc,
                     δ_eq=δ_eq, α_eq=α_eq, α=α, Jr=[Jr], Jc=[Jc])
end

"""
    verify_mc(se_r_vec, math_r, se_c_vec, math_c;
              δ_eq=$(DELTA_EQ_DEFAULT), α_eq=$(ALPHA_EQ_DEFAULT), α=$(ALPHA_OP_DEFAULT),
              k_r=nothing, k_c=nothing, M=nothing, U_nom=NaN,
              vv_lo=0.9, vv_hi=1.1, nev_thresh=10.0, vuf_thresh=2.0, pvur_thresh=2.0)
        → NamedTuple

Monte-Carlo digital-twin verification. `se_r_vec[j]` and `se_c_vec[j]` are the
paired reference / candidate DSSE result dicts for realisation `j`; both math
dicts are held fixed across realisations.

The statistical axis forms the per-measurement gap `d_j`, its mean `d̄` and
spread `s_d`, and the one-sided equivalence test at `(δ_eq, α_eq)`. The
operational axis pools the alarm confusion counts over every item and
realisation, then rates them into `η_D = 1 − ACC_D` and `η_∞ = max_D η_D`.
`MAPE_U` is aggregated as the mean of the per-realisation values (worst-bus
`MaxAPE_U` as the max). The joint verdict combines the two axes.

`M = |M|` defaults to `n_measurements(math_r)`; pass it explicitly when the
math dict has no `meas` block.
"""
function verify_mc(
    se_r_vec::AbstractVector, math_r::AbstractDict,
    se_c_vec::AbstractVector, math_c::AbstractDict;
    δ_eq::Real = DELTA_EQ_DEFAULT,
    α_eq::Real = ALPHA_EQ_DEFAULT,
    α::Real = ALPHA_OP_DEFAULT,
    k_r::Union{Integer, Nothing} = nothing,
    k_c::Union{Integer, Nothing} = nothing,
    M::Union{Real, Nothing} = nothing,
    U_nom::Float64 = NaN,
    vv_lo::Float64 = 0.9,
    vv_hi::Float64 = 1.1,
    nev_thresh::Float64 = 10.0,
    vuf_thresh::Float64 = 2.0,
    pvur_thresh::Float64 = 2.0,
)
    n = length(se_r_vec)
    if n != length(se_c_vec)
        error("DTVerify.verify_mc: se_r_vec ($n) and se_c_vec " *
              "($(length(se_c_vec))) must be the same length.")
    end
    n >= 1 || error("DTVerify.verify_mc: empty input vectors.")

    kr = isnothing(k_r) ? n_state_variables(math_r) : Int(k_r)
    kc = isnothing(k_c) ? n_state_variables(math_c) : Int(k_c)

    Jr = [J_star(se) for se in se_r_vec]
    Jc = [J_star(se) for se in se_c_vec]

    Mval = isnothing(M) ? n_measurements(math_r) : M
    (Mval isa Real && Mval > 0) ||
        @warn "verify_mc: measurement count |M| = $Mval; the fidelity gap will " *
              "be NaN. Pass `M = <n_measurements>` explicitly."

    fg = fidelity_gap(Jr, Jc, Mval)
    eq = equivalence_test(fg.d_bar, fg.s_d, n; δ_eq=δ_eq, α_eq=α_eq)

    # Operational axis: pool per-realisation confusion counts, then rate.
    raws = Vector{Any}(undef, n)
    mape_sum = 0.0; mape_n = 0; maxape_u = NaN
    for j in 1:n
        d = decision_disagreement(
            _sol(se_r_vec[j]), math_r, _sol(se_c_vec[j]), math_c;
            U_nom=U_nom, vv_lo=vv_lo, vv_hi=vv_hi,
            nev_thresh=nev_thresh, vuf_thresh=vuf_thresh, pvur_thresh=pvur_thresh,
        )
        raws[j] = d.raw
        if !isnan(d.mape_u)
            mape_sum += d.mape_u; mape_n += 1
        end
        if !isnan(d.maxape_u)
            maxape_u = isnan(maxape_u) ? d.maxape_u : max(maxape_u, d.maxape_u)
        end
    end
    raw = pool_confusion(raws)
    op = operational_axis(raw)
    mape_u = mape_n == 0 ? NaN : mape_sum / mape_n

    verdict = joint_verdict(eq.holds, op.eta_inf; α=α)

    return _result(; verdict=verdict, fg=fg, eq=eq, op=op, raw=raw,
                     mape_u=mape_u, maxape_u=maxape_u,
                     n=n, n_meas=Mval, kr=kr, kc=kc,
                     δ_eq=δ_eq, α_eq=α_eq, α=α, Jr=Jr, Jc=Jc)
end
