## ============================================================
## decisions.jl — Operational decision evaluators
## ============================================================
##
## Voltage-violation (VV), neutral-to-earth-voltage (NEV), voltage
## unbalance factor (VUF), and phase voltage unbalance rate (PVUR)
## alarms — adapted from `KPIUtils.jl` (lines 123–312, 611–667),
## with the Pliers / PMDUtils dependencies removed.
##
## References:
##   • EN 50160:2010 — Voltage characteristics at LV
##   • IEEE Std 1100-2005 — touch-voltage limit 10 V rms
##   • IEEE Std 141-1993 — PVUR magnitude-deviation definition
##   • Kersting (2012) — symmetrical components
## ============================================================

# ── Symmetrical-component transform ────────────────────────────────────────

"""
    sequence_components(Va, Vb, Vc) → (V0, V1, V2)

Zero-, positive-, and negative-sequence phasors from three phase phasors,
via `α = exp(j·2π/3)` (Kersting 2012, §3.3).
"""
function sequence_components(Va::ComplexF64, Vb::ComplexF64, Vc::ComplexF64)
    α = exp(im * 2π / 3)
    V0 = (Va        + Vb            + Vc        ) / 3
    V1 = (Va        + α * Vb        + α^2 * Vc  ) / 3
    V2 = (Va        + α^2 * Vb      + α * Vc    ) / 3
    return V0, V1, V2
end

# ── Voltage violation (VV) ─────────────────────────────────────────────────

"""
    voltage_violations(sol_dict, math; U_nom=NaN, lo=0.9, hi=1.1)
        → Dict{String, Dict{String, Bool}}

Raise a violation alarm per (bus, phase) when the phase-to-neutral
voltage magnitude lies outside `[lo, hi] · U_nom`.

  • EN model (terminal 4 present): `|V_phase - V_neutral|`
  • KRN model (no terminal 4):     `|V_phase|` (phase-to-ground)

`U_nom` defaults to `vbase_V(math)`. Only phases 1–3 are checked.
"""
function voltage_violations(
    sol_dict::AbstractDict,
    math::AbstractDict;
    U_nom::Float64 = NaN,
    lo::Float64 = 0.9,
    hi::Float64 = 1.1,
)::Dict{String, Dict{String, Bool}}

    if isnan(U_nom)
        U_nom = vbase_V(math)
    end

    bv = bus_voltages(sol_dict, math)
    alarms = Dict{String, Dict{String, Bool}}()

    for (b, volts) in bv
        phase_alarms = Dict{String, Bool}()
        vn = get(volts, "4", nothing)
        for phase in ("1", "2", "3")
            haskey(volts, phase) || continue
            v_phase = volts[phase]
            v_pn = isnothing(vn) ? v_phase : (v_phase - vn)
            mag = abs(v_pn) * U_nom
            alarm = (mag < lo * U_nom) || (mag > hi * U_nom)
            phase_alarms[phase] = alarm
        end
        if !isempty(phase_alarms)
            alarms[b] = phase_alarms
        end
    end
    return alarms
end

# ── Neutral-to-earth voltage (NEV) ─────────────────────────────────────────

"""
    nev_alarms(sol_dict, math; U_thresh=10.0)
        → Dict{String, NamedTuple{(:U_n_V, :alarm), Tuple{Float64, Bool}}}

Per-bus neutral-to-earth voltage magnitude (Volts) and alarm flag
(IEEE Std 1100-2005: |V_n| > 10 V).

A Kron-reduced (KRN) candidate has **no explicit neutral conductor**, so its
neutral-to-earth voltage is *assumed* to be 0 V and the alarm is always
`false`. This per-bus assumed-zero dict is still returned (rather than an
empty dict) so that comparing it against a reference EN model via
`compare_decisions` correctly produces a **false negative** wherever the
reference fires a real NEV alarm that the KRN model structurally misses.
"""
function nev_alarms(
    sol_dict::AbstractDict,
    math::AbstractDict;
    U_thresh::Float64 = 10.0,
)
    bv = bus_voltages(sol_dict, math)
    result = Dict{String, NamedTuple{(:U_n_V, :alarm), Tuple{Float64, Bool}}}()

    if is_krn(math)
        # KRN: neutral assumed 0 V everywhere → alarm=false per bus. Emitted
        # (not skipped) so the candidate-vs-reference comparison can register
        # false negatives where the reference's real |V_n| exceeds U_thresh.
        for b in keys(bv)
            result[b] = (U_n_V=0.0, alarm=false)
        end
        return result
    end

    vbase = vbase_V(math)
    for (b, volts) in bv
        vn = get(volts, "4", nothing)
        if isnothing(vn)
            result[b] = (U_n_V=0.0, alarm=false)
        else
            U_n_V = abs(vn) * vbase
            result[b] = (U_n_V=U_n_V, alarm = U_n_V > U_thresh)
        end
    end
    return result
end

# ── VUF ────────────────────────────────────────────────────────────────────

"""
    vuf_alarms(sol_dict, math; thresh=2.0)
        → Dict{String, NamedTuple{(:vuf_pct, :alarm), Tuple{Float64, Bool}}}

Per-bus voltage unbalance factor (negative-to-positive sequence ratio,
EN 50160 §2.5). Alarms when VUF > `thresh` (default 2 %). Buses without
three active phases are skipped.
"""
function vuf_alarms(
    sol_dict::AbstractDict,
    math::AbstractDict;
    thresh::Float64 = 2.0,
)
    bv = bus_voltages(sol_dict, math)
    result = Dict{String, NamedTuple{(:vuf_pct, :alarm), Tuple{Float64, Bool}}}()

    for (b, volts) in bv
        (haskey(volts, "1") && haskey(volts, "2") && haskey(volts, "3")) || continue
        vn = get(volts, "4", nothing)

        Va = isnothing(vn) ? volts["1"] : (volts["1"] - vn)
        Vb = isnothing(vn) ? volts["2"] : (volts["2"] - vn)
        Vc = isnothing(vn) ? volts["3"] : (volts["3"] - vn)

        _, V1, V2 = sequence_components(Va, Vb, Vc)
        vuf_pct = abs(V1) < 1e-12 ? 0.0 : abs(V2) / abs(V1) * 100.0
        result[b] = (vuf_pct=vuf_pct, alarm = vuf_pct > thresh)
    end
    return result
end

# ── PVUR ───────────────────────────────────────────────────────────────────

"""
    pvur_alarms(sol_dict, math; thresh=2.0)
        → Dict{String, NamedTuple{(:pvur_pct, :alarm), Tuple{Float64, Bool}}}

Per-bus phase-voltage-unbalance rate (IEEE Std 141-1993, magnitude
deviation). Alarms when PVUR > `thresh` (default 2 %). Buses without
three active phases are skipped.
"""
function pvur_alarms(
    sol_dict::AbstractDict,
    math::AbstractDict;
    thresh::Float64 = 2.0,
)
    bv = bus_voltages(sol_dict, math)
    result = Dict{String, NamedTuple{(:pvur_pct, :alarm), Tuple{Float64, Bool}}}()

    for (b, volts) in bv
        (haskey(volts, "1") && haskey(volts, "2") && haskey(volts, "3")) || continue
        vn = get(volts, "4", nothing)

        Va_c = isnothing(vn) ? volts["1"] : (volts["1"] - vn)
        Vb_c = isnothing(vn) ? volts["2"] : (volts["2"] - vn)
        Vc_c = isnothing(vn) ? volts["3"] : (volts["3"] - vn)

        Va = abs(Va_c)
        Vb = abs(Vb_c)
        Vc = abs(Vc_c)

        V_avg = (Va + Vb + Vc) / 3.0
        pvur_pct = if V_avg < 1e-12
            0.0
        else
            max_dev = max(abs(Va - V_avg), abs(Vb - V_avg), abs(Vc - V_avg))
            (max_dev / V_avg) * 100.0
        end

        result[b] = (pvur_pct=pvur_pct, alarm = pvur_pct > thresh)
    end
    return result
end

# ── Decision comparison primitives ─────────────────────────────────────────

"""
    _extract_alarm(x) → Bool | nothing

Pull the alarm flag out of either a raw `Bool`, or a `NamedTuple` that has
an `:alarm` field, or `nothing` if no flag is encoded.
"""
function _extract_alarm(x)
    isnothing(x) && return nothing
    x isa Bool && return x
    x isa NamedTuple && haskey(x, :alarm) && return x.alarm
    return nothing
end

"""
    compare_decisions(d_est, d_true) → NamedTuple{(:fp,:fn,:tp,:tn), …}

Compare two alarm dicts whose values are `Bool` or `NamedTuple` with an
`:alarm` field. Keyed by bus id. Missing keys on the `est` side are
treated as a false negative when `true_alarm` is true.
"""
function compare_decisions(d_est, d_true)
    fp = fn = tp = tn = 0
    all_keys = union(keys(d_est), keys(d_true))
    for k in all_keys
        est_alarm = _extract_alarm(get(d_est, k, nothing))
        true_alarm = _extract_alarm(get(d_true, k, nothing))

        if isnothing(est_alarm)
            true_alarm === true ? (fn += 1) : (tn += 1)
            continue
        end
        if isnothing(true_alarm)
            est_alarm ? (fp += 1) : (tn += 1)
            continue
        end

        if est_alarm && true_alarm;      tp += 1
        elseif est_alarm && !true_alarm; fp += 1
        elseif !est_alarm && true_alarm; fn += 1
        else;                            tn += 1
        end
    end
    return (fp=fp, fn=fn, tp=tp, tn=tn)
end

"""
    compare_decisions_vv(d_est, d_true) → NamedTuple{(:fp,:fn,:tp,:tn), …}

Variant for VV dicts that are `bus_id → Dict(phase => Bool)`. Per
(bus, phase) confusion counts.
"""
function compare_decisions_vv(d_est, d_true)
    fp = fn = tp = tn = 0
    all_buses = union(keys(d_est), keys(d_true))
    for b in all_buses
        est_phases  = get(d_est,  b, Dict{String,Bool}())
        true_phases = get(d_true, b, Dict{String,Bool}())
        all_phases  = union(keys(est_phases), keys(true_phases))
        for ph in all_phases
            est_alarm  = get(est_phases,  ph, false)
            true_alarm = get(true_phases, ph, false)
            if est_alarm && true_alarm;      tp += 1
            elseif est_alarm && !true_alarm; fp += 1
            elseif !est_alarm && true_alarm; fn += 1
            else;                            tn += 1
            end
        end
    end
    return (fp=fp, fn=fn, tp=tp, tn=tn)
end

# ── Continuous state-space fidelity (MAPE_U) ────────────────────────────────

"""
    voltage_mape(sol_r, math_r, sol_c, math_c; vfloor=1e-9)
        → (mape_u, maxape_u)

Mean and maximum absolute percent error of the candidate's phase-to-neutral
voltage **magnitudes** versus the reference, evaluated over the set of
`(bus, phase)` keys common to both solutions:

    MAPE_U   = 100 · mean_k |‖V_c(k)‖ − ‖V_r(k)‖| / ‖V_r(k)‖
    MaxAPE_U = 100 · max_k  |‖V_c(k)‖ − ‖V_r(k)‖| / ‖V_r(k)‖

Phase-to-neutral voltages come from `flatten_phase_voltages`, which subtracts
the neutral terminal for EN models and uses phase-to-ground for KRN, so EN
and KRN are compared on a like-for-like basis.

Keys whose reference magnitude is below `vfloor` (Volts/p.u.) are skipped to
avoid a divide-by-(near-)zero. If no usable common key remains, both metrics
are returned as `NaN`.
"""
function voltage_mape(
    sol_r::AbstractDict, math_r::AbstractDict,
    sol_c::AbstractDict, math_c::AbstractDict;
    vfloor::Float64 = 1e-9,
)
    vr = flatten_phase_voltages(sol_r, math_r)
    vc = flatten_phase_voltages(sol_c, math_c)

    sum_ape = 0.0
    max_ape = 0.0
    n = 0
    for (k, vr_k) in vr
        haskey(vc, k) || continue
        denom = abs(vr_k)
        denom < vfloor && continue        # guard: skip ~0 reference magnitude
        ape = abs(abs(vc[k]) - denom) / denom
        sum_ape += ape
        ape > max_ape && (max_ape = ape)
        n += 1
    end

    n == 0 && return (NaN, NaN)
    return (100.0 * sum_ape / n, 100.0 * max_ape)
end

# ── Disagreement aggregator ────────────────────────────────────────────────

_total_mismatch(c) = c.fp + c.fn   # "disagreements" = est ≠ truth

"""
    decision_disagreement(sol_r, math_r, sol_c, math_c;
                          U_nom=NaN, vv_lo=0.9, vv_hi=1.1,
                          nev_thresh=10.0, vuf_thresh=2.0,
                          pvur_thresh=2.0) → NamedTuple

Per-decision disagreement counts between reference and candidate
solutions. The reference side is `sol_r`; the candidate side is
`sol_c`. Returns:

  (vv, nev, vuf, pvur, total, raw)

where each scalar field is the number of disagreements (FP + FN), and
`raw` is a NamedTuple of the underlying (tp, tn, fp, fn) per decision.

Also reports the continuous voltage-magnitude fidelity metrics `mape_u`
and `maxape_u` (see [`voltage_mape`](@ref)).
"""
function decision_disagreement(
    sol_r::AbstractDict, math_r::AbstractDict,
    sol_c::AbstractDict, math_c::AbstractDict;
    U_nom::Float64 = NaN,
    vv_lo::Float64 = 0.9,
    vv_hi::Float64 = 1.1,
    nev_thresh::Float64 = 10.0,
    vuf_thresh::Float64 = 2.0,
    pvur_thresh::Float64 = 2.0,
)
    if isnan(U_nom)
        U_nom = vbase_V(math_r)
    end

    vv_r  = voltage_violations(sol_r, math_r; U_nom=U_nom, lo=vv_lo, hi=vv_hi)
    vv_c  = voltage_violations(sol_c, math_c; U_nom=U_nom, lo=vv_lo, hi=vv_hi)
    nev_r = nev_alarms(sol_r, math_r; U_thresh=nev_thresh)
    nev_c = nev_alarms(sol_c, math_c; U_thresh=nev_thresh)
    vuf_r = vuf_alarms(sol_r, math_r; thresh=vuf_thresh)
    vuf_c = vuf_alarms(sol_c, math_c; thresh=vuf_thresh)
    pvur_r = pvur_alarms(sol_r, math_r; thresh=pvur_thresh)
    pvur_c = pvur_alarms(sol_c, math_c; thresh=pvur_thresh)

    # Defensive bus-key alignment for NEV: a KRN candidate has no explicit
    # neutral, so its assumed-zero (alarm=false) dict must cover every bus the
    # reference reports. Any reference bus missing on the candidate side is
    # back-filled with an assumed-0 false alarm so that a reference NEV alarm
    # is registered as a false negative rather than silently dropped.
    if is_krn(math_c)
        for b in keys(nev_r)
            haskey(nev_c, b) || (nev_c[b] = (U_n_V=0.0, alarm=false))
        end
    end

    vv_raw   = compare_decisions_vv(vv_c, vv_r)
    nev_raw  = compare_decisions(nev_c, nev_r)
    vuf_raw  = compare_decisions(vuf_c, vuf_r)
    pvur_raw = compare_decisions(pvur_c, pvur_r)

    vv   = _total_mismatch(vv_raw)
    nev  = _total_mismatch(nev_raw)
    vuf  = _total_mismatch(vuf_raw)
    pvur = _total_mismatch(pvur_raw)

    mape_u, maxape_u = voltage_mape(sol_r, math_r, sol_c, math_c)

    return (
        vv = vv,
        nev = nev,
        vuf = vuf,
        pvur = pvur,
        total = vv + nev + vuf + pvur,
        mape_u = mape_u,
        maxape_u = maxape_u,
        raw = (vv=vv_raw, nev=nev_raw, vuf=vuf_raw, pvur=pvur_raw),
    )
end

# ── Operational axis: accuracy / recall / disagreement rates ─────────────────

_acc(c) = ((t = c.tp + c.tn + c.fp + c.fn); t == 0 ? NaN : (c.tp + c.tn) / t)
_rec(c) = ((d = c.tp + c.fn);              d == 0 ? NaN :  c.tp        / d)

"""
    operational_axis(raw) → NamedTuple

Collapse pooled confusion counts into the operational coordinate of
`3_framework.tex` §III ("Operational Decision Agreement"). `raw` is a
NamedTuple `(vv, nev, vuf, pvur)`, each a `(tp, tn, fp, fn)` NamedTuple of
counts *pooled over all items and Monte-Carlo realisations* (the paper pools
before rating). Returns per-decision accuracy `ACC_D`, recall `REC_D`, and
disagreement `η_D = 1 − ACC_D`, together with the worst-case coordinate

    η_∞ = max_D η_D                         (eq. for η_∞).

Decision categories with no observed items (`ACC_D = NaN`) are excluded from
the maximum; `η_∞ = 0` when no category is observed.
"""
function operational_axis(raw)
    cats = (raw.vv, raw.nev, raw.vuf, raw.pvur)
    accs = map(_acc, cats)
    recs = map(_rec, cats)
    etas = map(a -> isnan(a) ? NaN : 1 - a, accs)
    finite = filter(!isnan, collect(etas))
    η_inf = isempty(finite) ? 0.0 : maximum(finite)
    return (
        acc = (vv = accs[1], nev = accs[2], vuf = accs[3], pvur = accs[4]),
        rec = (vv = recs[1], nev = recs[2], vuf = recs[3], pvur = recs[4]),
        eta = (vv = etas[1], nev = etas[2], vuf = etas[3], pvur = etas[4]),
        eta_inf = η_inf,
    )
end

"""
    pool_confusion(raws) → NamedTuple

Sum a collection of per-realisation `raw` NamedTuples (each
`(vv, nev, vuf, pvur)` of `(fp, fn, tp, tn)` counts) into a single pooled
`raw` NamedTuple, as required by [`operational_axis`](@ref).
"""
function pool_confusion(raws)
    z() = (fp = 0, fn = 0, tp = 0, tn = 0)
    acc = (vv = z(), nev = z(), vuf = z(), pvur = z())
    _add(a, b) = (fp = a.fp + b.fp, fn = a.fn + b.fn, tp = a.tp + b.tp, tn = a.tn + b.tn)
    for r in raws
        r === nothing && continue
        acc = (vv = _add(acc.vv, r.vv), nev = _add(acc.nev, r.nev),
               vuf = _add(acc.vuf, r.vuf), pvur = _add(acc.pvur, r.pvur))
    end
    return acc
end
