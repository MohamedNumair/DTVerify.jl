## ============================================================
## extract.jl - PMDSE dict accessors 
## ============================================================
##
## Pure dict walkers. All functions operate on `Dict{String,Any}`
## structures that look like PMDSE result and math dicts. We never
## construct math or solution objects - we only read what's there.
## ============================================================

const _N_IDX = (4,)   # neutral terminal index (PMD convention)

"""
    J_star(se_result) → Float64

Return the converged M-estimator objective `J*` from a PMDSE result dict.
Looks up `se_result["objective"]`, then `se_result["solution"]["objective"]`.
Errors if neither is present.
"""
function J_star(se_result::AbstractDict)::Float64
    if haskey(se_result, "objective") && se_result["objective"] !== nothing
        return Float64(se_result["objective"])
    end
    if haskey(se_result, "solution") &&
       haskey(se_result["solution"], "objective") &&
       se_result["solution"]["objective"] !== nothing
        return Float64(se_result["solution"]["objective"])
    end
    error("DTVerify.J_star: no `objective` field found in se_result " *
          "(checked top-level and `solution`).")
end

"""
    n_measurements(math) → Int

Count total scalar measurements stored under `math["meas"][*]["dst"]`.
Returns 0 if no `meas` block is present.
"""
function n_measurements(math::AbstractDict)::Int
    haskey(math, "meas") || return 0
    n = 0
    for (_, meas) in math["meas"]
        if haskey(meas, "dst")
            n += length(meas["dst"])
        end
    end
    return n
end

"""
    n_state_variables(math; convention=:rigorous) → Int

Count the number of free state variables of the underlying \\gls{ivren}/\\gls{ivru}
state-estimation optimisation problem. In PMDSE, the optimisation declares many
JuMP variables (`vr, vi, cr, ci, crg, cig, crd, cid, crt, cit, res`, plus
measurement-conversion auxiliaries), but only **bus voltages** `(vr, vi)` are
independent decision variables; every other variable is removed by an equality
constraint (Ohm's law, KCL, load/gen current models, residual definitions,
measurement-conversion equalities).

The free-parameter count therefore equals

    n_state = 2 · Σ_i |ungrounded_terminals_i|  −  n_ref_equality_pins

where the second term accounts for equality constraints applied at the
reference bus.

Reference-bus convention is auto-detected from the bus keys:

    • `vm` + `va` present: `constraint_mc_voltage_fixed` applies to each
        ungrounded terminal (`vr = vm*cos(va)`, `vi = vm*sin(va)`) and removes
        **two** degrees of freedom per ungrounded terminal.
    • exactly one of `vm` / `va` present: one scalar equality per ungrounded
        terminal (`|V|` fixed or angle relation fixed), removing **one** degree
        of freedom per ungrounded terminal.
    • neither present: no reference-pin subtraction.

`convention`:
  • `:rigorous` (default) - formula above; counts neutral state variables
    at non-ref buses and properly distinguishes EN from KRN.
  • `:phase_only` - backwards-compatible mode that mirrors the original
    `KPIUtils.summarise_case` formula (excludes terminal 4 everywhere and
    only counts load/gen buses; pretends EN and KRN have the same `k`).
    Use only when you need bit-for-bit agreement with the old KPI pipeline.

The reference-pin convention can be overridden via `ref_pins` (total pins
subtracted across the whole network):

    n_state_variables(math; ref_pins = 6)   # force 3-phase fully fixed ref
"""
function n_state_variables(math::AbstractDict;
                           convention::Symbol = :rigorous,
                           ref_pins::Union{Nothing, Int} = nothing)::Int
    if convention === :phase_only
        return _n_state_phase_only(math)
    elseif convention !== :rigorous
        error("DTVerify.n_state_variables: unknown convention $(convention); " *
              "use :rigorous or :phase_only.")
    end

    haskey(math, "bus") || return 0

    n_terminals = 0
    n_ref_pins_auto = 0

    for (_, bus) in math["bus"]
        terms = get(bus, "terminals", Int[])
        grounded = get(bus, "grounded", fill(false, length(terms)))
        n_free = (length(grounded) == length(terms)) ? count(!, grounded) : length(terms)
        n_terminals += n_free

        if get(bus, "bus_type", 1) == 3
            has_vm = haskey(bus, "vm")
            has_va = haskey(bus, "va")

            if has_vm && has_va
                n_ref_pins_auto += 2 * n_free
            elseif has_vm || has_va
                n_ref_pins_auto += n_free
            end
        end
    end

    n_ref_pins = if ref_pins !== nothing
        ref_pins
    else
        n_ref_pins_auto
    end

    return 2 * n_terminals - n_ref_pins
end

# Backwards-compatible KPIUtils.summarise_case logic.
function _n_state_phase_only(math::AbstractDict)::Int
    load_buses = String[]
    if haskey(math, "load")
        for (_, lo) in math["load"]
            haskey(lo, "load_bus") && push!(load_buses, string(lo["load_bus"]))
        end
    end
    gen_buses = String[]
    if haskey(math, "gen")
        for (_, g) in math["gen"]
            haskey(g, "gen_bus") && push!(gen_buses, string(g["gen_bus"]))
        end
    end
    non_zero = unique(vcat(load_buses, gen_buses))

    ref_terms_list = Vector{Int}[]
    if haskey(math, "bus")
        for (_, bus) in math["bus"]
            if get(bus, "bus_type", 1) == 3 && haskey(bus, "terminals")
                push!(ref_terms_list, collect(bus["terminals"]))
            end
        end
    end
    n_ref_phases = isempty(ref_terms_list) ?
                   0 :
                   length(setdiff(first(ref_terms_list), _N_IDX))

    n_state = 0
    if haskey(math, "bus")
        for b in non_zero
            haskey(math["bus"], b) || continue
            terms = get(math["bus"][b], "terminals", Int[])
            n_state += length(setdiff(terms, _N_IDX)) * 2
        end
    end
    return n_state - n_ref_phases
end

"""
    is_krn(math) → Bool

`true` when no bus declares terminal 4 - i.e. the math is Kron-reduced
(no explicit neutral conductor).
"""
function is_krn(math::AbstractDict)::Bool
    haskey(math, "bus") || return true
    for (_, bus) in math["bus"]
        terms = get(bus, "terminals", Int[])
        if 4 in terms
            return false
        end
    end
    return true
end

"""
    vbase_V(math) → Float64

Voltage base in **true Volts** (≈ 230 V for a standard LV feeder). Walks the
math dict directly - no Pliers `calc_bases_from_dict` dependency:

  1. `math["bus"]["sourcebus"]["vbase"]`
  2. `math["settings"]["vbases_default"]` (first value)
  3. fallback `230.0`

PMD math dicts store `vbase` (and `vbases_default`) in **scaled** units -
typically kV (e.g. 0.23094) - together with a `math["settings"]
["voltage_scale_factor"]` (= 1000.0) that converts the scaled per-unit base
back to Volts. The resolved per-unit base is multiplied by this scale factor
**only when it is present**, so that:

  • PMD data (vbase = 0.23094 kV, scale = 1000) → 230.94 V, and
  • data already expressed in Volts (no scale factor) keeps its value.
"""
function vbase_V(math::AbstractDict)::Float64
    # Scale factor that converts the stored (scaled) p.u. voltage base to
    # Volts. Defaults to 1.0 when absent so already-in-Volts data is intact.
    scale = 1.0
    if haskey(math, "settings") &&
       haskey(math["settings"], "voltage_scale_factor")
        sf = math["settings"]["voltage_scale_factor"]
        if sf isa Number && sf > 0
            scale = Float64(sf)
        end
    end

    if haskey(math, "bus") &&
       haskey(math["bus"], "sourcebus") &&
       haskey(math["bus"]["sourcebus"], "vbase")
        v = math["bus"]["sourcebus"]["vbase"]
        v isa Number && v > 0 && return Float64(v) * scale
    end
    if haskey(math, "settings") &&
       haskey(math["settings"], "vbases_default")
        vbs = math["settings"]["vbases_default"]
        if vbs isa AbstractDict && !isempty(vbs)
            v = first(values(vbs))
            v isa Number && v > 0 && return Float64(v) * scale
        elseif vbs isa AbstractVector && !isempty(vbs)
            v = first(vbs)
            v isa Number && v > 0 && return Float64(v) * scale
        end
    end
    return 230.0
end

"""
    bus_voltages(sol_dict, math) → Dict{String, Dict{String, ComplexF64}}

Extract bus voltage phasors, keyed `bus_id → terminal_str → V (ComplexF64)`.
Handles both the dictified `bus["voltage"]` form and the raw `vr`/`vi`
vector form (which we zip against `math["bus"][b]["terminals"]`).

Accepts dicts that look like `sol_dict["bus"]…` or
`sol_dict["solution"]["bus"]…`.
"""
function bus_voltages(sol_dict::AbstractDict, math::AbstractDict)
    buses = if haskey(sol_dict, "bus")
        sol_dict["bus"]
    elseif haskey(sol_dict, "solution") && haskey(sol_dict["solution"], "bus")
        sol_dict["solution"]["bus"]
    else
        return Dict{String, Dict{String, ComplexF64}}()
    end

    result = Dict{String, Dict{String, ComplexF64}}()
    for (b, bus_data) in buses
        if haskey(bus_data, "voltage")
            result[b] = Dict{String, ComplexF64}(
                string(k) => ComplexF64(v) for (k, v) in bus_data["voltage"]
            )
        elseif haskey(bus_data, "vr") && haskey(bus_data, "vi")
            terms = haskey(math, "bus") && haskey(math["bus"], b) ?
                    get(math["bus"][b], "terminals", Int[]) :
                    Int[]
            vr = bus_data["vr"]
            vi = bus_data["vi"]
            volt_dict = Dict{String, ComplexF64}()
            for (i, t) in enumerate(terms)
                if i <= length(vr) && i <= length(vi)
                    volt_dict[string(t)] = ComplexF64(vr[i], vi[i])
                end
            end
            result[b] = volt_dict
        end
    end
    return result
end

"""
    flatten_phase_voltages(sol_dict, math) → Dict{Tuple{String,String}, ComplexF64}

Phase-to-neutral voltage phasors for phases 1–3 only, keyed `(bus, phase)`.
For EN models we subtract terminal 4; for KRN we use phase-to-ground.
"""
function flatten_phase_voltages(sol_dict::AbstractDict, math::AbstractDict)
    bv = bus_voltages(sol_dict, math)
    out = Dict{Tuple{String,String}, ComplexF64}()
    for (b, volts) in bv
        vn = get(volts, "4", nothing)
        for ph in ("1", "2", "3")
            haskey(volts, ph) || continue
            v = isnothing(vn) ? volts[ph] : (volts[ph] - vn)
            out[(b, ph)] = v
        end
    end
    return out
end

"""
    flatten_neutral_voltages(sol_dict, math) → Dict{String, ComplexF64}

Neutral-terminal (terminal 4) voltage phasors per bus. Returns an empty
dict for Kron-reduced models.
"""
function flatten_neutral_voltages(sol_dict::AbstractDict, math::AbstractDict)
    out = Dict{String, ComplexF64}()
    is_krn(math) && return out
    bv = bus_voltages(sol_dict, math)
    for (b, volts) in bv
        haskey(volts, "4") || continue
        out[b] = volts["4"]
    end
    return out
end
