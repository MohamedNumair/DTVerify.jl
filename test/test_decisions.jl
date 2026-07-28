## Synthetic-decision tests  

# Build a four-wire math + sol with controllable per-bus phase magnitudes.
# Voltages are in p.u. and the alarm code multiplies by U_nom internally.
function make_math_en(bus_ids::Vector{String})
    bus = Dict{String,Any}(
        "sourcebus" => Dict{String,Any}(
            "terminals" => [1, 2, 3, 4],
            "bus_type" => 3,
            "vbase" => 230.0,
        ),
    )
    load = Dict{String,Any}()
    for (i, b) in enumerate(bus_ids)
        bus[b] = Dict{String,Any}("terminals" => [1, 2, 3, 4], "bus_type" => 1)
        load[string(i)] = Dict{String,Any}("load_bus" => i)
    end
    return Dict{String,Any}(
        "bus"  => bus,
        "load" => load,
        "gen"  => Dict{String,Any}(),
    )
end

# Build a sol where each bus has a configurable phase magnitude in p.u.
function make_sol_from_mags(math, mags_per_bus::Dict{String, NTuple{3,Float64}}; vn::Float64=0.0)
    sol = Dict{String,Any}("bus" => Dict{String,Any}())
    α = exp(im * 2π/3)
    for (b, mags) in mags_per_bus
        m1, m2, m3 = mags
        # phase-to-ground = phase mag × balanced angle, then add neutral offset
        v1 = m1 * 1.0   + vn
        v2 = m2 * α^2   + vn
        v3 = m3 * α     + vn
        v4 = ComplexF64(vn, 0.0)
        sol["bus"][b] = Dict{String,Any}(
            "vr" => [real(v1), real(v2), real(v3), real(v4)],
            "vi" => [imag(v1), imag(v2), imag(v3), imag(v4)],
        )
    end
    return sol
end

# EN math whose vbase is stored in kV with a 1000× voltage_scale_factor,
# i.e. a PMD-style dict where vbase_V must resolve to ~230.94 V.
function make_math_en_scaled(bus_ids::Vector{String})
    m = make_math_en(bus_ids)
    m["bus"]["sourcebus"]["vbase"] = 0.23094
    m["settings"] = Dict{String,Any}("voltage_scale_factor" => 1000.0)
    return m
end

# KRN counterpart (phases only, no terminal 4), same buses + scaled base.
function make_math_krn_scaled(bus_ids::Vector{String})
    bus = Dict{String,Any}(
        "sourcebus" => Dict{String,Any}(
            "terminals" => [1, 2, 3], "bus_type" => 3, "vbase" => 0.23094,
        ),
    )
    load = Dict{String,Any}()
    for (i, b) in enumerate(bus_ids)
        bus[b] = Dict{String,Any}("terminals" => [1, 2, 3], "bus_type" => 1)
        load[string(i)] = Dict{String,Any}("load_bus" => i)
    end
    return Dict{String,Any}(
        "bus" => bus, "load" => load, "gen" => Dict{String,Any}(),
        "settings" => Dict{String,Any}("voltage_scale_factor" => 1000.0),
    )
end

# KRN solution: phase-to-ground voltages only (no neutral terminal).
function make_krn_sol_from_mags(mags_per_bus::Dict{String, NTuple{3,Float64}})
    sol = Dict{String,Any}("bus" => Dict{String,Any}())
    α = exp(im * 2π/3)
    for (b, mags) in mags_per_bus
        m1, m2, m3 = mags
        v1 = m1 * 1.0; v2 = m2 * α^2; v3 = m3 * α
        sol["bus"][b] = Dict{String,Any}(
            "vr" => [real(v1), real(v2), real(v3)],
            "vi" => [imag(v1), imag(v2), imag(v3)],
        )
    end
    return sol
end

@testset "nev_alarms: scaled base - ~13 V fires, ~9 V silent" begin
    math = make_math_en_scaled(["A"])
    @test vbase_V(math) ≈ 230.94 atol = 1e-6
    # ~13 V neutral: vn_pu = 13 / 230.94 ≈ 0.0563
    sol_fire = make_sol_from_mags(math, Dict("A" => (1.0,1.0,1.0)); vn = 13.0/230.94)
    nev_fire = nev_alarms(sol_fire, math; U_thresh=10.0)
    @test nev_fire["A"].U_n_V ≈ 13.0 atol = 0.1
    @test nev_fire["A"].alarm == true
    # ~9 V neutral: below the 10 V threshold ⇒ silent
    sol_ok = make_sol_from_mags(math, Dict("A" => (1.0,1.0,1.0)); vn = 9.0/230.94)
    nev_ok = nev_alarms(sol_ok, math; U_thresh=10.0)
    @test nev_ok["A"].U_n_V ≈ 9.0 atol = 0.1
    @test nev_ok["A"].alarm == false
end

@testset "nev: KRN candidate vs EN reference ⇒ FN > 0" begin
    bus_ids = ["A", "B"]
    en_math  = make_math_en_scaled(bus_ids)
    krn_math = make_math_krn_scaled(bus_ids)
    # Reference (EN) has a real ~13 V neutral at both buses ⇒ fires NEV.
    sol_r = make_sol_from_mags(en_math,
        Dict("A" => (1.0,1.0,1.0), "B" => (1.0,1.0,1.0)); vn = 13.0/230.94)
    # KRN candidate: phase-to-ground, structurally cannot see neutral.
    sol_c = make_krn_sol_from_mags(Dict("A" => (1.0,1.0,1.0), "B" => (1.0,1.0,1.0)))

    nev_r = nev_alarms(sol_r, en_math; U_thresh=10.0)
    nev_c = nev_alarms(sol_c, krn_math; U_thresh=10.0)
    raw = compare_decisions(nev_c, nev_r)
    @test raw.fn >= 2          # both reference NEV alarms missed by KRN
    @test raw.tp == 0

    # End-to-end via decision_disagreement (with key-alignment back-fill).
    d = decision_disagreement(sol_r, en_math, sol_c, krn_math)
    @test d.raw.nev.fn >= 2
    @test d.nev >= 2
end

@testset "VV confusion counts invariant to vbase scale" begin
    bus_ids = ["A", "B"]
    math_v  = make_math_en(bus_ids)            # vbase 230, no scale factor
    math_kv = make_math_en_scaled(bus_ids)     # vbase 0.23094 kV, scale 1000
    @test vbase_V(math_v)  ≈ 230.0
    @test vbase_V(math_kv) ≈ 230.94 atol = 1e-6

    mags_r = Dict("A" => (1.00,1.00,1.00), "B" => (1.00,1.00,1.00))
    mags_c = Dict("A" => (0.85,1.00,1.15), "B" => (1.00,0.80,1.00))
    sol_r = make_sol_from_mags(math_v, mags_r)
    sol_c = make_sol_from_mags(math_v, mags_c)

    # Compare VV confusion using each math's resolved base (U_nom default).
    vv_r_v  = voltage_violations(sol_r, math_v)
    vv_c_v  = voltage_violations(sol_c, math_v)
    vv_r_kv = voltage_violations(sol_r, math_kv)
    vv_c_kv = voltage_violations(sol_c, math_kv)
    c_v  = compare_decisions_vv(vv_c_v,  vv_r_v)
    c_kv = compare_decisions_vv(vv_c_kv, vv_r_kv)
    @test c_v == c_kv          # TP/TN/FP/FN identical regardless of base
end

@testset "voltage_violations: in-band ⇒ no alarms" begin
    math = make_math_en(["A"])
    sol  = make_sol_from_mags(math, Dict("A" => (1.00, 1.00, 1.00)))
    vv = voltage_violations(sol, math; U_nom=230.0, lo=0.9, hi=1.1)
    for (_, phases) in vv
        for (_, alarm) in phases
            @test alarm == false
        end
    end
end

@testset "voltage_violations: out-of-band ⇒ alarms" begin
    math = make_math_en(["A"])
    sol  = make_sol_from_mags(math, Dict("A" => (0.85, 1.00, 1.15)))
    vv = voltage_violations(sol, math; U_nom=230.0, lo=0.9, hi=1.1)
    @test vv["A"]["1"] == true   # 0.85 < 0.9
    @test vv["A"]["2"] == false
    @test vv["A"]["3"] == true   # 1.15 > 1.1
end

@testset "nev_alarms: KRN returns all-false" begin
    krn_math = Dict{String,Any}(
        "bus" => Dict{String,Any}(
            "sourcebus" => Dict{String,Any}("terminals" => [1,2,3], "bus_type"=>3, "vbase"=>230.0),
            "A" => Dict{String,Any}("terminals" => [1,2,3], "bus_type"=>1),
        ),
        "load" => Dict{String,Any}("1" => Dict("load_bus" => "A")),
        "gen"  => Dict{String,Any}(),
    )
    sol = Dict{String,Any}("bus" => Dict{String,Any}(
        "A" => Dict{String,Any}("vr" => [1.0, -0.5, -0.5], "vi" => [0.0, -0.87, 0.87]),
        "sourcebus" => Dict{String,Any}("vr" => [1.0, -0.5, -0.5], "vi" => [0.0, -0.87, 0.87]),
    ))
    nev = nev_alarms(sol, krn_math)
    for (_, v) in nev
        @test v.alarm == false
        @test v.U_n_V == 0.0
    end
end

@testset "nev_alarms: EN flags large neutral excursion" begin
    math = make_math_en(["A"])
    # neutral offset 0.1 p.u. × 230 V = 23 V > 10 V threshold
    sol = make_sol_from_mags(math, Dict("A" => (1.0, 1.0, 1.0)); vn=0.1)
    nev = nev_alarms(sol, math; U_thresh=10.0)
    @test nev["A"].alarm == true
    @test nev["A"].U_n_V > 10.0
end

@testset "vuf_alarms: balanced ⇒ low VUF" begin
    math = make_math_en(["A"])
    sol = make_sol_from_mags(math, Dict("A" => (1.0, 1.0, 1.0)))
    vuf = vuf_alarms(sol, math; thresh=2.0)
    @test vuf["A"].vuf_pct < 0.5
    @test vuf["A"].alarm == false
end

@testset "vuf_alarms: unbalanced ⇒ alarm" begin
    math = make_math_en(["A"])
    sol = make_sol_from_mags(math, Dict("A" => (1.0, 0.7, 1.2)))
    vuf = vuf_alarms(sol, math; thresh=2.0)
    @test vuf["A"].alarm == true
end

@testset "pvur_alarms: unbalanced ⇒ alarm" begin
    math = make_math_en(["A"])
    sol = make_sol_from_mags(math, Dict("A" => (1.0, 0.7, 1.2)))
    pvur = pvur_alarms(sol, math; thresh=2.0)
    @test pvur["A"].alarm == true
end

@testset "decision_disagreement: identical ⇒ zero" begin
    math = make_math_en(["A", "B"])
    sol  = make_sol_from_mags(math, Dict("A" => (1.0,1.0,1.0), "B" => (1.0,1.0,1.0)))
    d = decision_disagreement(sol, math, sol, math)
    @test d.total == 0
    @test d.vv == 0
    @test d.nev == 0
    @test d.vuf == 0
    @test d.pvur == 0
end

@testset "decision_disagreement: counts grow with mismatch" begin
    math = make_math_en(["A"])
    sol_g = make_sol_from_mags(math, Dict("A" => (1.0, 1.0, 1.0)))                # truth: no alarms
    sol_c = make_sol_from_mags(math, Dict("A" => (0.80, 0.80, 1.20)); vn=0.1)     # candidate alarms VV / NEV / unbalance
    d = decision_disagreement(sol_g, math, sol_c, math)
    @test d.total > 0
    @test d.vv >= 2   # phase 1 and 3 should disagree (and possibly 2)
    @test d.nev == 1
end

@testset "voltage_mape: identical ⇒ 0, ×1.05 ⇒ ≈5%" begin
    math = make_math_en(["A", "B"])
    sol_r = make_sol_from_mags(math,
        Dict("A" => (1.00, 0.97, 1.03), "B" => (0.99, 1.01, 1.00)))
    # Identical candidate ⇒ MAPE_U == 0
    m0, mx0 = voltage_mape(sol_r, math, sol_r, math)
    @test m0 == 0.0
    @test mx0 == 0.0

    # Candidate scaled by 1.05 on every phase magnitude ⇒ MAPE_U ≈ 5.0
    sol_c = make_sol_from_mags(math,
        Dict("A" => (1.05, 0.97*1.05, 1.03*1.05), "B" => (0.99*1.05, 1.01*1.05, 1.00*1.05)))
    m1, mx1 = voltage_mape(sol_r, math, sol_c, math)
    @test m1 ≈ 5.0 atol = 1e-6
    @test mx1 ≈ 5.0 atol = 1e-6

    # decision_disagreement surfaces the same metrics.
    d = decision_disagreement(sol_r, math, sol_c, math)
    @test d.mape_u ≈ 5.0 atol = 1e-6
    @test d.maxape_u ≈ 5.0 atol = 1e-6
end

@testset "verify_mc: MAPE_U aggregation (mean / max)" begin
    math = make_math_en(["A", "B"])
    sol_r = make_sol_from_mags(math,
        Dict("A" => (1.00, 0.97, 1.03), "B" => (0.99, 1.01, 1.00)))
    sol_c2 = make_sol_from_mags(math,
        Dict("A" => (1.02, 0.97*1.02, 1.03*1.02), "B" => (0.99*1.02, 1.01*1.02, 1.00*1.02)))
    sol_c8 = make_sol_from_mags(math,
        Dict("A" => (1.08, 0.97*1.08, 1.03*1.08), "B" => (0.99*1.08, 1.01*1.08, 1.00*1.08)))
    se_r = [merge(sol_r, Dict("objective" => 1.0)), merge(sol_r, Dict("objective" => 1.0))]
    se_c = [merge(sol_c2, Dict("objective" => 2.0)), merge(sol_c8, Dict("objective" => 2.0))]
    res = verify_mc(se_r, math, se_c, math)
    # mean of 2% and 8% per-realisation MAPE = 5%; max of maxAPE = 8%
    @test res.mape_u ≈ 5.0 atol = 1e-6
    @test res.maxape_u ≈ 8.0 atol = 1e-6
end

@testset "compare_decisions: missing keys count correctly" begin
    d_truth = Dict("A" => (U_n_V=20.0, alarm=true), "B" => (U_n_V=2.0, alarm=false))
    d_est   = Dict("A" => (U_n_V=20.0, alarm=true))
    c = compare_decisions(d_est, d_truth)
    @test c.tp == 1
    @test c.tn == 1   # B truth=false missing on est ⇒ tn
end

@testset "operational_axis: η_D = 1 − ACC_D and η_∞" begin
    raw = (
        vv   = (fp=0, fn=0, tp=5, tn=5),    # ACC 1.00 ⇒ η 0
        nev  = (fp=1, fn=1, tp=2, tn=6),    # ACC 8/10 = 0.80 ⇒ η 0.20
        vuf  = (fp=0, fn=0, tp=0, tn=10),   # ACC 1.00 ⇒ η 0
        pvur = (fp=0, fn=0, tp=0, tn=10),   # ACC 1.00 ⇒ η 0
    )
    op = operational_axis(raw)
    @test op.acc.vv  ≈ 1.0
    @test op.acc.nev ≈ 0.8
    @test op.eta.nev ≈ 0.2
    @test op.rec.nev ≈ 2 / 3                # tp/(tp+fn)
    @test op.eta_inf ≈ 0.2                  # worst-affected category
end

@testset "operational_axis: unobserved category excluded from η_∞" begin
    raw = (
        vv   = (fp=0, fn=0, tp=0, tn=0),    # no items ⇒ ACC NaN ⇒ excluded
        nev  = (fp=0, fn=0, tp=1, tn=9),    # ACC 1.0 ⇒ η 0
        vuf  = (fp=0, fn=0, tp=0, tn=0),
        pvur = (fp=0, fn=0, tp=0, tn=0),
    )
    op = operational_axis(raw)
    @test isnan(op.acc.vv)
    @test op.eta_inf ≈ 0.0
end

@testset "pool_confusion: sums per-realisation counts" begin
    r1 = (vv=(fp=1,fn=0,tp=2,tn=3), nev=(fp=0,fn=1,tp=0,tn=4),
          vuf=(fp=0,fn=0,tp=0,tn=5), pvur=(fp=0,fn=0,tp=0,tn=5))
    r2 = (vv=(fp=0,fn=1,tp=1,tn=4), nev=(fp=1,fn=0,tp=1,tn=3),
          vuf=(fp=0,fn=0,tp=0,tn=5), pvur=(fp=0,fn=0,tp=0,tn=5))
    p = pool_confusion([r1, r2])
    @test p.vv  == (fp=1, fn=1, tp=3, tn=7)
    @test p.nev == (fp=1, fn=1, tp=1, tn=7)
end

@testset "verify_mc: verdict schema is present" begin
    math = make_math_en(["A", "B"])
    sol  = make_sol_from_mags(math, Dict("A" => (1.0,1.0,1.0), "B" => (1.0,1.0,1.0)))
    se_r = [merge(sol, Dict("objective" => 1.0)) for _ in 1:5]
    se_c = [merge(sol, Dict("objective" => 1.0)) for _ in 1:5]
    res = verify_mc(se_r, math, se_c, math; M = 100)
    @test res.verdict in (:verified, :tolerable, :consequential)
    @test res.eta_inf ≈ 0.0                 # identical solutions ⇒ no disagreement
    @test res.d_bar ≈ 0.0                   # identical objectives ⇒ zero gap
    @test res.verdict == :verified          # equivalent + no alarm disagreement
    @test haskey(res, :s_d) && haskey(res, :equivalence_holds)
end
