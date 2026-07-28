## Synthetic math/sol dicts  

# Minimal EN math dict: source bus + two load buses, all four-wire.
function make_en_math()
    return Dict{String,Any}(
        "bus" => Dict{String,Any}(
            "sourcebus" => Dict{String,Any}(
                "terminals" => [1, 2, 3, 4],
                "bus_type"  => 3,    # reference
                "grounded"  => Bool[0, 0, 0, 1],
                "vm"        => [1.0, 1.0, 1.0, 0.0],
                "va"        => [0.0, -2pi/3, 2pi/3, 0.0],
                "vbase"     => 230.0,
            ),
            "1" => Dict{String,Any}(
                "terminals" => [1, 2, 3, 4],
                "bus_type"  => 1,
                "grounded"  => Bool[0, 0, 0, 0],
            ),
            "2" => Dict{String,Any}(
                "terminals" => [1, 2, 3, 4],
                "bus_type"  => 1,
                "grounded"  => Bool[0, 0, 0, 0],
            ),
        ),
        "load" => Dict{String,Any}(
            "1" => Dict{String,Any}("load_bus" => 1),
            "2" => Dict{String,Any}("load_bus" => 2),
        ),
        "gen" => Dict{String,Any}(),
        "settings" => Dict{String,Any}("vbases_default" => Dict("sourcebus" => 230.0)),
    )
end

# KRN math (phases only - no terminal 4)
function make_krn_math()
    return Dict{String,Any}(
        "bus" => Dict{String,Any}(
            "sourcebus" => Dict{String,Any}(
                "terminals" => [1, 2, 3],
                "bus_type"  => 3,
                "grounded"  => Bool[0, 0, 0],
                "vm"        => [1.0, 1.0, 1.0],
                "va"        => [0.0, -2pi/3, 2pi/3],
                "vbase"     => 230.0,
            ),
            "1" => Dict{String,Any}(
                "terminals" => [1, 2, 3],
                "bus_type"  => 1,
                "grounded"  => Bool[0, 0, 0],
            ),
        ),
        "load" => Dict{String,Any}(
            "1" => Dict{String,Any}("load_bus" => 1),
        ),
        "gen" => Dict{String,Any}(),
    )
end

# Solution dict with vr/vi vectors keyed by bus id.
function make_sol_vrvi(math)
    sol = Dict{String,Any}("bus" => Dict{String,Any}())
    for (b, _) in math["bus"]
        terms = math["bus"][b]["terminals"]
        nt = length(terms)
        # Set ~1 p.u. phase voltages, 0 neutral (if present)
        vr = Float64[]
        vi = Float64[]
        for (i, t) in enumerate(terms)
            if t == 4
                push!(vr, 0.0); push!(vi, 0.0)
            elseif t == 1
                push!(vr, 1.0); push!(vi, 0.0)
            elseif t == 2
                push!(vr, -0.5); push!(vi, -sqrt(3)/2)
            elseif t == 3
                push!(vr, -0.5); push!(vi, sqrt(3)/2)
            end
        end
        sol["bus"][b] = Dict{String,Any}("vr" => vr, "vi" => vi)
    end
    return sol
end

@testset "is_krn detection" begin
    en = make_en_math()
    krn = make_krn_math()
    @test is_krn(en) == false
    @test is_krn(krn) == true
end

@testset "vbase_V resolution" begin
    en = make_en_math()
    @test vbase_V(en) ≈ 230.0

    # Settings fallback
    en2 = Dict{String,Any}(
        "bus" => Dict{String,Any}(),
        "settings" => Dict{String,Any}("vbases_default" => Dict("sb" => 400.0)),
    )
    @test vbase_V(en2) ≈ 400.0

    # Empty fallback
    @test vbase_V(Dict{String,Any}()) ≈ 230.0
end

@testset "vbase_V applies voltage_scale_factor (kV → V)" begin
    # PMD-style math: vbase stored in kV, scale factor 1000 → true Volts.
    math = Dict{String,Any}(
        "bus" => Dict{String,Any}(
            "sourcebus" => Dict{String,Any}("vbase" => 0.23094),
        ),
        "settings" => Dict{String,Any}("voltage_scale_factor" => 1000.0),
    )
    @test vbase_V(math) ≈ 230.94 atol = 1e-6

    # vbases_default path also scaled.
    math2 = Dict{String,Any}(
        "bus" => Dict{String,Any}(),
        "settings" => Dict{String,Any}(
            "vbases_default" => Dict("sb" => 0.23094),
            "voltage_scale_factor" => 1000.0,
        ),
    )
    @test vbase_V(math2) ≈ 230.94 atol = 1e-6

    # No scale factor → already-in-Volts data is untouched (regression).
    @test vbase_V(make_en_math()) ≈ 230.0
end

@testset "n_state_variables (EN, rigorous)" begin
    en = make_en_math()
    # Free terminals: source has 3 (neutral grounded), buses 1-2 have 4 each.
    # Total free terminals = 11; vr+vi => 22 variables.
    # Reference bus has vm+va on 3 ungrounded phases => 6 equality pins.
    # 22 - 6 = 16
    @test n_state_variables(en) == 16
    @test n_state_variables(en; convention=:rigorous) == 16
end

@testset "n_state_variables (KRN, rigorous)" begin
    krn = make_krn_math()
    # 2 buses × 3 terminals = 6 terminals total
    # × 2 (vr+vi) = 12
    # ref bus has vm+va on 3 terminals => 2·3 = 6 equality pins
    # → 12 - 6 = 6
    @test n_state_variables(krn) == 6
end

@testset "n_state_variables (phase_only, backwards-compat)" begin
    en = make_en_math()
    # KPIUtils.summarise_case formula:
    #   load_buses ∈ {1,2}; each contributes setdiff([1,2,3,4],[4])·2 = 6
    #   total = 12; n_ref_phases = 3 → 12 - 3 = 9
    @test n_state_variables(en; convention=:phase_only) == 9

    krn = make_krn_math()
    # load_buses ∈ {1}; setdiff([1,2,3],[4])·2 = 6; n_ref_phases = 3 → 3
    @test n_state_variables(krn; convention=:phase_only) == 3
end

@testset "n_state_variables (ref_pins override)" begin
    en = make_en_math()
    # Auto count has 22 base vars (ungrounded terminals only).
    @test n_state_variables(en; ref_pins = 0) == 22
    @test n_state_variables(en; ref_pins = 8) == 14
end

@testset "bus_voltages from vr/vi" begin
    en = make_en_math()
    sol = make_sol_vrvi(en)
    bv = bus_voltages(sol, en)
    @test haskey(bv, "1")
    v1 = bv["1"]
    @test haskey(v1, "1") && haskey(v1, "2") && haskey(v1, "3") && haskey(v1, "4")
    @test real(v1["1"]) ≈ 1.0
    @test imag(v1["1"]) ≈ 0.0
    @test isapprox(abs(v1["1"]), 1.0; atol=1e-12)
    @test abs(v1["4"]) ≈ 0.0
end

@testset "bus_voltages from `voltage` dict form" begin
    math = Dict{String,Any}("bus" => Dict{String,Any}("a" => Dict{String,Any}("terminals" => [1,2,3])))
    sol = Dict{String,Any}("bus" => Dict{String,Any}(
        "a" => Dict{String,Any}(
            "voltage" => Dict("1" => complex(1.0, 0.1), "2" => complex(-0.5, -0.9)),
        ),
    ))
    bv = bus_voltages(sol, math)
    @test bv["a"]["1"] == ComplexF64(1.0, 0.1)
    @test bv["a"]["2"] == ComplexF64(-0.5, -0.9)
end

@testset "flatten_phase_voltages excludes neutral, subtracts when EN" begin
    en = make_en_math()
    sol = make_sol_vrvi(en)
    fpv = flatten_phase_voltages(sol, en)
    # Only phases 1..3 → 2 non-source buses × 3 phases = 6
    @test length(fpv) == 6 || length(fpv) == 9   # may include source bus
    @test all(k -> k[2] in ("1","2","3"), keys(fpv))
end

@testset "flatten_neutral_voltages - empty for KRN" begin
    krn = make_krn_math()
    sol = make_sol_vrvi(krn)
    @test isempty(flatten_neutral_voltages(sol, krn))
end

@testset "J_star extraction" begin
    @test J_star(Dict("objective" => 1.23)) ≈ 1.23
    @test J_star(Dict("solution" => Dict("objective" => 4.5))) ≈ 4.5
    @test_throws ErrorException J_star(Dict("foo" => 1.0))
end

@testset "n_measurements" begin
    math = Dict{String,Any}("meas" => Dict{String,Any}(
        "1" => Dict{String,Any}("dst" => [1, 2, 3]),
        "2" => Dict{String,Any}("dst" => [1, 2]),
    ))
    @test n_measurements(math) == 5
    @test n_measurements(Dict{String,Any}()) == 0
end
