## Synthetic-data tests for the statistical axis: fidelity gap, one-sided
## equivalence test, and the three-region joint verdict (no PMDSE dependency).

using Statistics
using Distributions: TDist, quantile

@testset "fidelity_gap: scalars and vectors" begin
    # Identical objectives ⇒ zero gap, zero spread (n ≥ 2).
    fg = fidelity_gap([1.0, 2.0, 3.0], [1.0, 2.0, 3.0], 10.0)
    @test fg.d_bar ≈ 0.0
    @test fg.s_d ≈ 0.0

    # Uniform per-measurement shift of +0.5 loss over M = 10 ⇒ d̄ = 0.05, s_d = 0.
    fg2 = fidelity_gap([1.0, 2.0, 3.0], [1.5, 2.5, 3.5], 10.0)
    @test fg2.d_bar ≈ 0.05
    @test fg2.s_d ≈ 0.0 atol = 1e-12

    # Single realisation ⇒ s_d undefined (NaN); d̄ = (Jc − Jr)/M.
    fg1 = fidelity_gap(1.0, 5.0, 4.0)
    @test fg1.d_bar ≈ 1.0
    @test isnan(fg1.s_d)
    @test fg1.d_j == [1.0]

    # Per-realisation M vector.
    fgv = fidelity_gap([0.0, 0.0], [1.0, 2.0], [10.0, 20.0])
    @test fgv.d_bar ≈ 0.1
    @test all(fgv.d_j .≈ 0.1)

    # Non-positive M ⇒ NaN gap, not an error.
    @test isnan(fidelity_gap(0.0, 1.0, 0.0).d_bar)

    # Length mismatch ⇒ error.
    @test_throws ErrorException fidelity_gap([1.0, 2.0], [1.0], 10.0)
end

@testset "equivalence_test: one-sided upper bound vs δ_eq" begin
    # d̄ = 0, s_d = 0, n = 30 ⇒ bound = 0 < δ_eq ⇒ equivalence holds.
    eq0 = equivalence_test(0.0, 0.0, 30)
    @test eq0.holds
    @test eq0.bound ≈ 0.0

    # d̄ well above the margin ⇒ fails.
    @test !equivalence_test(1.0, 0.01, 30).holds

    # Bound formula: d̄ + t_{0.95,n−1}·s_d/√n.
    d̄, s_d, n = 0.02, 0.05, 30
    eq = equivalence_test(d̄, s_d, n)
    @test eq.bound ≈ d̄ + eq.t_crit * s_d / sqrt(n)
    @test eq.t_crit ≈ quantile(TDist(n - 1), 0.95)

    # A tight sample just under the margin holds; widening s_d pushes it over.
    @test equivalence_test(0.045, 0.001, 30).holds
    @test !equivalence_test(0.045, 0.05, 30).holds

    # n < 2 ⇒ not testable (holds = false, bound NaN).
    eq1 = equivalence_test(0.0, NaN, 1)
    @test !eq1.holds
    @test isnan(eq1.bound)

    # Custom knobs propagate.
    eqc = equivalence_test(0.08, 0.0, 30; δ_eq = 0.10)
    @test eqc.holds && eqc.δ_eq == 0.10
end

@testset "joint_verdict: three regions (eq. joint-verdict)" begin
    @test joint_verdict(true,  0.00) == :verified
    @test joint_verdict(false, 0.00) == :tolerable
    @test joint_verdict(true,  0.50) == :consequential   # η_∞ > α overrides
    @test joint_verdict(false, 0.20) == :consequential
    # η_∞ exactly at α is within tolerance (η_∞ ≤ α).
    @test joint_verdict(true, 0.10) == :verified
    # Custom operational tolerance.
    @test joint_verdict(false, 0.04; α = 0.02) == :consequential
    @test joint_verdict(false, 0.04; α = 0.05) == :tolerable
end

@testset "statistical axis end-to-end" begin
    # Effectively equivalent candidate: per-measurement gap ~1e-6.
    Jr = collect(1.0:1.0:30.0)
    Jc_equiv = Jr .+ 1e-4                     # ΔJ = 1e-4 over M = 100 ⇒ d ≈ 1e-6
    fg_e = fidelity_gap(Jr, Jc_equiv, 100.0)
    eq_e = equivalence_test(fg_e.d_bar, fg_e.s_d, length(Jr))
    @test eq_e.holds
    @test joint_verdict(eq_e.holds, 0.0) == :verified

    # Distinguishable candidate: gap ~0.5 per measurement ⇒ fails equivalence,
    # but with no alarm disagreement it is tolerable, not consequential.
    Jc_dist = Jr .+ 50.0                      # ΔJ = 50 over M = 100 ⇒ d = 0.5
    fg_d = fidelity_gap(Jr, Jc_dist, 100.0)
    eq_d = equivalence_test(fg_d.d_bar, fg_d.s_d, length(Jr))
    @test !eq_d.holds
    @test fg_d.d_bar ≈ 0.5
    @test joint_verdict(eq_d.holds, 0.0)  == :tolerable
    @test joint_verdict(eq_d.holds, 0.25) == :consequential
end
