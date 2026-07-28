## kron_reduction_case.jl - Reproduce the paper's Kron-reduction
## Run:
##   cd DTVerify.jl/examples
##   julia --project=. -e 'using Pkg; Pkg.instantiate()'   # first time only
##   julia --project=. kron_reduction_case.jl

using Pkg
Pkg.activate(@__DIR__)

using CSV, DataFrames, Printf
using DTVerify

const DATA = joinpath(@__DIR__, "..", "data", "kron_case")

## Declared deployment knobs (paper defaults)
const δ_EQ = 0.05    # equivalence margin, loss units per measurement
const α_EQ = 0.05    # one-sided test level (95 % confidence)
const α_OP = 0.10    # operational tolerance on the worst-alarm disagreement

## ---- Statistical axis ------------------------------------------------------
# Per-realisation paired objectives -> fidelity gap and equivalence test.

dj = CSV.read(joinpath(DATA, "kron_dj_n30.csv"), DataFrame)
n  = nrow(dj)

gap = fidelity_gap(dj.J_star_r, dj.J_star_c, dj.n_meas)
eq  = equivalence_test(gap.d_bar, gap.s_d, n; δ_eq = δ_EQ, α_eq = α_EQ)

## ---- Operational axis ------------------------------------------------------
# Pooled confusion counts over all items and realisations -> per-decision
# accuracy/recall, worst-alarm disagreement, and the joint verdict.

kpi = CSV.read(joinpath(DATA, "block3_en_vs_krn.csv"), DataFrame)[1, :]

acc(a) = (kpi["TP_$a"] + kpi["TN_$a"]) /
         (kpi["TP_$a"] + kpi["TN_$a"] + kpi["FP_$a"] + kpi["FN_$a"])
rec(a) = (kpi["TP_$a"] + kpi["FN_$a"]) == 0 ? 1.0 :
         kpi["TP_$a"] / (kpi["TP_$a"] + kpi["FN_$a"])

const ALARMS = ("VV", "NEV", "VUF", "PVUR")
η_D   = Dict(a => 1 - acc(a) for a in ALARMS)
η_inf = maximum(values(η_D))

verdict = joint_verdict(eq.holds, η_inf; α = α_OP)
 