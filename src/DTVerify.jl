"""
    DTVerify

Decision-aware Digital Twin Verification (DTV) for distribution-system state
estimation (DSSE), implementing the framework of the companion paper,
Section 2 (*The Digital Twin Verification Framework*).

Given paired reference (`M_r`) and candidate (`M_c`) DSSE result dicts and
their math dicts, the package scores each comparison on two axes:

  • Statistical - the per-measurement fidelity gap `d̄` and a one-sided paired
    equivalence test at `(δ_eq, α_eq)`.
  • Operational - four standard-based alarms (VV, NEV, VUF, PVUR) rated into
    per-decision disagreement `η_D = 1 − ACC_D` and its worst case
    `η_∞ = max_D η_D`, plus the diagnostic `MAPE_U`.

The two axes combine into the joint verdict `:verified`, `:tolerable`, or
`:consequential`.

The statistical axis is estimator-agnostic: the M-estimator-specific
normalisation constant `C_ρ` is identical for `M_r` and `M_c` and cancels in
the paired fidelity gap (companion paper, "State Estimation as the
Likelihood Generator"), so WLS, WLAV, rWLAV, etc. share one code path.

Public API:
- `verify(se_r, math_r, se_c, math_c; kwargs...)` - single realisation
- `verify_mc(se_r_vec, math_r, se_c_vec, math_c; kwargs...)` - Monte-Carlo
- `fidelity_gap`, `equivalence_test`, `joint_verdict`
- `operational_axis`, `decision_disagreement`, `voltage_mape`
- `voltage_violations`, `nev_alarms`, `vuf_alarms`, `pvur_alarms`
- `J_star`, `n_state_variables`, `n_measurements`, `is_krn`, `bus_voltages`, `vbase_V`
"""
module DTVerify

using Statistics
using LinearAlgebra
using Distributions
using DataFrames

include("extract.jl")
include("tests.jl")
include("decisions.jl")
include("api.jl")

# Extraction
export J_star,
       n_state_variables,
       n_measurements,
       bus_voltages,
       vbase_V,
       is_krn,
       flatten_phase_voltages,
       flatten_neutral_voltages

# Statistical axis
export fidelity_gap,
       equivalence_test,
       joint_verdict

# Operational axis / decisions
export voltage_violations,
       nev_alarms,
       vuf_alarms,
       pvur_alarms,
       compare_decisions,
       compare_decisions_vv,
       voltage_mape,
       decision_disagreement,
       operational_axis,
       pool_confusion

# Main user-facing
export verify, verify_mc

end # module DTVerify
