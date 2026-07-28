# DTVerify.jl

Decision-aware **Digital Twin Verification (DTV)** for distribution-system
state estimation (DSSE). Reference implementation of the framework in the
companion paper (Section 2, *The Digital Twin Verification Framework*),
released under the MIT license together with the cable and network data of
the paper's case studies (see [`data/`](data/README.md)).

Given paired **reference** (`M_r`) and **candidate** (`M_c`) DSSE result dicts
and their math dicts, `DTVerify` scores each comparison on two axes and returns
a joint verdict:

- **Statistical axis** - the per-measurement fidelity gap `d̄` and a one-sided
  paired-*t* **equivalence test** at declared knobs `(δ_eq, α_eq)`.
- **Operational axis** - four standard-based alarms (VV, NEV, VUF, PVUR) rated
  into per-decision disagreement `η_D = 1 − ACC_D` and its worst case
  `η_∞ = max_D η_D`, plus the diagnostic `MAPE_U`.

The two axes combine into the verdict `:verified`, `:tolerable`, or
`:consequential`.

The statistical axis is **estimator-agnostic**: the M-estimator-specific
normalisation constant cancels in the paired fidelity gap, so WLS, WLAV, rWLAV,
etc. share one code path.

## Install

Not yet registered. Clone, then either `dev` it into another project:

```julia
pkg> dev /path/to/DTVerify.jl
```

or work in the package directory directly:

```julia
pkg> activate .
pkg> instantiate
```

## Quickstart

```julia
using DTVerify

# se_r_vec[j], se_c_vec[j]   paired reference/candidate DSSE result dicts for
# Monte-Carlo realisation j (output of solve_mc_se + dictify_solution!).
# math_r, math_c   the corresponding math dicts (output of transform_data_model).
result = verify_mc(se_r_vec, math_r, se_c_vec, math_c;
                   δ_eq = 0.05, α_eq = 0.05, α = 0.10)

result.verdict            # :verified | :tolerable | :consequential
result.d_bar              # per-measurement fidelity gap  d̄
result.s_d                # across-realisation spread of the gap
result.equivalence_holds  # one-sided equivalence test outcome (Bool)
result.eq_bound           # upper (1−α_eq) confidence bound on d̄
result.eta_inf            # worst-alarm operational disagreement  η_∞
result.eta                # (vv, nev, vuf, pvur)  per-decision η_D = 1 − ACC_D
result.acc, result.rec    # (vv, nev, vuf, pvur)  accuracy / recall
result.mape_u             # phase-voltage-magnitude MAPE diagnostic (%)
result.maxape_u           # worst-bus APE (%)
result.counts             # (vv, nev, vuf, pvur, total)  disagreement counts (FP+FN)
result.n_realizations     # n
result.k_r, result.k_c    # state-variable counts
```

`verify(se_r, math_r, se_c, math_c; …)` runs a single realisation and returns
the same schema. A single realisation carries no variance, so the equivalence
test cannot run (`equivalence_holds = false`) and the verdict is therefore
never `:verified` - use `verify_mc` with `n ≥ 2` for a full statistical verdict.

Pass `M = <n_measurements>` explicitly if the math dict carries no `meas` block
(otherwise `|M|` defaults to `n_measurements(math_r)`).

## Interpreting the verdict

The joint verdict reads the two axes together:

| Verdict           | Condition                              | Meaning |
|-------------------|----------------------------------------|---------|
| `:verified`       | `η_∞ ≤ α` **and** equivalence holds    | Statistically equivalent to the reference within the margin, and decision-preserving. |
| `:tolerable`      | `η_∞ ≤ α` **and** equivalence fails    | Statistically distinguishable, but does not change any operational decision beyond tolerance. |
| `:consequential`  | `η_∞ > α`                              | Changes operator decisions beyond the operational tolerance. |

Separating *detectable* from *decision-changing* model error is the central
function of the joint verdict. The declared knobs are deployment parameters;
the paper's results use `δ_eq = 0.05`, `α_eq = 0.05` (one-sided, 95 %
confidence), and `α = 0.10`.

## Building blocks

The public API decomposes into reusable pieces:

- `fidelity_gap(J_r, J_c, M) → (d_bar, s_d, d_j)`
- `equivalence_test(d_bar, s_d, n; δ_eq, α_eq) → (holds, bound, …)`
- `joint_verdict(equivalence_holds, η_inf; α) → Symbol`
- `voltage_violations`, `nev_alarms`, `vuf_alarms`, `pvur_alarms` - the alarms
- `decision_disagreement`, `operational_axis` - confusion counts → `η_D`, `η_∞`
- `voltage_mape` - the `MAPE_U` diagnostic
- `J_star`, `n_state_variables`, `n_measurements`, `is_krn`, `bus_voltages`,
  `vbase_V` - PMDSE dict accessors

See `examples/kron_reduction_case.jl` for the paper-reproduction example
and `test/` for synthetic-data unit tests.

## Data

The paper's cable and network data ship in [`data/`](data/README.md):

- `data/network_1/Feeder_1/` - ENWL LV test Network 1, Feeder 1 (OpenDSS),
  with the finite-element linecodes of the trunk and lateral cables.
- `data/zprim/` - primitive series-impedance matrices (JLD2) for the
  NAYCWY 3×95 / NAYY cables in both finite-element (`fem`) and analytical
  (`full_carson`, `simple_carson`, `deri`) formulations.
- `data/kron_case/` - archived per-realisation objectives and pooled alarm
  counts of the paper's n = 30 Kron-reduction run (Table 1, case 1).

## Reproducing the paper's Kron-reduction verdict

The single example replays the paper's headline case from the archived
n = 30 run and reproduces Table 1, case 1 exactly (fidelity gap
d̄ = 3.304, worst-alarm disagreement η∞ = 0.62 through the
neutral-to-earth-voltage alarm, verdict Consequential):

```
cd DTVerify.jl/examples
julia --project=. -e 'using Pkg; Pkg.instantiate()'   # first time only
julia --project=. kron_reduction_case.jl
```

It needs only this package plus CSV/DataFrames (Julia ≥ 1.11), and runs
in seconds.

## Reference

The framework, notation, and equations are defined in the companion paper,
Section 2 (*The Digital Twin Verification Framework*): the fidelity gap and
equivalence test in *Verification as a Likelihood Equivalence Test*, the alarm
disagreement in *Operational Decision Agreement*, and the three-region verdict
in *Joint Verdict*.

## License

MIT - see [LICENSE](LICENSE).
