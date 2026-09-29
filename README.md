# UV/AOP Kintecus Pipeline — 1,4-Dioxane Degradation

An R harness around the [Kintecus](http://www.kintecus.com/) chemical-kinetics
engine for simulating 1,4-dioxane destruction in UV/H₂O₂/chloramine advanced
oxidation. It models 254 nm photolysis with time-varying irradiance,
breakpoint-chlorination chemistry with concentration-dependent rate laws, and
~150 radical reactions, across a batch of experimental conditions read from a
spreadsheet.

Kintecus itself integrates a fixed reaction network with fixed rate constants.
Several rate constants in this system are not fixed — they depend on pH,
carbonate, hypochlorite, and on how much UV is currently reaching the water,
all of which change as the reaction proceeds. This harness solves that by
**operator splitting**: it runs Kintecus in 1-second substeps, recomputing
every concentration-dependent rate constant between them.

---

## Contents

| File | Role |
|---|---|
| `run_model.R` | Batch driver. Seeds concentrations, runs the substep loop, writes per-experiment CSVs. |
| `photolysis.R` | All chemistry constants and helper functions: speciation, photolysis rates, breakpoint rate laws. |
| `MODEL_template.DAT` | Kintecus reaction network, with `_K_*_` placeholders for the dynamic rate constants. |
| `SPECIES_template.DAT` | Bare species registry — one name per line. Nothing else. |
| `PARM_template.DAT` | Kintecus run parameters, with `_HRT_`, `_SAMP_INT_`, `_TEMP_K_` placeholders. |
| `Orchestra.qmd` | Quarto report: runs the batch, summarizes, plots, compares to the target model. |
| `xlsx_parsing.R` | One-time converter: source reaction spreadsheet → `MODEL.DAT`. Not part of the run loop. |
| `target_model.xlsx` | Reference trajectory from the model being validated against. |
---

## Requirements

- **Windows.** The driver shells out to `C:/KINTECUS/kintecus.exe` and calls
  `taskkill`. Paths are hardcoded at the top of `run_model.R`.
- **Kintecus** installed at `C:/KINTECUS/`, with a valid `KINKEY.TXT`. The
  driver copies the key into the working directory automatically on first run.
- **R** with `readxl`. `Orchestra.qmd` additionally needs `dplyr`, `ggplot2`,
  `purrr`, and Quarto.

---

## Running

Put `kintecus_experimental_batch.xlsx` in the working directory alongside the
scripts, then either:

```r
source("run_model.R")          # batch only
```

or render the full report:

```bash
quarto render Orchestra.qmd    # batch + summary + plots + target comparison
```

Each experiment row produces `CONC_Exp_<Exp_ID>.csv`. `Orchestra.qmd` also
writes `Master_Degradation_Summary.csv`.

---

## Input workbook schema

One row per experiment. Missing columns silently default to 0, so check
spelling — a typo'd column name is not an error, it's a zero dose.

### Reactor and UV

| Column | Units | Notes |
|---|---|---|
| `Exp_ID` | — | Used in the output filename |
| `HRT` | s | Total simulated time. Defaults to 1. |
| `E_avg` | mW/cm² | Incident average irradiance at 254 nm |
| `path_length` | cm | Optical path. Must be > 0. |
| `pH` | — | Held constant for the whole run (see *pH pinning*) |
| `UVT` or `UVT_pct` | % | 1 cm transmittance at 254 nm. Defaults to 100. |
| `Temp_K` **or** `Temp_C` | K / °C | `Temp_K` wins if both present. Defaults to 298.15 K. |

### Doses

Every dose column is in mass units and converted internally.

| Column | Units | Converted with |
|---|---|---|
| `H2O2_mgL` | mg/L as H₂O₂ | 34.0147 g/mol |
| `DIOXANE_ugL` | **µg/L** as C₄H₈O₂ | 88.106 g/mol; also the LRV baseline |
| `FAC_mgL` | mg/L as Cl₂ | 70.90 g/mol, then split HOCl / OCl⁻ by pKa 7.53 |
| `NH2Cl_mgL` | mg/L as Cl₂ | 70.90 g/mol (one Cl per molecule) |
| `NHCl2_mgL` | mg/L as Cl₂ | **141.80 g/mol** (two Cl per molecule) |
| `NH3_mgL_N` | mg/L as N | 14.01 g/mol, then split NH₃ / NH₄⁺ by pKa 9.30 |
| `Alk_mgL` | mg/L as CaCO₃ | 100.09 g/mol, ÷2 for equivalents → full carbonate speciation |
| `Cl_mgL` | mg/L as Cl⁻ | 35.453 g/mol |
| `TOC_mgL` | mg/L as C | 12.011 g/mol |
| `O2_mgL` | mg/L | 31.998 g/mol |
| `NO3_mgL` | mg/L **as N** | 14.01 g/mol |
| `NO2_mgL` | mg/L **as N** | 14.01 g/mol |

Note that dioxane is **µg/L**, not mg/L — it's the one column with a 10⁶
divisor. The nitrogen species are **as N**, not as the ion.

#### On "as Cl₂"

Online chloramine analyzers report oxidizing equivalents, not molecular mass.
Each chlorine atom counts as one Cl₂ equivalent, so the divisor is
70.90 × (number of Cl atoms):

- NH₂Cl → 1 Cl → **70.90 g/mol**
- NHCl₂ → 2 Cl → **141.80 g/mol**

`1.62 mg/L NHCl₂ ÷ 1000 ÷ 141.80 = 1.14 × 10⁻⁵ M`

---

## How it works

### The substep loop

For each experiment:

1. Seed all "n" species. Everything starts at zero, then the dosed species are
   filled in from the speciation helpers.
2. Compute fixed background absorbance `A_bg` from measured UVT, minus what
   the model's own photolyzing species absorb at t = 0.
3. Loop, in `DT = 1.0 s` substeps:
   - `k_photo()` — photolysis rate constants from the **current**
     concentrations
   - `k_breakpoint()` — breakpoint rate constants from the current
     concentrations
   - Write `PARM.DAT`, `MODEL.DAT` (placeholder substitution) and `SPECIES.DAT`
     (written from scratch)
   - Run `kintecus.exe -ig:mass`, read `CONC.TXT`
   - Carry **every** species forward into the next substep
4. Concatenate, compute fluence and LRV, write the CSV.

Kintecus samples at `dt/10`, so each substep contributes 10 output rows. The
duplicate row at each substep boundary is dropped.

### pH pinning

`PINNED <- c("H+", "OH-")` writes a concentration into Kintecus's
`Constant File?` column for those two species, which holds them fixed for the
entire run. pH is therefore an input, not a result. `OH⁻` is derived as
`10^(pH − pKw)`.

This is important for breakpoint reactions where `k207`, `k211`, `k212`, `k213` and `k215` are
all functions of `[OH⁻]`, so pinning makes those rate constants stable.

### Photolysis

Optically-dilute (Taylor-expanded) Beer-Lambert:

```
k = ln(10) · QY · I₀ · ε · L
I₀ = E_avg · I_ratio / (470697.6 · L)
```

`470697.6` is the molar photon energy at 254 nm in J/einstein; the mW→W and
mL→L factors cancel, so `I₀` comes out in einstein/(L·s). Note that **path
length cancels** between the two expressions — `L` reaches the rate constant
only through `I_ratio`.

At 254 nm:

| Species | ε (L/mol/cm) | Quantum yield |
|---|---|---|
| H₂O₂ | 18.6 | 0.50 |
| HOCl | 58.0 | 0.55 |
| OCl⁻ | 60.0 | 0.55 |
| NH₂Cl | 371.0 | 0.294 |
| NHCl₂ | 126.0 | 0.82 |

### The irradiance ratio

As absorbers photolyze away the water clears and the cell brightens, so the
rate constants rise. `I_ratio` tracks that, built from the **path-averaged
fluence rate** — what a well-mixed cell actually experiences — rather than the
intensity at the exit plane:

```
Ebar(A, L) = (1 − 10^(−A·L)) / (A·L)
I_ratio    = Ebar(A_sp, L) / Ebar(A_sp0, L)
```

- **Background absorbance is excluded from this ratio on purpose.** Only the
  five photolyzing species enter it. Under the old exit-plane form `A_bg`
  cancelled algebraically; under the averaging form it would not, so it is now
  kept out explicitly and used only for the reported `UVT_pct` / `A_path`
  columns. This matches the reference model. It is a simplification: with
  substantial real background absorbance, both models slightly overstate how
  much the cell brightens.

### Fluence

```r
run$Fluence_mJ_cm2 <- E_avg * run$I_ratio * run[[tcol]]
```

A **product**, not a cumulative integral — current irradiance × elapsed time.
This matches the reference model, whose dose column divided by `E_avg · t`
equals its `I ratio` on every row. 

### Breakpoint chemistry

`k_breakpoint()` returns the concentration- and temperature-dependent rate
constants substituted into `_K_201_` … `_K_215_`. Formula numbers follow the
source spreadsheet. Arrhenius terms use `T_K`; the rest are functions of the
current `[H⁺]`, `[OH⁻]`, `[HCO₃⁻]`, `[H₂CO₃]`, `[CO₃²⁻]`, `[HOCl]`, `[OCl⁻]`.

`k214` deliberately omits `[NHCl₂]` from the rate expression because Kintecus
supplies it as the reaction's own reactant.

### Carbonate speciation

```
Alk_eq = Ct·(α₁ + 2α₂) + [OH⁻] − [H⁺]
Ct     = (Alk_eq + [H⁺] − [OH⁻]) / (α₁ + 2α₂)
```

The `+[H⁺]` proton correction is not optional. At pH 5.6 with low-alkalinity RO
permeate it is ~4% of the alkalinity; dropping it put total carbonate 4.3%
below the reference. With it, t₀ carbonate matches to <0.1%.

Note that `pK1 = 6.32` here is for H₂CO₃\* (CO₂(aq) + H₂CO₃), while the
`H2CO3 ==> HCO3- + H+` pair in `MODEL_template.DAT` (5.00E+05 against
1.00E+10) is a kinetic pKa of 4.30, correct for *true* H₂CO₃. Both are
intentional and match the reference model, which re-equilibrates to
H₂CO₃/HCO₃⁻ = 0.048 within the first timestep while conserving total carbonate.

---

## Output schema

`CONC_Exp_<id>.csv` — one row per Kintecus sample point.

| Column(s) | Meaning |
|---|---|
| `Time(s)` | Cumulative time across substeps |
| 62 species columns | Molar concentrations |
| `I_ratio` | Path-averaged irradiance relative to t₀ |
| `E_now` | `E_avg × I_ratio`, mW/cm² |
| `A_path` | Total absorbance across the path (background included) |
| `UVT_pct` | `100 × 10^(−A_path)` |
| `k_photo_*` | Photolysis rate constants for the five species, s⁻¹ |
| `Fluence_mJ_cm2` | `E_avg × I_ratio × t` |
| `LRV` | `log10(DIOXANE_0 / DIOXANE)` |

`I_ratio`, `E_now`, `A_path` and `UVT_pct` are recomputed from each row's own
concentrations, so they track within a substep.

`k_photo_*` are the rate constants Kintecus **actually integrated** over that
substep, and are constant across its 10 rows. That is the operator splitting,
not a reporting shortcut: the solver genuinely used one value for the whole
second. They therefore lag `I_ratio` slightly within a substep. Shrink `DT` to
reduce the lag; it cannot be post-corrected.

---
---

## Known issues and open items

- **Formula 207 is wired but inert.** `k207 = 167·[OH⁻]` is computed in
  `photolysis.R` and substituted into `_K_207_`, but the reaction line in
  `MODEL_template.DAT` is still commented out:
  `#_K_207_  NHCl2 ==> NOH + 2H+ + 2Cl-`. Substitution rewrites the comment,
  so the reaction never fires. Uncomment the line to activate it.
- **Formula 214 stoichiometry is unresolved.** The rate law uses `[OCl⁻]` but
  the written reaction shows `2HOCl`. Flagged in `MODEL_template.DAT`; confirm
  against the source spreadsheet before trusting that pathway.
- **Reverse disproportionation** (`NHCl2 + NH3 ==> 2NH2Cl`, 6.11E+04) has no
  formula in the source and is kept as a static literature value.
- **`pKw = 14.0`** here; the reference model uses 13.967 (Kw = 1.079 × 10⁻¹⁴).
  Change in `CONST` if exact agreement is wanted.
- **Windows-only.** `C:/KINTECUS/` paths and `taskkill` are hardcoded.

---

## Extending the model

**To add a species:** add its name to `SPECIES_template.DAT` and use it in
`MODEL_template.DAT`. Nothing else — the driver validates that every seeded
species exists in the registry and errors loudly if not.

**To add a reaction with a fixed rate constant:** add a
`<rate>\t<reactants> ==> <products>` line to `MODEL_template.DAT`.

**To add a reaction with a dynamic rate constant:** use a `_K_<name>_`
placeholder in `MODEL_template.DAT` and return `k<name>` from
`k_breakpoint()`. The name-based mapping picks it up automatically.

**To add a photolyzing species:** add a row to `PHOTO` in `photolysis.R` with
its species name, its `_K_*_` placeholder, ε and quantum yield, and put the
matching placeholder on a reaction line in `MODEL_template.DAT`. The species
then contributes to the absorbance and to `I_ratio` automatically. Row order
does not matter — the placeholder travels on the same row as the species name,
and the driver pairs them by name.

**To change a constant:** everything tunable lives in `CONST` and `PHOTO` at
the top of `photolysis.R`. Nothing else in the pipeline hardcodes one.
