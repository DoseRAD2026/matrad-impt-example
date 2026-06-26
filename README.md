# matrad-impt-example

IMPT treatment planning pipeline for thorax (`1TH*`) and abdomen (`1AB*`) patients: DICOM → matRad MAT → IMPT optimization (PTV/geometry derived from the CT → fluence optimization → spot filtering).

---

## Requirements

- Python 3.8+ (standard library only — no pip packages required)
- MATLAB R2021a
- matRad (see [Setup](#setup) below)

---

## Setup

**1. Clone with submodule**

```bash
git clone --recurse-submodules https://github.com/DoseRAD2026/matrad-impt-example
```

Or if already cloned:

```bash
git submodule update --init
```

**2. Configure MATLAB path**

Both scripts hardcode the MATLAB binary path. Edit if yours differs:

```python
# script/01_dicom2mat.py and script/02_proton_optimization.py
MATLAB_BIN = "/usr/local/MATLAB/link/R2021a/matlab"
```

**3. Place DICOM input data**

```
data/input/<patient_id>/    ← CT + RTSTRUCT .dcm files
```

Patient IDs must contain `TH` (thorax/lung) or `AB` (abdomen) — the prescription is chosen accordingly.

---

## Structure

```
matrad-impt-example/
├── data/
│   ├── input/<patient_id>/       # CT + RTSTRUCT DICOM files
│   └── output/<patient_id>/      # all outputs (generated)
├── matRad-dev_VMAT_merge/        # matRad submodule (pinned commit)
└── script/
    ├── 01_dicom2mat.py
    └── 02_proton_optimization.py
```

---

## Usage

```bash
# Step 1 — DICOM → MAT
python script/01_dicom2mat.py                   # all patients
python script/01_dicom2mat.py 1ABB027           # single patient

# Step 2 — IMPT optimization
python script/02_proton_optimization.py             # all patients
python script/02_proton_optimization.py 1ABB027     # single patient
python script/02_proton_optimization.py 1ABB027 --force-recalc-dij   # skip dij cache
python script/02_proton_optimization.py 1ABB027 -c pln.json          # custom plan params
```

---

## Outputs

Written to `data/output/<patient_id>/`:

| File | Content |
|------|---------|
| `<id>.mat` | matRad `ct` + `cst` structs |
| `<id>_proton_stf_dij_cache.mat` | Beam geometry + dose influence matrix (cached) |
| `<id>_proton.mat` | Full result (`resultGUI`, `pln`, `stf`, `cst`, `ct`, opt history) |
| `<id>_proton.json` | Per-spot export (beam → ray → beamlet energies & weights) |
| `<id>_proton_dose.mha` | Optimized RBE-weighted dose cube, per-fraction (single precision) |
| `<id>_proton_dvh_report.txt` | DVH metrics + constraint PASS/FAIL |
| `<id>_proton_dvh.png` | DVH plot |

> The stf/dij calculation is slow and cached automatically. Re-runs reuse the cache unless parameters change or `--force-recalc-dij` is passed.

---

## Planning Protocol

**Prescription** (by site, inferred from patient ID):

| Site | ID pattern | Prescription |
|------|-----------|--------------|
| Lung | `1TH*` | 70 Gy(RBE) / 35 fx |
| Abdomen | `1AB*` | 60 Gy(RBE) / 30 fx |

| Parameter | Value |
|-----------|-------|
| Radiation mode | protons |
| Biological model | constRBE (RBE = 1.1), optimize on RBExDose |
| Dose engine | HongPB (analytical pencil beam) |
| Dose grid | 3 mm |
| Optimizer | IPOPT |
| Spot filter | drop spots with weight < 0.1 |

**Self-contained geometry derivation** (no external plan files needed):

- **PTV** — selected automatically by priority: `PTV` > `PTVOPT` > `PTVUNION` > `PTV1` > `PTV_HAUPTPLAN` > any remaining `PTV*`.

> **On the optimization target.** The structure sets follow conventional (photon-based) clinical delineation. `PTV_OPT` is an auxiliary structure created for *photon* plan optimization (e.g. cropped from the skin/OARs), not a proton-specific target. This pipeline therefore optimizes directly on the nominal clinical `PTV` — the prescription/evaluation target — and falls back to `PTV_OPT` (or the other variants) only when no nominal `PTV` exists. No robust optimization is applied.
- **Spot grid** — adapts to PTV volume: < 15 cc → bixel 3 mm / spacing 3 mm; ≥ 300 cc → bixel 8 mm / spacing 7 mm; otherwise 6 mm / 5 mm.
- **Gantry geometry** — an ipsilateral 3-field template chosen by PTV laterality relative to the body midline: PTV right → `[210°, 270°, 320°]`, PTV left → `[40°, 90°, 150°]` (the two are mirror images about the AP axis), couch 0°.

> The gantry selection is a laterality-based heuristic that reproduces the clinical beam arrangement in most cases but is not a substitute for case-by-case clinical beam-angle optimization.

A ring structure (PTV+5 mm to PTV+12 mm, clipped to body, GTV-excluded) is added to control dose fall-off.

| OAR (lung) | Constraint |
|-----|-----------|
| Spinal cord | Dmax < 50.5 Gy |
| Lung | V60Gy < 33%, V5Gy < 33% |
| Oesophagus | Dmax < 73.5 Gy, V60Gy < 15.3% |
| Brachial plexus | Dmax < 66 Gy |

| OAR (abdomen) | Constraint |
|-----|-----------|
| Spinal cord | Dmax < 50.5 Gy |
| Liver | Dmean < 30 Gy |
| Kidney | Dmean < 18 Gy |
| Stomach | Dmax < 45 Gy |
| Bladder | Dmax < 65 Gy |
