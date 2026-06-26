#!/usr/bin/env python3
"""
DICOM CT + RTSTRUCT → matRad MAT file converter.

Directory layout (relative to this script):
    ../data/input/<patient_id>/   ← DICOM files (.dcm)
    ../data/output/<patient_id>/  ← output MAT file

Usage:
    python 01_dicom2mat.py                  # convert all patients
    python 01_dicom2mat.py 1ABB027          # convert one patient
"""

import argparse
import subprocess
import tempfile
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR  = Path(__file__).parent
DICOM_ROOT  = SCRIPT_DIR.parent / "data" / "input"   # one sub-folder per patient
OUTPUT_ROOT = SCRIPT_DIR.parent / "data" / "output"  # mirrored structure

# ---------------------------------------------------------------------------
# External tools
# ---------------------------------------------------------------------------
MATLAB_BIN  = "/usr/local/MATLAB/link/R2021a/matlab"
MATRAD_PATH = SCRIPT_DIR.parent / "matRad-dev_VMAT_merge"


# ---------------------------------------------------------------------------
# MATLAB script template
# ---------------------------------------------------------------------------

def build_matlab_script(dicom_dir: Path, output_mat: Path) -> str:
    """Return a self-contained MATLAB script that imports DICOM and saves MAT."""
    return f"""\
try
    % ---- setup ----
    addpath(genpath('{MATRAD_PATH}'));
    cfg = MatRad_Config.instance();
    cfg.disableGUI = true;
    cfg.logLevel   = 3;

    % ---- import ----
    importer = matRad_DicomImporter('{dicom_dir}');
    if size(importer.importFiles.ct, 1) == 0
        error('No CT slices found in: {dicom_dir}');
    end
    importer.matRad_importDicom();

    ct  = evalin('base', 'ct');
    cst = evalin('base', 'cst');

    % ---- save ----
    save('{output_mat}', 'ct', 'cst', '-v7.3');
    fprintf('Saved → %s\\n', '{output_mat}');

catch err
    fprintf('ERROR: %s\\n', err.message);
    exit(1);
end
"""


# ---------------------------------------------------------------------------
# Run helpers
# ---------------------------------------------------------------------------

def run_matlab(script: str) -> bool:
    """Write *script* to a temp .m file and execute it with MATLAB."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.m', delete=False) as fh:
        fh.write(script)
        tmp = Path(fh.name)

    try:
        cmd = [
            MATLAB_BIN,
            '-nodisplay', '-nosplash', '-nodesktop',
            '-r', f"cd('{tmp.parent}'); {tmp.stem}; exit;",
        ]
        return subprocess.run(cmd, text=True, timeout=600).returncode == 0

    except subprocess.TimeoutExpired:
        print("ERROR: MATLAB timed out (>10 min)")
        return False
    except Exception as exc:
        print(f"ERROR: could not launch MATLAB — {exc}")
        return False
    finally:
        tmp.unlink(missing_ok=True)


def convert_patient(patient_id: str) -> bool:
    """Convert one patient's DICOM folder to a MAT file."""
    dicom_dir  = DICOM_ROOT  / patient_id
    output_mat = OUTPUT_ROOT / patient_id / f"{patient_id}.mat"

    # --- validate input ---
    if not dicom_dir.exists():
        print(f"  [SKIP] Input folder not found: {dicom_dir}")
        return False
    if not list(dicom_dir.glob("*.dcm")):
        print(f"  [SKIP] No .dcm files in: {dicom_dir}")
        return False

    # --- convert ---
    output_mat.parent.mkdir(parents=True, exist_ok=True)
    print(f"  {dicom_dir}  →  {output_mat}")
    return run_matlab(build_matlab_script(dicom_dir, output_mat))


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> bool:
    parser = argparse.ArgumentParser(
        description="Convert DICOM (CT + RTSTRUCT) to matRad MAT format."
    )
    parser.add_argument(
        "patient", nargs="?", default=None,
        metavar="PATIENT_ID",
        help="Patient folder name inside data/input/ (omit to process all)",
    )
    args = parser.parse_args()

    # Collect patients to process
    if args.patient:
        patients = [args.patient]
    else:
        patients = sorted(p.name for p in DICOM_ROOT.iterdir() if p.is_dir())

    print(f"Processing {len(patients)} patient(s)...\n")

    succeeded, failed = 0, []
    for pid in patients:
        print(f"[{pid}]")
        if convert_patient(pid):
            succeeded += 1
        else:
            failed.append(pid)

    print(f"\n{'─'*50}")
    print(f"Done:  {succeeded} succeeded,  {len(failed)} failed")
    if failed:
        print("Failed:", ", ".join(failed))

    return len(failed) == 0


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
