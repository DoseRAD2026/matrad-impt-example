#!/usr/bin/env python3
"""
IMPT treatment planning & optimization using matRad (ipsilateral 3-field, self-contained geometry).

Directory layout (relative to this script):
    ../data/output/<patient_id>/<patient_id>.mat  ← input from 01_dicom2mat.py
    ../data/output/<patient_id>/                  ← MHA, JSON, DVH report, PNG outputs

Usage:
    python 02_proton_optimization.py                  # optimize all patients
    python 02_proton_optimization.py 1ABB027          # optimize one patient
    python 02_proton_optimization.py 1ABB027 -t 3600  # custom timeout (s)
    python 02_proton_optimization.py 1ABB027 -c pln.json  # custom plan params
"""

import argparse
import json
import subprocess
import tempfile
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR  = Path(__file__).parent
MAT_ROOT    = SCRIPT_DIR.parent / "data" / "output"

# ---------------------------------------------------------------------------
# External tools
# ---------------------------------------------------------------------------
MATLAB_BIN  = "/usr/local/MATLAB/link/R2021a/matlab"
MATRAD_PATH = SCRIPT_DIR.parent / "matRad-dev_VMAT_merge"

# ---------------------------------------------------------------------------
# Default proton planning parameters
# ---------------------------------------------------------------------------
DEFAULT_PLN = {
    "radiationMode":            "protons",
    "machine":                  "Generic",
    "bioModel":                 "constRBE",
    "multScen":                 "nomScen",
    # Beam geometry (default spot grid; adjusted by PTV volume, see below)
    "bixelWidth":               6,
    "longitudinalSpotSpacing":  5,
    # Dose calculation
    "doseGridResolution":       3,
    "doseEngine":               "HongPB",
    "calcLET":                  0,
    # Optimization
    "quantityOpt":              "RBExDose",
    "optimizer":                "IPOPT",
    "maxIter":                  500,
    "minAbsWeight":             0.1,   # drop spots with w < 0.1 (absolute threshold)
}

# ---------------------------------------------------------------------------
# PTV-volume-dependent spot grid (cc thresholds -> [bixelWidth, longSpacing] mm)
# ---------------------------------------------------------------------------
SMALL_PTV_CC,  SMALL_BIXEL,  SMALL_SPACING  = 15.0,  3, 3   # small  PTV -> finer grid
XLARGE_PTV_CC, XLARGE_BIXEL, XLARGE_SPACING = 300.0, 8, 7   # large  PTV -> coarser grid

# ---------------------------------------------------------------------------
# Ipsilateral 3-field gantry templates (mirror images about the AP axis)
# ---------------------------------------------------------------------------
GANTRY_PTV_RIGHT = [210, 270, 320]   # PTV on patient's right → beams from the right
GANTRY_PTV_LEFT  = [40, 90, 150]     # PTV on patient's left  → beams from the left


# ---------------------------------------------------------------------------
# MATLAB script template
# ---------------------------------------------------------------------------

def build_matlab_script(
    patient_id: str,
    mat_file: Path,
    pln_params: dict,
    force_recalc_dij: bool,
) -> str:
    """Return a self-contained MATLAB script for IMPT optimization."""
    pln = {**DEFAULT_PLN, **pln_params}

    return f"""\
try
    % ---- setup ----
    addpath(genpath('{MATRAD_PATH}'));
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.disableGUI = true;
    matRad_cfg.logLevel = 1;
    matRad_rc;

    passFail = {{'(FAIL)', '(PASS)'}};

    % ---- load patient data ----
    load('{mat_file}');

    % ---- output paths ----
    outDir      = '{mat_file.parent}';
    resultMat   = fullfile(outDir, '{patient_id}_proton.mat');
    spotJson    = fullfile(outDir, '{patient_id}_proton.json');
    doseMhaPath = fullfile(outDir, '{patient_id}_proton_dose.mha');
    reportPath  = fullfile(outDir, '{patient_id}_proton_dvh_report.txt');
    dvhPngPath  = fullfile(outDir, '{patient_id}_proton_dvh.png');

    patientID = '{patient_id}';
    if contains(patientID, 'TH')
        patientType = 'lung';
        prescriptionDose = 70;
        numFractions = 35;
    elseif contains(patientID, 'AB')
        patientType = 'abdomen';
        prescriptionDose = 60;
        numFractions = 30;
    else
        error('Unknown patient type: %s', patientID);
    end
    fprintf('Patient: %s (%s) %dGy(RBE)/%dfx\\n', patientID, patientType, prescriptionDose, numFractions);

    cleanName = @(s) lower(regexprep(s, '[_ ]s2$', ''));

    % ---- select TARGET (priority: PTV > PTVOPT > PTVUNION > PTV1 > PTV_HAUPTPLAN > any PTV*) ----
    ptvPatterns = {{
        '^ptv$',
        'ptv[ _-]*opt',
        '^ptv[ _-]*union$|^ptvunion$',
        '^ptv1$',
        '^ptv[ _-]*hauptplan$|ptv[ _-]*hauptplan',
        '^ptv'
    }};
    selectedTargetIdx = 0;
    for p = 1:numel(ptvPatterns)
        for i = 1:size(cst, 1)
            if isequal(cst{{i,3}},'TARGET') && ~isempty(regexpi(cleanName(cst{{i,2}}), ptvPatterns{{p}}))
                selectedTargetIdx = i; break;
            end
        end
        if selectedTargetIdx > 0, break; end
    end
    if selectedTargetIdx == 0
        error('No suitable PTV found in cst for %s', patientID);
    end
    selectedPTVName = cst{{selectedTargetIdx, 2}};
    fprintf('Selected TARGET: %s (index %d)\\n', selectedPTVName, selectedTargetIdx);

    % Lung structure indices (for lung OAR objectives)
    lungStructureIndices = [];
    if strcmp(patientType, 'lung')
        lungTotalIdx = 0; lungLeftIdx = 0; lungRightIdx = 0;
        for i = 1:size(cst, 1)
            name = cleanName(cst{{i, 2}});
            if ~isempty(regexpi(name, '^lung_?total$|^lungtotal$'))
                lungTotalIdx = i;
            elseif ~isempty(regexpi(name, '^lung_?left$|^lungleft$'))
                lungLeftIdx = i;
            elseif ~isempty(regexpi(name, '^lung_?right$|^lungright$'))
                lungRightIdx = i;
            end
        end
        if lungTotalIdx > 0
            lungStructureIndices = lungTotalIdx;
        elseif lungLeftIdx > 0 || lungRightIdx > 0
            lungStructureIndices = [lungLeftIdx, lungRightIdx];
            lungStructureIndices(lungStructureIndices == 0) = [];
        end
    end

    % ---- PTV mask, volume and body midline ----
    if isfield(ct, 'cubeHU') && ~isempty(ct.cubeHU)
        baseCube = ct.cubeHU; if iscell(baseCube); baseCube = baseCube{{1}}; end
    elseif isfield(ct, 'cube') && ~isempty(ct.cube)
        baseCube = ct.cube; if iscell(baseCube); baseCube = baseCube{{1}}; end
    else
        cubeDim = round(ct.cubeDim(:)');
        baseCube = zeros(cubeDim(1), cubeDim(2), cubeDim(3));
    end
    ptvMask = false(size(baseCube));
    ptvIdx  = unique(vertcat(cst{{selectedTargetIdx, 4}}{{:}}));
    ptvIdx  = ptvIdx(ptvIdx >= 1 & ptvIdx <= numel(ptvMask));
    if isempty(ptvIdx); error('Selected PTV has no valid voxels'); end
    ptvMask(ptvIdx) = true;

    % Spot grid: finer for small PTVs, coarser for very large PTVs
    ptvVolumeCC = numel(ptvIdx) * ct.resolution.x * ct.resolution.y * ct.resolution.z / 1000;
    bixelWidth = {pln["bixelWidth"]}; longSpacing = {pln["longitudinalSpotSpacing"]};
    if ptvVolumeCC < {SMALL_PTV_CC}
        bixelWidth = {SMALL_BIXEL}; longSpacing = {SMALL_SPACING};
    elseif ptvVolumeCC >= {XLARGE_PTV_CC}
        bixelWidth = {XLARGE_BIXEL}; longSpacing = {XLARGE_SPACING};
    end
    fprintf('PTV volume: %.2f cc -> bixel %g mm, longSpacing %g mm\\n', ptvVolumeCC, bixelWidth, longSpacing);

    % Body mask (patient midline reference + ring clipping)
    bodyMask = false(size(ptvMask));
    for idx = 1:size(cst, 1)
        if ~isempty(regexpi(cleanName(cst{{idx, 2}}), 'body|external|skin|aussenhaut'))
            bIdx = unique(vertcat(cst{{idx, 4}}{{:}}));
            bodyMask(bIdx(bIdx >= 1 & bIdx <= numel(bodyMask))) = true;
        end
    end

    % Gantry geometry: ipsilateral 3-field, side chosen by PTV laterality vs body midline
    [~, ptvCols, ~] = ind2sub(size(ptvMask), ptvIdx);
    ptvCx = mean(ptvCols);
    if any(bodyMask(:))
        [~, bodyCols, ~] = ind2sub(size(bodyMask), find(bodyMask));
        midlineCx = mean(bodyCols);
    else
        midlineCx = size(ptvMask, 2) / 2;
    end
    if ptvCx < midlineCx   % smaller column = smaller x = patient's right (DICOM +x = left)
        gantryAngles = {GANTRY_PTV_RIGHT};
    else                   % PTV on the patient's left
        gantryAngles = {GANTRY_PTV_LEFT};
    end
    couchAngles = zeros(1, numel(gantryAngles));
    fprintf('Gantry: %s (PTV col %.1f vs midline %.1f)\\n', mat2str(gantryAngles), ptvCx, midlineCx);

    % Dose-fall-off ring: PTV+5..PTV+12 mm, clipped to body, GTV-excluded
    ringIdx = 0;
    ptvPlusInner = matRad_addMargin(ptvMask, cst, ct.resolution, struct('x',5,'y',5,'z',5), true);
    ptvPlusOuter = matRad_addMargin(ptvMask, cst, ct.resolution, struct('x',12,'y',12,'z',12), true);
    ringMask = ptvPlusOuter & ~ptvPlusInner & ~ptvMask;

    gtvMask = false(size(ptvMask));
    for idx = 1:size(cst, 1)
        if isequal(cst{{idx, 3}}, 'TARGET') && ~isempty(regexpi(cleanName(cst{{idx, 2}}), '\\<gtv\\>'))
            gtvIdx = unique(vertcat(cst{{idx, 4}}{{:}}));
            gtvMask(gtvIdx(gtvIdx >= 1 & gtvIdx <= numel(gtvMask))) = true;
        end
    end
    ringMask = ringMask & ~gtvMask;
    if any(bodyMask(:)); ringMask = ringMask & bodyMask; end

    ringVoxels = find(ringMask);
    if ~isempty(ringVoxels)
        ringIdx = size(cst, 1) + 1;
        structIds = cell2mat(cst(cellfun(@(x) isnumeric(x) && isscalar(x), cst(:,1)), 1));
        ringStructId = max([structIds(:); 0]) + 1;
        cst{{ringIdx, 1}} = ringStructId;
        cst{{ringIdx, 2}} = 'PTV_RING_5_12MM';
        cst{{ringIdx, 3}} = 'OAR';
        cst{{ringIdx, 4}} = {{ringVoxels}};
        cst{{ringIdx, 5}} = cst{{selectedTargetIdx, 5}};
        if isempty(cst{{ringIdx, 5}}) || ~isstruct(cst{{ringIdx, 5}}); cst{{ringIdx, 5}} = struct(); end
        cst{{ringIdx, 5}}.Priority = 2;
        cst{{ringIdx, 5}}.alphaX = 0.1; cst{{ringIdx, 5}}.betaX = 0.05;
        cst{{ringIdx, 5}}.Visible = 1; cst{{ringIdx, 5}}.visibleColor = [0.1 0.4 0.9];
        cst{{ringIdx, 6}} = {{struct(DoseObjectives.matRad_SquaredOverdosing(850, 1.05*prescriptionDose))}};
    end

    % Assign objectives to all structures
    for i = 1:size(cst, 1)
        name = cleanName(cst{{i, 2}});
        if isequal(cst{{i, 3}}, 'TARGET')
            if i == selectedTargetIdx
                cst{{i, 6}} = {{
                    struct(DoseObjectives.matRad_SquaredDeviation(2600, prescriptionDose));
                    struct(DoseObjectives.matRad_MinDVH(3000, 0.95*prescriptionDose, 95));
                    struct(DoseObjectives.matRad_MaxDVH(100, 1.12*prescriptionDose, 2));
                    struct(DoseObjectives.matRad_SquaredUnderdosing(1000, 0.95*prescriptionDose))
                }};
            else
                cst{{i, 6}} = {{}};
            end

        elseif isequal(cst{{i, 3}}, 'OAR')
            if ringIdx > 0 && i == ringIdx; continue; end
            cst{{i, 6}} = {{}};
            if ~isempty(regexpi(name, 'spinal|spinalkanal'))
                cst{{i, 6}}{{1}} = struct(DoseObjectives.matRad_SquaredOverdosing(400, 50.5));
            elseif strcmp(patientType, 'lung')
                if ismember(i, lungStructureIndices)
                    cst{{i, 6}} = {{struct(DoseObjectives.matRad_MaxDVH(80, 60, 33));
                                   struct(DoseObjectives.matRad_MaxDVH(20, 5, 33));
                                   struct(DoseObjectives.matRad_MeanDose(1, 0))}};
                elseif ~isempty(regexpi(name, 'brachial|plexus'))
                    cst{{i, 6}}{{1}} = struct(DoseObjectives.matRad_SquaredOverdosing(200, 66));
                elseif ~isempty(regexpi(name, 'oesophagus|esophagus'))
                    cst{{i, 6}} = {{struct(DoseObjectives.matRad_SquaredOverdosing(200, 73.5));
                                   struct(DoseObjectives.matRad_MaxDVH(100, 60, 15.3));
                                   struct(DoseObjectives.matRad_MeanDose(2, 0))}};
                elseif ~isempty(regexpi(name, 'skin|body|external|boundary|aussenhaut'))
                    cst{{i, 6}} = {{struct(DoseObjectives.matRad_SquaredOverdosing(120, 1.10*prescriptionDose));
                                   struct(DoseObjectives.matRad_MeanDose(1, 0))}};
                end
            else  % abdomen
                if ~isempty(regexpi(name, 'liver|leber'))
                    cst{{i, 6}}{{1}} = struct(DoseObjectives.matRad_SquaredOverdosing(25, 30));
                elseif ~isempty(regexpi(name, 'kidney|niere'))
                    cst{{i, 6}}{{1}} = struct(DoseObjectives.matRad_SquaredOverdosing(30, 18));
                elseif ~isempty(regexpi(name, 'stomach|magen'))
                    cst{{i, 6}}{{1}} = struct(DoseObjectives.matRad_SquaredOverdosing(220, 45));
                elseif ~isempty(regexpi(name, 'bladder|blase'))
                    cst{{i, 6}}{{1}} = struct(DoseObjectives.matRad_SquaredOverdosing(160, 65));
                elseif ~isempty(regexpi(name, 'skin|body|external|boundary|aussenhaut'))
                    cst{{i, 6}} = {{struct(DoseObjectives.matRad_SquaredOverdosing(80, 1.10*prescriptionDose));
                                   struct(DoseObjectives.matRad_MeanDose(1, 0))}};
                end
            end
        end
    end

    % ---- configure plan (pln) ----
    pln.numOfFractions = numFractions;
    pln.radiationMode  = '{pln["radiationMode"]}';
    pln.machine        = '{pln["machine"]}';
    pln.bioModel       = '{pln["bioModel"]}';
    pln.multScen       = '{pln["multScen"]}';
    pln.propStf.gantryAngles            = gantryAngles;
    pln.propStf.couchAngles             = couchAngles;
    pln.propStf.bixelWidth              = bixelWidth;
    pln.propStf.longitudinalSpotSpacing = longSpacing;
    pln.propStf.numOfBeams              = numel(pln.propStf.gantryAngles);
    pln.propStf.isoCenter               = ones(pln.propStf.numOfBeams, 1) * matRad_getIsoCenter(cst(selectedTargetIdx,:), ct, 0);
    pln.propDoseCalc.engine             = '{pln["doseEngine"]}';
    pln.propDoseCalc.calcLET            = {pln["calcLET"]};
    pln.propDoseCalc.doseGrid.resolution.x = {pln["doseGridResolution"]};
    pln.propDoseCalc.doseGrid.resolution.y = {pln["doseGridResolution"]};
    pln.propDoseCalc.doseGrid.resolution.z = {pln["doseGridResolution"]};
    pln.propOpt.quantityOpt   = '{pln["quantityOpt"]}';
    pln.propOpt.optimizer     = '{pln["optimizer"]}';
    pln.propOpt.runDAO        = 0;
    pln.propOpt.maxIter       = {pln["maxIter"]};
    pln.propSeq.runSequencing = 0;

    % ---- stf / dij (with cache) ----
    forceRecalcDij = {"true" if force_recalc_dij else "false"};
    dijCachePath   = fullfile('{mat_file.parent}', sprintf('%s_proton_stf_dij_cache.mat', patientID));

    expectedCacheMeta = struct('patientID', patientID, 'machine', pln.machine, ...
        'radiationMode', pln.radiationMode, 'bioModel', pln.bioModel, ...
        'engine', pln.propDoseCalc.engine, 'bixelWidth', pln.propStf.bixelWidth, ...
        'longitudinalSpotSpacing', pln.propStf.longitudinalSpotSpacing, ...
        'gantryAngles', pln.propStf.gantryAngles, 'couchAngles', pln.propStf.couchAngles, ...
        'isoCenter', pln.propStf.isoCenter, ...
        'doseGrid', [pln.propDoseCalc.doseGrid.resolution.x, ...
                     pln.propDoseCalc.doseGrid.resolution.y, ...
                     pln.propDoseCalc.doseGrid.resolution.z]);

    cacheLoaded = false;
    if ~forceRecalcDij && exist(dijCachePath, 'file')
        try
            cacheData = load(dijCachePath, 'stf', 'dij', 'cacheMeta');
            if all(isfield(cacheData, {{'stf','dij','cacheMeta'}})) && isequaln(cacheData.cacheMeta, expectedCacheMeta)
                stf = cacheData.stf; dij = cacheData.dij; cacheLoaded = true;
                fprintf('Loaded stf/dij cache.\\n');
            else
                fprintf('Cache mismatch; recalculating.\\n');
            end
        catch e
            fprintf('Cache load failed (%s); recalculating.\\n', e.message);
        end
    end

    if ~cacheLoaded
        stf = matRad_generateStf(ct, cst, pln);
        dij = matRad_calcDoseInfluence(ct, cst, stf, pln);
        cacheMeta = expectedCacheMeta;
        try
            save(dijCachePath, 'stf', 'dij', 'cacheMeta', '-v7.3');
            fprintf('Saved stf/dij cache.\\n');
        catch e
            fprintf('WARNING: cache save failed (%s)\\n', e.message);
        end
    end

    % ---- optimization ----
    resultGUI = matRad_fluenceOptimization(dij, cst, pln);
    optHistory = [];  % objective convergence history (best-effort; private in some matRad versions)
    try
        if isfield(resultGUI, 'usedOptimizer') && ~isempty(resultGUI.usedOptimizer.allObjectiveFunctionValues)
            optHistory = resultGUI.usedOptimizer.allObjectiveFunctionValues;
        end
    catch
        optHistory = [];
    end

    % Filter spots with w < minAbsWeight
    minAbsWeight = {pln["minAbsWeight"]};
    w = resultGUI.w;
    spotsBeforeFilter = nnz(w > 0);
    spotsDropped = 0; droppedWeightFraction = 0;
    if minAbsWeight > 0
        droppedMask = (w > 0) & (w < minAbsWeight);
        spotsDropped = nnz(droppedMask);
        if spotsDropped > 0
            droppedWeightFraction = sum(w(droppedMask)) / max(sum(w), eps);
            w(droppedMask) = 0;
            resultGUI = matRad_calcDoseForward(ct, cst, stf, pln, w);
            resultGUI.w = w;
        end
    end
    resultGUI.numSpotsTotal         = numel(resultGUI.w);
    resultGUI.numSpotsActive        = nnz(resultGUI.w > 0);
    resultGUI.numSpotsBeforeFilter  = spotsBeforeFilter;
    resultGUI.numSpotsDropped       = spotsDropped;
    resultGUI.minAbsWeight          = minAbsWeight;
    resultGUI.droppedWeightFraction = droppedWeightFraction;
    fprintf('Spots: %d active / %d total (dropped %d)\\n', resultGUI.numSpotsActive, resultGUI.numSpotsTotal, spotsDropped);

    % ---- export dose cube (MHA) ----
    if isfield(resultGUI, 'RBExDose') && ~isempty(resultGUI.RBExDose)
        doseCube = resultGUI.RBExDose;
    else
        doseCube = resultGUI.physicalDose;
    end
    doseSize = size(doseCube);
    ctSize   = [length(ct.y), length(ct.x), length(ct.z)];
    if isequal(doseSize, ctSize)
        gridX = ct.x; gridY = ct.y; gridZ = ct.z;
        doseRes = [ct.resolution.x, ct.resolution.y, ct.resolution.z];
    else
        gridX = dij.doseGrid.x; gridY = dij.doseGrid.y; gridZ = dij.doseGrid.z;
        doseRes = [dij.doseGrid.resolution.x, dij.doseGrid.resolution.y, dij.doseGrid.resolution.z];
    end
    doseMetadata = struct('resolution', doseRes, 'imageOrigin', [gridX(1), gridY(1), gridZ(1)], 'datatype', 'single');
    matRad_writeCube(doseMhaPath, doseCube, 'single', doseMetadata);
    fprintf('Saved dose MHA.\\n');

    % ---- export spots (JSON) ----
    exportData = struct('patient_id', patientID, 'patient_type', patientType, ...
        'num_fractions', pln.numOfFractions, 'prescription_dose', prescriptionDose, ...
        'bio_model', pln.bioModel, 'quantity_opt', pln.propOpt.quantityOpt, ...
        'dose_engine', pln.propDoseCalc.engine, 'bixel_width', pln.propStf.bixelWidth, ...
        'longitudinal_spot_spacing', pln.propStf.longitudinalSpotSpacing, ...
        'selected_ptv', selectedPTVName, ...
        'num_spots_total', resultGUI.numSpotsTotal, 'num_spots_active', resultGUI.numSpotsActive, ...
        'num_spots_before_filter', resultGUI.numSpotsBeforeFilter, ...
        'num_spots_dropped_filter', resultGUI.numSpotsDropped, ...
        'min_abs_weight', resultGUI.minAbsWeight, ...
        'dropped_weight_fraction', resultGUI.droppedWeightFraction);

    beamCells = cell(1, numel(stf));
    beamletCounter = 0; keptCounter = 0; perBeamKept = zeros(1, numel(stf));
    for b = 1:numel(stf)
        beam = struct('beam_idx', b-1, 'gantry_angle', stf(b).gantryAngle, ...
            'couch_angle', stf(b).couchAngle, 'SAD', stf(b).SAD, 'iso_center', stf(b).isoCenter);
        rayCells = {{}};
        for r = 1:numel(stf(b).ray)
            nBixels = numel(stf(b).ray(r).energy);
            beamletCells = {{}};
            for k = 1:nBixels
                beamletCounter = beamletCounter + 1;
                wk = resultGUI.w(beamletCounter);
                if wk <= 0; continue; end
                bl = struct('beamlet_idx', k-1, 'energy', stf(b).ray(r).energy(k), 'weight', wk);
                if isfield(stf(b).ray(r), 'focusIx') && numel(stf(b).ray(r).focusIx) >= k
                    bl.focusIx = stf(b).ray(r).focusIx(k);
                end
                beamletCells{{end+1}} = bl;
                keptCounter = keptCounter + 1; perBeamKept(b) = perBeamKept(b) + 1;
            end
            if ~isempty(beamletCells)
                rayCells{{end+1}} = struct('ray_idx', r-1, 'rayPos_bev', stf(b).ray(r).rayPos_bev, 'beamlets', {{beamletCells}});
            end
        end
        beam.rays = rayCells; beam.num_kept_spots = perBeamKept(b);
        beamCells{{b}} = beam;
    end
    exportData.beams = beamCells;
    exportData.total_beamlet_count = keptCounter;
    exportData.per_beam_kept_spots = perBeamKept;
    fid = fopen(spotJson, 'w');
    if fid ~= -1; fprintf(fid, '%s', jsonencode(exportData, 'PrettyPrint', true)); fclose(fid); end
    fprintf('Saved JSON (%d spots).\\n', keptCounter);

    % ---- DVH structures ----
    targetPlotNames = {{selectedPTVName}};
    clinicalPTVName = '';
    for j = 1:size(cst, 1)
        if ~isequal(cst{{j, 3}}, 'TARGET'); continue; end
        nm = cst{{j, 2}};
        if ~strcmp(nm, selectedPTVName) && ~isempty(regexpi(cleanName(nm), '^ptv$'))
            clinicalPTVName = nm; targetPlotNames{{end+1}} = nm; break;
        end
    end
    selectedGTVName = ''; bestGtvScore = -inf;
    for j = 1:size(cst, 1)
        if ~isequal(cst{{j, 3}}, 'TARGET'); continue; end
        nm = cst{{j, 2}}; nmClean = cleanName(nm);
        if ~contains(nmClean, 'gtv'); continue; end
        gv = numel(unique(vertcat(cst{{j, 4}}{{:}}))); sc = gv;
        if ~isempty(regexpi(nmClean, '^gtv$')); sc = 1e12+gv;
        elseif ~isempty(regexpi(nmClean, '^gtv[ _-]*union$|^gtvunion$')); sc = 1e11+gv; end
        if sc > bestGtvScore; bestGtvScore = sc; selectedGTVName = nm; end
    end
    if ~isempty(selectedGTVName) && ~any(strcmp(targetPlotNames, selectedGTVName))
        targetPlotNames{{end+1}} = selectedGTVName;
    end

    if strcmp(patientType, 'lung')
        oarPriorityPatterns = {{'lung','spinal|spinalkanal','oesophagus|esophagus','brachial|plexus','heart','skin|aussenhaut','body|external'}};
    else
        oarPriorityPatterns = {{'stomach|magen','liver|leber','kidney|niere','spinal|spinalkanal','duodenum|bowel|darm','bladder|blase','skin|aussenhaut','body|external'}};
    end
    selectedOars = {{}};
    for p = 1:length(oarPriorityPatterns)
        if length(selectedOars) >= 5; break; end
        for j = 1:size(cst, 1)
            if ~isequal(cst{{j, 3}}, 'OAR'); continue; end
            nm = cst{{j, 2}}; nmLower = lower(nm);
            if contains(nmLower,'ptv')||contains(nmLower,'gtv')||contains(nmLower,'ctv')||contains(nmLower,'target'); continue; end
            if ~isempty(regexpi(nmLower, oarPriorityPatterns{{p}})) && ~any(strcmp(selectedOars, nm))
                selectedOars{{end+1}} = nm;
                if length(selectedOars) >= 5; break; end
            end
        end
    end

    % ---- DVH report (TXT) ----
    refVol = [2 5 50 95 98];
    refGyTotal = [45 50 60]; if strcmp(patientType, 'lung'); refGyTotal = [5 60 70]; end
    qi = matRad_calcQualityIndicators(cst, pln, doseCube, refGyTotal/numFractions, refVol);
    qiByName = containers.Map(); for i = 1:length(qi); qiByName(qi(i).name) = i; end
    reportList = [targetPlotNames, selectedOars];

    fid = fopen(reportPath, 'w');
    if fid ~= -1
        fprintf(fid, '=== Proton DVH Report ===\\n');
        fprintf(fid, 'Patient ID:       %s\\n', patientID);
        fprintf(fid, 'Patient Type:     %s\\n', patientType);
        fprintf(fid, 'Prescription:     %d Gy(RBE) / %d fractions\\n', prescriptionDose, numFractions);
        fprintf(fid, 'Selected PTV:     %s  *** OPTIMIZATION TARGET ***\\n', selectedPTVName);
        if ~isempty(clinicalPTVName)
            fprintf(fid, 'Clinical PTV:     %s  *** PLAN ACCEPTANCE TARGET ***\\n', clinicalPTVName);
        end
        fprintf(fid, 'Bio Model:        %s (quantityOpt=%s)\\n', pln.bioModel, pln.propOpt.quantityOpt);
        fprintf(fid, 'Dose Engine:      %s\\n', pln.propDoseCalc.engine);
        fprintf(fid, 'Bixel Width:      %.1f mm,  LongSpacing: %.1f mm\\n', pln.propStf.bixelWidth, pln.propStf.longitudinalSpotSpacing);
        fprintf(fid, 'Dose Grid:        %.1f x %.1f x %.1f mm\\n', doseRes(1), doseRes(2), doseRes(3));
        fprintf(fid, '\\n--- Beam Geometry ---\\n');
        for b = 1:numel(stf)
            fprintf(fid, '  Beam %d: gantry=%.1f deg, couch=%.1f deg, kept_spots=%d\\n', b, stf(b).gantryAngle, stf(b).couchAngle, perBeamKept(b));
        end
        fprintf(fid, '\\n--- Pencil Beam Stats ---\\n');
        fprintf(fid, 'Total spots generated:    %d\\n', resultGUI.numSpotsTotal);
        fprintf(fid, 'Active before filter:     %d  (post-IPOPT, w > 0)\\n', resultGUI.numSpotsBeforeFilter);
        fprintf(fid, 'Dropped by filter:        %d  (w < %.4g)\\n', resultGUI.numSpotsDropped, resultGUI.minAbsWeight);
        fprintf(fid, 'Final active spots:       %d\\n', resultGUI.numSpotsActive);
        fprintf(fid, 'Dropped weight fraction:  %.6f%% of total\\n', 100*resultGUI.droppedWeightFraction);
        fprintf(fid, '\\n--- Key Structures (total Gy(RBE)) ---\\n');
        for k = 1:length(reportList)
            structName = reportList{{k}};
            if ~qiByName.isKey(structName); continue; end
            q = qi(qiByName(structName));
            roleTag = '';
            if strcmp(structName, selectedPTVName); roleTag = '  >>> SELECTED PTV (optimized) <<<';
            elseif strcmp(structName, clinicalPTVName); roleTag = '  >>> CLINICAL PTV (acceptance) <<<';
            elseif any(strcmp(targetPlotNames, structName)); roleTag = '  (GTV — informational)'; end
            fprintf(fid, '\\n[%s]%s\\n', structName, roleTag);
            fprintf(fid, '  Dmean: %.2f Gy(RBE)\\n', q.mean*numFractions);
            fprintf(fid, '  Dmax:  %.2f Gy(RBE)\\n', q.max*numFractions);
            fprintf(fid, '  Dmin:  %.2f Gy(RBE)\\n', q.min*numFractions);
            for fld = {{'D_2','D_5','D_50','D_95','D_98'}}
                if isfield(q, fld{{1}}) && ~isempty(q.(fld{{1}}))
                    fprintf(fid, '  %-5s: %.2f Gy(RBE)\\n', strrep(fld{{1}},'_','%'), q.(fld{{1}})*numFractions);
                end
            end
            for g = 1:numel(refGyTotal)
                fnm = strrep(sprintf('V_%g', refGyTotal(g)/numFractions), '.', '_');
                if isfield(q, fnm) && ~isempty(q.(fnm))
                    fprintf(fid, '  V%dGy: %.2f%%\\n', refGyTotal(g), q.(fnm)*100);
                end
            end
            nmLow = lower(structName);
            if contains(nmLow,'spinal')
                fprintf(fid, '  Check Dmax<50.5 Gy %s\\n', passFail{{(q.max*numFractions<50.5)+1}});
            elseif strcmp(patientType,'lung')
                if contains(nmLow,'lung'); fprintf(fid, '  Check Dmean<18 Gy %s\\n', passFail{{(q.mean*numFractions<18)+1}});
                elseif contains(nmLow,'oesophagus')||contains(nmLow,'esophagus'); fprintf(fid, '  Check Dmax<73.5 Gy %s\\n', passFail{{(q.max*numFractions<73.5)+1}});
                elseif contains(nmLow,'brachial')||contains(nmLow,'plexus'); fprintf(fid, '  Check Dmax<66 Gy %s\\n', passFail{{(q.max*numFractions<66)+1}}); end
            else
                if contains(nmLow,'liver')||contains(nmLow,'leber'); fprintf(fid, '  Check Dmean<30 Gy %s\\n', passFail{{(q.mean*numFractions<30)+1}});
                elseif contains(nmLow,'kidney')||contains(nmLow,'niere'); fprintf(fid, '  Check Dmean<18 Gy %s\\n', passFail{{(q.mean*numFractions<18)+1}});
                elseif contains(nmLow,'stomach')||contains(nmLow,'magen'); fprintf(fid, '  Check Dmax<45 Gy %s\\n', passFail{{(q.max*numFractions<45)+1}});
                elseif contains(nmLow,'bladder')||contains(nmLow,'blase'); fprintf(fid, '  Check Dmax<65 Gy %s\\n', passFail{{(q.max*numFractions<65)+1}}); end
            end
        end
        fclose(fid);
    end

    % ---- DVH plot (PNG) ----
    try
        dvh = matRad_calcDVH(cst, doseCube);
        dvhNames = {{}}; for i = 1:length(dvh); dvhNames{{i}} = dvh(i).name; end
        oarColors = [0.00 0.45 0.74; 0.20 0.65 0.20; 0.00 0.70 0.70; 0.45 0.35 0.75; 0.80 0.60 0.00];
        dvhFig = figure('visible','off','Position',[100 100 1000 600]); hold on;
        legendHandles = []; legendNames = {{}}; oarIdx = 0;
        for k = 1:length([selectedOars, targetPlotNames])
            plotOrder = [selectedOars, targetPlotNames];
            structName = plotOrder{{k}};
            di = find(strcmp(dvhNames, structName), 1);
            if isempty(di) || isempty(dvh(di).volumePoints) || all(dvh(di).volumePoints==0); continue; end
            if strcmp(structName, selectedPTVName); col=[1 0 0]; lw=3.0;
            elseif strcmp(structName, clinicalPTVName); col=[0.55 0.05 0.05]; lw=2.6;
            elseif any(strcmp(targetPlotNames, structName)); col=[0.3 0 0]; lw=2.2;
            else; oarIdx=oarIdx+1; col=oarColors(mod(oarIdx-1,size(oarColors,1))+1,:); lw=1.6; end
            h = plot(dvh(di).doseGrid*numFractions, dvh(di).volumePoints, '-', 'Color', col, 'LineWidth', lw);
            legendHandles(end+1) = h; legendNames{{end+1}} = structName;
        end
        hold off;
        xlabel('Dose [Gy(RBE)]','FontSize',12); ylabel('Volume [%]','FontSize',12);
        title(sprintf('DVH - %s (%s, %dGy/%dfx, %d spots)', patientID, patientType, prescriptionDose, numFractions, resultGUI.numSpotsActive),'FontSize',13);
        grid on; xlim([0, prescriptionDose*1.15]); ylim([0,105]);
        xline(prescriptionDose,'--k',sprintf('Rx: %dGy',prescriptionDose),'LineWidth',1.5,'LabelHorizontalAlignment','left');
        if ~isempty(legendHandles); legend(legendHandles, legendNames,'Location','eastoutside','Interpreter','none','FontSize',9); end
        saveas(dvhFig, dvhPngPath); close(dvhFig);
    catch e
        fprintf('WARNING: DVH plot failed (%s)\\n', e.message);
    end

    save(resultMat, 'resultGUI', 'pln', 'stf', 'cst', 'ct', 'optHistory', 'selectedPTVName', '-v7.3');
    fprintf('Done: %s\\n', resultMat);

catch err
    fprintf('ERROR: %s\\n', err.message);
    for k = 1:length(err.stack)
        fprintf('  %s line %d (%s)\\n', err.stack(k).file, err.stack(k).line, err.stack(k).name);
    end
    exit(1);
end
"""


# ---------------------------------------------------------------------------
# Run helpers
# ---------------------------------------------------------------------------

def run_matlab(script: str, timeout: int = 7200) -> bool:
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
        return subprocess.run(cmd, text=True, timeout=timeout).returncode == 0

    except subprocess.TimeoutExpired:
        print(f"ERROR: MATLAB timed out (>{timeout}s)")
        return False
    except Exception as exc:
        print(f"ERROR: could not launch MATLAB — {exc}")
        return False
    finally:
        tmp.unlink(missing_ok=True)


def optimize_patient(patient_id: str, pln_params: dict, timeout: int, force_recalc_dij: bool) -> bool:
    """Run IMPT optimization for one patient."""
    mat_file = MAT_ROOT / patient_id / f"{patient_id}.mat"

    if not mat_file.exists():
        print(f"  [SKIP] MAT file not found: {mat_file}  (run 01_dicom2mat.py first)")
        return False

    print(f"  {mat_file}")
    return run_matlab(
        build_matlab_script(patient_id, mat_file, pln_params, force_recalc_dij),
        timeout=timeout,
    )


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> bool:
    parser = argparse.ArgumentParser(
        description="IMPT treatment planning & optimization using matRad."
    )
    parser.add_argument(
        "patient", nargs="?", default=None,
        metavar="PATIENT_ID",
        help="Patient folder name inside data/output/ (omit to process all)",
    )
    parser.add_argument(
        "-t", "--timeout", type=int, default=7200,
        metavar="SECONDS",
        help="MATLAB timeout in seconds (default: 7200)",
    )
    parser.add_argument(
        "-c", "--config", type=Path, default=None,
        metavar="FILE",
        help="JSON file with custom planning parameters",
    )
    parser.add_argument(
        "--force-recalc-dij", action="store_true",
        help="Force stf/dij recalculation even if a valid cache exists",
    )
    args = parser.parse_args()

    pln_params = json.loads(args.config.read_text()) if args.config and args.config.exists() else {}

    if args.patient:
        patients = [args.patient]
    else:
        patients = sorted(
            p.name for p in MAT_ROOT.iterdir()
            if p.is_dir() and (p / f"{p.name}.mat").exists()
        )

    print(f"Processing {len(patients)} patient(s)...\n")

    succeeded, failed = 0, []
    for pid in patients:
        print(f"[{pid}]")
        if optimize_patient(pid, pln_params, args.timeout, args.force_recalc_dij):
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
