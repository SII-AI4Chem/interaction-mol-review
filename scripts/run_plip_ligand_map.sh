#!/usr/bin/env bash

set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <input.cif|pdb> [hetid] [chain] [resid] [scale] [output_dir]" >&2
  exit 1
fi

INPUT=$1
HETID=${2:-auto}
CHAIN=${3:-auto}
RESID=${4:-auto}
SCALE=${5:-55}
OUTPUT_ROOT=${6:-$(pwd)}

mkdir -p "$OUTPUT_ROOT"
OUTPUT_ROOT=$(cd "$OUTPUT_ROOT" && pwd)

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

if [[ "$HETID" == auto || "$CHAIN" == auto || "$RESID" == auto ]]; then
  PLIP_TAG="auto"
else
  PLIP_TAG="${HETID}_${CHAIN}_${RESID}"
fi

PLIP_OUT_DIR="${OUTPUT_ROOT}/plip_output_${PLIP_TAG}"
mkdir -p "$PLIP_OUT_DIR"

STRUCT_BASE=${INPUT##*/}
STRUCT_BASE=${STRUCT_BASE%.*}
PDB_FILE="${OUTPUT_ROOT}/${STRUCT_BASE}.pdb"

TOTAL_STEPS=6

echo "[1/${TOTAL_STEPS}] Converting ${INPUT} -> ${PDB_FILE}" >&2
obabel "$INPUT" -O "$PDB_FILE" --addtotitle "$STRUCT_BASE" >/dev/null

echo "[2/${TOTAL_STEPS}] Running PLIP (output: $PLIP_OUT_DIR)" >&2
python -m plip.plipcmd -f "$PDB_FILE" -o "$PLIP_OUT_DIR" -t -x -y >/dev/null

if [[ "$HETID" == auto || "$CHAIN" == auto || "$RESID" == auto ]]; then
  REPORT_XML="${PLIP_OUT_DIR}/report.xml"
  read HETID CHAIN RESID < <(REPORT_XML="$REPORT_XML" python - <<'PY'
import os
import xml.etree.ElementTree as ET
from pathlib import Path

report_path = Path(os.environ["REPORT_XML"])
if not report_path.exists():
    raise SystemExit(f"Missing report XML: {report_path}")

root = ET.parse(report_path).getroot()
for bs in root.findall('bindingsite'):
    if bs.get('has_interactions', 'True') != 'False':
        identifiers = bs.find('identifiers')
        hetid = identifiers.findtext('hetid', default='') if identifiers is not None else ''
        chain = identifiers.findtext('chain', default='') if identifiers is not None else ''
        position = identifiers.findtext('position', default='') if identifiers is not None else ''
        print(hetid, chain, position)
        break
else:
    bs = root.find('bindingsite')
    identifiers = bs.find('identifiers') if bs is not None else None
    hetid = identifiers.findtext('hetid', default='UNK') if identifiers is not None else 'UNK'
    chain = identifiers.findtext('chain', default='') if identifiers is not None else ''
    position = identifiers.findtext('position', default='') if identifiers is not None else ''
    print(hetid, chain, position)
PY
)
  HETID=$(printf '%s' "$HETID" | tr '[:lower:]' '[:upper:]')
  CHAIN=$(printf '%s' "$CHAIN" | tr '[:lower:]' '[:upper:]')
  RESID=${RESID}
  if [[ -z "$HETID" || -z "$CHAIN" || -z "$RESID" ]]; then
    echo "Failed to auto-detect ligand parameters." >&2
    exit 2
  fi
  echo "[3/${TOTAL_STEPS}] Auto-detected ligand ${HETID} chain ${CHAIN} resid ${RESID}" >&2
else
  HETID=$(printf '%s' "$HETID" | tr '[:lower:]' '[:upper:]')
  CHAIN=$(printf '%s' "$CHAIN" | tr '[:lower:]' '[:upper:]')
  echo "[3/${TOTAL_STEPS}] Using ligand ${HETID} chain ${CHAIN} resid ${RESID}" >&2
fi

NEW_TAG="${HETID}_${CHAIN}_${RESID}"
NEW_PLIP_DIR="${OUTPUT_ROOT}/plip_output_${NEW_TAG}"
if [[ "$PLIP_OUT_DIR" != "$NEW_PLIP_DIR" ]]; then
  rm -rf "$NEW_PLIP_DIR"
  mv "$PLIP_OUT_DIR" "$NEW_PLIP_DIR"
  PLIP_OUT_DIR="$NEW_PLIP_DIR"
fi

LIG_PDB="${OUTPUT_ROOT}/ligand_${HETID}_${CHAIN}_${RESID}.pdb"
LIG_SDF="${OUTPUT_ROOT}/ligand_${HETID}_${CHAIN}_${RESID}_2d.sdf"

echo "[4/${TOTAL_STEPS}] Extracting ligand ${HETID} chain ${CHAIN} resid ${RESID} -> ${LIG_PDB}" >&2
python "${SCRIPT_DIR}/extract_ligand.py" "$PDB_FILE" "$LIG_PDB" "$HETID" "$CHAIN" "$RESID"

echo "[5/${TOTAL_STEPS}] Generating 2D SDF ${LIG_SDF}" >&2
obabel "$LIG_PDB" -O "$LIG_SDF" --gen2d >/dev/null

MAP_OUTPUT="${PLIP_OUT_DIR}/mol_${HETID}_${CHAIN}_${RESID}.svg"
echo "[6/${TOTAL_STEPS}] Rendering interaction map -> ${MAP_OUTPUT}" >&2
python "${SCRIPT_DIR}/mol_map6.py" \
  --report "${PLIP_OUT_DIR}/report.xml" \
  --ligand-pdb "$LIG_PDB" \
  --ligand-sdf "$LIG_SDF" \
  --style basic \
  --scale "$SCALE" \
  --output "$MAP_OUTPUT"

echo "Done. Files written to ${PLIP_OUT_DIR} and ${OUTPUT_ROOT}." >&2
