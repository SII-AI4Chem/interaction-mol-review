#!/usr/bin/env python3
"""Extract a single ligand residue from a PDB file."""

from __future__ import annotations

import argparse
import pathlib


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_pdb", type=pathlib.Path, help="Source PDB file")
    parser.add_argument("output_pdb", type=pathlib.Path, help="Destination PDB file")
    parser.add_argument("hetid", help="Ligand residue name / hetid (e.g. A1I)")
    parser.add_argument("chain", help="Chain identifier")
    parser.add_argument("resid", help="Residue number (as it appears in the PDB, e.g. 0 or 169)")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    hetid = args.hetid.upper()
    chain = args.chain.upper()
    resid = args.resid.strip()

    matching_lines: list[str] = []

    with args.input_pdb.open() as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if line[17:20].strip().upper() != hetid:
                continue
            if line[21].strip().upper() != chain:
                continue
            if line[22:26].strip() != resid:
                continue
            matching_lines.append(line.rstrip("\n"))

    if not matching_lines:
        raise SystemExit(f"No atoms found for {hetid} chain {chain} resid {resid} in {args.input_pdb}")

    with args.output_pdb.open("w") as out:
        for line in matching_lines:
            out.write(line + "\n")
        out.write("TER\nEND\n")

    print(f"Extracted {len(matching_lines)} atoms to {args.output_pdb}")


if __name__ == "__main__":
    main()
