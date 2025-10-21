#!/usr/bin/env python3
"""
Variant of mol_map2.py that renders atoms using coloured element letters
instead of coloured filled circles.
"""

from __future__ import annotations

import argparse
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
import xml.etree.ElementTree as ET


ELEMENT_COLORS: Dict[str, str] = {
    "H": "#95a5a6",
    "C": "#2c3e50",
    "N": "#1a4d7a",
    "O": "#c0392b",
    "S": "#f39c12",
    "P": "#d35400",
    "F": "#27ae60",
    "Cl": "#16a085",
    "Br": "#8e44ad",
    "I": "#5e35b1",
    "B": "#0d47a1",
    "Si": "#6d4c41",
}


@dataclass
class Atom:
    serial: int
    name: str
    element: str
    coords3d: Tuple[float, float, float]
    x2d: float = 0.0
    y2d: float = 0.0
    depth: float = 0.0


@dataclass
class Interaction:
    kind: str
    residue_label: str
    residue_display: str
    ligand_serial: int
    ligand_point: Tuple[float, float]
    metadata: Dict[str, float]
    count: int = 1


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--ligand-pdb", required=True, type=Path)
    parser.add_argument("--ligand-sdf", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--scale", type=float, default=65.0)
    parser.add_argument("--margin", type=float, default=90.0)
    parser.add_argument("--style", choices=["basic", "shaded"], default="basic")
    return parser.parse_args()


def load_ligand_atoms(pdb_path: Path, sdf_path: Path) -> Tuple[List[Atom], List[Tuple[int, int, int]]]:
    pdb_atoms: List[Dict[str, object]] = []
    with pdb_path.open() as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")):
                pdb_atoms.append(
                    {
                        "serial": int(line[6:11]),
                        "name": line[12:16].strip(),
                        "element": (line[76:78].strip() or line[12:16].strip()[0]).capitalize(),
                        "coords3d": (
                            float(line[30:38]),
                            float(line[38:46]),
                            float(line[46:54]),
                        ),
                    }
                )

    sdf_lines = sdf_path.read_text().splitlines()
    counts_line = sdf_lines[3]
    num_atoms = int(counts_line[:3])
    num_bonds = int(counts_line[3:6])
    if num_atoms != len(pdb_atoms):
        raise ValueError("Atom count mismatch between PDB and SDF inputs.")

    atoms: List[Atom] = []
    for idx in range(num_atoms):
        line = sdf_lines[4 + idx]
        pdb_atom = pdb_atoms[idx]
        atoms.append(
            Atom(
                serial=pdb_atom["serial"],  # type: ignore[arg-type]
                name=pdb_atom["name"],  # type: ignore[arg-type]
                element=pdb_atom["element"],  # type: ignore[arg-type]
                coords3d=pdb_atom["coords3d"],  # type: ignore[arg-type]
                x2d=float(line[0:10]),
                y2d=float(line[10:20]),
            )
        )

    bonds: List[Tuple[int, int, int]] = []
    for idx in range(num_bonds):
        line = sdf_lines[4 + num_atoms + idx]
        a1 = int(line[0:3]) - 1
        a2 = int(line[3:6]) - 1
        order = int(line[6:9])
        bonds.append((a1, a2, order))

    return atoms, bonds


def _make_label(resname: str, resnr: str, chain: str) -> Tuple[str, str]:
    compact = f"{resname}{resnr}{chain}"
    display = f"{resname} {resnr} ({chain})"
    return compact, display


def extract_interactions(root: ET.Element, atoms: Sequence[Atom]) -> List[Interaction]:
    serial_lookup = {atom.serial: atom for atom in atoms}
    interactions: List[Interaction] = []

    for hb in root.findall(".//hydrogen_bond"):
        resname = hb.findtext("restype", default="UNK")
        resnr = hb.findtext("resnr", default="?")
        chain = hb.findtext("reschain", default="")
        prot_is_donor = hb.findtext("protisdon", default="False") == "True"
        donor_idx = int(hb.findtext("donoridx"))
        acceptor_idx = int(hb.findtext("acceptoridx"))
        ligand_serial = acceptor_idx if prot_is_donor else donor_idx
        atom = serial_lookup.get(ligand_serial)
        if atom is None:
            continue
        dist = float(hb.findtext("dist_d-a", default="0"))
        label, display = _make_label(resname, resnr, chain)
        interactions.append(
            Interaction(
                kind="hydrogen",
                residue_label=label,
                residue_display=display,
                ligand_serial=ligand_serial,
                ligand_point=(atom.x2d, atom.y2d),
                metadata={"distances": [dist]},
            )
        )

    for hyd in root.findall(".//hydrophobic_interaction"):
        resname = hyd.findtext("restype", default="UNK")
        resnr = hyd.findtext("resnr", default="?")
        chain = hyd.findtext("reschain", default="")
        ligand_serial = int(hyd.findtext("ligcarbonidx"))
        atom = serial_lookup.get(ligand_serial)
        if atom is None:
            continue
        dist = float(hyd.findtext("dist", default="0"))
        label, display = _make_label(resname, resnr, chain)
        interactions.append(
            Interaction(
                kind="hydrophobic",
                residue_label=label,
                residue_display=display,
                ligand_serial=ligand_serial,
                ligand_point=(atom.x2d, atom.y2d),
                metadata={"distances": [dist]},
            )
        )

    return interactions


def aggregate_interactions(interactions: Sequence[Interaction]) -> List[Interaction]:
    grouped: Dict[Tuple[str, str, int], Interaction] = {}
    order: List[Tuple[str, str, int]] = []

    for interaction in interactions:
        key = (interaction.kind, interaction.residue_label, interaction.ligand_serial)
        if key not in grouped:
            grouped[key] = Interaction(
                kind=interaction.kind,
                residue_label=interaction.residue_label,
                residue_display=interaction.residue_display,
                ligand_serial=interaction.ligand_serial,
                ligand_point=interaction.ligand_point,
                metadata={"distances": list(interaction.metadata["distances"])},
                count=interaction.count,
            )
            order.append(key)
        else:
            grouped[key].metadata["distances"].extend(interaction.metadata["distances"])
            grouped[key].count += interaction.count

    aggregated = []
    for key in order:
        aggregated.append(grouped[key])
    return aggregated


def format_distance_block(interaction: Interaction) -> str:
    dists = interaction.metadata.get("distances", [])
    if not dists:
        return ""
    formatted = "/".join(f"{dist:.2f}" for dist in dists)
    label = "H-bond" if interaction.kind == "hydrogen" else "Hydrophobic"
    if len(dists) > 1:
        return f"{label} {formatted} Å (x{len(dists)})"
    return f"{label} {formatted} Å"


def assign_projected_coordinates(atoms: Sequence[Atom]) -> None:
    coords = np.array([atom.coords3d for atom in atoms])
    centroid = coords.mean(axis=0)
    centered = coords - centroid
    if np.allclose(centered, 0):
        return
    _, _, vh = np.linalg.svd(centered, full_matrices=False)
    basis = vh[:3]
    projected = centered @ basis.T
    for atom, proj in zip(atoms, projected):
        atom.x2d = float(proj[0])
        atom.y2d = float(proj[1])
        atom.depth = float(proj[2])


def _clamp(value: float, lower: float = 0.0, upper: float = 255.0) -> int:
    return int(max(lower, min(upper, value)))


def _adjust_color(hex_color: str, factor: float) -> str:
    hex_color = hex_color.lstrip("#")
    r = int(hex_color[0:2], 16)
    g = int(hex_color[2:4], 16)
    b = int(hex_color[4:6], 16)
    if factor >= 0:
        r = _clamp(r + (255 - r) * factor)
        g = _clamp(g + (255 - g) * factor)
        b = _clamp(b + (255 - b) * factor)
    else:
        r = _clamp(r * (1 + factor))
        g = _clamp(g * (1 + factor))
        b = _clamp(b * (1 + factor))
    return f"#{r:02x}{g:02x}{b:02x}"


def _svg_coords(
    x: float,
    y: float,
    min_x: float,
    min_y: float,
    scale: float,
    margin: float,
    height: float,
) -> Tuple[float, float]:
    sx = (x - min_x) * scale + margin
    sy = height - ((y - min_y) * scale + margin)
    return sx, sy


def _bond_segments(
    start: Tuple[float, float],
    end: Tuple[float, float],
    order: int,
    gap: float,
) -> Iterable[Tuple[Tuple[float, float], Tuple[float, float]]]:
    x1, y1 = start
    x2, y2 = end
    dx = x2 - x1
    dy = y2 - y1
    length = math.hypot(dx, dy)
    if length == 0:
        return []
    ox = -dy / length
    oy = dx / length
    if order == 1:
        yield (start, end)
    elif order == 2:
        yield ((x1 + ox * gap, y1 + oy * gap), (x2 + ox * gap, y2 + oy * gap))
        yield ((x1 - ox * gap, y1 - oy * gap), (x2 - ox * gap, y2 - oy * gap))
    elif order == 3:
        yield (start, end)
        yield ((x1 + ox * gap, y1 + oy * gap), (x2 + ox * gap, y2 + oy * gap))
        yield ((x1 - ox * gap, y1 - oy * gap), (x2 - ox * gap, y2 - oy * gap))
    else:
        yield (start, end)


def _interaction_style(kind: str) -> Tuple[str, str, float]:
    if kind == "hydrogen":
        return "#1f77b4", "6 6", 2.8
    if kind == "hydrophobic":
        return "#e65100", "3 10", 2.8
    return "#424242", "6 6", 2.5


def _rotate(vx: float, vy: float, angle: float) -> Tuple[float, float]:
    cos_a = math.cos(angle)
    sin_a = math.sin(angle)
    return vx * cos_a - vy * sin_a, vx * sin_a + vy * cos_a


def _find_label_position(
    atom_x: float,
    atom_y: float,
    base_angle: float,
    base_radius: float,
    rect_w: float,
    rect_h: float,
    bounds: Tuple[float, float, float, float],
    placed: List[Tuple[float, float, float, float]],
    radial_step: float = 32.0,
    angle_step_deg: float = 15.0,
    max_radius_steps: int = 12,
) -> Tuple[float, float, Tuple[float, float, float, float]]:
    min_x_bound, min_y_bound, max_x_bound, max_y_bound = bounds
    angle_offsets = [0]
    for k in range(1, 8):
        angle_offsets.extend([k * angle_step_deg, -k * angle_step_deg])

    for r_idx in range(max_radius_steps):
        radius = base_radius + radial_step * r_idx
        for offset in angle_offsets:
            angle = base_angle + math.radians(offset)
            ux = math.cos(angle)
            uy = math.sin(angle)
            svg_ux, svg_uy = ux, -uy

            center_x = atom_x + svg_ux * radius
            center_y = atom_y + svg_uy * radius

            center_x = max(min_x_bound + rect_w / 2, min(max_x_bound - rect_w / 2, center_x))
            center_y = max(min_y_bound + rect_h / 2, min(max_y_bound - rect_h / 2, center_y))

            box = (center_x - rect_w / 2, center_y - rect_h / 2, center_x + rect_w / 2, center_y + rect_h / 2)

            if not any(
                box[0] < existing[2] and box[2] > existing[0] and box[1] < existing[3] and box[3] > existing[1]
                for existing in placed
            ):
                placed.append(box)
                return center_x, center_y, box

    placed.append(box)
    return center_x, center_y, box


def draw_svg(
    atoms: Sequence[Atom],
    bonds: Sequence[Tuple[int, int, int]],
    interactions: Sequence[Interaction],
    output_path: Path,
    scale: float,
    margin: float,
    style: str,
) -> None:
    xs = [atom.x2d for atom in atoms]
    ys = [atom.y2d for atom in atoms]
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    width = (max_x - min_x) * scale + margin * 2.0
    height = (max_y - min_y) * scale + margin * 2.0

    centroid_x = sum(xs) / len(xs)
    centroid_y = sum(ys) / len(ys)

    defs: List[str] = [
        '<linearGradient id="bg-grad" x1="0%" y1="0%" x2="0%" y2="100%">'
        '<stop offset="0%" stop-color="#ffffff"/><stop offset="100%" stop-color="#f0f3f8"/></linearGradient>',
        '<marker id="arrow-blue" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto" markerUnits="strokeWidth">'
        '<path d="M0,0 L6,3 L0,6 z" fill="#1f77b4"/></marker>',
        '<marker id="arrow-orange" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto" markerUnits="strokeWidth">'
        '<path d="M0,0 L6,3 L0,6 z" fill="#e65100"/></marker>',
    ]

    svg_lines: List[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width:.0f}" height="{height:.0f}" viewBox="0 0 {width:.1f} {height:.1f}">',
        "<defs>",
        *defs,
        "</defs>",
        '<rect x="0" y="0" width="100%" height="100%" fill="url(#bg-grad)"/>',
    ]

    for idx1, idx2, order in bonds:
        start = _svg_coords(atoms[idx1].x2d, atoms[idx1].y2d, min_x, min_y, scale, margin, height)
        end = _svg_coords(atoms[idx2].x2d, atoms[idx2].y2d, min_x, min_y, scale, margin, height)
        for (x1, y1), (x2, y2) in _bond_segments(start, end, order, gap=5.0):
            svg_lines.append(
                f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="#404040" '
                'stroke-width="2.1" stroke-linecap="round"/>'
            )

    highlighted_serials = {interaction.ligand_serial for interaction in interactions}

    for atom in atoms:
        x_svg, y_svg = _svg_coords(atom.x2d, atom.y2d, min_x, min_y, scale, margin, height)
        element = atom.element
        colour = ELEMENT_COLORS.get(element, '#34495e')
        highlighted = atom.serial in highlighted_serials
        outer_radius = 12 if highlighted else 9
        inner_radius = outer_radius - (3.0 if highlighted else 2.4)
        stroke_colour = '#111111' if highlighted else '#2c3e50'
        svg_lines.append(
            f'<circle cx="{x_svg:.1f}" cy="{y_svg:.1f}" r="{outer_radius:.1f}" fill="#fefefe" '            f'stroke="{stroke_colour}" stroke-width="2.6"/>'
        )
        svg_lines.append(
            f'<circle cx="{x_svg:.1f}" cy="{y_svg:.1f}" r="{inner_radius:.1f}" fill="#ffffff" '            'stroke="#d9d9d9" stroke-width="0.6" opacity="0.9"/>'
        )
        font_size = 13 if highlighted else 12
        svg_lines.append(
            f'<text x="{x_svg:.1f}" y="{y_svg + 4:.1f}" font-family="Helvetica" '            f'font-size="{font_size}" text-anchor="middle" fill="{colour}" font-weight="bold">{element}</text>'
        )
    marker_map = {"hydrogen": "url(#arrow-blue)", "hydrophobic": "url(#arrow-orange)"}
    color_map = {"hydrogen": "#1f77b4", "hydrophobic": "#e65100"}
    dash_map = {"hydrogen": "6 6", "hydrophobic": "3 10"}

    bins: Dict[Tuple[int, int], int] = defaultdict(int)
    placed_boxes: List[Tuple[float, float, float, float]] = []
    bounds = (20.0, 20.0, width - 20.0, height - 20.0)
    rect_w, rect_h = 120.0, 48.0
    base_label_offset = 115.0
    arrow_len = 60.0

    for interaction in interactions:
        atom = next(a for a in atoms if a.serial == interaction.ligand_serial)
        ax, ay = _svg_coords(atom.x2d, atom.y2d, min_x, min_y, scale, margin, height)
        vx = atom.x2d - centroid_x
        vy = atom.y2d - centroid_y
        if abs(vx) < 1e-5 and abs(vy) < 1e-5:
            vx, vy = 1.0, 0.0
        norm = math.hypot(vx, vy)
        vx /= norm
        vy /= norm

        base_angle = math.atan2(vy, vx)
        bin_phi = int(round(base_angle / (math.pi / 18)))
        key = (bin_phi, interaction.ligand_serial)
        duplicate_index = bins[key]
        bins[key] += 1
        angle_shift = (duplicate_index - (bins[key] - 1) / 2.0) * math.radians(10)
        shifted_vx, shifted_vy = _rotate(vx, vy, angle_shift)
        shifted_vy *= -1

        ex = ax + shifted_vx * arrow_len
        ey = ay + shifted_vy * arrow_len
        label_center_x, label_center_y, box = _find_label_position(
            ax,
            ay,
            base_angle,
            base_label_offset,
            rect_w,
            rect_h,
            bounds,
            placed_boxes,
        )

        color, dash, stroke_width = _interaction_style(interaction.kind)
        marker = marker_map.get(interaction.kind, "url(#arrow-blue)")

        svg_lines.append(
            f'<line x1="{ax:.1f}" y1="{ay:.1f}" x2="{ex:.1f}" y2="{ey:.1f}" stroke="{color}" '
            f'stroke-width="{stroke_width}" stroke-dasharray="{dash}" marker-end="{marker}"/>'
        )
        svg_lines.append(
            f'<rect x="{box[0]:.1f}" y="{box[1]:.1f}" width="{rect_w:.1f}" height="{rect_h:.1f}" rx="8" ry="8" '
            f'fill="{color}" opacity="0.12"/>'
        )
        svg_lines.append(
            f'<text x="{label_center_x:.1f}" y="{label_center_y - 8:.1f}" font-family="Helvetica" font-size="12" '
            f'text-anchor="middle" fill="{color}" font-weight="bold">{interaction.residue_label}</text>'
        )
        svg_lines.append(
            f'<text x="{label_center_x:.1f}" y="{label_center_y + 10:.1f}" font-family="Helvetica" font-size="11" '
            f'text-anchor="middle" fill="{color}">{format_distance_block(interaction)}</text>'
        )

    title_x = width - 28
    title_y = height - 34
    subtitle_y = height - 12
    svg_lines.append(
        f'<text x="{title_x:.1f}" y="{title_y:.1f}" font-family="Helvetica" font-size="20" '
        'fill="#212121" font-weight="bold" text-anchor="end">'
        "Protein–Ligand Interaction Map</text>"
    )
    svg_lines.append(
        f'<text x="{title_x:.1f}" y="{subtitle_y:.1f}" font-family="Helvetica" font-size="12" '
        'fill="#616161" text-anchor="end">'
        "Highlights: blue = H-bonds, orange = hydrophobic contacts</text>"
    )
    svg_lines.append("</svg>")

    output_path.write_text("\n".join(svg_lines))


def main() -> None:
    args = parse_args()
    atoms, bonds = load_ligand_atoms(args.ligand_pdb, args.ligand_sdf)
    if args.style == "shaded":
        assign_projected_coordinates(atoms)
    root = ET.parse(args.report).getroot()
    interactions = aggregate_interactions(extract_interactions(root, atoms))
    draw_svg(
        atoms,
        bonds,
        interactions,
        args.output,
        scale=args.scale,
        margin=args.margin,
        style=args.style,
    )
    print(f"SVG written to {args.output}")


if __name__ == "__main__":
    main()
