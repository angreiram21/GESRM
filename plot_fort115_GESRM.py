#!/usr/bin/env python3
"""Plot GESRM rest-frame energy and baryon-density plane cuts.

Example
-------
python3 plot_fort115_GESRM.py \
  printforinhydro/fort_115_b0_0.25.dat \
  printforinhydro/fort_115_b0_0.50.dat \
  printforinhydro/fort_115_b0_0.75.dat \
  --output-dir figuras

The expected fort.115 columns are

    i  j  k  T00  N0  pressure  vz  T0z

where T00 and N0 are defined in the GESRM calculation frame.  The local
rest-frame fields are reconstructed with the same relations used by GESRM:

    e_rest = T00 / ((1 + cs2) / (1 - vz**2) - cs2)
    n_rest = N0 * sqrt(1 - vz**2)

The transverse plane is the cut at z=0.  The reaction plane is the cut at the
y-cell center nearest to y=0.  No integration is performed: plotted values
remain volume densities (GeV/fm^3 and 1/fm^3).
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Colormap, ListedColormap


@dataclass
class Grid:
    dx: float = 0.7
    dy: float = 0.7
    dz: float = 0.007
    xmin: float = -14.9262
    ymin: float = -14.9262
    nz: int = 4001


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Draw rest-frame energy and baryon-density plane cuts from "
            "GESRM fort_115.dat."
        )
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        type=Path,
        help="one or more fort_115_*.dat input files",
    )
    parser.add_argument(
        "-m",
        "--metadata",
        type=Path,
        default=Path("output/output1"),
        help="GESRM output1 file used to read the grid (default: output/output1)",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=Path("."),
        help="directory for generated figures (default: current directory)",
    )
    parser.add_argument(
        "--save-data",
        action="store_true",
        help="also save one compressed <input>_cuts.npz file per input",
    )
    parser.add_argument("--cs2", type=float, default=1.0 / 3.0, help="sound speed squared")
    parser.add_argument("--dx", type=float, help="override x cell size [fm]")
    parser.add_argument("--dy", type=float, help="override y cell size [fm]")
    parser.add_argument("--dz", type=float, help="override z cell size [fm]")
    parser.add_argument("--xmin", type=float, help="override minimum x node [fm]")
    parser.add_argument("--ymin", type=float, help="override minimum y node [fm]")
    parser.add_argument("--nz", type=int, help="override number of z nodes")
    parser.add_argument(
        "--clip-negative",
        action="store_true",
        help="replace negative reconstructed energy/density with zero",
    )
    parser.add_argument(
        "--vmax-percentile",
        type=float,
        default=100.0,
        help="upper color limit percentile, in (0,100] (default: 100)",
    )
    parser.add_argument(
        "--colormap",
        type=Path,
        default=Path("magcolor"),
        help="three-column RGB palette file (default: magcolor)",
    )
    parser.add_argument("--dpi", type=int, default=180, help="figure resolution")
    parser.add_argument("--show", action="store_true", help="open an interactive window")
    return parser.parse_args()


def read_grid(path: Path) -> Grid:
    grid = Grid()
    if not path.exists():
        print(f"Warning: metadata file {path} not found; using built-in grid defaults.")
        return grid

    text = path.read_text(encoding="utf-8", errors="replace")
    patterns = {
        "dx": r"\bdx=\s*([-+0-9.Ee]+)",
        "dy": r"\bdy=\s*([-+0-9.Ee]+)",
        "dz": r"\bdz=\s*([-+0-9.Ee]+)",
        "xmin": r"\bxmin=\s*([-+0-9.Ee]+)",
        "ymin": r"\bymin=\s*([-+0-9.Ee]+)",
        "nz": r"\bnz=\s*([0-9]+)",
    }
    for name, pattern in patterns.items():
        match = re.search(pattern, text)
        if match:
            value = int(match.group(1)) if name == "nz" else float(match.group(1))
            setattr(grid, name, value)
    return grid


def apply_overrides(grid: Grid, args: argparse.Namespace) -> Grid:
    for name in ("dx", "dy", "dz", "xmin", "ymin", "nz"):
        value = getattr(args, name)
        if value is not None:
            setattr(grid, name, value)
    if min(grid.dx, grid.dy, grid.dz) <= 0 or grid.nz < 2:
        raise ValueError("Grid spacings must be positive and nz must be at least 2.")
    return grid


def load_fort115(path: Path) -> np.ndarray:
    data = np.loadtxt(path, dtype=np.float64)
    if data.ndim != 2 or data.shape[1] != 8:
        raise ValueError(f"Expected 8 columns in {path}; found shape {data.shape}.")
    if not np.isfinite(data).all():
        raise ValueError("Input contains NaN or infinite values.")
    indices = data[:, :3]
    if not np.equal(indices, np.rint(indices)).all():
        raise ValueError("The first three columns must contain integer grid indices.")
    return data


def reconstruct_rest_frame(
    t00: np.ndarray, n0: np.ndarray, vz: np.ndarray, cs2: float
) -> tuple[np.ndarray, np.ndarray]:
    if not 0.0 < cs2 < 1.0:
        raise ValueError("--cs2 must be between 0 and 1.")
    mask = np.abs(vz) >= 1.0
    if np.any(mask):
        vz = vz.copy()
        vz[mask] = np.sign(vz[mask]) * 0.999
    one_minus_v2 = 1.0 - vz * vz
    factor = (1.0 + cs2) / one_minus_v2 - cs2
    return t00 / factor, n0 * np.sqrt(one_minus_v2)


def make_plane_cuts(
    i: np.ndarray,
    j: np.ndarray,
    k: np.ndarray,
    values: np.ndarray,
    grid: Grid,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int, int]:
    imin, imax = int(i.min()), int(i.max())
    jmin, jmax = int(j.min()), int(j.max())
    kmin, kmax = int(k.min()), int(k.max())
    ni, nj, nk = imax - imin + 1, jmax - jmin + 1, kmax - kmin + 1

    # i,j label transverse cells; k follows the z-node convention in GESRM.f.
    x = grid.xmin + (np.arange(imin, imax + 1) - 0.5) * grid.dx
    y = grid.ymin + (np.arange(jmin, jmax + 1) - 0.5) * grid.dy
    z0 = -0.5 * (grid.nz - 1) * grid.dz
    z = z0 + (np.arange(kmin, kmax + 1) - 1.0) * grid.dz

    target_k = int(round(1.0 - z0 / grid.dz))
    target_j = jmin + int(np.argmin(np.abs(y)))
    if not kmin <= target_k <= kmax:
        raise ValueError(
            f"The z=0 index k={target_k} is outside the stored range {kmin}..{kmax}."
        )

    xy = np.zeros((ni, nj), dtype=np.float64)
    mask_xy = k == target_k
    xy[i[mask_xy] - imin, j[mask_xy] - jmin] = values[mask_xy]

    xz = np.zeros((ni, nk), dtype=np.float64)
    mask_xz = j == target_j
    xz[i[mask_xz] - imin, k[mask_xz] - kmin] = values[mask_xz]
    return x, y, z, xy, xz, target_k, target_j


def color_limits(array: np.ndarray, percentile: float) -> tuple[float, float]:
    finite = array[np.isfinite(array)]
    if finite.size == 0:
        return 0.0, 1.0
    vmin = min(0.0, float(finite.min()))
    vmax = float(np.percentile(finite, percentile))
    if vmax <= vmin:
        vmax = vmin + 1.0
    return vmin, vmax


def load_colormap(path: Path) -> ListedColormap:
    colors = np.loadtxt(path, dtype=np.float64)
    if colors.ndim != 2 or colors.shape[1] != 3:
        raise ValueError(f"Expected three RGB columns in {path}; found shape {colors.shape}.")
    if not np.isfinite(colors).all() or np.any((colors < 0.0) | (colors > 1.0)):
        raise ValueError(f"RGB values in {path} must be finite and lie between 0 and 1.")
    return ListedColormap(colors, name=path.stem or "magcolor")


def draw_panel(
    ax: plt.Axes,
    image: np.ndarray,
    extent: tuple[float, float, float, float],
    xlabel: str,
    ylabel: str,
    cbar_label: str,
    cmap: Colormap,
    percentile: float,
) -> None:
    vmin, vmax = color_limits(image, percentile)
    artist = ax.imshow(
        image,
        origin="lower",
        extent=extent,
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
    )
    ax.set(xlabel=xlabel, ylabel=ylabel)
    ax.text(
        0.03,
        0.96,
        "GESRM",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=18,
        color="black",
    )
    plt.colorbar(artist, ax=ax, label=cbar_label, pad=0.02)


def save_figure(
    image: np.ndarray,
    extent: tuple[float, float, float, float],
    xlabel: str,
    ylabel: str,
    cbar_label: str,
    colormap: Colormap,
    percentile: float,
    output: Path,
    dpi: int,
    xlim: tuple[float, float] | None = None,
) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(6, 5), constrained_layout=True)
    draw_panel(
        ax, image, extent, xlabel, ylabel, cbar_label, colormap, percentile
    )
    if xlim is not None:
        ax.set_xlim(*xlim)
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    print(f"Saved figure: {output}")
    return fig


def process_input(
    input_path: Path,
    args: argparse.Namespace,
    grid: Grid,
    colormap: Colormap,
) -> list[plt.Figure]:
    print(f"Processing: {input_path}")
    data = load_fort115(input_path)
    i, j, k = (data[:, column].astype(np.int64) for column in range(3))
    t00, n0, vz = data[:, 3], data[:, 4], data[:, 6]
    energy, baryon = reconstruct_rest_frame(t00, n0, vz, args.cs2)

    negative_energy = int(np.count_nonzero(energy < 0.0))
    negative_baryon = int(np.count_nonzero(baryon < 0.0))
    if negative_energy or negative_baryon:
        print(
            "Warning: reconstructed negative cells: "
            f"energy={negative_energy}, baryon={negative_baryon}."
        )
    if args.clip_negative:
        energy = np.maximum(energy, 0.0)
        baryon = np.maximum(baryon, 0.0)

    x, y, z, energy_xy, energy_xz, target_k, target_j = make_plane_cuts(
        i, j, k, energy, grid
    )
    _, _, _, baryon_xy, baryon_xz, _, _ = make_plane_cuts(
        i, j, k, baryon, grid
    )
    selected_y = grid.ymin + (target_j - 0.5) * grid.dy
    print(f"Transverse cut: k={target_k}, z=0 fm")
    print(f"Reaction-plane cut: j={target_j}, y={selected_y:.6g} fm (nearest to zero)")

    x_extent = (x[0] - grid.dx / 2.0, x[-1] + grid.dx / 2.0)
    y_extent = (y[0] - grid.dy / 2.0, y[-1] + grid.dy / 2.0)
    z_extent = (z[0] - grid.dz / 2.0, z[-1] + grid.dz / 2.0)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    stem = input_path.stem
    figures = [
        save_figure(
            energy_xy,
            (*y_extent, *x_extent),
            r"$y$ [fm]",
            r"$x$ [fm]",
            r"$e$ [GeV/fm$^3$]",
            colormap,
            args.vmax_percentile,
            args.output_dir / f"GESRM_{stem}_energy_xy.png",
            args.dpi,
        ),
        save_figure(
            baryon_xy,
            (*y_extent, *x_extent),
            r"$y$ [fm]",
            r"$x$ [fm]",
            r"$n_B$ [fm$^{-3}$]",
            colormap,
            args.vmax_percentile,
            args.output_dir / f"GESRM_{stem}_baryon_xy.png",
            args.dpi,
        ),
        save_figure(
            energy_xz,
            (*z_extent, *x_extent),
            r"$z$ [fm]",
            r"$x$ [fm]",
            r"$e$ [GeV/fm$^3$]",
            colormap,
            args.vmax_percentile,
            args.output_dir / f"GESRM_{stem}_energy_xz.png",
            args.dpi,
            xlim=(-7.5, 7.5),
        ),
        save_figure(
            baryon_xz,
            (*z_extent, *x_extent),
            r"$z$ [fm]",
            r"$x$ [fm]",
            r"$n_B$ [fm$^{-3}$]",
            colormap,
            args.vmax_percentile,
            args.output_dir / f"GESRM_{stem}_baryon_xz.png",
            args.dpi,
            xlim=(-7.5, 7.5),
        ),
    ]

    if args.save_data:
        data_output = args.output_dir / f"{stem}_cuts.npz"
        np.savez_compressed(
            data_output,
            x=x,
            y=y,
            z=z,
            energy_z0=energy_xy,
            baryon_z0=baryon_xy,
            energy_y_nearest_zero=energy_xz,
            baryon_y_nearest_zero=baryon_xz,
            selected_k=target_k,
            selected_j=target_j,
            selected_z=0.0,
            selected_y=selected_y,
            dx=grid.dx,
            dy=grid.dy,
            dz=grid.dz,
            cs2=args.cs2,
        )
        print(f"Saved cut data: {data_output}")
    return figures


def main() -> None:
    args = parse_args()
    if not 0.0 < args.vmax_percentile <= 100.0:
        raise ValueError("--vmax-percentile must lie in (0, 100].")

    grid = apply_overrides(read_grid(args.metadata), args)
    colormap = load_colormap(args.colormap)
    plt.rcParams.update({"font.size": 18})

    figures: list[plt.Figure] = []
    for input_path in args.inputs:
        figures.extend(process_input(input_path, args, grid, colormap))

    if args.show:
        plt.show()
    else:
        for figure in figures:
            plt.close(figure)


if __name__ == "__main__":
    main()
