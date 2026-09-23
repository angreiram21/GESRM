#!/usr/bin/env python3

import argparse
from pathlib import Path
import math

import numpy as np
import matplotlib.pyplot as plt


# =============================================================================
# GRID PARAMETERS
# =============================================================================

DX = 0.638  # fm
DY = 0.638  # fm
DZ = 0.638  # fm


# =============================================================================
# FUNCTIONS
# =============================================================================

def gamma_factor(v):
    """Lorentz factor."""
    return 1.0 / math.sqrt(1.0 - v * v)


def centered_coordinates(index_min, index_max, spacing):
    """
    Same coordinate convention used in GESRM.m:

        x = (i - (imax/2 + imin/2)) * dx
    """
    indices = np.arange(index_min, index_max + 1, dtype=float)

    return (
        indices
        - 0.5 * (index_max + index_min)
    ) * spacing


def find_zero_layer(coordinates, name):
    """
    Find explicitly the layer whose physical coordinate is zero.
    """

    idx = np.where(
        np.isclose(
            coordinates,
            0.0,
            atol=1.0e-12,
            rtol=0.0
        )
    )[0]

    if len(idx) != 1:
        raise RuntimeError(
            f"Expected exactly one {name}=0 layer, "
            f"but found {len(idx)}.\n"
            f"{name} coordinates:\n{coordinates}"
        )

    return int(idx[0])


# =============================================================================
# READ fort_100.dat
# =============================================================================

def read_input(filename):

    data = np.loadtxt(filename)

    if data.ndim == 1:
        data = data[np.newaxis, :]

    if data.shape[1] < 7:
        raise ValueError(
            "fort_100.dat must contain at least seven columns."
        )

    return data


# =============================================================================
# BUILD THE SAME ARRAYS AS IN GESRM.m
# =============================================================================

def build_fields(data):

    # Original GESRM grid indices
    i_input = data[:, 0].astype(int)
    j_input = data[:, 1].astype(int)
    k_input = data[:, 2].astype(int)

    imin = i_input.min()
    imax = i_input.max()

    jmin = j_input.min()
    jmax = j_input.max()

    kmin = k_input.min()
    kmax = k_input.max()

    # Extra space is kept because the MATLAB code stores the
    # original cell at (i+1,k+1,j+1).
    shape = (
        imax + 2,
        kmax + 2,
        jmax + 2
    )

    EDCF = np.zeros(shape)
    Vz   = np.zeros(shape)

    for row in data:

        i0 = int(row[0])
        j0 = int(row[1])
        k0 = int(row[2])

        # Same columns as GESRM.m
        T00 = row[3]   # column 4
        vz  = row[6]   # column 7

        # Same prescription as MATLAB
        if vz == 1.0:
            vz = 0.99

        if vz == -1.0:
            vz = -0.99

        if abs(vz) >= 1.0:
            raise ValueError(
                f"|vz| >= 1 at cell "
                f"(i,j,k)=({i0},{j0},{k0})"
            )

        # Same +1 shift used in GESRM.m
        ii = i0 + 1
        jj = j0 + 1
        kk = k0 + 1

        EDCF[ii, kk, jj] = T00
        Vz[ii, kk, jj]   = vz

    bounds = {
        "imin": int(imin),
        "imax": int(imax),
        "jmin": int(jmin),
        "jmax": int(jmax),
        "kmin": int(kmin),
        "kmax": int(kmax),
    }

    return EDCF, Vz, bounds


# =============================================================================
# VORTICITY CALCULATION
# =============================================================================

def calculate_vorticity(EDCF, Vz, bounds):

    imin = bounds["imin"]
    imax = bounds["imax"]

    jmin = bounds["jmin"]
    jmax = bounds["jmax"]

    kmin = bounds["kmin"]
    kmax = bounds["kmax"]

    classvort = np.zeros_like(Vz)
    relvort   = np.zeros_like(Vz)

    # Same total quantities as MATLAB
    tcv = 0.0
    trv = 0.0

    for iii in range(imin, imax + 1):

        for iy in range(jmin, jmax + 1):

            for kkk in range(kmin, kmax + 1):

                # MATLAB:
                #
                # if (Vz(iii,kkk,iy) ~= 0)

                if Vz[iii, kkk, iy] == 0.0:
                    continue


                # =============================================================
                # DIAGONAL NEIGHBOURS
                # =============================================================

                # -------------------------------------------------------------
                # v12
                # -------------------------------------------------------------

                if EDCF[iii - 1, kkk + 1, iy] == 0.0:

                    v12 = Vz[iii, kkk + 1, iy]
                    gamma12 = gamma_factor(v12)
                    w12 = 0

                else:

                    v12 = Vz[iii - 1, kkk + 1, iy]
                    gamma12 = gamma_factor(v12)
                    w12 = 1


                # -------------------------------------------------------------
                # v21
                # -------------------------------------------------------------

                if EDCF[iii + 1, kkk - 1, iy] == 0.0:

                    v21 = Vz[iii, kkk - 1, iy]
                    gamma21 = gamma_factor(v21)
                    w21 = 0

                else:

                    v21 = Vz[iii + 1, kkk - 1, iy]
                    gamma21 = gamma_factor(v21)
                    w21 = 1


                # -------------------------------------------------------------
                # v22
                # -------------------------------------------------------------

                if EDCF[iii + 1, kkk + 1, iy] == 0.0:

                    v22 = Vz[iii, kkk + 1, iy]
                    gamma22 = gamma_factor(v22)
                    w22 = 0

                else:

                    v22 = Vz[iii + 1, kkk + 1, iy]
                    gamma22 = gamma_factor(v22)
                    w22 = 1


                # -------------------------------------------------------------
                # v11
                # -------------------------------------------------------------

                if EDCF[iii - 1, kkk - 1, iy] == 0.0:

                    v11 = Vz[iii, kkk - 1, iy]
                    gamma11 = gamma_factor(v11)
                    w11 = 0

                else:

                    v11 = Vz[iii - 1, kkk - 1, iy]
                    gamma11 = gamma_factor(v11)
                    w11 = 1


                # =============================================================
                # NEAREST NEIGHBOURS IN x
                # =============================================================

                # -------------------------------------------------------------
                # v10
                # -------------------------------------------------------------

                if EDCF[iii - 1, kkk, iy] == 0.0:

                    v10 = Vz[iii, kkk, iy]
                    gamma10 = gamma_factor(v10)
                    w10 = 0

                else:

                    v10 = Vz[iii - 1, kkk, iy]
                    gamma10 = gamma_factor(v10)
                    w10 = 1


                # -------------------------------------------------------------
                # v20
                # -------------------------------------------------------------

                if EDCF[iii + 1, kkk, iy] == 0.0:

                    v20 = Vz[iii, kkk, iy]
                    gamma20 = gamma_factor(v20)
                    w20 = 0

                else:

                    v20 = Vz[iii + 1, kkk, iy]
                    gamma20 = gamma_factor(v20)
                    w20 = 1


                # =============================================================
                # BOUNDARY CONDITIONS
                #
                # Exact logic used in GESRM.m
                # =============================================================

                if (v22 == 0.0) and (w22 == 0):

                    v22 = Vz[iii - 1, kkk + 1, iy]
                    gamma22 = gamma_factor(v22)
                    w22 = -1


                if (v12 == 0.0) and (w12 == 0):

                    v12 = Vz[iii + 1, kkk + 1, iy]
                    gamma12 = gamma_factor(v12)
                    w12 = -1


                if (v21 == 0.0) and (w21 == 0):

                    v21 = Vz[iii - 1, kkk - 1, iy]
                    gamma21 = gamma_factor(v21)
                    w21 = -1


                if (v11 == 0.0) and (w11 == 0):

                    v11 = Vz[iii + 1, kkk - 1, iy]
                    gamma11 = gamma_factor(v11)
                    w11 = -1


                # =============================================================
                # WEIGHTS
                # =============================================================

                q1 = w20 + w10
                q2 = w22 + w11 + w21 + w12

                if q1 == 0:
                    q1 = 1

                if q2 == 0:
                    q2 = 1


                # =============================================================
                # CENTRAL CELL
                # =============================================================

                v00 = Vz[iii, kkk, iy]

                gamma00 = gamma_factor(v00)


                # =============================================================
                # DISCRETE DERIVATIVES
                #
                # Exactly the combinations entering GESRM.m
                # =============================================================

                dvz_dx = (
                    (v20 - v10) / (q1 * DX)
                    +
                    (
                        v22
                        - v11
                        + v21
                        - v12
                    ) / (q2 * DX)
                )


                dgamma_dx = (
                    (gamma20 - gamma10) / (q1 * DX)
                    +
                    (
                        gamma22
                        - gamma11
                        + gamma21
                        - gamma12
                    ) / (q2 * DX)
                )


                # =============================================================
                # CLASSICAL VORTICITY
                # =============================================================

                classvort[iii, kkk, iy] = (
                    -0.5 * dvz_dx
                )


                # =============================================================
                # RELATIVISTIC VORTICITY
                #
                # Sign convention adopted in the thesis and in GESRM.m:
                #
                # omega_y^rel =
                #
                # -1/2 [ gamma d(vz)/dx + vz d(gamma)/dx ]
                # =============================================================

                relvort[iii, kkk, iy] = (
                    -0.5
                    *
                    (
                        gamma00 * dvz_dx
                        +
                        v00 * dgamma_dx
                    )
                )


                # =============================================================
                # TOTAL VORTICITY
                # =============================================================

                tcv += classvort[iii, kkk, iy]
                trv += relvort[iii, kkk, iy]


    # Same final operation appearing in GESRM.m
    tcv_matlab = 100.0 * tcv
    trv_matlab = 100.0 * trv

    return (
        classvort,
        relvort,
        tcv,
        trv,
        tcv_matlab,
        trv_matlab,
    )


# =============================================================================
# EXTRACT DISTRIBUTIONS
# =============================================================================

def extract_maps(classvort, relvort, bounds):

    imin = bounds["imin"]
    imax = bounds["imax"]

    jmin = bounds["jmin"]
    jmax = bounds["jmax"]

    kmin = bounds["kmin"]
    kmax = bounds["kmax"]


    # Physical coordinates
    x = centered_coordinates(
        imin,
        imax,
        DX
    )

    y = centered_coordinates(
        jmin,
        jmax,
        DY
    )

    z = centered_coordinates(
        kmin,
        kmax,
        DZ
    )


    # Compact 3D arrays:
    #
    # axis 0 -> x
    # axis 1 -> z
    # axis 2 -> y

    class_xzy = classvort[
        imin:imax + 1,
        kmin:kmax + 1,
        jmin:jmax + 1
    ]

    rel_xzy = relvort[
        imin:imax + 1,
        kmin:kmax + 1,
        jmin:jmax + 1
    ]


    # =========================================================================
    # REACTION PLANE: y = 0
    # =========================================================================

    iy0 = find_zero_layer(
        y,
        "y"
    )

    omega_class_xz = class_xzy[:, :, iy0]
    omega_rel_xz   = rel_xzy[:, :, iy0]


    # =========================================================================
    # SUM OVER ALL y-LAYERS
    #
    # These correspond to the Omega_y distributions of Figs. 30 and 32.
    #
    # IMPORTANT:
    # this is a discrete sum, not an integral.
    # Therefore there is NO multiplication by DY.
    # =========================================================================

    Omega_class_xz = np.sum(
        class_xzy,
        axis=2
    )

    Omega_rel_xz = np.sum(
        rel_xzy,
        axis=2
    )


    return {
        "x": x,
        "y": y,
        "z": z,

        "iy0": iy0,

        "class_xzy": class_xzy,
        "rel_xzy": rel_xzy,

        "omega_class_xz": omega_class_xz,
        "omega_rel_xz": omega_rel_xz,

        "Omega_class_xz": Omega_class_xz,
        "Omega_rel_xz": Omega_rel_xz,
    }


# =============================================================================
# PLOTTING
# =============================================================================

def plot_xz(
    x,
    z,
    values,
    title,
    colorbar_label,
    filename
):

    # Symmetric scale around zero
    vmax = np.nanmax(
        np.abs(values)
    )

    if vmax == 0.0:
        vmax = 1.0


    fig, ax = plt.subplots(
        figsize=(7, 6)
    )


    # Paper orientation:
    #
    # horizontal -> z
    # vertical   -> x

    mesh = ax.pcolormesh(
        z,
        x,
        values,
        shading="auto",
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax
    )


    ax.set_xlabel(
        r"$z$ (fm)",
        fontsize=14
    )

    ax.set_ylabel(
        r"$x$ (fm)",
        fontsize=14
    )

    ax.set_title(
        title,
        fontsize=14
    )

    ax.set_aspect(
        "equal",
        adjustable="box"
    )


    cbar = fig.colorbar(
        mesh,
        ax=ax
    )

    cbar.set_label(
        colorbar_label,
        fontsize=13
    )


    ax.tick_params(
        labelsize=12
    )

    fig.tight_layout()

    fig.savefig(
        filename,
        bbox_inches="tight"
    )

    plt.close(fig)


# =============================================================================
# MAIN
# =============================================================================

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "input",
        nargs="?",
        default="fort_100.dat"
    )

    parser.add_argument(
        "--output-dir",
        default="vorticity_output"
    )

    args = parser.parse_args()


    input_file = Path(
        args.input
    )

    output_dir = Path(
        args.output_dir
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )


    # =========================================================================
    # READ INPUT
    # =========================================================================

    data = read_input(
        input_file
    )


    # =========================================================================
    # BUILD FIELDS
    # =========================================================================

    EDCF, Vz, bounds = build_fields(
        data
    )


    # =========================================================================
    # CALCULATE VORTICITY
    # =========================================================================

    (
        classvort,
        relvort,
        tcv,
        trv,
        tcv_matlab,
        trv_matlab
    ) = calculate_vorticity(
        EDCF,
        Vz,
        bounds
    )


    # =========================================================================
    # EXTRACT MAPS
    # =========================================================================

    maps = extract_maps(
        classvort,
        relvort,
        bounds
    )


    x = maps["x"]
    y = maps["y"]
    z = maps["z"]


    omega_class_xz = maps[
        "omega_class_xz"
    ]

    omega_rel_xz = maps[
        "omega_rel_xz"
    ]

    Omega_class_xz = maps[
        "Omega_class_xz"
    ]

    Omega_rel_xz = maps[
        "Omega_rel_xz"
    ]


    # =========================================================================
    # SAVE NUMERICAL DATA
    # =========================================================================

    np.savez_compressed(
        output_dir / "GESRM_vorticity_maps.npz",

        x=x,
        y=y,
        z=z,

        class_xzy=maps["class_xzy"],
        rel_xzy=maps["rel_xzy"],

        omega_class_xz=omega_class_xz,
        omega_rel_xz=omega_rel_xz,

        Omega_class_xz=Omega_class_xz,
        Omega_rel_xz=Omega_rel_xz
    )


    np.savetxt(
        output_dir / "omega_class_xz_y0.dat",
        omega_class_xz
    )

    np.savetxt(
        output_dir / "omega_rel_xz_y0.dat",
        omega_rel_xz
    )

    np.savetxt(
        output_dir / "Omega_class_xz_sum_y.dat",
        Omega_class_xz
    )

    np.savetxt(
        output_dir / "Omega_rel_xz_sum_y.dat",
        Omega_rel_xz
    )


    # =========================================================================
    # FIG. 29 TYPE
    # =========================================================================

    plot_xz(
        x,
        z,
        omega_class_xz,

        r"$\omega_y^{\mathrm{class}}(x,z)$ at $y=0$",

        r"$\omega_y^{\mathrm{class}}$ (fm$^{-1}$)",

        output_dir
        / "omega_class_xz_y0.pdf"
    )


    # =========================================================================
    # FIG. 31 TYPE
    # =========================================================================

    plot_xz(
        x,
        z,
        omega_rel_xz,

        r"$\omega_y^{\mathrm{rel}}(x,z)$ at $y=0$",

        r"$\omega_y^{\mathrm{rel}}$ (fm$^{-1}$)",

        output_dir
        / "omega_rel_xz_y0.pdf"
    )


    # =========================================================================
    # FIG. 30 TYPE
    # =========================================================================

    plot_xz(
        x,
        z,
        Omega_class_xz,

        (
            r"$\Omega_y^{\mathrm{class}}(x,z)"
            r"=\sum_y\omega_y^{\mathrm{class}}$"
        ),

        r"$\Omega_y^{\mathrm{class}}$ (fm$^{-1}$)",

        output_dir
        / "Omega_class_xz_sum_y.pdf"
    )


    # =========================================================================
    # FIG. 32 TYPE
    # =========================================================================

    plot_xz(
        x,
        z,
        Omega_rel_xz,

        (
            r"$\Omega_y^{\mathrm{rel}}(x,z)"
            r"=\sum_y\omega_y^{\mathrm{rel}}$"
        ),

        r"$\Omega_y^{\mathrm{rel}}$ (fm$^{-1}$)",

        output_dir
        / "Omega_rel_xz_sum_y.pdf"
    )


    # =========================================================================
    # DIAGNOSTICS
    # =========================================================================

    iy0 = maps["iy0"]

    print()
    print("GESRM vorticity calculation completed.")
    print()

    print(
        "Input file:",
        input_file
    )

    print(
        "y=0 layer:",
        iy0
    )

    print(
        "Physical y coordinate:",
        y[iy0],
        "fm"
    )

    print()

    print(
        "Total classical vorticity:",
        tcv
    )

    print(
        "Total relativistic vorticity:",
        trv
    )

    print()

    print(
        "MATLAB tcv = 100 * total:",
        tcv_matlab
    )

    print(
        "MATLAB trv = 100 * total:",
        trv_matlab
    )

    print()

    print(
        "max |omega_class(y=0)| =",
        np.max(
            np.abs(
                omega_class_xz
            )
        )
    )

    print(
        "max |omega_rel(y=0)| =",
        np.max(
            np.abs(
                omega_rel_xz
            )
        )
    )

    print(
        "max |Omega_class| =",
        np.max(
            np.abs(
                Omega_class_xz
            )
        )
    )

    print(
        "max |Omega_rel| =",
        np.max(
            np.abs(
                Omega_rel_xz
            )
        )
    )


if __name__ == "__main__":
    main()
