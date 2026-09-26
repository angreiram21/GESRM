# GESRM

**Generalized Effective String Rope Model**

Fortran implementation of the **Generalized Effective String Rope Model (GESRM)** for event-by-event simulations of high-energy nuclear collisions.

---

## Authors

**V. Magas & A. Reina**

Barcelona, **April 2022**

---

## Description

The **Generalized Effective String Rope Model (GESRM)** is based on the **Effective String Rope Model (ESRM)** and has been generalized to perform **event-by-event simulations of high-energy nuclear collisions with initial-state fluctuations generated using a Glauber Monte Carlo approach**.

This repository contains the Fortran implementation of the GESRM model. The code calculates the resulting initial state on a three-dimensional spatial grid, providing physical quantities such as:

* Energy density
* Baryon density
* Longitudinal momentum
* Rapidity
* Pressure
* Vorticity

The resulting initial conditions can be used as input for subsequent hydrodynamical calculations.

The main program is contained in:

```text
GESRM.f
```

The code can also generate a `fort.100` file containing initial-condition information for hydrodynamical calculations.

---

## Scientific model

The **Generalized Effective String Rope Model (GESRM)** is based on the Effective String Rope Model (ESRM). It generalizes the ESRM framework to allow for **event-by-event simulations and initial-state fluctuations**.

The nuclear geometry and collision configuration are generated using a **Glauber Monte Carlo approach**, allowing different initial configurations to be obtained for individual collision events.

The GESRM model and its application to the initial stages of ultrarelativistic heavy-ion collisions are described in:

> **A. Reina Ramirez, V. K. Magas, L. P. Csernai, and D. Strottman**,
> *Generalized Effective String Rope Model for the Initial Stages of Ultrarelativistic Heavy Ion Collisions*,
> **Physical Review C 107, 034915 (2023)**.
> DOI: `10.1103/PhysRevC.107.034915`
> arXiv: `2205.15797`

---

## Repository structure

The main files and directories are:

```text
GESRM/
│
├── GESRM.f
├── GESRM_vorticity.py
├── plot_fort115_GESRM.py
│
├── input/
│   └── global_parameters_GESRM_code
│
├── output/
│   └── output1
│
├── printforinhydro/
│   └── fort_115_*.dat
│
└── magcolor/
```

The contents of the output directories depend on the simulations performed.

---

## Requirements

### Fortran

The main program is written in **Fortran** and uses fixed-format source code together with legacy Fortran features.

A Fortran compiler such as **GNU Fortran (`gfortran`)** is required.

On Debian/Ubuntu-based systems:

```bash
sudo apt install gfortran
```

### Python

The repository also contains Python scripts for analysis and visualisation.

The plotting utilities require:

* Python 3
* NumPy
* Matplotlib

Install them with:

```bash
pip install numpy matplotlib
```

---

## Compilation

From the root directory of the repository, compile the main program with:

```bash
gfortran -O2 -std=legacy -ffixed-line-length-none GESRM.f -o GESRM
```

This creates the executable:

```text
GESRM
```

The `-std=legacy` option allows compilation of legacy Fortran constructs, while `-ffixed-line-length-none` avoids problems with long lines in fixed-format source code.

---

## Running the program

The program uses relative paths to access the input and output directories. Therefore, it should be executed from the **root directory of the repository**.

Run:

```bash
./GESRM
```

The main input file is:

```text
input/global_parameters_GESRM_code
```

The main program output is written to:

```text
output/output1
```

Depending on the configuration, additional files required for hydrodynamical calculations are generated.

---

## Input parameters

The main input file,

```text
input/global_parameters_GESRM_code
```

contains the physical and numerical parameters used in the simulation.

The main parameters read by the program include:

| Parameter   | Description                                  |
| ----------- | -------------------------------------------- |
| `epsilon00` | Initial energy per nucleon                   |
| `b0`        | Reduced impact parameter                     |
| `A1`        | Mass number of projectile nucleus            |
| `A2`        | Mass number of target nucleus                |
| `r1`        | Projectile nuclear radius                    |
| `r2`        | Target nuclear radius                        |
| `aws1`      | Projectile surface/diffuseness parameter     |
| `aws2`      | Target surface/diffuseness parameter         |
| `Mn`        | Nucleon mass                                 |
| `rn`        | Nucleon radius                               |
| `A`         | Nuclear mass-number parameter                |
| `rho00`     | Initial nuclear density                      |
| `tnc`       | Number of transverse grid cells              |
| `nteventos` | Number of simulated events                   |
| `vortin`    | Enables/disables the vorticity calculation   |
| `for100`    | Controls generation of the `fort.100` output |
| `check`     | Check/debug option                           |
| `end_time`  | End-time parameter                           |

The program reads these parameters at the beginning of the calculation and reports them in the output.

---

## Event-by-event simulation

The number of events is controlled by:

```text
nteventos
```

For each event, the code performs a Glauber Monte Carlo calculation to generate the corresponding initial nuclear configuration.

The resulting configuration is then used to calculate the physical properties of the initial state on a three-dimensional spatial grid.

Among the quantities calculated are:

* Energy density
* Baryon density
* Momentum
* Rapidity
* Pressure
* Vorticity

The transverse grid spacing is determined from the nuclear radii and the value of `tnc`.

The longitudinal grid spacing is smaller than the transverse spacing, according to the definitions implemented in `GESRM.f`.

---

## Vorticity

GESRM includes an option for calculating vorticity.

The parameter

```text
vortin
```

controls whether the vorticity calculation is enabled.

The selected configuration is reported when the program starts.

Additional Python functionality related to vorticity is provided in:

```text
GESRM_vorticity.py
```

---

## Output

One of the main output files produced by the program is:

```text
fort.100
```

This file contains initial-state information that can be used as input for subsequent hydrodynamical calculations.

The program also writes simulation information and parameters to:

```text
output/output1
```

Additional files containing the calculated fields can be generated for further analysis.

In particular, `fort_115` data files can be used with the included plotting script.

---

## Visualisation

The repository contains:

```text
plot_fort115_GESRM.py
```

This Python script can be used to visualise the information contained in `fort_115` files.

The script processes the fields stored in the data and reconstructs rest-frame quantities such as the energy density and baryon density.

It produces four types of plots:

1. Energy density in the transverse `xy` plane
2. Baryon density in the transverse `xy` plane
3. Energy density in the reaction `xz` plane
4. Baryon density in the reaction `xz` plane

For example:

```bash
python3 plot_fort115_GESRM.py \
    printforinhydro/fort_115_b0_0.25.dat \
    --output-dir figuras
```

The generated figures are saved in:

```text
figuras/
```

The script provides additional options for controlling the calculation and visualisation.

To see all available options:

```bash
python3 plot_fort115_GESRM.py --help
```

---

## Typical workflow

A typical calculation can be performed as follows.

### 1. Clone the repository

```bash
git clone https://github.com/angreiram21/GESRM.git
cd GESRM
```

### 2. Configure the simulation

Edit the main input file:

```text
input/global_parameters_GESRM_code
```

and specify the desired collision and simulation parameters.

### 3. Compile the code

```bash
gfortran -O2 -std=legacy -ffixed-line-length-none GESRM.f -o GESRM
```

### 4. Run the simulation

```bash
./GESRM
```

### 5. Inspect the output

Check:

```text
output/output1
```

and the generated initial-condition files.

### 6. Visualise the results

If `fort_115` data are available:

```bash
python3 plot_fort115_GESRM.py \
    <fort_115_file> \
    --output-dir figuras
```

---

## Notes

### Relative paths

The program expects the `input/` and `output/` directories to be available relative to the directory from which the executable is launched.

For this reason, it is recommended to run the executable from the repository root:

```bash
./GESRM
```

### Fortran source format

`GESRM.f` is written in fixed-format Fortran and uses legacy language features.

Depending on the compiler and version, the compilation options may need to be adjusted.

### Reproducibility

For reproducible calculations, it is recommended to record:

* GESRM version or Git commit
* Input parameters
* Collision system
* Collision energy
* Impact parameter
* Number of simulated events
* Compiler and compiler version
* Compilation options

---

## References

### GESRM

**A. Reina Ramirez, V. K. Magas, L. P. Csernai, and D. Strottman**,
*Generalized Effective String Rope Model for the Initial Stages of Ultrarelativistic Heavy Ion Collisions*,
**Physical Review C 107, 034915 (2023)**.

DOI: `10.1103/PhysRevC.107.034915`
arXiv: `2205.15797`

### Effective String Rope Model

**V. K. Magas, L. P. Csernai, and D. D. Strottman**,
*The Effective String Rope Model*,
**Physical Review C 64, 014901 (2001)**.

### Glauber Monte Carlo

**M. L. Miller, K. Reygers, S. J. Sanders, and P. Steinberg**,
*Glauber Modeling in High Energy Nuclear Collisions*,
**Annual Review of Nuclear and Particle Science 57, 205–243 (2007)**.

---

## Citation

If you use **GESRM** in a scientific publication, please cite both the scientific paper describing the model and the software itself.

### Scientific paper

```bibtex
@article{ReinaRamirez2023,
  author        = {Reina Ramirez, Angel and Magas, V. K. and Csernai, L. P. and Strottman, D.},
  title         = {Generalized Effective String Rope Model for the Initial Stages of Ultrarelativistic Heavy Ion Collisions},
  journal       = {Physical Review C},
  volume        = {107},
  pages         = {034915},
  year          = {2023},
  doi           = {10.1103/PhysRevC.107.034915},
  eprint        = {2205.15797},
  archivePrefix = {arXiv}
}
```

### GESRM software

```bibtex
@software{GESRM2022,
  author    = {Magas, V. and Reina, Angel},
  title     = {GESRM: Generalized Effective String Rope Model},
  year      = {2022},
  month     = {4},
  publisher = {GitHub},
  url       = {https://github.com/angreiram21/GESRM},
  note      = {Fortran implementation of the Generalized Effective String Rope Model}
}
```

The software citation refers to the version of the GESRM code developed in **Barcelona in April 2022**, as indicated in the source code.

For reproducibility, users are encouraged to cite the specific Git commit or release used for their calculations whenever possible.

---

## License

GESRM is free software: you can redistribute it and/or modify it under the terms of the **GNU General Public License v3.0 (GPL-3.0)**.

See the `LICENSE` file for the full license text.

---

## Repository

The GESRM source code is available at:

**https://github.com/angreiram21/GESRM**
