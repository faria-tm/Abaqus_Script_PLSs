# Model-Creator: Parametric Lattice Unit Cell Generator (Abaqus/CAE)

An Abaqus/CAE Python script that generates and simulates a combinatorial family of 3D-printable lattice unit cells. Each cell is defined by a binary "fingerprint" that encodes which of 20 wall positions are populated, enabling batch generation of structurally diverse designs — including drainage holes for resin removal after vat photopolymerization (SLA/DLP) printing.

## Overview

For each fingerprint listed in the input file, the script:

1. Constructs the base wall geometries (square, diagonal, and triangular panels), each with a drainage hole to allow uncured resin to escape after printing.
2. Assembles the walls indicated by the fingerprint's binary digits into a unit cell, merges them into a single part, and mirrors it to complete the cell.
3. Adds a simple steel load-cell fixture, material properties, contact interactions, an implicit dynamics compression step, and boundary conditions.
4. Meshes the unit cell and load cell.
5. Computes the unit cell's mass and appends it to a results file.
6. Saves the completed model as its own `.cae` file.

This tooling was developed to support combinatorial design-space exploration of lattice structures, as described in the accompanying publication.

## Requirements

- Abaqus/CAE (developed and tested against a recent Abaqus release; requires the `main`/`secondary` contact-pair API, so older versions using `master`/`slave` may need adjustment).
- No external Python packages beyond what ships with Abaqus/CAE's scripting kernel.

## Repository Structure

```
.
├── Model-creator-module-loadcellless_0.25.py   # Main generation + simulation script
├── fingerprint.txt                             # Example input: one fingerprint per line
├── output/                                      # Created automatically; holds results
│   ├── mass_results.txt                        # fingerprint, mass_kg for each cell
│   └── <fingerprint>.cae                       # One CAE file per generated structure
└── README.md
```

## Usage

1. Create a `fingerprint.txt` file in the working directory with one 20-digit binary string per line, e.g.:

   ```
   11000000000000000000
   10100000000000000000
   01100000000000000000
   ```

   Each digit corresponds to one of the 20 possible wall positions in the unit cell (see comments in the script for the position-to-index mapping).

2. Run the script through Abaqus/CAE:

   ```bash
   abaqus cae noGUI=Model-creator-module-loadcellless_0.25.py
   ```

   Or open it from within the Abaqus/CAE Python kernel (File → Run Script).

3. Results appear in `output/`:
   - `mass_results.txt` — a CSV-style log of `fingerprint, mass_kg` for every successfully built structure.
   - One `.cae` file per fingerprint, containing the full model (geometry, mesh, materials, interactions, and step definition) ready for job submission.

Structures that fail to build (e.g. due to an invalid or geometrically infeasible fingerprint) are skipped, logged to the console with the error message, and do not appear in the results file.

## Configuration

Key parameters can be adjusted near the top of the script:

| Variable | Description | Default |
|---|---|---|
| `FINGERPRINT_FILE` | Path to the input file of fingerprints | `fingerprint.txt` |
| `OUTPUT_DIR` | Directory for `.cae` files and the mass log | `./output` |
| `width` | Unit cell width (mm) | `10` |
| `thickness` | Wall thickness (mm) | `0.25` |
| `hole_radius` | Drainage hole radius (mm) | `0.7` |

Material properties (elastic modulus, Poisson's ratio, density) for both the lattice material and the steel load-cell fixture are defined inline in the script and can be edited directly.

## Units

Abaqus is unit-agnostic; this script follows the mm–tonne–N–MPa–s convention (densities are given in tonne/mm³, e.g. `2.85e-09`).

## Citing

If you use this script in academic work, please cite the associated publication:

> https://doi.org/10.1016/j.giant.2024.100282

