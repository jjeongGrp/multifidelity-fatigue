# Multi-fidelity fatigue

**A multi-fidelity simulation workflow for defect-driven fatigue life prediction in additively manufactured 316L stainless steel.**

This repository accompanies the study:

> **An integrative multi-fidelity framework for fatigue life prediction in additively manufactured 316L stainless steel**

The code provides an integrated computational workflow for XCT-informed microstructure generation, large-scale crystal-elasticity simulation, fatigue-hotspot identification, PRISMS-Plasticity crystal-plasticity subdomain preparation, and CPFE-informed fatigue-life prediction for laser powder bed fused (LPBF) 316L stainless steel.

The framework is designed for process-aware fatigue assessment across optimized EOS, lack-of-fusion, LOF, and keyholing, KH, process windows.



## Overview

This work develops and validates a multi-fidelity framework to predict defect-driven fatigue in LPBF 316L stainless steel. The workflow combines experimental characterization and simulation across multiple fidelity levels:

- **Experimental characterization**
  - EBSD-based grain morphology and texture analysis
  - high-resolution XCT-based pore and defect characterization
  - strain-controlled fatigue testing
  - SEM fractography for crack-initiation validation

- **Large-scale crystal-elasticity screening**
  - voxel-based polycrystalline microstructures
  - XCT-resolved pore distributions
  - pore voxels treated as inactive material
  - full-field elastic stress simulation using an OpenMP/HDF5 Fortran solver
  - fatigue-hotspot identification using Theory-of-Critical-Distances-inspired stress averaging, pore-surface bands, and morphology-aware severity metrics

- **Defect-resolved crystal-plasticity simulation**
  - CPFEM representative volume elements, RVEs, constructed around top-ranked elastic hotspots
  - XCT-resolved pore geometry embedded in synthetic grain structures
  - EBSD-informed upper-tail grain-size and aspect-ratio statistics
  - process-specific crystallographic textures

- **Fatigue-life prediction**
  - stabilized hotspot strain amplitudes extracted from CPFEM
  - Coffin-Manson-Basquin, CMB, strain-life relation calibrated to EOS fatigue data
  - life prediction for EOS, LOF, and KH specimens
  - sensitivity studies for pore morphology and microstructure/texture fidelity

The repository supports the following simulation chain:

```text
XCT pore data
      │
      ▼
Synthetic microstructure generation
      │
      ▼
Large-scale ELAS3D-Xtal elastic simulation
      │
      ▼
Pore-surface hotspot identification
      │
      ▼
CPFEM subdomain generation for PRISMS-Plasticity
      │
      ▼
Hotspot strain extraction
      │
      ▼
Coffin-Manson-Basquin fatigue-life prediction
```

---

## Repository Structure

```text
elas3d-xtal/
├── src/                              # Reusable source files for microstructure generation, ELAS3D-Xtal simulation, postprocessing, PRISMS export, and fatigue prediction
│   ├── MicroGen_elas3dxtal_input.m   # Generates XCT-informed polycrystalline microstructures and writes input_structure_poly.h5
│   ├── elas3dxtal_pcg.f90            # OpenMP/HDF5 preconditioned conjugate-gradient crystal-elasticity solver
│   ├── elas3dxtal_postprocessing.m   # Postprocesses full-field stress results and identifies pore-associated fatigue hotspots
│   ├── prismsplasticity_input.m      # Exports selected hotspot subdomains to PRISMS-Plasticity input files
│   └── fatigue_life_prediction.m     # Predicts fatigue life from CPFE hotspot equivalent strain data
│
├── applications/                     # Case-specific files used to reproduce the paper results
│   ├── sec4_1/                       # Section 4.1: Validation of synthetic morphology, XCT pores, and texture
│   │   ├── Set11/                    # Set11 synthetic microstructure generation and validation case
│   │   ├── Set61/                    # Set61 synthetic microstructure generation and validation case
│   │   └── Set91/                    # Set91 synthetic microstructure generation and validation case
│   │
│   ├── sec4_2/                       # Section 4.2: Large-scale ELAS3D-Xtal elastic simulation and hotspot screening
│   │   ├── Set11/                    # Set11 ELAS3D-Xtal input generation, elastic simulation, hotspot analysis, and PRISMS export
│   │   ├── Set61/                    # Set61 ELAS3D-Xtal input generation, elastic simulation, hotspot analysis, and PRISMS export
│   │   └── Set91/                    # Set91 ELAS3D-Xtal input generation, elastic simulation, hotspot analysis, and PRISMS export
│   │
│   ├── sec4_3/                       # Section 4.3: Crystal-plasticity fatigue simulations for selected hotspot subdomains
│   │   ├── Set11/                    # Set11 PRISMS-Plasticity inputs and CPFE-informed fatigue-life prediction
│   │   ├── Set61/                    # Set61 PRISMS-Plasticity inputs and CPFE-informed fatigue-life prediction
│   │   └── Set91/                    # Set91 PRISMS-Plasticity inputs and CPFE-informed fatigue-life prediction
│   │
│   ├── sec4_4/                       # Section 4.4: Sensitivity analysis for pore morphology
│   │   ├── xct_set61/                # Set61 CPFE case using XCT-resolved pore morphology
│   │   ├── xct_set62/                # Set62 CPFE case using XCT-resolved pore morphology
│   │   ├── xct_set63/                # Set63 CPFE case using XCT-resolved pore morphology
│   │   ├── xct_set64/                # Set64 CPFE case using XCT-resolved pore morphology
│   │   ├── sphere_set61/             # Set61 CPFE case using volume-equivalent spherical pore morphology
│   │   ├── sphere_set62/             # Set62 CPFE case using volume-equivalent spherical pore morphology
│   │   ├── sphere_set63/             # Set63 CPFE case using volume-equivalent spherical pore morphology
│   │   ├── sphere_set64/             # Set64 CPFE case using volume-equivalent spherical pore morphology
│   │   ├── ellip_set61/              # Set61 CPFE case using volume-equivalent ellipsoidal pore morphology
│   │   ├── ellip_set62/              # Set62 CPFE case using volume-equivalent ellipsoidal pore morphology
│   │   ├── ellip_set63/              # Set63 CPFE case using volume-equivalent ellipsoidal pore morphology
│   │   └── ellip_set64/              # Set64 CPFE case using volume-equivalent ellipsoidal pore morphology
│   │
│   └── sec4_5/                       # Section 4.5: Sensitivity analysis for microstructure and texture
│       ├── xct_set92_textured/       # Set92 CPFE case with XCT pore morphology and process-specific textured grains
│       ├── xct_set92_random/         # Set92 CPFE case with XCT pore morphology and random grain orientations
│       └── xct_set92_grain/          # Set92 CPFE case with XCT pore morphology and modified grain morphology/statistics
│
├── exp/                              # XCT pore-domain extraction scripts and process-specific experimental pore data
│   ├── EOS/                          # Optimized EOS condition: XCT pore subdomain extraction for gas/equiaxed pore cases
│   ├── LOF/                          # Lack-of-fusion condition: XCT pore subdomain extraction for flat/irregular LOF defects
│   └── KH/                           # Keyhole condition: XCT pore subdomain extraction for keyhole pore cases
│
├── LICENSE                           # License and attribution information
└── README.md                         # Repository documentation

```

---

## 🛠️ Software Requirements

To compile and run the ELAS3D-Xtal workflow on Windows, please install the following software in the order listed below.

The workflow uses:

- **MATLAB** for XCT pore-domain processing, microstructure generation, postprocessing, PRISMS-Plasticity input generation, and fatigue-life prediction.
- **Intel Fortran + HDF5** for compiling and running the ELAS3D-Xtal crystal-elasticity solver.
- **PRISMS-Plasticity** for crystal-plasticity finite element, CPFE, simulations.
- **ParaView** for extracting CPFE equivalent strain data from simulation output files.

---

### Visual Studio 2022 Community, Professional, or Enterprise

* **Required For:** The Microsoft linker, `link.exe`, and build tools required by the Intel Fortran Compiler on Windows.
* **Download:** [Visual Studio 2022](https://visualstudio.microsoft.com/downloads/)
* **Installation Note:** During installation, select the following **Workloads**:
    * **Desktop development with C++**
    * **.NET desktop development**
    * **Python development**, optional but useful for scientific workflows

---

### Intel Fortran Compiler, `ifx`

* **Required For:** Compiling the Fortran source code, `elas3dxtal_pcg.f90`, with OpenMP support.
* **Download:** [Intel oneAPI HPC Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html)
* **Installation Notes:**
    1. Run the installer and select **Custom Installation**.
    2. When prompted, **"Would you like to integrate with an IDE?"**, ensure that the checkbox for **Microsoft Visual Studio 2022** is selected.
    3. After installation, compile the solver from the **Intel oneAPI Command Prompt for Intel 64**, not from the standard Windows `cmd.exe`.

---

### HDF5 Library, Version 2.0.0

* **Required For:** High-performance binary input/output for large voxel grids and full-field simulation results.
* **Download:** [HDF5 2.0.0 Release on GitHub](https://github.com/HDFGroup/hdf5/releases/tag/hdf5_2.0.0)
* **Recommended Windows Installer:** `hdf5-2.0.0-win-vs2022_intel.msi`
* **SHA256 Hash:** `8c625b68cb9b429208391f33b4ae675513c65c226ad3fe2a10226d93b16e2d35`
* **Default Installation Path Assumed by the Build Commands:**

    ```text
    C:\Program Files\HDF_Group\HDF5\2.0.0\
    ```

* **Installation Note:** If HDF5 is installed in a different location, update the include and library paths in the Fortran compilation command.

---

### MATLAB
* **Required For:** Running the MATLAB workflow scripts for XCT pore-domain processing, microstructure generation, ELAS3D-Xtal input generation, stress postprocessing, PRISMS-Plasticity input generation, and fatigue-life prediction.
* **Recommended Version:** R2021b or newer.

**A. Required MathWorks Toolboxes:**
You must have the following toolboxes installed (check via `ver` command in MATLAB):
* **Statistics and Machine Learning Toolbox** (required for `lognrnd` and grain statistics).
* **Image Processing Toolbox** (required for voxel/pore-mask manipulation, connected-component analysis, morphological operations, filtering, and region-property calculations via 'bwconncomp', 'regionprops3', 'imdilate', 'bwareaopen', 'imgaussfilt3', 'smooth3', 'isosurface').
* **Parallel Computing Toolbox** (required for spatial filtering acceleration and fast block processing via `parfor`).

**B. External Libraries (MTEX):**
The code relies on **MTEX** for crystallographic texture analysis, orientation generation, Rodrigues-vector handling, IPF coloring, and pole-figure visualization.
1.  **Download:** [MTEX Toolbox Website](https://mtex-toolbox.github.io/)
2.  **Install:**
    * Extract the folder (e.g., to `C:\Matlab_Toolboxes\mtex`).
    * Open MATLAB and run:
        ```matlab
        cd 'C:\Matlab_Toolboxes\mtex'
        startup_mtex
        ```
    * Run the ELAS3D-Xtal MATLAB scripts only *after* MTEX is initialized.

---

### PRISMS-Plasticity
* **Required For:** Crystal-plasticity finite element, CPFE, simulations of selected hotspot subdomains.
* **Download:** [PRISMS-Plasticity 1.5.0 Release on GitHub](https://github.com/prisms-center/plasticity)

---

### ParaView
* **Required For:** Extracting equivalent strain data from PRISMS-Plasticity CPFE simulation output files for fatigue-life prediction.
* **Download:** [ParaView Homepage](https://www.paraview.org/download/)
* **Recommended Version:** ParaView 6.1.1.

---

## ⚙️ Compilation & Usage

This section describes how to compile the ELAS3D-Xtal Fortran solver and run the main workflow.

The typical workflow is:

```text
MATLAB microstructure generation
        │
        ▼
input_structure_poly.h5
        │
        ▼
ELAS3D-Xtal Fortran solver
        │
        ▼
fullfield_poly.h5
        │
        ▼
MATLAB stress postprocessing and hotspot identification
        │
        ▼
PRISMS-Plasticity input generation
        │
        ▼
CPFEM simulation
        │
        ▼
ParaView strain extraction
        │
        ▼
MATLAB fatigue-life prediction
```

### 1. Compile the ELAS3D-Xtal Solver

The main Fortran solver is:

```text
elas3dxtal_pcg.f90
```

This file is located in:

```text
src/
```

and also appears as case-specific copies inside some `applications/sec4_2/` folders.

---

#### Windows Compilation with Intel Fortran and HDF5

Open the:

```text
Intel oneAPI Command Prompt for Intel 64
```

Do **not** use the standard Windows `cmd.exe`, because the standard command prompt may not have the Intel compiler and linker environment variables loaded.

Navigate to the folder containing `elas3dxtal_pcg.f90`. For example:

```cmd
cd path\to\elas3d-xtal\src
```

Compile using:

```cmd
ifx /O3 /QxHost /Qunroll /Qopenmp /Qipo /I"C:\Program Files\HDF_Group\HDF5\2.0.0\mod\shared" elas3dxtal_pcg.f90 /Fe:elas3dxtal_pcg.exe /link /LIBPATH:"C:\Program Files\HDF_Group\HDF5\2.0.0\lib" hdf5_fortran.lib hdf5.lib
```

This creates:

```text
elas3dxtal_pcg.exe
```

---
#### Notes on HDF5 Paths

The compile command above assumes that HDF5 is installed at:

```text
C:\Program Files\HDF_Group\HDF5\2.0.0\
```

If HDF5 is installed elsewhere, update these paths:

```text
/I"C:\Program Files\HDF_Group\HDF5\2.0.0\mod\shared"
```

and:

```text
/LIBPATH:"C:\Program Files\HDF_Group\HDF5\2.0.0\lib"
```

to match your local HDF5 installation.

---

#### Case-Specific Compilation

For application cases in:

```text
applications/sec4_2/
```

Each set folder may contain a case-specific copy of:

```text
elas3dxtal_pcg.f90
```

For example:

```text
applications/sec4_2/Set11/elas3dxtal_pcg.f90
applications/sec4_2/Set61/elas3dxtal_pcg.f90
applications/sec4_2/Set91/elas3dxtal_pcg.f90
```

To reproduce a specific case, navigate to the corresponding set folder and compile the local solver file there.

Example:

```cmd
cd path\to\elas3d-xtal\applications\sec4_2\Set11
```

Then run the same compile command:

```cmd
ifx /O3 /QxHost /Qunroll /Qopenmp /Qipo ^
/I"C:\Program Files\HDF_Group\HDF5\2.0.0\mod\shared" ^
elas3dxtal_pcg.f90 ^
/Fe:elas3dxtal_pcg.exe ^
/link /LIBPATH:"C:\Program Files\HDF_Group\HDF5\2.0.0\lib" ^
hdf5_fortran.lib hdf5.lib
```

---


### 2. Generate the ELAS3D-Xtal Input File

The ELAS3D-Xtal solver requires the HDF5 input file:

```text
input_structure_poly.h5
```

This file is generated using MATLAB.

For the generic source workflow, run:

```matlab
elas3dxtal_input
```

For case-specific workflows in `applications/sec4_2/`, run the corresponding set-specific script.

Examples:

```matlab
elas3dxtal_input_set11
```

```matlab
elas3dxtal_input_set61
```

```matlab
elas3dxtal_input_set91
```

The input-generation script creates files such as:

```text
input_structure_poly.h5
xct_poly_params_SS316L.mat
xct_poly_params_SS316L_full.mat
voxel_volume_summary.txt
```

The file required by the Fortran solver is:

```text
input_structure_poly.h5
```

---
### 3. Run the ELAS3D-Xtal Solver

After `input_structure_poly.h5` has been generated, run the compiled solver from the same folder:

```cmd
elas3dxtal_pcg.exe
```

The solver reads:

```text
input_structure_poly.h5
```

and writes:

```text
fullfield_poly.h5
```

The output file `fullfield_poly.h5` contains the full-field elastic simulation results, including stress fields used for hotspot analysis.

---
### 4. Postprocess ELAS3D-Xtal Results

After the solver finishes, run the MATLAB postprocessing script.

For the generic source workflow:

```matlab
elas3dxtal_postprocessing
```

For case-specific workflows, run the corresponding set-specific script.

Examples:

```matlab
elas3dxtal_postprocessing_set11
```

```matlab
elas3dxtal_postprocessing_set61
```

```matlab
elas3dxtal_postprocessing_set91
```

The postprocessing script reads:

```text
fullfield_poly.h5
input_structure_poly.h5
xct_poly_params_SS316L.mat
```

and generates outputs such as:

```text
hotspot_analysis_log.txt
xct_poly_post_SS316L_full.mat
```

The script identifies high-stress hotspot regions near pores and provides recommended local domains for subsequent crystal-plasticity finite element simulations.

---
### 5. Generate PRISMS-Plasticity Input Files

After hotspot identification, generate a local PRISMS-Plasticity subdomain around the selected hotspot.

For the generic source workflow:

```matlab
prismsplasticity_input
```

For case-specific workflows, run the corresponding set-specific script.

Examples:

```matlab
prismsplasticity_input_set11
```

```matlab
prismsplasticity_input_set61
```

```matlab
prismsplasticity_input_set91
```

This generates:

```text
micro.msh
ori.txt
grain_id_map_subdomain.txt
```

The main files required by PRISMS-Plasticity are:

```text
micro.msh
ori.txt
```

For the Section 4.3 CPFE simulations, copy these files from the corresponding `sec4_2` folder to the matching `sec4_3` folder.

Example for Set11:

```text
applications/sec4_2/Set11/micro.msh  →  applications/sec4_3/Set11/micro.msh
applications/sec4_2/Set11/ori.txt    →  applications/sec4_3/Set11/ori.txt
```

---
### 6. Run PRISMS-Plasticity CPFEM Simulations

Each PRISMS-Plasticity simulation folder contains input files such as:

```text
BCinfoTable.txt
LatentHardeningRatio.txt
prm.prm
slipDirections.txt
slipNormals.txt
micro.msh
ori.txt
```

A typical command for running PRISMS-Plasticity is:

```bash
mpirun -n <number_of_cpus> $PLAS_DIR/applications/crystalPlasticity/main prm.prm
```

where:

```text
<number_of_cpus>
```

is the number of MPI processes, and:

```text
$PLAS_DIR
```

is the PRISMS-Plasticity installation directory.

Please refer to the official PRISMS-Plasticity documentation for installation, compilation, and detailed solver usage:

```text
https://github.com/prisms-center/plasticity
```

---
### 7. Extract CPFE Equivalent Strain Data Using ParaView

After the PRISMS-Plasticity simulation finishes, open the output file in ParaView.

Typical output file:

```text
projectedFields-0000.pvtu
```

Equivalent strain data should be extracted at:

```text
+0.003 strain  → maximum tensile strain state
-0.003 strain  → maximum compressive strain state
```

For the loading/output schedule used in this repository, the corresponding time frames are typically:

```text
time frame 30  → +0.003 strain
time frame 90  → -0.003 strain
```

Export CSV files using ParaView:

```text
Filters → Domain → Data Analysis → Probe Location → Apply → Export Spreadsheet
```

Example output filenames:

```text
data_export_emax_plus_full_set1_1.csv
data_export_emax_minus_full_set1_1.csv
```

Revise the filenames for each set or sensitivity-analysis case.

The exported CSV files should include columns such as:

```text
Eqv_strain
Points_0
Points_1
Points_2
```

---
### 8. Run Fatigue-Life Prediction

After exporting the equivalent strain CSV files, run the corresponding MATLAB fatigue-life prediction script.

For the generic source workflow:

```matlab
fatigue_life_prediction
```

For Section 4.3 case-specific workflows:

```matlab
fatigue_life_set11
```

```matlab
fatigue_life_set61
```

```matlab
fatigue_life_set91
```

For Section 4.4 pore-morphology sensitivity cases, use the corresponding script, for example:

```matlab
fatiguelife_xct_set61
```

For Section 4.5 microstructure-sensitivity cases, use the corresponding script, for example:

```matlab
fatiguelife_set92_textured
```

The fatigue-life scripts read the positive and negative CPFE strain CSV files, compute a hotspot equivalent strain amplitude, and generate a CPFE-informed fatigue-life prediction curve.

---

## 📊 Reproducing Paper Results

The `applications/` folder contains case-specific scripts and input files used to reproduce the main results and sensitivity studies.

```text
applications/
├── sec4_1/   # Section 4.1: Validation of synthetic morphology, pores, and texture
├── sec4_2/   # Section 4.2: Large-scale elastic simulation results
├── sec4_3/   # Section 4.3: Crystal plasticity fatigue simulation results
├── sec4_4/   # Section 4.4: Sensitivity analysis: pore morphology
└── sec4_5/   # Section 4.5: Sensitivity analysis: microstructure
```

---
	
### 1. Validation of Synthetic Morphology, Pores, and Texture (Section 4.1)
Location:

```text
applications/sec4_1/
```

Subfolders:

```text
Set11/
Set61/
Set91/
```

This section validates the synthetic microstructure generation procedure by importing XCT-derived pore distributions and generating voxelized synthetic polycrystalline microstructures with process-specific morphology and texture.

Each folder contains case-specific MATLAB files to generate and plot the synthetic microstructure. The primary script follows the naming pattern:

```text
microgen_elas3dxtal_setXX.m
```

For example:

```text
microgen_elas3dxtal_set11.m
microgen_elas3dxtal_set61.m
microgen_elas3dxtal_set91.m
```

### Purpose

The scripts in this section:

- import XCT-derived pore distributions,
- generate 3D voxel-based synthetic polycrystalline microstructures,
- reproduce process-specific grain morphology and texture,
- embed XCT pore distributions into the synthetic microstructure,
- generate validation plots for morphology, pores, and texture.

### How to Run

Open MATLAB, navigate to the desired case folder, initialize MTEX if needed, and run the set-specific script.

Example for Set11:

```matlab
cd applications/sec4_1/Set11
startup_mtex
microgen_elas3dxtal_set11
```

Example for Set61:

```matlab
cd applications/sec4_1/Set61
startup_mtex
microgen_elas3dxtal_set61
```

Example for Set91:

```matlab
cd applications/sec4_1/Set91
startup_mtex
microgen_elas3dxtal_set91
```

### Expected Outputs

Depending on the case-specific script, outputs may include:

```text
microstructure visualization figures
pore-distribution plots
IPF orientation maps
pole figures
MATLAB workspace or parameter files
```

---
		
### 2. Large-Scale Elastic Simulation Results (Section 4.2)
Location:

```text
applications/sec4_2/
```

Subfolders:

```text
Set11/
Set61/
Set91/
```

This section performs large-scale ELAS3D-Xtal elastic simulations using XCT-informed synthetic microstructures. Each set folder contains the files needed to generate the ELAS3D-Xtal input, run the Fortran solver, postprocess the elastic stress field, identify hotspots, and export a local PRISMS-Plasticity subdomain.

Each folder contains files such as:

```text
*.mat
elas3dxtal_input_setXX.m
elas3dxtal_pcg.f90
elas3dxtal_postprocessing_setXX.m
prismsplasticity_input_setXX.m
```

where `XX` corresponds to the set number, for example, `11`, `61`, or `91`.

### Workflow

#### Step 1: Generate ELAS3D-Xtal Input File

Run the set-specific MATLAB input-generation script.

Example for Set11:

```matlab
cd applications/sec4_2/Set11
startup_mtex
elas3dxtal_input_set11
```

This generates:

```text
input_structure_poly.h5
xct_poly_params_SS316L.mat
xct_poly_params_SS316L_full.mat
voxel_volume_summary.txt
```

The required file for the ELAS3D-Xtal solver is:

```text
input_structure_poly.h5
```

---

#### Step 2: Compile and Run ELAS3D-Xtal

Compile the local `elas3dxtal_pcg.f90` file in the same set folder.

Example for Set11 on Windows:

```cmd
cd path\to\elas3d-xtal\applications\sec4_2\Set11
```

```cmd
ifx /O3 /QxHost /Qunroll /Qopenmp /Qipo ^
/I"C:\Program Files\HDF_Group\HDF5\2.0.0\mod\shared" ^
elas3dxtal_pcg.f90 ^
/Fe:elas3dxtal_pcg.exe ^
/link /LIBPATH:"C:\Program Files\HDF_Group\HDF5\2.0.0\lib" ^
hdf5_fortran.lib hdf5.lib
```

Run the solver:

```cmd
elas3dxtal_pcg.exe
```

The solver reads:

```text
input_structure_poly.h5
```

and writes:

```text
fullfield_poly.h5
```

---

#### Step 3: Postprocess Elastic Stress and Identify Hotspots

After the ELAS3D-Xtal simulation finishes, run the set-specific postprocessing script.

Example for Set11:

```matlab
elas3dxtal_postprocessing_set11
```

This script:

- reads `fullfield_poly.h5`,
- reconstructs the padded simulation grid,
- computes von Mises stress maps,
- performs pore-surface hotspot detection,
- ranks hotspots using morphology-aware severity metrics,
- suggests a local CPFE simulation domain.

Typical outputs include:

```text
hotspot_analysis_log.txt
xct_poly_post_SS316L_full.mat
hotspot visualization figures
```

---

#### Step 4: Generate PRISMS-Plasticity Input Files

Run the set-specific PRISMS-Plasticity input-generation script.

Example for Set11:

```matlab
prismsplasticity_input_set11
```

This generates:

```text
micro.msh
ori.txt
grain_id_map_subdomain.txt
```

The two main files needed for PRISMS-Plasticity are:

```text
micro.msh
ori.txt
```

Copy these files to the matching folder in `applications/sec4_3/`.

Example for Set11:

```text
applications/sec4_2/Set11/micro.msh  →  applications/sec4_3/Set11/micro.msh
applications/sec4_2/Set11/ori.txt    →  applications/sec4_3/Set11/ori.txt
```

---
	
	
### 3. Crystal Plasticity Fatigue Simulation Results (Section 4.3)
Location:

```text
applications/sec4_3/
```

Subfolders:

```text
Set11/
Set61/
Set91/
```

This section contains PRISMS-Plasticity simulation inputs and fatigue-life postprocessing scripts for the hotspot subdomains identified in Section 4.2.

Each set folder contains files such as:

```text
BCinfoTable.txt
LatentHardeningRatio.txt
prm.prm
slipDirections.txt
slipNormals.txt
fatigue_life_setXX.m
```

The following files should be copied from the corresponding `applications/sec4_2/SetXX/` folder:

```text
micro.msh
ori.txt
```

### Step 1: Copy PRISMS-Plasticity Mesh and Orientation Files

Example for Set11:

```text
applications/sec4_2/Set11/micro.msh  →  applications/sec4_3/Set11/micro.msh
applications/sec4_2/Set11/ori.txt    →  applications/sec4_3/Set11/ori.txt
```

---

### Step 2: Run PRISMS-Plasticity

Refer to the PRISMS-Plasticity documentation for installation and full usage instructions:

```text
https://github.com/prisms-center/plasticity
```

A typical command is:

```bash
mpirun -n <number_of_cpus> $PLAS_DIR/applications/crystalPlasticity/main prm.prm
```

where:

```text
<number_of_cpus>
```

is the number of MPI processes.

---

### Step 3: Extract Equivalent Strain CSV Files

After the CPFE simulation finishes, use ParaView to extract equivalent strain data from:

```text
projectedFields-0000.pvtu
```

Extract data at:

```text
time frame 20  → +0.003 tensile strain
time frame 60  → -0.003 compressive strain
```

Export CSV files using filenames expected by the set-specific fatigue-life script.

Example for Set11:

```text
data_export_emax_plus_full_set1_1.csv
data_export_emax_minus_full_set1_1.csv
```

The CSV files should contain:

```text
Eqv_strain
Points_0
Points_1
Points_2
```

---

### Step 4: Run Fatigue-Life Prediction

Run the set-specific fatigue-life script.

Example for Set11:

```matlab
fatigue_life_set11
```

For other sets:

```matlab
fatigue_life_set61
fatigue_life_set91
```

The script reads the positive- and negative-equivalent strain CSV files, computes hotspot strain amplitudes, and generates the CPFE-informed fatigue-life prediction curve.

---

## 4.4 Sensitivity Analysis: Pore Morphology

Location:

```text
applications/sec4_4/
```

This section studies how simplified pore morphology affects CPFE-predicted fatigue life relative to XCT-resolved pore geometry.

The folder contains 12 subfolders:

```text
xct_set61/
xct_set62/
xct_set63/
xct_set64/

sphere_set61/
sphere_set62/
sphere_set63/
sphere_set64/

ellip_set61/
ellip_set62/
ellip_set63/
ellip_set64/
```

### Morphology Classes

| Folder prefix | Description |
|---|---|
| `xct_` | XCT-resolved pore morphology |
| `sphere_` | Volume-equivalent spherical pore approximation |
| `ellip_` | Volume-equivalent ellipsoidal pore approximation |

The goal is to compare fatigue-life predictions obtained from:

- the original XCT pore geometry,
- simplified spherical pore geometry,
- simplified ellipsoidal pore geometry.

---

### Files in Each Subfolder

Each subfolder contains PRISMS-Plasticity input files and case-specific MATLAB scripts.

Typical PRISMS-Plasticity files include:

```text
BCinfoTable.txt
LatentHardeningRatio.txt
prm.prm
slipDirections.txt
slipNormals.txt
```

The PRISMS-Plasticity mesh and orientation files are generated using morphology-specific MATLAB scripts.

Examples:

```text
prismsplasticity_input_xctset61.m
prismsplasticity_input_sphset61.m
prismsplasticity_input_ellset61.m
```

These scripts generate:

```text
micro.msh
ori.txt
grain_id_map_subdomain.txt
```

Each subfolder also contains a fatigue-life prediction script, for example:

```text
fatiguelife_xct_set61.m
fatiguelife_sph_set61.m
fatiguelife_ell_set61.m
```

The exact filename should be checked in the corresponding folder.

---

### Workflow

The workflow for each pore-morphology case is the same as the CPFE workflow in Section 4.3.

#### Step 1: Generate PRISMS-Plasticity Input Files

Navigate to the desired folder and run the morphology-specific MATLAB script.

Example for XCT Set61:

```matlab
cd applications/sec4_4/xct_set61
prismsplasticity_input_xctset61
```

Example for spherical Set61:

```matlab
cd applications/sec4_4/sphere_set61
prismsplasticity_input_sphset61
```

Example for ellipsoidal Set61:

```matlab
cd applications/sec4_4/ellip_set61
prismsplasticity_input_ellset61
```

Expected outputs:

```text
micro.msh
ori.txt
grain_id_map_subdomain.txt
```

---

#### Step 2: Run PRISMS-Plasticity

Run the CPFE simulation using PRISMS-Plasticity.

Typical command:

```bash
mpirun -n <number_of_cpus> $PLAS_DIR/applications/crystalPlasticity/main prm.prm
```

Refer to the PRISMS-Plasticity documentation for installation and detailed usage:

```text
https://github.com/prisms-center/plasticity
```

---

#### Step 3: Extract Equivalent Strain CSV Files

After the CPFE simulation finishes, open the PRISMS-Plasticity output in ParaView.

Typical output file:

```text
projectedFields-0000.pvtu
```

Extract equivalent strain data at:

```text
time frame 20  → +0.003 tensile strain
time frame 60  → -0.003 compressive strain
```

Export the data as CSV files. Example filenames for XCT Set61:

```text
data_export_emax_plus_full_xct_set61.csv
data_export_emax_minus_full_xct_set61.csv
```

The exact filenames should match the `csvFiles` variable in the corresponding fatigue-life script.

---

#### Step 4: Run Fatigue-Life Prediction

Run the case-specific fatigue-life script.

Example for XCT Set61:

```matlab
fatiguelife_xct_set61
```

Example for spherical Set61:

```matlab
fatiguelife_sph_set61
```

Example for ellipsoidal Set61:

```matlab
fatiguelife_ell_set61
```

The script reads the positive- and negative-equivalent strain CSV files and generates the fatigue-life prediction curve for that pore-morphology case.

---
## 4.5 Sensitivity Analysis: Microstructure

Location:

```text
applications/sec4_5/
```

This section examines how assumptions about microstructure and texture affect CPFE-predicted fatigue life.

The folder contains three subfolders:

```text
xct_set92_textured/
xct_set92_random/
xct_set92_grain/
```

### Microstructure Cases

| Folder | Description |
|---|---|
| `xct_set92_textured/` | XCT pore morphology with process-specific textured grain orientations |
| `xct_set92_random/` | XCT pore morphology with random grain orientations |
| `xct_set92_grain/` | XCT pore morphology with modified grain morphology/statistics |

Each folder contains PRISMS-Plasticity input files, including:

```text
micro.msh
ori.txt
BCinfoTable.txt
LatentHardeningRatio.txt
prm.prm
slipDirections.txt
slipNormals.txt
```

Each folder also contains a fatigue-life prediction script:

```text
fatiguelife_set92_textured.m
fatiguelife_set92_random.m
fatiguelife_set92_grain.m
```

---

### Workflow

The CPFE workflow is identical to Section 4.3.

#### Step 1: Run PRISMS-Plasticity

Navigate to the desired folder and run:

```bash
mpirun -n <number_of_cpus> $PLAS_DIR/applications/crystalPlasticity/main prm.prm
```

Example:

```bash
cd applications/sec4_5/xct_set92_textured
mpirun -n <number_of_cpus> $PLAS_DIR/applications/crystalPlasticity/main prm.prm
```

---

#### Step 2: Extract Equivalent Strain CSV Files

After the simulation completes, open the output in ParaView:

```text
projectedFields-0000.pvtu
```

Extract equivalent strain data at:

```text
time frame 20  → +0.003 tensile strain
time frame 60  → -0.003 compressive strain
```

Export the CSV files using filenames expected by the corresponding fatigue-life script.

Example:

```text
data_export_emax_plus_full_set92_textured.csv
data_export_emax_minus_full_set92_textured.csv
```

---

#### Step 3: Run Fatigue-Life Prediction

Run the corresponding MATLAB script.

For textured microstructure:

```matlab
fatiguelife_set92_textured
```

For random texture:

```matlab
fatiguelife_set92_random
```

For modified grain morphology/statistics:

```matlab
fatiguelife_set92_grain
```

The scripts generate CPFE-informed fatigue-life prediction curves for the microstructure-sensitivity study.

---

## 🧪 Experimental XCT Pore-Domain Preparation

Location:

```text
exp/
```

Subfolders:

```text
EOS/
LOF/
KH/
```

The `exp/` folder contains scripts for extracting user-defined pore subdomains from XCT scan data for each LPBF process/defect class.

| Folder | Process / Defect Class | Description |
|---|---|---|
| `EOS/` | Optimized EOS condition | Gas/equiaxed pore dominated XCT pore distributions |
| `LOF/` | Lack-of-fusion condition | Large, flat, irregular lack-of-fusion defect distributions |
| `KH/` | Keyhole condition | Keyhole pore distributions |

---

### Purpose

The scripts in `exp/` convert physical XCT pore data into cubic pore-domain MAT files for use in the ELAS3D-Xtal microstructure-generation workflow.

These extracted pore-domain files can be used as input data in:

```text
applications/sec4_2/
```

for large-scale ELAS3D-Xtal elastic simulations.

---

### Typical Scripts

Each process folder contains a process-specific MATLAB script to select and export a user-defined XCT pore subdomain.

Examples:

```text
exp/EOS/phy2cube_EOS_subdomain.m
exp/LOF/phy2cube_LOF_subdomain.m
exp/KH/phy2cube_KH_subdomain.m
```

---

### General Workflow

1. Navigate to the desired process folder.

   Example:

   ```matlab
   cd exp/EOS
   ```

2. Run the corresponding pore-domain extraction script.

   Example for EOS:

   ```matlab
   phy2cube_EOS_subdomain
   ```

   Example for LOF:

   ```matlab
   phy2cube_LOF_subdomain
   ```

   Example for KH:

   ```matlab
   phy2cube_KH_subdomain
   ```

3. The script extracts a user-defined pore subdomain from the XCT scan data and writes a MAT file containing the pore distribution.

4. Use the generated MAT file as input for the corresponding `elas3dxtal_input_*.m` script in:

   ```text
   applications/sec4_2/
   ```

---

### Connection to ELAS3D-Xtal Simulations

The XCT-derived pore-domain MAT files from `exp/` are used to create microstructures for large-scale elastic simulation.

For example:

```text
exp/EOS/
      │
      ▼
applications/sec4_2/Set11/
      │
      ▼
elas3dxtal_input_set11.m
      │
      ▼
input_structure_poly.h5
      │
      ▼
elas3dxtal_pcg.exe
```

---


## 📄 Citation

If you use this code in your research, please cite **both** the new methodology and the original NIST algorithm:

**1. Multi-fidelity fatigue (This Work):**

> Jeong, J., Bidar, A., Jame. M.S.R., Andani, M.T., Shao, S., Shamsaei, N., & Sundararaghavan, V. (2026). *An integrative multi-fidelity framework for fatigue life prediction in additively manufactured 316L stainless steel*. *Under review* (2026)

**2. Original NIST Solver:**

> Garboczi, E. (1998). *Finite Element and Finite Difference Programs for Computing the Linear Electric and Elastic Properties of Digital Images of Random Materials*. NIST Interagency/Internal Report (NISTIR) 6269.

### BibTeX

```bibtex
@article{garboczi1998finite,
  title={Finite element and finite difference programs for computing the linear electric and elastic properties of digital images of random materials},
  author={Garboczi, Edward J},
  year={1998},
  journal={Report NISTIR},
  publisher={Edward J. Garboczi}
}

```

---

## ⚖️ License & Attribution

This project uses a **hybrid license** model to respect the original government work:

1. **New contributions:** Licensed under the **MIT License**.
2. **Original Solver:** Derived from NISTIR 6269. This software is a derivative work of NISTIR 6269, republished courtesy of the National Institute of Standards and Technology.

**Full details:** See the [`LICENSE`] file.

---

## 🤝 Acknowledgements

This work was supported by the Defense Advanced Research Projects Agency (DARPA) SURGE program under Cooperative Agreement No. HR0011-25-2-0009, "Predictive Real-time Intelligence for Metallic Endurance (PRIME)."

*Disclaimer: Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of DARPA or the United States Government.*

---

## 📫 Contact

For questions regarding the code, methodology, or the associated publication, please feel free to reach out:

* **Juyoung Jeong** - jjuyoung@umich.edu
* **Veera Sundararaghavan** - veeras@umich.edu
* **Affiliation:** Department of Aerospace Engineering, University of Michigan, Ann Arbor, MI 48109, USA

**Bug Reports & Feature Requests:** If you encounter any issues while compiling or running the solver, please use the [GitHub Issues](https://github.com/jjeongGrp/multifidelity-fatigue/issues) page to report them.

**Project Link:** [https://github.com/jjeongGrp/multifidelity-fatigue](https://github.com/jjeongGrp/multifidelity-fatigue)
