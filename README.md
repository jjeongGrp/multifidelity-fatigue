# Multifidelity Fatigue

A multifidelity framework for fatigue assessment of laser powder bed fusion (LPBF) 316L stainless steel.

This repository accompanies:

> **An integrative multi-fidelity framework for fatigue life prediction in additively manufactured 316L stainless steel**

The framework combines high-resolution X-ray computed tomography (XCT), electron backscatter diffraction (EBSD), melt-pool-informed synthetic microstructure generation, full-field crystal-elasticity hotspot screening, local crystal-plasticity finite element method (CPFEM) simulations, and Coffin–Manson–Basquin (CMB) fatigue-life mapping.

The repository supports optimized (`OPT`), lack-of-fusion (`LOF`), and keyholing (`KH`) LPBF process conditions.

---

## Workflow

```text
XCT pore data
    │
    ▼
MicroGen: synthetic microstructure generation
    │
    ▼
ELAS3D-XTAL: full-field crystal-elasticity simulation
    │
    ▼
Pore-aware hotspot identification and ranking
    │
    ▼
CPFEM-Fatigue: local CPFEM simulation and strain–life mapping
    │
    ▼
Fatigue-life prediction
```

XCT provides pore morphology and spatial distribution. EBSD informs the grain-size, aspect-ratio, and texture statistics used to generate the surrounding synthetic grain structure. ELAS3D-XTAL screens the large voxelized domain to identify fatigue-critical pore or grain neighborhoods, and CPFEM resolves stabilized local cyclic strain fields within selected hotspot domains.

---

## Repository Structure

```text
multifidelity-fatigue/
├── applications/
│   ├── CPFEM-Fatigue/                # Local CPFEM and fatigue-life post-processing
│   │   ├── OPT/
│   │   ├── LOF/
│   │   └── KH/
│   │
│   ├── ELAS3D-XTAL/                  # Crystal-elasticity screening and hotspot extraction
│   │   ├── OPT/
│   │   ├── LOF/
│   │   └── KH/
│   │
│   ├── MicroGen/                     # XCT-informed synthetic microstructure generation
│   │   ├── OPT/
│   │   ├── LOF/
│   │   └── KH/
│   │
│   ├── Sensitivity_Microstructure/   # Texture and grain-morphology sensitivity
│   │   ├── KH_101_Texture/
│   │   ├── KH_Grain/
│   │   └── KH_Random_Texture/
│   │
│   └── Sensitivity_Pore_Morphology/  # Pore-shape sensitivity for LOF specimens
│       ├── LOF_XCT_Set61/ ... LOF_XCT_Set64/
│       ├── LOF_Sphere_Set61/ ... LOF_Sphere_Set64/
│       └── LOF_Ellip_Set61/ ... LOF_Ellip_Set64/
│
├── exp/                              # XCT pore-domain preparation
│   ├── OPT/
│   ├── LOF/
│   └── KH/
│
├── src/                              # Reusable source files
│   ├── MicroGen_elas3dxtal_input.m
│   ├── elas3dxtal_pcg.f90
│   ├── elas3dxtal_postprocessing.m
│   ├── fatigue_life_prediction.m
│   └── prismsplasticity_input.m
│
├── LICENSE
└── README.md
```
---

## Main Components

| Component | Function |
|---|---|
| `exp/` | Prepares pore subdomains from XCT data |
| `MicroGen/` | Creates XCT-informed synthetic grain structures using EBSD and melt-pool information |
| `ELAS3D-XTAL/` | Performs voxel-based crystal-elasticity simulations and ranks fatigue hotspots |
| `CPFEM-Fatigue/` | Runs local cyclic CPFEM analyses and maps stabilized strain fields to fatigue life |
| `Sensitivity_Pore_Morphology/` | Compares XCT-resolved, spherical, and ellipsoidal LOF pore geometries |
| `Sensitivity_Microstructure/` | Evaluates KH texture and grain-morphology effects |
| `src/` | Contains reusable implementations of the workflow stages |

---
## Main Capabilities

- XCT-resolved pore geometry and spatial distributions
- EBSD-informed grain size, grain aspect ratio, and crystallographic texture
- Melt-pool-informed synthetic microstructure generation
- Voxel-based anisotropic crystal-elasticity simulation using ELAS3D-XTAL
- Pore-aware hotspot detection and morphology-sensitive ranking
- PRISMS-Plasticity input generation for selected hotspot subdomains
- CPFEM-informed CMB strain–life prediction
- Pore-morphology and microstructure-sensitivity analyses

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

## Running a Case

The standard workflow for a process condition is:

1. Extract a pore domain from `exp/<condition>/`.
2. Generate an XCT-informed microstructure in `applications/MicroGen/<condition>/`.
3. Run ELAS3D-XTAL input generation, solver execution, and hotspot post-processing in `applications/ELAS3D-XTAL/<condition>/`.
4. Export the selected hotspot domain to PRISMS-Plasticity.
5. Run the corresponding CPFEM case in `applications/CPFEM-Fatigue/<condition>/`.
6. Export tensile and compressive equivalent-strain fields using ParaView.
7. Run the MATLAB fatigue-life post-processing script.

Here, `<condition>` is `OPT`, `LOF`, or `KH`.

### ELAS3D-XTAL compilation

Compile `elas3dxtal_pcg.f90` from an Intel oneAPI command prompt. An example Windows command is:

```cmd
ifx /O3 /QxHost /Qunroll /Qopenmp /Qipo ^
/I"C:\Program Files\HDF_Group\HDF5\2.0.0\mod\shared" ^
elas3dxtal_pcg.f90 ^
/Fe:elas3dxtal_pcg.exe ^
/link /LIBPATH:"C:\Program Files\HDF_Group\HDF5\2.0.0\lib" ^
hdf5_fortran.lib hdf5.lib
```

Update the HDF5 include and library paths for your installation.

The solver reads:

```text
input_structure_poly.h5
```

and writes:

```text
fullfield_poly.h5
```

### CPFEM execution

PRISMS-Plasticity runs use the case-specific `prm.prm`, mesh, orientation, and boundary-condition files. A typical command is:

```bash
mpirun -n <number_of_processes> \
$PLAS_DIR/applications/crystalPlasticity/main prm.prm
```

See the [PRISMS-Plasticity documentation](https://github.com/prisms-center/plasticity) for installation and solver details.

---

## Sensitivity Studies

### Pore morphology

`applications/Sensitivity_Pore_Morphology/` compares three pore representations for LOF specimens:

| Prefix | Pore representation |
|---|---|
| `LOF_XCT_` | Original XCT-resolved pore morphology |
| `LOF_Sphere_` | Volume-equivalent spherical pore representation |
| `LOF_Ellip_` | Volume- and second-moment-equivalent ellipsoidal representation |

The cases `Set61`–`Set64` correspond to individual LOF specimens.

### Microstructure

`applications/Sensitivity_Microstructure/` compares three KH microstructure descriptions:

| Folder | Description |
|---|---|
| `KH_101_Texture/` | Process-specific KH texture and reference grain morphology |
| `KH_Random_Texture/` | Identical grain morphology with random orientations |
| `KH_Grain/` | Modified grain morphology with KH texture retained |

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
