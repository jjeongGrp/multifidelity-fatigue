# Multifidelity Fatigue

A multi-fidelity computational framework for fatigue assessment of laser powder bed fusion (LPBF) 316L stainless steel.

This repository accompanies the manuscript:

> **An integrative multi-fidelity framework for fatigue life prediction in additively manufactured 316L stainless steel**

The framework combines X-ray computed tomography (XCT), electron backscatter diffraction (EBSD), melt-pool-informed synthetic microstructure generation, full-field crystal-elasticity hotspot screening, local crystal-plasticity finite-element method (CPFEM) simulations, and Coffin--Manson--Basquin (CMB) fatigue-life mapping.

The repository includes case files for three LPBF process conditions:

- optimized (`OPT`);
- lack of fusion (`LOF`);
- keyholing (`KH`).

---

## Workflow

```text
Processed XCT pore mask and EBSD-derived statistics
                         │
                         ▼
       MicroGen: synthetic microstructure generation
                         │
                         ▼
   ELAS3D-XTAL: full-field crystal-elasticity screening
                         │
                         ▼
       Pore-aware hotspot identification and ranking
                         │
                         ▼
  Hotspot-domain extraction and local CPFEM simulation
                         │
                         ▼
        CMB-based fatigue-life post-processing
                         │
                         ▼
                Fatigue-life prediction
```

XCT provides the three-dimensional pore morphology and spatial distribution. EBSD supplies section-level grain-size, aspect-ratio, and crystallographic-texture statistics used to generate the surrounding synthetic grain structures.

ELAS3D-XTAL screens large voxelized domains to identify candidate fatigue-critical pore or grain neighborhoods. Local PRISMS-Plasticity simulations then resolve the stabilized cyclic response within selected hotspot domains. The resulting strain fields are mapped to fatigue life using the experimentally calibrated CMB relation.

---

## Repository Structure

```text
multifidelity-fatigue/
├── applications/
│   ├── CPFEM-Fatigue/                # Local PRISMS-Plasticity cases
│   │   ├── OPT/                      # OPT case and PRISMS input files
│   │   ├── LOF/                      # LOF case and PRISMS input files
│   │   └── KH/                       # KH case and PRISMS input files
│   │
│   ├── ELAS3D-XTAL/                  # Elastic screening and hotspot extraction
│   │   ├── OPT/
│   │   ├── LOF/
│   │   └── KH/
│   │
│   ├── MicroGen/                     # Synthetic-microstructure case files
│   │   ├── OPT/
│   │   ├── LOF/
│   │   └── KH/
│   │
│   ├── Sensitivity_Microstructure/   # Texture and morphology sensitivity
│   │   ├── KH_101_Texture/
│   │   ├── KH_Grain/
│   │   └── KH_Random_Texture/
│   │
│   └── Sensitivity_Pore_Morphology/  # LOF pore-shape sensitivity cases
│       ├── LOF_XCT_Set61/
│       ├── LOF_XCT_Set62/
│       ├── LOF_XCT_Set63/
│       ├── LOF_XCT_Set64/
│       ├── LOF_Sphere_Set61/
│       ├── LOF_Sphere_Set62/
│       ├── LOF_Sphere_Set63/
│       ├── LOF_Sphere_Set64/
│       ├── LOF_Ellip_Set61/
│       ├── LOF_Ellip_Set62/
│       ├── LOF_Ellip_Set63/
│       └── LOF_Ellip_Set64/
│
├── src/                              # Reusable workflow source files
│   ├── MicroGen_elas3dxtal_input.m
│   ├── elas3dxtal_pcg.f90
│   ├── elas3dxtal_postprocessing.m
│   └── fatigue_life_prediction.m
│
├── LICENSE
└── README.md
```
---

## Main Components

| Component | Function |
|---|---|
| `applications/MicroGen/` | Creates XCT-informed synthetic grain structures using EBSD and melt-pool information |
| `applications/ELAS3D-XTAL/` | Performs voxel-based crystal-elasticity simulations and ranks fatigue hotspots |
| `applications/CPFEM-Fatigue/` | Runs local cyclic CPFEM analyses and maps stabilized strain fields to fatigue life |
| `applications/Sensitivity_Pore_Morphology/` | Compares XCT-resolved, spherical, and ellipsoidal LOF pore geometries |
| `applications/Sensitivity_Microstructure/` | Evaluates KH texture and grain-morphology effects |
| `src/` | Contains reusable implementations of the principal workflow stages |

---
## Main Capabilities

- XCT-resolved pore geometry and spatial distributions
- EBSD-informed grain size, grain aspect ratio, and crystallographic texture
- Melt-pool-informed synthetic microstructure generation
- Voxel-based anisotropic crystal-elasticity simulation using ELAS3D-XTAL
- Pore-aware hotspot detection and morphology-sensitive ranking
- Local cyclic CPFEM simulation using PRISMS-Plasticity;
- CPFEM-informed CMB strain–life prediction
- Pore-morphology and microstructure-sensitivity analyses

---

## 🛠️ Software Requirements

To compile and run the ELAS3D-Xtal workflow on Windows, please install the following software in the order listed below.

The workflow uses:

- **MATLAB** for XCT pore-domain processing, statistical microstructure generation, ELAS3D-XTAL input generation, full-field post-processing, hotspot identification, and fatigue-life prediction;
- **Intel Fortran + HDF5** for compiling and running the ELAS3D-XTAL crystal-elasticity solver;
- **PRISMS-Plasticity** for local crystal-plasticity finite-element simulations;
- **ParaView** for visualizing PRISMS-Plasticity output and exporting equivalent-strain fields used in fatigue-life post-processing.

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
* **Required For:** Running the MATLAB workflow scripts for XCT pore-mask processing, synthetic microstructure generation, ELAS3D-Xtal input generation, full-field stress postprocessing, hotspot identification and extraction, crystallographic-texture processing, and CMB fatigue-life prediction.
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

The standard workflow for an included process condition is:

1. Select the desired condition:

   ```text
   OPT
   LOF
   KH
   ```

2. Review or run the corresponding statistical microstructure case in:

   ```text
   applications/MicroGen/<condition>/
   ```

3. Prepare and run the corresponding ELAS3D-XTAL case in:

   ```text
   applications/ELAS3D-XTAL/<condition>/
   ```

4. Post-process the crystal-elasticity result to identify and rank candidate fatigue hotspots.

5. Use the included PRISMS-Plasticity input files in:

   ```text
   applications/CPFEM-Fatigue/<condition>/
   ```

   to run the selected local CPFEM case.

6. Visualize the PRISMS-Plasticity output and export the tensile- and compressive-extreme equivalent-strain fields using ParaView.

7. Run the MATLAB fatigue-life post-processing workflow using:

   ```text
   src/fatigue_life_prediction.m
   ```

---

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

## ELAS3D-XTAL Execution and Post-Processing

A typical ELAS3D-XTAL workflow is:

1. Generate `input_structure_poly.h5` using the appropriate MicroGen or case-specific MATLAB workflow.
2. Copy or compile the ELAS3D-XTAL executable in the desired case directory.
3. Run:

   ```cmd
   elas3dxtal_pcg.exe
   ```

4. Post-process `fullfield_poly.h5` using:

   ```text
   src/elas3dxtal_postprocessing.m
   ```

The post-processing stage performs pore-aware stress averaging, hotspot identification, ranking, and local-domain extraction according to the case configuration.

---

### CPFEM execution

PRISMS-Plasticity runs use the case-specific `prm.prm`, mesh, orientation, and boundary-condition files. 

Navigate to the desired case directory:

```bash
cd applications/CPFEM-Fatigue/<condition>/
```

A typical PRISMS-Plasticity execution command is:

```bash
mpirun -n <number_of_processes> \
$PLAS_DIR/applications/crystalPlasticity/main prm.prm
```
where:

- `<number_of_processes>` is the desired MPI process count;
- `$PLAS_DIR` is the local PRISMS-Plasticity installation directory;
- `prm.prm` is the principal case parameter file.

The exact command and required filenames may depend on the PRISMS-Plasticity version and the organization of the included case. Consult:

- the files in the selected case directory;
- the [PRISMS-Plasticity documentation](https://github.com/prisms-center/plasticity).

---

## Fatigue-Life Post-Processing

After the local CPFEM simulation:

1. Open the PRISMS-Plasticity output in ParaView.
2. Extract or export the equivalent-strain fields at the tensile and compressive extrema of the stabilized loading cycle.
3. Provide the exported data to:

   ```text
   src/fatigue_life_prediction.m
   ```

4. Verify the case-specific:
   - nominal strain amplitude;
   - ROI definition;
   - CMB parameters;
   - output paths;
   - percentile- or damage-equivalent aggregation settings.

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
