!------------------------------------------------------------------------------
! 3D Periodic Microstructure Elasticity Solver
!------------------------------------------------------------------------------
!
! This program computes the linear elastic response of a 3D voxel/nodal
! microstructure using a regular-grid finite element formulation with periodic
! boundary conditions. The solver is designed for digital microstructures such
! as XCT-derived porous metals and synthetic polycrystals.
!
! The code supports:
!   - isotropic matrix/void calculations, and
!   - anisotropic polycrystalline calculations using grain-wise orientations.
!
! For the polycrystalline case, each solid grain is assigned a crystallographic
! orientation using a Rodrigues vector. The corresponding cubic single-crystal
! stiffness tensor is rotated into the global/sample frame and assigned to that
! grain phase. The pore/void phase is modeled as a zero-stiffness phase.
!
! The equilibrium displacement field is obtained by minimizing the total elastic
! energy using a Preconditioned Conjugate Gradient (PCG) algorithm. The matrix
! operations are performed in a matrix-free/stencil-based form using local
! finite-element stiffness matrices and a 27-point periodic neighbor table.
! OpenMP is used to accelerate assembly, energy evaluation, PCG iterations, and
! post-processing.
!
! The imposed loading is strain-controlled. User-prescribed macroscopic strain
! components are applied through periodic displacement jumps across the domain
! boundaries. After convergence, the code computes volume-averaged stress/strain
! quantities and full-field stress/von Mises stress values.
!
! Numerical formulation:
!   - Element type      : 8-node trilinear hexahedral element
!   - Grid type         : regular 3D periodic grid
!   - Neighbor stencil  : 27-point periodic stencil
!   - Solver            : Preconditioned Conjugate Gradient minimization
!   - Preconditioner    : 3 x 3 block Jacobi preconditioner
!   - Parallelization   : OpenMP
!   - Input/output      : HDF5 and text files
!
! Primary input:
!   - input_structure_poly.h5
!       HDF5 file containing:
!         * /pix
!             1D phase-label array for the simulation grid.
!             Phase IDs 1:(nphase-1) correspond to solid grains.
!             Phase ID nphase corresponds to the pore/void phase.
!
!         * /orientation
!             Rodrigues orientation vectors for the solid grain phases.
!             Stored as:
!               [r1_x, r1_y, r1_z, r2_x, r2_y, r2_z, ...]
!
! Main user-defined parameters:
!   - Grid dimensions and padding:
!       md, pad, nx, ny, nz, ns
!
!   - Microstructure/phase parameters:
!       n_grains, nphase, PORE_PHASE
!
!   - Single-crystal elastic constants:
!       C11_local, C12_local, C44_local
!
!   - Material-model flag:
!       flag_m = 0  isotropic solid phases
!       flag_m = 1  anisotropic polycrystalline solid phases
!
!   - Applied macroscopic strain:
!       aml_exx, aml_eyy, aml_ezz, aml_exz, aml_eyz, aml_exy
!
!   - PCG solver controls:
!       kmax, ldemb, tol, block_size
!
! Outputs:
!   - output_SS316L_polycrystal.txt
!       Text summary of grid information, phase statistics, applied strains,
!       solver progress, final volume-averaged stress/strain values, and timing.
!
!   - cg_history_polycrystal.txt
!       High-level convergence history written from the main program.
!       Contains iteration count, total energy, and relative residual norm.
!
!   - cgitr_SS316L_polycrystal.txt
!       Detailed PCG iteration history written by dembx_OpenMP.
!       Contains local CG iteration number and relative residual norm.
!
!   - fullfield_poly.h5
!       HDF5 file containing full-field output:
!         * /pix
!             Phase-label field.
!         * /vm
!             von Mises equivalent stress field.
!         * /stress_tensor
!             Full stress tensor field stored as an ns x 6 array using Voigt
!             ordering:
!               [sxx, syy, szz, sxz, syz, sxy]
!
! Major workflow:
!   1. Define grid, phase, material, loading, and solver parameters.
!   2. Allocate the main solver data structure.
!   3. Read phase labels and grain orientations from input_structure_poly.h5.
!   4. Build the periodic 27-neighbor connectivity table.
!   5. Build an active-node mask to exclude nodes connected only to pore/void
!      regions.
!   6. Compute phase volume fractions.
!   7. Assemble phase stiffness matrices and local FE stiffness tensors.
!   8. Initialize the displacement field from the applied macroscopic strain.
!   9. Evaluate initial energy and residual.
!  10. Solve for the equilibrium displacement field using PCG.
!  11. Compute volume-averaged stress and strain.
!  12. Compute and export full-field stress and von Mises stress.
!  13. Deallocate memory and report runtime.
!
! Notes:
!   - The Fortran solver operates on grid points/nodes. If the input was created
!     from a voxel-centered MATLAB phase map, the MATLAB preprocessing step should
!     convert voxel phases to nodal phases before writing /pix.
!
!   - The shear strain components exz, eyz, and exy used internally are tensorial
!     shear strains. Engineering shear strains are twice these values.
!
!   - PORE_PHASE is assigned zero stiffness. Nodes connected only to pore/void
!     elements are marked inactive and constrained to zero displacement/residual
!     during the solve.
!
!   - Stress and strain tensor components use Voigt ordering:
!       1 = xx, 2 = yy, 3 = zz, 4 = xz, 5 = yz, 6 = xy
!
! Based on:
!   Garboczi, E. J. (1998). Finite element and finite difference programs for
!   computing the linear electric and elastic properties of digital images of
!   random materials. NISTIR 6269. National Institute of Standards and Technology.
!
! Authors:
!   Juyoung Jeong, Veera Sundararaghavan
!
! Affiliation:
!   Department of Aerospace Engineering, University of Michigan,
!   Ann Arbor, MI 48109, USA.
!
! Repository:
!   https://github.com/jjeongGrp/multifidelity-fatigue
!
! Acknowledgment:
!   This work was supported by the Defense Advanced Research Projects Agency
!   (DARPA) SURGE program under Cooperative Agreement No. HR0011-25-2-0009,
!   "Predictive Real-time Intelligence for Metallic Endurance (PRIME)".
!
! License:
!   MIT License. See LICENSE file in the repository.
!                              
!------------------------------------------------------------------------------  
      module elas3d_mod
      ! ==============================================================================
      ! MODULE: elas3d_mod
      ! ------------------------------------------------------------------------------
      ! PURPOSE:
      !   Defines global parameters, material constants, solver controls, and the
      !   main derived data type used by the 3D periodic linear-elasticity solver.
      !
      !   The module centralizes:
      !     - grid dimensions and finite-element discretization parameters,
      !     - phase and material definitions,
      !     - prescribed macroscopic strain components,
      !     - conjugate-gradient solver controls,
      !     - dynamically allocated arrays grouped in elas3d_data_type.
      !
      !   The derived type elas3d_data_type stores the large field arrays used
      !   throughout the FEM assembly, energy evaluation, PCG solve, and
      !   post-processing steps. Passing this type between subroutines keeps the
      !   solver interface organized and avoids very long argument lists.
      !
      ! NOTES:
      !   - The simulation uses a regular 3D grid with periodic boundary conditions.
      !   - The phase map is stored as a 1D array using Fortran-style linear indexing.
      !   - Grid points/nodes are used by the Fortran solver. If the input is created
      !     from a voxel/cell-centered MATLAB volume, the nodal grid has one more
      !     point than the voxel grid in each direction.
      !   - PORE_PHASE identifies the zero-stiffness pore/void phase.
      !
      ! MAIN DATA STRUCTURE:
      !   type(elas3d_data_type)
      !
      !   Primary components:
      !     u(ns,ndof)
      !       Nodal displacement field.
      !
      !     gb(ns,ndof)
      !       Energy gradient / residual vector used by the CG solver.
      !
      !     h(ns,ndof)
      !       Conjugate-gradient search direction.
      !
      !     b(ns,ndof)
      !       Global right-hand-side/equivalent force vector associated with
      !       imposed macroscopic strain jumps.
      !
      !     cmod(nphmax,6,6)
      !       6x6 constitutive stiffness matrix for each phase.
      !
      !     dk(nphmax,nnode_fe,ndof,nnode_fe,ndof)
      !       Integrated 8-node hexahedral element stiffness matrices for each phase.
      !
      !     ib(ns,nfaces)
      !       Periodic neighbor table. For each grid site, stores the 27 neighboring
      !       site indices used by the finite-element stencil.
      !
      !     pix(ns)
      !       Phase label at each grid site/element location.
      !
      !     orientation(3*(nphase-1))
      !       Rodrigues orientation vectors for the solid grain phases.
      !
      !     rotatedstiffness(36*(nphase-1))
      !       Flattened 6x6 rotated stiffness tensor for each grain phase.
      !
      !     stress_field(ns,6), vm(ns)
      !       Full-field stress tensor and von Mises stress output arrays.
      !
      !     is_active(ns)
      !       Logical mask identifying nodes connected to at least one solid element.
      !       Inactive nodes belong only to pore/void regions and are excluded from
      !       the mechanical solve.
      ! ==============================================================================
       use, intrinsic :: iso_fortran_env, dp=>real64 
       implicit none
       
       ! ======================================================================
       ! ---- 1. Domain and Grid Resolution -----------------------------------
       ! Grid definition for the periodic finite-element solve.
       !
       ! md:
       !   Original XCT/MATLAB voxel count per direction before adding nodal
       !   conversion and padding.
       !
       ! pad:
       !   Number of padding layers added on each side of the domain.
       !
       ! nx, ny, nz:
       !   Number of grid points/nodes in each direction used by the Fortran
       !   solver. For a voxelized input, the corresponding nodal dimension is
       !   voxel_count + 1, plus padding on both sides.
       !
       ! ns:
       !   Total number of grid sites/nodes in the periodic domain.
       ! ======================================================================
       integer, parameter :: md = 546
       integer, parameter :: pad = 11
       integer, parameter :: nx = md+1+pad*2
       integer, parameter :: ny = md+1+pad*2
       integer, parameter :: nz = md+1+pad*2
       integer, parameter :: ns = nx * ny * nz

       ! ======================================================================
       ! ---- 2. Microstructure and FE Topology -------------------------------
       ! Defines the number of material phases and local finite-element topology.
       !
       ! Phase convention:
       !   phases 1:(nphase-1) = solid grain phases
       !   phase nphase        = pore/void phase
       !
       ! nfaces:
       !   Number of entries in the periodic neighbor stencil. Despite the name,
       !   this is a 27-point 3x3x3 neighbor table, not a count of geometric faces.
       !
       ! nnode_fe:
       !   Number of nodes in each trilinear hexahedral finite element.
       !
       ! ngauss:
       !   Number of integration points per coordinate direction for the
       !   3x3x3 element quadrature rule.
       ! ======================================================================
       integer, parameter :: n_grains = 230566
       integer, parameter :: nphase = n_grains+1
       integer, parameter :: nphmax = nphase    
       integer, parameter :: nfaces = 27
       integer, parameter :: ndof = 3
	   integer, parameter :: nnode_fe = 8
       integer, parameter :: ngauss   = 3      
       
       ! ======================================================================
       ! ---- 3. Material Properties: Single-Crystal Stiffness ----------------
       ! Independent elastic constants for the cubic single-crystal stiffness
       ! tensor of the solid material. Units are Pascals.
       !
       ! Reference:
       !   Ledbetter, H. M. (1981). Predicted single-crystal elastic constants
       !   of stainless-steel 316. British Journal of Non-Destructive Testing,
       !   23(6).
       ! ======================================================================   
       real(dp), parameter :: C11_local = 204.60d9
       real(dp), parameter :: C12_local = 137.70d9
       real(dp), parameter :: C44_local = 126.20d9
	   
       ! ======================================================================
       ! ---- 4. Isotropic Equivalents and Phase Flags ------------------------
       ! Computes isotropic Voigt-type equivalent properties from the cubic
       ! single-crystal stiffness constants. These values are used when the
       ! matrix is treated isotropically and may also serve as baseline material
       ! properties.
       !
       ! E1, nu1:
       !   Young's modulus and Poisson's ratio assigned to the pore/void phase.
       !   Here, both are zero to represent a zero-stiffness defect phase.
       !
       ! flag_m:
       !   0 = isotropic matrix/grain phases using E0 and nu0
       !   1 = polycrystalline matrix with orientation-dependent anisotropic grains
       ! ======================================================================
	   real(dp), parameter :: K0 = (C11_local + 2.0d0*C12_local)/3.0d0
       real(dp), parameter :: G0 = (C11_local - C12_local + 3.0d0*C44_local)/5.0d0	   
       real(dp), parameter :: E0 = 9.0d0*K0*G0/(3.0d0*K0 + G0)
       real(dp), parameter :: nu0 = (3.0d0*K0 - 2.0d0*G0) / (2.0d0*(3.0d0*K0+G0))
	   
      ! Pore/void phase properties.  
       real(dp), parameter :: E1 = 0.0d0
       real(dp), parameter :: nu1 = 0.0d0       
	   
       ! Microstructure material model flag.			   
       integer, parameter :: flag_m = 1   
    
       
       ! ======================================================================
       ! ---- 5. Macroscopic Loading ------------------------------------------
       ! Prescribed average strain components applied to the periodic domain.
       !
       ! Normal strain components:
       !   aml_exx, aml_eyy, aml_ezz
       !
       ! Tensorial shear strain components:
       !   aml_exz, aml_eyz, aml_exy
       !
       ! Note:
       !   Engineering shear strain is twice the tensorial shear strain.
       ! ======================================================================
       real(dp), parameter :: aml     =  3.0e-3
       real(dp), parameter :: aml_exx =  0.00d0
       real(dp), parameter :: aml_eyy =  0.00d0       
       real(dp), parameter :: aml_ezz =  aml      
       real(dp), parameter :: aml_exz =  0.00d0
       real(dp), parameter :: aml_eyz =  0.00d0
       real(dp), parameter :: aml_exy =  0.00d0       
       
       ! ======================================================================
       ! ---- 6. PCG Solver Controls ------------------------------------------
       ! Controls the conjugate-gradient minimization iterations.
       !
       ! kmax:
       !   Number of outer restart/load-step calls to dembx_OpenMP.
       !   For the current linear-elastic formulation, kmax = 1 is preferred
       !   because restarting discards conjugate search directions and can slow
       !   convergence.
       !
       ! ldemb:
       !   Maximum number of CG iterations inside one dembx_OpenMP call.
       !
       ! n_iter:
       !   Maximum stored history length across all outer calls.
       !
       ! block_size:
       !   Block length used in the matrix-vector product for cache-friendly
       !   processing.
       ! ======================================================================
       integer, parameter :: kmax = 1
       integer, parameter :: ldemb = 10000
       integer, parameter :: n_iter = kmax * ldemb
       integer, parameter :: block_size = 4096
       
       ! ======================================================================
       ! ---- 7. Convergence and Phase Constants ------------------------------
       ! tol:
       !   Relative residual tolerance for the PCG solve.
       !
       ! PORE_PHASE:
       !   Phase ID assigned to pores/voids. This is the final phase index.
       ! ======================================================================
       real, parameter :: tol = 1.0d-8
       integer, parameter :: PORE_PHASE = nphase   

       ! ======================================================================
       ! ---- 8. Main Solver Data Structure -----------------------------------
       ! Groups all major dynamically allocated arrays and scalar state variables
       ! used by the FEM assembly, CG solver, and post-processing routines.
       ! ======================================================================
       type :: elas3d_data_type
	   
         ! Nodal displacement, residual/gradient, right-hand side, and
         ! conjugate-gradient search direction.
         real(dp), allocatable  :: u(:,:), gb(:,:), b(:,:)
		 real(dp), allocatable  :: h(:,:)
		 
         ! Constitutive stiffness matrices and integrated local element
         ! stiffness matrices for each phase.		 
         real(dp), allocatable  :: cmod(:,:,:), dk(:,:,:,:,:)
		 
         ! Periodic neighbor table and phase-label map.		 
         integer, allocatable   :: ib(:,:), pix(:)
		 
         ! Applied macroscopic strain components.	 
         real(dp) :: exx, eyy, ezz, exz, eyz, exy

         ! Constant energy offset from imposed periodic strain jumps.		 
         real(dp) :: C

         ! Volume-averaged stress components.
         real(dp) :: strxx, stryy, strzz, strxz, stryz, strxy
		 
         ! Volume-averaged strain components.
         real(dp) :: sxx, syy, szz, sxz, syz, sxy
       
         ! Per-site/full-field stress outputs for post-processing.
         real(dp), allocatable  :: vm(:), stress_field(:,:)
       
         ! Grain orientation data and rotated anisotropic stiffness tensors.
         real(dp), allocatable  :: orientation(:)
         real(dp), allocatable  :: rotatedstiffness(:)
		 
         ! CG convergence history arrays.
         real(dp), allocatable  :: energies(:), ggs(:)
         integer, allocatable  :: cgiter(:)
		 
         ! Active-node mask.
         ! .true.  = node is connected to at least one solid element
         ! .false. = node belongs only to pore/void surroundings and is excluded
         !           from displacement updates and residual calculations.
         logical, allocatable :: is_active(:)		 

       end type elas3d_data_type
       
      end module elas3d_mod
    
!=========================================================================    
      program elas3dxtal
      ! ==============================================================================
      ! PROGRAM: elas3d
      ! ------------------------------------------------------------------------------
      ! PURPOSE:
      !   Main driver for the 3D periodic voxel/nodal finite-element elasticity solver.
      !
      !   This program initializes the simulation, reads the microstructure and grain
      !   orientation data, constructs the periodic neighbor table, assembles local
      !   finite-element stiffness data, solves for the equilibrium displacement field
      !   using a preconditioned conjugate-gradient method, and post-processes both
      !   volume-averaged and full-field stress quantities.
      !
      ! OVERVIEW:
      !   The solver operates on a regular 3D periodic grid. Phase labels and grain
      !   orientations are read from an HDF5 input file by ppixel_hdf5. The material
      !   response may be isotropic or polycrystalline anisotropic depending on flag_m.
      !   A prescribed macroscopic strain is applied through periodic displacement
      !   jumps. The resulting equilibrium displacement field is found by minimizing
      !   the total elastic energy.
      !
      ! MAIN WORKFLOW:
      !   1. Allocate all major solver arrays.
      !   2. Initialize strain, stress, energy, and convergence-history variables.
      !   3. Open output files and start the wall-clock timer.
      !   4. Define isotropic baseline material properties for solid and/or pore phases.
      !   5. Build the 27-point periodic neighbor table.
      !   6. Read phase labels and grain orientations from input_structure_poly.h5.
      !   7. Build an active-node mask for nodes connected to solid material.
      !   8. Compute phase volume fractions.
      !   9. Assign prescribed macroscopic strain components.
      !  10. Assemble constitutive matrices, local FE stiffness tensors, force vector,
      !      and energy offset using femat_OpenMP.
      !  11. Initialize the displacement field from the homogeneous macroscopic strain.
      !  12. Evaluate the initial energy and residual.
      !  13. Run the PCG minimization solver dembx_OpenMP.
      !  14. Compute final volume-averaged stress and strain.
      !  15. Compute and write full-field stress and von Mises stress.
      !  16. Deallocate all dynamically allocated memory.
      !
      ! IMPORTANT DATA:
      !   e3d
      !     Main derived type containing solver arrays, phase labels, displacement
      !     fields, stiffness matrices, stress outputs, and convergence histories.
      !
      !   phasemod(:,1)
      !     Initially stores Young's modulus E for isotropic phases, then is converted
      !     to bulk modulus K.
      !
      !   phasemod(:,2)
      !     Initially stores Poisson's ratio nu, then is converted to shear modulus G.
      !
      !   prob(:)
      !     Phase volume fractions computed from e3d%pix.
      !
      !   e3d%is_active(:)
      !     Logical mask identifying nodes connected to at least one solid element.
      !     Inactive pore-only nodes are constrained to zero displacement/residual.
      !
      ! OUTPUT FILES OPENED/WRITTEN HERE:
      !   - output_SS316L_polycrystal.txt
      !       Main text summary, material information, solver progress, averaged
      !       stress/strain, and timing.
      !
      !   - cg_history_polycrystal.txt
      !       High-level convergence history from the main program.
      !
      ! Additional output:
      !   - cgitr_SS316L_polycrystal.txt
      !       Detailed CG iteration history written inside dembx_OpenMP.
      !
      !   - fullfield_poly.h5
      !       Full-field stress and von Mises stress written inside
      !       stress_fullfield_OpenMP.
      ! ==============================================================================	  
      use elas3d_mod
      use omp_lib   
      implicit none

      ! ======================================================================
      ! ---- 1. Variable Declarations ----------------------------------------
      ! Loop indices, file/status variables, timing variables, material-property
      ! work arrays, convergence measures, and the main solver data structure.
      ! ======================================================================
      integer :: i1, j1, k1, m1
      integer :: i, j, k, m, n
      integer :: micro, npoints, ltot, kkk, Lstep
      integer :: ios, nxy
      integer :: im(nfaces), jm(nfaces), km(nfaces)
      real(dp) :: gg, utot, youngs
      real(dp) :: x, y, z, elapsed_time, gg0, phmod(2,2)
	  real(dp) :: rel_res
      real(dp), allocatable :: phasemod(:,:), prob(:)
      integer :: t1, t2, tc 
      integer :: lhist
	  
      ! Main solver state container defined in elas3d_mod.
      ! Stores phase labels, displacements, stiffness data, residuals,
      ! stress outputs, orientation data, and convergence histories. 
      type(elas3d_data_type) :: e3d
      
      
      ! ======================================================================
      ! ---- 2. Memory Allocation --------------------------------------------
      ! Allocate all large arrays required by the solver. Array dimensions are
      ! determined from module parameters such as ns, ndof, nphase, and nphmax.
      !
      ! Major allocated fields:
      !   e3d%u, e3d%gb, e3d%b, e3d%h
      !     Displacement, gradient/residual, right-hand side, and CG search
      !     direction arrays.
      !
      !   e3d%cmod, e3d%dk
      !     Per-phase constitutive stiffness matrices and integrated local
      !     finite-element stiffness tensors.
      !
      !   e3d%ib, e3d%pix
      !     Periodic neighbor table and phase-label field.
      !
      !   e3d%orientation, e3d%rotatedstiffness
      !     Grain orientation data and rotated stiffness tensors for the
      !     polycrystalline case.
      !
      !   e3d%energies, e3d%ggs, e3d%cgiter
      !     Solver convergence-history arrays.
      ! ======================================================================
      allocate(phasemod(nphase,2), prob(nphmax))
      allocate(e3d%u(ns,ndof), e3d%gb(ns,ndof))
      allocate(e3d%b(ns,ndof), e3d%h(ns,ndof))

      
      allocate(e3d%cmod(nphmax,6,6))
      allocate(e3d%dk(nphmax,nnode_fe,ndof,nnode_fe,ndof))
      allocate(e3d%vm(ns), e3d%stress_field(ns,6))

      allocate(e3d%ib(ns,nfaces), e3d%pix(ns))
      allocate(e3d%orientation(3*(nphase-1)))
      allocate(e3d%rotatedstiffness(36*(nphase-1)))    
      allocate(e3d%energies(n_iter+1))
      allocate(e3d%ggs(n_iter+1))
      allocate(e3d%cgiter(n_iter+1))

      ! ======================================================================
      ! ---- 3. Variable Initialization --------------------------------------
      ! Initialize applied strain components, averaged stress/strain outputs,
      ! energy offset, and CG convergence-history arrays.
      !
      ! The actual prescribed macroscopic strains are assigned later from the
      ! module parameters aml_exx, aml_eyy, aml_ezz, aml_exz, aml_eyz, aml_exy.
      ! ======================================================================
      e3d%exx = 0.d0;      e3d%eyy = 0.d0;      e3d%ezz = 0.d0;
      e3d%exz = 0.d0;      e3d%eyz = 0.d0;      e3d%exy = 0.d0;
      e3d%C   = 0.d0
      e3d%strxx = 0.d0;    e3d%stryy = 0.d0;    e3d%strzz = 0.d0;
      e3d%strxz = 0.d0;    e3d%stryz = 0.d0;    e3d%strxy = 0.d0;
      e3d%sxx = 0.d0;      e3d%syy = 0.d0;      e3d%szz = 0.d0;
      e3d%sxz = 0.d0;      e3d%syz = 0.d0;      e3d%sxy = 0.d0;
      
      ! Initialize CG convergence-history arrays.
      e3d%energies = 0.d0
      e3d%ggs      = 0.d0
      e3d%cgiter   = 0
      
      ! ======================================================================
      ! ---- 4. Open Main Output File and Start Timer -------------------------
      ! Open the main text output file and start the wall-clock timer using
      ! system_clock.
      ! ======================================================================  
      open(unit=7, file='output_SS316L_polycrystal.txt', status='replace', iostat=ios)
      
      ! Start wall-clock timing.
      call system_clock(t1, tc)  
      
      ! Report grid dimensions to both the output file and console.
      write(7, '(a,i4,a,i4,a,i4,a,i12)') 'nx = ', nx, ', ny =', ny, ', nz =', nz, ', ns =', ns
      print '(a,i4,a,i4,a,i4,a,i12)', 'nx = ', nx, ', ny =', ny, ', nz =', nz, ', ns =', ns


      ! ======================================================================
      ! ---- 5. Initialize Isotropic Material Properties ----------------------
      ! Set isotropic Young's modulus and Poisson's ratio values for phases
      ! that are treated isotropically.
      !
      ! If flag_m == 0:
      !   All solid grain phases are assigned the isotropic equivalent matrix
      !   properties E0 and nu0.
      !
      ! If flag_m == 1:
      !   Solid grain phases are handled later using rotated anisotropic
      !   single-crystal stiffness tensors. Only the pore/void phase properties
      !   are initialized here.
      !
      ! phasemod(:,1):
      !   Initially stores Young's modulus E. It is converted below to bulk
      !   modulus K.
      !
      ! phasemod(:,2):
      !   Initially stores Poisson's ratio nu. It is converted below to shear
      !   modulus G.
      ! ======================================================================
      ! phasemod(i,1) temporarily stores E; phasemod(i,2) stores nu.
      if (flag_m == 0) then
        do i = 1, nphase-1
           phasemod(i,1) = E0   ! matrix
           phasemod(i,2) = nu0
        end do
      else
	    phasemod(nphase,1) = E1  ! Inclusions
		phasemod(nphase,2) = nu1    
      end if
      
      ! Convert E and nu to bulk modulus K and shear modulus G.
      ! After this block:
      !   phasemod(:,1) = K
      !   phasemod(:,2) = G
	  if (flag_m == 0) then
        do i = 1, nphase-1
          youngs = phasemod(i,1)
          phasemod(i,1) = phasemod(i,1) / 3.d0 / (1.d0 - 2.d0*phasemod(i,2)) ! K = E / 3(1-2v)
          phasemod(i,2) = youngs / 2.d0 / (1.d0 + phasemod(i,2))             ! G = E / 2(1+v) 
        end do
      else
	    youngs = phasemod(nphase,1)
		phasemod(nphase,1) = phasemod(nphase,1) / 3.d0 / (1.d0 - 2.d0*phasemod(nphase,2)) ! K = E/(3*(1-2*nu))
		phasemod(nphase,2) = youngs / 2.d0 / (1.d0 + phasemod(nphase,2))     ! G = E/(2*(1+nu))     
      end if
		
      ! ======================================================================
      ! ---- 6. Build Periodic Neighbor Table --------------------------------
      ! Construct the 27-point periodic neighbor stencil for every grid site.
      !
      ! e3d%ib(m,n) stores the 1D index of neighbor n associated with grid
      ! site m. Boundary indices are wrapped periodically in x, y, and z.
      !
      ! The 27 entries correspond to the local 3x3x3 neighborhood used by the
      ! finite-element stencil and matrix-free operations.
      ! ======================================================================
      ! Define relative coordinate shifts for the 27-neighbor stencil.
      im(1:8) = [0, 1, 1, 1, 0, -1, -1, -1]
      jm(1:8) = [1, 1, 0,-1,-1, -1,  0,  1]

      do n = 1, 8
        km(n)    = 0
        km(n+8)  = -1
        km(n+16) = 1
        im(n+8)  = im(n)
        im(n+16) = im(n)
        jm(n+8)  = jm(n)
        jm(n+16) = jm(n)
      end do

      im(25:27) = 0
      jm(25:27) = 0
      km(25) = -1; km(26) = 1; km(27) = 0      

      ! Construct the neighbor table using 1D linear indices.
      ! e3d%ib(m,n) gives the 1D index of the nth periodic neighbor of site m.
	  
      ! Number of grid sites in one x-y plane, used for 3D-to-1D indexing.
      nxy = nx * ny

      !$omp parallel do collapse(3) schedule(guided) default(shared) &
      !$omp private(i,j,k,n,m,m1,i1,j1,k1)	  
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
		  
            ! Linear index of the current grid site.
            m = nxy * (k-1) + nx * (j-1) + i
			
            do n = 1, 27
              ! Neighbor coordinates with periodic wrapping
              i1 = i + im(n)
              j1 = j + jm(n)
              k1 = k + km(n)
              ! Apply periodic wrapping				  
              if (i1 < 1)  i1 = i1 + nx
              if (i1 > nx) i1 = i1 - nx
              if (j1 < 1)  j1 = j1 + ny
              if (j1 > ny) j1 = j1 - ny
              if (k1 < 1)  k1 = k1 + nz
              if (k1 > nz) k1 = k1 - nz
			  
              ! Compute 1D index of the wrapped neighbor	
              m1 = nxy * (k1-1) + nx*(j1-1) + i1
			  
			  ! Store in neighbor table
              e3d%ib(m,n) = m1
            end do
          end do
        end do
      end do
      !$omp end parallel do
      
      ! ======================================================================
      ! ---- 7. Microstructure Pre-processing and FEM Setup -------------------
      ! Process each microstructure realization. In the current configuration,
      ! npoints = 1, so only one HDF5 input microstructure is analyzed.
      !
      ! For each microstructure:
      !   - read phase labels and orientations,
      !   - build the active-node mask,
      !   - compute phase volume fractions,
      !   - assign macroscopic strain loading,
      !   - assemble stiffness data and force terms,
      !   - initialize the displacement field.
      ! ======================================================================
      npoints = 1  ! Number of distinct microstructures to process
      do micro = 1, npoints
        ! Read phase labels and grain orientations from the HDF5 input file.
        ! ppixel_hdf5 populates e3d%pix and e3d%orientation.
        call ppixel_hdf5(e3d)          
      
        print *, 'Unique phases in pix:', minval(e3d%pix), maxval(e3d%pix)
        do i = 1, nphase
          write(7, '(A,I7,A,F20.6,A,F20.6)') 'Phase ',i,' bulk=',phasemod(i,1),' shear=',phasemod(i,2)
        end do        

        ! Allocate and build active-node mask.
        ! Active nodes are connected to at least one non-pore element.
       allocate(e3d%is_active(ns))
       call build_active_mask(e3d)

       print *, 'Active nodes  :', count(e3d%is_active)
       print *, 'Inactive nodes:', count(.not. e3d%is_active)
        
        ! Compute and write phase volume fractions.
		call assig(e3d, prob)
        do i = 1, nphase
          write(7, '(A,I7,A,F20.6)') 'Volume fraction of phase ',i,' is ',prob(i)
        end do


        ! Assign prescribed macroscopic strain components.
        !
        ! exx, eyy, and ezz are normal strain components.
        ! exz, eyz, and exy are tensorial shear strain components.
        ! Engineering shear strains are therefore:
        !   gamma_xz = 2*exz
        !   gamma_yz = 2*eyz
        !   gamma_xy = 2*exy
        e3d%exx = aml_exx
        e3d%eyy = aml_eyy
        e3d%ezz = aml_ezz
        e3d%exz = aml_exz
        e3d%eyz = aml_eyz
        e3d%exy = aml_exy   
     
        write(7, *) 'Applied strains: exx eyy ezz gamma_xz gamma_yz gamma_xy'
        write(7, *) e3d%exx, e3d%eyy, e3d%ezz, 2.d0*e3d%exz, 2.d0*e3d%eyz, 2.d0*e3d%exy
        
        print *, 'Applied strains: exx eyy ezz gamma_xz gamma_yz gamma_xy'
        print '(6E15.5)', e3d%exx, e3d%eyy, e3d%ezz, 2.d0*e3d%exz, 2.d0*e3d%eyz, 2.d0*e3d%exy     


        ! Assemble material stiffness data, local finite-element stiffness
        ! matrices, the equivalent force vector b, and the energy offset C
        ! required by the energy functional:
        !
        !   U = 0.5*u^T*K*u + b^T*u + C
        !
        ! For flag_m == 0, phmod contains isotropic solid/void phase moduli.
        ! For flag_m == 1, phmod contains only the isotropic pore/void moduli;
        ! solid grain stiffnesses are assigned from rotated single-crystal
        ! tensors inside femat_OpenMP.
		if (flag_m == 0) then
		  phmod(1:2,1:2) = phasemod(1:2,1:2)
		else
          phmod(1:1,1:2) = phasemod(nphase:nphase,1:2)    
        end if 
        ! OpenMP-parallelized finite-element assembly routine.		
        call femat_OpenMP(e3d, phmod)
        ! Orientation data and rotated stiffness tensors are no longer needed
        ! after e3d%cmod and e3d%dk have been assembled.	
        deallocate(e3d%orientation, e3d%rotatedstiffness)
      
        ! Initialize the displacement field using the homogeneous displacement
        ! associated with the prescribed macroscopic strain tensor:
        !
        !   u_x = exx*x + exy*y + exz*z
        !   u_y = exy*x + eyy*y + eyz*z
        !   u_z = exz*x + eyz*y + ezz*z
        !
        ! Inactive pore-only nodes are initialized to zero displacement.       
        !$omp parallel do collapse(3) schedule(guided) default(shared) private(i,j,k,m,x,y,z)		
        do k = 1, nz
		  do j = 1, ny
            do i = 1, nx
			  m = nxy*(k-1)+nx*(j-1)+i
              if (.not. e3d%is_active(m)) then
                e3d%u(m,:) = 0.d0   ! inactive pore-only node
              else			  
                x = dble(i-1)
			    y = dble(j-1)
			    z = dble(k-1)
		   
                e3d%u(m,1) = x * e3d%exx + y * e3d%exy + z * e3d%exz
                e3d%u(m,2) = x * e3d%exy + y * e3d%eyy + z * e3d%eyz
                e3d%u(m,3) = x * e3d%exz + y * e3d%eyz + z * e3d%ezz
			  end if	
            end do
          end do
        end do
		!$omp end parallel do
     
        
      ! ======================================================================
      ! ---- 8. PCG Energy-Minimization Loop ---------------------------------
      ! Compute the initial energy/residual and solve for the equilibrium
      ! displacement field using preconditioned conjugate-gradient minimization.
      !
      ! The residual metric gg is the preconditioned gradient inner product
      ! used by the solver. The relative residual norm is sqrt(gg/gg0).
      ! ====================================================================== 
        ltot = 0
        lhist = 1
        
        ! Compute initial energy and gradient/residual for the homogeneous
        ! displacement-field guess.
        call energy_OpenMP(e3d, utot)
        
        
        ! Compute the initial residual norm measure used for relative convergence.
        gg = sum(e3d%gb(:,:)**2)
        gg0 = gg    ! Store initial norm as baseline
		
        write(7,*) 'Initial energy=', utot, 'gg=', gg
        write(7,*) '   '
        call flush(7)
        print '(A,E12.5,A,E12.5)', 'Initial energy=', utot, ', Initial gg=', gg
        print *, '   '
        
        e3d%energies(lhist) = utot
        e3d%ggs(lhist)      = sqrt(gg/gg0) ! Store initial convergence history entry.
        e3d%cgiter(lhist)   = 0
        
        ! Open high-level convergence-history file.
        open(unit=99, file='cg_history_polycrystal.txt', status='replace')
        write(99, '(I8,1X,1PE18.10,1X,1PE18.10)') e3d%cgiter(lhist), e3d%energies(lhist), e3d%ggs(lhist)
        
        ! Outer PCG restart/load-step loop.
        ! For this linear problem, kmax is normally set to 1.	
        do kkk=1, kmax
           ! Run the OpenMP-parallelized PCG solver.
           call dembx_OpenMP(e3d, Lstep, gg, gg0, kkk) 
                       

           ltot = ltot + Lstep
           lhist = lhist + 1
           
           ! Recompute total energy after the PCG call. If convergence was reached,
           ! this is the final energy; otherwise, it provides an intermediate
           ! diagnostic value. 
           call energy_OpenMP(e3d, utot)
         

           write(7,*) 'Energy=', utot, 'gg=', gg
           write(7,*) 'Number of conjugate steps=', ltot
           call flush(7) 
           print '(A,E12.5,A,E12.5,A,E12.5)', 'Energy=', utot, ', gg=', gg, ', rel. residual norm=', sqrt(gg/gg0)
           print *, 'Number of conjugate steps=', ltot
           
           ! Save high-level CG convergence history.
           e3d%energies(lhist) = utot
           e3d%ggs(lhist)      = sqrt(gg/gg0)
           e3d%cgiter(lhist)   = ltot           
           
           write(99, '(I8,1X,1PE18.10,1X,1PE18.10)') e3d%cgiter(lhist), e3d%energies(lhist), e3d%ggs(lhist)

                      
           ! Exit the outer loop if the relative residual tolerance is met.
		   rel_res = sqrt(gg/gg0)
           if (rel_res < tol) exit
           
           
           ! Compute and output intermediate volume-averaged stress/strain
           ! if another outer iteration would be required.
           call stress_OpenMP(e3d)
           

           write(7,*) 'stresses: xx, yy, zz, xz, yz, xy'
           write(7,'(6E15.6)') e3d%strxx, e3d%stryy, e3d%strzz, e3d%strxz, e3d%stryz, e3d%strxy
           write(7,*) 'strains: xx, yy, zz, xz, yz, xy'
           write(7,'(6E15.6)') e3d%sxx, e3d%syy, e3d%szz, e3d%sxz, e3d%syz, e3d%sxy
           write(7,*) '     '
           call flush(7)
           
        end do
        
        ! Close file for CG history
        close(99)

        ! Compute final volume-averaged stresses and strains using the
        ! converged displacement field.
        call stress_OpenMP(e3d)


        write(7,*) 'stresses: xx, yy, zz, xz, yz, xy'
        write(7,'(6E15.6)') e3d%strxx, e3d%stryy, e3d%strzz, e3d%strxz, e3d%stryz, e3d%strxy
        write(7,*) 'strains: xx, yy, zz, xz, yz, xy'
        write(7,'(6E15.6)') e3d%sxx, e3d%syy, e3d%szz, e3d%sxz, e3d%syz, e3d%sxy

        
      end do
      
      print *, 'End of CG'
      
      ! ======================================================================
      ! ---- 9. Post-processing and Teardown ----------------------------------
      ! Stop the timer, report runtime, compute full-field stress/von Mises
      ! output, and deallocate all dynamically allocated arrays.
      ! ======================================================================
      call system_clock(t2)
      elapsed_time = dble(t2-t1)/dble(tc)
      print *, ''
      print *, 'Calculation complete!'
      print *, 'Total calculation time: ', elapsed_time, ' seconds'   
	  write(7,*) ''
	  write(7,*) 'Calculation complete!'
      write(7,*) 'Total calculation time: ', elapsed_time, ' seconds'   	  
      
      ! Compute and write full-field stress and von Mises stress after
      ! convergence using the final displacement field.
      call stress_fullfield_OpenMP(e3d)    

      
      ! Deallocate all dynamically allocated arrays before program exit.
      deallocate(phasemod, prob)
      deallocate(e3d%u, e3d%gb, e3d%b)
      deallocate(e3d%h)
      deallocate(e3d%cmod, e3d%dk)
      deallocate(e3d%vm, e3d%stress_field )
      deallocate(e3d%ib, e3d%pix)
      deallocate(e3d%energies, e3d%ggs, e3d%cgiter)
      deallocate(e3d%is_active)
      
      end program elas3dxtal    
    
!=========================================================================
      subroutine build_active_mask(e3d)
      ! ==============================================================================
      ! SUBROUTINE: build_active_mask
      ! ------------------------------------------------------------------------------
      ! PURPOSE:
      !   Builds a logical active-node mask used to exclude mechanically inactive
      !   pore-only nodes from the displacement solve.
      !
      !   A node is marked active if it belongs to at least one neighboring solid
      !   finite element. A node is marked inactive only if all eight hexahedral
      !   elements sharing that node are pore/void elements.
      !
      ! WHY THIS IS NEEDED:
      !   The pore phase has zero stiffness. Nodes surrounded entirely by pore/void
      !   elements have no mechanical stiffness contribution and can create singular
      !   or poorly conditioned degrees of freedom in the PCG solve. These nodes are
      !   therefore constrained to zero displacement/residual in later routines.
      !
      ! INPUT/OUTPUT:
      !   e3d%ib
      !     Periodic 27-neighbor table. Used to locate the eight element origins
      !     surrounding each node.
      !
      !   e3d%pix
      !     Phase label array. e3d%pix(m) gives the phase ID associated with grid
      !     site/element m.
      !
      !   e3d%is_active
      !     Logical output mask. Must be allocated before calling this subroutine.
      !       .true.  = node is connected to at least one non-pore element
      !       .false. = node is connected only to pore/void elements
      !
      ! PHASE CONVENTION:
      !   PORE_PHASE is the phase ID assigned to pore/void material.
      !
      ! PERIODICITY:
      !   The neighbor lookup uses e3d%ib, so the active-node check is consistent
      !   with the periodic boundary conditions used throughout the solver.
      ! ==============================================================================	  
      use elas3d_mod
      implicit none
      type(elas3d_data_type), intent(inout) :: e3d

      ! Neighbor-table slots corresponding to the eight hexahedral elements
      ! that share the current node m.
      !
      ! In this code, each grid site can be interpreted as the origin/corner of
      ! an element, and the 27-neighbor table stores nearby element origins using
      ! periodic indexing. The following eight slots identify the elements whose
      ! local 8-node connectivity includes node m.
      !
      ! This ordering is consistent with the element-node mapping used in the
      ! stiffness assembly and matrix-free stencil operations.
      integer, parameter :: elem_nbrs(8) = [27, 7, 6, 5, 25, 15, 14, 13]
      integer :: m, q, nbr_elem
      logical :: any_solid

      ! Loop over all grid nodes. Each node is active if at least one of the
      ! eight surrounding elements is not the pore/void phase.
      !$omp parallel do schedule(guided) private(m, q, nbr_elem, any_solid)
      do m = 1, ns
        any_solid = .false.
		
        ! Check the eight periodically wrapped elements that share node m.		
        do q = 1, 8
          nbr_elem = e3d%ib(m, elem_nbrs(q))
		  
          ! If any adjacent element is solid, retain this node in the solve.		  
          if (e3d%pix(nbr_elem) /= PORE_PHASE) then
            any_solid = .true.
            exit
          end if
        end do
		
        ! Store active/inactive status for later use in assembly, energy,
        ! residual, and PCG update routines.		
        e3d%is_active(m) = any_solid
      end do
      !$omp end parallel do

      end subroutine build_active_mask
!-------------------------------------------------------------------------
      subroutine femat_OpenMP(e3d, phasemod)    
      ! ==============================================================================
      ! SUBROUTINE: femat_OpenMP
      ! ------------------------------------------------------------------------------
      ! PURPOSE:
      !   Builds the phase stiffness matrices and integrated finite-element stiffness
      !   tensors required by the energy-based PCG solver.
      !
      !   This routine:
      !     1. Constructs 6x6 constitutive stiffness matrices for each phase.
      !     2. Computes the 8-node hexahedral element stiffness matrix for each phase.
      !     3. Applies periodic macroscopic-strain jump corrections on boundary
      !        faces, edges, and corners.
      !     4. Assembles the equivalent force vector b and scalar energy offset C
      !        used in the total energy functional:
      !
      !          U(u) = 0.5*u^T*K*u + b^T*u + C
      !
      !   The global stiffness matrix K is not explicitly assembled. Instead, the
      !   local phase stiffness tensors e3d%dk are later used by matrix-free stencil
      !   operations in energy_OpenMP and dembx_OpenMP.
      !
      ! ARGUMENTS:
      !   e3d
      !     Main solver data structure.
      !
      !     Input components:
      !       e3d%pix
      !         Phase label at each grid site/element.
      !
      !       e3d%ib
      !         Periodic 27-neighbor connectivity table.
      !
      !       e3d%exx, e3d%eyy, e3d%ezz, e3d%exz, e3d%eyz, e3d%exy
      !         Prescribed macroscopic strain components. Shear components are
      !         tensorial shear strains.
      !
      !       e3d%orientation
      !         Rodrigues orientation vectors for solid grain phases. Used only
      !         when flag_m == 1.
      !
      !       e3d%rotatedstiffness
      !         Workspace for rotated grain stiffness tensors. Populated by
      !         rotatestiffnessby when flag_m == 1.
      !
      !       e3d%is_active
      !         Active-node mask. Inactive pore-only nodes receive zero forcing.
      !
      !     Output components:
      !       e3d%cmod
      !         6x6 constitutive stiffness matrix for each phase.
      !
      !       e3d%dk
      !         Integrated local finite-element stiffness tensor for each phase.
      !
      !       e3d%b
      !         Equivalent force vector associated with imposed periodic strain jumps.
      !
      !       e3d%C
      !         Constant energy offset associated with imposed periodic strain jumps.
      !
      !   phasemod
      !     Isotropic material-property work array.
      !
      !     Convention:
      !       phasemod(:,1) = bulk modulus K
      !       phasemod(:,2) = shear modulus G
      !
      !     In the current polycrystal case, phasemod(1,:) is used for the
      !     pore/void phase. Solid grain stiffnesses are generated from rotated
      !     anisotropic single-crystal tensors.
      !
      ! NUMERICAL FORMULATION:
      !   - Element type:
      !       8-node trilinear hexahedral element.
      !
      !   - Integration:
      !       3 x 3 x 3 Simpson/Newton-Cotes integration over the unit voxel.
      !       The weights are 1, 4, 16, and 64 depending on the number of midpoint
      !       coordinates, with final normalization by 216.
      !
      !   - Constitutive notation:
      !       Voigt ordering is:
      !         1 = xx, 2 = yy, 3 = zz, 4 = xz, 5 = yz, 6 = xy
      !
      !   - Boundary conditions:
      !       Periodic displacement jumps are applied on upper faces, edges, and
      !       the upper corner to represent the prescribed macroscopic strain.
      !
      ! PARALLELIZATION:
      !   OpenMP is used for phase stiffness integration and boundary-force
      !   assembly. Atomic updates are used when multiple threads may contribute
      !   to the same global force-vector entry.
      ! ==============================================================================
      use elas3d_mod
      implicit none

      ! ======================================================================
      ! ---- 1. Subroutine Arguments and Local Variables ----------------------
      ! ======================================================================

      ! Main solver state structure.
      type(elas3d_data_type), intent(inout) :: e3d
	  
      ! Isotropic phase properties:
      !   phasemod(:,1) = bulk modulus K
      !   phasemod(:,2) = shear modulus G
      !
      ! In polycrystal mode, only phasemod(1,:) is used here for the pore/void
      ! phase; solid grain stiffnesses are assigned from rotated anisotropic
      ! stiffness tensors.	  
      real(dp), intent(in) :: phasemod(2,2)

      ! Loop counters and indexing variables.
      integer :: i, j, k, l, m, n
      integer :: ijk, mm, nn, ii, jj, kk, ll
      integer :: nxy, i3, i8, j3, k3, n8, m3, m8
	  
      ! Local finite-element arrays.
      ! dndx, dndy, dndz:
      !   Shape-function derivatives at the current integration point.
      !
      ! ck, cmu:
      !   Volumetric and deviatoric basis matrices for isotropic stiffness.
      !
      ! g:
      !   3D Simpson integration weights.
      !
      ! es:
      !   Strain-displacement matrix B stored as
      !   es(strain_component, local_node, displacement_component).
      !
      ! delta:
      !   Prescribed periodic displacement jump for local element nodes.	  
      integer :: is(nnode_fe), ind
      real(dp) :: dndx(nnode_fe), dndy(nnode_fe), dndz(nnode_fe)
      real(dp) :: ck(6,6), cmu(6,6), g(ngauss,ngauss,ngauss)
      real(dp) :: es(6,nnode_fe,ndof), delta(nnode_fe,ndof)
      real(dp) :: left(6), right(6)
      real(dp) :: sumval, C
      real(dp) :: x, y, z
      real(dp) :: exx, eyy, ezz, exz, eyz, exy
	  
      ! Local work arrays assembled in this subroutine and copied into e3d.
      ! dk stores phase-wise local element stiffness matrices.
      ! b stores equivalent force-vector contributions from periodic strain jumps.
      real(dp), allocatable :: dk(:,:,:,:,:)
      real(dp), allocatable :: b(:,:)

      ! ======================================================================
      ! ---- 2. Memory Allocation and Initialization --------------------------
      ! Allocate local work arrays and zero all temporary/global output arrays
      ! before assembly.
      !
      ! dk and b are local to this subroutine. They are copied back to e3d%dk
      ! and e3d%b after assembly is complete.
      ! ======================================================================
      allocate(dk(nphmax,nnode_fe,ndof,nnode_fe,ndof))
	  allocate(b(ns,ndof))

      g = 0.d0;   es = 0.d0;   delta = 0.d0
      dk = 0.d0;  ck = 0.d0;   cmu = 0.d0
      left = 0.d0; right = 0.d0; b = 0.d0; C = 0.d0
      dndx = 0.d0; dndy = 0.d0; dndz = 0.d0

      e3d%cmod = 0.d0;   e3d%dk = 0.d0;   e3d%b = 0.d0
      e3d%C = 0.d0

      ! ======================================================================
      ! ---- 3. Finite-Element Node Mapping -----------------------------------
      ! Map the eight local nodes of a hexahedral voxel element to entries in
      ! the 27-point periodic neighbor table e3d%ib.
      !
      ! The ordering in is(:) must remain consistent with the local stiffness
      ! tensor dk and with the stencil expressions used in energy_OpenMP and
      ! dembx_OpenMP.
      ! ======================================================================
      is = (/27, 3, 2, 1, 26, 19, 18, 17/)
	  
      nxy = nx * ny


      ! ======================================================================
      ! ---- 4. Isotropic Stiffness Basis Matrices ----------------------------
      ! Construct the volumetric and deviatoric basis matrices used to form an
      ! isotropic 6x6 stiffness matrix from bulk modulus K and shear modulus G:
      !
      !   C = K*ck + G*cmu
      !
      ! ck:
      !   Volumetric contribution. The upper-left 3x3 block is all ones.
      !
      ! cmu:
      !   Deviatoric/shear contribution. The upper-left 3x3 block contains
      !   4/3 on the diagonal and -2/3 off diagonal. The shear entries
      !   cmu(4,4), cmu(5,5), and cmu(6,6) are one.
      !
      ! Voigt ordering:
      !   1=xx, 2=yy, 3=zz, 4=xz, 5=yz, 6=xy
      ! ======================================================================
      ck(1:3,1:3) = 1.d0
      do i = 1, 3
        cmu(i,i) = 4.d0/3.d0
        do j = 1, 3
          if (i /= j) cmu(i,j) = -2.d0/3.d0
        end do
      end do
      cmu(4,4) = 1.d0
      cmu(5,5) = 1.d0
      cmu(6,6) = 1.d0      

      ! ======================================================================
      ! ---- 5. Constitutive Matrix Assembly: e3d%cmod ------------------------
      ! Build the 6x6 stiffness matrix for each material phase.
      !
      ! If flag_m == 0:
      !   The solid matrix phases are treated isotropically using bulk and shear
      !   moduli from phasemod.
      !
      ! If flag_m == 1:
      !   Solid grain phases are treated as anisotropic crystals. The cubic
      !   single-crystal stiffness tensor is rotated into the sample frame for
      !   each grain using its Rodrigues orientation vector.
      !
      !   The final phase, nphase, is the pore/void phase and is assigned an
      !   isotropic zero-stiffness tensor from phasemod(1,:).
      !
      ! Isotropic stiffness construction:
      !   C = K*ck + G*cmu
      ! ======================================================================
	  if (flag_m == 0) then
        ! Isotropic mode.
        ! Construct isotropic stiffness matrices from K and G.
	    do k = 1,2; do j = 1,6; do i = 1,6;
		  e3d%cmod(k,i,j) = phasemod(k,1)*ck(i,j) + phasemod(k,2)*cmu(i,j)
        end do; end do; end do
		
      else
        ! Polycrystal mode.
        ! Rotate the single-crystal stiffness tensor into the global/sample
        ! frame for each grain phase.
		call rotatestiffnessby(e3d)	
      
        ! Assign the transformed anisotropic stiffness to each grain phase	
        !$omp parallel do collapse(3) schedule(guided) default(shared)
        do k = 1, nphase-1
		  do i = 1,6
		    do j = 1,6
              ! Map flattened rotated stiffness storage into e3d%cmod.
              ! Flattening convention:
              !   rotatedstiffness(36*(grain-1) + 6*(i-1) + j)		
			  e3d%cmod(k,i,j) = e3d%rotatedstiffness(36*(k-1)+(i-1)*6+j)
			end do
          end do
        end do
		!$omp end parallel do
		
		
        ! Assign isotropic pore/void stiffness to the final phase.
        ! For the present model, E1 = 0 and nu1 = 0, so this phase has
        ! zero stiffness.
        do j = 1, 6
		  do i = 1, 6
            e3d%cmod(nphase,i,j) = phasemod(1,1)*ck(i,j) + phasemod(1,2)*cmu(i,j)
          end do
		end do
      end if

      ! ======================================================================
      ! ---- 6. Numerical Integration Weights ---------------------------------
      ! Set up the 3 x 3 x 3 Simpson/Newton-Cotes integration weights over the
      ! unit voxel.
      !
      ! Each coordinate direction uses weights [1, 4, 1]. The 3D weight is the
      ! product of the three 1D weights, giving possible values:
      !   1, 4, 16, and 64.
      !
      ! The final normalization by 216 appears in the stiffness accumulation:
      !   216 = 6^3
      ! ======================================================================
      do k3 = 1, ngauss; do j3 = 1, ngauss; do i3 = 1, ngauss;
        n = 0
        if(i3==2) n = n+1
        if(j3==2) n = n+1
        if(k3==2) n = n+1
        g(i3,j3,k3) = 4.0 ** n
      end do; end do; end do

      ! ======================================================================
      ! ---- 7. Element Stiffness Integration ---------------------------------
      ! For each phase, integrate the local element stiffness matrix:
      !
      !   k_e = integral( B^T * C * B dV )
      !
      ! over the unit hexahedral voxel using the 3 x 3 x 3 Simpson integration
      ! rule defined above.
      !
      ! The resulting tensor is stored as:
      !   dk(phase, local_node_i, dof_i, local_node_j, dof_j)
      !
      ! Parallelization note:
      !   The outer loop is over phase index ijk. Each thread writes to a
      !   different dk(ijk,:,:,:,:), so no atomic updates are needed here.
      ! ======================================================================
      !$omp parallel default(shared) private(k3,j3,i3,x,y,z,dndx,dndy,dndz,es,n8,mm,nn,ii,jj,left,right,sumval)
      !$omp do schedule(guided)
      do ijk = 1, nphase
	    do k3 = 1, ngauss
          do j3 = 1, ngauss
            do i3 = 1, ngauss;
              ! Local coordinates of the current integration point in [0,1]^3.			
              x = dble(i3-1)/2.d0
              y = dble(j3-1)/2.d0
              z = dble(k3-1)/2.d0

              ! Derivatives of the 8 trilinear shape functions with respect to
              ! the local unit-cell coordinates x, y, and z.
              dndx(1) = -(1.d0-y)*(1.d0-z)
              dndx(2) =  (1.d0-y)*(1.d0-z)
              dndx(3) =         y*(1.d0-z)
              dndx(4) =        -y*(1.d0-z)
              dndx(5) = -(1.d0-y)*z
              dndx(6) =  (1.d0-y)*z
              dndx(7) =         y*z
              dndx(8) =        -y*z

              dndy(1) = -(1.d0-x)*(1.d0-z)
              dndy(2) =        -x*(1.d0-z)
              dndy(3) =         x*(1.d0-z)
              dndy(4) =  (1.d0-x)*(1.d0-z)
              dndy(5) = -(1.d0-x)*z
              dndy(6) =        -x*z
              dndy(7) =         x*z
              dndy(8) =  (1.d0-x)*z

              dndz(1) = -(1.d0-x)*(1.d0-y)
              dndz(2) =        -x*(1.d0-y)
              dndz(3) =        -x*y
              dndz(4) = -(1.d0-x)*y
              dndz(5) =  (1.d0-x)*(1.d0-y)
              dndz(6) =         x*(1.d0-y)
              dndz(7) =         x*y
              dndz(8) =  (1.d0-x)*y

              ! Construct the strain-displacement matrix B.
              ! Voigt strain ordering:
              !   1=xx, 2=yy, 3=zz, 4=xz, 5=yz, 6=xy
              es = 0.d0
              do n8 = 1, nnode_fe
                es(1,n8,1) = dndx(n8)
                es(2,n8,2) = dndy(n8)
                es(3,n8,3) = dndz(n8)
                es(4,n8,1) = dndz(n8)
                es(4,n8,3) = dndx(n8)
                es(5,n8,2) = dndz(n8)
                es(5,n8,3) = dndy(n8)
                es(6,n8,1) = dndy(n8)
                es(6,n8,2) = dndx(n8)
              end do

              ! Perform matrix multiplication: B^T * D * B
              ! Accumulate into local stiffness matrix `dk` using Gauss weights
              do mm = 1, ndof
                do nn = 1, ndof
                  do ii = 1, nnode_fe
                    left = es(:, ii, mm)
                    do jj = 1, nnode_fe
                      right = es(:, jj, nn)
                      sumval = dot_product(left, matmul(e3d%cmod(ijk,:,:), right))
                      ! Accumulate B_i^T * C * B_j weighted by the integration
                      ! weight and normalized by the Simpson-rule factor 216.					  
                      dk(ijk,ii,mm,jj,nn) = dk(ijk,ii,mm,jj,nn) + g(i3,j3,k3)*sumval/216.d0
                    end do
                  end do
                end do
              end do

            end do
          end do
        end do
      end do
	  !$omp end do
	  !$omp end parallel

      ! Copy the completed phase-wise local stiffness tensors into the solver state.
      e3d%dk = dk

      ! ======================================================================
      ! ---- 8. Periodic Macroscopic-Strain Corrections -----------------------
      ! Apply the contribution of prescribed macroscopic strain jumps to the
      ! energy functional.
      !
      ! For periodic boundary conditions with an imposed average strain, nodes
      ! on opposite sides of the domain differ by a known displacement jump:
      !
      !   delta_u = E_macro * L
      !
      ! Rather than explicitly duplicating boundary degrees of freedom, this
      ! routine accounts for those jumps by assembling:
      !
      !   b : equivalent force vector
      !   C : constant energy offset
      !
      ! Boundary corrections are applied on:
      !   - upper x, y, and z faces,
      !   - upper edges where two periodic jumps combine,
      !   - upper corner where all three periodic jumps combine.
      ! ======================================================================
      ! Work copies of equivalent force vector and energy offset.	  
      b = e3d%b
      C = e3d%C
	  
      ! Unpack prescribed macroscopic strain components.
      exx = e3d%exx; eyy = e3d%eyy; ezz = e3d%ezz
      exz = e3d%exz; eyz = e3d%eyz; exy = e3d%exy

      ! -------- Face: x = nx ------------------------------------------------
      ! Define displacement jumps for local element nodes crossing the periodic
      ! x-boundary.
      delta = 0.d0
      do i8 = 1, nnode_fe;
        if (i8==2 .or. i8==3 .or. i8==6 .or. i8==7) then
          delta(i8,1) = exx*nx
          delta(i8,2) = exy*nx
          delta(i8,3) = exz*nx
        end if
      end do

      ! Apply force corrections due to displacement jump
      !$omp parallel do collapse(2) schedule(guided) private(j,k,m,nn,mm,sumval,m3,m8) reduction(+:C)	  
      do j = 1, ny-1
        do k = 1, nz-1
		  ! 1. Calculate Element Index
          m = nxy*(k-1) + j*nx
		  
          do nn = 1, ndof
            do mm = 1, nnode_fe
			
			  ! 2. Compute Contribution (sumval)
              ! This is purely local work, no race conditions here.
              sumval = 0.d0
              do m3 = 1, ndof
                do m8 = 1, nnode_fe
				  ! NOTE: Ensure 'dk' dimensions are (m8, m3, mm, nn, pix) 
                  ! for fastest Fortran access (Column-Major).
                  sumval = sumval + delta(m8,m3) * dk(e3d%pix(m),m8,m3,mm,nn)
				  
				  ! Scalar reduction on C is safe and fast
                  C = C + 0.5d0 * delta(m8,m3) * dk(e3d%pix(m),m8,m3,mm,nn) * delta(mm,nn)
                end do
              end do
			  
			  ! 3. Pre-calculate Index to keep Atomic block tight
              ind = e3d%ib(m,is(mm))

              ! 4. Critical Update			  
              !$omp atomic update
			  b(ind,nn) = b(ind,nn) + sumval

            end do
          end do
        end do
      end do
      !$omp end parallel do

      ! -------- Face: y = ny --------
      delta = 0.d0
      do i8 = 1, nnode_fe;
        if (i8==3 .or. i8==4 .or. i8==7 .or. i8==8) then
          delta(i8,1) = exy*ny
          delta(i8,2) = eyy*ny
          delta(i8,3) = eyz*ny
        end if
      end do;

      !$omp parallel do collapse(2) schedule(guided) private(i,k,m,nn,mm,sumval,m3,m8) reduction(+:C)	  
      do i = 1, nx-1
        do k = 1, nz-1
          m = nxy*(k-1) + nx*(ny-1) + i
          do nn = 1, ndof
            do mm = 1, nnode_fe
              sumval = 0.d0
              do m3 = 1, ndof
                do m8 = 1, nnode_fe
                  sumval = sumval + delta(m8,m3) * dk(e3d%pix(m),m8,m3,mm,nn)
                  C = C + 0.5d0 * delta(m8,m3) * dk(e3d%pix(m),m8,m3,mm,nn) * delta(mm,nn)
                end do
              end do
			  ind = e3d%ib(m,is(mm))
              !$omp atomic update
			  b(ind,nn) = b(ind,nn) + sumval
            end do
          end do
        end do
      end do
      !$omp end parallel do

      ! -------- Face: z = nz --------
      delta = 0.d0
      do i8 = 1, nnode_fe;
        if (i8==5 .or. i8==6 .or. i8==7 .or. i8==8) then
          delta(i8,1) = exz*nz
          delta(i8,2) = eyz*nz
          delta(i8,3) = ezz*nz
        end if
      end do;

      !$omp parallel do collapse(2) schedule(guided) private(i,j,m,nn,mm,sumval,m3,m8) reduction(+:C)	  
      do i = 1, nx-1
        do j = 1, ny-1
          m = nxy*(nz-1) + nx*(j-1) + i
          do nn = 1, ndof
            do mm = 1, nnode_fe
              sumval = 0.d0
              do m3 = 1, ndof
                do m8 = 1, nnode_fe
                  sumval = sumval + delta(m8,m3) * dk(e3d%pix(m),m8,m3,mm,nn)
                  C = C + 0.5d0 * delta(m8,m3) * dk(e3d%pix(m),m8,m3,mm,nn) * delta(mm,nn)
                end do
              end do
              ind = e3d%ib(m,is(mm))
              !$omp atomic update
			  b(ind,nn) = b(ind,nn) + sumval
            end do
          end do
        end do
      end do
      !$omp end parallel do

      ! -------- Edge: x = nx, y = ny ----------------------------------------
      ! Nodes on this edge receive the combined x- and y-periodic displacement
      ! jumps.
      delta = 0.d0
      do i8 = 1, nnode_fe;
        if (i8==2 .or. i8==6) then
          delta(i8,1) = exx*nx
          delta(i8,2) = exy*nx
          delta(i8,3) = exz*nx
        end if
        if (i8==4 .or. i8==8) then
          delta(i8,1) = exy*ny
          delta(i8,2) = eyy*ny
          delta(i8,3) = eyz*ny
        end if
        if (i8==3 .or. i8==7) then
          delta(i8,1) = exy*ny + exx*nx
          delta(i8,2) = eyy*ny + exy*nx
          delta(i8,3) = eyz*ny + exz*nx
        end if
      end do;

      !$omp parallel do schedule(guided) private(k,m,nn,mm,sumval,m3,m8) reduction(+:C)	  
      do k = 1, nz-1
        m = nxy*k
        do nn = 1, ndof
          do mm = 1, nnode_fe
            sumval = 0.d0
            do m3 = 1, ndof
              do m8 = 1, nnode_fe
                sumval = sumval + delta(m8,m3) * dk(e3d%pix(m),m8,m3,mm,nn)
                C = C + 0.5d0 * delta(m8,m3) * dk(e3d%pix(m),m8,m3,mm,nn)*delta(mm,nn)
              end do
            end do
            ind = e3d%ib(m,is(mm))
            !$omp atomic update
            b(ind,nn) = b(ind,nn) + sumval
          end do
        end do
      end do
      !$omp end parallel do

      ! -------- Edge: x=nx, z=nz --------
      delta = 0.d0
      do i8 = 1, nnode_fe;
        if(i8==2 .or. i8==3) then
          delta(i8,1) = exx*nx
          delta(i8,2) = exy*nx
          delta(i8,3) = exz*nx
        end if
        if(i8==5 .or. i8==8) then
          delta(i8,1) = exz*nz
          delta(i8,2) = eyz*nz
          delta(i8,3) = ezz*nz
        end if
        if(i8==6 .or. i8==7) then
          delta(i8,1) = exz*nz + exx*nx
          delta(i8,2) = eyz*nz + exy*nx
          delta(i8,3) = ezz*nz + exz*nx
        end if
      end do;

      !$omp parallel do schedule(guided) private(j,m,nn,mm,sumval,m3,m8) reduction(+:C)	  
      do j = 1, ny-1
        m = nxy*(nz-1) + nx*(j-1) + nx
        do nn = 1, ndof
          do mm = 1, nnode_fe
            sumval = 0.d0
            do m3 = 1, ndof
              do m8 = 1, nnode_fe
                sumval = sumval + delta(m8,m3)*dk(e3d%pix(m),m8,m3,mm,nn)
                C = C + 0.5d0 * delta(m8,m3)*dk(e3d%pix(m),m8,m3,mm,nn)*delta(mm,nn)
              end do
            end do
            ind = e3d%ib(m,is(mm))
            !$omp atomic update
            b(ind,nn) = b(ind,nn) + sumval
          end do
        end do
      end do
      !$omp end parallel do

      ! -------- Edge: y=ny, z=nz --------
      delta = 0.d0
      do i8 = 1, nnode_fe;
        if(i8==5 .or. i8==6) then
          delta(i8,1) = exz*nz
          delta(i8,2) = eyz*nz
          delta(i8,3) = ezz*nz
        end if
        if(i8==3 .or. i8==4) then
          delta(i8,1) = exy*ny
          delta(i8,2) = eyy*ny
          delta(i8,3) = eyz*ny
        end if
        if(i8==7 .or. i8==8) then
          delta(i8,1) = exy*ny + exz*nz
          delta(i8,2) = eyy*ny + eyz*nz
          delta(i8,3) = eyz*ny + ezz*nz
        end if
      end do;

      !$omp parallel do schedule(guided) private(i,m,nn,mm,sumval,m3,m8) reduction(+:C)	  
      do i = 1, nx-1
        m = nxy*(nz-1) + nx*(ny-1) + i
        do nn = 1, ndof
          do mm = 1, nnode_fe
            sumval = 0.d0
            do m3 = 1, ndof
              do m8 = 1, nnode_fe
                sumval = sumval + delta(m8,m3)*dk(e3d%pix(m),m8,m3,mm,nn)
                C = C + 0.5d0 * delta(m8,m3)*dk(e3d%pix(m),m8,m3,mm,nn)*delta(mm,nn)
              end do
            end do
            ind = e3d%ib(m,is(mm))
            !$omp atomic update
            b(ind,nn) = b(ind,nn) + sumval
          end do
        end do
      end do
      !$omp end parallel do

      ! -------- Corner: x = nx, y = ny, z = nz ------------------------------
      ! The upper periodic corner receives the combined displacement jumps from
      ! all three coordinate directions.
      delta = 0.d0
      do i8 = 1, nnode_fe;
        if(i8==2) then
          delta(i8,1) = exx*nx
          delta(i8,2) = exy*nx
          delta(i8,3) = exz*nx
        end if
        if(i8==4) then
          delta(i8,1) = exy*ny
          delta(i8,2) = eyy*ny
          delta(i8,3) = eyz*ny
        end if
        if(i8==5) then
          delta(i8,1) = exz*nz
          delta(i8,2) = eyz*nz
          delta(i8,3) = ezz*nz
        end if
        if(i8==8) then
          delta(i8,1) = exy*ny + exz*nz
          delta(i8,2) = eyy*ny + eyz*nz
          delta(i8,3) = eyz*ny + ezz*nz
        end if
        if(i8==6) then
          delta(i8,1) = exx*nx + exz*nz
          delta(i8,2) = exy*nx + eyz*nz
          delta(i8,3) = exz*nx + ezz*nz
        end if
        if(i8==3) then
          delta(i8,1) = exx*nx + exy*ny
          delta(i8,2) = exy*nx + eyy*ny
          delta(i8,3) = exz*nx + eyz*ny
        end if
        if(i8==7) then
          delta(i8,1) = exx*nx + exy*ny + exz*nz
          delta(i8,2) = exy*nx + eyy*ny + eyz*nz
          delta(i8,3) = exz*nx + eyz*ny + ezz*nz
        end if
      end do;

      m = nx*ny*nz
      do nn = 1, ndof
        do mm = 1, nnode_fe
          sumval = 0.d0
          do m3 = 1, ndof
            do m8 = 1, nnode_fe
              sumval = sumval + delta(m8,m3)*dk(e3d%pix(m),m8,m3,mm,nn)
              C = C + 0.5d0 * delta(m8,m3)*dk(e3d%pix(m),m8,m3,mm,nn)*delta(mm,nn)
            end do
          end do
          b(e3d%ib(m,is(mm)),nn) = b(e3d%ib(m,is(mm)),nn) + sumval
        end do
      end do
	  ! -- End boundary and periodic terms section --

      ! Enforce zero equivalent force on inactive pore-only nodes.
      ! These nodes are excluded from the mechanical solve because they are not
      ! connected to any solid element.
      !$omp parallel do collapse(2) schedule(guided)
      do j = 1, ndof
        do m = 1, ns
          if (.not. e3d%is_active(m)) b(m,j) = 0.d0
        end do
      end do
      !$omp end parallel do

      ! ======================================================================
      ! ---- 9. Final Assignment and Cleanup ----------------------------------
      ! Store the assembled force vector and energy offset in the solver state,
      ! then release local work arrays.
      ! ======================================================================
      ! Copy assembled periodic-boundary terms back to the global solver state.	  
      e3d%b = b
      e3d%C = C

      ! Release local work arrays.
      deallocate(dk, b)

      end subroutine femat_OpenMP

!-------------------------------------------------------------------------   
      subroutine energy_OpenMP(e3d, utot)
      ! ==============================================================================
      ! SUBROUTINE: energy_OpenMP
      ! ------------------------------------------------------------------------------
      ! PURPOSE:
      !   Evaluate the total elastic energy for the current displacement field and
      !   compute the corresponding energy gradient/residual vector.
      !
      !   This routine performs the matrix-free operation:
      !
      !      gb = A*u
      !
      !   where A is the global stiffness operator implied by the local element
      !   stiffness tensors e3d%dk and the periodic 27-neighbor stencil. The global
      !   stiffness matrix A is not explicitly assembled.
      !
      !   After the matrix-free multiplication, the routine computes:
      !
      !      U(u) = C + 0.5*u^T*A*u + b^T*u
      !
      !   and finalizes the gradient:
      !
      !      grad(U) = A*u + b
      !
      !   The finalized gradient is stored in e3d%gb and used by the PCG solver.
      !
      ! ARGUMENTS:
      !   e3d
      !     Main solver data structure.
      !
      !     Input components:
      !       e3d%u
      !         Current displacement field, dimensioned ns x ndof.
      !
      !       e3d%b
      !         Equivalent force vector from imposed periodic macroscopic strain
      !         jumps.
      !
      !       e3d%C
      !         Constant energy offset from imposed periodic macroscopic strain
      !         jumps.
      !
      !       e3d%dk
      !         Local integrated element stiffness tensors for each phase.
      !
      !       e3d%pix
      !         Phase label field used to select the appropriate local stiffness.
      !
      !       e3d%ib
      !         27-point periodic neighbor table.
      !
      !       e3d%is_active
      !         Active-node mask. Inactive pore-only nodes are assigned zero
      !         gradient after assembly.
      !
      !     Output component:
      !       e3d%gb
      !         Final energy gradient/residual vector, grad(U) = A*u + b.
      !
      !   utot
      !     Total elastic energy of the current displacement field.
      !
      ! NUMERICAL FORMULATION:
      !   Energy:
      !      U = C + 0.5*u^T*A*u + b^T*u
      !
      !   Gradient:
      !      g = A*u + b
      !
      !   Matrix-free stiffness action:
      !      A*u is evaluated using the precomputed local element stiffness
      !      tensors e3d%dk and the periodic neighbor stencil e3d%ib.
      !
      ! PARALLELIZATION:
      !   OpenMP is used over grid sites and displacement components. Scalar energy
      !   accumulation uses reduction clauses.
      ! ==============================================================================	  
      use elas3d_mod
      implicit none

      ! ======================================================================
      ! ---- 1. Subroutine Arguments and Local Variables ----------------------
      ! ======================================================================

      ! Main solver state structure. Provides displacements, stiffness data,
      ! phase labels, periodic connectivity, force vector, and energy offset.
      type(elas3d_data_type), intent(inout) :: e3d
	  
      ! Total elastic energy for the current displacement field.	  
      real(dp), intent(out) :: utot

      ! Loop counters.
      integer :: m, m3, j
	  
      ! Local copy of the constant energy offset and temporary stencil sum.	  
      real(dp) :: C, gbsum
	  
      ! Local gradient/work array.
      ! Before adding e3d%b, gb stores A*u.
      ! After adding e3d%b, gb stores grad(U) = A*u + b.	  
      real(dp), allocatable :: gb(:,:)

      ! ======================================================================
      ! ---- 2. Memory Allocation and Initialization --------------------------
      ! Allocate the local work array gb for all grid sites and displacement
      ! components.
      !
      ! gb is first used to store the matrix-free stiffness product A*u.
      ! It is later converted to the full energy gradient A*u + b.
      ! ======================================================================
      allocate( gb(ns,ndof) )

      ! Local copy of the constant energy offset from periodic strain jumps.
      C   = e3d%C

      ! Zero the work array before accumulating stencil contributions.
      gb(:,:) = 0.d0

      ! ======================================================================
      ! ---- 3. Matrix-Free Stiffness Application: gb = A*u -------------------
      ! Compute the stiffness action A*u using the 27-point periodic stencil.
      !
      ! The outer loop over j selects the output displacement/gradient component:
      !   j = 1 : x-component of A*u
      !   j = 2 : y-component of A*u
      !   j = 3 : z-component of A*u
      !
      ! For each grid site m, the expanded stencil expression collects
      ! contributions from the x-, y-, and z-displacement components of the
      ! 27 periodic neighboring sites.
      !
      ! The long stencil expressions below are algebraically equivalent to
      ! multiplying the global stiffness operator A by the displacement vector u,
      ! but avoid explicitly storing A.
      ! ======================================================================	  
      do j = 1, ndof

        !$omp parallel do private(m, gbsum) schedule(guided)		
        do m = 1, ns
          gbsum = 0.0d0

          ! ==================================================================
          ! Contributions from x-displacement components, u(:,1)
          ! ------------------------------------------------------------------
          ! e3d%ib(m,idx) gives the periodic neighbor index.
          ! e3d%pix(...) selects the material phase associated with the
          ! neighboring element contribution.
          ! e3d%dk(...) provides the corresponding local stiffness entry.
          !
          ! This block contributes to output component j of A*u at grid site m.
          ! ==================================================================
          ! === input displacement component n = 1, ux ===
          gbsum = gbsum                                 &
        + e3d%u(e3d%ib(m,1),1 )*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,4,1) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,3,1) + e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,8,1) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,7,1) ) &
        + e3d%u(e3d%ib(m,2),1 )*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,3,1) + e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,7,1) )                                            &
        + e3d%u(e3d%ib(m,3),1 )*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,2,1) + e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,3,1) + e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,7,1) + e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,6,1) ) &
        + e3d%u(e3d%ib(m,4),1 )*( e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,2,1) + e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,6,1) )                                             &
        + e3d%u(e3d%ib(m,5),1 )*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,2,1) + e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,1,1) + e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,6,1) + e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,5,1) ) &
        + e3d%u(e3d%ib(m,6),1 )*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,1,1) + e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,5,1) )                                              &
        + e3d%u(e3d%ib(m,7),1 )*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,4,1) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,1,1) + e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,8,1) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,5,1) ) &
        + e3d%u(e3d%ib(m,8),1 )*( e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,4,1) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,8,1) )                                              &
        + e3d%u(e3d%ib(m,9),1 )*( e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,4,1) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,3,1) )                                             &
        + e3d%u(e3d%ib(m,10),1)*( e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,3,1) )                                               &
        + e3d%u(e3d%ib(m,11),1)*( e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,3,1) + e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,2,1) )                                             &
        + e3d%u(e3d%ib(m,12),1)*( e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,2,1) )                                               &
        + e3d%u(e3d%ib(m,13),1)*( e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,1,1) + e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,2,1) )                                             &
        + e3d%u(e3d%ib(m,14),1)*( e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,1,1) )                                               &
        + e3d%u(e3d%ib(m,15),1)*( e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,4,1) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,1,1) )                                             &
        + e3d%u(e3d%ib(m,16),1)*( e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,4,1) )                                               &
        + e3d%u(e3d%ib(m,17),1)*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,8,1) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,7,1) )                                              &
        + e3d%u(e3d%ib(m,18),1)*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,7,1) )                                               &
        + e3d%u(e3d%ib(m,19),1)*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,6,1) + e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,7,1) )                                              &
        + e3d%u(e3d%ib(m,20),1)*( e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,6,1) )                                                &
        + e3d%u(e3d%ib(m,21),1)*( e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,5,1) + e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,6,1) )                                               &
        + e3d%u(e3d%ib(m,22),1)*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,5,1) )                                                &
        + e3d%u(e3d%ib(m,23),1)*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,8,1) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,5,1) )                                              &
        + e3d%u(e3d%ib(m,24),1)*( e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,8,1) )                                                &
        + e3d%u(e3d%ib(m,25),1)*( e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,3,1) + e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,4,1) +  &
                                      e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,2,1) + e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,1,1) ) &
        + e3d%u(e3d%ib(m,26),1)*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,7,1) + e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,8,1) +  &
                                      e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,5,1) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,6,1) ) &
        + e3d%u(e3d%ib(m,27),1)*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,1,1) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,2,1) +  &
                                      e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,3,1) + e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,4,1) +  &
                                      e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,5,1) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,6,1) + &
                                      e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,7,1) + e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,8,1) )

          ! ==================================================================
          ! Contributions from y-displacement components, u(:,2)
          ! ------------------------------------------------------------------
          ! Adds the coupling between the y-displacement field and output
          ! component j of A*u at grid site m.
          ! ==================================================================
          ! === input displacement component n = 2, uy ===
          gbsum = gbsum                                 &
        + e3d%u(e3d%ib(m,1),2 )*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,4,2) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,3,2) + e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,8,2) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,7,2) ) &
        + e3d%u(e3d%ib(m,2),2 )*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,3,2) + e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,7,2) )                                            &
        + e3d%u(e3d%ib(m,3),2 )*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,2,2) + e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,3,2) + e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,7,2) + e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,6,2) ) &
        + e3d%u(e3d%ib(m,4),2 )*( e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,2,2) + e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,6,2) )                                             &
        + e3d%u(e3d%ib(m,5),2 )*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,2,2) + e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,1,2) + e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,6,2) + e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,5,2) ) &
        + e3d%u(e3d%ib(m,6),2 )*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,1,2) + e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,5,2) )                                              &
        + e3d%u(e3d%ib(m,7),2 )*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,4,2) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,1,2) + e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,8,2) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,5,2) ) &
        + e3d%u(e3d%ib(m,8),2 )*( e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,4,2) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,8,2) )                                              &
        + e3d%u(e3d%ib(m,9),2 )*( e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,4,2) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,3,2) )                                             &
        + e3d%u(e3d%ib(m,10),2)*( e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,3,2) )                                               &
        + e3d%u(e3d%ib(m,11),2)*( e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,3,2) + e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,2,2) )                                             &
        + e3d%u(e3d%ib(m,12),2)*( e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,2,2) )                                               &
        + e3d%u(e3d%ib(m,13),2)*( e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,1,2) + e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,2,2) )                                             &
        + e3d%u(e3d%ib(m,14),2)*( e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,1,2) )                                               &
        + e3d%u(e3d%ib(m,15),2)*( e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,4,2) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,1,2) )                                             &
        + e3d%u(e3d%ib(m,16),2)*( e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,4,2) )                                               &
        + e3d%u(e3d%ib(m,17),2)*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,8,2) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,7,2) )                                              &
        + e3d%u(e3d%ib(m,18),2)*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,7,2) )                                               &
        + e3d%u(e3d%ib(m,19),2)*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,6,2) + e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,7,2) )                                              &
        + e3d%u(e3d%ib(m,20),2)*( e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,6,2) )                                                &
        + e3d%u(e3d%ib(m,21),2)*( e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,5,2) + e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,6,2) )                                               &
        + e3d%u(e3d%ib(m,22),2)*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,5,2) )                                                &
        + e3d%u(e3d%ib(m,23),2)*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,8,2) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,5,2) )                                              &
        + e3d%u(e3d%ib(m,24),2)*( e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,8,2) )                                                &
        + e3d%u(e3d%ib(m,25),2)*( e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,3,2) + e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,4,2) +  &
                                      e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,2,2) + e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,1,2) ) &
        + e3d%u(e3d%ib(m,26),2)*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,7,2) + e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,8,2) +  &
                                      e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,5,2) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,6,2) ) &
        + e3d%u(e3d%ib(m,27),2)*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,1,2) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,2,2) +  &
                                      e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,3,2) + e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,4,2) +  &
                                      e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,5,2) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,6,2) + &
                                      e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,7,2) + e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,8,2) )

          ! ==================================================================
          ! Contributions from z-displacement components, u(:,3)
          ! ------------------------------------------------------------------
          ! Adds the coupling between the z-displacement field and output
          ! component j of A*u at grid site m.
          ! ==================================================================
          ! === input displacement component n = 3, uz ===
          gbsum = gbsum                                 &
        + e3d%u(e3d%ib(m,1),3 )*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,4,3) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,3,3) + e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,8,3) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,7,3) ) &
        + e3d%u(e3d%ib(m,2),3 )*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,3,3) + e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,7,3) )                                            &
        + e3d%u(e3d%ib(m,3),3 )*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,2,3) + e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,3,3) + e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,7,3) + e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,6,3) ) &
        + e3d%u(e3d%ib(m,4),3 )*( e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,2,3) + e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,6,3) )                                             &
        + e3d%u(e3d%ib(m,5),3 )*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,2,3) + e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,1,3) + e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,6,3) + e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,5,3) ) &
        + e3d%u(e3d%ib(m,6),3 )*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,1,3) + e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,5,3) )                                              &
        + e3d%u(e3d%ib(m,7),3 )*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,4,3) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,1,3) + e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,8,3) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,5,3) ) &
        + e3d%u(e3d%ib(m,8),3 )*( e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,4,3) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,8,3) )                                              &
        + e3d%u(e3d%ib(m,9),3 )*( e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,4,3) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,3,3) )                                             &
        + e3d%u(e3d%ib(m,10),3)*( e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,3,3) )                                               &
        + e3d%u(e3d%ib(m,11),3)*( e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,3,3) + e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,2,3) )                                             &
        + e3d%u(e3d%ib(m,12),3)*( e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,2,3) )                                               &
        + e3d%u(e3d%ib(m,13),3)*( e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,1,3) + e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,2,3) )                                             &
        + e3d%u(e3d%ib(m,14),3)*( e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,1,3) )                                               &
        + e3d%u(e3d%ib(m,15),3)*( e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,4,3) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,1,3) )                                             &
        + e3d%u(e3d%ib(m,16),3)*( e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,4,3) )                                               &
        + e3d%u(e3d%ib(m,17),3)*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,8,3) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,7,3) )                                              &
        + e3d%u(e3d%ib(m,18),3)*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,7,3) )                                               &
        + e3d%u(e3d%ib(m,19),3)*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,6,3) + e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,7,3) )                                              &
        + e3d%u(e3d%ib(m,20),3)*( e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,6,3) )                                                &
        + e3d%u(e3d%ib(m,21),3)*( e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,5,3) + e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,6,3) )                                               &
        + e3d%u(e3d%ib(m,22),3)*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,5,3) )                                                &
        + e3d%u(e3d%ib(m,23),3)*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,8,3) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,5,3) )                                              &
        + e3d%u(e3d%ib(m,24),3)*( e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,8,3) )                                                &
        + e3d%u(e3d%ib(m,25),3)*( e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,3,3) + e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,4,3) +  &
                                      e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,2,3) + e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,1,3) ) &
        + e3d%u(e3d%ib(m,26),3)*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,7,3) + e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,8,3) +  &
                                      e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,5,3) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,6,3) ) &
        + e3d%u(e3d%ib(m,27),3)*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,1,3) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,2,3) +  &
                                      e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,3,3) + e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,4,3) +  &
                                      e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,5,3) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,6,3) + &
                                      e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,7,3) + e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,8,3) )

          ! Store output component j of the matrix-free product A*u.
          gb(m, j) = gbsum
        end do
        !$omp end parallel do
      end do

      ! ======================================================================
      ! ---- 4. Total Energy and Gradient Finalization ------------------------
      ! At this point, gb contains the matrix-free stiffness product A*u.
      !
      ! Compute:
      !   U = C + 0.5*u^T*A*u + b^T*u
      !
      ! Then convert:
      !   gb = A*u
      !
      ! into the full gradient:
      !   gb = A*u + b
      !
      ! The resulting e3d%gb is used by the PCG solver as the residual/gradient
      ! direction.
      ! ======================================================================
      ! Start energy accumulation with the constant offset from periodic jumps.	  
      utot = C
      do m3 = 1, ndof
        !$omp parallel do schedule(guided) reduction(+:utot)		
        do m = 1, ns
          ! Add this DOF contribution to the total energy and finalize the
          ! gradient by adding the linear force term b.		
          utot = utot + 0.5d0 * e3d%u(m, m3) * gb(m, m3) + e3d%b(m, m3) * e3d%u(m, m3)
          gb(m, m3) = gb(m, m3) + e3d%b(m, m3)
        end do
        !$omp end parallel do
      end do
	  
      ! ======================================================================
      ! ---- 5. Store Gradient and Enforce Inactive-Node Constraints ----------
      ! Copy the finalized gradient into the global solver state. Then enforce
      ! zero gradient on inactive pore-only nodes so they do not participate in
      ! the PCG update.
      ! ======================================================================
      e3d%gb = gb
	  
      ! Inactive nodes are connected only to pore/void elements. Set their
      ! gradient entries to zero to remove them from the mechanical solve.
      !$omp parallel do collapse(2) schedule(guided)
      do j = 1, ndof
        do m = 1, ns
          if (.not. e3d%is_active(m)) e3d%gb(m,j) = 0.d0
        end do
      end do
      !$omp end parallel do 

      ! Release local work array.
      deallocate( gb )

      end subroutine energy_OpenMP
  
!-------------------------------------------------------------------------
      subroutine dembx_OpenMP(e3d, Lstep, gg, gg0, kkk)
      ! ==============================================================================
      ! SUBROUTINE: dembx_OpenMP
      ! ------------------------------------------------------------------------------
      ! PURPOSE:
      !   Solves for the equilibrium displacement field by minimizing the total
      !   elastic energy using a Preconditioned Conjugate Gradient (PCG) algorithm.
      !
      !   The solver operates on the energy gradient/residual stored in e3d%gb.
      !   The global stiffness matrix is not assembled explicitly. Instead, the
      !   matrix-vector product A*h is evaluated using the same 27-point periodic
      !   stencil structure used by energy_OpenMP.
      !
      ! MATHEMATICAL FORM:
      !   The energy functional has the form:
      !
      !      U(u) = 0.5*u^T*A*u + b^T*u + C
      !
      !   with gradient:
      !
      !      g = A*u + b
      !
      !   This routine iteratively drives g toward zero.
      !
      !   The preconditioned residual/gradient is:
      !
      !      z = M^{-1} g
      !
      !   and the scalar convergence measure used internally is:
      !
      !      gg = z^T g
      !
      !   The relative residual reported by the code is:
      !
      !      rel_res = sqrt(gg/gg0)
      !
      ! ARGUMENTS:
      !   e3d
      !     Main solver data structure.
      !
      !     Input components:
      !       e3d%u
      !         Current displacement field.
      !
      !       e3d%gb
      !         Current energy gradient/residual, usually computed by energy_OpenMP.
      !
      !       e3d%h
      !         Previous or initial conjugate search direction.
      !
      !       e3d%dk
      !         Phase-wise local finite-element stiffness tensors.
      !
      !       e3d%pix
      !         Phase-label field.
      !
      !       e3d%ib
      !         27-point periodic neighbor table.
      !
      !       e3d%is_active
      !         Active-node mask. Inactive pore-only nodes are constrained to zero
      !         displacement and zero residual.
      !
      !     Output/updated components:
      !       e3d%u
      !         Updated displacement field after PCG iterations.
      !
      !       e3d%gb
      !         Updated gradient/residual.
      !
      !       e3d%h
      !         Updated conjugate search direction.
      !
      !   Lstep
      !     Number of PCG iterations completed in this call.
      !
      !   gg
      !     Final preconditioned residual inner product, z^T g. The relative
      !     residual norm is computed as sqrt(gg/gg0).
      !
      !   gg0
      !     Initial reference value for the residual measure. Used to compute
      !     the relative residual norm.
      !
      !   kkk
      !     Outer solver/restart index. In the current configuration, kmax = 1,
      !     so kkk is normally 1. The block-Jacobi preconditioner is assembled
      !     when kkk == 1.
      !
      ! PRECONDITIONER:
      !   A 3 x 3 block-Jacobi preconditioner is built for each active grid node.
      !   The block is formed from the local diagonal stiffness contributions of
      !   the eight hexahedral elements sharing that node. Each 3 x 3 block is
      !   inverted using invert3x3.
      !
      !   Note:
      !     M_block_inv is local to this subroutine and is assembled only when
      !     kkk == 1. This is consistent with the current module setting kmax = 1.
      !     If multiple outer calls are enabled, the preconditioner should either
      !     be rebuilt every call or stored persistently.
      !
      ! MATRIX-FREE OPERATOR:
      !   The product Ah = A*h is evaluated using an expanded 27-point stencil.
      !   Each grid site receives contributions from x-, y-, and z-components of
      !   the neighboring displacement/search-direction field.
      !
      ! CONVERGENCE AND STABILITY:
      !   The loop exits when:
      !
      !      sqrt(gg/gg0) < tol
      !
      !   A residual-spike safeguard resets the search direction to the
      !   preconditioned gradient direction if the residual increases sharply.
      !
      ! OUTPUT:
      !   Writes detailed local PCG residual history to:
      !
      !      cgitr_SS316L_polycrystal.txt
      !
      ! PARALLELIZATION:
      !   OpenMP is used for vector operations, preconditioner construction,
      !   matrix-free stencil evaluation, and reductions.
      ! ==============================================================================
      use elas3d_mod
      implicit none

      ! ======================================================================
      ! ---- 1. Subroutine Arguments -----------------------------------------
      ! PCG solver for the equilibrium displacement field.
      !
      ! The residual/gradient convention is:
      !   gb = A*u + b
      !
      ! The preconditioned gradient is:
      !   z = M^{-1}*gb
      !
      ! The code updates:
      !   u  <- u  - lambda*h
      !   gb <- gb - lambda*A*h
      !
      ! where h is the current preconditioned conjugate search direction.
      ! ======================================================================
      type(elas3d_data_type), intent(inout) :: e3d
	  
      ! Number of PCG iterations performed in this call.	  
      integer, intent(out) :: Lstep
	  
      ! Final preconditioned residual inner product, gg = z^T*gb.
      ! The relative residual norm is sqrt(gg/gg0).
      real(dp), intent(out) :: gg
	  
      ! Reference residual measure from the initial state.
      real(dp), intent(in)  :: gg0
	  
      ! Outer solver/restart index. The current code expects kkk == 1 because
      ! the preconditioner is local and assembled only for kkk == 1.
      integer, intent(in)   :: kkk

      ! ======================================================================
      ! ---- 2. Local Array Declarations -------------------------------------
      ! Local work arrays used during PCG iterations.
      !
      ! u, gb, h:
      !   Local copies of displacement, gradient/residual, and conjugate search
      !   direction. These are copied back to e3d at the end of the solve.
      !
      ! Ah:
      !   Matrix-free product A*h.
      !
      ! M_block_inv:
      !   Inverse 3 x 3 block-Jacobi preconditioner for each grid site.
      !
      ! z:
      !   Preconditioned gradient/residual, z = M^{-1}*gb.
      !
      ! history_step, history_res:
      !   Local convergence-history arrays written to cgitr_SS316L_polycrystal.txt.
      ! ======================================================================
      real(dp), allocatable :: u(:,:), gb(:,:), h(:,:)
      real(dp), allocatable :: Ah(:,:)
      real(dp), allocatable :: M_block_inv(:,:,:) ! Precomputed inverse
      real(dp), allocatable :: z(:,:)      ! Preconditioned residual  
      integer, allocatable :: history_step(:)
      real(dp), allocatable :: history_res(:)
    
      ! ======================================================================
      ! ---- 3. Local Scalars and PCG Variables -------------------------------
      ! lambda:
      !   PCG step length, equivalent to alpha in standard notation.
      !
      ! gamma:
      !   Conjugate-direction update factor, equivalent to beta.
      !
      ! hAh:
      !   Curvature along the current search direction, h^T*A*h.
      !
      ! gglast:
      !   Previous preconditioned residual inner product.
      !
      ! rel_res:
      !   Relative residual norm, sqrt(gg/gg0).
      !
      ! reset_count:
      !   Number of times the search direction has been reset due to a residual
      !   spike.
      ! ======================================================================
      integer :: ns_local, L1, L2
      integer :: m, ijk, i, j, k
      integer :: ios	  
      real(dp) :: lambda, gamma, hAh, gglast, rel_res
      real(dp) :: rel_res_prev, rel_res_min
      integer :: reset_count	 

      ! Residual-spike safeguard parameters.
      ! If rel_res increases by RESET_RATIO and is still above RESET_ABS_TOL,
      ! the search direction is reset to the preconditioned gradient.
      real(dp), parameter :: RESET_RATIO    = 100.0d0
      real(dp), parameter :: RESET_ABS_TOL = 1.0d-4
	  
      ! ======================================================================
      ! ---- 4. Thread-Private Stencil Work Variables -------------------------
      ! ibl:
      !   Local copy of the 27 periodic neighbor indices for one grid site.
      !
      ! pl:
      !   Phase labels associated with those 27 neighboring sites.
      !
      ! ahsum1, ahsum2, ahsum3:
      !   Accumulators for the three components of A*h at one grid site.
      !
      ! M_block:
      !   Local 3 x 3 diagonal stiffness block used to build the block-Jacobi
      !   preconditioner.
      ! ======================================================================
      integer :: ibl(27), pl(27)
      real(dp) :: ahsum1, ahsum2, ahsum3 
	  real(dp) :: M_block(3,3)

      ! Open detailed PCG iteration-history file.	  
	  open(unit=100, file='cgitr_SS316L_polycrystal.txt', status='replace', iostat=ios)

      ns_local = ns

      ! ======================================================================
      ! ---- 5. Memory Allocation --------------------------------------------
      ! Allocate local PCG work arrays using the global number of grid sites and
      ! displacement components.
      !
      ! Note:
      !   M_block_inv is allocated locally and must be assembled before it is used.
      !   With the current setting kmax = 1, this occurs when kkk == 1.
      ! ======================================================================
      allocate( u(ns,ndof), gb(ns,ndof), h(ns,ndof) )
      allocate( Ah(ns,ndof) )
	  allocate( M_block_inv(ns,ndof,ndof) )
      allocate( z(ns,ndof) )
      allocate(history_step(ldemb))
      allocate(history_res(ldemb))	  

      ! ======================================================================
      ! ---- 6. Copy Solver State to Local Work Arrays ------------------------
      ! Copy displacement, gradient/residual, and search direction from e3d into
      ! local arrays. Local arrays improve readability and can reduce repeated
      ! derived-type component access inside tight loops.
      ! ======================================================================
	  !$omp parallel do schedule(guided)
      do m = 1, ns
        u(m,:)  = e3d%u(m,:)
        gb(m,:) = e3d%gb(m,:)
        h(m,:)  = e3d%h(m,:)
      end do
      !$omp end parallel do
	  
      ! Enforce zero displacement, gradient, and search direction on inactive
      ! pore-only nodes before starting PCG iterations.
      !$omp parallel do collapse(2) schedule(guided)
      do j = 1, ndof
        do m = 1, ns
          if (.not. e3d%is_active(m)) then
            u(m,j)  = 0.d0
            gb(m,j) = 0.d0
            h(m,j)  = 0.d0
          end if
        end do
      end do
      !$omp end parallel do	  
	  
  
      ! ======================================================================
      ! ---- 7. Block-Jacobi Preconditioner Assembly --------------------------
      ! Build and invert a 3 x 3 block-Jacobi preconditioner for each grid site.
      !
      ! For each node m, the local diagonal block is assembled by summing the
      ! stiffness contributions from the eight hexahedral elements sharing that
      ! node. The relevant element phase labels are obtained from the 27-point
      ! periodic neighbor table.
      !
      ! The inverse block M_block_inv(m,:,:) is later applied to the gradient:
      !
      !   z(m,:) = M_block_inv(m,:,:) * gb(m,:)
      !
      ! This block is currently assembled only when kkk == 1, consistent with
      ! kmax = 1 in the module.
      ! ======================================================================
	  if (kkk == 1) then
	    ! Build Preconditioner 
        !$omp parallel do private(m,ibl,pl,i,j,M_block) schedule(guided)
        do m = 1, ns
          ! Cache periodic neighbor indices and their phase labels.		
          ibl = e3d%ib(m,1:27)
          pl  = e3d%pix(ibl(1:27))
		  
          ! Assemble the full 3 x 3 diagonal stiffness block for node m.
          do i = 1, 3
            do j = 1, 3
              M_block(i,j) = e3d%dk(pl(27),1,i,1,j) + e3d%dk(pl(7),2,i,2,j) +  &
                             e3d%dk(pl(6),3,i,3,j)  + e3d%dk(pl(5),4,i,4,j) +  &
                             e3d%dk(pl(25),5,i,5,j) + e3d%dk(pl(15),6,i,6,j) + &
                             e3d%dk(pl(14),7,i,7,j) + e3d%dk(pl(13),8,i,8,j)
            end do
          end do
          ! Invert the local block and store it for repeated use during PCG.
          ! If M_block is singular, invert3x3 returns a zero matrix.
          ! This is acceptable for inactive nodes because their residuals are
          ! constrained to zero, but active singular blocks should be checked
          ! carefully if convergence problems occur.		  
          call invert3x3(M_block, M_block_inv(m,:,:))
        end do
        !$omp end parallel do
      end if
  
      ! ======================================================================
      ! ---- 8. PCG Initialization -------------------------------------------
      ! Apply the block-Jacobi preconditioner to the initial gradient:
      !
      !   z = M^{-1}*gb
      !
      ! Then compute the preconditioned residual inner product:
      !
      !   gg = z^T*gb
      !
      ! For the first outer call, initialize the search direction h with z.
      ! Since the displacement update uses u <- u - lambda*h, this corresponds
      ! to moving in the negative preconditioned-gradient direction.
      ! ======================================================================
      !$omp parallel do collapse(2) schedule(guided)
      do i = 1, ndof 
	    do m = 1, ns 
        z(m,i) = M_block_inv(m,i,1) * gb(m,1) + &
		         M_block_inv(m,i,2) * gb(m,2) + &
				 M_block_inv(m,i,3) * gb(m,3)
        end do
	  end do
      !$omp end parallel do
  
  
      ! ======================================================================
      ! Compute the preconditioned residual/gradient inner product gg = z^T*gb.
      ! ======================================================================
      gg = 0.d0 
      !$omp parallel do collapse(2) schedule(guided) reduction(+:gg)	  
      do j = 1, ndof 
	    do i = 1, ns
          gg = gg + z(i,j)*gb(i,j)
        end do
	  end do
      !$omp end parallel do    

      ! ======================================================================
      ! Initialize the search direction with the preconditioned gradient.
      ! The descent sign is applied in the update u = u - lambda*h.
      ! ======================================================================
      if (kkk == 1) then
        !$omp parallel do collapse(2) schedule(guided)
        do i = 1, ndof 
		  do m = 1, ns
            h(m,i) = z(m,i)
          end do
        end do
        !$omp end parallel do
      end if
	  
      ! Initialize residual-spike tracking variables.
      rel_res_prev = sqrt(gg/gg0)
      rel_res_min = rel_res_prev
      reset_count = 0	  


      ! ======================================================================
      ! ---- 9. Main Preconditioned Conjugate-Gradient Loop -------------------
      ! Each iteration performs:
      !   1. Matrix-free product Ah = A*h.
      !   2. Step-size calculation lambda = (z^T*g)/(h^T*A*h).
      !   3. Displacement and gradient update.
      !   4. Preconditioner application z = M^{-1}*g.
      !   5. Conjugate-direction update or residual-spike reset.
      !   6. Convergence check using sqrt(gg/gg0).
      ! ======================================================================
      Lstep = 0
      do ijk = 1, ldemb
        Lstep = Lstep + 1
        ! Clear matrix-vector product array before accumulating Ah = A*h.	
        Ah = 0.d0

        ! ====================================================================
        ! ---- 9a. Matrix-Free Product: Ah = A*h ------------------------------
        ! Evaluate the action of the global stiffness operator A on the current
        ! search direction h using the expanded 27-point periodic stencil.
        !
        ! The calculation is blocked to improve cache locality and reduce memory
        ! bandwidth pressure for large 3D grids.
        !
        ! For each grid site m:
        !   Ah(m,1) receives x-component contributions,
        !   Ah(m,2) receives y-component contributions,
        !   Ah(m,3) receives z-component contributions.
        ! ====================================================================
        do L1 = 1, ns_local, block_size
          L2 = min(L1 + block_size - 1, ns_local)

          ! Matrix multiply: Ah = A * h
          !$omp parallel do private(m, ibl, pl, ahsum1, ahsum2, ahsum3) schedule(guided)
          do m = L1, L2
            ! Cache the 27 periodic neighbor indices and their phase labels
            ! for the current grid site.
            ibl = e3d%ib(m,1:27)      ! All 27 neighbors
            pl  = e3d%pix(ibl(1:27))  ! Phase index for all neighbors

            ! Output DOF 1: x-component of Ah = A*h.
            ahsum1 = 0.d0
            ! Contributions from input displacement/search component n = 1 (x).
            ahsum1 = ahsum1 &
          + h(ibl(1),1 )*( e3d%dk(pl(27),1,1,4,1) + e3d%dk(pl(7),2,1,3,1) + e3d%dk(pl(25),5,1,8,1) + e3d%dk(pl(15),6,1,7,1) ) &
          + h(ibl(2),1)*( e3d%dk(pl(27),1,1,3,1) + e3d%dk(pl(25),5,1,7,1) ) &
          + h(ibl(3),1)*( e3d%dk(pl(27),1,1,2,1) + e3d%dk(pl(5),4,1,3,1) + e3d%dk(pl(13),8,1,7,1) + e3d%dk(pl(25),5,1,6,1) ) &
          + h(ibl(4),1)*( e3d%dk(pl(5),4,1,2,1) + e3d%dk(pl(13),8,1,6,1) ) &
          + h(ibl(5),1)*( e3d%dk(pl(6),3,1,2,1) + e3d%dk(pl(5),4,1,1,1) + e3d%dk(pl(14),7,1,6,1) + e3d%dk(pl(13),8,1,5,1) ) &
          + h(ibl(6),1)*( e3d%dk(pl(6),3,1,1,1) + e3d%dk(pl(14),7,1,5,1) ) &
          + h(ibl(7),1)*( e3d%dk(pl(6),3,1,4,1) + e3d%dk(pl(7),2,1,1,1) + e3d%dk(pl(14),7,1,8,1) + e3d%dk(pl(15),6,1,5,1) ) &
          + h(ibl(8),1)*( e3d%dk(pl(7),2,1,4,1) + e3d%dk(pl(15),6,1,8,1) ) &
          + h(ibl(9),1)*( e3d%dk(pl(25),5,1,4,1) + e3d%dk(pl(15),6,1,3,1) ) &
          + h(ibl(10),1)*( e3d%dk(pl(25),5,1,3,1) ) &
          + h(ibl(11),1)*( e3d%dk(pl(13),8,1,3,1) + e3d%dk(pl(25),5,1,2,1) ) &
          + h(ibl(12),1)*( e3d%dk(pl(13),8,1,2,1) ) &
          + h(ibl(13),1)*( e3d%dk(pl(13),8,1,1,1) + e3d%dk(pl(14),7,1,2,1) ) &
          + h(ibl(14),1)*( e3d%dk(pl(14),7,1,1,1) ) &
          + h(ibl(15),1)*( e3d%dk(pl(14),7,1,4,1) + e3d%dk(pl(15),6,1,1,1) ) &
          + h(ibl(16),1)*( e3d%dk(pl(15),6,1,4,1) ) &
          + h(ibl(17),1)*( e3d%dk(pl(27),1,1,8,1) + e3d%dk(pl(7),2,1,7,1) ) &
          + h(ibl(18),1)*( e3d%dk(pl(27),1,1,7,1) ) &
          + h(ibl(19),1)*( e3d%dk(pl(27),1,1,6,1) + e3d%dk(pl(5),4,1,7,1) ) &
          + h(ibl(20),1)*( e3d%dk(pl(5),4,1,6,1) ) &
          + h(ibl(21),1)*( e3d%dk(pl(5),4,1,5,1) + e3d%dk(pl(6),3,1,6,1) ) &
          + h(ibl(22),1)*( e3d%dk(pl(6),3,1,5,1) ) &
          + h(ibl(23),1)*( e3d%dk(pl(6),3,1,8,1) + e3d%dk(pl(7),2,1,5,1) ) &
          + h(ibl(24),1)*( e3d%dk(pl(7),2,1,8,1) ) &
          + h(ibl(25),1)*( e3d%dk(pl(14),7,1,3,1) + e3d%dk(pl(13),8,1,4,1) + e3d%dk(pl(15),6,1,2,1) + e3d%dk(pl(25),5,1,1,1) ) &
          + h(ibl(26),1)*( e3d%dk(pl(6),3,1,7,1) + e3d%dk(pl(5),4,1,8,1) + e3d%dk(pl(27),1,1,5,1) + e3d%dk(pl(7),2,1,6,1) ) &
          + h(ibl(27),1)*( e3d%dk(pl(27),1,1,1,1) + e3d%dk(pl(7),2,1,2,1) + e3d%dk(pl(6),3,1,3,1) + e3d%dk(pl(5),4,1,4,1) &
                          + e3d%dk(pl(25),5,1,5,1) + e3d%dk(pl(15),6,1,6,1) + e3d%dk(pl(14),7,1,7,1) + e3d%dk(pl(13),8,1,8,1) )
            ! Contributions from input displacement/search component n = 2 (y).
            ahsum1 = ahsum1 &
          + h(ibl(1),2 )*( e3d%dk(pl(27),1,1,4,2) + e3d%dk(pl(7),2,1,3,2) + e3d%dk(pl(25),5,1,8,2) + e3d%dk(pl(15),6,1,7,2) ) &
          + h(ibl(2),2)*( e3d%dk(pl(27),1,1,3,2) + e3d%dk(pl(25),5,1,7,2) ) &
          + h(ibl(3),2)*( e3d%dk(pl(27),1,1,2,2) + e3d%dk(pl(5),4,1,3,2) + e3d%dk(pl(13),8,1,7,2) + e3d%dk(pl(25),5,1,6,2) ) &
          + h(ibl(4),2)*( e3d%dk(pl(5),4,1,2,2) + e3d%dk(pl(13),8,1,6,2) ) &
          + h(ibl(5),2)*( e3d%dk(pl(6),3,1,2,2) + e3d%dk(pl(5),4,1,1,2) + e3d%dk(pl(14),7,1,6,2) + e3d%dk(pl(13),8,1,5,2) ) &
          + h(ibl(6),2)*( e3d%dk(pl(6),3,1,1,2) + e3d%dk(pl(14),7,1,5,2) ) &
          + h(ibl(7),2)*( e3d%dk(pl(6),3,1,4,2) + e3d%dk(pl(7),2,1,1,2) + e3d%dk(pl(14),7,1,8,2) + e3d%dk(pl(15),6,1,5,2) ) &
          + h(ibl(8),2)*( e3d%dk(pl(7),2,1,4,2) + e3d%dk(pl(15),6,1,8,2) ) &
          + h(ibl(9),2)*( e3d%dk(pl(25),5,1,4,2) + e3d%dk(pl(15),6,1,3,2) ) &
          + h(ibl(10),2)*( e3d%dk(pl(25),5,1,3,2) ) &
          + h(ibl(11),2)*( e3d%dk(pl(13),8,1,3,2) + e3d%dk(pl(25),5,1,2,2) ) &
          + h(ibl(12),2)*( e3d%dk(pl(13),8,1,2,2) ) &
          + h(ibl(13),2)*( e3d%dk(pl(13),8,1,1,2) + e3d%dk(pl(14),7,1,2,2) ) &
          + h(ibl(14),2)*( e3d%dk(pl(14),7,1,1,2) ) &
          + h(ibl(15),2)*( e3d%dk(pl(14),7,1,4,2) + e3d%dk(pl(15),6,1,1,2) ) &
          + h(ibl(16),2)*( e3d%dk(pl(15),6,1,4,2) ) &
          + h(ibl(17),2)*( e3d%dk(pl(27),1,1,8,2) + e3d%dk(pl(7),2,1,7,2) ) &
          + h(ibl(18),2)*( e3d%dk(pl(27),1,1,7,2) ) &
          + h(ibl(19),2)*( e3d%dk(pl(27),1,1,6,2) + e3d%dk(pl(5),4,1,7,2) ) &
          + h(ibl(20),2)*( e3d%dk(pl(5),4,1,6,2) ) &
          + h(ibl(21),2)*( e3d%dk(pl(5),4,1,5,2) + e3d%dk(pl(6),3,1,6,2) ) &
          + h(ibl(22),2)*( e3d%dk(pl(6),3,1,5,2) ) &
          + h(ibl(23),2)*( e3d%dk(pl(6),3,1,8,2) + e3d%dk(pl(7),2,1,5,2) ) &
          + h(ibl(24),2)*( e3d%dk(pl(7),2,1,8,2) ) &
          + h(ibl(25),2)*( e3d%dk(pl(14),7,1,3,2) + e3d%dk(pl(13),8,1,4,2) + e3d%dk(pl(15),6,1,2,2) + e3d%dk(pl(25),5,1,1,2) ) &
          + h(ibl(26),2)*( e3d%dk(pl(6),3,1,7,2) + e3d%dk(pl(5),4,1,8,2) + e3d%dk(pl(27),1,1,5,2) + e3d%dk(pl(7),2,1,6,2) ) &
          + h(ibl(27),2)*( e3d%dk(pl(27),1,1,1,2) + e3d%dk(pl(7),2,1,2,2) + e3d%dk(pl(6),3,1,3,2) + e3d%dk(pl(5),4,1,4,2) &
                          + e3d%dk(pl(25),5,1,5,2) + e3d%dk(pl(15),6,1,6,2) + e3d%dk(pl(14),7,1,7,2) + e3d%dk(pl(13),8,1,8,2) )
            ! Contributions from input displacement/search component n = 3 (z).
            ahsum1 = ahsum1 &
          + h(ibl(1),3 )*( e3d%dk(pl(27),1,1,4,3) + e3d%dk(pl(7),2,1,3,3) + e3d%dk(pl(25),5,1,8,3) + e3d%dk(pl(15),6,1,7,3) ) &
          + h(ibl(2),3)*( e3d%dk(pl(27),1,1,3,3) + e3d%dk(pl(25),5,1,7,3) ) &
          + h(ibl(3),3)*( e3d%dk(pl(27),1,1,2,3) + e3d%dk(pl(5),4,1,3,3) + e3d%dk(pl(13),8,1,7,3) + e3d%dk(pl(25),5,1,6,3) ) &
          + h(ibl(4),3)*( e3d%dk(pl(5),4,1,2,3) + e3d%dk(pl(13),8,1,6,3) ) &
          + h(ibl(5),3)*( e3d%dk(pl(6),3,1,2,3) + e3d%dk(pl(5),4,1,1,3) + e3d%dk(pl(14),7,1,6,3) + e3d%dk(pl(13),8,1,5,3) ) &
          + h(ibl(6),3)*( e3d%dk(pl(6),3,1,1,3) + e3d%dk(pl(14),7,1,5,3) ) &
          + h(ibl(7),3)*( e3d%dk(pl(6),3,1,4,3) + e3d%dk(pl(7),2,1,1,3) + e3d%dk(pl(14),7,1,8,3) + e3d%dk(pl(15),6,1,5,3) ) &
          + h(ibl(8),3)*( e3d%dk(pl(7),2,1,4,3) + e3d%dk(pl(15),6,1,8,3) ) &
          + h(ibl(9),3)*( e3d%dk(pl(25),5,1,4,3) + e3d%dk(pl(15),6,1,3,3) ) &
          + h(ibl(10),3)*( e3d%dk(pl(25),5,1,3,3) ) &
          + h(ibl(11),3)*( e3d%dk(pl(13),8,1,3,3) + e3d%dk(pl(25),5,1,2,3) ) &
          + h(ibl(12),3)*( e3d%dk(pl(13),8,1,2,3) ) &
          + h(ibl(13),3)*( e3d%dk(pl(13),8,1,1,3) + e3d%dk(pl(14),7,1,2,3) ) &
          + h(ibl(14),3)*( e3d%dk(pl(14),7,1,1,3) ) &
          + h(ibl(15),3)*( e3d%dk(pl(14),7,1,4,3) + e3d%dk(pl(15),6,1,1,3) ) &
          + h(ibl(16),3)*( e3d%dk(pl(15),6,1,4,3) ) &
          + h(ibl(17),3)*( e3d%dk(pl(27),1,1,8,3) + e3d%dk(pl(7),2,1,7,3) ) &
          + h(ibl(18),3)*( e3d%dk(pl(27),1,1,7,3) ) &
          + h(ibl(19),3)*( e3d%dk(pl(27),1,1,6,3) + e3d%dk(pl(5),4,1,7,3) ) &
          + h(ibl(20),3)*( e3d%dk(pl(5),4,1,6,3) ) &
          + h(ibl(21),3)*( e3d%dk(pl(5),4,1,5,3) + e3d%dk(pl(6),3,1,6,3) ) &
          + h(ibl(22),3)*( e3d%dk(pl(6),3,1,5,3) ) &
          + h(ibl(23),3)*( e3d%dk(pl(6),3,1,8,3) + e3d%dk(pl(7),2,1,5,3) ) &
          + h(ibl(24),3)*( e3d%dk(pl(7),2,1,8,3) ) &
          + h(ibl(25),3)*( e3d%dk(pl(14),7,1,3,3) + e3d%dk(pl(13),8,1,4,3) + e3d%dk(pl(15),6,1,2,3) + e3d%dk(pl(25),5,1,1,3) ) &
          + h(ibl(26),3)*( e3d%dk(pl(6),3,1,7,3) + e3d%dk(pl(5),4,1,8,3) + e3d%dk(pl(27),1,1,5,3) + e3d%dk(pl(7),2,1,6,3) ) &
          + h(ibl(27),3)*( e3d%dk(pl(27),1,1,1,3) + e3d%dk(pl(7),2,1,2,3) + e3d%dk(pl(6),3,1,3,3) + e3d%dk(pl(5),4,1,4,3) &
                          + e3d%dk(pl(25),5,1,5,3) + e3d%dk(pl(15),6,1,6,3) + e3d%dk(pl(14),7,1,7,3) + e3d%dk(pl(13),8,1,8,3) )

            Ah(m,1) = ahsum1

            ! Output DOF 2: y-component of Ah = A*h.
            ahsum2 = 0.d0
            ! Contributions from input displacement/search component n = 1 (x).
            ahsum2 = ahsum2 &
          + h(ibl(1),1 )*( e3d%dk(pl(27),1,2,4,1) + e3d%dk(pl(7),2,2,3,1) + e3d%dk(pl(25),5,2,8,1) + e3d%dk(pl(15),6,2,7,1) ) &
          + h(ibl(2),1)*( e3d%dk(pl(27),1,2,3,1) + e3d%dk(pl(25),5,2,7,1) ) &
          + h(ibl(3),1)*( e3d%dk(pl(27),1,2,2,1) + e3d%dk(pl(5),4,2,3,1) + e3d%dk(pl(13),8,2,7,1) + e3d%dk(pl(25),5,2,6,1) ) &
          + h(ibl(4),1)*( e3d%dk(pl(5),4,2,2,1) + e3d%dk(pl(13),8,2,6,1) ) &
          + h(ibl(5),1)*( e3d%dk(pl(6),3,2,2,1) + e3d%dk(pl(5),4,2,1,1) + e3d%dk(pl(14),7,2,6,1) + e3d%dk(pl(13),8,2,5,1) ) &
          + h(ibl(6),1)*( e3d%dk(pl(6),3,2,1,1) + e3d%dk(pl(14),7,2,5,1) ) &
          + h(ibl(7),1)*( e3d%dk(pl(6),3,2,4,1) + e3d%dk(pl(7),2,2,1,1) + e3d%dk(pl(14),7,2,8,1) + e3d%dk(pl(15),6,2,5,1) ) &
          + h(ibl(8),1)*( e3d%dk(pl(7),2,2,4,1) + e3d%dk(pl(15),6,2,8,1) ) &
          + h(ibl(9),1)*( e3d%dk(pl(25),5,2,4,1) + e3d%dk(pl(15),6,2,3,1) ) &
          + h(ibl(10),1)*( e3d%dk(pl(25),5,2,3,1) ) &
          + h(ibl(11),1)*( e3d%dk(pl(13),8,2,3,1) + e3d%dk(pl(25),5,2,2,1) ) &
          + h(ibl(12),1)*( e3d%dk(pl(13),8,2,2,1) ) &
          + h(ibl(13),1)*( e3d%dk(pl(13),8,2,1,1) + e3d%dk(pl(14),7,2,2,1) ) &
          + h(ibl(14),1)*( e3d%dk(pl(14),7,2,1,1) ) &
          + h(ibl(15),1)*( e3d%dk(pl(14),7,2,4,1) + e3d%dk(pl(15),6,2,1,1) ) &
          + h(ibl(16),1)*( e3d%dk(pl(15),6,2,4,1) ) &
          + h(ibl(17),1)*( e3d%dk(pl(27),1,2,8,1) + e3d%dk(pl(7),2,2,7,1) ) &
          + h(ibl(18),1)*( e3d%dk(pl(27),1,2,7,1) ) &
          + h(ibl(19),1)*( e3d%dk(pl(27),1,2,6,1) + e3d%dk(pl(5),4,2,7,1) ) &
          + h(ibl(20),1)*( e3d%dk(pl(5),4,2,6,1) ) &
          + h(ibl(21),1)*( e3d%dk(pl(5),4,2,5,1) + e3d%dk(pl(6),3,2,6,1) ) &
          + h(ibl(22),1)*( e3d%dk(pl(6),3,2,5,1) ) &
          + h(ibl(23),1)*( e3d%dk(pl(6),3,2,8,1) + e3d%dk(pl(7),2,2,5,1) ) &
          + h(ibl(24),1)*( e3d%dk(pl(7),2,2,8,1) ) &
          + h(ibl(25),1)*( e3d%dk(pl(14),7,2,3,1) + e3d%dk(pl(13),8,2,4,1) + e3d%dk(pl(15),6,2,2,1) + e3d%dk(pl(25),5,2,1,1) ) &
          + h(ibl(26),1)*( e3d%dk(pl(6),3,2,7,1) + e3d%dk(pl(5),4,2,8,1) + e3d%dk(pl(27),1,2,5,1) + e3d%dk(pl(7),2,2,6,1) ) &
          + h(ibl(27),1)*( e3d%dk(pl(27),1,2,1,1) + e3d%dk(pl(7),2,2,2,1) + e3d%dk(pl(6),3,2,3,1) + e3d%dk(pl(5),4,2,4,1) &
                          + e3d%dk(pl(25),5,2,5,1) + e3d%dk(pl(15),6,2,6,1) + e3d%dk(pl(14),7,2,7,1) + e3d%dk(pl(13),8,2,8,1) )
            ! Contributions from input displacement/search component n = 2 (y).
            ahsum2 = ahsum2 &
          + h(ibl(1),2 )*( e3d%dk(pl(27),1,2,4,2) + e3d%dk(pl(7),2,2,3,2) + e3d%dk(pl(25),5,2,8,2) + e3d%dk(pl(15),6,2,7,2) ) &
          + h(ibl(2),2)*( e3d%dk(pl(27),1,2,3,2) + e3d%dk(pl(25),5,2,7,2) ) &
          + h(ibl(3),2)*( e3d%dk(pl(27),1,2,2,2) + e3d%dk(pl(5),4,2,3,2) + e3d%dk(pl(13),8,2,7,2) + e3d%dk(pl(25),5,2,6,2) ) &
          + h(ibl(4),2)*( e3d%dk(pl(5),4,2,2,2) + e3d%dk(pl(13),8,2,6,2) ) &
          + h(ibl(5),2)*( e3d%dk(pl(6),3,2,2,2) + e3d%dk(pl(5),4,2,1,2) + e3d%dk(pl(14),7,2,6,2) + e3d%dk(pl(13),8,2,5,2) ) &
          + h(ibl(6),2)*( e3d%dk(pl(6),3,2,1,2) + e3d%dk(pl(14),7,2,5,2) ) &
          + h(ibl(7),2)*( e3d%dk(pl(6),3,2,4,2) + e3d%dk(pl(7),2,2,1,2) + e3d%dk(pl(14),7,2,8,2) + e3d%dk(pl(15),6,2,5,2) ) &
          + h(ibl(8),2)*( e3d%dk(pl(7),2,2,4,2) + e3d%dk(pl(15),6,2,8,2) ) &
          + h(ibl(9),2)*( e3d%dk(pl(25),5,2,4,2) + e3d%dk(pl(15),6,2,3,2) ) &
          + h(ibl(10),2)*( e3d%dk(pl(25),5,2,3,2) ) &
          + h(ibl(11),2)*( e3d%dk(pl(13),8,2,3,2) + e3d%dk(pl(25),5,2,2,2) ) &
          + h(ibl(12),2)*( e3d%dk(pl(13),8,2,2,2) ) &
          + h(ibl(13),2)*( e3d%dk(pl(13),8,2,1,2) + e3d%dk(pl(14),7,2,2,2) ) &
          + h(ibl(14),2)*( e3d%dk(pl(14),7,2,1,2) ) &
          + h(ibl(15),2)*( e3d%dk(pl(14),7,2,4,2) + e3d%dk(pl(15),6,2,1,2) ) &
          + h(ibl(16),2)*( e3d%dk(pl(15),6,2,4,2) ) &
          + h(ibl(17),2)*( e3d%dk(pl(27),1,2,8,2) + e3d%dk(pl(7),2,2,7,2) ) &
          + h(ibl(18),2)*( e3d%dk(pl(27),1,2,7,2) ) &
          + h(ibl(19),2)*( e3d%dk(pl(27),1,2,6,2) + e3d%dk(pl(5),4,2,7,2) ) &
          + h(ibl(20),2)*( e3d%dk(pl(5),4,2,6,2) ) &
          + h(ibl(21),2)*( e3d%dk(pl(5),4,2,5,2) + e3d%dk(pl(6),3,2,6,2) ) &
          + h(ibl(22),2)*( e3d%dk(pl(6),3,2,5,2) ) &
          + h(ibl(23),2)*( e3d%dk(pl(6),3,2,8,2) + e3d%dk(pl(7),2,2,5,2) ) &
          + h(ibl(24),2)*( e3d%dk(pl(7),2,2,8,2) ) &
          + h(ibl(25),2)*( e3d%dk(pl(14),7,2,3,2) + e3d%dk(pl(13),8,2,4,2) + e3d%dk(pl(15),6,2,2,2) + e3d%dk(pl(25),5,2,1,2) ) &
          + h(ibl(26),2)*( e3d%dk(pl(6),3,2,7,2) + e3d%dk(pl(5),4,2,8,2) + e3d%dk(pl(27),1,2,5,2) + e3d%dk(pl(7),2,2,6,2) ) &
          + h(ibl(27),2)*( e3d%dk(pl(27),1,2,1,2) + e3d%dk(pl(7),2,2,2,2) + e3d%dk(pl(6),3,2,3,2) + e3d%dk(pl(5),4,2,4,2) &
                          + e3d%dk(pl(25),5,2,5,2) + e3d%dk(pl(15),6,2,6,2) + e3d%dk(pl(14),7,2,7,2) + e3d%dk(pl(13),8,2,8,2) )
            ! Contributions from input displacement/search component n = 3 (z).
            ahsum2 = ahsum2 &
          + h(ibl(1),3 )*( e3d%dk(pl(27),1,2,4,3) + e3d%dk(pl(7),2,2,3,3) + e3d%dk(pl(25),5,2,8,3) + e3d%dk(pl(15),6,2,7,3) ) &
          + h(ibl(2),3)*( e3d%dk(pl(27),1,2,3,3) + e3d%dk(pl(25),5,2,7,3) ) &
          + h(ibl(3),3)*( e3d%dk(pl(27),1,2,2,3) + e3d%dk(pl(5),4,2,3,3) + e3d%dk(pl(13),8,2,7,3) + e3d%dk(pl(25),5,2,6,3) ) &
          + h(ibl(4),3)*( e3d%dk(pl(5),4,2,2,3) + e3d%dk(pl(13),8,2,6,3) ) &
          + h(ibl(5),3)*( e3d%dk(pl(6),3,2,2,3) + e3d%dk(pl(5),4,2,1,3) + e3d%dk(pl(14),7,2,6,3) + e3d%dk(pl(13),8,2,5,3) ) &
          + h(ibl(6),3)*( e3d%dk(pl(6),3,2,1,3) + e3d%dk(pl(14),7,2,5,3) ) &
          + h(ibl(7),3)*( e3d%dk(pl(6),3,2,4,3) + e3d%dk(pl(7),2,2,1,3) + e3d%dk(pl(14),7,2,8,3) + e3d%dk(pl(15),6,2,5,3) ) &
          + h(ibl(8),3)*( e3d%dk(pl(7),2,2,4,3) + e3d%dk(pl(15),6,2,8,3) ) &
          + h(ibl(9),3)*( e3d%dk(pl(25),5,2,4,3) + e3d%dk(pl(15),6,2,3,3) ) &
          + h(ibl(10),3)*( e3d%dk(pl(25),5,2,3,3) ) &
          + h(ibl(11),3)*( e3d%dk(pl(13),8,2,3,3) + e3d%dk(pl(25),5,2,2,3) ) &
          + h(ibl(12),3)*( e3d%dk(pl(13),8,2,2,3) ) &
          + h(ibl(13),3)*( e3d%dk(pl(13),8,2,1,3) + e3d%dk(pl(14),7,2,2,3) ) &
          + h(ibl(14),3)*( e3d%dk(pl(14),7,2,1,3) ) &
          + h(ibl(15),3)*( e3d%dk(pl(14),7,2,4,3) + e3d%dk(pl(15),6,2,1,3) ) &
          + h(ibl(16),3)*( e3d%dk(pl(15),6,2,4,3) ) &
          + h(ibl(17),3)*( e3d%dk(pl(27),1,2,8,3) + e3d%dk(pl(7),2,2,7,3) ) &
          + h(ibl(18),3)*( e3d%dk(pl(27),1,2,7,3) ) &
          + h(ibl(19),3)*( e3d%dk(pl(27),1,2,6,3) + e3d%dk(pl(5),4,2,7,3) ) &
          + h(ibl(20),3)*( e3d%dk(pl(5),4,2,6,3) ) &
          + h(ibl(21),3)*( e3d%dk(pl(5),4,2,5,3) + e3d%dk(pl(6),3,2,6,3) ) &
          + h(ibl(22),3)*( e3d%dk(pl(6),3,2,5,3) ) &
          + h(ibl(23),3)*( e3d%dk(pl(6),3,2,8,3) + e3d%dk(pl(7),2,2,5,3) ) &
          + h(ibl(24),3)*( e3d%dk(pl(7),2,2,8,3) ) &
          + h(ibl(25),3)*( e3d%dk(pl(14),7,2,3,3) + e3d%dk(pl(13),8,2,4,3) + e3d%dk(pl(15),6,2,2,3) + e3d%dk(pl(25),5,2,1,3) ) &
          + h(ibl(26),3)*( e3d%dk(pl(6),3,2,7,3) + e3d%dk(pl(5),4,2,8,3) + e3d%dk(pl(27),1,2,5,3) + e3d%dk(pl(7),2,2,6,3) ) &
          + h(ibl(27),3)*( e3d%dk(pl(27),1,2,1,3) + e3d%dk(pl(7),2,2,2,3) + e3d%dk(pl(6),3,2,3,3) + e3d%dk(pl(5),4,2,4,3) &
                          + e3d%dk(pl(25),5,2,5,3) + e3d%dk(pl(15),6,2,6,3) + e3d%dk(pl(14),7,2,7,3) + e3d%dk(pl(13),8,2,8,3) )

            Ah(m,2) = ahsum2

            ! Output DOF 3: z-component of Ah = A*h.
            ahsum3 = 0.d0
            ! Contributions from input displacement/search component n = 1 (x).
            ahsum3 = ahsum3 &
          + h(ibl(1),1 )*( e3d%dk(pl(27),1,3,4,1) + e3d%dk(pl(7),2,3,3,1) + e3d%dk(pl(25),5,3,8,1) + e3d%dk(pl(15),6,3,7,1) ) &
          + h(ibl(2),1)*( e3d%dk(pl(27),1,3,3,1) + e3d%dk(pl(25),5,3,7,1) ) &
          + h(ibl(3),1)*( e3d%dk(pl(27),1,3,2,1) + e3d%dk(pl(5),4,3,3,1) + e3d%dk(pl(13),8,3,7,1) + e3d%dk(pl(25),5,3,6,1) ) &
          + h(ibl(4),1)*( e3d%dk(pl(5),4,3,2,1) + e3d%dk(pl(13),8,3,6,1) ) &
          + h(ibl(5),1)*( e3d%dk(pl(6),3,3,2,1) + e3d%dk(pl(5),4,3,1,1) + e3d%dk(pl(14),7,3,6,1) + e3d%dk(pl(13),8,3,5,1) ) &
          + h(ibl(6),1)*( e3d%dk(pl(6),3,3,1,1) + e3d%dk(pl(14),7,3,5,1) ) &
          + h(ibl(7),1)*( e3d%dk(pl(6),3,3,4,1) + e3d%dk(pl(7),2,3,1,1) + e3d%dk(pl(14),7,3,8,1) + e3d%dk(pl(15),6,3,5,1) ) &
          + h(ibl(8),1)*( e3d%dk(pl(7),2,3,4,1) + e3d%dk(pl(15),6,3,8,1) ) &
          + h(ibl(9),1)*( e3d%dk(pl(25),5,3,4,1) + e3d%dk(pl(15),6,3,3,1) ) &
          + h(ibl(10),1)*( e3d%dk(pl(25),5,3,3,1) ) &
          + h(ibl(11),1)*( e3d%dk(pl(13),8,3,3,1) + e3d%dk(pl(25),5,3,2,1) ) &
          + h(ibl(12),1)*( e3d%dk(pl(13),8,3,2,1) ) &
          + h(ibl(13),1)*( e3d%dk(pl(13),8,3,1,1) + e3d%dk(pl(14),7,3,2,1) ) &
          + h(ibl(14),1)*( e3d%dk(pl(14),7,3,1,1) ) &
          + h(ibl(15),1)*( e3d%dk(pl(14),7,3,4,1) + e3d%dk(pl(15),6,3,1,1) ) &
          + h(ibl(16),1)*( e3d%dk(pl(15),6,3,4,1) ) &
          + h(ibl(17),1)*( e3d%dk(pl(27),1,3,8,1) + e3d%dk(pl(7),2,3,7,1) ) &
          + h(ibl(18),1)*( e3d%dk(pl(27),1,3,7,1) ) &
          + h(ibl(19),1)*( e3d%dk(pl(27),1,3,6,1) + e3d%dk(pl(5),4,3,7,1) ) &
          + h(ibl(20),1)*( e3d%dk(pl(5),4,3,6,1) ) &
          + h(ibl(21),1)*( e3d%dk(pl(5),4,3,5,1) + e3d%dk(pl(6),3,3,6,1) ) &
          + h(ibl(22),1)*( e3d%dk(pl(6),3,3,5,1) ) &
          + h(ibl(23),1)*( e3d%dk(pl(6),3,3,8,1) + e3d%dk(pl(7),2,3,5,1) ) &
          + h(ibl(24),1)*( e3d%dk(pl(7),2,3,8,1) ) &
          + h(ibl(25),1)*( e3d%dk(pl(14),7,3,3,1) + e3d%dk(pl(13),8,3,4,1) + e3d%dk(pl(15),6,3,2,1) + e3d%dk(pl(25),5,3,1,1) ) &
          + h(ibl(26),1)*( e3d%dk(pl(6),3,3,7,1) + e3d%dk(pl(5),4,3,8,1) + e3d%dk(pl(27),1,3,5,1) + e3d%dk(pl(7),2,3,6,1) ) &
          + h(ibl(27),1)*( e3d%dk(pl(27),1,3,1,1) + e3d%dk(pl(7),2,3,2,1) + e3d%dk(pl(6),3,3,3,1) + e3d%dk(pl(5),4,3,4,1) &
                          + e3d%dk(pl(25),5,3,5,1) + e3d%dk(pl(15),6,3,6,1) + e3d%dk(pl(14),7,3,7,1) + e3d%dk(pl(13),8,3,8,1) )
            ! Contributions from input displacement/search component n = 2 (y).
            ahsum3 = ahsum3 &
          + h(ibl(1),2 )*( e3d%dk(pl(27),1,3,4,2) + e3d%dk(pl(7),2,3,3,2) + e3d%dk(pl(25),5,3,8,2) + e3d%dk(pl(15),6,3,7,2) ) &
          + h(ibl(2),2)*( e3d%dk(pl(27),1,3,3,2) + e3d%dk(pl(25),5,3,7,2) ) &
          + h(ibl(3),2)*( e3d%dk(pl(27),1,3,2,2) + e3d%dk(pl(5),4,3,3,2) + e3d%dk(pl(13),8,3,7,2) + e3d%dk(pl(25),5,3,6,2) ) &
          + h(ibl(4),2)*( e3d%dk(pl(5),4,3,2,2) + e3d%dk(pl(13),8,3,6,2) ) &
          + h(ibl(5),2)*( e3d%dk(pl(6),3,3,2,2) + e3d%dk(pl(5),4,3,1,2) + e3d%dk(pl(14),7,3,6,2) + e3d%dk(pl(13),8,3,5,2) ) &
          + h(ibl(6),2)*( e3d%dk(pl(6),3,3,1,2) + e3d%dk(pl(14),7,3,5,2) ) &
          + h(ibl(7),2)*( e3d%dk(pl(6),3,3,4,2) + e3d%dk(pl(7),2,3,1,2) + e3d%dk(pl(14),7,3,8,2) + e3d%dk(pl(15),6,3,5,2) ) &
          + h(ibl(8),2)*( e3d%dk(pl(7),2,3,4,2) + e3d%dk(pl(15),6,3,8,2) ) &
          + h(ibl(9),2)*( e3d%dk(pl(25),5,3,4,2) + e3d%dk(pl(15),6,3,3,2) ) &
          + h(ibl(10),2)*( e3d%dk(pl(25),5,3,3,2) ) &
          + h(ibl(11),2)*( e3d%dk(pl(13),8,3,3,2) + e3d%dk(pl(25),5,3,2,2) ) &
          + h(ibl(12),2)*( e3d%dk(pl(13),8,3,2,2) ) &
          + h(ibl(13),2)*( e3d%dk(pl(13),8,3,1,2) + e3d%dk(pl(14),7,3,2,2) ) &
          + h(ibl(14),2)*( e3d%dk(pl(14),7,3,1,2) ) &
          + h(ibl(15),2)*( e3d%dk(pl(14),7,3,4,2) + e3d%dk(pl(15),6,3,1,2) ) &
          + h(ibl(16),2)*( e3d%dk(pl(15),6,3,4,2) ) &
          + h(ibl(17),2)*( e3d%dk(pl(27),1,3,8,2) + e3d%dk(pl(7),2,3,7,2) ) &
          + h(ibl(18),2)*( e3d%dk(pl(27),1,3,7,2) ) &
          + h(ibl(19),2)*( e3d%dk(pl(27),1,3,6,2) + e3d%dk(pl(5),4,3,7,2) ) &
          + h(ibl(20),2)*( e3d%dk(pl(5),4,3,6,2) ) &
          + h(ibl(21),2)*( e3d%dk(pl(5),4,3,5,2) + e3d%dk(pl(6),3,3,6,2) ) &
          + h(ibl(22),2)*( e3d%dk(pl(6),3,3,5,2) ) &
          + h(ibl(23),2)*( e3d%dk(pl(6),3,3,8,2) + e3d%dk(pl(7),2,3,5,2) ) &
          + h(ibl(24),2)*( e3d%dk(pl(7),2,3,8,2) ) &
          + h(ibl(25),2)*( e3d%dk(pl(14),7,3,3,2) + e3d%dk(pl(13),8,3,4,2) + e3d%dk(pl(15),6,3,2,2) + e3d%dk(pl(25),5,3,1,2) ) &
          + h(ibl(26),2)*( e3d%dk(pl(6),3,3,7,2) + e3d%dk(pl(5),4,3,8,2) + e3d%dk(pl(27),1,3,5,2) + e3d%dk(pl(7),2,3,6,2) ) &
          + h(ibl(27),2)*( e3d%dk(pl(27),1,3,1,2) + e3d%dk(pl(7),2,3,2,2) + e3d%dk(pl(6),3,3,3,2) + e3d%dk(pl(5),4,3,4,2) &
                          + e3d%dk(pl(25),5,3,5,2) + e3d%dk(pl(15),6,3,6,2) + e3d%dk(pl(14),7,3,7,2) + e3d%dk(pl(13),8,3,8,2) )
            ! Contributions from input displacement/search component n = 3 (z).
            ahsum3 = ahsum3 &
          + h(ibl(1),3 )*( e3d%dk(pl(27),1,3,4,3) + e3d%dk(pl(7),2,3,3,3) + e3d%dk(pl(25),5,3,8,3) + e3d%dk(pl(15),6,3,7,3) ) &
          + h(ibl(2),3)*( e3d%dk(pl(27),1,3,3,3) + e3d%dk(pl(25),5,3,7,3) ) &
          + h(ibl(3),3)*( e3d%dk(pl(27),1,3,2,3) + e3d%dk(pl(5),4,3,3,3) + e3d%dk(pl(13),8,3,7,3) + e3d%dk(pl(25),5,3,6,3) ) &
          + h(ibl(4),3)*( e3d%dk(pl(5),4,3,2,3) + e3d%dk(pl(13),8,3,6,3) ) &
          + h(ibl(5),3)*( e3d%dk(pl(6),3,3,2,3) + e3d%dk(pl(5),4,3,1,3) + e3d%dk(pl(14),7,3,6,3) + e3d%dk(pl(13),8,3,5,3) ) &
          + h(ibl(6),3)*( e3d%dk(pl(6),3,3,1,3) + e3d%dk(pl(14),7,3,5,3) ) &
          + h(ibl(7),3)*( e3d%dk(pl(6),3,3,4,3) + e3d%dk(pl(7),2,3,1,3) + e3d%dk(pl(14),7,3,8,3) + e3d%dk(pl(15),6,3,5,3) ) &
          + h(ibl(8),3)*( e3d%dk(pl(7),2,3,4,3) + e3d%dk(pl(15),6,3,8,3) ) &
          + h(ibl(9),3)*( e3d%dk(pl(25),5,3,4,3) + e3d%dk(pl(15),6,3,3,3) ) &
          + h(ibl(10),3)*( e3d%dk(pl(25),5,3,3,3) ) &
          + h(ibl(11),3)*( e3d%dk(pl(13),8,3,3,3) + e3d%dk(pl(25),5,3,2,3) ) &
          + h(ibl(12),3)*( e3d%dk(pl(13),8,3,2,3) ) &
          + h(ibl(13),3)*( e3d%dk(pl(13),8,3,1,3) + e3d%dk(pl(14),7,3,2,3) ) &
          + h(ibl(14),3)*( e3d%dk(pl(14),7,3,1,3) ) &
          + h(ibl(15),3)*( e3d%dk(pl(14),7,3,4,3) + e3d%dk(pl(15),6,3,1,3) ) &
          + h(ibl(16),3)*( e3d%dk(pl(15),6,3,4,3) ) &
          + h(ibl(17),3)*( e3d%dk(pl(27),1,3,8,3) + e3d%dk(pl(7),2,3,7,3) ) &
          + h(ibl(18),3)*( e3d%dk(pl(27),1,3,7,3) ) &
          + h(ibl(19),3)*( e3d%dk(pl(27),1,3,6,3) + e3d%dk(pl(5),4,3,7,3) ) &
          + h(ibl(20),3)*( e3d%dk(pl(5),4,3,6,3) ) &
          + h(ibl(21),3)*( e3d%dk(pl(5),4,3,5,3) + e3d%dk(pl(6),3,3,6,3) ) &
          + h(ibl(22),3)*( e3d%dk(pl(6),3,3,5,3) ) &
          + h(ibl(23),3)*( e3d%dk(pl(6),3,3,8,3) + e3d%dk(pl(7),2,3,5,3) ) &
          + h(ibl(24),3)*( e3d%dk(pl(7),2,3,8,3) ) &
          + h(ibl(25),3)*( e3d%dk(pl(14),7,3,3,3) + e3d%dk(pl(13),8,3,4,3) + e3d%dk(pl(15),6,3,2,3) + e3d%dk(pl(25),5,3,1,3) ) &
          + h(ibl(26),3)*( e3d%dk(pl(6),3,3,7,3) + e3d%dk(pl(5),4,3,8,3) + e3d%dk(pl(27),1,3,5,3) + e3d%dk(pl(7),2,3,6,3) ) &
          + h(ibl(27),3)*( e3d%dk(pl(27),1,3,1,3) + e3d%dk(pl(7),2,3,2,3) + e3d%dk(pl(6),3,3,3,3) + e3d%dk(pl(5),4,3,4,3) &
                          + e3d%dk(pl(25),5,3,5,3) + e3d%dk(pl(15),6,3,6,3) + e3d%dk(pl(14),7,3,7,3) + e3d%dk(pl(13),8,3,8,3) )

            Ah(m,3) = ahsum3

          end do
          !$omp end parallel do

        end do ! blocks

        ! ====================================================================
        ! ---- 9b. Curvature and Step Size ------------------------------------
        ! Compute the curvature along the current search direction:
        !
        !   hAh = h^T*A*h
        !
        ! The PCG step size is:
        !
        !   lambda = gg / hAh
        !
        ! where:
        !
        !   gg = z^T*gb
        !
        ! Because the code stores the energy gradient, the displacement update
        ! below uses:
        !
        !   u <- u - lambda*h
        ! ====================================================================

        ! Remove any inactive-node entries from Ah. These nodes are excluded
        ! from the solve and should not contribute to inner products.
        !$omp parallel do collapse(2) schedule(guided)
        do i = 1, ndof
          do m = 1, ns
            if (.not. e3d%is_active(m)) Ah(m,i) = 0.d0
          end do
        end do
        !$omp end parallel do
				
        ! Compute hAh = sum(h .* Ah)
        hAh = 0.d0
        !$omp parallel do collapse(2) schedule(guided) reduction(+:hAh)		
        do i = 1, ndof
		  do m = 1, ns
            hAh = hAh + h(m,i) * Ah(m,i)
          end do
        end do
        !$omp end parallel do

        ! Protect against division by zero or loss of positive curvature.
        if(dabs(hAh) < 1d-14) then
          print *,"WARNING: hAh is nearly zero, breaking CG step."
          exit
        end if
		
		
        ! Compute PCG step size, lambda, often denoted alpha:
        !   lambda = (z^T*gb)/(h^T*A*h)
        lambda = gg / hAh

        ! ====================================================================
        ! ---- 9c. Displacement and Gradient Update ---------------------------
        ! Update displacement and gradient along the current search direction:
        !
        !   u  <- u  - lambda*h
        !   gb <- gb - lambda*A*h
        !
        ! The minus sign appears because gb is the energy gradient.
        ! ====================================================================
        !$omp parallel do collapse(2) schedule(guided)		
        do i = 1, ndof
		  do m = 1, ns
            u(m,i) = u(m,i) - lambda * h(m,i)
            gb(m,i) = gb(m,i) - lambda * Ah(m,i)
		  end do
        end do
        !$omp end parallel do
		
        ! Re-enforce zero displacement and zero gradient on inactive pore-only
        ! nodes after the update.
        !$omp parallel do collapse(2) schedule(guided)
        do i = 1, ndof
          do m = 1, ns
            if (.not. e3d%is_active(m)) then
              u(m,i)  = 0.d0
              gb(m,i) = 0.d0
            end if
          end do
        end do
        !$omp end parallel do

        ! ====================================================================
        ! ---- 9d. Apply Preconditioner and Update Residual Measure -----------
        ! Apply the block-Jacobi preconditioner to the updated gradient:
        !
        !   z = M^{-1}*gb
        !
        ! Then compute the new preconditioned residual inner product:
        !
        !   gg = z^T*gb
        ! ====================================================================
        !$omp parallel do collapse(2) schedule(guided)		
        do i = 1, ndof
		  do m = 1, ns
            z(m,i) = M_block_inv(m,i,1)*gb(m,1) + &
		             M_block_inv(m,i,2)*gb(m,2) + &
			         M_block_inv(m,i,3)*gb(m,3)
          end do
        end do
        !$omp end parallel do
    
        ! Store previous residual inner product and compute the updated value.
        gglast = gg
        gg = 0.d0
        !$omp parallel do collapse(2) schedule(guided) reduction(+:gg)		
        do i = 1, ndof
		  do m = 1, ns
            gg = gg + z(m,i) * gb(m,i)
          end do
        end do
        !$omp end parallel do    

        ! Compute relative residual norm for convergence and reset checks.
		rel_res = sqrt(gg/gg0)
        ! Store local iteration history for output.
        history_step(Lstep) = Lstep
        history_res(Lstep)  = rel_res		

        ! ====================================================================
        ! ---- 9e. Residual-Spike Safeguard -----------------------------------
        ! If the relative residual increases sharply, the conjugate search
        ! direction may have lost effectiveness due to numerical roundoff,
        ! ill-conditioning, or pore-induced stiffness contrast.
        !
        ! In that case, reset the search direction to the preconditioned gradient:
        !
        !   h = z
        !
        ! Because the displacement update is u <- u - lambda*h, this reset
        ! corresponds to a preconditioned steepest-descent step.
        ! ====================================================================
        if (rel_res > RESET_RATIO * rel_res_prev .and. rel_res > RESET_ABS_TOL .and. Lstep > 5) then
          print *, 'WARNING: Residual increased!'
          print *, '  Previous rel_res =', rel_res_prev
          print *, '  Current rel_res  =', rel_res
          print *, '  Resetting search direction to steepest descent.'
          
          ! Reset direction to the preconditioned gradient.
          ! With u <- u - lambda*h, this is preconditioned steepest descent.
          !$omp parallel do collapse(2) schedule(guided)		
          do i = 1, ndof
            do m = 1, ns
              h(m,i) = z(m,i)
            end do
          end do
          !$omp end parallel do
          
          ! Reset residual tracking after direction reset.
          rel_res_min = rel_res
          reset_count = reset_count + 1
          
          print *, '  Direction reset count:', reset_count
		  
        else
          ! ==================================================================
          ! ---- 9f. Conjugate Direction Update -------------------------------
          ! Update the PCG search direction using the Fletcher-Reeves form:
          !
          !   gamma = gg_new / gg_old
          !   h     = z + gamma*h
          !
          ! Here gg = z^T*gb is the preconditioned residual inner product.
          ! ==================================================================
          gamma = gg / gglast
          !$omp parallel do collapse(2) schedule(guided)		
          do i = 1, ndof
		    do m = 1, ns
		      h(m,i) = z(m,i) + gamma * h(m,i)
		    end do
          end do
          !$omp end parallel do 

          ! Track the best residual observed after the most recent reset.
          if (rel_res < rel_res_min) then
            rel_res_min = rel_res
          end if		
        end if
        
        ! Store residual for the next spike-detection check.
        rel_res_prev = rel_res        
        
        ! Exit early if the relative residual tolerance is met.	        
        if (rel_res < tol) exit
    
        print *, 'Number of conjugate steps_local iteration =', Lstep
        print *, 'rel. residual norm=', rel_res


      end do ! CG steps
	  
      ! ======================================================================
      ! ---- 10. Write History and Copy Local State Back to e3d ---------------
      ! Write detailed PCG convergence history to disk and copy the updated
      ! displacement, gradient, and search direction fields back to the global
      ! solver data structure.
      ! ======================================================================
      do i = 1, Lstep
         write(100, '(I8, 4X, E20.8)') history_step(i), history_res(i)
      end do

      ! Close detailed PCG history file.
      close(100)
	  
      ! Persist final local solver state back to e3d.
      !$omp parallel do collapse(2) schedule(guided)  
      do i = 1, ndof
	    do m = 1, ns
          e3d%u(m,i)  = u(m,i)
          e3d%gb(m,i) = gb(m,i)
          e3d%h(m,i)  = h(m,i)
        end do
      end do
      !$omp end parallel do	  


      ! ======================================================================
      ! ---- 11. Memory Cleanup ----------------------------------------------
      ! Release local PCG work arrays before returning to the main program.
      ! ======================================================================
      deallocate( u, gb, h, Ah )
      deallocate( M_block_inv, z )  
      deallocate( history_step, history_res )	      

      end subroutine dembx_OpenMP
	  
!-------------------------------------------------------------------------  
      pure subroutine invert3x3(a, ainv)
      ! ==============================================================================
      ! SUBROUTINE: invert3x3
      ! ------------------------------------------------------------------------------
      ! PURPOSE:
      !   Compute the inverse of a 3 x 3 matrix using an explicit adjugate/cofactor
      !   formula.
      !
      !   In this solver, the routine is used primarily to invert the local
      !   3 x 3 block-Jacobi preconditioner matrix at each grid node in the PCG
      !   solver.
      !
      ! ARGUMENTS:
      !   a
      !     Input 3 x 3 matrix.
      !
      !   ainv
      !     Output inverse matrix. If the determinant is below the singularity
      !     threshold, ainv is returned as the zero matrix.
      !
      ! METHOD:
      !   The inverse is computed as:
      !
      !      A^{-1} = adj(A) / det(A)
      !
      !   where adj(A) is the adjugate matrix. The entries are evaluated
      !   explicitly to avoid the overhead of a general-purpose matrix inversion
      !   routine for this small fixed-size case.
      !
      ! SINGULARITY HANDLING:
      !   If abs(det(A)) < 1.0d-14, the matrix is treated as singular or nearly
      !   singular and the routine returns a zero matrix.
      !
      !   This behavior is useful for inactive pore-only nodes, whose stiffness
      !   blocks may be zero. For active solid-connected nodes, frequent singular
      !   blocks may indicate a conditioning or connectivity issue.
      !
      ! THREAD SAFETY:
      !   The routine is declared PURE. It has no side effects, performs no I/O,
      !   and modifies only its explicit output argument. It is therefore safe to
      !   call inside OpenMP-parallel loops.
      !
      ! NOTE:
      !   This direct formula is efficient for 3 x 3 matrices but is not intended
      !   as a general-purpose numerically robust matrix inversion routine.
      ! ==============================================================================	  
      use elas3d_mod
      implicit none

      ! Input matrix and output inverse.	  
      real(dp), intent(in)  :: a(3,3)
      real(dp), intent(out) :: ainv(3,3)
	  
      ! Determinant and reciprocal determinant.	  
      real(dp) :: det, inv_det
	  
      ! Entries of the adjugate matrix.
      ! The variable names correspond to their final location in ainv:
      !   ainv(i,j) = cij / det
      real(dp) :: c11, c12, c13
      real(dp) :: c21, c22, c23
      real(dp) :: c31, c32, c33

      ! ======================================================================
      ! ---- 1. Compute First Column of the Adjugate --------------------------
      ! Compute the first column of adj(A), which will become the first column
      ! of A^{-1} after scaling by 1/det(A).
      !
      ! These three entries are also sufficient to compute det(A) by expanding
      ! along the first row of A.
      ! ======================================================================
      c11 =  a(2,2)*a(3,3) - a(2,3)*a(3,2)
      c21 = -a(2,1)*a(3,3) + a(2,3)*a(3,1)  ! Note: Sign absorbed
      c31 =  a(2,1)*a(3,2) - a(2,2)*a(3,1)

      ! ======================================================================
      ! ---- 2. Compute Determinant ------------------------------------------
      ! Use the first-row expansion:
      !
      !   det(A) = a11*adj11 + a12*adj21 + a13*adj31
      !
      ! where c11, c21, and c31 are the first-column entries of adj(A).
      ! ======================================================================
      det = a(1,1)*c11 + a(1,2)*c21 + a(1,3)*c31
			

      ! ======================================================================
      ! ---- 3. Singularity Check --------------------------------------------
      ! If the determinant is too small, return a zero inverse. This avoids
      ! division by a near-zero number.
      !
      ! In the PCG preconditioner, zero blocks can occur for inactive pore-only
      ! nodes. Those nodes are constrained elsewhere by the active-node mask.
      ! ======================================================================
      if (abs(det) < 1d-14) then
        ainv(:,:) = 0.d0
        return
      end if
	  
      ! Precompute reciprocal determinant so the final scaling uses
      ! multiplication instead of repeated division.
      inv_det = 1.0d0 / det	  
      
      ! ======================================================================
      ! ---- 4. Compute Remaining Adjugate Entries ----------------------------
      ! Compute the remaining entries of adj(A). The naming convention follows
      ! the output inverse location:
      !
      !   ainv(i,j) = cij * inv_det
      !
      ! ======================================================================
      ! Column 2 of Inverse (Row 2 cofactors of A)
      c12 = -a(1,2)*a(3,3) + a(1,3)*a(3,2)
      c22 =  a(1,1)*a(3,3) - a(1,3)*a(3,1)
      c32 = -a(1,1)*a(3,2) + a(1,2)*a(3,1)

      ! Column 3 of Inverse (Row 3 cofactors of A)
      c13 =  a(1,2)*a(2,3) - a(1,3)*a(2,2)
      c23 = -a(1,1)*a(2,3) + a(1,3)*a(2,1)
      c33 =  a(1,1)*a(2,2) - a(1,2)*a(2,1)	  
		
      ! ======================================================================
      ! ---- 5. Assemble Inverse Matrix ---------------------------------------
      ! Scale the adjugate entries by 1/det(A) to obtain A^{-1}.
      !
      ! Assignments are grouped by column, matching Fortran's column-major
      ! storage order.
      ! ======================================================================
      ainv(1,1) = c11 * inv_det
      ainv(2,1) = c21 * inv_det
      ainv(3,1) = c31 * inv_det

      ainv(1,2) = c12 * inv_det
      ainv(2,2) = c22 * inv_det
      ainv(3,2) = c32 * inv_det

      ainv(1,3) = c13 * inv_det
      ainv(2,3) = c23 * inv_det
      ainv(3,3) = c33 * inv_det

      end subroutine invert3x3		  
!-------------------------------------------------------------------------
      subroutine stress_OpenMP(e3d)
      ! ==============================================================================
      ! SUBROUTINE: stress_OpenMP
      ! ------------------------------------------------------------------------------
      ! PURPOSE:
      !   Compute volume-averaged stress and strain components for the current
      !   displacement field.
      !
      !   This routine post-processes the converged or intermediate displacement
      !   field e3d%u. For each grid site/element, it gathers the eight local nodal
      !   displacements, applies periodic macroscopic-strain jump corrections where
      !   needed, computes the local strain at the voxel/element center, computes
      !   stress using the phase stiffness tensor, and accumulates global averages.
      !
      ! ARGUMENTS:
      !   e3d
      !     Main solver data structure.
      !
      !     Input components:
      !       e3d%u
      !         Current nodal displacement field.
      !
      !       e3d%cmod
      !         6 x 6 constitutive stiffness matrix for each phase.
      !
      !       e3d%pix
      !         Phase label at each grid site/element.
      !
      !       e3d%ib
      !         Periodic 27-neighbor connectivity table.
      !
      !       e3d%exx, e3d%eyy, e3d%ezz, e3d%exz, e3d%eyz, e3d%exy
      !         Prescribed macroscopic strain components. The shear quantities
      !         stored here are tensorial shear components used to generate
      !         displacement jumps.
      !
      !     Output components:
      !       e3d%strxx, e3d%stryy, e3d%strzz, e3d%strxz, e3d%stryz, e3d%strxy
      !         Volume-averaged stress components.
      !
      !       e3d%sxx, e3d%syy, e3d%szz, e3d%sxz, e3d%syz, e3d%sxy
      !         Volume-averaged strain components.
      !
      ! VOIGT ORDERING:
      !   Stress:
      !     1 = xx, 2 = yy, 3 = zz, 4 = xz, 5 = yz, 6 = xy
      !
      !   Strain:
      !     1 = eps_xx
      !     2 = eps_yy
      !     3 = eps_zz
      !     4 = gamma_xz = du_x/dz + du_z/dx
      !     5 = gamma_yz = du_y/dz + du_z/dy
      !     6 = gamma_xy = du_x/dy + du_y/dx
      !
      !   Thus, the shear strain components accumulated here are engineering
      !   shear strains.
      !
      ! PORE HANDLING:
      !   Sites with phase ID PORE_PHASE are skipped and contribute zero to the
      !   accumulated stress and strain sums.
      !
      ! AVERAGING CONVENTION:
      !   Final averages are divided by ns, the total number of grid sites. This
      !   means pores are included in the domain average as zero-contribution
      !   regions. If a solid-only average is desired, the normalization should be
      !   changed to the number of non-pore sites.
      !
      ! NUMERICAL METHOD:
      !   - Uses 8-node trilinear hexahedral element kinematics.
      !   - Evaluates strain at the element/voxel center.
      !   - Applies periodic boundary displacement jumps before strain evaluation.
      !   - Uses OpenMP reductions to accumulate volume averages.
      ! ==============================================================================
      use elas3d_mod
      implicit none

      ! Main solver state structure.
      type(elas3d_data_type), intent(inout) :: e3d

      ! Loop indices and local element/grid-site indices.
      integer :: i, j, k, m,  mm, n8, n3, n
	  integer :: m_local, k_local, j_local, i_local, pixm

      ! Periodic connectivity and indexing.
      integer :: ibl(27)
      integer :: nxy
	  
      ! Shape-function derivatives at the element center.
      real(dp) :: dndx(nnode_fe), dndy(nnode_fe), dndz(nnode_fe)
	  
      ! Strain-displacement matrix:
      !   es(strain_component, local_node, displacement_component)	  
	  real(dp) :: es(6,nnode_fe,ndof)
	  
      ! Local nodal displacement array for one element.
	  real(dp) :: uu(nnode_fe,ndof)
	  
      ! Optional/local scalar placeholders.
      ! Note: str11, s11, etc. are declared here but the current implementation
      ! uses the *_t accumulators below for local stress/strain contributions.	  
      real(dp) :: str11, str22, str33, str13, str23, str12
      real(dp) :: s11, s22, s33, s13, s23, s12
	  
      ! Local per-site stress and strain accumulators.
      real(dp) :: str11_t, str22_t, str33_t, str13_t, str23_t, str12_t
      real(dp) :: s11_t, s22_t, s33_t, s13_t, s23_t, s12_t      
	  
      ! Global stress and strain sums accumulated using OpenMP reductions.
      real(dp) :: strxx, stryy, strzz, strxz, stryz, strxy
      real(dp) :: sxx, syy, szz, sxz, syz, sxy
	  
      ! Prescribed macroscopic strain components unpacked from e3d.
	  real(dp) :: exx, eyy, ezz, exz, eyz, exy



      ! ======================================================================
      ! ---- 1. Initialization and B-Matrix Setup -----------------------------
      ! Unpack prescribed macroscopic strain components and construct the
      ! strain-displacement matrix evaluated at the center of an 8-node
      ! trilinear hexahedral element.
      !
      ! The derivatives below correspond to the element-center point of the
      ! unit voxel. Rows 4-6 of es produce engineering shear strains.
      ! ======================================================================
      exx = e3d%exx; eyy = e3d%eyy; ezz = e3d%ezz
      exz = e3d%exz; eyz = e3d%eyz; exy = e3d%exy

      nxy = nx * ny

      ! Shape-function derivatives at the element center.
      dndx = [-0.25d0, 0.25d0, 0.25d0, -0.25d0, -0.25d0, 0.25d0, 0.25d0, -0.25d0]
      dndy = [-0.25d0, -0.25d0, 0.25d0, 0.25d0, -0.25d0, -0.25d0, 0.25d0, 0.25d0]
      dndz = [-0.25d0, -0.25d0, -0.25d0, -0.25d0, 0.25d0, 0.25d0, 0.25d0, 0.25d0]

      ! Build strain-displacement matrix B.
      !
      ! Voigt strain ordering:
      !   1 = eps_xx
      !   2 = eps_yy
      !   3 = eps_zz
      !   4 = gamma_xz
      !   5 = gamma_yz
      !   6 = gamma_xy
      es = 0.d0
      do n = 1, nnode_fe
        es(1,n,1) = dndx(n)
        es(2,n,2) = dndy(n)
        es(3,n,3) = dndz(n)
		
        es(4,n,1) = dndz(n);   es(4,n,3) = dndx(n); !  xz
        es(5,n,2) = dndz(n);   es(5,n,3) = dndy(n); !  yz
        es(6,n,1) = dndy(n);   es(6,n,2) = dndx(n); !  xy
      end do

      ! Initialize global stress and strain accumulators.
      ! These are summed over all non-pore sites and normalized after the loop.
      strxx = 0.d0; stryy = 0.d0; strzz = 0.d0
      strxz = 0.d0; stryz = 0.d0; strxy = 0.d0
      sxx = 0.d0;   syy = 0.d0;   szz = 0.d0
      sxz = 0.d0;   syz = 0.d0;   sxy = 0.d0
      ! ======================================================================
      ! ---- 2. Parallel Loop over Grid Sites ---------------------------------
      ! For each grid site/element:
      !   1. Gather the eight local nodal displacements.
      !   2. Apply periodic macroscopic displacement jumps on upper boundaries.
      !   3. Skip pore/void sites.
      !   4. Compute local strain at the element center.
      !   5. Compute local stress using the phase stiffness tensor.
      !   6. Accumulate stress and strain into OpenMP reduction variables.
      ! ======================================================================
      ! The reduction clause accumulates global stress/strain sums safely across
      ! threads. Local arrays such as uu and ibl are private to each thread.	  
      !$omp parallel do collapse(3) default(shared) private(m_local,mm,n8,n3,n,ibl,uu,str11_t,str22_t,str33_t,str13_t,str23_t,str12_t, s11_t,s22_t,s33_t,s13_t,s23_t,s12_t, pixm) reduction(+:strxx,stryy,strzz,strxz,stryz,strxy,sxx,syy,szz,sxz,syz,sxy) schedule(guided)
      do k_local = 1, nz
        do j_local = 1, ny
          do i_local = 1, nx
		  
            m_local = (k_local-1)*nxy + (j_local-1)*nx + i_local

            ! Gather only the neighbor-table entries needed to form the eight
            ! local nodes of the hexahedral element associated with m_local.
            ibl( 1) = e3d%ib(m_local,  1)
			ibl( 2) = e3d%ib(m_local,  2)
            ibl( 3) = e3d%ib(m_local,  3)
            ibl(17) = e3d%ib(m_local, 17)
			ibl(18) = e3d%ib(m_local, 18)
            ibl(19) = e3d%ib(m_local, 19)
			ibl(26) = e3d%ib(m_local, 26)

            ! Load local nodal displacements into element-node ordering.
            ! The neighbor indices correspond to the element connectivity used
            ! consistently throughout the stiffness and stencil routines.
            uu = 0.d0
            do mm = 1, ndof
              uu(1,mm) = e3d%u(m_local,mm)
              uu(2,mm) = e3d%u(ibl(3 ),mm)
              uu(3,mm) = e3d%u(ibl(2 ),mm)
              uu(4,mm) = e3d%u(ibl(1 ),mm)
              uu(5,mm) = e3d%u(ibl(26),mm)
              uu(6,mm) = e3d%u(ibl(19),mm)
              uu(7,mm) = e3d%u(ibl(18),mm)
              uu(8,mm) = e3d%u(ibl(17),mm)
            end do

            ! Apply periodic macroscopic displacement jumps for elements that
            ! cross the upper x, y, or z periodic boundaries.
            !
            ! The stored shear strain components exz, eyz, and exy are tensorial
            ! shear components. Because displacements are applied symmetrically,
            ! the resulting B*u shear strains are engineering shear strains.
            if(i_local == nx) then
              uu(2,1) = uu(2,1) + exx * dble(nx)
              uu(3,1) = uu(3,1) + exx * dble(nx)
              uu(6,1) = uu(6,1) + exx * dble(nx)
              uu(7,1) = uu(7,1) + exx * dble(nx)
              uu(2,2) = uu(2,2) + exy * dble(nx)
              uu(3,2) = uu(3,2) + exy * dble(nx)
              uu(6,2) = uu(6,2) + exy * dble(nx)
              uu(7,2) = uu(7,2) + exy * dble(nx)
              uu(2,3) = uu(2,3) + exz * dble(nx)
              uu(3,3) = uu(3,3) + exz * dble(nx)
              uu(6,3) = uu(6,3) + exz * dble(nx)
              uu(7,3) = uu(7,3) + exz * dble(nx)
            end if
            if(j_local == ny) then
              uu(3,1) = uu(3,1) + exy * dble(ny)
              uu(4,1) = uu(4,1) + exy * dble(ny)
              uu(7,1) = uu(7,1) + exy * dble(ny)
              uu(8,1) = uu(8,1) + exy * dble(ny)
              uu(3,2) = uu(3,2) + eyy * dble(ny)
              uu(4,2) = uu(4,2) + eyy * dble(ny)
              uu(7,2) = uu(7,2) + eyy * dble(ny)
              uu(8,2) = uu(8,2) + eyy * dble(ny)
              uu(3,3) = uu(3,3) + eyz * dble(ny)
              uu(4,3) = uu(4,3) + eyz * dble(ny)
              uu(7,3) = uu(7,3) + eyz * dble(ny)
              uu(8,3) = uu(8,3) + eyz * dble(ny)
            end if
            if(k_local == nz) then
              uu(5,1) = uu(5,1) + exz * dble(nz)
              uu(6,1) = uu(6,1) + exz * dble(nz)
              uu(7,1) = uu(7,1) + exz * dble(nz)
              uu(8,1) = uu(8,1) + exz * dble(nz)
              uu(5,2) = uu(5,2) + eyz * dble(nz)
              uu(6,2) = uu(6,2) + eyz * dble(nz)
              uu(7,2) = uu(7,2) + eyz * dble(nz)
              uu(8,2) = uu(8,2) + eyz * dble(nz)
              uu(5,3) = uu(5,3) + ezz * dble(nz)
              uu(6,3) = uu(6,3) + ezz * dble(nz)
              uu(7,3) = uu(7,3) + ezz * dble(nz)
              uu(8,3) = uu(8,3) + ezz * dble(nz)
            end if

            ! Pore/void sites have zero stiffness and are excluded from the
            ! stress/strain sums. Since final normalization uses ns, pores are
            ! included in the total-domain average as zero-contribution regions.
            if (e3d%pix(m_local) == PORE_PHASE) then
              ! Pore/void: zero contribution to volume averages.
              cycle
            end if

            ! Reset local stress and strain accumulators for this grid site.
            str11_t = 0.d0; str22_t = 0.d0; str33_t = 0.d0
            str13_t = 0.d0; str23_t = 0.d0; str12_t = 0.d0
            s11_t = 0.d0;   s22_t = 0.d0;   s33_t = 0.d0
            s13_t = 0.d0;   s23_t = 0.d0;   s12_t = 0.d0

            pixm = e3d%pix(m_local)
			
            ! Compute local strain eps = B*u and stress sigma = C*eps.
            !
            ! The strain accumulators s11_t, ..., s12_t store the six Voigt
            ! strain components at the element center. The stress accumulators
            ! str11_t, ..., str12_t store the corresponding stress components.
            !$omp simd
            do n3 = 1, ndof
              do n8 = 1, nnode_fe
                s11_t = s11_t + es(1,n8,n3)*uu(n8,n3)
                s22_t = s22_t + es(2,n8,n3)*uu(n8,n3)
                s33_t = s33_t + es(3,n8,n3)*uu(n8,n3)
                s13_t = s13_t + es(4,n8,n3)*uu(n8,n3)
                s23_t = s23_t + es(5,n8,n3)*uu(n8,n3)
                s12_t = s12_t + es(6,n8,n3)*uu(n8,n3)
                do n = 1, 6
                  str11_t = str11_t + e3d%cmod(pixm,1,n)*es(n,n8,n3)*uu(n8,n3)
                  str22_t = str22_t + e3d%cmod(pixm,2,n)*es(n,n8,n3)*uu(n8,n3)
                  str33_t = str33_t + e3d%cmod(pixm,3,n)*es(n,n8,n3)*uu(n8,n3)
                  str13_t = str13_t + e3d%cmod(pixm,4,n)*es(n,n8,n3)*uu(n8,n3)
                  str23_t = str23_t + e3d%cmod(pixm,5,n)*es(n,n8,n3)*uu(n8,n3)
                  str12_t = str12_t + e3d%cmod(pixm,6,n)*es(n,n8,n3)*uu(n8,n3)
                end do
              end do
            end do

            ! Accumulate local stress and strain into global reduction sums.
            strxx = strxx + str11_t
            stryy = stryy + str22_t
            strzz = strzz + str33_t
            strxz = strxz + str13_t
            stryz = stryz + str23_t
            strxy = strxy + str12_t
            sxx   = sxx   + s11_t
            syy   = syy   + s22_t
            szz   = szz   + s33_t
            sxz   = sxz   + s13_t
            syz   = syz   + s23_t
            sxy   = sxy   + s12_t

          end do
        end do
      end do
      !$omp end parallel do

      ! ======================================================================
      ! ---- 3. Volume Averaging and Storage ----------------------------------
      ! Normalize accumulated sums by the total number of grid sites, ns, and
      ! store the resulting volume-averaged stress and strain components in e3d.
      !
      ! Because pore/void sites were skipped in the loop but ns is used in the
      ! denominator, pores are treated as zero-stress/zero-strain volume in the
      ! reported domain averages.
      !
      ! If a solid-only average is desired, replace ns with the number of
      ! non-pore sites used in the accumulation.
      ! ======================================================================
      strxx = strxx / dble(ns)
      stryy = stryy / dble(ns)
      strzz = strzz / dble(ns)
      strxz = strxz / dble(ns)
      stryz = stryz / dble(ns)
      strxy = strxy / dble(ns)
	  
      sxx = sxx / dble(ns)
      syy = syy / dble(ns)
      szz = szz / dble(ns)
      sxz = sxz / dble(ns)
      syz = syz / dble(ns)
      sxy = sxy / dble(ns)

      ! Store volume-averaged stress components.
      e3d%strxx = strxx
      e3d%stryy = stryy
      e3d%strzz = strzz
      e3d%strxz = strxz
      e3d%stryz = stryz
      e3d%strxy = strxy
      ! Store volume-averaged strain components.	  
      e3d%sxx   = sxx
      e3d%syy   = syy
      e3d%szz   = szz
      e3d%sxz   = sxz
      e3d%syz   = syz
      e3d%sxy   = sxy

   
      
      end subroutine stress_OpenMP
!-------------------------------------------------------------------------     
      subroutine assig(e3d, prob)
      ! ==============================================================================
      ! SUBROUTINE: assig
      ! ------------------------------------------------------------------------------
      ! PURPOSE:
      !   Compute the grid-site fraction of each phase in the phase-label array
      !   e3d%pix.
      !
      !   The routine counts how many entries of e3d%pix are assigned to each
      !   valid phase ID and normalizes the counts by ns, the total number of
      !   grid sites.
      !
      ! ARGUMENTS:
      !   e3d
      !     Main solver data structure.
      !
      !     Input component:
      !       e3d%pix
      !         Phase-label array of length ns.
      !
      !   prob
      !     Output phase-fraction array of length nphase.
      !
      !     prob(i) is the fraction of grid sites with:
      !
      !       e3d%pix(m) == i
      !
      !     for phase i.
      !
      ! PHASE CONVENTION:
      !   Valid phase IDs are:
      !
      !      1, 2, ..., nphase
      !
      !   The pore/void phase is:
      !
      !      PORE_PHASE = nphase
      !
      ! AVERAGING CONVENTION:
      !   Fractions are normalized by ns. Because the grid is regular, these
      !   fractions correspond to volume fractions when each grid site represents
      !   an equal-volume cell/site.
      !
      ! INVALID LABEL HANDLING:
      !   Entries with phase labels outside [1,nphase] are ignored by this
      !   routine. The input reader ppixel_hdf5 performs a stricter global
      !   bounds check after reading the phase map.
      ! ==============================================================================  
      use elas3d_mod
      implicit none
      
      ! Main solver state containing the phase-label field.
      type(elas3d_data_type), intent(in) :: e3d
	  
      ! Phase fractions for phases 1:nphase.	  
      real(dp), intent(out) :: prob(nphase)
      
      ! Linear grid-site index.
      integer :: m

      ! ======================================================================
      ! ---- 1. Initialize Phase Counts ---------------------------------------
      ! Start from zero counts/fractions for every phase.
      ! ======================================================================
      prob(:) = 0.d0

      ! ======================================================================
      ! ---- 2. Count Phase Labels --------------------------------------------
      ! Traverse the phase-label array and increment the count for each valid
      ! phase ID. Out-of-range labels are ignored here.
      ! ======================================================================
      do m = 1, ns
        if(e3d%pix(m) >= 1 .and. e3d%pix(m) <= nphase) then
          prob(e3d%pix(m)) = prob(e3d%pix(m)) + 1.d0
        end if
      end do

      ! ======================================================================
      ! ---- 3. Normalize Counts ----------------------------------------------
      ! Convert counts to fractions by dividing by the total number of grid sites.
      ! ======================================================================
      prob(:) = prob(:) / dble(ns)


      end subroutine assig	  
!-------------------------------------------------------------------------
      subroutine rotatestiffnessby(e3d)
      ! ==============================================================================
      ! SUBROUTINE: rotatestiffnessby
      ! ------------------------------------------------------------------------------
      ! PURPOSE:
      !   Rotate the cubic single-crystal stiffness tensor into the global/sample
      !   coordinate frame for each solid grain phase.
      !
      !   Each grain phase has an associated Rodrigues orientation vector stored in
      !   e3d%orientation. This routine converts that Rodrigues vector into a 3 x 3
      !   rotation matrix, builds the corresponding 6 x 6 Bond transformation matrix
      !   in Voigt notation, and computes the rotated stiffness tensor.
      !
      ! ARGUMENTS:
      !   e3d
      !     Main solver data structure.
      !
      !     Input component:
      !       e3d%orientation
      !         Rodrigues orientation vectors for the solid grain phases.
      !
      !     Output component:
      !       e3d%rotatedstiffness
      !         Flattened rotated 6 x 6 stiffness tensor for each solid grain.
      !
      ! PHASE CONVENTION:
      !   Rotates solid grain phases 1:(nphase-1). The final phase, nphase, is
      !   the pore/void phase and is not handled here.
      !
      ! VOIGT ORDERING:
      !   1=xx, 2=yy, 3=zz, 4=xz, 5=yz, 6=xy
      !
      ! TRANSFORMATION:
      !   A rotation matrix Q is constructed from each Rodrigues vector. A 6 x 6
      !   Bond transformation matrix newr is then formed, and the rotated stiffness
      !   is computed as:
      !
      !      C_global = newr^T * C_crystal * newr
      !
      ! NOTES:
      !   The orientation convention must match the preprocessing code that writes
      !   the Rodrigues vectors.
      ! ==============================================================================
      use elas3d_mod
      implicit none

      ! Main solver state containing orientation input and rotated stiffness output.
      type(elas3d_data_type), intent(inout) :: e3d
      
      ! Loop indices.
      integer :: i, j, k, l, g
	  
      ! Rodrigues-vector quantities.
      real(dp) :: term1, term2, r(3), rdotr
	  
      ! Cubic stiffness tensor in the local crystal frame.
	  real(dp) :: StiffnessInCrystalFrame(6,6)
	  
      ! Rotation matrix and Bond transformation matrix.
	  real(dp) :: OrientationMatrix(3,3)
      real(dp) :: newr(6,6)
	  
      ! Rotated stiffness tensor for the current grain.
      real(dp) :: rotated(6,6)
	  
      ! Alias for the rotation matrix used in Bond matrix construction.
	  real(dp) :: Q(3,3)

	  
      ! ======================================================================
      ! ---- 1. Initialize Cubic Stiffness in the Crystal Frame ---------------
      ! Assemble the unrotated single-crystal stiffness matrix using cubic
      ! elastic constants.
      ! ======================================================================
      StiffnessInCrystalFrame(:,:) = 0.d0
      StiffnessInCrystalFrame(1,1) = C11_local
      StiffnessInCrystalFrame(2,2) = C11_local
      StiffnessInCrystalFrame(3,3) = C11_local
      StiffnessInCrystalFrame(1,2) = C12_local
      StiffnessInCrystalFrame(1,3) = C12_local
      StiffnessInCrystalFrame(2,3) = C12_local
      StiffnessInCrystalFrame(2,1) = C12_local
      StiffnessInCrystalFrame(3,1) = C12_local
      StiffnessInCrystalFrame(3,2) = C12_local
      StiffnessInCrystalFrame(4,4) = C44_local
      StiffnessInCrystalFrame(5,5) = C44_local
      StiffnessInCrystalFrame(6,6) = C44_local

      ! ======================================================================
      ! ---- 2. Loop Over Solid Grain Phases ----------------------------------
      ! Rotate the crystal-frame stiffness tensor for each solid grain phase.
      ! ======================================================================
      do g = 1, nphase-1
	  
        ! Load Rodrigues orientation vector for grain phase g.
        r(1:3) = e3d%orientation((g-1)*3+1:(g-1)*3+3)    
		              
        ! Build the 3 x 3 rotation matrix from the Rodrigues vector.
        rdotr = sum(r(:)**2)
        term1 = 1.d0 - rdotr
        term2 = 1.d0 + rdotr
		
        OrientationMatrix(:,:) = 0.d0
        OrientationMatrix(1,1) = term1
        OrientationMatrix(2,2) = term1
        OrientationMatrix(3,3) = term1
        
        ! Rodrigues/Cayley formula: diagonal and symmetric r*r^T contribution.
        OrientationMatrix = OrientationMatrix + 2.d0 * spread(r, 2, 3) * spread(r, 1, 3)
       
        ! Add skew-symmetric terms. The signs define the rotation convention.
        OrientationMatrix(1,2) = OrientationMatrix(1,2) - 2.d0*r(3)
        OrientationMatrix(1,3) = OrientationMatrix(1,3) + 2.d0*r(2)
        OrientationMatrix(2,3) = OrientationMatrix(2,3) - 2.d0*r(1)
        OrientationMatrix(2,1) = OrientationMatrix(2,1) + 2.d0*r(3)
        OrientationMatrix(3,1) = OrientationMatrix(3,1) - 2.d0*r(2)
        OrientationMatrix(3,2) = OrientationMatrix(3,2) + 2.d0*r(1)
        
        ! Normalize by 1 + r.r. The fallback protects against pathological input.
        if (term2 > 1.0d-12) then
           term2 = 1.0d0 / term2
           OrientationMatrix(:,:) = OrientationMatrix(:,:) * term2
        else
           OrientationMatrix = 0.d0
           OrientationMatrix(1,1) = 1.0d0
           OrientationMatrix(2,2) = 1.0d0
           OrientationMatrix(3,3) = 1.0d0
        end if
		
        ! Short alias for readability.
        Q = OrientationMatrix
		
        ! ====================================================================
        ! ---- 3. Build 6 x 6 Bond Transformation Matrix ----------------------
        ! Construct the Voigt-space transformation matrix associated with Q.
        !
        ! Voigt ordering:
        !   [xx, yy, zz, xz, yz, xy]
        !
        ! Factors of 2 in the lower-left block account for the engineering
        ! shear-strain convention used in the finite-element B-matrix.
        ! ====================================================================
		
        ! Row 1: transformed xx component.
        newr(1,1) = Q(1,1)**2
        newr(1,2) = Q(1,2)**2
        newr(1,3) = Q(1,3)**2
        newr(1,4) = Q(1,2)*Q(1,3)
        newr(1,5) = Q(1,3)*Q(1,1)
        newr(1,6) = Q(1,1)*Q(1,2)
    
        ! Row 2: transformed yy component.
        newr(2,1) = Q(2,1)**2
        newr(2,2) = Q(2,2)**2
        newr(2,3) = Q(2,3)**2
        newr(2,4) = Q(2,2)*Q(2,3)
        newr(2,5) = Q(2,3)*Q(2,1)
        newr(2,6) = Q(2,1)*Q(2,2)
    
        ! Row 3: transformed zz component.
        newr(3,1) = Q(3,1)**2
        newr(3,2) = Q(3,2)**2
        newr(3,3) = Q(3,3)**2
        newr(3,4) = Q(3,2)*Q(3,3)
        newr(3,5) = Q(3,3)*Q(3,1)
        newr(3,6) = Q(3,1)*Q(3,2)
    
        ! Row 4: transformed xz shear component.
        newr(4,1) = 2.d0 * Q(3,1)*Q(1,1)
        newr(4,2) = 2.d0 * Q(3,2)*Q(1,2)
        newr(4,3) = 2.d0 * Q(3,3)*Q(1,3)
        newr(4,4) = Q(1,2)*Q(3,3) + Q(1,3)*Q(3,2)
        newr(4,5) = Q(1,3)*Q(3,1) + Q(1,1)*Q(3,3)
        newr(4,6) = Q(1,1)*Q(3,2) + Q(1,2)*Q(3,1)
    
        ! Row 5: transformed yz shear component.
        newr(5,1) = 2.d0 * Q(2,1)*Q(3,1)
        newr(5,2) = 2.d0 * Q(2,2)*Q(3,2)
        newr(5,3) = 2.d0 * Q(2,3)*Q(3,3)
        newr(5,4) = Q(2,2)*Q(3,3) + Q(2,3)*Q(3,2)
        newr(5,5) = Q(2,1)*Q(3,3) + Q(2,3)*Q(3,1)
        newr(5,6) = Q(2,2)*Q(3,1) + Q(2,1)*Q(3,2)
    
        ! Row 6: transformed xy shear component.
        newr(6,1) = 2.d0 * Q(1,1)*Q(2,1)
        newr(6,2) = 2.d0 * Q(1,2)*Q(2,2)
        newr(6,3) = 2.d0 * Q(1,3)*Q(2,3)
        newr(6,4) = Q(1,2)*Q(2,3) + Q(1,3)*Q(2,2)
        newr(6,5) = Q(1,3)*Q(2,1) + Q(1,1)*Q(2,3)
        newr(6,6) = Q(1,1)*Q(2,2) + Q(1,2)*Q(2,1)      

        ! ====================================================================
        ! ---- 4. Rotate Stiffness Tensor and Store ---------------------------
        ! Rotate the crystal-frame stiffness tensor into the global/sample frame:
        !
        !   C_global = newr^T * C_crystal * newr
        !
        ! Store the result in flattened grain-block order.
        ! ====================================================================
        rotated = matmul(matmul(transpose(newr), StiffnessInCrystalFrame), newr)	

		do i = 1, 6
		  do l = 1, 6
            e3d%rotatedstiffness((g-1)*36 + (i-1)*6 + l) = rotated(i, l)
          end do
        end do
		
      end do
      
      end subroutine rotatestiffnessby	  
!-------------------------------------------------------------------------
      subroutine stress_fullfield_OpenMP(e3d) 
      ! ==============================================================================
      ! SUBROUTINE: stress_fullfield_OpenMP
      ! ------------------------------------------------------------------------------
      ! PURPOSE:
      !   Compute and export the full-field stress tensor and von Mises equivalent
      !   stress over the entire simulation grid.
      !
      !   This routine post-processes the current displacement field e3d%u after
      !   convergence. For each grid site/element, it gathers the eight local nodal
      !   displacements, applies periodic macroscopic displacement jumps, computes
      !   the local strain at the element center, evaluates stress using the
      !   phase stiffness tensor, and stores the result in e3d%stress_field.
      !
      !   The von Mises equivalent stress is then computed from the stored stress
      !   tensor field and written, together with phase labels and full stress
      !   tensors, to an HDF5 output file.
      !
      ! ARGUMENTS:
      !   e3d
      !     Main solver data structure.
      !
      !     Input components:
      !       e3d%u
      !         Current nodal displacement field.
      !
      !       e3d%cmod
      !         6 x 6 constitutive stiffness matrix for each phase.
      !
      !       e3d%pix
      !         Phase-label array over the simulation grid.
      !
      !       e3d%ib
      !         Periodic 27-neighbor connectivity table.
      !
      !       e3d%exx, e3d%eyy, e3d%ezz, e3d%exz, e3d%eyz, e3d%exy
      !         Prescribed macroscopic strain components used to reconstruct
      !         periodic displacement jumps.
      !
      !     Output components:
      !       e3d%stress_field(ns,6)
      !         Full-field stress tensor in Voigt order:
      !
      !           [sxx, syy, szz, sxz, syz, sxy]
      !
      !       e3d%vm(ns)
      !         von Mises equivalent stress at each grid site.
      !
      ! VOIGT ORDERING:
      !   Stress:
      !     1 = sxx
      !     2 = syy
      !     3 = szz
      !     4 = sxz
      !     5 = syz
      !     6 = sxy
      !
      !   Strain:
      !     1 = eps_xx
      !     2 = eps_yy
      !     3 = eps_zz
      !     4 = gamma_xz = du_x/dz + du_z/dx
      !     5 = gamma_yz = du_y/dz + du_z/dy
      !     6 = gamma_xy = du_x/dy + du_y/dx
      !
      !   The shear strain components are engineering shear strains. The shear
      !   stress components are tensor stress components.
      !
      ! PORE HANDLING:
      !   Sites with phase ID PORE_PHASE are assigned:
      !
      !      stress_field(m,:) = 0
      !      vm(m)             = 0
      !
      !   and are skipped from the constitutive stress calculation.
      !
      ! OUTPUT FILE:
      !   fullfield_poly.h5
      !
      !   Datasets written:
      !     /pix
      !       Phase-label array, length ns.
      !
      !     /vm
      !       von Mises equivalent stress array, length ns.
      !
      !     /stress_tensor
      !       Full stress tensor array, dimension ns x 6.
      !
      ! NUMERICAL METHOD:
      !   - Uses 8-node trilinear hexahedral element kinematics.
      !   - Evaluates strain at the element/grid-site center.
      !   - Applies periodic macroscopic displacement jumps before strain evaluation.
      !   - Uses OpenMP for the full-field stress calculation and von Mises loop.
      ! ==============================================================================
      use elas3d_mod
      use hdf5
      implicit none
  
      ! Main solver state structure.
      type(elas3d_data_type), intent(inout) :: e3d

      ! Loop indices and grid-site indexing variables.
      integer :: i, j, k, m, nxy, mm, n8, n3, n
	  integer :: m_l, k_l, j_l, i_l 
	  
      ! Shape-function derivatives at the element center.
      real(dp) :: dndx(nnode_fe), dndy(nnode_fe), dndz(nnode_fe)
	  
      ! Strain-displacement matrix:
      !   es(strain_component, local_node, displacement_component)
      real(dp) :: es(6,nnode_fe,ndof)         

      ! Local element nodal displacement array.
      real(dp) :: uu(nnode_fe,ndof)
	  
      ! Local strain vector at one grid site.
      ! Voigt order:
      !   [eps_xx, eps_yy, eps_zz, gamma_xz, gamma_yz, gamma_xy]
      real(dp) :: eps(6)

      ! Periodic neighbor indices needed to build the local element connectivity.
      integer :: ibl(27)	  

      ! Prescribed macroscopic strain components.
      real(dp) :: exx, eyy, ezz, exz, eyz, exy  

      ! HDF5 identifiers and dimension arrays.
      integer(HID_T)    :: file_id, dset_id, space_id
      integer(HSIZE_T)  :: dims(2)
      integer           :: hdferr   

      ! ======================================================================
      ! ---- 1. Initialization and B-Matrix Setup -----------------------------
      ! Unpack prescribed macroscopic strain components and construct the
      ! strain-displacement matrix evaluated at the center of an 8-node
      ! trilinear hexahedral element.
      !
      ! Rows 1-3 of es compute normal strain components.
      ! Rows 4-6 compute engineering shear strain components.
      ! ======================================================================
      exx = e3d%exx; eyy = e3d%eyy; ezz = e3d%ezz
      exz = e3d%exz; eyz = e3d%eyz; exy = e3d%exy  
  
      ! Number of grid sites in one x-y plane, used for 3D-to-1D indexing.
      nxy = nx * ny  

      ! Trilinear shape-function derivatives at the element center.
      dndx = [-0.25d0, 0.25d0, 0.25d0,-0.25d0,-0.25d0, 0.25d0, 0.25d0,-0.25d0]
      dndy = [-0.25d0,-0.25d0, 0.25d0, 0.25d0,-0.25d0,-0.25d0, 0.25d0, 0.25d0]
      dndz = [-0.25d0,-0.25d0,-0.25d0,-0.25d0, 0.25d0, 0.25d0, 0.25d0, 0.25d0]
      
      ! Build strain-displacement matrix B.
      !
      ! Voigt strain ordering:
      !   1 = eps_xx
      !   2 = eps_yy
      !   3 = eps_zz
      !   4 = gamma_xz
      !   5 = gamma_yz
      !   6 = gamma_xy
      es = 0.d0
      do n = 1,8
        es(1,n,1) = dndx(n)
        es(2,n,2) = dndy(n)
        es(3,n,3) = dndz(n)
        es(4,n,1) = dndz(n);    es(4,n,3) = dndx(n)
        es(5,n,2) = dndz(n);    es(5,n,3) = dndy(n)
        es(6,n,1) = dndy(n);    es(6,n,2) = dndx(n)
      end do


      ! ======================================================================
      ! ---- 2. Parallel Full-Field Stress Calculation ------------------------
      ! For each grid site/element:
      !   1. Gather the eight local nodal displacements.
      !   2. Apply periodic macroscopic displacement jumps on upper boundaries.
      !   3. Compute the local strain vector eps = B*u.
      !   4. If the site is a pore, assign zero stress and skip Hooke's law.
      !   5. Otherwise compute stress = C*eps and store it in stress_field.
      !
      ! The loop is parallelized over the full 3D grid.
      ! ======================================================================
      !$omp parallel do collapse(3) schedule(guided) default(shared) private(k_l,j_l,i_l,m_l,mm,n8,n3,ibl,uu,eps) 
      do k_l = 1, nz
	    do j_l = 1, ny
		  do i_l = 1, nx
        
		    m_l = (k_l-1)*nxy + (j_l-1)*nx + i_l


            ! Gather the neighbor-table entries needed to form the eight local
            ! nodes of the hexahedral element associated with m_l.
            ibl( 1) = e3d%ib(m_l,  1); 
		    ibl( 2) = e3d%ib(m_l,  2); 
		    ibl( 3) = e3d%ib(m_l,  3); 
            ibl(17) = e3d%ib(m_l, 17); 
		    ibl(18) = e3d%ib(m_l, 18); 
		    ibl(19) = e3d%ib(m_l, 19); 
		    ibl(26) = e3d%ib(m_l, 26);
        
            ! Load local nodal displacements into the element-node ordering used
            ! by the finite-element stiffness and strain calculations.
            uu = 0.d0
            do mm = 1, ndof
              uu(1,mm) = e3d%u(m_l,    mm)
              uu(2,mm) = e3d%u(ibl(3 ),mm)
              uu(3,mm) = e3d%u(ibl(2 ),mm)
              uu(4,mm) = e3d%u(ibl(1 ),mm)
              uu(5,mm) = e3d%u(ibl(26),mm)
              uu(6,mm) = e3d%u(ibl(19),mm)
              uu(7,mm) = e3d%u(ibl(18),mm)
              uu(8,mm) = e3d%u(ibl(17),mm)
            end do
    
            ! Apply periodic macroscopic displacement jumps for elements crossing
            ! the upper x, y, or z periodic boundaries.
            !
            ! The stored shear strain components exz, eyz, and exy are tensorial
            ! shear strain components used to construct the displacement field.
            if (i_l == nx) then
              uu(2,1) = uu(2,1) + exx * dble(nx)
              uu(3,1) = uu(3,1) + exx * dble(nx)
              uu(6,1) = uu(6,1) + exx * dble(nx)
              uu(7,1) = uu(7,1) + exx * dble(nx)
              uu(2,2) = uu(2,2) + exy * dble(nx)
              uu(3,2) = uu(3,2) + exy * dble(nx)
              uu(6,2) = uu(6,2) + exy * dble(nx)
              uu(7,2) = uu(7,2) + exy * dble(nx)
              uu(2,3) = uu(2,3) + exz * dble(nx)
              uu(3,3) = uu(3,3) + exz * dble(nx)
              uu(6,3) = uu(6,3) + exz * dble(nx)
              uu(7,3) = uu(7,3) + exz * dble(nx)
            end if
            if (j_l == ny) then
              uu(3,1) = uu(3,1) + exy * dble(ny)
              uu(4,1) = uu(4,1) + exy * dble(ny)
              uu(7,1) = uu(7,1) + exy * dble(ny)
              uu(8,1) = uu(8,1) + exy * dble(ny)
              uu(3,2) = uu(3,2) + eyy * dble(ny)
              uu(4,2) = uu(4,2) + eyy * dble(ny)
              uu(7,2) = uu(7,2) + eyy * dble(ny)
              uu(8,2) = uu(8,2) + eyy * dble(ny)
              uu(3,3) = uu(3,3) + eyz * dble(ny)
              uu(4,3) = uu(4,3) + eyz * dble(ny)
              uu(7,3) = uu(7,3) + eyz * dble(ny)
              uu(8,3) = uu(8,3) + eyz * dble(ny)
            end if
            if (k_l == nz) then
              uu(5,1) = uu(5,1) + exz * dble(nz)
              uu(6,1) = uu(6,1) + exz * dble(nz)
              uu(7,1) = uu(7,1) + exz * dble(nz)
              uu(8,1) = uu(8,1) + exz * dble(nz)
              uu(5,2) = uu(5,2) + eyz * dble(nz)
              uu(6,2) = uu(6,2) + eyz * dble(nz)
              uu(7,2) = uu(7,2) + eyz * dble(nz)
              uu(8,2) = uu(8,2) + eyz * dble(nz)
              uu(5,3) = uu(5,3) + ezz * dble(nz)
              uu(6,3) = uu(6,3) + ezz * dble(nz)
              uu(7,3) = uu(7,3) + ezz * dble(nz)
              uu(8,3) = uu(8,3) + ezz * dble(nz)
            end if         

            ! Compute local strain vector eps = B*u at the element center.
            !
            ! Voigt strain order:
            !   [eps_xx, eps_yy, eps_zz, gamma_xz, gamma_yz, gamma_xy]
            !
            ! Shear entries are engineering shear strains.
            eps = 0.d0
            do n3 = 1, ndof
              do n8 = 1, nnode_fe
                eps(1) = eps(1) + es(1,n8,n3)*uu(n8,n3)
                eps(2) = eps(2) + es(2,n8,n3)*uu(n8,n3)
                eps(3) = eps(3) + es(3,n8,n3)*uu(n8,n3)
                eps(4) = eps(4) + es(4,n8,n3)*uu(n8,n3)
                eps(5) = eps(5) + es(5,n8,n3)*uu(n8,n3)
                eps(6) = eps(6) + es(6,n8,n3)*uu(n8,n3)
              end do
            end do

            ! Pore/void sites have zero stiffness. Store zero stress and zero
            ! von Mises stress, then skip the constitutive calculation.
            if (e3d%pix(m_l) == PORE_PHASE) then
              e3d%stress_field(m_l,:) = 0.d0
              e3d%vm(m_l)             = 0.d0
              cycle
            end if

            ! Compute and store local stress:
            !
            !   sigma = C_phase * eps
            !
            ! Stress is stored in Voigt order:
            !   [sxx, syy, szz, sxz, syz, sxy]
            e3d%stress_field(m_l,:) = matmul(e3d%cmod(e3d%pix(m_l),1:6,1:6), eps)
        
          end do
		end do
      end do
      !$omp end parallel do

      
      ! ======================================================================
      ! ---- 3. von Mises Equivalent Stress Calculation -----------------------
      ! Compute von Mises equivalent stress from the stored stress tensor:
      !
      !   sigma_vm = sqrt( 0.5*((sxx-syy)^2 + (syy-szz)^2 + (szz-sxx)^2
      !                    + 6*(sxz^2 + syz^2 + sxy^2)) )
      !
      ! Pore entries were already assigned zero stress, so their von Mises stress
      ! remains zero.
      ! ====================================================================== 
      !$omp parallel do schedule(guided)
      do m = 1, ns
        e3d%vm(m) = sqrt(0.5d0*((e3d%stress_field(m,1)-e3d%stress_field(m,2))**2 + &
                                (e3d%stress_field(m,2)-e3d%stress_field(m,3))**2 + &
                                (e3d%stress_field(m,3)-e3d%stress_field(m,1))**2 + &
                         6.0d0*(e3d%stress_field(m,4)**2 + &
                                e3d%stress_field(m,5)**2 + &
						        e3d%stress_field(m,6)**2)))
      end do
      !$omp end parallel do


      ! ======================================================================
      ! ---- 4. HDF5 Data Export ---------------------------------------------
      ! Write phase labels, von Mises stress, and full stress tensors to an HDF5
      ! file for post-processing and visualization.
      !
      ! Output file:
      !   fullfield_poly.h5
      !
      ! Datasets:
      !   /pix
      !     Integer phase-label array of length ns.
      !
      !   /vm
      !     Double-precision von Mises stress array of length ns.
      !
      !   /stress_tensor
      !     Double-precision stress tensor array with dimensions ns x 6.
      !     Voigt order:
      !       [sxx, syy, szz, sxz, syz, sxy]
      ! ======================================================================
      call h5open_f(hdferr)
      ! Create or overwrite the full-field output file.
      call h5fcreate_f('fullfield_poly.h5', H5F_ACC_TRUNC_F, file_id, hdferr)
	  
      ! Create 1D dataspace for arrays of length ns.
      dims(1) = ns 
      call h5screate_simple_f(1, dims, space_id, hdferr)
	  
      ! Dataset /pix: phase labels.
      call h5dcreate_f(file_id, "pix", H5T_NATIVE_INTEGER, space_id, dset_id, hdferr)
	  call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, e3d%pix, dims, hdferr)
      call h5dclose_f(dset_id, hdferr)
	  
      ! Dataset /vm: von Mises equivalent stress.
      call h5dcreate_f(file_id, "vm", H5T_NATIVE_DOUBLE, space_id, dset_id, hdferr)
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, e3d%vm, dims, hdferr)
      call h5dclose_f(dset_id, hdferr)
      
      ! Close 1D dataspace.
      call h5sclose_f(space_id, hdferr)
	  
      ! Dataset /stress_tensor: full stress tensor field, dimension ns x 6.
      dims(1) = ns
      dims(2) = 6
      call h5screate_simple_f(2, dims, space_id, hdferr)

      call h5dcreate_f(file_id, "stress_tensor", H5T_NATIVE_DOUBLE, space_id, dset_id, hdferr)
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, e3d%stress_field, dims, hdferr)
      call h5dclose_f(dset_id, hdferr)
	  
      ! Close HDF5 dataset, dataspace, file, and library handles.
	  call h5sclose_f(space_id, hdferr)	  
      call h5fclose_f(file_id, hdferr)
      call h5close_f(hdferr)

      
      end subroutine stress_fullfield_OpenMP    
!-------------------------------------------------------------------------    
      subroutine ppixel_hdf5(e3d)
      ! ==============================================================================
      ! SUBROUTINE: ppixel_hdf5
      ! ------------------------------------------------------------------------------
      ! PURPOSE:
      !   Read the simulation microstructure and grain orientation data from the
      !   HDF5 input file:
      !
      !      input_structure_poly.h5
      !
      !   The routine populates:
      !
      !      e3d%pix
      !      e3d%orientation
      !
      !   These arrays are then used to assign phase-dependent stiffness tensors
      !   and construct the periodic finite-element problem.
      !
      ! ARGUMENTS:
      !   e3d
      !     Main solver data structure.
      !
      !     Output components:
      !       e3d%pix(ns)
      !         Phase-label array for the simulation grid.
      !
      !       e3d%orientation(3*(nphase-1))
      !         Rodrigues orientation vectors for the solid grain phases.
      !
      ! INPUT FILE FORMAT:
      !   File:
      !      input_structure_poly.h5
      !
      !   Required datasets:
      !
      !      /pix
      !        Integer phase-label dataset.
      !        Total number of entries must be ns = nx*ny*nz.
      !        Valid phase labels are 1:nphase, where nphase is the pore/void phase.
      !
      !      /orientation
      !        Double-precision Rodrigues-vector dataset for solid grain phases.
      !        Total number of entries must be 3*(nphase-1).
      !
      ! VALIDATION:
      !   The routine checks:
      !     - /pix total size equals ns,
      !     - /orientation total size equals size(e3d%orientation),
      !     - phase labels are within [1,nphase].
      !
      ! NOTES:
      !   - This routine checks total dataset size, not optional HDF5 attributes.
      !   - The /pix ordering must match the Fortran solver's expected indexing.
      !   - The orientation convention must be consistent with rotatestiffnessby.
      ! ==============================================================================
      use elas3d_mod
      use hdf5
      implicit none

      ! Main solver state structure to be populated from HDF5.
      type(elas3d_data_type), intent(inout) :: e3d

      ! HDF5 file, dataset, and dataspace identifiers.
      integer(HID_T)   :: file_id, dset_id, space_id
	  
      ! Dataset dimensions and maximum dimensions returned by HDF5.
      integer(HSIZE_T), allocatable :: dims(:), maxdims(:)
	  
      ! Dataset rank and HDF5 error/status code.
      integer           :: rank, hdferr
	  
      ! Total number of entries in the current dataset.
	  integer(HSIZE_T)  :: total_elems
	  
      ! One-dimensional dimension array used when reading orientation data.
      integer(HSIZE_T)  :: odims(1)

      ! ======================================================================
      ! ---- 1. Open HDF5 Library and Input File ------------------------------
      ! Initialize the HDF5 Fortran interface and open the input file.
      ! ======================================================================
      call h5open_f(hdferr)
      call h5fopen_f('input_structure_poly.h5', H5F_ACC_RDONLY_F, file_id, hdferr)
	  if (hdferr /= 0) stop 'Error opening HDF5 file'
	  
      ! ======================================================================
      ! ---- 2. Read /pix Dataset: Phase Labels -------------------------------
      ! Read the integer phase-label array. The total number of entries must
      ! match ns, the number of grid sites used by the solver.
      ! ======================================================================
      call h5dopen_f(file_id, 'pix', dset_id, hdferr)
      call h5dget_space_f(dset_id, space_id, hdferr)

      ! Query dataset rank and dimensions.
      call h5sget_simple_extent_ndims_f(space_id, rank, hdferr)
      allocate(dims(rank), maxdims(rank))
      call h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdferr)


      ! Verify that the total number of phase labels matches the solver grid.
      ! Note: for very large datasets, care is needed because product(dims)
      ! may depend on integer kind behavior.
	  total_elems = product(dims)
	  
      if (total_elems /= ns) then
        print *, 'ERROR: HDF5 pix size (', total_elems, ') does not match simulation ns (', ns, ')'
        stop
      end if

      ! Read phase labels into the solver state.
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, e3d%pix, dims, hdferr)
	  
      ! Close /pix dataset handles and release dimension arrays.
      call h5dclose_f(dset_id, hdferr)
      call h5sclose_f(space_id, hdferr)
      deallocate(dims, maxdims)
	  
      ! ======================================================================
      ! ---- 3. Read /orientation Dataset: Rodrigues Vectors ------------------
      ! Read Rodrigues orientation vectors for the solid grain phases.
      ! Expected total number of entries is 3*(nphase-1), matching
      ! size(e3d%orientation).
      ! ======================================================================
      call h5dopen_f(file_id, 'orientation', dset_id, hdferr)
      call h5dget_space_f(dset_id, space_id, hdferr)
	  
      ! Query rank and dimensions of the orientation dataset.
      call h5sget_simple_extent_ndims_f(space_id, rank, hdferr)
      allocate(dims(rank), maxdims(rank))
      call h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdferr)

      ! Verify that the dataset contains exactly one Rodrigues vector per solid
      ! grain phase.
	  total_elems = product(dims)
      if (total_elems /= size(e3d%orientation)) then
          print *, 'ERROR: HDF5 orientation size (', total_elems, ') mismatch with allocated (', size(e3d%orientation), ')'
          stop
      end if
	  
      ! Read orientation data as a 1D Rodrigues-vector array.
      odims(1) = size(e3d%orientation)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, e3d%orientation, odims, hdferr)
	  
      ! Close /orientation dataset handles and release dimension arrays.
      call h5dclose_f(dset_id, hdferr)
      call h5sclose_f(space_id, hdferr)
      deallocate(dims, maxdims)


      ! Close HDF5 file and library handles.
      call h5fclose_f(file_id, hdferr)
      call h5close_f(hdferr)

      ! Final phase-label bounds check.
      ! All phase IDs must lie in [1,nphase]. Phase nphase is the pore/void phase.
      if (minval(e3d%pix) < 1 .or. maxval(e3d%pix) > nphase) then
        print *, "ERROR: Phase label in pix out of bounds"
        stop
      end if

      end subroutine ppixel_hdf5