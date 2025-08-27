! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_param
use mod_types
use mod_initgrid, only: pos_array
#if defined(_OPENACC)
use cudecomp
#endif
implicit none
public
!
! parameters
!
real(rp), parameter :: pi = acos(-1._rp)
#if !defined(_EPS_EXACT_ZERO) /* recommended */
real(rp), parameter :: eps = epsilon(1._rp)
#else
real(rp), parameter :: eps = 0._rp
#endif
real(rp), parameter :: small = epsilon(1._rp)*10**(precision(1._rp)/2)
character(len=100), parameter :: datadir = 'data/'
real(rp), parameter, dimension(2,3) :: rkcoeff = reshape([32._rp/60._rp,  0._rp        , &
                                                          25._rp/60._rp, -17._rp/60._rp, &
                                                          45._rp/60._rp, -25._rp/60._rp], shape(rkcoeff))
real(rp), parameter, dimension(3)   :: rkcoeff12 = rkcoeff(1,:)+rkcoeff(2,:)
type xyz
  real(rp), allocatable :: x(:)
  real(rp), allocatable :: y(:)
  real(rp), allocatable :: z(:)
end type
!
! variables to be determined from the input file
!
integer , dimension(3) :: ng ! not protected -> could be changed if non-uniform meshes are used
real(rp), protected, dimension(3) :: l
integer , protected, dimension(3) :: gtype
real(rp), protected, dimension(3) :: gr
logical , protected :: is_cell_length
real(rp), protected :: cfl,dtmax,dt_f
real(rp), protected :: visci
!
character(len=100), protected :: inivel
logical, protected :: is_wallturb
!
integer , protected :: nstep
real(rp), protected :: time_max,tw_max
logical , protected, dimension(3) :: stop_type
logical , protected :: restart,is_overwrite_save
integer , protected :: nsaves_max
integer , protected :: icheck,iout0d,iout1d,iout2d,iout3d,isave
!
integer , dimension(2) :: dims ! not protected -> overwritten when autotuning
integer , protected :: ipencil_axis = 1
!
integer, dimension(0:1,3) :: nb
logical, dimension(0:1,3) :: is_bound
character(len=1), protected, dimension(0:1,3,3) ::  cbcvel
real(rp)        , protected, dimension(0:1,3,3) ::   bcvel
character(len=1), protected, dimension(0:1,3)   ::  cbcpre
real(rp)        , protected, dimension(0:1,3)   ::   bcpre
!
real(rp), protected, dimension(3) :: bforce
logical , protected, dimension(3) :: is_forced
real(rp), protected, dimension(3) :: velf
!
real(rp), protected :: visc
logical, protected, dimension(3) :: is_user_mesh
character(len=1) :: dim_char
real(rp) :: val
integer, dimension(3) :: npoints
!
! scalar input parameters
!
integer , protected :: nscal ! number of transported scalars
real(rp), protected, allocatable, dimension(:) :: alphai
real(rp), protected :: beta
real(rp), protected, dimension(3) :: gacc
character(len=100), protected, allocatable, dimension(:)     :: iniscal
character(len=1)  , protected, allocatable, dimension(:,:,:) :: cbcscal ! size (0:1,3,nscal)
real(rp)          , protected, allocatable, dimension(:,:,:) ::  bcscal ! size (0:1,3,nscal)
real(rp), protected, allocatable, dimension(:) :: ssource
logical , protected, allocatable, dimension(:) :: is_sforced
real(rp), protected, allocatable, dimension(:) :: scalf
logical , protected :: is_boussinesq_buoyancy = .false.
real(rp), protected :: alpha_max
!
#if defined(_OPENACC)
!
! cuDecomp input parameters
!
integer, protected :: cudecomp_t_comm_backend,cudecomp_h_comm_backend
logical, protected :: cudecomp_is_t_comm_autotune ,cudecomp_is_h_comm_autotune , &
                      cudecomp_is_t_enable_nccl   ,cudecomp_is_h_enable_nccl   , &
                      cudecomp_is_t_enable_nvshmem,cudecomp_is_h_enable_nvshmem, &
                      cudecomp_is_t_in_place
#endif
!
! other options: debugging/benchmarking
!
logical, protected :: is_debug = .true., is_debug_poisson = .false., &
                      is_timing = .true., &
                      is_mask_divergence_check = .false.
!
! other options: numerics
!
logical, protected :: is_impdiff = .false., is_impdiff_1d = .false., &
                      is_poisson_pcr_tdma = .false., &
                      is_fast_mom_kernels = .true., &
                      is_gridpoint_natural_channel = .false.
contains
  subroutine read_input(myid, mesh_output, dl, dli)
    use, intrinsic :: iso_fortran_env, only: iostat_end
    use mpi
    implicit none
    character(len=*), parameter :: input_file = 'input.nml'
    character(len=*), parameter :: mesh_file = 'mesh.csv'
    integer, intent(in) :: myid
    type(xyz), allocatable, intent(out) :: mesh_output
    ! dl and dli are declared after reading the input file, but they are both intent(out) 
    integer :: iunit,ierr,iunit_mesh_file
    character(len=1024) :: c_iomsg
    integer :: i_dim
    integer :: x_counter = 1
    integer :: y_counter = 1
    integer :: z_counter = 1
    namelist /dns/ &
                  ng, &
                  l, &
                  gtype, &
                  gr, &
                  is_cell_length, &
                  cfl,dtmax,dt_f, &
                  visci, &
                  inivel, &
                  is_wallturb, &
                  nstep,time_max,tw_max, &
                  stop_type, &
                  restart,is_overwrite_save,nsaves_max, &
                  icheck,iout0d,iout1d,iout2d,iout3d,isave, &
                  cbcvel,cbcpre,bcvel,bcpre, &
                  bforce, &
                  is_forced, &
                  velf, &
                  gacc, &
                  nscal, &
                  dims,ipencil_axis
#if defined(_OPENACC)
    namelist /cudecomp/ &
                       cudecomp_t_comm_backend,cudecomp_is_t_enable_nccl,cudecomp_is_t_enable_nvshmem, &
                       cudecomp_h_comm_backend,cudecomp_is_h_enable_nccl,cudecomp_is_h_enable_nvshmem
#endif
    namelist /scalar/ &
                     alphai,beta, &
                     iniscal, &
                     cbcscal,bcscal, &
                     ssource, &
                     is_sforced, &
                     scalf, &
                     is_boussinesq_buoyancy
    namelist /numerics/ &
                       is_impdiff,is_impdiff_1d, &
                       is_poisson_pcr_tdma, &
                       is_gridpoint_natural_channel
    namelist /other_options/ &
                            is_debug,is_debug_poisson, &
                            is_timing, &
                            is_mask_divergence_check
    !
    ! defaults
    !
    dt_f = -1.
    gacc(:) = 0.
    nscal = 0
    is_cell_length = .false.
    is_user_mesh(:) = .false.
    open(newunit=iunit,file=input_file,status='old',action='read',iostat=ierr,iomsg=c_iomsg)
      if(ierr /= 0) then
        if(myid == 0) print*, 'Error reading the input file: ', trim(c_iomsg)
        if(myid == 0) print*, 'Aborting...'
        call MPI_FINALIZE(ierr)
        close(iunit)
        error stop
      end if
      read(iunit,nml=dns,iostat=ierr,iomsg=c_iomsg)
      if(ierr /= 0) then
        if(myid == 0) print*, 'Error reading `dns` namelist: ', trim(c_iomsg)
        if(myid == 0) print*, 'Aborting...'
        call MPI_FINALIZE(ierr)
        close(iunit)
        error stop
      end if
      !
      ! read mesh.csv file - this section can be throughly optimised, it is repeating exactly the same loop twice once for
      ! allocation and once for assignment. there has to be a better way of handling this
      !
      if(any(ng(:) == 0)) then
        type(pos_array), intent(out) :: dl(3), dli(3)
        ncell_user(:) = 0
        do i_dim = 1,3
          if(ng(i_dim) == 0) then
            is_user_mesh(i_dim) = .true.
          end if
        end do
        open(newunit=iunit_mesh_file,file=mesh_file,status='old',action='read',iostat=ierr,iomsg=c_iomsg)
          if(ierr /= 0) then
            if(myid == 0) print*, 'Error opening the mesh file: ', trim(c_iomsg)
            if(myid == 0) print*, 'Aborting...'
            call MPI_FINALIZE(ierr)
            close(iunit_mesh_file)
            error stop
          end if
          do ! determine size of each user given direction in the mesh.csv
            read(iunit_mesh_file, *, iostat=ierr, iomsg=c_iomsg) dim_char, val
            if(ierr /= 0) then
              if(myid == 0) print*, 'Error reading the mesh file: ', trim(c_iomsg)
              if(myid == 0) print*, 'Aborting...'
              call MPI_FINALIZE(ierr)
              close(iunit_mesh_file)
              error stop
            end if
            !
            select case (dim_char)
            case ('x', 'X')
              npoints(1) = npoints(1) + 1
            case ('y', 'Y')
              npoints(2) = npoints(2) + 1
            case ('z', 'Z')
              npoints(3) = npoints(3) + 1
            case default
              if(myid == 0) print *, 'Warning: unknown direction present in mesh file: ', dim_char
              if(myid == 0) print *, 'Aborting...'
              call MPI_FINALIZE(ierr)
              close(iunit_mesh_file)
              error stop
            end select
          end do
          if(npoints(1) /= 0) allocate(mesh_output%x(1:npoints(1)))
          if(npoints(2) /= 0) allocate(mesh_output%y(1:npoints(2)))
          if(npoints(3) /= 0) allocate(mesh_output%z(1:npoints(3)))
          do
            read(iunit_mesh_file, *, iostat=ierr, iomsg=c_iomsg) dim_char, val
            select case(dim_char)
            case ('x', 'X')
              mesh_output%x(x_counter) = val
              x_counter = x_counter + 1
            case ('y', 'Y')
              mesh_output%y(y_counter) = val
              y_counter = y_counter + 1
            case ('z', 'Z')
              mesh_output%z(z_counter) = val
              z_counter = z_counter + 1
          end do
      !
      else
        real(rp), protected, intent(out) :: dl(3), dli(3)
        dl(:) = l(:)/(1.*ng(:))
        dli(:) = dl(:)**(-1)
          !
          ! section of code where the contents of the mesh file are read, resulting in three separate vectors which are the z%f and
          ! z%c arrays - note that for the directions which correspond to non-uniform mesh but are
          ! generated using the grid-stretching functions (igrid directions), z(igrid)%f and z(igrid)%c are still to go through the
          ! initgrid submodule 
          !
          if is_cell_length==.false. then
            ! z%f is obtained automatically
            ! z%c is obtained by (z(i)%f(j) + z(i)%f(j+1))/2
          else
            ! z%c is obtained automatically
            ! z%f is also obtained
          end if
      end if
      !
      visc = visci**(-1)
      !
      if(all([1,2,3] /= ipencil_axis)) then
        ipencil_axis = 1 ! default to one
        if(myid == 0) print*, 'Warning: prescribed value of `ipencil_axis` different than 1/2/3.', trim(c_iomsg)
        if(myid == 0) print*, 'Defaulting to `ipencil_axis = 1` (x-aligned pencils)...'
      end if
#if defined(_OPENACC)
      !
      ! reading cuDecomp parameters, if these are set
      !
      ! defaults
      !
      cudecomp_is_t_comm_autotune  = .true.
      cudecomp_is_h_comm_autotune  = .true.
      cudecomp_is_t_enable_nccl    = .true.
      cudecomp_is_h_enable_nccl    = .true.
      cudecomp_is_t_enable_nvshmem = .true.
      cudecomp_is_h_enable_nvshmem = .true.
      rewind(iunit)
      read(iunit,nml=cudecomp,iostat=ierr,iomsg=c_iomsg)
      if(ierr /= 0) then
        if(myid == 0) print*, 'Error reading cudecomp namelist: ', trim(c_iomsg)
        if(myid == 0) print*, 'Aborting...'
        call MPI_FINALIZE(ierr)
        close(iunit)
        error stop
      end if
      !
      if(cudecomp_t_comm_backend >= 1 .and. cudecomp_t_comm_backend <= 7) then
        cudecomp_is_t_comm_autotune = .false. ! do not autotune if backend is prescribed
        select case(cudecomp_t_comm_backend)
        case(1)
          cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_MPI_P2P
        case(2)
          cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_MPI_P2P_PL
        case(3)
          cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_MPI_A2A
        case(4)
          cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_NCCL
        case(5)
          cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_NCCL_PL
        case(6)
          cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_NVSHMEM
        case(7)
          cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_NVSHMEM_PL
        case default
          cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_MPI_P2P
        end select
      end if
      if(cudecomp_h_comm_backend >= 1 .and. cudecomp_h_comm_backend <= 4) then
        cudecomp_is_h_comm_autotune = .false. ! do not autotune if backend is prescribed
        select case(cudecomp_h_comm_backend)
        case(1)
          cudecomp_h_comm_backend = CUDECOMP_HALO_COMM_MPI
        case(2)
          cudecomp_h_comm_backend = CUDECOMP_HALO_COMM_MPI_BLOCKING
        case(3)
          cudecomp_h_comm_backend = CUDECOMP_HALO_COMM_NCCL
        case(4)
          cudecomp_h_comm_backend = CUDECOMP_HALO_COMM_NVSHMEM
        case(5)
          cudecomp_h_comm_backend = CUDECOMP_HALO_COMM_NVSHMEM_BLOCKING
        case default
          cudecomp_h_comm_backend = CUDECOMP_HALO_COMM_MPI
        end select
      end if
      !
      ! manually set cuDecomp out-of-place transposes by default
      !
      cudecomp_is_t_in_place = .false.
#endif
      !
      ! reading scalar transport parameters, if these are set
      !
      if(nscal > 0) then
        !
        ! allocate memory
        !
        allocate(alphai(nscal),iniscal(nscal), &
                 cbcscal(0:1,3,nscal),bcscal(0:1,3,nscal), &
                 ssource(nscal),is_sforced(nscal),scalf(nscal))
        !
        ! set default values
        !
        beta          = 0.
        ssource(:)    = 0.
        is_sforced(:) = .false.
        scalf(:)      = 0.
        !
        ! read `scalar` namelist
        !
        rewind(iunit)
        read(iunit,nml=scalar,iostat=ierr,iomsg=c_iomsg)
        if(ierr /= 0) then
          if(myid == 0) print*, 'Error reading scalar namelist: ', trim(c_iomsg)
          if(myid == 0) print*, 'Aborting...'
          call MPI_FINALIZE(ierr)
          close(iunit)
          error stop
        end if
      else
        nscal = 0 ! negative values equivalent to nscal = 0
      end if
      alpha_max = huge(1._rp)
      alpha_max = minval(alphai(1:nscal))
      alpha_max = alpha_max**(-1)
      if(is_boussinesq_buoyancy) then
        if(nscal == 0) then
          if(myid == 0) print*, 'error reading the input file: `is_boussinesq_buoyancy = T` requires `nscal > 0`.'
          if(myid == 0) print*, 'Aborting...'
          call MPI_FINALIZE(ierr)
          error stop
        end if
      end if
      !
      ! read `numerics` namelist
      !
      rewind(iunit)
      read(iunit,nml=numerics,iostat=ierr,iomsg=c_iomsg)
      if(ierr /= 0 .and. ierr /= iostat_end) then
        if(myid == 0) print*, 'Error reading `numerics` namelist: ', trim(c_iomsg)
        if(myid == 0) print*, 'Aborting...'
        call MPI_FINALIZE(ierr)
        close(iunit)
        error stop
      end if
      if(is_impdiff_1d) is_impdiff = .true.
      !
      ! read `other_options` namelist
      !
      rewind(iunit)
      read(iunit,nml=other_options,iostat=ierr,iomsg=c_iomsg)
      if(ierr /= 0 .and. ierr /= iostat_end) then
        if(myid == 0) print*, 'Error reading `other_options` namelist: ', trim(c_iomsg)
        if(myid == 0) print*, 'Aborting...'
        call MPI_FINALIZE(ierr)
        close(iunit)
        error stop
      end if
    close(iunit)
  end subroutine read_input
end module mod_param
