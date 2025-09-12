! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2025 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_initsolver
  use, intrinsic :: iso_c_binding, only: C_PTR
  use mod_fft  , only: fftini
  use mod_types
  implicit none
  private
  public initsolver
  contains
  subroutine initsolver(ng,n_x_fft,n_y_fft,lo_z,hi_z,dxci_g,dxfi_g,dyci_g,dyfi_g,dzci_g,dzfi_g,cbc,bc,lambdaxy,dxdy, &
                        c_or_f,a,b,c,arrplan,normfft,rhsbx,rhsby,rhsbz)
    !
    ! initializes the Poisson/Helmholtz solver
    !
    implicit none
    integer , intent(in), dimension(3) :: ng,n_x_fft,n_y_fft,lo_z,hi_z
    real(rp), intent(in), dimension(0:) :: dxci_g,dxfi_g,dyci_g,dyfi_g,dzci_g,dzfi_g
    character(len=1), intent(in), dimension(0:1,3) :: cbc
    real(rp)        , intent(in), dimension(0:1,3) :: bc
    real(rp), intent(out), dimension(lo_z(1):,lo_z(2):) :: lambdaxy,dxdy
    character(len=1), intent(in), dimension(3) :: c_or_f
    real(rp), intent(out), dimension(lo_z(3):) :: a,b,c
#if !defined(_OPENACC) || defined(_USE_HIP)
    type(C_PTR), intent(out), dimension(2,2) :: arrplan
#else
    integer    , intent(out), dimension(2,2) :: arrplan
#endif
    real(rp), intent(out), dimension(:,:,0:) :: rhsbx
    real(rp), intent(out), dimension(:,:,0:) :: rhsby
    real(rp), intent(out), dimension(:,:,0:) :: rhsbz
    real(rp), intent(out) :: normfft
    real(rp), dimension(3)         :: dl
    real(rp), dimension(0:ng(1)+1) :: dxc_g,dxf_g
    real(rp), dimension(0:ng(2)+1) :: dyc_g,dyf_g
    real(rp), dimension(0:ng(3)+1) :: dzc_g,dzf_g
    integer :: i,j,info
    real(rp), dimension(ng(1)) :: lambdax
    real(rp), dimension(ng(2)) :: lambday
    real(rp), dimension(ng(1)-1) :: ax_g, cx_g
    real(rp), dimension(ng(2)-1) :: ay_g, cy_g
    real(rp), dimension(ng(3)-1) :: az_g, cz_g
    real(rp), dimension(ng(1))   :: bx_g
    real(rp), dimension(ng(2))   :: by_g
    real(rp), dimension(ng(3))   :: bz_g
    logical,  dimension(2) :: non_uniform_grid = .false.
    !
    if(any(dxc_g)/=dxc_g(lbound(dxc_g))) then
      non_uniform_grid(1) = .true.
    end if 
    if(any(dyc_g)/=dyc_g(lbound(dyc_g))) then
      non_uniform_grid(2) = .true.
    end if
    !
    dxc_g(:) = dxci_g(:)**(-1)
    dxf_g(:) = dxfi_g(:)**(-1)
    dyc_g(:) = dyci_g(:)**(-1)
    dyf_g(:) = dyfi_g(:)**(-1)
    dzc_g(:) = dzci_g(:)**(-1)
    dzf_g(:) = dzfi_g(:)**(-1)
    !
    ! compute tridiagonal matrix for non uniform grid
    !
    if(non_uniform_grid(1)) then
      call tridmatrix(cbc(:,1),ng(1),dxci_g,dxfi_g,c_or_f(1),ax_g,bx_g,cx_g)
    end if
    if(non_uniform_grid(2)) then
      call tridmatrix(cbc(:,2),ng(2),dyci_g,dxfi_g,c_or_f(2),ay_g,by_g,cy_g) 
    end if
    !
    ! compute and distribute coefficients for tridiagonal solver
    !
    call tridmatrix(cbc(:,3),ng(3),dzci_g,dzfi_g,c_or_f(3),az_g,bz_g,cz_g)
    !
    ! z tri-diagonal system is not symmetrized -- r.h.s. requires additional scaling
    ! deleted scaling for z since it is no longer required for any case
    !
    !
    ! generating eigenvalues (and eigenvectors for non uniform grids)
    !
    dl(1) = dxf_g(0) ! == dxc_g(0)
    dl(2) = dyf_g(0) ! == dyc_g(0)
    !
    if(.not. non_uniform_grid(1)) then
      call eigenvalues(ng(1),cbc(:,1),c_or_f(1),lambdax)
      lambdax(:) = lambdax(:)/dl(1)
    else
      call dstevd('V', )
    end if
    !
    if(.not. non_uniform_grid(2)) then
      call eigenvalues(ng(2),cbc(:,2),c_or_f(2),lambday)
      lambday(:) = lambday(:)/dl(2)
    else
      call dstevd('V', )
    end if
    !
    ! add scaled eigenvalues and store cell areas
    !
    select case(c_or_f(1)//c_or_f(2)//c_or_f(3))
    case('ccc','ccf')
      do j=lo_z(2),hi_z(2)
        do i=lo_z(1),hi_z(1)
          lambdaxy(i,j) = lambdax(i)+lambday(j)
          dxdy(i,j) = dxf_g(i)*dyf_g(j)
        end do
      end do
    case('fcc')
      do j=lo_z(2),hi_z(2)
        do i=lo_z(1),hi_z(1)
          lambdaxy(i,j) = lambdax(i)+lambday(j)
          dxdy(i,j) = dxc_g(i)*dyf_g(j)
        end do
      end do
    case('cfc')
      do j=lo_z(2),hi_z(2)
        do i=lo_z(1),hi_z(1)
          lambdaxy(i,j) = lambdax(i)+lambday(j)
          dxdy(i,j) = dxf_g(i)*dyc_g(j)
        end do
      end do
    case('fcc')
      do j=lo_z(2),hi_z(2)
        do i=lo_z(1),hi_z(1)
          lambdaxy(i,j) = lambdax(i)+lambday(j)
          dxdy(i,j) = dxc_g(i)*dyf_g(j)
        end do
      end do
    case('cfc')
      do j=lo_z(2),hi_z(2)
        do i=lo_z(1),hi_z(1)
          lambdaxy(i,j) = lambdax(i)+lambday(j)
          dxdy(i,j) = dxf_g(i)*dyc_g(j)
        end do
      end do
    end select
    !
    ! compute values to be added to the right hand side
    !
    if(     c_or_f(1) == 'c') then
      call bc_rhs(cbc(:,1),bc(:,1),[dxc_g(0),dxc_g(ng(1)  )],[dxf_g(1),dxf_g(ng(1))],c_or_f(1),rhsbx)
    else if(c_or_f(1) == 'f') then
      call bc_rhs(cbc(:,1),bc(:,1),[dxc_g(1),dxc_g(ng(1)-1)],[dxf_g(1),dxf_g(ng(1))],c_or_f(1),rhsbx)
    end if
    if(     c_or_f(2) == 'c') then
      call bc_rhs(cbc(:,2),bc(:,2),[dyc_g(0),dyc_g(ng(2)  )],[dyf_g(1),dyf_g(ng(2))],c_or_f(2),rhsby)
    else if(c_or_f(2) == 'f') then
      call bc_rhs(cbc(:,2),bc(:,2),[dyc_g(1),dyc_g(ng(2)-1)],[dyf_g(1),dyf_g(ng(2))],c_or_f(2),rhsby)
    end if
    !
    ! z tri-diagonal system is not symmetrized -- r.h.s. requires additional scaling
    !
    if(     c_or_f(3) == 'c') then
      call bc_rhs(cbc(:,3),bc(:,3),[dzc_g(0),dzc_g(ng(3)  )],[dzf_g(1),dzf_g(ng(3))],c_or_f(3),rhsbz)
      rhsbz(:,:,0) = rhsbz(:,:,0)/dzf_g(1)
      rhsbz(:,:,1) = rhsbz(:,:,1)/dzf_g(ng(3))
    else if(c_or_f(3) == 'f') then
      call bc_rhs(cbc(:,3),bc(:,3),[dzc_g(1),dzc_g(ng(3)-1)],[dzf_g(1),dzf_g(ng(3))],c_or_f(3),rhsbz)
      rhsbz(:,:,0) = rhsbz(:,:,0)/dzc_g(1)
      rhsbz(:,:,1) = rhsbz(:,:,1)/dzc_g(ng(3)-1)
    end if
    !
    ! prepare ffts
    !
    call fftini(ng,n_x_fft,n_y_fft,cbc(:,1:2),c_or_f(1:2),arrplan,normfft)
  end subroutine initsolver
  !
  subroutine eigenvalues(n,bc,c_or_f,lambda)
    use mod_param, only: pi
    implicit none
    integer , intent(in ) :: n
    character(len=1), intent(in), dimension(0:1) :: bc
    character(len=1), intent(in) :: c_or_f ! c -> cell-centered; f -> face-centered
    real(rp), intent(out), dimension(n) :: lambda
    integer :: l
    select case(bc(0)//bc(1))
    case('PP')
      do l=1,n
        lambda(l  )   = -2.*(1.-cos((2*(l-1))*pi/(1.*n)))
      end do
#if defined(_OPENACC)
      block
        !
        ! new format: (r[0],r[n],r[1],i[1],...,r[n-1],i[n-1])
        ! note that i[0] = i[n] = 0 in a R2C DFT
        !
        integer :: nh,iswap(n)
        nh = (n+1)/2
        iswap(1) = 1
        iswap(2) = nh+(1-mod(n,2))
        do l=2,n-1
          if(l <= nh) then ! real eigenvalue
            iswap(2*l-1                  ) = l
          else             ! imaginary eigenvalue
            iswap(n-2*(l-(nh+1))-mod(n,2)) = l+1
          end if
        end do
        lambda(:) = lambda(iswap(:))
      end block
#endif
    case('NN')
      if(     c_or_f == 'c') then
        do l=1,n
          lambda(l)   = -2.*(1.-cos((l-1  )*pi/(1.*n)))
        end do
      else if(c_or_f == 'f') then
        do l=1,n
          lambda(l)   = -2.*(1.-cos((l-1  )*pi/(1.*(n-1+1))))
        end do
      end if
    case('DD')
      if(     c_or_f == 'c') then
        do l=1,n
          lambda(l)   = -2.*(1.-cos((l    )*pi/(1.*n)))
        end do
      else if(c_or_f == 'f') then
        do l=1,n-1 ! point at n is a boundary and is excluded here
          lambda(l)   = -2.*(1.-cos((l    )*pi/(1.*(n+1-1))))
        end do
        lambda(n) = 0.
      end if
    case('ND','DN')
      do l=1,n
        lambda(l)     = -2.*(1.-cos((2*l-1)*pi/(2.*n)))
      end do
    end select
  end subroutine eigenvalues
        a(k) = dzci(k-1)
        c(k) = dzci(k)
      end do
    case('f')
      do k = 1,n
        a(k) = dzfi(k)
        c(k) = dzfi(k+1)
      end do
    end select
    b(:) = -(a(:)+c(:))
    do ibound = 0,1
      select case(bc(ibound))
      case('P')
        factor(ibound) = 0.
      case('D')
        factor(ibound) = -1.
      case('N')
        factor(ibound) = 1.
      end select
    end do
    select case(c_or_f)
    case('c')
      b(1) = b(1) + factor(0)*a(1)
      b(n) = b(n) + factor(1)*c(n)
    case('f')
      if(bc(0) == 'N') b(1) = b(1) + factor(0)*a(1)
      if(bc(1) == 'N') b(n) = b(n) + factor(1)*c(n)
    end select
  end subroutine tridmatrix
  !
  subroutine bc_rhs(cbc,bc,dlc,dlf,c_or_f,rhs)
    implicit none
    character(len=1), intent(in), dimension(0:1) :: cbc
    real(rp), intent(in), dimension(0:1) :: bc
    real(rp), intent(in), dimension(0:1) :: dlc,dlf
    real(rp), intent(out), dimension(:,:,0:) :: rhs
    character(len=1), intent(in) :: c_or_f ! c -> cell-centered; f -> face-centered
    real(rp) :: factor
    real(rp) :: sgn
    integer :: ibound
    !
    select case(c_or_f)
    case('c')
      do ibound = 0,1
        select case(cbc(ibound))
        case('P')
          factor = 0.
        case('D')
          factor = -2.*bc(ibound)
        case('N')
          if(ibound == 0) sgn =  1.
          if(ibound == 1) sgn = -1.
          factor = sgn*dlc(ibound)*bc(ibound)
        end select
        rhs(:,:,ibound) = factor/dlc(ibound)
      end do
    case('f')
      do ibound = 0,1
        select case(cbc(ibound))
        case('P')
          factor = 0.
        case('D')
          factor = -bc(ibound)
        case('N')
          if(ibound == 0) sgn =  1.
          if(ibound == 1) sgn = -1.
          factor = sgn*dlf(ibound)*bc(ibound)
        end select
        rhs(:,:,ibound) = factor/dlf(ibound)
      end do
    end select
  end subroutine bc_rhs
end module mod_initsolver
