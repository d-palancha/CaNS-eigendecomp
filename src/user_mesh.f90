! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_user_mesh
  use mod_types
  use mod_param
  use mod_initgrid, only: pos_array
  implicit none
  private
  public user_mesh
  contains
  subroutine user_mesh(mesh_input, is_cell_length, ng, l, dl, dli, z_g, dz_g, dzi_g)
    type(xyz),       intent(in) :: mesh_input
    logical,         intent(in) :: is_cell_length
    integer,         intent(inout), dimension(3) :: ng
    real(rp),        intent(in),    dimension(3) :: l
    type(pos_array), intent(out),   dimension(3) :: dl, dli
    type(pos_array), intent(out),   dimension(3) :: z_g, dz_g, dzi_g
    real(rp) :: dl_temp
    integer :: k, i
    !
    if(is_cell_length == .false.) then
      do i = 1, 3
        if(ng(i) == 0) then
          select case i ! THIS IS NOT CORRECT AS THE SIZE OF z(i)%f(:) == z(i)%c(:) IN main.f90
          case (1)
            z(i)%f(:) = mesh_input%x(:)
          case (2)
            z(i)%f(:) = mesh_input%y(:)
          case (3)
            z(i)%f(:) = mesh_input%z(:)
          end select
        else
          dl_temp = l(i)/(1.*ng(i))
          do concurrent k = 0:ng(i)
            z(i)%f(k) = dl_temp*k ! NOTICE HOW z(i)%f(0:ng(i)), THIS DOES NOT MATCH MAIN
          end do
        end if
        do concurrent k = 1:ng(i)
          z(i)%c(k) = (z(i)%f(k-1)+z(i)%f(k))/2.
        end do
      end do
    elseif(is_cell_length == .true.) then
      do i = 1, 3
        if(ng(i) == 0) then
          select case i ! THIS IS NOT CORRECT AS THE SIZE OF z(i)%f(:) == z(i)%c(:) IN main.f90
          case (1)
            z(i)%c(:) = mesh_input%x(:)
          case (2)
            z(i)%c(:) = mesh_input%y(:)
          case (3)
            z(i)%c(:) = mesh_input%z(:)
          end select
        else
          dl_temp = l(i)/(1.*ng(i))
          do concurrent k = 1:ng(i) ! NOTICE HOW z(i)%c(1:ng(i)), THIS DOES NOT MATCH MAIN
            z(i)%c(k) = dl_temp*k 
          end do
        end if
        ! STILL NEED TO CONVERT z(i)%c(:) TO z(i)%f(:)
      end do
    end if
    
  end subroutine user_mesh
end module mod_user_mesh
