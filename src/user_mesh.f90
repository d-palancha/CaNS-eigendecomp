module mod_user_mesh
  use mod_types
  use mod_param
  use mod_initgrid, only: pos_array
  implicit none
  private
  public user_mesh
  contains
  subroutine user_mesh(mesh_input, is_cell_length, ng, l, dl, dli, z_g, dz_g, dzi_g, is_cell_length)
    type(xyz),       intent(in) :: mesh_input
    logical,         intent(in) :: is_cell_length
    integer,         intent(inout), dimension(3) :: ng
    real(rp),        intent(in),    dimension(3) :: l
    type(pos_array), intent(out),   dimension(3) :: dl, dli
    type(pos_array), intent(out),   dimension(3) :: z_g, dz_g, dzi_g
    integer :: i
    !
    do i = 1, 3
      if(ng(i) == 0) then
        if(is_cell_length == .false.) then
          select case i
          case (1)
            z(i)%f(:) = mesh_input%x(:)
          case (2)
            z(i)%f(:) = mesh_input%y(:)
          case (3)
            z(i)%f(:) = mesh_input%z(:)
          end select
        elseif(is_cell_length == .true.) then
          
        end if
      else
        
      end if
    end do
    mesh_input%x
    mesh_input%y
    mesh_input%z
    
  end subroutine user_mesh
