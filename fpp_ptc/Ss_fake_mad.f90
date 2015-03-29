module madx_ptc_module
  !  use madx_keywords
  use S_fitting_new
  implicit none
  public
  TYPE(INTERNAL_STATE),POINTER :: my_state
  TYPE(layout),POINTER :: my_ring,bmadl
  type(mad_universe), pointer :: m_u,m_t

contains

  subroutine ptc_INI()
    implicit none

    allocate(m_u)
    allocate(m_t)
    allocate(bmadl)
    call set_up_universe(m_u)
    call append_empty_layout(m_u)
    call set_up_universe(m_t)
    call append_empty_layout(m_t)
    call set_up(bmadl)
    bmadl%NAME='BMAD REUSED FIBRE LAYOUT'

    call point_m_u(m_u,m_t)
  END subroutine ptc_ini

  subroutine ptc_ini_no_append()
    implicit none

    allocate(m_u)
    call set_up_universe(m_u)
    allocate(m_t)
    call set_up_universe(m_t)
    allocate(bmadl)
    call set_up(bmadl)
    bmadl%NAME='BMAD REUSED FIBRE LAYOUT'
    call point_m_u(m_u,m_t)
  END subroutine ptc_ini_no_append

  subroutine ptc_end(graphics_maybe,flat_file)
    implicit none
    integer i,i_layout
    character(120) filename
    logical, optional :: graphics_maybe,flat_file
    type(layout), pointer :: mring
    if(present(flat_file)) then
      if(flat_file) then
             mring=>m_u%start
       do i=1,m_u%n
        write(filename,*) "flat",i,".txt"
        call context(filename)
          call print_new_flat(mring,filename)
          mring=>mring%next
       enddo
      endif
    endif
    if(present(graphics_maybe)) then
      call open_gino_graphics
    endif
    if(present(graphics_maybe)) then
     call close_gino_graphics
    endif
    call kill_universe(m_t)
    call kill_universe(m_u)
    call kill_tpsa
    call kill(bmadl)
    do i=1,size(s_b)
       call nul_coef(s_b(i))
    enddo
    deallocate(s_b)



    firsttime_coef=.true.

  end subroutine ptc_end


end module madx_ptc_module

character * 48 function charconv(tint)
  !----------------------------------------------------------------------*
  ! purpose:                                                             *
  !   converts integer array to string (based on ascii)                  *
  ! input:                                                               *
  !   tint  (int array)  1 = length, rest = string                       *
  !----------------------------------------------------------------------*
  implicit none
  integer tint(*)
  integer i, j, m, n
  parameter (m = 128)
  character *(m) letter
  data letter /                                                     &
       &'                                !"#$%&''()*+,-./0123456789:;<=>?@&
       &ABCDEFGHIJKLMNOPQRSTUVWXYZ[ ]^_`abcdefghijklmnopqrstuvwxyz{|}~'/
  charconv = ' '
  n = tint(1)
  do i = 1, n
     j = tint(i+1)
     if (j .lt. m)  charconv(i:i) = letter(j:j)
  enddo
end function charconv
