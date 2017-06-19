!  map_alone_duan_zhe.f90 
!
!  FUNCTIONS:
!  map_alone_duan_zhe - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: map_alone_duan_zhe
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program map_alone_duan_zhe
    use duan_zhe_map
    implicit none
    TYPE(TREE_ELEMENT) t(3)
    character(120) file
    type(probe) xs0
    real(dp) x(6)
    type(internal_state) state
    logical dofix0,dofix
    dofix=.true.
    dofix0=.true.
    file="one_turn_map.txt"
    file="one_turn_map5.txt"

    call zhe_ini
    call zhe_read_tree(t,file)

    state=only_4d0
    x=0.001d0
    xs0=x
    call track_TREE_probe_complex(T,xs0,dofix0,dofix,state)  !,jump)

   call print(xs0,6)

    end program map_alone_duan_zhe

