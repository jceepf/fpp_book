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
    integer i

    dofix=.true.
    dofix0=.true.
    file="one_turn_map.txt"



    call zhe_ini



    file="map_rad7.txt"
   
        call read_tree_zhe(t,file)
    x=0.0001d0
  
    xs0=x
        call track_TREE_probe_complex_zhe(T,xs0,.false.)
   call print(xs0,6)
   call kill_tree_zhe(t)

    file="one_turn_map5.txt"
    call read_tree_zhe(t,file)

    state=only_4d0
    x=0.001d0
    xs0=x
    call track_TREE_probe_complex_ptc(T,xs0,dofix0,dofix,state)  

   call print(xs0,6)




    end program map_alone_duan_zhe

