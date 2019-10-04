
    !left just of max 16 char string (+1 for C null termination)
    function ch16lft(in)
     implicit none
     character(*) :: in
     character(len=17) :: ch16lft

     write(ch16lft,'(a16)') in
     ch16lft = adjustl(ch16lft)
    end function ch16lft

    subroutine putGnormaltable(gen,maxorder)
    !gets generating function that are the linear part of A_t
      !use madx_ptc_module
      use pointer_lattice
      use c_TPSA
      implicit none
      type(c_taylor) :: gen
      integer     :: order
      integer     :: ind(10), i
      character(len=17):: nn, nick
      logical skew
      integer     	:: r, myn1,myn2,mynres,o
      complex(dp)   :: c_val
      real(dp)    :: im_val, re_val, d_val,  piotreps=1e-6
      integer     :: maxorder
      character(len=250)   :: fmt
      character(len=17) :: ch16lft
      

      call print(gen,6)
      

      ind(:) = 0
      myn1 = 0
      myn2 = 0
      mynres = 0
      i=1
      call c_taylor_cycle(gen,size=mynres)


      do o=1,maxorder !print order by order, I don't know how to sort c_taylor (piotr)

        do r=1,mynres

          call c_taylor_cycle(gen,ii=r,value=c_val,j=ind(1:c_%nv))
          
          
          order = sum(ind(1:6))

          if ( order .ne. o) then
            cycle
          endif

          !print*,"GNFU ",ind(1:6)

          im_val = imag(c_val)
          re_val = real(c_val)
          d_val  = hypot(re_val, im_val)

          ! if amplitude is close to zero then it is not worth to output
          if (d_val .lt. piotreps) then
            print*,"putGnormaltable idx=",r," ",d_val," smaller then piotreps=",piotreps, " skipping "
            cycle
          endif


          write(nn,'(a4,6(a1,i1))') 'gnfa','_',ind(1),'_',ind(2),'_',ind(3), &
                                          '_',ind(4),'_',ind(5),'_',ind(6)

         ! write(nick,'(a2,6(i1))') 'f_',ind(1),ind(2),ind(3), &
         !           	ind(4),ind(5),ind(6)

          write (fmt,'(a,i1,a)')  '(a2,2(a16,1x),ES16.8,',7,'(1x,i16))'
          write(6,fmt) '  ',ch16lft(nn),  ch16lft(nn), &
                         d_val, order, ind(1:6)




          write(nn,'(a4,6(a1,i1))') 'gnfs','_',ind(1),'_',ind(2),'_',ind(3), &
                                          '_',ind(4),'_',ind(5),'_',ind(6)
         ! write(nick,'(a2,6(i1),a3)') 'f_',ind(1),ind(2),ind(3), &
         !                                  ind(4),ind(5),ind(6),'_im'
          write(6,fmt) '  ',ch16lft(nn),  ch16lft(nn), &
                         im_val, order, ind(1:6)

          write(nn,'(a4,6(a1,i1))') 'gnfc','_',ind(1),'_',ind(2),'_',ind(3), &
                                          '_',ind(4),'_',ind(5),'_',ind(6)
         ! write(nick,'(a2,6(i1),a3)') 'f_',ind(1),ind(2),ind(3), &
         !                                  ind(4),ind(5),ind(6),'_re'
          write(6,fmt) '  ',ch16lft(nn),  ch16lft(nn), &
                         re_val, order, ind(1:6)

        enddo
      enddo

      myn1 = 0
      myn2 = 0
      mynres = 0


    end subroutine putGnormaltable

program modulated_map
use madx_ptc_module
use pointer_lattice
use c_TPSA
implicit none

    interface
       subroutine build_lattice_als(ALS,MIS,error,exact,sl,thin,onecell)
         use madx_ptc_module
         use pointer_lattice
         implicit none
         type(layout), target :: ALS
         logical(lp) mis
         real(dp),optional :: error(6)
         logical, optional :: exact,sl,thin,onecell
       end subroutine build_lattice_als
    end interface
	
type(layout), pointer:: ALS
real(dp) prec,closed_orbit(6) 
real(dp) L,Kq,k0,mu_mod,beta,dmu,mu_x,circ,energy,deltap
 
type(internal_state),target :: state 
logical(lp) :: mis=.false. ,rf
type(c_damap)  one_turn_map, quasi_diagonal, diagonal,a_ac, id
type(c_normal_form) normal_form
integer :: pos =1 
integer i,map_order,mf1,mfmap
type(probe) ray_closed
type(probe_8) ray
type(real_8) y(6)
complex(dp) g1,g2,al1,al2
type(fibre),pointer :: p
type(work) werk
!piotr
type(c_damap)  :: c_Map, c_Map2, q_Map, a_cs, a_cs_1
type(c_taylor)  :: nrmlzdPseudoHam, g_io
type(c_vector_field) vf, vf_kernel

!!!!!!!!!!!!!!!!!!!!!

use_quaternion=.true.


c_verbose=.false.
prec=1.d-10 ! for printing
longprint=.false. 

call ptc_ini_no_append
call append_empty_layout(m_u)
ALS=>m_u%start


 
call kanalnummer(mfmap,"maps.txt") 


state= nocavity0 + modulation0 !
!state= nocavity0 

map_order=3
call init_all(state,map_order,0)


call alloc(y)
call alloc(one_turn_map,quasi_diagonal, diagonal,id,a_ac)
call alloc(normal_form)
call alloc(ray)


call build_lattice_als(ALS,mis,exact=.false.) 
 


!!!! circ is the circumference of the ring !!!! 
call get_length(als,circ)
!!!! AC_modulate.txt sets the magnet QF1 as a modulated magnet !!!! 
rf=.true.
rf=.false.

call kanalnummer(mf1,file="AC_modulation.txt")


! setup of rf modulation variables (quantities) and their amplitudes
  if(rf) then
    ! RF cav
    mu_mod=twopi*0.0212345d0; 

    write(mf1,*) "select layout"                  
    write(mf1,*) 1
    write(mf1,*) " MODULATERF"               
    write(mf1,*) " CAV  1" ! name and number of frequencies              
    write(mf1,*) "1.d0 0 0       !DC_ac,A_ac,theta_ac"
    write(mf1,*) "1.d0   0       ! D_ac,n_ac  "
    write(mf1,*) " 0.000d0 , 0.1d0 !  d_volt , d_phas"
    write(mf1,*) " return "
  else
   !dipole
    write(mf1,*) "select layout"                  
    write(mf1,*) 1
    write(mf1,*) " MODULATE"               
    write(mf1,*) " BEND1 1" ! name and number of frequencies              
    write(mf1,*) "1.d0 0 0       !DC_ac,A_ac,theta_ac"
    write(mf1,*) "1.d0   1       ! D_ac,n_ac  "
    write(mf1,*) "1 0.1d0 0      ! n d_bn(n) d_an(n)  "  ! (A)
    write(mf1,*) "0  0 0 " 
    write(mf1,*) " return "
    mu_mod=twopi*0.0212345d0; 
    
endif
close(mf1)
 
 call read_ptc_command77("AC_modulation.txt")
 
!  junk_e=.false.
 !call read_ptc_command77("compare.txt")

! stop
 p=>als%start
 call move_to(als,p,"CAV")
 write(6,*) "name: ", p%mag%name
 write(6,*) "nmul and phase: ", p%mag%p%nmul,p%mag%phas
 

 
!!!! set a modulation clock !!!!!!

ray_closed%ac(1)%om=mu_mod/circ ! (B1) differs from the first edition 
ray_closed%ac(1)%x=0.d0 ;       ! (B2) differs from the first edition 
write(6,*) " Modulation tune in radians =",circ*ray_closed%ac(1)%om

! Definition of clock
probe_graphical%nac=1
probe_graphical%ac(1)%om=mu_mod/circ 
probe_graphical%ac(1)%x(1)=1.d0
probe_graphical%ac(1)%x(2)=0.d0
closed_orbit=0.d0;                                                   ! (C)

call make_node_layout(als)

closed_orbit=0
ray_closed=closed_orbit  

     werk=als%start
     write(6,*) "gamma ", 1.d0/werk%gamma0I,werk%beta0
     write(6,*) "mass ",  werk%mass
     
 p=>als%start
 call move_to(als,p,"BEND")
 write(6,*) p%mag%name


! call propagate(ray_closed,state,fibre1=p,fibre2=p%next) 
! write(6,*) ray_closed%x(5)*werk%p0c
! stop

call find_orbit_x(als,closed_orbit(1:6),STATE,1.e-8_dp,fibre1=pos)   ! (D)
 write(6,*) "Closed orbit"
 write(6,*) closed_orbit(1:3)
 write(6,*) closed_orbit(4:6) 
     call GET_loss(als,energy,deltap)
     
     
     
     werk=als%start
     write(6,*) werk%beta0,1.d0/werk%gamma0I
     write(6,*) energy,deltap

ray_closed=closed_orbit     ! (E)

id=1;    
! ray= closed orbit + identity map  

ray=id+ray_closed;          ! (F)
                  
call propagate(als,RAY,state,fibre1=pos)  ! (G)
 
one_turn_map=ray                         ! (H)
write(mfmap,*); write(mfmap,*) " Map produced by code " ; write(mfmap,*);  
call print(one_turn_map,mfmap)


id=normal_form%a_t*from_phasor()
call  c_normal(one_turn_map,normal_form)    ! (J2)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

    call alloc(a_CS)

    ! raw transformation are subject to random rotations due to numerical instabilities
    ! c_canonise fixes them straight to fit the Courant Snyder format
    call c_canonise(normal_form%atot,a_CS)
    
    call alloc(vf)
    call alloc(g_io)
    call alloc(a_CS_1)



    a_CS=to_phasor()*a_CS*from_phasor()
    call c_factor_map(a_CS,a_CS_1,vf,0)

    g_io = cgetpb(vf)

    call putGnormaltable(g_io,map_order)



call kill(y)
call kill(one_turn_map,quasi_diagonal, diagonal,id,a_ac)
call kill(normal_form)
call kill(ray)
 
 close(mfmap)
 
 
call ptc_end(graphics_maybe=1,flat_file=.false.)

 
end program modulated_map

