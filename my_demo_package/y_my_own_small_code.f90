module my_own_little_code_utilities
use my_own_da
use my_analysis
integer, parameter :: nmul=3
integer :: ptc_def=0

private add_r_m,add_m_r,equal_ray_map,equal_map_ray,equal_real_ray,equal_ray_ray
 

real(dp) :: e1 = 0


type ray
 type(my_taylor) z(3)
! real(dp) stochastic_envelope(2,2) !@  &nbsp; not used yet
end type ray 

type ray8
 real(dp) z(3)
! real(dp) stochastic_envelope(2,2) !@  &nbsp; not used yet
end type ray8 

type magnet
 real(dp) L
 type(my_taylor) bn(0:nmul) 
 real(dp) :: bnr(0:nmul) = 0.0_dp
 real(dp) h 
 integer n
 character(8) name
end type magnet

  INTERFACE assignment (=)
     MODULE PROCEDURE equal_ray_map
     MODULE PROCEDURE equal_map_ray
     MODULE PROCEDURE equal_real_ray
     MODULE PROCEDURE equal_ray_ray
  end  INTERFACE  
  
    INTERFACE OPERATOR (+)
     MODULE PROCEDURE add_r_m   !@1 &nbsp;  ray+map  
     MODULE PROCEDURE add_m_r   !@1 &nbsp;  map+ray  
  END INTERFACE





contains 
!!!
  FUNCTION add_r_m( S1, S2 )
    implicit none
    TYPE (ray) add_r_m 
    real(dp), INTENT (IN) :: S1(3) 
    TYPE (my_map), INTENT (IN) :: s2
    TYPE (my_map) a
    a=s2
    a%v(1)%a(0)=0
    a%v(2)%a(0)=0

      add_r_m%z(1)=a%v(1)+s1(1)
      add_r_m%z(2)=a%v(2)+s1(2)
     if(delta_is_3rd_parameter) then
      add_r_m%z(3)=s1(3)+(1.0_dp.mon.3)    
     else
      add_r_m%z(3)=s1(3)
     endif
  end FUNCTION add_r_m

  FUNCTION add_m_r( S2, S1 )
    implicit none
    TYPE (ray) add_m_r 
    real(dp), INTENT (IN) :: S1(3) 
    TYPE (my_map), INTENT (IN) :: s2
    TYPE (my_map) a
    a=s2
    a%v(1)%a(0)=0
    a%v(2)%a(0)=0

     add_m_r%z(1)=a%v(1)+s1(1)
     add_m_r%z(2)=a%v(2)+s1(2)
     if(delta_is_3rd_parameter) then
      add_m_r%z(3)=s1(3)+(1.0_dp.mon.3)    
     else
      add_m_r%z(3)=s1(3)
     endif
     
  end FUNCTION add_m_r


 subroutine equal_ray_map(r,m)
 implicit none
 type(my_map), intent(in) :: m
 type(ray), intent(inout) :: r
  r%z(1)=m%v(1)
  r%z(2)=m%v(2)
 end subroutine equal_ray_map

 subroutine equal_ray_ray(r,m)
 implicit none
 type(ray), intent(in) :: m
 type(ray), intent(inout) :: r
  r%z(1)=m%z(1)
  r%z(2)=m%z(2)
  r%z(3)=m%z(3)
 end subroutine equal_ray_ray

 subroutine equal_map_ray(m,r)
 implicit none
 type(my_map), intent(inout) :: m
 type(ray), intent(in) :: r

  m%v(1)=r%z(1)
  m%v(2)=r%z(2)

 end subroutine equal_map_ray

 subroutine equal_real_ray(f,r)
 implicit none
 real(dp), intent(inout) :: f(3)
 type(ray), intent(in) :: r

  f(1)=r%z(1)
  f(2)=r%z(2)

 end subroutine equal_real_ray





  subroutine CANONISE_ray( z,PHASE_ADVANCE, a_cs,disp,A_l,A_nl )    !@1 &nbsp; Puts the ray into a special form
    implicit none
    TYPE (MY_TAYLOR) , INTENT (inout) ::   z(2)
    TYPE(MY_TAYLOR), INTENT (inout) ::  PHASE_ADVANCE
    TYPE(my_map), INTENT (inout) ::  a_cs,disp,A_l,A_nl
    TYPE (my_map) a,r
    a=z
     call CANONISE( A,A_cs,disp,A_l,A_nl,R,PHASE_ADVANCE )
    z=a_cs 
    
  end subroutine CANONISE_ray

  subroutine track_magnet8( r,mag )
    implicit none
    real(dp) z(3),rad
    TYPE (ray8) , INTENT (inout) ::   r
    TYPE (magnet) , INTENT (in) ::   mag
    integer n,i,j
    real(dp) dl,fac


    if(mag%l/=0) then  
    do i=1,mag%n
     call track_step8( r,mag )
    enddo
   else
        z=r%z
     fac=1.0_dp
     do j=0,nmul
      z(2)=z(2)-mag%bnr(j)*z(1)**(j)/fac
      fac=fac*(j+1)
    enddo   
       r%z=z 
   endif
   

   
  end subroutine track_magnet8

  subroutine track_step8( r,mag )
    implicit none
    real(dp) z(3),rad
    TYPE (ray8) , INTENT (inout) ::   r
    TYPE (magnet) , INTENT (in) ::   mag
    integer n,i,j
    real(dp) dl,fac

     dl=mag%L/mag%n  
        z=r%z
     z(1)=z(1)+DL/2.d0*z(2)/(1.d0+z(3))*(1.0_dp+e1*mag%bnr(1))
     fac=1.0_dp
     do j=0,nmul
      z(2)=z(2)-mag%bnr(j)*dl*z(1)**(j)/fac
       fac=fac*(j+1)
     enddo
 
! Ideal sector bend focussing and removal of design orbit  b(0)=h
      z(2)=z(2)-mag%bnr(0)*mag%h*dl*z(1)+(1.d0+z(3))*mag%h*dl
! End of sector adjustments
! fake radiation
   if(pseudo_radiation)  then
     rad=0.0_dp
     fac=1.0_dp
     do j=0,nmul
      rad=rad+(mag%bnr(j)*z(1)**(j)/fac)**2
       fac=fac*(j+1)
     enddo
     z(1)=z(1)-dl*crad*rad*z(1)
   endif
!End fake radiation
     z(1)=z(1)+DL/2.d0*z(2)/(1.d0+z(3))*(1.0_dp+e1*mag%bnr(1))


   
   r%z=z
   
  end subroutine track_step8

    subroutine track_lattice8( r,L,loc1, loc2)
    implicit none
    TYPE (magnet) , INTENT (in) ::   L(:)
    TYPE (ray8) , INTENT (inout) ::  r
    integer n,loc1, loc2,i
    

    n=size(L)
    if(loc1> n) stop 1
    if(loc2> n) stop 2
    
    if(loc2>loc1) then
     do i=loc1,loc2-1
      call track_magnet8( r,L(i) )
     enddo
    elseif(loc1>=loc2) then
     do i=loc2,n
      call track_magnet8( r,L(i) )
     enddo
     do i=1,loc1-1
      call track_magnet8( r,L(i) )
     enddo
    endif         
    end subroutine track_lattice8



  subroutine track_magnet( r,mag )
    implicit none
    TYPE (MY_TAYLOR) z(3),rad
    TYPE (ray) , INTENT (inout) ::   r
    TYPE (magnet) , INTENT (in) ::   mag
    integer n,i,j
    real(dp) dl,fac


    if(mag%l/=0) then  
    do i=1,mag%n
     call track_step( r,mag )
    enddo
   else
        z=r%z
     fac=1.0_dp
     do j=0,nmul
      z(2)=z(2)-mag%bn(ptc_def+j)*z(1)**(j)/fac
      fac=fac*(j+1)
    enddo   
       r%z=z 
   endif
   

   
  end subroutine track_magnet
  
  subroutine track_step( r,mag )
    implicit none
    TYPE (MY_TAYLOR) z(3),rad
    TYPE (ray) , INTENT (inout) ::   r
    TYPE (magnet) , INTENT (in) ::   mag
    integer n,i,j
    real(dp) dl,fac

     dl=mag%L/mag%n  
        z=r%z
     z(1)=z(1)+DL/2.d0*z(2)/(1.d0+z(3))*(1.0_dp+e1*mag%bn(ptc_def+1))
 
     fac=1.0_dp
     do j=0,nmul
      z(2)=z(2)-mag%bn(ptc_def+j)*dl*z(1)**(j)/fac
       fac=fac*(j+1)
     enddo
! Ideal sector bend focussing and removal of design orbit  b(0)=h
      z(2)=z(2)-mag%bn(ptc_def+0)*mag%h*dl*z(1)+(1.d0+z(3))*mag%h*dl
! End of sector adjustments
! fake radiation
   if(pseudo_radiation)  then
     rad=0.0_dp
     fac=1.0_dp
     do j=0,nmul
      rad=rad+(mag%bn(ptc_def+j)*z(1)**(j)/fac)**2
       fac=fac*(j+1)
     enddo
     z(1)=z(1)-dl*crad*rad*z(1)
   endif
!End fake radiation
     z(1)=z(1)+DL/2.d0*z(2)/(1.d0+z(3))*(1.0_dp+e1*mag%bn(ptc_def+1))


   
   r%z=z
   
  end subroutine track_step
  
    subroutine track_lattice( r,L,loc1, loc2)
    implicit none
    TYPE (magnet) , INTENT (in) ::   L(:)
    TYPE (ray) , INTENT (inout) ::  r
    integer n,loc1, loc2,i
    

    n=size(L)
    if(loc1> n) stop 1
    if(loc2> n) stop 2
    
    if(loc2>loc1) then
     do i=loc1,loc2-1
      call track_magnet( r,L(i) )
     enddo
    elseif(loc1>=loc2) then
     do i=loc2,n
      call track_magnet( r,L(i) )
     enddo
     do i=1,loc1-1
      call track_magnet( r,L(i) )
     enddo
    endif         
    end subroutine track_lattice
  
    subroutine find_closed_orbit( fix,L,loc)
    implicit none
    real(dp) fix(3),dfix(2)
    TYPE (magnet) , INTENT (in) ::   L(:)
    type(my_map) m,id
    type(ray) r
    integer loc,i,my_old_taylor
    real(dp) dnorm1,dnorm2
    
    id=1
  !  fix=0.0_dp
    my_old_taylor=my_order
    my_order=1
    dnorm1=-1
    
    
    do i=1,100
       m=1
       r=fix+m
        
       call track_lattice(r,L,loc,loc)
       
       
       m=r%z
       
       m%v(1)=m%v(1)-id%v(1)-fix(1)
       m%v(2)=m%v(2)-id%v(2)-fix(2)
       m=m.oo.(-1)
       dfix=m
       fix(1:2)=fix(1:2)+dfix
        dnorm2=abs(dfix(1))+abs(dfix(2))
       if(i>10) then
        if(dnorm2>=dnorm1) exit
       endif
        dnorm1=dnorm2
    enddo
     if(i>99) write(6,*) " Find_closed_orbit did not converged "
    
    my_order=my_old_taylor

    end subroutine find_closed_orbit
end module my_own_little_code_utilities

