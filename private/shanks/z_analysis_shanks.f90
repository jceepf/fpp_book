module shanks_crittenden
use madx_ptc_module
use pointer_lattice
use c_TPSA
implicit none
private bfieldr,fxr,fxp,bfieldp,feval,rk4_crittendenr,rk4_crittendenp,rk4
private integrate_crittendenr,integrate_crittendenp,conv_to_xpr,conv_to_xpp,conv_to_pxr,conv_to_pxp
 type(internal_state) :: state_correct
 type(node_array) , allocatable :: sverige(:) 
! This is a beam to be tracked
integer, private :: dirb =1
  real(dp), private, allocatable :: bf(:,:,:), zp(:),dx0(:)

  INTERFACE bfield
     MODULE PROCEDURE bfieldr
     MODULE PROCEDURE bfieldp
  END INTERFACE

  INTERFACE feval
     MODULE PROCEDURE fxr
     MODULE PROCEDURE fxp
  END INTERFACE

  INTERFACE rk4 
     MODULE PROCEDURE rk4_crittendenr
     MODULE PROCEDURE rk4_crittendenp
  END INTERFACE

  INTERFACE  integrate_crittenden
     MODULE PROCEDURE integrate_crittendenr
     MODULE PROCEDURE integrate_crittendenp
  END INTERFACE

  INTERFACE conv_to_xp_crittenden
     MODULE PROCEDURE conv_to_xpr
     MODULE PROCEDURE conv_to_xpp
  END INTERFACE

  INTERFACE conv_to_px_crittenden
     MODULE PROCEDURE conv_to_pxr
     MODULE PROCEDURE conv_to_pxp
  END INTERFACE


contains



! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$   Jim Shanks Stuff $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 

 

 subroutine locate_shanks(ring,ds,n)
IMPLICIT NONE
type(layout), pointer :: ring
type(integration_node), pointer :: t
real(dp) s,ds,s1
type(internal_state) state_no
TYPE(C_TAYLOR) PHASE(3)
type(c_damap) id,m,r
type(c_normal_form) normal
type(real_8) y(6)
integer n,i,k
integer j(6)
 
s=0
s1=0
k=1
t=>ring%t%start
do i=1,ring%t%n

s1=t%s(1)-s 

if(s1>=ds) then
k=k+1
s=t%s(1)
endif


t=>t%next
enddo
  
k=k+1

call alloc_node_array(sverige,k,n)
s=0
s1=0
k=1
t=>ring%t%start
sverige(1)%t=>t 

do i=1,ring%t%n

s1=t%s(1)-s 

if(s1>=ds) then
k=k+1
s=t%s(1)
sverige(k)%t=>t
endif


t=>t%next
enddo
  sverige(k+1)%t=>ring%t%start
end subroutine locate_shanks

 subroutine analysis_shanks_3(k)
IMPLICIT NONE
integer, optional :: k
integer i,i1,i2

i1=1
i2=size(sverige)-1 
if(present(k)) then
 i1=k
 i2=k
endif

do i=i1,i2
! 3nux
 sverige(i)%s(1)=(sverige(i)%f%v(2).sub.'2000') !/3
!  nux
 sverige(i)%s(2)=(sverige(i)%f%v(2).sub.'11') !/2
!  2nux+nuy
 sverige(i)%s(3)=(sverige(i)%f%v(2).sub.'101') !/2
!  2nux-nuy
 sverige(i)%s(4)=(sverige(i)%f%v(2).sub.'1001') !/2
!  nux+2nuy 
 sverige(i)%s(5)=(sverige(i)%f%v(2).sub.'002')
!  nux-2nuy
 sverige(i)%s(6)=(sverige(i)%f%v(2).sub.'0002')
!  3nuy
 sverige(i)%s(7)=(sverige(i)%f%v(4).sub.'0020') !/3
!  nuy
 sverige(i)%s(8)=(sverige(i)%f%v(4).sub.'0011') !/2
!  chromx'
 sverige(i)%s(9)=(sverige(i)%f%v(2).sub.'01001')
!  chromy'
sverige(i)%s(10)=(sverige(i)%f%v(4).sub.'00011')
!  2nux'
sverige(i)%s(11)=(sverige(i)%f%v(2).sub.'01001') !/2
!  2nuy'
sverige(i)%s(12)=(sverige(i)%f%v(4).sub.'10001') !/2
!  (nux+nuy)'
sverige(i)%s(13)=(sverige(i)%f%v(2).sub.'00101')
!  (nux-nuy)'
sverige(i)%s(14)=(sverige(i)%f%v(2).sub.'00011')
!  etax delta**2
sverige(i)%s(15)=(sverige(i)%f%v(2).sub.'00002')
!  etay delta**2
sverige(i)%s(16)=(sverige(i)%f%v(4).sub.'00002')
enddo

end subroutine analysis_shanks_3



 subroutine analysis_shanks_fast(ring,state_no)
IMPLICIT NONE
type(layout), pointer :: ring
real(dp) x(6),xm,prec,c(2,2),d(2,2)
type(internal_state) state_no
TYPE(C_TAYLOR) PHASE(2)
type(c_damap) id,L,r,U,a,u_c,one_turn,a0,a1,a2
type(c_normal_form) normal
type(real_8) y(6),yb(6)
integer n,i,j,km,mf,sh
type(c_vector_field) logN,f
type(c_taylor) x2
logical :: dbeta=.true.
prec=1.d-8
 j=size(sverige)
!call kanalnummer(mf,"x2.txt")

!!!!!!!!!!!!!!!!!!!!!!   third order calculation !!!!!!!!!!!!!!!!!!!
n=3
 call init_all(state_no,n,0)
 call alloc(id,a, a0, a1, a2,L); call alloc(normal);call alloc(y);call alloc(sverige(j)%f);

 x=0.d0
 call find_orbit(ring,x,1,state_no,1e-5_dp)
 id=1
 y=id+x
 call propagate(ring,y,state_no,fibre1=1)
 id=y
 call c_normal(id,normal)

 

normal%ker%f(1)=normal%ker   !ker%f(1:3)    rot=exp(ker%f(1).grad)exp(ker%f(2).grad)exp(ker%f(3).grad)
normal%ker%f(1)=(1.d0/twopi)*normal%ker%f(1)   !tunes



a=normal%a_t
call c_full_canonise(a,L,a0=a0,a1=a1,a2=a2)  ! a= L*rot= a0*a1*a2*rot   L=a0*a1*a2 courant-snyder
a1=a0*(a1.sub.1)  
 
id=a1**(-1)*id*a1
call print(id%v(1),6)

pause
 L=to_phasor()*id*from_phasor() 

 call c_factor_map(L,L,sverige(j)%f,1)  

 
call analysis_shanks_3(j)


sh=24
!  nux_total
 sverige(j)%s(sh+1)=(normal%ker%f(1)%v(2).sub.'01000')
!  nux_total'
 sverige(j)%s(sh+2)=(normal%ker%f(1)%v(2).sub.'01001')
!  nux_total [2]
 sverige(j)%s(sh+3)=(normal%ker%f(1)%v(2).sub.'01002')
!  nuy_total
 sverige(j)%s(sh+4)=(normal%ker%f(1)%v(4).sub.'00010')
!  nuy_total'
 sverige(j)%s(sh+5)=(normal%ker%f(1)%v(4).sub.'00011')
!  nuy_total[2]
 sverige(j)%s(sh+6)=(normal%ker%f(1)%v(4).sub.'00012')
!  nux_total function of Jx
 sverige(j)%s(sh+7)=(normal%ker%f(1)%v(2).sub.'12000')
!  nux_total function of Jy
 sverige(j)%s(sh+8)=(normal%ker%f(1)%v(2).sub.'01110')
!  nuy_total function of Jy
 sverige(j)%s(sh+9)=(normal%ker%f(1)%v(4).sub.'00120')
! write(6,*) sverige(j)%s(1:9)


 call kill(id,a, a0, a1, a2,L); call kill(normal);call kill(y);call kill(sverige(j)%f);

n=2
 
 call init_all(state_no,n,0)
 call alloc(id,U,a,u_c,r,L,one_turn,a0,a1,a2); call alloc(normal);call alloc(y);CALL ALLOC(PHASE)
call alloc(logn);call alloc(f); call alloc(x2);call alloc(yb)
call alloc_node_array_tpsa(sverige)
 x=0.d0

 call find_orbit(ring,x,1,state_no,1e-5_dp)
 id=1


 y=id+x

 call propagate(ring,y,state_no,fibre1=1)

 id=y
 one_turn=id
 call c_normal(id,normal)






r=1
if(dbeta) then 
 yb=normal%A_t+x
endif

call  c_canonise_shanks(normal%A_t ,U_c,a)  
one_turn=u_c**(-1)*one_turn*U_c
one_turn=to_phasor()*one_turn*from_phasor() 

 y=U_c+x

do i=1,size(sverige)-1
 call propagate(y,state_no,node1=sverige(i)%t,node2=sverige(i+1)%t)
 call propagate(yb,state_no,node1=sverige(i)%t,node2=sverige(i+1)%t)

  U=y ! copying in map  ! (15b)

 if(dbeta) then

  a=yb
 
call c_full_canonise(a,id,a0=a0,a1=a1,a2=a2)  ! call c_canonise(a,id,a0,a1,a2) 

!  write(mf,*) sverige(i+1)%t%parent_fibre%pos,sverige(i+1)%t%parent_fibre%mag%name 

c(1,1)=a1%v(1).sub.'10' 
c(1,2)=a1%v(1).sub.'01'
c(2,1)=a1%v(2).sub.'10'
c(2,2)=a1%v(2).sub.'01'
d(1,1)=a1%v(1).sub.'10001' 
d(1,2)=a1%v(1).sub.'01001'
d(2,1)=a1%v(2).sub.'10001'
d(2,2)=a1%v(2).sub.'01001'

sverige(i)%s(17)=2*(c(1,1)*d(1,1)+c(1,2)*d(1,2))   ! beta_x '
sverige(i)%s(18)=-(c(1,1)*d(2,1)+d(1,1)*c(2,1)+c(1,2)*d(2,2)+d(1,2)*d(2,2)) ! alpha_x '

c(1,1)=a1%v(3).sub.'0010' 
c(1,2)=a1%v(3).sub.'0001'
c(2,1)=a1%v(4).sub.'0010'
c(2,2)=a1%v(4).sub.'0001'
d(1,1)=a1%v(3).sub.'00101' 
d(1,2)=a1%v(3).sub.'00011'
d(2,1)=a1%v(4).sub.'00101'
d(2,2)=a1%v(4).sub.'00011'
sverige(i)%s(19)=2*(c(1,1)*d(1,1)+c(1,2)*d(1,2)) ! beta_y '
sverige(i)%s(20)=-(c(1,1)*d(2,1)+d(1,1)*c(2,1)+c(1,2)*d(2,2)+d(1,2)*d(2,2)) ! alpha_y '
 
     

  sverige(i)%s(21)=(a0%v(1).sub.'00002')   ! etax [2]
  sverige(i)%s(22)=(a0%v(2).sub.'00002') ! etapx [2]
  sverige(i)%s(23)=(a0%v(3).sub.'00002') ! etay [2]
  sverige(i)%s(24)=(a0%v(4).sub.'00002') ! etapy [2]


!  x2=2.d0.cmono.'20'
!  call C_AVERAGE(x2,A1,x2)  
! write(mf,*) sverige(i)%s(17)
!   write(mf,*) "2 <x**2>" 
!  call print(x2,mf)

!  x2=(((-1.d0).cmono.'1')*  (1.d0.cmono.'01'))*2.d0
!  call C_AVERAGE(x2,A1,x2)
! write(mf,*) sverige(i)%s(18) 
!   write(mf,*) " -<x px>" 
!  call print(x2,mf)


!  x2=2.d0.cmono.'002'
!  call C_AVERAGE(x2,A1,x2) 
! write(mf,*) sverige(i)%s(19)  
! write(mf,*) "2 <y**2>" 
!  call print(x2,mf)

!  x2=(((-1.d0).cmono.'001')*  (1.d0.cmono.'0001')*2.d0)
!  call C_AVERAGE(x2,A1,x2) 
! write(mf,*) sverige(i)%s(20)
!   write(mf,*) " -<y py>" 
!  call print(x2,mf)

!  x2=  (1.d0.cmono.'0001')
!  call C_AVERAGE(x2,A0,x2) 
! write(mf,*) sverige(i)%s(24)
!   write(mf,*) " y' " 
!  call print(x2,mf)

 endif
 
  call c_canonise_shanks(U,U_c,a)

 
 ! U=U_c**(-1)*U                
 

 L=to_phasor()*U_c**(-1)*U*from_phasor() 

 call c_factor_map(L,sverige(i)%m,sverige(i)%f,1)  
 

 

 x=y
 y=x+U_c 
 if(mod(i,100)==1) write(6,*) i
 
enddo
 
  sverige(size(sverige))%f=0
  sverige(size(sverige))%m=1

!!! rewriting   the vector fields as  exp(f1)...exp(f_n) R
j=size(sverige)
call alloc(r);call alloc(f);
r=1
do i=1,size(sverige)-1
 r=sverige(i)%m*r
 sverige(i)%f=r*sverige(i)%f
 sverige(j)%f=sverige(j)%f+sverige(i)%f
enddo

 normal%ker%f(1)=normal%ker
normal%ker%f(1)=(1.d0/twopi)*normal%ker%f(1)


 

 call kill(id,U,a,u_c,r,L,one_turn,a0,a1,a2); call kill(normal);call kill(y);CALL kill(PHASE);call kill(x2)
call kill(logn);call kill(yb)
!close(mf)
end subroutine analysis_shanks_fast

subroutine c_canonise_shanks(at,a_cs,rot) 
    implicit none                             
    type(c_damap) , intent(inout) :: at,a_cs,rot 
    real(dp)  a(4,4),r(4,4), a1,a2
    r=0
    a=at
     
    a1=-atan2(a(1,2),a(1,1))
    a2=-atan2(a(3,4),a(3,3))
   
    r(1,1)=cos(a1);r(2,2)=r(1,1);r(1,2)=sin(a1);r(2,1)=-r(1,2);
    r(3,3)=cos(a2);r(4,4)=r(3,3);r(3,4)=sin(a2);r(4,3)=-r(3,4);
    
     rot=r

    a_cs=at*rot

    a=a_cs
    a1=0;a2=0;
    
    if(a(1,1)<0) then
     a1=pi
    endif
    if(a(3,3)<0) then
     a2=pi
    endif
    if(a1>0.or.a2>0) then
     write(6,*) "carp"
     r(1,1)=cos(a1);r(2,2)=r(1,1);r(1,2)=sin(a1);r(2,1)=-r(1,2);
     r(3,3)=cos(a2);r(4,4)=r(3,3);r(3,4)=sin(a2);r(4,3)=-r(3,4);   
    rot=r
    a_cs=at*rot
    endif
    a_cs=a_cs.sub.1
!    a=a_cs
!    a_cs=a
end subroutine c_canonise_shanks

!!!!!!!!!!!!!!!!!! Jim Crittenden

subroutine read_crittenden(w)
implicit None
integer ns,nfit,i,mf,j
real(dp) brho
type(work) w
ns=133
nfit=17
call kanalnummer(mf,"C:\document\etienne_programs_new\read_m_u_m_t\crittenden\sarc_16_fitpars.txt")


allocate(bf(ns,3,nfit),zp(ns),dx0(ns))

bf=0.d0
zp=0.d0
dx0=0.d0

do i=1,ns
do j=1,3
read(mf,*) zp(i),bf(i,j,1:nfit)

enddo
enddo
zp=zp/100.d0
bf=bf/w%brho

close(mf) 
end subroutine read_crittenden



 subroutine fxr(i,x,k,f,hcurv,w)
    implicit none

    real(dp)  d(3),c(6),BETA0,GAMMA0I
    real(dp)  b(3),hcurv
    type(work)  w
    real(dp) ,intent(inout) :: x(6)
    real(dp), intent(out):: f(6)
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K
    integer i
    if(k%time) then
       beta0=w%beta0;GAMMA0I=w%GAMMA0I;
    else
       beta0=1.0_dp;GAMMA0I=0.0_dp;
    endif
    
    x(1)=x(1)-dx0(i)
 call bfield(i,x(1),x(3),b)
    x(1)=x(1)+dx0(i)

    d(1)=root(x(2)**2+x(4)**2+(1.0_dp+hcurv*x(1))**2)
    d(2)=(d(1)**3)/root(1.0_dp+2*x(5)/beta0+x(5)**2)
    d(3)=1.0_dp+hcurv*x(1)

    c(1)=d(1)**2-x(2)**2
    c(2)=-x(2)*x(4)
    c(3)= x(2)*x(4)
    c(4)=-d(1)**2+x(4)**2
    c(5)=d(2)*(x(4)*b(3)-d(3)*b(2)) +hcurv*d(3)*(d(1)**2+x(2)**2)
    c(6)=d(2)*(x(2)*b(3)-d(3)*b(1)) -hcurv*d(3)*c(3)

    d(3)=c(1)*c(4)-c(2)*c(3)
    f(1)=x(2)
    f(2)=(c(4)*c(5)-c(2)*c(6))/d(3)
    f(3)=x(4)
    f(4)=(c(1)*c(6)-c(3)*c(5))/d(3)
    d(2)=1.0_dp+2.0_dp*x(5)/beta0+x(5)**2
    d(2)=gamma0I/beta0/d(2)
    f(6)=root((1+d(2)**2))*d(1)  ! (time)-prime = dt/dz
 
    f(5)=0.0_dp
 

  end subroutine fxr
 

  subroutine fxp(i,x,k,f,hcurv,w)
    implicit none

    type(real_8)  d(3),c(6)
    type(real_8) ,intent(inout) :: x(6)
    type(real_8)  b(3)
    real(dp)   BETA0,GAMMA0I,hcurv
    type(real_8), intent(out):: f(6)
    type(work) w
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K
    integer i
    call alloc(d,3)
    call alloc(c,6)
    call alloc(b,3)
 
    if(k%time) then
       beta0=w%beta0;GAMMA0I=w%GAMMA0I;
    else
       beta0=1.0_dp;GAMMA0I=0.0_dp;
    endif
    x(1)=x(1)-dx0(i)
 call bfield(i,x(1),x(3),b)
    x(1)=x(1)+dx0(i)

    d(1)=SQRT(x(2)**2+x(4)**2+(1.0_dp+hcurv*x(1))**2)
    d(2)=(d(1)**3)/SQRT(1.0_dp+2*x(5)/beta0+x(5)**2)
    d(3)=1.0_dp+hcurv*x(1)

    c(1)=d(1)**2-x(2)**2
    c(2)=-x(2)*x(4)
    c(3)= x(2)*x(4)
    c(4)=-d(1)**2+x(4)**2
    c(5)=d(2)*(x(4)*b(3)-d(3)*b(2)) +hcurv*d(3)*(d(1)**2+x(2)**2)
    c(6)=d(2)*(x(2)*b(3)-d(3)*b(1)) -hcurv*d(3)*c(3)

    d(3)=c(1)*c(4)-c(2)*c(3)
    f(1)=x(2)
    f(2)=(c(4)*c(5)-c(2)*c(6))/d(3)
    f(3)=x(4)
    f(4)=(c(1)*c(6)-c(3)*c(5))/d(3)
    f(5)=0.0_dp
    d(2)=1.0_dp+2.0_dp*x(5)/beta0+x(5)**2
 

    d(2)=gamma0I/beta0/d(2)
    f(6)=SQRT((1+d(2)**2))*d(1)  ! (time)-prime = dt/dz

    f(5)=0.0_dp

    call kill(d,3)
    call kill(c,6)
    call kill(b,3)
  end subroutine fxp

  subroutine bfieldr(j,x,y,b)
    implicit none
    integer i,j
    real(dp),intent(inout):: b(3),x,y
do i=1,3
 b(j)=bf(i,j,1)+bf(i,j,2)*x+bf(i,j,3)*y
 b(j)=b(j)+bf(i,j,4)*x*y+bf(i,j,5)*x**2+bf(i,j,6)*y**2
 b(j)=b(j)+bf(i,j,7)*y*x**2+bf(i,j,8)*x*y**2+ bf(i,j,9)*x**3+bf(i,j,10)*y*x**3+bf(i,j,11)*y**2*x**3
 b(j)=b(j)+bf(i,j,12)*x**4+bf(i,j,13)*y*x**4+bf(i,j,14)*y**2*x**4
 b(j)=b(j)+bf(i,j,15)*x**5+bf(i,j,16)*y*x**5+bf(i,j,17)*y**2*x**5
enddo
  end subroutine bfieldr


  subroutine bfieldp(j,x,y,b)
    implicit none
    integer i,j
    type(real_8),intent(inout):: b(3),x,y

do i=1,3
 b(j)=bf(i,j,1)+bf(i,j,2)*x+bf(i,j,3)*y
 b(j)=b(j)+bf(i,j,4)*x*y+bf(i,j,5)*x**2+bf(i,j,6)*y**2
 b(j)=b(j)+bf(i,j,7)*y*x**2+bf(i,j,8)*x*y**2+ bf(i,j,9)*x**3+bf(i,j,10)*y*x**3+bf(i,j,11)*y**2*x**3
 b(j)=b(j)+bf(i,j,12)*x**4+bf(i,j,13)*y*x**4+bf(i,j,14)*y**2*x**4
 b(j)=b(j)+bf(i,j,15)*x**5+bf(i,j,16)*y*x**5+bf(i,j,17)*y**2*x**5
enddo

  end subroutine bfieldp

 subroutine rk4_crittendenr(ti,h,hcurv,w,y,k)
    IMPLICIT none

    integer ne
    parameter (ne=6)
    real(dp), INTENT(INOUT)::  y(ne)
    real(dp)  yt(ne),f(ne),a(ne),b(ne),c(ne),d(ne)
    integer j
    real(dp), intent(inout) :: h,hcurv
    integer, intent(inout) :: ti
    type(work) w
    INTEGER TT
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K

    call feval(tI,y,k,f,hcurv,w)
    do  j=1,ne
       a(j)=h*f(j)
    enddo
    do  j=1,ne
       yt(j)=y(j)+a(j)/2.0_dp
    enddo

    tt=tI+dirb
    call feval(tt,yt,k,f,hcurv,w)
    do  j=1,ne
       b(j)=h*f(j)
    enddo
    do   j=1,ne
       yt(j)=y(j) + b(j)/2.0_dp
    enddo


    tt=tI+dirb
    call feval(tt,yt,k,f,hcurv,w)
    do  j=1,ne
       c(j)=h*f(j)
    enddo
    do  j=1,ne
       yt(j)=y(j)+c(j)
    enddo

    tt=tI+2*dirb
    call feval(tt,yt,k,f,hcurv,w)
    do  j=1,ne
       d(j)=h*f(j)
    enddo

    do  j=1,ne
       y(j) = y(j)+(a(j)+2.0_dp*b(j)+2.0_dp*c(j)+d(j))/6.0_dp
    enddo
    tI=ti+2*dirb

 !   if(k%TIME) then
 !      Y(6)=Y(6)-(1-k%TOTALPATH)*GR%P%LD/GR%P%beta0/GR%P%nst
 !   else
 !      Y(6)=Y(6)-(1-k%TOTALPATH)*GR%P%LD/GR%P%nst
 !   endif

    return
  end  subroutine rk4_crittendenr

 
 

 

 subroutine rk4_crittendenp(ti,h,hcurv,w,y,k)
    IMPLICIT none

    type(work) w
    integer ne
    parameter (ne=6)
    TYPE(REAL_8), INTENT(INOUT)::  y(ne)
    TYPE(REAL_8)  yt(ne),f(ne),a(ne),b(ne),c(ne),d(ne)
    integer j
    real(dp), intent(inout) :: h,hcurv
    integer, intent(inout) :: ti
    INTEGER TT 
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K
    call alloc(yt,ne)
    call alloc(f,ne)
    call alloc(a,ne)
    call alloc(b,ne)
    call alloc(c,ne)
    call alloc(d,ne)


    call feval(tI,y,k,f,hcurv,w)
    do  j=1,ne
       a(j)=h*f(j)
    enddo
    do  j=1,ne
       yt(j)=y(j)+a(j)/2.0_dp
    enddo

    tt=tI+dirb
    call feval(tt,yt,k,f,hcurv,w)
    do  j=1,ne
       b(j)=h*f(j)
    enddo
    do   j=1,ne
       yt(j)=y(j) + b(j)/2.0_dp
    enddo


    tt=tI+dirb
    call feval(tt,yt,k,f,hcurv,w)
    do  j=1,ne
       c(j)=h*f(j)
    enddo
    do  j=1,ne
       yt(j)=y(j)+c(j)
    enddo

    tt=tI+2*dirb
    call feval(tt,yt,k,f,hcurv,w)
    do  j=1,ne
       d(j)=h*f(j)
    enddo

    do  j=1,ne
       y(j) = y(j)+(a(j)+2.0_dp*b(j)+2.0_dp*c(j)+d(j))/6.0_dp
    enddo
    tI=ti+2*dirb

 !   if(k%TIME) then
 !      Y(6)=Y(6)-(1-k%TOTALPATH)*GR%P%LD/GR%P%beta0/GR%P%nst
 !   else
 !      Y(6)=Y(6)-(1-k%TOTALPATH)*GR%P%LD/GR%P%nst
 !   endif

    call KILL(yt,ne)
    call KILL(f,ne)
    call KILL(a,ne)
    call KILL(b,ne)
    call KILL(c,ne)
    call KILL(d,ne)

    return
  end  subroutine rk4_crittendenp

subroutine integrate_crittendenr(y,w,k)
implicit none
integer nt,is,i
real(dp) h,hcurv
    real(dp), INTENT(INOUT)::  y(6)
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K
    type(work) w

h=zp(2)-zp(1)
hcurv=0.d0
          IS=1
          nt=size(bf,1)
          DO I=1,nt
             call rk4(is,h,hcurv,w,y,k)  
          ENDDO

end subroutine integrate_crittendenr

subroutine integrate_crittendenp(y,w,k)
implicit none
integer nt,is,i
real(dp) h,hcurv
    TYPE(REAL_8), INTENT(INOUT)::  y(6)
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K
    type(work) w

h=zp(2)-zp(1)
hcurv=0.d0
          IS=1
          nt=size(bf,1)
          DO I=1,nt
             call rk4(is,h,hcurv,w,y,k)  
          ENDDO

end subroutine integrate_crittendenp


  SUBROUTINE conv_to_xpr(hcurv,beta0,X,k)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    real(dp) hcurv,beta0
    real(dp) ti
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K
    if(k%TIME) then
       ti=ROOT(1.0_dp+2.0_dp*X(5)/beta0+x(5)**2-X(2)**2-X(4)**2)
       x(2)=(1.0_dp+hcurv*X(1))*x(2)/ti
       x(4)=(1.0_dp+hcurv*X(1))*x(4)/ti
    else
       ti=ROOT((1.0_dp+x(5))**2-X(2)**2-X(4)**2)
       x(2)=(1.0_dp+hcurv*X(1))*x(2)/ti
       x(4)=(1.0_dp+hcurv*X(1))*x(4)/ti
    endif

  end SUBROUTINE conv_to_xpr

  SUBROUTINE conv_to_xpp( hcurv,beta0,X,k)
    IMPLICIT NONE
    type(real_8),INTENT(INOUT):: X(6)
    real(dp) hcurv,beta0
    type(real_8) ti
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K
    call alloc(ti)
    if(k%TIME) then
       ti=sqrt(1.0_dp+2.0_dp*X(5)/beta0+x(5)**2-X(2)**2-X(4)**2)
       x(2)=(1.0_dp+hcurv*X(1))*x(2)/ti
       x(4)=(1.0_dp+hcurv*X(1))*x(4)/ti
    else
       ti=sqrt((1.0_dp+x(5))**2-X(2)**2-X(4)**2)
       x(2)=(1.0_dp+hcurv*X(1))*x(2)/ti
       x(4)=(1.0_dp+hcurv*X(1))*x(4)/ti
    endif
    call kill(ti)

  end SUBROUTINE conv_to_xpp

 SUBROUTINE conv_to_pxr(hcurv,beta0,X,k)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    real(dp) hcurv,beta0
    real(dp) ti
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K
    ti=ROOT((1.0_dp+hcurv*X(1))**2+X(2)**2+X(4)**2)
    if(k%TIME) then
       x(2)=x(2)*ROOT(1.0_dp+2.0_dp*X(5)/beta0+x(5)**2)/ti
       x(4)=x(4)*ROOT(1.0_dp+2.0_dp*X(5)/beta0+x(5)**2)/ti
    else
       x(2)=x(2)*(1.0_dp+x(5))/ti
       x(4)=x(4)*(1.0_dp+x(5))/ti
    endif
  end SUBROUTINE conv_to_pxr

  SUBROUTINE conv_to_pxp(hcurv,beta0,X,k)
    IMPLICIT NONE
    type(real_8),INTENT(INOUT):: X(6)
    real(dp) hcurv,beta0
    type(real_8) ti
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K
    call alloc(ti)
    ti=SQRT((1.0_dp+hcurv*X(1))**2+X(2)**2+X(4)**2)
    if(k%TIME) then
       x(2)=x(2)*sqrt(1.0_dp+2.0_dp*X(5)/beta0+x(5)**2)/ti
       x(4)=x(4)*sqrt(1.0_dp+2.0_dp*X(5)/beta0+x(5)**2)/ti
    else
       x(2)=x(2)*(1.0_dp+x(5))/ti
       x(4)=x(4)*(1.0_dp+x(5))/ti
    endif
    call kill(ti)

  end SUBROUTINE conv_to_pxp



end module shanks_crittenden