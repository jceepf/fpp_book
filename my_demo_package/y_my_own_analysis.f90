module my_analysis
use my_own_da
implicit none
private concat,concatda,real_in_map,map_in_real,concatreal,map_zero,mat_map,map_mat
private unaryminus,POWMAP,POWMAPda,der,input_map_in_vec,input_vec_in_map,subpartmap,subparttaylor
private  comp_in_map,map_in_comp,print_my_vector_field,vector_field_zero,gradt,concatintaylor,texpt,input_map_in_my_taylor
private inverse_dragt_finn_zero,input_idf_map,canonise_map,normal_zero
integer, allocatable :: mres(:) 
real(dp) :: epsresonance = -1.d-10 
logical ::  delta_is_3rd_parameter=.false.
logical ::  pseudo_radiation=.false.
real(dp) :: crad=1.e-4_dp,cfluc=1.e-5_dp  !@ not used yet

 TYPE my_map
       type(my_taylor) v(2)   !@1 &nbsp; good old Taylor map of optics code
 end TYPE my_map
 
type(my_map) to_phasor,from_phasor,c_phasor,ci_phasor !@1 &nbsp; x+ip x-ip basis transformation

 TYPE vector_field         !@1 &nbsp;  map= exp(f.grad)I  
       integer n          !@1 &nbsp;  map= exp(f.grad)I expanded up to order f^n  
       type(my_taylor)v(2)    
 end TYPE vector_field
 
 TYPE inverse_dragt_finn
       type(my_map) m    
       type(vector_field) f    
 end TYPE inverse_dragt_finn
 

 type normalform
       type(my_map) a_t   !@1 &nbsp; Canonical transformation such that Map= at r at^-1
       type(my_map) r    !@1 &nbsp; normal form : rotation usually
       type(my_map) disp    !@1 &nbsp; dispersion transformation
       type(my_map) a_l    !@1 &nbsp; linear part of transformation
       type(my_map) a_nl    !@1 &nbsp; nonlinear part of transformation
       type(my_taylor) total_tune    !@1 &nbsp; total tune
       real(dp) tune,damping   !@1 &nbsp; angle of rotation and damping if nonsymplectic
       real(dp) dtune_dA     !@1 &nbsp; Tune shift with amplitude A=x^2+p^2=2J 
       real(dp) dtune_dk     !@1 &nbsp; Tune shift with parameter 
 end type normalform

  INTERFACE assignment (=)
     MODULE PROCEDURE input_my_taylor_in_map   !@1 &nbsp; map=Taylors
     MODULE PROCEDURE input_map_in_my_taylor   !@1 &nbsp; Taylors=map
     MODULE PROCEDURE input_my_map_in_map   !@1 &nbsp; map=map
     MODULE PROCEDURE input_vec_in_map   !@1 &nbsp; map=vector_field
     MODULE PROCEDURE input_map_in_vec   !@1 &nbsp; vector_field=map
     MODULE PROCEDURE input_idf_map   !@1 &nbsp; inverse dragt finn=map
     MODULE PROCEDURE real_in_map   !@1 &nbsp; map=real
     MODULE PROCEDURE map_in_real   !@1 &nbsp; real=map
     MODULE PROCEDURE comp_in_map  !@1 &nbsp; map=complex
     MODULE PROCEDURE map_in_comp  !@1 &nbsp; complex=map
     MODULE PROCEDURE map_zero   !@1 &nbsp; map=0  or map=identity  (syntax map=0 or map=1)
     MODULE PROCEDURE inverse_dragt_finn_zero !@1 &nbsp; idf=0  or idf=identity  (syntax idf=0 or idf=1)
     MODULE PROCEDURE vector_field_zero   !@1 &nbsp; vector_field=0     
     MODULE PROCEDURE mat_map   !@1 &nbsp;  mat=map
     MODULE PROCEDURE map_mat   !@1 &nbsp;  map=mat
     MODULE PROCEDURE mat_mapc   !@1 &nbsp;  mat=map  : complex matrix
     MODULE PROCEDURE map_matc   !@1 &nbsp;  map=mat  : complex matrix
     MODULE PROCEDURE normal_zero !@1 &nbsp;  zeroes normal form (bug found with gfortran)
     
  end  INTERFACE  
  
    INTERFACE OPERATOR (.o.)
!@3  <b><font color="#FF0000">map=map.o.map is interesting </font></b> click below </br>
     MODULE PROCEDURE concat   !@1 &nbsp; map 0 map   
     MODULE PROCEDURE concatreal   !@1 &nbsp; map 0 real(2)   
  END INTERFACE 
  
      INTERFACE OPERATOR (**) !@2 <font color="#FF0000"><b>Constant part is ignored in a &quot;DA&quot; calculation</b></font>
        MODULE PROCEDURE  POWMAPDA  !@1 &nbsp; map is raised to a power; can be negative (inverse)
      END INTERFACE 
      INTERFACE OPERATOR (.oo.) !@2 <font color="#FF0000"><b>Constant part is used in a &quot;TPSA&quot; calculation</b></font>
        MODULE PROCEDURE  POWMAP  !@1 &nbsp; map is raised to a power; can be negative (inverse)
      END INTERFACE 

    INTERFACE OPERATOR (*)
!@3  <b><font color="#FF0000">map=map*map is interesting </font></b> click below </br>
     MODULE PROCEDURE concatda   !@1 &nbsp; map * map  
     MODULE PROCEDURE concatintaylor   !@1 &nbsp; taylor o map
     MODULE PROCEDURE concatrealda   !@1 &nbsp; map 0 real(2) 
  END INTERFACE

  INTERFACE OPERATOR (.d.)
     MODULE PROCEDURE der   !@1 &nbsp;  derivative  
  END INTERFACE

  INTERFACE OPERATOR (.grad.)
     MODULE PROCEDURE gradt   !@1 &nbsp;  <font color="#FF0000"><b>F<font face="Symbol">Å~Å˜</font>M&nbsp; where F is type vector_field and M is type my_map</b></font>
  END INTERFACE

  INTERFACE OPERATOR (.subpart.)
     MODULE PROCEDURE subparttaylor   !@1 &nbsp;  gets the <=ith part of the taylor  
     MODULE PROCEDURE subpartmap   !@1 &nbsp;  gets the <=ith part of the map  
  END INTERFACE

  INTERFACE OPERATOR (.i.)
     MODULE PROCEDURE subtaylor   !@1 &nbsp;  gets the =ith part of the taylor  
     MODULE PROCEDURE submap   !@1 &nbsp;  gets the =ith part of the map  
  END INTERFACE

  INTERFACE OPERATOR (.vf.)
     MODULE PROCEDURE extract_leading_vecfield  !@1 &nbsp;  gets the ith part of the map and put it in a vector field
  END INTERFACE

  INTERFACE OPERATOR (-)
     MODULE PROCEDURE unaryminus   !@1 &nbsp;  -map  
  END INTERFACE

  INTERFACE texp
     MODULE PROCEDURE texpt   !@1 &nbsp; <font color="#FF0000"><b>exp(F<font face="Symbol">Å~Å˜</font>)M where F is type vector_field and M is type my_map</b></font></font>  
  END INTERFACE

  INTERFACE exp
     MODULE PROCEDURE texpt !@1 &nbsp; <font color="#FF0000"><b>exp(F<font face="Symbol">Å~Å˜</font>)M where F is type vector_field and M is type my_map</b></font></font>  
  END INTERFACE

  INTERFACE canonise
     MODULE PROCEDURE canonise_map   !@1 &nbsp; Puts the canonical transformation into a special form
  END INTERFACE


  
  INTERFACE print
     MODULE PROCEDURE print_my_map   !@1 &nbsp; prints map for human eye
     module procedure print_my_vector_field !@1 &nbsp; prints a vector field for human eye
  END INTERFACE

  INTERFACE print_for_human
     MODULE PROCEDURE print_for_human_map   !@1 &nbsp; prints map for human eye
     module procedure print_for_human_my_vector_field !@1 &nbsp; prints a vector field for human eye
  END INTERFACE



contains

  subroutine map_mat( S2, S1 )
    implicit none
    real(dp), INTENT (IN) :: S1(3,3)
    TYPE (my_map), INTENT (inout) :: S2
    integer i,j 
      s2=0
      
      do i=1,2
      do j=1,3
       S2%v(i)%a(j)=s1(i,j)
      enddo
      enddo
  END subroutine map_mat
  
  subroutine mat_map( S2, S1 )
    implicit none
    real(dp), INTENT (INout) :: S2(3,3)
    TYPE (my_map), INTENT (in) :: S1
    integer i,j 
      
      s2=0.d0
      s2(3,3)=1
      do i=1,2
      do j=1,2
       s2(i,j)=S1%v(i)%a(j)
      enddo
      enddo

  END subroutine mat_map
  
  subroutine map_matc( S2, S1 )
    implicit none
    complex(dp), INTENT (IN) :: S1(2,2)
    TYPE (my_map), INTENT (inout) :: S2
    integer i,j 
      s2=0
      
      do i=1,2
      do j=1,2
       S2%v(i)%a(j)=s1(i,j)
      enddo
      enddo

  END subroutine map_matc
  
  subroutine mat_mapc( S2, S1 )
    implicit none
    complex(dp), INTENT (INout) :: S2(2,2)
    TYPE (my_map), INTENT (in) :: S1
    integer i,j 
      
      s2=0.d0
      
      do i=1,2
      do j=1,2
       s2(i,j)=S1%v(i)%a(j)
      enddo
      enddo

  END subroutine mat_mapc


  function texpt(s2, S1 )
    implicit none
    TYPE (my_map) texpt
    TYPE (my_map),INTENT (IN) :: S1
    TYPE (vector_field),INTENT (IN) :: S2
    TYPE (my_map) t,sf
    integer i,j
    
    t=s1
    sf=s1

    do i=1,s2%n
     t=s2.grad.t
     do j=1,2
      t%v(j)=t%v(j)/REAL(i,kind=DP)   ! creating the exponential
      sf%v(j)=sf%v(j)+t%v(j)
     enddo
    enddo
    
    texpt=sf

  END function texpt

  function gradt(s2, S1 )
    implicit none
    TYPE (my_map) gradt
    TYPE (my_map),INTENT (IN) :: S1
    TYPE (vector_field),INTENT (IN) :: S2
    integer i,j
    TYPE (my_map) t
    
    t=0
    
    do j=1,2
    do i=1,2
     t%v(j)= s2%v(i)*(s1%v(j).d.i) + t%v(j)
    enddo
    enddo
    
    gradt=t
    
  END function gradt


  
  function der(s2, S1 )
    implicit none
    TYPE (my_taylor) der
    integer, INTENT (IN) :: S1
    TYPE (my_taylor),INTENT (IN) :: S2
    integer i,k

     if(first_mul) then
      call multiplication_table
      first_mul=.false.
     endif

      der=0.0_dp 
      
      if(s1==1) then
       do i=1,n_mono
        k=der_table(i,1)
        if(k>=0.and.jorder(k)<=my_order) then
         der%a(k)=jexp1(i)*s2%a(i)
        endif
       enddo
      endif

      if(s1==2) then
       do i=1,n_mono
        k=der_table(i,2)
        if(k>=0.and.jorder(k)<=my_order) then
         der%a(k)=jexp2(i)*s2%a(i)
        endif
       enddo
      endif
      
      if(s1==3) then
       do i=1,n_mono
        k=der_table(i,3)
        if(k>=0.and.jorder(k)<=my_order) then
         der%a(k)=jexp3(i)*s2%a(i)
        endif
       enddo
      endif

  END function der

  function subparttaylor(s2, S1 )
    implicit none
    TYPE (my_taylor) subparttaylor
    integer, INTENT (IN) :: S1
    TYPE (my_taylor),INTENT (IN) :: S2
    integer i
    
    
      subparttaylor=0.0_dp
      
      do i=0,n_mono
       if(jorder(i)<=s1) subparttaylor%a(i)=s2%a(i)
      enddo 


  END function subparttaylor

  function subpartmap(s2, S1 )
    implicit none
    TYPE (my_map) subpartmap
    integer, INTENT (IN) :: S1
    TYPE (my_map),INTENT (IN) :: S2
    
    subpartmap%v(1)=s2%v(1).subpart.s1
    subpartmap%v(2)=s2%v(2).subpart.s1
    
  END function subpartmap

  function subtaylor(s2, S1 )
    implicit none
    TYPE (my_taylor) subtaylor
    integer, INTENT (IN) :: S1
    TYPE (my_taylor),INTENT (IN) :: S2
    integer i
    
    
      subtaylor=0.0_dp
      
      do i=0,n_mono
       if(jorder(i)==s1) subtaylor%a(i)=s2%a(i)
      enddo 


  END function subtaylor

  function submap(s2, S1 )
    implicit none
    TYPE (my_map) submap
    integer, INTENT (IN) :: S1
    TYPE (my_map),INTENT (IN) :: S2
    
    submap%v(1)=s2%v(1).i.s1
    submap%v(2)=s2%v(2).i.s1
    
  END function submap

  function extract_leading_vecfield(s2, S1 )
    implicit none
    TYPE (vector_field) extract_leading_vecfield
    integer, INTENT (IN) :: S1
    TYPE (my_map),INTENT (IN) :: S2
   
    extract_leading_vecfield%v(1)=s2%v(1).i.s1
    extract_leading_vecfield%v(2)=s2%v(2).i.s1
    
  END function extract_leading_vecfield


  subroutine map_zero( S2, S1 )
    implicit none
    integer, INTENT (IN) :: S1
    TYPE (my_map), INTENT (inout) :: S2
     
      S2%v(1)=0.d0
      S2%v(2)=0.d0
      
      if(s1==1) then
       S2%v(1)%a(1)=1.d0
       S2%v(2)%a(2)=1.d0
      endif

  END subroutine map_zero

  subroutine normal_zero( S2, S1 )
    implicit none
    integer, INTENT (IN) :: S1
    TYPE (normalform), INTENT (inout) :: S2
     
      s2%a_t=0  
      s2%r=0    
      s2%disp=0 
      s2%a_l=0  
      s2%a_nl=0 
      s2%tune=0.d0
      s2%damping=0.d0
       s2%dtune_dA=0.d0    
       s2%dtune_dk=0.d0    
       s2%total_tune=0.d0

  END subroutine normal_zero
  

    subroutine inverse_dragt_finn_zero( S2, S1 )
    implicit none
    integer, INTENT (IN) :: S1
    TYPE (inverse_dragt_finn), INTENT (inout) :: S2
     
      s2%m=s1
      s2%f=s1

  END subroutine inverse_dragt_finn_zero

  subroutine vector_field_zero( S2, S1 )
    implicit none
    integer, INTENT (IN) :: S1
    TYPE (vector_field), INTENT (inout) :: S2
    integer power
      S2%v(1)=0.0_dp
      S2%v(2)=0.0_dp
      
      power=s1
      if(power<=my_order) power=my_order+1
      s2%n=power

  END subroutine vector_field_zero

  subroutine input_idf_map( S2, S1 )
    implicit none
    TYPE (my_map), INTENT (IN) :: S1
    TYPE (inverse_dragt_finn), INTENT (inout) :: S2
    TYPE (my_map) m,DF,m0
    integer i
    
    s2=0
    s2%m=s1.subpart.1
    
    m=s2%m**(-1)*s1
   m0=m
    
    do i=1,my_order
    
     DF=m.subpart.1
     s2%f%V(1)=-(M%V(1)-DF%V(1))+s2%f%V(1)
     s2%f%V(2)=-(M%V(2)-DF%V(2))+s2%f%V(2)
       m=exp(s2%f,m0)
    
    enddo
     
     s2%f%V(1)=-s2%f%V(1)
     s2%f%V(2)=-s2%f%V(2)
    

    
  END subroutine input_idf_map



  
  subroutine input_my_taylor_in_map( S2, S1 )
    implicit none
    TYPE (my_taylor), INTENT (IN) :: S1(2)
    TYPE (my_map), INTENT (inout) :: S2

     S2%v(1)=S1(1)
     S2%v(2)=S1(2)
     
  END subroutine input_my_taylor_in_map
  
  subroutine input_map_in_my_taylor( S2, S1 )
    implicit none
    TYPE (my_taylor), INTENT (inout) :: S2(2)
    TYPE (my_map), INTENT (in) :: S1

     S2(1)=S1%v(1)
     S2(2)=S1%v(2)
     
  END subroutine input_map_in_my_taylor
  
  subroutine input_my_map_in_map( S2, S1 )
    implicit none
    TYPE (my_map), INTENT (IN) :: S1 
    TYPE (my_map), INTENT (inout) :: S2

     S2%v(1)=S1%v(1)
     S2%v(2)=S1%v(2)
     
  END subroutine input_my_map_in_map

  subroutine input_vec_in_map( S2, S1 )
    implicit none
    TYPE (vector_field), INTENT (IN) :: S1 
    TYPE (my_map), INTENT (inout) :: S2

     S2%v(1)=S1%v(1)
     S2%v(2)=S1%v(2)
     
  END subroutine input_vec_in_map

  subroutine input_map_in_vec( S2, S1 )
    implicit none
    TYPE (my_map), INTENT (IN) :: S1 
    TYPE (vector_field), INTENT (inout) :: S2

     S2%v(1)=S1%v(1)
     S2%v(2)=S1%v(2)
     
  END subroutine input_map_in_vec

  
   subroutine real_in_map( S2, S1 )
    implicit none
    real(dp), INTENT (IN) :: S1(2)
    TYPE (my_map), INTENT (inout) :: S2

     S2%v(1)=S1(1)
     S2%v(2)=S1(2)
     
  END subroutine real_in_map
 
  subroutine map_in_real( S2, S1 )
    implicit none
    real(dp), INTENT (INOUT) :: S2(2)
    TYPE (my_map), INTENT (in) :: S1

     S2(1)=S1%v(1)
     S2(2)=S1%v(2)
     
  END subroutine map_in_real
 
    subroutine comp_in_map( S2, S1 )
    implicit none
    complex(dp), INTENT (IN) :: S1(2)
    TYPE (my_map), INTENT (inout) :: S2

     S2%v(1)=S1(1)
     S2%v(2)=S1(2)
     
  END subroutine comp_in_map

  subroutine map_in_comp( S2, S1 )
    implicit none
    complex(dp), INTENT (INOUT) :: S2(2)
    TYPE (my_map), INTENT (in) :: S1

     S2(1)=S1%v(1)
     S2(2)=S1%v(2)
     
  END subroutine map_in_comp
 
  subroutine print_my_map( S1, mfi,title )
    implicit none
    TYPE (my_map), INTENT (in) :: S1
    TYPE (my_map) s2
    integer mf
	logical fi
     character*3 fii
	integer,optional :: mfi
    character(*) , optional :: title
     mf=6
     if(present(mfi)) mf=mfi


    s2=s1
    write(mf,*) "  "
    write(mf,*) " variable 1 of map"
     call print(s2%v(1),mf,title)
    write(mf,*) "  "
    write(mf,*) " variable 2 of map"
     call print(s2%v(2),mf,title)

  END subroutine print_my_map
  
  subroutine print_my_map2( S1,s1p, mfi,title )
    implicit none
    TYPE (my_map), INTENT (in) :: S1
    TYPE (my_map),optional, INTENT (in) :: s1p
    TYPE (my_map) s2,s2p
    integer mf
	logical fi
     character*3 fii
	integer,optional :: mfi
    character(*) , optional :: title
     mf=6
     if(present(mfi)) mf=mfi


    s2=s1
    write(mf,*) "  "
     if(present(s1p)) then 
      s2p=s1p
     call print_my_taylor4(s2%v(1),s2%v(2),s2p%v(1),s2p%v(2),mf,title)
      else
     call print_my_taylor2(s2%v(1),s2%v(2),mf,title)
     endif

  END subroutine print_my_map2



  subroutine print_my_vector_field( S2, mfi,title )
    implicit none
    TYPE (vector_field), INTENT (inout) :: S2
    integer mf
	logical fi
     character*3 fii
	integer,optional :: mfi
    character(*) , optional :: title
     mf=6
     if(present(mfi)) mf=mfi
    write(mf,*) "  "

    write(mf,*) " variable 1 of vector field"
     call print(s2%v(1),mf,title)
    write(mf,*) "  "
     write(mf,*) " variable 2 of vector field"
     call print(s2%v(2),mf,title)

  END subroutine print_my_vector_field
  

  subroutine print_for_human_map( S1, mfi,title )
    implicit none
    TYPE (my_map), INTENT (in) :: S1
    TYPE (my_map) s2
    integer mf
	logical fi
     character*3 fii
	integer,optional :: mfi
    character(*) , optional :: title
     mf=6
     if(present(mfi)) mf=mfi


    s2=s1
    write(mf,*) "  "
    write(mf,*) " variable 1 of map"
     call print_for_human(s2%v(1),mf,title)
    write(mf,*) "  "
    write(mf,*) " variable 2 of map"
     call print_for_human(s2%v(2),mf,title)

  END subroutine print_for_human_map

  subroutine print_for_human_my_vector_field( S2, mfi,title )
    implicit none
    TYPE (vector_field), INTENT (inout) :: S2
    integer mf
	logical fi
     character*3 fii
	integer,optional :: mfi
    character(*) , optional :: title
     mf=6
     if(present(mfi)) mf=mfi
    write(mf,*) "  "

    write(mf,*) " variable 1 of vector field"
     call print_for_human(s2%v(1),mf,title)
    write(mf,*) "  "
     write(mf,*) " variable 2 of vector field"
     call print_for_human(s2%v(2),mf,title)

  END subroutine print_for_human_my_vector_field

  FUNCTION unaryminus( S2 )
    implicit none
     TYPE (my_map) unaryminus
     TYPE (my_map), INTENT (in) :: S2
     integer i
    
     unaryminus%v(1)= -S2%v(1)
     unaryminus%v(2)= -S2%v(2)
     
  END FUNCTION unaryminus
 


  FUNCTION concat( S1, S2 )
    implicit none
    TYPE (my_map) concat
    TYPE (my_map), INTENT (IN) :: S1, S2
    TYPE (my_taylor)  x(3),y(2),z(2)
    integer i,j,k
    x(1)=s2%v(1)
    x(2)=s2%v(2)
    x(3)=1.d0.mon.3
    y(1)=s1%v(1)
    y(2)=s1%v(2)
    
    z(1)=0.d0
    z(2)=0.d0
!@3 <p><font color="#FF0000"><b>New Code using tables to 4<sup>th </sup>order</b></font></p>
!@3 <p><font color="#FF0000"><b> This is ABSOLUTELY necessary for speed</b></font></p>
    
    table(0,0,0)=1.0_dp
    table(1,0,0)=x(1)
    table(0,1,0)=x(2)
    table(0,0,1)=x(3)
    table(2,0,0)=x(1)**2
    table(0,2,0)=x(2)**2
    table(0,0,2)=x(3)**2
    table(1,1,0)=x(1)*x(2)
    table(1,0,1)=x(1)*x(3)
    table(0,1,1)=x(2)*x(3)
    table(3,0,0)=table(1,0,0)*table(2,0,0)
    table(0,3,0)=table(0,1,0)*table(0,2,0)
    table(0,0,3)=table(0,0,1)*table(0,0,2)
    table(2,1,0)=table(2,0,0)*table(0,1,0)
    table(2,0,1)=table(2,0,0)*table(0,0,1)
    table(0,2,1)=table(0,2,0)*table(0,0,1)
    table(1,2,0)=table(0,2,0)*table(1,0,0)
    table(1,0,2)=table(0,0,2)*table(1,0,0)
    table(0,1,2)=table(0,0,2)*table(0,1,0)
    table(1,1,1)=x(1)*table(0,1,1)
    
    table(4,0,0)=table(2,0,0)*table(2,0,0)
    table(0,4,0)=table(0,2,0)*table(0,2,0)
    table(0,0,4)=table(0,0,2)*table(0,0,2)
    table(3,1,0)=table(3,0,0)*table(0,1,0)
    table(3,0,1)=table(3,0,0)*table(0,0,1)
    table(0,3,1)=table(0,3,0)*table(0,0,1)
    table(1,3,0)=table(0,3,0)*table(1,0,0)
    table(0,1,3)=table(0,0,3)*table(0,1,0)
    table(1,0,3)=table(0,0,3)*table(1,0,0)
    table(1,1,2)=table(0,0,2)*table(1,1,0)
    table(1,2,1)=table(0,2,0)*table(1,0,1)
    table(2,1,1)=table(2,0,0)*table(0,1,1)
    table(2,2,0)=table(2,0,0)*table(0,2,0)
    table(0,2,2)=table(0,2,0)*table(0,0,2)
    table(2,0,2)=table(2,0,0)*table(0,0,2)

    do i=1,2
      z(i)=0.d0  ! y(i)%a(0) bug found by Lingyun Yang 2013.04.26
     do j=0,n_mono
      if(jorder(j)>my_order) cycle
       z(i)=z(i)+ y(i)%a(j)*table(jexp1(j),jexp2(j),jexp3(j))
     enddo 
    enddo
 !@3 <p><font color="#FF0000"><b> Switch the above lines and notice the speed degradation</b></font></p>
 
    concat%v(1)=z(1)
    concat%v(2)=z(2)
     
  END FUNCTION concat


  FUNCTION concatda( S1, S2 )
    implicit none
    TYPE (my_map) concatda
    TYPE (my_map), INTENT (IN) :: S1, S2
    TYPE (my_map)   S20
      
    s20=s2
    s20%v(1)%a(0)=0.d0
    s20%v(2)%a(0)=0.d0
 

      concatda=S1.O.S20

 
  END FUNCTION concatda

  FUNCTION concatintaylor( T, S2 )
    implicit none
    TYPE (my_taylor) concatintaylor
    TYPE (my_taylor), INTENT (IN) :: t
    TYPE (my_map), INTENT (IN) :: S2
    TYPE (my_map)   S20
      
    s20=0
    
    s20%v(1)=t
     
     s20=S20*S2

      concatintaylor=s20%v(1)

 
  END FUNCTION concatintaylor



  FUNCTION concatreal( S1, z )
    implicit none
    real(dp)  concatreal(2)
    TYPE (my_map), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: z(2)
    TYPE (my_map) zm
    
    zm=z

    zm=s1.o.zm
 
    concatreal=zm
     
  END FUNCTION concatreal


  FUNCTION concatrealda( S1, z )
    implicit none
    real(dp)  concatrealda(2)
    TYPE (my_map), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: z(2)
    TYPE (my_map) zm,s10
    
    s10=s1
    s10%v(1)%a(0)=0.d0
    s10%v(2)%a(0)=0.d0
    
    zm=z

    zm=s10.o.zm
 
    concatrealda=zm
     
  END FUNCTION concatrealda



  subroutine invert_map( S1, S2 )
    implicit none
    TYPE (my_map), INTENT (IN) :: S1 
    TYPE (my_map), INTENT (inout) :: S2
    TYPE (my_map) at,s10,a
    complex(dp) det,mi(2,2),d(2),w(2),v(2)
    integer i,ex(3)
    
    mi=s1
    call invmatc( mi,mi )
!    mi(1,1)=s1%v(2)%a(2)
!    mi(2,2)=s1%v(1)%a(1)
!    mi(1,2)=-s1%v(1)%a(2)
!    mi(2,1)=-s1%v(2)%a(1)
!    det=mi(1,1)*mi(2,2)-mi(1,2)*mi(2,1)
!    mi=mi/det
    
    
    s10=s1
    s10%v(1)%a(0)=0.d0
    s10%v(2)%a(0)=0.d0
    
    at=0   
    at=mi
    
    ex=0;ex(3)=1;

    v(1)=s1%v(1).sub.ex
    v(2)=s1%v(2).sub.ex
    w(1)=-mi(1,1)*v(1)-mi(1,2)*v(2)
    w(2)=-mi(2,1)*v(1)-mi(2,2)*v(2)
    at%v(1)=at%v(1)+w(1)*(1.0_dp.mon.3)
    at%v(2)=at%v(2)+w(2)*(1.0_dp.mon.3)
    
    
    a=at
    s10=at*s10
    do i=2, my_order
    at=s10
     at%v(1)%a(1:2)= -at%v(1)%a(1:2)
     at%v(2)%a(1:2)= -at%v(2)%a(1:2)
     at=-at
     s10=at*s10
     a=at*a    
    enddo
    d=S1
    at=1
    at%v(1)=at%v(1)-d(1)
    at%v(2)=at%v(2)-d(2)
    s2=a.o.at
  END subroutine invert_map

  FUNCTION POWMAP( S1, R2 )
    implicit none
    TYPE (my_map) POWMAP
    TYPE (my_map), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: R2
    TYPE (my_map) S11
    INTEGER I,R22

    s11=1


    R22=IABS(R2)
    DO I=1,R22
       s11=s1*s11
    ENDDO

    IF(R2.LT.0) THEN
       CALL invert_map( S11, S11 ) 
    ENDIF

    powmap=s11

  END FUNCTION POWMAP

  FUNCTION POWMAPda( S1, R2 )
    implicit none
    TYPE (my_map) POWMAPda
    TYPE (my_map), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: R2
    TYPE (my_map) S11,S10
    INTEGER I,R22

    s11=1
    s10=S1
    S10%V(1)%A(0)=0.D0
    S10%V(2)%A(0)=0.D0

    R22=IABS(R2)
    DO I=1,R22
       s11=s10*s11
    ENDDO

    IF(R2.LT.0) THEN
       CALL invert_map( S11, S11 ) 
    ENDIF
    POWMAPda=s11

  END FUNCTION POWMAPda

  subroutine compute_dispersion( A,disp )
    implicit none
    TYPE (my_map) , INTENT (inout) ::   A,disp 

    disp=0
    
    disp=a*disp    !@1 &nbsp; removes transverse parts
    
    disp%v(1)=(1.0_dp.mon.1)+disp%v(1) !@1 &nbsp; Adds identity
    disp%v(2)=(1.0_dp.mon.2)+disp%v(2) !@1 &nbsp; Adds identity
  
  end subroutine compute_dispersion
      
  subroutine canonise_map( Atot,A_cs,disp0,a_l0,a_nl0,R0,PHASE_ADVANCE0 )
    implicit none
    TYPE (my_map) , INTENT (inout) ::   Atot,A_cs
     TYPE (my_map) ,optional,  INTENT (inout) ::  R0,disp0,A_l0,A_nl0
    TYPE(MY_TAYLOR), OPTIONAL, INTENT (inout) ::  PHASE_ADVANCE0
    TYPE (my_map) RT,AN,A, R,disp,A_l,A_nl
    TYPE(MY_TAYLOR) a11,a12,D_MU,C,S,PHASE_ADVANCE
    type(inverse_dragt_finn) idf,idfn
    type(vector_field) f
    integer i1(3),i2(3),i

    PHASE_ADVANCE=0.0_DP
    IF(PRESENT(PHASE_ADVANCE0)) PHASE_ADVANCE=PHASE_ADVANCE0
    A=Atot
    !@2 Compute_dispersion gives us the parameter dependent fixed point.
    
    call compute_dispersion( A,disp ) !@1 &nbsp; Parameter dependent dispersion
        
    a_cs=disp**(-1)*A !@1 &nbsp; a_cs diagonolises a map that is already around its parameter dependent fixed point 
    
    call factor_linear(a_cs,A_l,A_nl ) !@1 &nbsp; a_cs = A_l o A_nl
    
    call  get_a11_a12( A_l,a11,a12 )  !@1 &nbsp; (A_11,A_12) are computed as functions of the parameter
    


    
    R=1
    D_MU=ATAN2(a12,a11)
    d_mu=d_mu.subpart.(my_order-1)

    C=COS(D_MU)
    S=SIN(D_MU)
    PHASE_ADVANCE=PHASE_ADVANCE+D_MU/TWOPI
    
    RT=0    
    RT%V(1)=C*(1.0_dp.mon.1)+S*(1.0_dp.mon.2)
    RT%V(2)=C*(1.0_dp.mon.2)-S*(1.0_dp.mon.1)
    
    A_l=A_l*RT**(-1)
    A_nl=RT*A_nl*RT**(-1)   !@1 &nbsp; A_l and A_nl  are rotated so that A_l_12=0

    R=RT*R

!@2 A canonical rule is added for the nonlinear part of the map; 
!@2 rotation generators are removed   

    call create_phasors
    
    idf=A_nl
    AN=1
    AN=exp(idf%f,AN)
    AN=to_phasor *AN*from_phasor
    idfn=an
    
    i1=0;i1(1)=2;i1(2)=1;
    i2=0;i2(1)=1;i2(2)=2;
    
    f=0
    do i=0,n_mono 
     if(jexp1(i)==2.and.jexp2(i)==1) then !@1 &nbsp; x_1^2*x_2 d/dx_2 generates rotation
      f%v(1)%a(i)=idfn%f%v(1)%a(i)
      PHASE_ADVANCE=PHASE_ADVANCE- &
     ((1.0_DP.mon.1)**2+(1.0_DP.mon.2)**2)*IMAG(idfn%f%v(1)%a(i)) &
     *(1.0_dp.mon.3)**jexp3(i)/TWOPI
     elseif(jexp1(i)==1.and.jexp2(i)==2) then !@1 &nbsp; x_1*x_2^2 d/dx_1 generates rotation
      f%v(2)%a(i)=idfn%f%v(2)%a(i)
     endif
    enddo
    RT=1

    RT=exp(f,RT)
    RT=from_phasor*RT*to_phasor
    A_nl=A_nl*rt**(-1)
    R=RT*R
    
    A_cs=disp*a_l*a_nl !@1 &nbsp; Final form of A_cs &nbsp; A = A_cs o R 
    
    if(present(disp0) ) disp0=disp
    if(present(r0) )    r0=r
    if(present(a_nl0) ) a_nl0=a_nl
    if(present(a_l0) ) a_l0=a_l
    IF(PRESENT(PHASE_ADVANCE0)) PHASE_ADVANCE0=PHASE_ADVANCE
  END subroutine canonise_map

  
  subroutine factor_linear( A,A_l,A_nl )
    implicit none
    TYPE (my_map) , INTENT (inout) ::   A,A_l,A_nl 
    integer i,j
     
    A_l=0
    do i=0,N_mono
      if(jexp1(i)+jexp2(i)==1) then !@1 &nbsp; x_1^i x_2^j &nbsp;  i+j=1 implies linear!
       do j=1,2
        A_l%v(j)%a(i)=a%v(j)%a(i)
       enddo
      endif
    enddo
    a_nl=a_l**(-1)*a !@1 &nbsp; a_nl is a fully nonlinear map including parameter dependent terms
    
  end subroutine factor_linear
   
  subroutine get_a11_a12( A_l,a11,a12 )
    implicit none
    TYPE (my_map) , INTENT (inout) ::   A_l
    type(my_taylor) , intent(inout):: a11,a12 
    integer i
     
    a11=0.0_dp
    a12=0.0_dp
    do i=0,N_mono
      if(jexp1(i)==1.and.jexp2(i)==0) then
       a11%a(i)=a11%a(i)+ A_l%v(1)%a(i)
      endif
      if(jexp1(i)==0.and.jexp2(i)==1) then
       a12%a(i)=a12%a(i)+ A_l%v(1)%a(i)
      endif
    enddo
    a11=a11.d.1
    a12=a12.d.2
    
  end subroutine get_a11_a12
   
  subroutine normalise( M,N )
    implicit none
    TYPE (normalform), INTENT (inout) ::N
    TYPE (my_map), INTENT (in) :: M
    real(dp) det,c,s,ang,alpha,beta
    TYPE (my_map) A0,A1,R,NL,nlc,b
    type(vector_field) f
    type(my_taylor) junk_phase
    integer k,j(3)

     N=0
    call create_phasors  !@3 &nbsp; transformation to a complex phasors basis is computed

!@3 !!!!  Linear Analysis starts here
    call find_disp(m,a0)  !@3 &nbsp; Find closed orbit to first order in parameter

    n%a_t=a0
!@1 !!!!  The map is put around the liner fixed point

    N%R=A0**(-1)*m*A0
    
    
  !@1 !!!!  The linear Courant-Snyder map's computation begins!
    call diag_mat(N%R,a1,n%tune,n%damping) 

  !@1 !!!!  The linear Courant-Snyder map's computation is now done!
  
    n%a_t=n%a_t*a1
    N%R=A1**(-1)*N%R*A1
   

    R=N%R.i.1  !@3 &nbsp; LINEAR PART

 !@3 !!!! Linear normalisation done
 
    NL=N%R*R**(-1)   !@3 &nbsp; NONLINEAR PART   N%R= NL * R 

        
    do k=1,my_order-1
  !@1 The order by order normalization into a rotation begins
     call analyse_kernel(k+1,n%tune,n%damping,nl,f,n%dtune_dA,n%dtune_dk)
     b=1    


     b=exp(f,b)

     b=from_phasor *b*to_phasor

     n%r=b**(-1)*n%r*b

     NL=N%R*R**(-1)   !@3 &nbsp; NONLINEAR PART   N%R= NL * R 
     n%a_t=n%a_t*b

    enddo
   !@1 The map n%a_t is put into a canonical form 
    call canonise( n%a_t,n%a_t,n%disp,n%a_l,n%a_nl,R,junk_phase )


    !@1 The map n%a_t*R is put into a canonical phase just to get the tunes!
     call canonise( n%R,a1,a0,nl,nlc,R,n%total_tune )

 

     
     det=n%total_tune
     if(det<0) n%total_tune%a(0)=det+1.0_dp
     
     j=0
     j(1)=2
     n%dtune_dA=n%total_tune.sub.j
     j=0
     j(3)=1
     n%dtune_dk=n%total_tune.sub.j
   END subroutine normalise
   
    
    
  subroutine analyse_kernel(order,tune,alpha,nl,g,dtune_dA,dtune_dk)
    implicit none
    TYPE (my_map), INTENT (inout) :: nl
    type(vector_field), INTENT (inout) ::  g
    TYPE (my_map) t
    integer i,j,k,order,i1,i2
    type(vector_field) f
    real(dp) dtune_dA,dtune_dk,tune,alpha
    complex(dp) lam(2),denominator

    lam(1)=exp(-i_*tune*twopi+alpha)
    lam(2)=exp(i_*tune*twopi+alpha)
    
    g=0
    
     t=to_phasor*nl*from_phasor 
     f=t.vf.order
     
 do i=1,2
 do j=0,n_mono    
     if(jorder(j)/=order) cycle
      denominator= 1.0_dp-lam(i)*exp(i_*(jexp1(j)-jexp2(j))*tune*twopi-(jexp1(j)+jexp2(j))*alpha)
      
      if(epsresonance>0) then    !@1 &nbsp; determined by closeness to a resonance
          if(abs(denominator)>epsresonance) then
              g%v(i)%a(j) = f%v(i)%a(j)/denominator 
          endif 
      else             !@1 &nbsp; leaves tune shifts and amplitude dependent damping
          if((jexp1(j)-jexp2(j)-1/=0.and.i==1).or.(jexp1(j)-jexp2(j)+1/=0.and.i==2)) then
              g%v(i)%a(j) = f%v(i)%a(j)/denominator 
          endif
      endif 
      
  !    if(allocated(mres)) then
  !    do k=1,size(mres)   !@ &nbsp; undoes the operation above and leaves the resonance
  !        if((abs(jexp1(j)-jexp2(j)-1)==mres(k).and.i==1).or.(abs(jexp1(j)-jexp2(j)+1)==mres(k).and.i==2)) then
  !            g%v(i)%a(j) = 0.0_dp 
  !        endif
  !    enddo
  !    endif
           
 enddo
 enddo

  END subroutine analyse_kernel
  
   subroutine find_disp( m,a0 )
    implicit none
    TYPE (my_map), INTENT (in) :: m
    TYPE (my_map), INTENT (inout) :: a0
    complex(dp) mat(2,2),det,eta(2),v(2)
    integer i,ex(3)
     ex=0
     ex(3)=1
   
    mat=m
    do i=1,2
     v(i)=m%v(i).sub.ex
    enddo
    
    do i=1,2
     mat(i,i)=mat(i,i)-1.0_dp
    enddo
    
    mat=-mat
    
    call invmatc( mat,mat )
    
    
    eta(1)=mat(1,1)*v(1)+mat(1,2)*v(2)
    eta(2)=mat(2,1)*v(1)+mat(2,2)*v(2)

    a0=1
    a0%v(1)=a0%v(1)+eta(1)*(1.0_dp.mon.3)
    a0%v(2)=a0%v(2)+eta(2)*(1.0_dp.mon.3)
    
    end subroutine find_disp 
    
   subroutine diag_mat( m,a1,tune,dam )
    implicit none
    TYPE (my_map), INTENT (inout) :: m
    TYPE (my_map), INTENT (inout) :: a1
    real(dp) mat(3,3)
    complex(dp) eta(2),v(2),lam,trace,del
    real(dp) det,x(2),px(2),norm,norm1,norm2,tune,dam
    
  !  type(my_map) A_cs,disp,a_l,a_nl,R
  !  type(my_taylor) PHASE_ADVANCE
!@1 first we find the eigenvectors which we assume here to be complex (stable) 

    mat=m
    
    det=mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1)    
    trace=mat(1,1)+mat(2,2)
    
    del=i_*sqrt(abs(trace**2-4.0_dp*det))
    lam=(trace+del)/2.0_dp
    
!@1 we need an eigenvector of the transpose of mat to construct A    
    v(1)=1.0_dp
    v(2)=(lam-mat(1,1))/mat(2,1)
    
!@1  Insuring that    Poisson bracket should be equal to one
    x(1)=real(v(1))   ;  x(2)=real(v(2))   ;
    px(1)=aimag(v(1)) ;  px(2)=aimag(v(2)) ;
    norm=x(1)*px(2)-x(2)*px(1)
    if(norm<0) then
     norm1= sqrt(abs(norm)); norm2 =-norm1;
    else
     norm1= sqrt(abs(norm)); norm2 = norm1;
    endif
    x=x/norm1
    px=px/norm2
!@1 Now the Poisson bracket should be equal to one    

    a1%v(1)= ( x(1).mon.1)+( x(2).mon.2)
    a1%v(2)= (px(1).mon.1)+(px(2).mon.2)
    
    
    a1=a1**(-1)
   !   call canonise( A1,A_cs,disp,a_l,a_nl,R,PHASE_ADVANCE )
   ! a1=a_cs
 !@1 computation of tune and damping   
   
    tune=atan2(imag(lam),real(lam))
    if(mat(1,2)<0) tune=-tune+twopi
    tune=tune/twopi
    dam=LOG(DET)/2.d0
    
    end subroutine diag_mat 
    
    
   subroutine invmatc( mat0,mati )
    implicit none
    complex(dp) mat(2,2),mati(2,2),det,mat0(2,2) 
   
    mat=mat0
    
    det=(mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1)) 
    mati(1,1)=mat(2,2)/det
    mati(2,2)=mat(1,1)/det
    mati(1,2)=-mat(1,2)/det
    mati(2,1)=-mat(2,1)/det
    
   end subroutine invmatc
   
   subroutine create_phasors(use_J) 
    implicit none
    logical, optional :: use_J
    real(dp) x
    
    x=1.0_dp
    if(present(use_J)) then
     if(use_J) x=sqrt(2.0_dp)
    endif
     !@1 &nbsp; complex transformation to phasors basis
     
    to_phasor=0
    to_phasor%v(1)%a(1)=1.0_dp/x
    to_phasor%v(1)%a(2)=i_/x
    to_phasor%v(2)%a(1)=1.0_dp/x
    to_phasor%v(2)%a(2)=-i_/x
    
    from_phasor=to_phasor**(-1)
    c_phasor=from_phasor
    ci_phasor=to_phasor

   end subroutine create_phasors 

!!! Some useful routines

   SUBROUTINE AVERAGE(F,A_CS,F_FLOQUET) 
    IMPLICIT NONE
    TYPE(MY_MAP) A_CS
    TYPE(MY_TAYLOR) F,F_FLOQUET
    INTEGER I
         
! 1) CREATES THE CHANGE OF BASIS FROM A ROTATION TO A FULLY DIAGONAL COMPLEX MATRIX
    CALL CREATE_PHASORS()  

!2) APPLIES A TO THE FUNCTION AND THEN GOES INTO THE COMPLEX PHASORS BASIS    
    F_FLOQUET=(F*A_CS)*C_PHASOR
    
    DO I=0,N_MONO
      IF(JORDER(I)>MY_ORDER) CYCLE
! 3) REJECTS TERMS OF UNEQUAL POWERS
      IF(JEXP1(I)/=JEXP2(I)) THEN 
       F_FLOQUET%A(I)=0.0_DP
      ENDIF       
    ENDDO
    
   END SUBROUTINE AVERAGE 


end module my_analysis