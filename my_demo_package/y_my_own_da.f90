
module my_own_da   
use precision_constants  
use file_handler 
implicit none
!@3 <font color="#FF0000" size="5"><p>Private routines only accessible through an interface</font></p>
 private input_my_taylor_in_real
 private input_real_in_my_taylor
 private input_comp_in_my_taylor
 private input_my_taylor_in_comp
 private input_my_taylor_in_my_taylor
 private mul   
 private dmulsc 
 private dscmul 
 private dmulscc
 private dscmulc
 private div     
 private ddivsc  
 private ddivscc 
 private dscdiv  
 private Idivsc  
 private add      
 private unaryADD 
 private daddsc   
 private dscadd   
 private daddscc  
 private dscaddc  
 private subs     
 private unarySUB 
 private dsubsc   
 private dscsub   
 private dsubscc  
 private dscsubc  
 private POW,powr
 private DEXPT  
 private DLOGT  
 private DSQRTT 
 private DCOST
 private DSINT
 private atan2t,dtant,atant
 private tpsa_expt
 private print_my_taylor
 private print_my_real
 private MONOT
 private SUBT,SUBTindex
 private alloc_my_taylor  
 private morpht,part,cutt
 private clean_taylor
 !@  <hr>
 integer :: my_order = 0    !@1 &nbsp; Order of the TPSA calculations 冪級数（べききゅうすう)
  complex   (dp), parameter :: I_=(0,1) !@1 &nbsp; complex number <font face="Times New Roman">&#8730;-1</font>
  real(dp) ::  epsclean=1.d-10   !@1 &nbsp; number used to insure that small numbers are set to zero
 integer :: n_tpsa_exp = 10      !@1 &nbsp; used in a "bad" way to do the exponential of a Taylor series
  real(dp) ::  epsprint=1.d-10
integer,parameter :: max_order =4
integer,parameter :: n_mono =34
integer, allocatable :: jexp1(:)
integer, allocatable :: jexp2(:)
integer, allocatable :: jexp3(:)
integer, allocatable :: jorder(:)
integer, allocatable :: mul_table(:,:)
integer, allocatable :: sub_index(:,:,:)
integer, allocatable :: der_table(:,:)
logical :: first_mul=.true.

! Definitions                
 TYPE my_taylor             
    complex(dp) a(0:n_mono)      !@2  &nbsp; Taylor series テイラー展開（テイラーてんかい) 浮動小数点数（ふどうしょうすうてんすう）
	!@3  Taylor = a<sub>0</sub>+<font color="#0000FF">a<sub>1</sub>x<sub>1</sub></font>+<font color="#0000FF">a
	!@3 <sub>2</sub>x<sub>2</sub></font>+<font color="#FF0000">a<sub>3</sub>x<sub>1</sub><sup>2</sup></font>+
	!@3 <font color="#FF0000">a<sub>4</sub>x<sub>2</sub><sup>2</sup></font>+<font color="#FF0000">a<sub>5</sub>x<sub>1</sub>
	!@3 x<sub>2</sub></font>+<font color="#CC00CC">a<sub>6</sub>x<sub>1</sub><sup>3</sup></font>
	!@3 +<font color="#CC00CC">a<sub>7</sub>x<sub>2</sub><sup>3</sup></font>+<font color="#CC00CC">a<sub>8</sub>x<sub>1</sub><sup>2</sup>x</font><sub><font color="#CC00CC">2</font></sub>+<font color="#CC00CC">a<sub>9</sub>x<sub>1</sub>x<sub>2</sub><sup>2</sup></font>

END TYPE my_taylor
type(my_taylor), allocatable :: table(:,:,:)
type(my_taylor) dx_1,dx_2,dx_3

 TYPE order 
  integer :: i=0
 end TYPE order 

 type(order) order_of_taylor

  INTERFACE assignment (=)
     MODULE PROCEDURE input_my_taylor_in_my_taylor   !@1 &nbsp; Taylor=Taylor
     MODULE PROCEDURE input_my_taylor_in_real   !@1 &nbsp; real=Taylor
     MODULE PROCEDURE input_real_in_my_taylor  !@1 &nbsp; Taylor=real 
     MODULE PROCEDURE input_real_in_my_taylors  !@1 &nbsp; Taylor=real 
     MODULE PROCEDURE input_comp_in_my_taylor  !@1 &nbsp; complex=Taylor
     MODULE PROCEDURE input_my_taylor_in_comp   !@1 &nbsp; complex=Taylor
     MODULE PROCEDURE input_order
  end  INTERFACE     
!@ <font color="#800000"><b>N.B. </b></font>In this toy DA, Taylor=Taylor is not 
!@ needed because Fortran90 does automatically the correct job. This is <b>
!@ <font color="#FF0000">not</font></b> true when we overload a pointer as in the 
!@ case of Berz's package.

  INTERFACE OPERATOR (*)  !@2  &nbsp; Operator overloading (or ad hoc polymorphism)  オペレータ(演算子)のオーバーロード
!@3  <b><font color="#FF0000">Taylor=Taylor*Taylor is interesting </font></b> click below </br>
     MODULE PROCEDURE mul   !@1 &nbsp; Taylor * Taylor   
     MODULE PROCEDURE dmulsc    !@1 &nbsp; Taylor * Real(dp)
     MODULE PROCEDURE dscmul     !@1 &nbsp;  Real(dp) * Taylor
     MODULE PROCEDURE dmulscc    !@1 &nbsp; Taylor * complex(dp)
     MODULE PROCEDURE dscmulc     !@1 &nbsp;  complex(dp) * Taylor
  END INTERFACE
  
  INTERFACE OPERATOR (/)  !@2  &nbsp; Operator overloading (or ad hoc polymorphism)  オペレータ(演算子)のオーバーロード
!@3  <b><font color="#FF0000">Taylor=Taylor/Taylor is also interesting because first case separating TPSA from DA part to create a nilpotent operator!</font></b> click below </br>
     MODULE PROCEDURE div    !@1 &nbsp; Taylor / Taylor
     MODULE PROCEDURE ddivsc  !@1 &nbsp; Taylor / Real(dp)
     MODULE PROCEDURE ddivscc  !@1 &nbsp; Taylor / complex(dp)
     MODULE PROCEDURE dscdiv  !@1 &nbsp; Real(dp) / Taylor
     MODULE PROCEDURE Idivsc  !@1 &nbsp; Taylor / Integer  <font color="#FF0000">&nbsp; &#8594; &nbsp; added because useful in example code</font>
  END INTERFACE

  INTERFACE OPERATOR (+)  !@2  &nbsp; Operator overloading (or ad hoc polymorphism)  オペレータ(演算子)のオーバーロード
     MODULE PROCEDURE add       !@1 &nbsp; Taylor + Taylor
     MODULE PROCEDURE unaryADD  !@1 &nbsp;  +Taylor
     MODULE PROCEDURE daddsc    !@1 &nbsp;  Taylor + Real(dp)
     MODULE PROCEDURE dscadd    !@1 &nbsp;  Real(dp) + Taylor
     MODULE PROCEDURE daddscc    !@1 &nbsp;  Taylor + complex(dp)
     MODULE PROCEDURE dscaddc    !@1 &nbsp;  complex(dp) + Taylor
  END INTERFACE

  INTERFACE OPERATOR (-)  !@2  &nbsp; Operator overloading (or ad hoc polymorphism)  オペレータ(演算子)のオーバーロード
     MODULE PROCEDURE subs      !@1 &nbsp; Taylor - Taylor
     MODULE PROCEDURE unarySUB     !@1 &nbsp;  -Taylor
     MODULE PROCEDURE dsubsc   !@1 &nbsp;  Taylor - Real(dp)
     MODULE PROCEDURE dscsub    !@1 &nbsp;  Real(dp) - Taylor
     MODULE PROCEDURE dsubscc   !@1 &nbsp;  Taylor - complex(dp)
     MODULE PROCEDURE dscsubc    !@1 &nbsp;  complex(dp) - Taylor
  END INTERFACE


  INTERFACE OPERATOR (**)  !@2  &nbsp; Operator overloading (or ad hoc polymorphism)  オペレータ(演算子)のオーバーロード
     MODULE PROCEDURE POW    !@1 &nbsp; Taylor ** Integer 
     MODULE PROCEDURE POWr    !@1 &nbsp; Taylor ** real(dp) 
  END INTERFACE


!@3  <p><i><font size="5">&nbsp;Overloading standard procedures </font></i></p>

  INTERFACE exp
     MODULE PROCEDURE DEXPT   !@1 &nbsp; exp(Taylor)
  END INTERFACE

  INTERFACE tpsa_exp
     MODULE PROCEDURE tpsa_expt   !@1 &nbsp; exp(Taylor) using a sum up to n_tpsa_exp
  END INTERFACE

  INTERFACE LOG
     MODULE PROCEDURE DLOGT   !@1 &nbsp; log(Taylor)
  END INTERFACE

  INTERFACE SQRT
     MODULE PROCEDURE DSQRTT    !@1 &nbsp; sqrt(Taylor)
  END INTERFACE

  INTERFACE COS
     MODULE PROCEDURE DCOST    !@1 &nbsp; cos(Taylor)
  END INTERFACE

  INTERFACE SIN
     MODULE PROCEDURE DSINT   !@1 &nbsp; sin(Taylor)
  END INTERFACE

  INTERFACE tan
     MODULE PROCEDURE dtant   !@1 &nbsp; tan(Taylor)
  END INTERFACE

  INTERFACE atan
     MODULE PROCEDURE atant   !@1 &nbsp; tan(Taylor)
  END INTERFACE

  INTERFACE atan2
     MODULE PROCEDURE atan2t   !@1 &nbsp; Atan2(Taylor,Taylor)
  END INTERFACE

 INTERFACE print_for_human
     MODULE PROCEDURE print_my_taylor_human
  END INTERFACE

  INTERFACE print
     MODULE PROCEDURE print_my_taylor   !@1 &nbsp; Fancy print for human eye
     MODULE PROCEDURE print_my_real
  END INTERFACE

!@3  <p><i><font size="5">User defined operator</font></i></p>

  INTERFACE OPERATOR (.mon.)
     MODULE PROCEDURE MONOT    !@1 <font color="#FF0000">Creates the monomial </font><font color="#0000FF">r x<sub>i</sub></font>
  END INTERFACE

  INTERFACE OPERATOR (.SUB.)
     MODULE PROCEDURE SUBT    !@1 <font color="#FF0000">Retrieves "r" from the monomial </font><font color="#0000FF">(r+i c) x<sub>i</sub></font>
  END INTERFACE

  INTERFACE OPERATOR (.par.)
     MODULE PROCEDURE part    
  END INTERFACE


  INTERFACE OPERATOR (.cut.)
     MODULE PROCEDURE cutt    
  END INTERFACE

  INTERFACE OPERATOR (.index.)
     MODULE PROCEDURE SUBTindex    
  END INTERFACE

! Destructors and Constructors for My_taylor

  INTERFACE alloc
     MODULE PROCEDURE alloc_my_taylor   !@1 &nbsp; for compatibility with FPP
  END INTERFACE

  INTERFACE kill
     MODULE PROCEDURE alloc_my_taylor   !@1 &nbsp; for compatibility with FPP
  END INTERFACE

  INTERFACE morph
     MODULE PROCEDURE morpht   !@1 &nbsp; for compatibility with FPP
  END INTERFACE

  INTERFACE clean
     MODULE PROCEDURE clean_taylor   !@1 &nbsp;  removes numbers smaller than epsclean
  END INTERFACE

contains
   subroutine input_order(o,i)
    implicit none
     integer, intent(in) :: i
     type(order), intent(out) :: o
     integer j,nmono

     if(i>max_order) then
      write(6,*) " Increase max_order"
      write(6,*) " The analysis package will not work beyond order=4"
       nmono=1
       do j=i+3,i+1,-1
        nmono=nmono*j
       enddo
       nmono=nmono/6-1
       write(6,*) "n_mono should be ",nmono
      stop
     endif
     o%i=i
     my_order=i

     dx_1=1.0_dp.mon.1
     dx_2=1.0_dp.mon.2
     dx_3=1.0_dp.mon.3

    
   end subroutine input_order

  subroutine clean_taylor( S2 )
    implicit none
    TYPE (my_taylor), INTENT (INout) :: S2
    integer i
    
    do i=0,n_mono
     if(abs (real(S2%a(i),kind=dp) ) <epsclean) S2%a(i)= i_*aimag(S2%a(i))
     if(abs(aimag(S2%a(i)))  <epsclean)         S2%a(i)=real(S2%a(i),kind=dp)
    enddo     
  END subroutine clean_taylor

  subroutine clean_taylor_PRINT( S2 )
    implicit none
    TYPE (my_taylor), INTENT (INout) :: S2
    integer i
    
    do i=0,n_mono
     if(abs (real(S2%a(i),kind=dp) ) <EPSPRINT) S2%a(i)= i_*aimag(S2%a(i))
     if(abs(aimag(S2%a(i)))  <EPSPRINT)         S2%a(i)=real(S2%a(i),kind=dp)
    enddo     
  END subroutine clean_taylor_PRINT


  subroutine alloc_my_taylor( S2 )
    implicit none
    TYPE (my_taylor), INTENT (INout) :: S2

     S2%a=0.0_dp
     
  END subroutine alloc_my_taylor



  FUNCTION morpht( t )
    implicit none
    TYPE (my_taylor) morpht
    TYPE (my_taylor), INTENT(IN) :: t
    
	morpht= t
  
  END FUNCTION morpht

  FUNCTION MONOT( R, I )
    implicit none
    TYPE (my_taylor) MONOT
    REAL(DP), INTENT(IN) :: R
    INTEGER, INTENT(IN) :: I
    integer j(3)
     if(first_mul) then
      call set_j_arrays 
      call multiplication_table
      first_mul=.false.
     endif

    j=0
    j(i)=1
    
     MONOT=0.0_DP
     MONOT%a( sub_index(j(1),j(2),j(3)) )= R

  END FUNCTION MONOT

  FUNCTION SUBTindex( T, I )
    implicit none
    real(dp) SUBTindex
    TYPE (my_taylor), INTENT(IN) :: T
    INTEGER, INTENT(IN) :: I
    INTEGER K,j(3)
	
	     SUBTindex=0.0_DP
    
     j=0
     j(i)=1
      
    SUBTindex=t%a(sub_index(j(1),j(2),j(3)))

  END FUNCTION SUBTindex

  FUNCTION SUBT( T, I )
    implicit none
    real(dp) SUBT
    TYPE (my_taylor), INTENT(IN) :: T
    INTEGER, INTENT(IN) :: I(:)
    INTEGER K,j(3)
	
	     SUBT=0.0_DP
    
    j=0
    do k=1,size(i)
     j(k)=i(k)
    enddo  
      
    SUBT=t%a(sub_index(j(1),j(2),j(3)))

    
  END FUNCTION SUBT

  FUNCTION part( T, I )
    implicit none
    TYPE (my_taylor) part,te
    TYPE (my_taylor), INTENT(IN) :: T
    INTEGER, INTENT(IN) :: I(2)
    INTEGER K 
	

         te=0.0_dp
       do k=0,n_mono
       if(jexp1(k)==i(1).and.jexp2(k)==i(2)) then
         te=t%a(sub_index(i(1),i(2),jexp3(k)))*dx_3**jexp3(k)+te
       endif
      enddo
      
      part=te

  END FUNCTION part

  FUNCTION cutt( T, I )
    implicit none
    TYPE (my_taylor) cutt,te
    TYPE (my_taylor), INTENT(IN) :: T
    INTEGER, INTENT(IN) :: I 
    INTEGER K 
	

         te=0.0_dp
       do k=0,n_mono
       if(jexp1(k)+jexp2(k)+jexp3(k)<i) then
         te%a(k) = t%a(k) 
       endif
      enddo
      
      cutt=te

  END FUNCTION cutt

  FUNCTION add( S1, S2 )
    implicit none
    TYPE (my_taylor) add
    TYPE (my_taylor), INTENT (IN) :: S1, S2

     add%a=S1%a + S2%a     

  END FUNCTION add

  FUNCTION daddsc( S1, sc )
    implicit none
    TYPE (my_taylor) daddsc
    TYPE (my_taylor), INTENT (IN) :: S1 
    real(dp), INTENT (IN) :: sc
     
     daddsc=s1
     daddsc%a(0)= s1%a(0) + sc    
     call clean(daddsc)
  END FUNCTION daddsc

  FUNCTION dscadd( sc ,  S1)
    implicit none
    TYPE (my_taylor) dscadd
    TYPE (my_taylor), INTENT (IN) :: S1 
    real(dp), INTENT (IN) :: sc
     
     dscadd=s1
     dscadd%a(0)= s1%a(0) + sc    
     call clean(dscadd)

  END FUNCTION dscadd

  FUNCTION daddscc( S1, sc )
    implicit none
    TYPE (my_taylor) daddscc
    TYPE (my_taylor), INTENT (IN) :: S1 
    complex(dp), INTENT (IN) :: sc
     
     daddscc=s1
     daddscc%a(0)= s1%a(0) + sc    
     call clean(daddscc)

  END FUNCTION daddscc

  FUNCTION dscaddc( sc ,  S1)
    implicit none
    TYPE (my_taylor) dscaddc
    TYPE (my_taylor), INTENT (IN) :: S1 
    complex(dp), INTENT (IN) :: sc
     
     dscaddc=s1
     dscaddc%a(0)= s1%a(0) + sc    
     call clean(dscaddc)
  END FUNCTION dscaddc

  FUNCTION unaryadd( S1 )
    implicit none
    TYPE (my_taylor) unaryadd
    TYPE (my_taylor), INTENT (IN) :: S1 
     
     unaryadd=s1

  END FUNCTION unaryadd


  FUNCTION subs( S1, S2 )
    implicit none
    TYPE (my_taylor) subs
    TYPE (my_taylor), INTENT (IN) :: S1, S2

     subs%a=S1%a - S2%a     
     call clean(subs)

  END FUNCTION subs

  FUNCTION unarySUB( S1 )
    implicit none
    TYPE (my_taylor) unarySUB
    TYPE (my_taylor), INTENT (IN) :: S1 
     
     unarySUB%a=-s1%a

  END FUNCTION unarySUB

  FUNCTION dsubsc( S1, sc )
    implicit none
    TYPE (my_taylor) dsubsc
    TYPE (my_taylor), INTENT (IN) :: S1 
    real(dp), INTENT (IN) :: sc
     
     dsubsc=s1
     dsubsc%a(0)= s1%a(0) - sc    
     call clean(dsubsc)

  END FUNCTION dsubsc


  FUNCTION dscsub( sc , S1 )
    implicit none
    TYPE (my_taylor) dscsub
    TYPE (my_taylor), INTENT (IN) :: S1 
    real(dp), INTENT (IN) :: sc
     
     dscsub=-s1   ! uses unary sub
     dscsub=sc + dscsub   ! uses add   
     call clean(dscsub)

  END FUNCTION dscsub

  FUNCTION dsubscc( S1, sc )
    implicit none
    TYPE (my_taylor) dsubscc
    TYPE (my_taylor), INTENT (IN) :: S1 
    complex(dp), INTENT (IN) :: sc
     
     dsubscc=s1
     dsubscc%a(0)= s1%a(0) - sc    
     call clean(dsubscc)
  END FUNCTION dsubscc


  FUNCTION dscsubc( sc , S1 )
    implicit none
    TYPE (my_taylor) dscsubc
    TYPE (my_taylor), INTENT (IN) :: S1 
    complex(dp), INTENT (IN) :: sc
     
     dscsubc=-s1   ! uses unary sub
     dscsubc=sc + dscsubc   ! uses add   
     call clean(dscsubc)

  END FUNCTION dscsubc


  FUNCTION mul( S1, S2 )
    implicit none
    TYPE (my_taylor) mul
    TYPE (my_taylor), INTENT (IN) :: S1, S2
    integer i,j,k
!@3  <hr color="#00FFFF" align="left" size="4">
!@3  T<sub>1</sub> =a<sub>0</sub>+<font color="#0000FF">a<sub>1</sub>x<sub>1</sub></font>+<font color="#0000FF">a<sub>2</sub>x<sub>2</sub></font>+<font color="#FF0000">a<sub>3</sub>x<sub>1</sub><sup>2</sup></font>+
!@3  	<font color="#FF0000">a<sub>4</sub>x<sub>2</sub><sup>2</sup></font>+<font color="#FF0000">a<sub>5</sub>x<sub>1</sub>
!@3  	x<sub>2</sub></font>+<font color="#CC00CC">a<sub>6</sub>x<sub>1</sub><sup>3</sup></font>
!@3  	+<font color="#CC00CC">a<sub>7</sub>x<sub>2</sub><sup>3</sup></font>+<font color="#CC00CC">a<sub>8</sub>x<sub>1</sub><sup>2</sup>x</font><sub><font color="#CC00CC">2</font></sub>+<font color="#CC00CC">a<sub>9</sub>x<sub>1</sub>x<sub>2</sub><sup>2</sup></font>&nbsp;&nbsp; 
!@3       = (a<sub>0</sub>,a<sub>1</sub>,a<sub>2</sub>,a<sub>3</sub>,a<sub>4</sub>,a<sub>5</sub>,a<sub>6</sub>,a<sub>7</sub>,a<sub>8</sub>,a<sub>9</sub>)</br>
!@3  T<sub>2</sub> =b<sub>0</sub>+<font color="#0000ff">b<sub>1</sub>x<sub>1</sub></font>+<font color="#0000ff">b<sub>2</sub>x<sub>2</sub></font>+<font color="#ff0000">b<sub>3</sub>x<sub>1</sub><sup>2</sup></font>+
!@3   <font color="#ff0000">b<sub>4</sub>x<sub>2</sub><sup>2</sup></font>+<font color="#ff0000">b<sub>5</sub>x<sub>1</sub> 
!@3   x<sub>2</sub></font>+<font color="#cc00cc">b<sub>6</sub>x<sub>1</sub><sup>3</sup></font> 
!@3   +<font color="#cc00cc">b<sub>7</sub>x<sub>2</sub><sup>3</sup></font>+<font color="#cc00cc">b<sub>8</sub>x<sub>1</sub><sup>2</sup>x</font><sub><font color="#cc00cc">2
!@3   </font></sub>+<font color="#cc00cc">b<sub>9</sub>x<sub>1</sub>x<sub>2</sub><sup>2</sup></font> 
!@3   = (b<sub>0</sub>,b<sub>1</sub>,b<sub>2</sub>,b<sub>3</sub>,b<sub>4</sub>,b<sub>5</sub>,b<sub>6</sub>,b<sub>7</sub>,b<sub>8</sub>,b<sub>9</sub>)</br>
!@3   </br> 
!@3  T<sub>1</sub> * T<sub>2</sub> =a<sub>0</sub>b<sub>0</sub>+{a<sub>0</sub><font color="#0000ff">b<sub>1</sub></font>+b<sub>0</sub>
!@3  <font color="#0000FF">a<sub>1</sub>}x<sub>1</sub></font>+{a<sub>0</sub><font color="#0000ff">b<sub>2</sub></font>+b<sub>0</sub>
!@3  <font color="#0000FF">a<sub>2</sub>}x<sub>2</sub></font>+{a<sub>0</sub><font color="#0000ff">b</font><sub><font color="#0000FF">3</font></sub>+b<sub>0</sub>
!@3  <font color="#0000FF">a<sub>3</sub></font>+<font color="#0000FF">a<sub>1</sub><sup>2</sup><sub> </sub>x<sub>1</sub><sup>2</sup>}x<sub>1</sub><sup>2</sup>
!@3  </font>+{a<sub>0</sub><font color="#0000ff">b</font><sub><font color="#0000FF">4</font></sub>+b<sub>0</sub><font color="#0000FF">a</font><sub><font color="#0000FF">4</font>
!@3  </sub>+<font color="#0000FF">a<sub>2</sub><sup>2</sup>}x<sub>2</sub><sup>2</sup></font>+{a<sub>0</sub><font color="#0000ff">b</font><sub><font color="#0000FF">5</font></sub>+<sub>
!@3  </sub>b<sub>0</sub><font color="#0000FF">a</font><sub><font color="#0000FF">5</font></sub>+a<sub>1</sub><font color="#0000ff">b</font><sub><font color="#0000FF">2</font></sub>+<sub>
!@3  </sub>b<sub>1</sub><font color="#0000FF">a<sub>2</sub>}x<sub>1</sub>x<sub>2</sub>+ 
!@3       ...</font></br> 
!@3  <font color="#0000FF">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; = </font><font size="4">(</font>a<sub>0</sub>b<sub>0</sub><font size="4">
!@3   <b>,</b> </font>a<sub>0</sub><font color="#0000ff">b<sub>1</sub></font>+b<sub>0</sub><font color="#0000FF">a<sub>1<font size="4">
!@3   </font></sub></font><font size="4"><b>,</b> </font>a<sub>0</sub><font color="#0000ff">b<sub>2</sub></font>+b<sub>0</sub><font color="#0000FF">a<sub>2 </sub></font>
!@3   <b><font size="4">,</font></b> a<sub>0</sub><font color="#0000ff">b</font><sub><font color="#0000FF">3</font></sub>+b<sub>0</sub><font color="#0000FF">a<sub>3</sub></font>+<font color="#0000FF">a<sub>1</sub><sup>2</sup><sub> </sub>x<sub>1</sub><sup>2 </sup></font>
!@3   <font size="4"><b>,</b></font> a<sub>0</sub><font color="#0000ff">b</font><sub><font color="#0000FF">4</font></sub>+b<sub>0</sub><font color="#0000FF">a</font><sub><font color="#0000FF">4</font></sub>+<font color="#0000FF">a<sub>2</sub><sup>2 </sup></font>
!@3   <font size="4"><b>,</b> </font>a<sub>0</sub><font color="#0000ff">b</font><sub><font color="#0000FF">5</font></sub>+b<sub>0</sub><font color="#0000FF">a</font><sub><font color="#0000FF">5</font></sub>+a<sub>1</sub><font color="#0000ff">b</font><sub><font color="#0000FF">2</font></sub>+b<sub>1</sub><font color="#0000FF">a<sub>2<font size="4">
!@3   </font></sub> </font><font size="4"><b>,...</b>)</font> </br>
!@3   <a href="#*">Back up</a>
!@3  <hr color="#00FFFF" align="left" size="4">	 
     if(first_mul) then
      call set_j_arrays 
      call multiplication_table
      first_mul=.false.
     endif
!@3 <p>Old Code done by hand to 3<sup>rd </sup>order</p>
!@3 <dl>
!@3<dt><br>
!@3<b>  mul%a(0)=S1%a(0)*S2%a(0)</b> </br> 
!@3 <br><b>  mul%a(1)=S1%a(1)*S2%a(0)+S1%a(0)*S2%a(1)</b></br>
!@3 <br><b>  mul%a(2)=S1%a(2)*S2%a(0)+S1%a(0)*S2%a(2)</b></br>
!@3 <br><b>  mul%a(3)=S1%a(3)*S2%a(0)+S1%a(0)*S2%a(3)+S1%a(1)*S2%a(1)</b></br>
!@3 <br><b>  mul%a(4)=S1%a(4)*S2%a(0)+S1%a(0)*S2%a(4)+S1%a(2)*S2%a(2) </b></br>
!@3 <br><b>  mul%a(5)=S1%a(5)*S2%a(0)+S1%a(0)*S2%a(5)+S1%a(1)*S2%a(2)+S1%a(2)*S2%a(1) </b></br>
!@3 <br><b>  mul%a(6)=S1%a(6)*S2%a(0)+S1%a(0)*S2%a(6)+S1%a(1)*S2%a(3)+S1%a(3)*S2%a(1) </b></br>
!@3 <br><b>  mul%a(7)=S1%a(7)*S2%a(0)+S1%a(0)*S2%a(7)+S1%a(2)*S2%a(4)+S1%a(4)*S2%a(2) </b></br>
!@3 <br><b>  mul%a(8)=S1%a(8)*S2%a(0)+S1%a(0)*S2%a(8)+S1%a(1)*S2%a(5)+S1%a(5)*S2%a(1)+S1%a(2)*S2%a(3)+S1%a(3)*S2%a(2) </b></br>
!@3 <br><b>  mul%a(9)=S1%a(9)*S2%a(0)+S1%a(0)*S2%a(9)+S1%a(1)*S2%a(4)+S1%a(4)*S2%a(1)+S1%a(2)*S2%a(5)+S1%a(5)*S2%a(2)</b></br>
!@3 </dl></dt>


!@3 <p><font color="#FF0000"><b>New Code using tables to 4<sup>th </sup>order</b></font></p>
!@3 <p><font color="#FF0000"><b> This is necessary for speed</b></font></p>

     
     mul%a=0.0_dp
     do i=0,n_mono
     do j=0,n_mono
      k=mul_table(i,j)
      if(k==-1.or.jorder(k)>my_order) cycle
      
      mul%a(k)=S1%a(i)*S2%a(j)+mul%a(k) 
     
     enddo
     enddo
     
         
     call clean(mul)

  END FUNCTION mul

  subroutine set_j_arrays()
  implicit none
  integer i,j,k,l,m_n
  
  m_n=0
  do i=0,max_order
  
   do l=0,i
    do k=0,i
     do j=0,i
      if(l+k+j/=i) cycle
      m_n=1+m_n
     enddo
    enddo
   enddo
  
  enddo

  allocate(jexp1(0:m_n-1))
  allocate(jexp2(0:m_n-1))
  allocate(jexp3(0:m_n-1))
  allocate(jorder(-1:m_n-1))
  allocate(mul_table(0:m_n-1,0:m_n-1))
  allocate(der_table(0:m_n-1,3))
  allocate(sub_index(0:max_order,0:max_order,0:max_order))
  allocate(table(0:max_order,0:max_order,0:max_order))
  
  sub_index=-1
  
  jorder=-1
  mul_table=-1
  der_table=-1
  m_n=0
  do i=0,max_order
  
   do l=0,i
    do k=0,i
     do j=0,i
      if(l+k+j/=i) cycle
      jexp1(m_n)=j
      jexp2(m_n)=k
      jexp3(m_n)=l
      jorder(m_n)=i
      m_n=1+m_n
     enddo
    enddo
   enddo
  
  enddo
  
  

  end subroutine set_j_arrays

  subroutine multiplication_table()
  implicit none
  integer i,j,k,n,i1,i2,i3
   do i=0,n_mono
    do j=0,n_mono
     i1=  jexp1(i)+jexp1(j) 
     i2=  jexp2(i)+ jexp2(j)
     i3=  jexp3(i)+ jexp3(j)
     n=i1+i2+i3
      if(n>max_order) cycle
     do k=0,n_mono
  
      if(i1==jexp1(k).and.i2==jexp2(k).and.i3==jexp3(k)) then
       mul_table(i,j)=k      
      exit
    endif
 
    enddo 
   enddo
   enddo
      do i=1,n_mono
       i1=  jexp1(i)-1
       i2=  jexp2(i)
       i3=  jexp3(i)
            do k=0,n_mono
  
              if(i1==jexp1(k).and.i2==jexp2(k).and.i3==jexp3(k)) then
                der_table(i,1)=k
              endif
            
           enddo
      enddo
           
      do i=1,n_mono
       i1=  jexp1(i)
       i2=  jexp2(i)-1
       i3=  jexp3(i)
            do k=0,n_mono
  
              if(i1==jexp1(k).and.i2==jexp2(k).and.i3==jexp3(k)) then
                der_table(i,2)=k
              endif
            
           enddo
      enddo
 
      do i=1,n_mono
       i1=  jexp1(i)
       i2=  jexp2(i)
       i3=  jexp3(i)-1
            do k=0,14
  
              if(i1==jexp1(k).and.i2==jexp2(k).and.i3==jexp3(k)) then
                der_table(i,3)=k
              endif
            
           enddo
      enddo
      
       do k=0,n_mono
       if(jexp1(k)+jexp2(k)+jexp3(k)<=max_order) then
         sub_index(jexp1(k),jexp2(k),jexp3(k))=k
       endif
      enddo

  end subroutine multiplication_table

  FUNCTION dmulsc( S1, sc )
    implicit none
    TYPE (my_taylor) dmulsc
    TYPE (my_taylor), INTENT (IN) :: S1 
    real(dp), INTENT (IN) :: sc

     dmulsc%a= s1%a*sc    
     call clean(dmulsc)

  END FUNCTION dmulsc

  FUNCTION dscmul( sc ,S1 )
    implicit none
    TYPE (my_taylor) dscmul
    TYPE (my_taylor), INTENT (IN) :: S1 
    real(dp), INTENT (IN) :: sc

     dscmul%a= s1%a*sc    
     call clean(dscmul)

  END FUNCTION dscmul

  FUNCTION dmulscc( S1, sc )
    implicit none
    TYPE (my_taylor) dmulscc
    TYPE (my_taylor), INTENT (IN) :: S1 
    complex(dp), INTENT (IN) :: sc

     dmulscc%a= s1%a*sc    
     call clean(dmulscc)

  END FUNCTION dmulscc

  FUNCTION dscmulc( sc ,S1 )
    implicit none
    TYPE (my_taylor) dscmulc
    TYPE (my_taylor), INTENT (IN) :: S1 
    complex(dp), INTENT (IN) :: sc

     dscmulc%a= s1%a*sc    
     call clean(dscmulc)
  END FUNCTION dscmulc

  FUNCTION ddivsc( S1, sc )
    implicit none
    TYPE (my_taylor) ddivsc
    TYPE (my_taylor), INTENT (IN) :: S1 
    real(dp), INTENT (IN) :: sc
     
     ddivsc%a= s1%a/sc    
     call clean(ddivsc)

  END FUNCTION ddivsc

  FUNCTION ddivscc( S1, sc )
    implicit none
    TYPE (my_taylor) ddivscc
    TYPE (my_taylor), INTENT (IN) :: S1 
    complex(dp), INTENT (IN) :: sc
     
     ddivscc%a= s1%a/sc    
     call clean(ddivscc)
  END FUNCTION ddivscc



  FUNCTION Idivsc( S1, sc )
    implicit none
    TYPE (my_taylor) Idivsc
    TYPE (my_taylor), INTENT (IN) :: S1 
    INTEGER, INTENT (IN) :: sc
     
     Idivsc%a= s1%a/sc    
     call clean(Idivsc)
  END FUNCTION Idivsc

  FUNCTION POW( S1,N )
    implicit none
    TYPE (my_taylor) POW , T
    TYPE (my_taylor), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: N
INTEGER I

    POW=1.D0

    IF(N<0) THEN
     T=1.D0/S1
     DO I=1,-N
 POW=POW*T
 ENDDO
ELSE
     DO I=1,N
 POW=POW*S1
 ENDDO

ENDIF
     call clean(POW)

    END FUNCTION POW

  FUNCTION POWr( S1,N )
    implicit none
    TYPE (my_taylor) POWr  
    TYPE (my_taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: N
 

       powr=log(s1)*n
       powr=exp(powr)
     call clean(POWr)

    END FUNCTION POWr

  FUNCTION inv( S1 )
    implicit none
    TYPE (my_taylor) inv,t,tt
    TYPE (my_taylor), INTENT (IN) :: S1
    integer i
!@3 T<sub>1</sub><sup>-1</sup> = (a<sub>0</sub>+<font color="#0000FF">a<sub>1</sub>x<sub>1</sub></font>+<font color="#0000FF">a<sub>2</sub>x<sub>2</sub>
!@3  </font>+<font color="#FF0000">a<sub>3</sub>x<sub>1</sub><sup>2</sup></font>+
!@3  	<font color="#FF0000">a<sub>4</sub>x<sub>2</sub><sup>2</sup></font>+<font color="#FF0000">a<sub>5</sub>x<sub>1</sub>x<sub>2</sub></font>+
!@3  	<font color="#CC00CC">a<sub>6</sub>x<sub>1</sub><sup>3</sup></font>
!@3  	+<font color="#CC00CC">a<sub>7</sub>x<sub>2</sub><sup>3</sup></font>+<font color="#CC00CC">a<sub>8</sub>x<sub>1</sub><sup>2</sup>x</font><sub>
!@3  	<font color="#CC00CC">2</font></sub>+<font color="#CC00CC">a<sub>9</sub>x<sub>1</sub>x<sub>2</sub><sup>2</sup>)<sup>-1</sup></font></br>
!@3  = a<sub>0</sub><sup>-1</sup>(1 + a<sub>1/</sub>a<sub>0 </sub>x<sub>1</sub>+a<sub>2/</sub>a<sub>0 </sub>x<sub>2</sub>+a<sub>3/</sub>a<sub>0 
!@3  </sub>x<sub>1</sub><sup>2</sup>+ a<sub>4/</sub>a<sub>0 </sub>x<sub>2</sub><sup>2</sup>+a<sub>5/</sub>a<sub>0 </sub>x<sub>1</sub>x<sub>2</sub>
!@3  +a<sub>6/</sub>a<sub>0 </sub>x<sub>1</sub><sup>3</sup>
!@3  	+a<sub>7/</sub>a<sub>0 </sub>x<sub>2</sub><sup>3</sup>+a<sub>8/</sub>a<sub>0 </sub>x<sub>1</sub><sup>2</sup>x<sub>2</sub>+a<sub>9/</sub>a<sub>0 
!@3  	</sub>x<sub>1</sub>x<sub>2</sub><sup>2</sup>)<sup>-1</sup></br>
!@3  = a<sub>0</sub><sup>-1</sup> (1 + T)<sup>-1</sup> <font color="#FF0000"><b>where T is 
!@3  nilpotent, i.e., T<sup>4</sup>=0!</b></font> </br>
!@3  <font color="#FF0000"><b>
!@3   Therefore (1+T)<sup>-1 </sup>= 1-T+T<sup>2</sup>-T<sup>3</sup> exactly in our algebra!</b></font></br>
!@3   <a href="#/">Back up</a>
    
T=S1/S1%A(0)
T%A(0)=0.D0
 INV=1.0_DP
 TT=1.0_DP
do i=1,my_order
 TT=-T*TT
 INV=inv+TT
enddo

    INV=INV/S1%A(0)
     call clean(INV)

  END FUNCTION inv

  FUNCTION DIV( S1, S2 )
    implicit none
    TYPE (my_taylor) DIV
    TYPE (my_taylor), INTENT (IN) :: S1, S2
!@3<b><font color="#FF0000">Go to function </font></b><a href="#INV">INV</a></br>
   DIV=INV(S2)
   DIV=S1*DIV
     call clean(DIV)
  END FUNCTION DIV

  FUNCTION dscdiv(  sc , S1)
    implicit none
    TYPE (my_taylor) dscdiv
    TYPE (my_taylor), INTENT (IN) :: S1 
    real(dp), INTENT (IN) :: sc
     
     dscdiv=INV(S1)
     dscdiv= SC *  dscdiv
     call clean(dscdiv)

  END FUNCTION dscdiv

!  DEFININING my_taylor=CONSTANT   
  subroutine input_real_in_my_taylor( S2, S1 )
    implicit none
    real(dp), INTENT (IN) :: S1
    TYPE (my_taylor), INTENT (inout) :: S2

     S2%a=0.0_dp 
     S2%a(0)=s1

  END subroutine input_real_in_my_taylor

!  DEFININING my_taylor=CONSTANT   
  subroutine input_real_in_my_taylors( S2, S1 )
    implicit none
    real(dp), INTENT (IN) :: S1
    TYPE (my_taylor), INTENT (inout) :: S2(:)
    integer i
     do i=1,size(s2)
      S2(i)%a=0.0_dp 
      S2(i)%a(0)=s1
     enddo
  END subroutine input_real_in_my_taylors

  subroutine input_comp_in_my_taylor( S2, S1 )
    implicit none
    complex(dp), INTENT (IN) :: S1
    TYPE (my_taylor), INTENT (inout) :: S2

     S2%a=0.0_dp 
     S2%a(0)=s1

  END subroutine input_comp_in_my_taylor

  subroutine input_my_taylor_in_comp( S2, S1 )
    implicit none
    TYPE (my_taylor), INTENT (IN) :: S1
    complex(dp), INTENT (inout) :: S2

     S2=S1%a(0) 
     
  END subroutine input_my_taylor_in_comp
  
    subroutine input_my_taylor_in_my_taylor( S2, S1 )
    implicit none
    TYPE (my_taylor), INTENT (IN) :: S1
    TYPE (my_taylor), INTENT (inout) :: S2

     S2%a=S1%a 
     
  END subroutine input_my_taylor_in_my_taylor

    subroutine input_my_taylor_in_real( S2, S1 )
    implicit none
    TYPE (my_taylor), INTENT (IN) :: S1
    real(dp), INTENT (inout) :: S2

     S2=S1%a(0) 
     
  END subroutine input_my_taylor_in_real


!  DEFININING EXP

  FUNCTION DEXPT( S1 )
    implicit none
    TYPE (my_taylor) DEXPT,t,tt
    TYPE (my_taylor), INTENT (IN) :: S1
    integer i,no
T=S1
T%A(0)=0.0_DP
    

DEXPT=1.0_dp
tt=1.0_dp

	do i=1,MY_ORDER
	 tt=tt*t/i
     DEXPT=DEXPT + tt
    enddo

    DEXPT=EXP(S1%A(0))*DEXPT
     call clean(DEXPT)

  END FUNCTION DEXPT 


  FUNCTION DLOGT( S1 )
    implicit none
    TYPE (my_taylor) DLOGT,T,TT
    TYPE (my_taylor), INTENT (IN) :: S1
    INTEGER I
T=S1/S1%A(0); T%A(0)=0.0_DP;
      DLOGT=0.0_dp
     TT=-1.0_DP
     DO I=1,MY_ORDER 
      TT=-T*TT
      DLOGT=DLOGT+TT/i      
     ENDDO
     
    DLOGT=DLOGT+ LOG(S1%A(0))
     call clean(DLOGT)

  END FUNCTION  DLOGT 

  FUNCTION DSQRTT( S1 )
    implicit none
    TYPE (my_taylor) DSQRTT
    TYPE (my_taylor), INTENT (IN) :: S1
!@3  <p><i><sup><font size="4" color="#0066FF"> </font></sup></i></p>

    DSQRTT= log(s1)/2.0_dp
    DSQRTT=exp(DSQRTT)

     call clean(DSQRTT)

  END FUNCTION  DSQRTT 

  FUNCTION DCOST(S1)
    implicit none
    TYPE (my_taylor) DCOST
    TYPE (my_taylor), INTENT (IN) :: S1
!@3  <p><i><sup><font size="4" color="#0066FF"> </font></sup></i></p>

DCOST=(exp(i_*s1)+exp(-i_*s1))/2.0_dp

     call clean(DCOST)

  END FUNCTION  DCOST 

  FUNCTION atant(S1)
    implicit none
    TYPE (my_taylor) atant,c
    TYPE (my_taylor), INTENT (IN) :: S1
 
!@3  <p><i><sup><font size="4" color="#0066FF"> </font></sup></i></p>
     c=sqrt(1.0_dp/(1.0_dp+s1**2))
     c=c+i_*s1*c  ! cos(atant)+i * sin(atant)

     atant=log(c)/i_

     call clean(atant)

  END FUNCTION  atant 

  FUNCTION DtanT(S1)
    implicit none
    TYPE (my_taylor) DtanT
    TYPE (my_taylor), INTENT (IN) :: S1
    integer i
!@3  <p><i><sup><font size="4" color="#0066FF"> </font></sup></i></p>

    DtanT=sin(s1)/cos(s1)

     call clean(DtanT)

  END FUNCTION  DtanT 

  FUNCTION atan2t(S1,s2)
    implicit none
    TYPE (my_taylor) atan2t
    TYPE (my_taylor), INTENT (IN) :: S1,s2


     atan2t=atan(s1/s2)
     atan2t%a(0)=0.0_dp

     atan2t= atan2t + atan2(real(s1%a(0),kind=dp),real(s2%a(0),kind=dp))
     call clean(atan2t)

  END FUNCTION  atan2t 

  FUNCTION DSINT(S1)
    implicit none
    TYPE (my_taylor) DSINT
    TYPE (my_taylor), INTENT (IN) :: S1
    
!@3  <p><i><sup><font size="4" color="#0066FF">Lazy implementation: should be  generalized, only good to 4th order.</font></sup></i></p>
 
DSINT=(exp(i_*s1)-exp(-i_*s1))/2.0_dp/i_
 
     call clean(DSINT)

  END FUNCTION  DSINT 

  FUNCTION tpsa_expt( S1 )
    implicit none
    TYPE (my_taylor) tpsa_expt,t
    TYPE (my_taylor), INTENT (IN) :: S1
    integer i

   
   T=1.0_DP
   tpsa_expt=1.0_DP

   do i=1,n_tpsa_exp
    T=S1*T/I
    tpsa_expt=tpsa_expt+T
   enddo
      call clean(tpsa_expt)
   
  END FUNCTION tpsa_expt 



  subroutine print_my_real( s , mfi)
    implicit none
    real(dp) s
	integer mf
	integer,optional :: mfi
 
101  FORMAT(3(3x,a7,E20.13))
100  FORMAT(3x,a7,E20.13)
     mf=6
     if(present(mfi)) mf=mfi
	 WRITE(mf,100) " (0,0) ", S 

   END subroutine print_my_real

  subroutine print_my_taylor( s , mfi,title)
    implicit none
    TYPE (my_taylor) s
	integer mf,i
	logical pr
	integer,optional :: mfi
    character(*) , optional :: title
     mf=6
     if(present(mfi)) mf=mfi
     write(mf,'(a)')
     if(present(title))  then
     write(mf,'(a)') title
     else
      write(mf,*) "Variable written as an array "
     endif

	 if(first_mul) then
      call set_j_arrays 
      call multiplication_table
      first_mul=.false.
     endif
	
	pr=.false.
    call clean_taylor_PRINT(s)
101  FORMAT((3x,a1,i1,a1,i1,a1,i1,a1,E20.13))
100  FORMAT((3x,a1,i1,a1,i1,a1,i1,a1,E20.13,a5,E20.13))
           call clean(s)
     do i=0,n_mono
     if(jorder(i)>my_order) cycle
      if(abs(s%a(i))<=epsprint) cycle
     	 if(abs(imag(S%A(i)))<epsprint) then
     	  pr=.true.
  	      WRITE(mf,101) "(",jexp1(i),",",jexp2(i),",",jexp3(i),")", real(S%A(i))
	     else
     	  pr=.true.
	      WRITE(mf,100) "(",jexp1(i),",",jexp2(i),",",jexp3(i),")",  real(S%a(i),kind=dp)," + i ",aimag(S%A(i))
	     endif
     enddo
     if(.not.pr) then
       WRITE(mf,101) "(",0,",",0,",",0,")", 0.0_dp
     endif

   END subroutine print_my_taylor

  subroutine print_my_taylor2( s ,t, mfi,title)
    implicit none
    TYPE (my_taylor) s,t
	integer mf,i
	logical pr
	integer,optional :: mfi
    character(*) , optional :: title
     mf=6
     if(present(mfi)) mf=mfi
     write(mf,'(a)')
     if(present(title))  then
     write(mf,'(a)') title
     else
      write(mf,*) "Variable written as an array "
     endif

	 if(first_mul) then
      call set_j_arrays 
      call multiplication_table
      first_mul=.false.
     endif
	
	pr=.false.
     call clean_taylor_PRINT(s)
     call clean_taylor_PRINT(t)
101  FORMAT(2(3x,a1,i1,a1,i1,a1,i1,a1,E11.4))
100  FORMAT(2(3x,a1,i1,a1,i1,a1,i1,a1,E11.4,a5,E11.4))
           call clean(s);          call clean(t);
     do i=0,n_mono
     if(jorder(i)>my_order) cycle
      if(abs(s%a(i))<=epsprint.and.abs(t%a(i))<=epsprint) cycle
     	 if(ABS(imag(S%A(i)))<epsprint.and.ABS(imag(t%A(i)))<epsprint) then
     	  pr=.true.
  	      WRITE(mf,101) "(",jexp1(i),",",jexp2(i),",",jexp3(i),")", real(S%A(i)) &
            ,"(",jexp1(i),",",jexp2(i),",",jexp3(i),")", real(t%A(i))
	     else
     	  pr=.true.
	      WRITE(mf,100) "(",jexp1(i),",",jexp2(i),",",jexp3(i),")",  real(S%a(i),kind=dp)," + i ",aimag(S%A(i)) &
             ,"(",jexp1(i),",",jexp2(i),",",jexp3(i),")",  real(t%a(i),kind=dp)," + i ",aimag(t%A(i))
	     endif
     enddo
     if(.not.pr) then
       WRITE(mf,101) "(",0,",",0,",",0,")", 0.0_dp,"(",0,",",0,",",0,")", 0.0_dp
     endif

   END subroutine print_my_taylor2

  subroutine print_my_taylor4( s ,t,u,v, mfi,title)
    implicit none
    TYPE (my_taylor) s,t,u,v
	integer mf,i
	logical pr
	integer,optional :: mfi
    character(*) , optional :: title
     mf=6
     if(present(mfi)) mf=mfi
     write(mf,'(a)')
     if(present(title))  then
     write(mf,'(a)') title
     else
      write(mf,*) "Variable written as an array "
     endif

	 if(first_mul) then
      call set_j_arrays 
      call multiplication_table
      first_mul=.false.
     endif
	
	pr=.false.
     call clean_taylor_PRINT(s)
     call clean_taylor_PRINT(t)
     call clean_taylor_PRINT(u)
     call clean_taylor_PRINT(v)

101  FORMAT(4(3x,a1,i1,a1,i1,a1,i1,a1,E11.4))
100  FORMAT(4(3x,a1,i1,a1,i1,a1,i1,a1,E11.4,a5,E11.4))
           call clean(s);          call clean(t);  call clean(u);          call clean(v);
     do i=0,n_mono
     if(jorder(i)>my_order) cycle
      if(abs(s%a(i))<=epsprint.and.abs(t%a(i))<=epsprint.and.abs(u%a(i))<=epsprint.and.abs(v%a(i))<=epsprint) cycle
  if(ABS(imag(S%A(i)))<epsprint.and.ABS(imag(t%A(i)))<epsprint.and.ABS(imag(u%A(i)))<epsprint.and.ABS(imag(v%A(i)))<epsprint) then
     	  pr=.true.
  	      WRITE(mf,101) "(",jexp1(i),",",jexp2(i),",",jexp3(i),")", real(S%A(i)) &
            ,"(",jexp1(i),",",jexp2(i),",",jexp3(i),")", real(t%A(i)) &
            ,"(",jexp1(i),",",jexp2(i),",",jexp3(i),")", real(u%A(i)) &
            ,"(",jexp1(i),",",jexp2(i),",",jexp3(i),")", real(v%A(i)) 
	     else
     	  pr=.true.
	      WRITE(mf,100) "(",jexp1(i),",",jexp2(i),",",jexp3(i),")",  real(S%a(i),kind=dp)," + i ",aimag(S%A(i)) &
             ,"(",jexp1(i),",",jexp2(i),",",jexp3(i),")",  real(t%a(i),kind=dp)," + i ",aimag(t%A(i)) &
             ,"(",jexp1(i),",",jexp2(i),",",jexp3(i),")",  real(u%a(i),kind=dp)," + i ",aimag(u%A(i)) &
             ,"(",jexp1(i),",",jexp2(i),",",jexp3(i),")",  real(v%a(i),kind=dp)," + i ",aimag(v%A(i))
	     endif
     enddo
     if(.not.pr) then
       WRITE(mf,101) "(",0,",",0,",",0,")", 0.0_dp,"(",0,",",0,",",0,")", 0.0_dp
     endif

   END subroutine print_my_taylor4

  subroutine print_my_taylor_human( s , mfi,title)
    implicit none
    TYPE (my_taylor) s
	integer mf,i
	logical pr,fi
     character*3 fii
	integer,optional :: mfi
    character(*) , optional :: title
     mf=6
     if(present(mfi)) mf=mfi
     write(mf,'(a)')
     if(present(title))  then
     write(mf,'(a)') title
     else
      write(mf,*) "Variable written as Taylor series "
     endif
	 if(first_mul) then
      call set_j_arrays 
      call multiplication_table
      first_mul=.false.
     endif
	fi=.true.
    fii=' '
	pr=.false.
 
101  FORMAT((a3,E20.13,a5,a6,i1,a10,i1,a10,i1))
100  FORMAT((a3,a1,E20.13,a5,E20.13,a1,a5,a6,i1,a10,i1,a10,i1))
!FORMAT((3x,a1,i1,a1,i1,a1,i1,a1,E20.13,a5,E20.13))
           call clean(s)
     do i=0,n_mono
     if(jorder(i)>my_order) cycle
      if(abs(s%a(i))==0.0_dp) cycle


     	 if(imag(S%A(i))==0) then
     	  pr=.true.
WRITE(mf,101) fii,real(S%A(i)),"  *  ", "dx_1 ^",jexp1(i),"  * dx_2 ^",jexp2(i)," * dx_3 ^",jexp3(i) 
	     else
     	  pr=.true.
WRITE(mf,100) fii,"(",real(S%a(i),kind=dp)," + i ",aimag(S%A(i)),")","  *  ", "dx_1 ^",jexp1(i),"  * dx_2 ^",jexp2(i),"  * dx_3 ^" &
    ,jexp3(i) 
	     endif
          if(fi) then
           fii=' + '
            fi=.false.
          endif
     enddo
     if(.not.pr) then
       WRITE(mf,*)  0.0_dp
     endif

   END subroutine print_my_taylor_human
end module my_own_da


