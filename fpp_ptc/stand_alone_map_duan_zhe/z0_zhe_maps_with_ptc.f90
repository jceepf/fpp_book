program spin_phase_advance_isf
use pointer_lattice
use duan_zhe_map, probe_zhe=>probe,tree_element_zhe=>tree_element,dp_zhe=>dp, & 
DEFAULT0_zhe=>DEFAULT0,TOTALPATH0_zhe=>TOTALPATH0,TIME0_zhe=>TIME0,ONLY_4d0_zhe=>ONLY_4d0,RADIATION0_zhe=>RADIATION0, &
NOCAVITY0_zhe=>NOCAVITY0,FRINGE0_zhe=>FRINGE0,STOCHASTIC0_zhe=>STOCHASTIC0,ENVELOPE0_zhe=>ENVELOPE0, &
DELTA0_zhe=>DELTA0,SPIN0_zhe=>SPIN0,MODULATION0_zhe=>MODULATION0,only_2d0_zhe=>only_2d0 , &
INTERNAL_STATE_zhe=>INTERNAL_STATE

implicit none
type(probe) xs0
type(layout), pointer :: als
integer mf,k
type(internal_state),target :: state
real(dp)  closed(6), x(6)
TYPE(tree_element_zhe)  T_zhe(3)
type(probe_zhe) xs0_zhe
type(internal_state_zhe) state_zhe
logical dofix0,dofix

    dofix=.true.
    dofix0=.true.

use_info = .true.;   ! not needed except for my windows interface


call ptc_ini_no_append


call append_empty_layout(m_u)
als=>m_u%start
call build_lattice(als,.false. ,exact=.false.,thin=.false.,onecell=.false.) 

state=radiation0+time0


allocate(my_estate)

my_estate=state
 
call read_ptc_command("map_for_zhe.txt")


x=0.0001d0
xs0=x
k=67

call propagate(als,xs0,state,fibre1=1,fibre2=k)

 call print(xs0,6)

     call zhe_ini


call read_tree_zhe(T_zhe,"map_rad7.txt")
 
    x=0.0001d0
 !x=[0.4449664625524d-09,  0.6802159036503d-07  , 0 ,      0  ,    0.3007906910141d-04 ,-0.8622664234125d-02]
    xs0_zhe=x
 
        call track_TREE_probe_complex_zhe(T_zhe,xs0_zhe,.false.)
   call print(xs0_zhe,6)


call kill_tree_zhe(t_zhe)

call read_tree_zhe(T_zhe,"one_turn_map5.txt")

    state_zhe=only_4d0_zhe
    x=0.001d0
    xs0_zhe=x
    call track_TREE_probe_complex_ptc(T_zhe,xs0_zhe,dofix0,dofix,state_zhe)  !,jump)

   call print(xs0_zhe,6)






111 call ptc_end(graphics_maybe=1,flat_file=.false.)
contains




subroutine  build_lattice(als,mis,error,exact,sl,thin,onecell)
use madx_ptc_module
use pointer_lattice
implicit none

type(layout), target :: als
real(dp),optional :: error(6)
logical, optional :: exact,sl,thin,onecell
real(dp) :: alpha,lbend, cut, ksd, ksf 
type(fibre)  l1,l2,l3,l4,l5,l6,l7,l8,l9,l10 
type(fibre)  l11,l12,l13,l14,l15,l16,l17,l18,l19,l20,cavm
type(fibre)  l21,l22,l23,l24,l25,l26,l27,l27a,l27b,l27c,l27d,ds
type(fibre)  qf1,qf2,qd1,qd2,qfa1,qfa2,sf,sd,cav,bend,vc5,bend1 
type(layout) :: sfline,sdline,sup1,supb, simple
logical(lp) :: mis,thi=.false.,oneperiod
!-----------------------------------
if(present(thin)) thi=thin

call make_states(.true.)
exact_model = .false.;oneperiod = .false.
if(present(exact)) exact_model=exact
if(present(onecell)) oneperiod=onecell
call update_states
madlength = .false.


call set_mad(energy = 1.5d0, method = 6, step = 4)

madkind2 = drift_kick_drift


l1  = drift("l1 ",  2.832695d0);l2  = drift("l2 ",  0.45698d0);
l3  = drift("l3 ",  0.08902d0);l4  = drift("l4 ",  0.2155d0);
l5  = drift("l5 ",  0.219d0);l6  = drift("l6 ",  0.107078d0);
l7  = drift("l7 ",  0.105716d0);l8  = drift("l8 ",  0.135904d0);
l9  = drift("l9 ",  0.2156993d0);l10 = drift("l10",  0.089084d0);
l11= drift("l11",  0.235416d0);l12= drift("l12",  0.1245d0);
l13= drift("l13",  0.511844d0);l14= drift("l14",  0.1788541d0);
l15= drift("l15",  0.1788483d0);l16= drift("l16",  0.511849d0);
l17= drift("l17",  0.1245d0);l18= drift("l18",  0.235405d0);
l19= drift("l19",  0.089095d0);l20= drift("l20",  0.2157007d0);
l21= drift("l21",  0.177716d0);l22= drift("l22",  0.170981d0);
l23= drift("l23",  0.218997d0);l24 = drift ("l24",  0.215503d0);
l25 = drift ("l25",  0.0890187d0);l26 = drift ("l26",  0.45698d0);
l27 = drift ("l27",  2.832696d0);l27a  = drift (" l27a",  0.8596d0);
l27b  = drift (" l27b",  0.1524d0);l27c  = drift (" l27c",  0.04445d0);
l27d  = drift (" l27d",  1.776246d0);ds  = drift (" ds  ", 0.1015d0);

qf1 = quadrupole(" qf1 ",0.344d0, k1= 2.2474d0+6.447435260914397d-03)
qf2 = quadrupole(" qf2 ",0.344d0, k1= 2.2474d0)
qd1 = quadrupole(" qd1 ",0.187d0, k1= -2.3368d0-2.593018157427161d-02); 
qd2 = quadrupole(" qd2 ",0.187d0, k1= -2.3368d0);  
qfa1= quadrupole(" qfa1",0.448d0, k1= 2.8856d0);  
qfa2= quadrupole(" qfa2",0.448d0, k1= 2.8856d0);  

!!! 1/2 mad-x value
ksf=-41.3355516397069748d0;
ksd=56.2564709584745489d0;

sf=sextupole ("sf",2.d0*0.1015d0, k2= ksf);
sd= sextupole("sd", 2.d0*0.1015d0, k2= ksd);

 vc5=marker("vc5");
alpha=0.17453292519943295769236907684886d0;
 
lbend=0.86621d0;
 
bend = rbend("bend", lbend, angle=alpha).q.(-0.778741d0)
bend1 = rbend("bend1", lbend, angle=alpha).q.(-0.778741d0)
 
cavm=mark("cavm");
cav=rfcavity("cav",l=0.0000d0,volt=-1.0d0,rev_freq=500.0d6)

if(thi) then
 sf=sextupole ("sf",0.d0, k2= ksf*0.203d0);
 sd= sextupole("sd", 0.d0, k2= ksd*0.203d0);
  sfline=(ds+sf+ds);
  sdline=(ds+sd+ds);
else
 sfline=1*sf;
 sdline=1*sd;
endif


sup1=l1+l2+l3+qf1+vc5+l4+l5+qd1+l6+l7+l8+vc5+bend+vc5+l9+sfline+l10+&
           l11+qfa1+l12+sdline+l13+ &
           l14+bend+l15+l16+sdline+l17+ &
           qfa2+l18+l19+sfline+l20+bend+l21+&
           l22+qd2+l23+l24+qf2+l25+ &
           l26+vc5+l27+cavm;


sup1=l1+l2+l3+qf1+vc5+l4+l5+qd1+l6+l7+l8+vc5+bend+vc5+l9+sfline+l10+&
           l11+qfa1+l12+sdline+l13+ &
           l14+bend+l15+l16+sdline+l17+ &
           qfa2+l18+l19+sfline+l20+bend+l21+&
           l22+qd2+l23+l24+qf2+l25+ &
           l26+vc5+l27+cavm;

supb=l1+l2+l3+qf1+vc5+l4+l5+qd1+l6+l7+l8+vc5+bend+vc5+l9+sfline+l10+&
           l11+qfa1+l12+sdline+l13+ &
           l14+bend+l15+l16+sdline+l17+ &
           qfa2+l18+l19+sfline+l20+bend1+l21+&
           l22+qd2+l23+l24+qf2+l25+ &
           l26+vc5+l27+cav;

simple=1*(qf2+l2+qd2+l2+qd2+l2); 

simple=vc5+vc5; 
!simple=15*(qf2+l2+qd2+l2+qd2+l2); 
call add(supb%start%next%next%next,4,1,1.d0)

if(oneperiod) then
 als = simple;  !11*sup1+supb;
else
 als = 11*sup1+supb;
! als = supb;
endif
if(present(sl)) then
l1  = drift("l1 ",  2.832695d0);
 if( sl ) then
  qf1 = quadrupole(" qf1 ",l=0.d0, k1= 0.01d0 ); l1  = drift("l1 ",l=0.1d0);
  als=l1+qf1;
 endif 
endif

als = .ring.als

call survey(als)


if(mis) then
 sig=1.d-5; cut=4.d0; 
 if(present(error)) sig=error
 call mess_up_alignment(als,sig,cut);
endif
end subroutine build_lattice

end program spin_phase_advance_isf




