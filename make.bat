@ECHO off 
rem SET gf=g95
SET gf=gfortran
@echo off
SET flags= -c 
rem SET flags= -c 
rem SET flags= -c -g
rem SET flags=-c -ftrace=full
SET FPP=..\fpp_ptc
SET PTC=..\fpp_ptc
fpp_ptc
SET main=book_examples
SET out=executables_%gf%\
type %g95%
REM
ECHO on
rem cls
@type intro.txt
@ECHO off
@mkdir object_files_%gf%
@mkdir executables_%gf%
@cd book_examples
dir /B /D *
@cd ..
@type blank.txt
: start
set /p mainprogram=  Write main program name = 
@type blank.txt
if '%mainprogram%'=='' goto start
set /p choice= type : full, common, main, link, carriage return to exit = 
@cd object_files_%gf%
if '%choice%'=='FULL' goto full
if '%choice%'=='full' goto full
if '%choice%'=='link' goto link
if '%choice%'=='LINK' goto link
if '%choice%'=='main' goto main
if '%choice%'=='MAIN' goto main
if '%choice%'=='common' goto common
if '%choice%'=='COMMON' goto common
if '%choice%'=='' goto end 
@echo on
GOTO START
: full
del *.mod
del *.obj
del *.o
ECHO ON
%gf% %flags% %FPP%\a_scratch_size.f90
%gf% %flags% %FPP%\b_da_arrays_all.f90
%gf% %flags% %FPP%\c_dabnew.f90
%gf% %flags% %FPP%\d_lielib.f90
rem Complex files
%gf% %flags% %FPP%\complex\cb_da_arrays_all.f90
%gf% %flags% %FPP%\complex\cc_dabnew.f90
rem
%gf% %flags% %FPP%\h_definition.f90
%gf% %flags% %FPP%\i_tpsa.f90
%gf% %flags% %FPP%\j_tpsalie.f90
%gf% %flags% %FPP%\k_tpsalie_analysis.f90
%gf% %flags% %FPP%\l_complex_taylor.f90
%gf% %flags% %FPP%\m_real_polymorph.f90
%gf% %flags% %FPP%\n_complex_polymorph.f90
%gf% %flags% %FPP%\o_tree_element.f90
rem Complex files
%gf% %flags% %FPP%\complex\Ci_tpsa.f90
rem
%gf% %flags% %PTC%\Sa_extend_poly.f90
%gf% %flags% %PTC%\Sb_sagan_pol_arbitrary.f90
%gf% %flags% %PTC%\Sc_euclidean.f90
%gf% %flags% %PTC%\Sd_frame.f90
%gf% %flags% %PTC%\Se_status.f90
%gf% %flags% %PTC%\Sf_def_all_kinds.f90
%gf% %flags% %PTC%\Sg_sagan_wiggler.f90
%gf% %flags% %PTC%\Sh_def_kind.f90
%gf% %flags% %PTC%\Si_def_element.f90
%gf% %flags% %PTC%\Sk_link_list.f90
%gf% %flags% %PTC%\Sl_family.f90
%gf% %flags% %PTC%\Sm_tracking.f90
%gf% %flags% %PTC%\Sma0_beam_beam_ptc.f90
%gf% %flags% %PTC%\Sma_multiparticle.f90
%gf% %flags% %PTC%\Sn_mad_like.f90
%gf% %flags% %PTC%\So_fitting.f90
%gf% %flags% %PTC%\Sp_keywords.f90
%gf% %flags% %PTC%\Spb_fake_gino_sub.f90
%gf% %flags% %PTC%\Sq_orbit_ptc.f90
%gf% %flags% %PTC%\Sr_spin.f90
%gf% %flags% %PTC%\Sra_fitting.f90
%gf% %flags% %PTC%\Ss_fake_mad.f90
%gf% %flags% %PTC%\St_pointers.f90
%gf% %flags% %PTC%\zzy_run_madx.f90
: common
%gf% %flags% ..\%main%\als\z_als_lattice.f90
: main
ECHO ON
%gf% %flags% ..\%main%\%mainprogram%.f90
:link
ECHO ON
@cd ..
@MKDIR %out%%mainprogram%
@COPY /Y terminal.lnk %out%%mainprogram%\terminal.lnk
@cd object_files_%gf%
%gf% -o ..\%out%%mainprogram%\%mainprogram% St_pointers.o ^
a_scratch_size.o b_da_arrays_all.o c_dabnew.o ^
d_lielib.o h_definition.o ^
i_tpsa.o j_tpsalie.o k_tpsalie_analysis.o l_complex_taylor.o m_real_polymorph.o ^
n_complex_polymorph.o o_tree_element.o zzy_run_madx.o ^
Sa_extend_poly.o Sb_sagan_pol_arbitrary.o Sc_euclidean.o Sd_frame.o Se_status.o Sra_fitting.o ^
Sf_def_all_kinds.o Sg_sagan_wiggler.o Sh_def_kind.o Si_def_element.o Spb_fake_gino_sub.o ^
Sk_link_list.o Sl_family.o Sm_tracking.o Sma_multiparticle.o Sma0_beam_beam_ptc.o ^
Sn_mad_like.o So_fitting.o Sp_keywords.o Sr_spin.o Ss_fake_mad.o Sq_orbit_ptc.o %mainprogram%.o ^
Ci_tpsa.o cb_da_arrays_all.o cc_dabnew.o z_als_lattice.o
: end
@cd ..



