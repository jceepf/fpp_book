# |
# o---------------------------------------------------------------------o
# |
# | MAD makefile - icl/icc compiler settings
# |
# o---------------------------------------------------------------------o
# |
# | Methodical Accelerator Design
# |
# | Copyright (c) 2011+ CERN, mad@cern.ch
# |
# | For more information, see http://cern.ch/mad
# |
# o---------------------------------------------------------------------o
# |
# | $Id$
# |

#####################
# ICC specific
#

#
# preprocessor flags
#

CPPFLAGS += -D_ICC -D_ICL

#
# command translator
#

ICL_CC1 := -D%  -I% /O0
ICL_CC2 := /D%  /I% /Od

###############
# C language
#

ifeq ($(CCNAME),icl)

#
# makedep
#

ifneq ($(SED),)
CDEP = $(CC) /nologo /Zs /QMM $(addprefix /I,$(CC_DIR))
CDEP_tr = | $(SED) -e "s/$(call f2bs,$(CURDIR)/)//gi" -e "s/\.obj:/\.o:/g"
endif

#
# compiler
#

CFLAGS = /Qstd=c99 /Wall /Wcheck /Wp64 /O$(NOPT) /c

#
# diagnostics
#

CFLAGS += /D_CRT_SECURE_NO_WARNINGS /Qdiag-disable:2259,1572,981 # /Qdiag-enable:sc2

#
# options flags
#

ifeq ($(DEBUG),yes)
CFLAGS += /debug:all
endif

ifeq ($(PROFILE),yes)
CFLAGS += /Qprof-use
endif

#
# extra flags
#

CFLAGS += /nologo /Qprec /fp:strict /EHc /Qrestrict $(addprefix /I,$(CC_DIR))

#
# command translator
#

CC_tr = $(strip $(subst $(SPACE)-o , /Fo,$(call trans,$(ICL_CC1),$(ICL_CC2),$1)))

endif

###############
# C++ language
#

ifeq ($(CXXNAME),icl)

#
# makedep
#

ifneq ($(SED),)
CXXDEP = $(CXX) /nologo /Zs /QMM $(addprefix /I,$(CXX_DIR))
CXXDEP_tr = | $(SED) -e "s/$(call f2bs,$(CURDIR)/)//gi" -e "s/\.obj:/\.o:/g"
endif

#
# compiler
#

CXXFLAGS = /Qstd=c++0x /Wall /Wcheck /Wp64 /O$(NOPT) /c

#
# diagnostics
#

CXXFLAGS += /D_CRT_SECURE_NO_WARNINGS /Qdiag-disable:2259,1572,981 # /Qdiag-enable:sc2

#
# options flags
#

ifeq ($(DEBUG),yes)
CXXFLAGS += /debug:all
endif

ifeq ($(PROFILE),yes)
CXXFLAGS += /Qprof-use
endif

#
# extra flags
#

CXXFLAGS += /nologo /Qprec /fp:strict /EHc /Qrestrict $(addprefix /I,$(CXX_DIR))

#
# command translator
#

CXX_tr = $(strip $(subst $(SPACE)-o , /Fo,$(call trans,$(ICL_CC1),$(ICL_CC2),$1)))

endif

# end of makefile
