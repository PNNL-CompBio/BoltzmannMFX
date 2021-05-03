EBASE     ?= bmx

USE_PARTICLES = TRUE
USE_MG        = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

#These are the directories in bmx
Bdirs 	:= src
Bdirs 	+= src/chemistry
Bdirs 	+= src/deposition
Bdirs 	+= src/des
Bdirs 	+= src/diffusion
Bdirs 	+= src/interpolation
Bdirs 	+= src/io
Bdirs 	+= src/leveldata
Bdirs 	+= src/mods
Bdirs 	+= src/setup
Bdirs 	+= src/timestepping
Bdirs 	+= src/util

Bpack	+= $(foreach dir, $(Bdirs), $(TOP)/$(dir)/Make.package)
Blocs	+= $(foreach dir, $(Bdirs), $(TOP)/$(dir))

include $(Bpack)
INCLUDE_LOCATIONS += $(Blocs)
VPATH_LOCATIONS   += $(Blocs)

Pdirs   := Base AmrCore Boundary Particle

Ppack	+= $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)
Plocs	+= $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir))

include $(Ppack)
INCLUDE_LOCATIONS += $(Plocs)
VPATH_LOCATIONS   += $(Plocs)

include $(AMREX_HOME)/Src/LinearSolvers/MLMG/Make.package
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/LinearSolvers/MLMG
VPATH_LOCATIONS   += $(AMREX_HOME)/Src/LinearSolvers/MLMG

#These are the directories in AMReX

all: $(executable)
	$(SILENT) $(RM) AMReX_buildInfo.cpp
	@echo SUCCESS

# job_info support
CEXE_sources += AMReX_buildInfo.cpp
CEXE_headers += $(AMREX_HOME)/Tools/C_scripts/AMReX_buildInfo.H
INCLUDE_LOCATIONS +=  $(AMREX_HOME)/Tools/C_scripts

CPPFLAGS += -DNEW_CHEM

AMReX_buildInfo.cpp:
	$(AMREX_HOME)/Tools/C_scripts/makebuildinfo_C.py \
          --amrex_home "$(AMREX_HOME)" \
          --COMP "$(COMP)" --COMP_VERSION "$(COMP_VERSION)" \
          --CXX_comp_name "$(CXX)" --CXX_flags "$(CXXFLAGS) $(CPPFLAGS) $(includes)" \
          --F_comp_name "$(F90)" --F_flags "$(F90FLAGS)" \
          --link_flags "$(LDFLAGS)" --libraries "$(libraries)" \
          --GIT "$(TOP) $(AMREX_HOME)"

vpath %.c   . $(VPATH_LOCATIONS)
vpath %.cpp . $(VPATH_LOCATIONS)
vpath %.h   . $(VPATH_LOCATIONS)
vpath %.H   . $(VPATH_LOCATIONS)
vpath %.F   . $(VPATH_LOCATIONS)
vpath %.f90 . $(VPATH_LOCATIONS)
vpath %.f   . $(VPATH_LOCATIONS)
vpath %.fi  . $(VPATH_LOCATIONS)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

clean::
	$(SILENT) $(RM) AMReX_buildInfo.cpp

