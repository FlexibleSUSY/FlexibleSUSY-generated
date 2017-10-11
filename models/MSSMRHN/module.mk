DIR          := models/MSSMRHN
MODNAME      := MSSMRHN
SARAH_MODEL  := MSSMRHN
WITH_$(MODNAME) := yes

MSSMRHN_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

MSSMRHN_MK     := \
		$(DIR)/module.mk

MSSMRHN_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

MSSMRHN_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

MSSMRHN_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

MSSMRHN_INCLUDE_MK := \
		$(MSSMRHN_SUSY_BETAS_MK) \
		$(MSSMRHN_SOFT_BETAS_MK)

MSSMRHN_SLHA_INPUT := \
		$(DIR)/LesHouches.in.MSSMRHN_generated \


MSSMRHN_GNUPLOT := \
		$(DIR)/MSSMRHN_plot_rgflow.gnuplot \
		$(DIR)/MSSMRHN_plot_spectrum.gnuplot

MSSMRHN_TARBALL := \
		$(MODNAME).tar.gz

LIBMSSMRHN_SRC := \
		$(DIR)/MSSMRHN_a_muon.cpp \
		$(DIR)/MSSMRHN_edm.cpp \
		$(DIR)/MSSMRHN_effective_couplings.cpp \
		$(DIR)/MSSMRHN_info.cpp \
		$(DIR)/MSSMRHN_input_parameters.cpp \
		$(DIR)/MSSMRHN_mass_eigenstates.cpp \
		$(DIR)/MSSMRHN_observables.cpp \
		$(DIR)/MSSMRHN_physical.cpp \
		$(DIR)/MSSMRHN_slha_io.cpp \
		$(DIR)/MSSMRHN_soft_parameters.cpp \
		$(DIR)/MSSMRHN_susy_parameters.cpp \
		$(DIR)/MSSMRHN_utilities.cpp \
		$(DIR)/MSSMRHN_weinberg_angle.cpp

EXEMSSMRHN_SRC := \
		$(DIR)/run_MSSMRHN.cpp \
		$(DIR)/run_cmd_line_MSSMRHN.cpp \
		$(DIR)/scan_MSSMRHN.cpp
LLMSSMRHN_LIB  :=
LLMSSMRHN_OBJ  :=
LLMSSMRHN_SRC  := \
		$(DIR)/MSSMRHN_librarylink.cpp

LLMSSMRHN_MMA  := \
		$(DIR)/MSSMRHN_librarylink.m \
		$(DIR)/run_MSSMRHN.m

LIBMSSMRHN_HDR := \
		$(DIR)/MSSMRHN_cxx_diagrams.hpp \
		$(DIR)/MSSMRHN_a_muon.hpp \
		$(DIR)/MSSMRHN_convergence_tester.hpp \
		$(DIR)/MSSMRHN_edm.hpp \
		$(DIR)/MSSMRHN_effective_couplings.hpp \
		$(DIR)/MSSMRHN_ewsb_solver.hpp \
		$(DIR)/MSSMRHN_ewsb_solver_interface.hpp \
		$(DIR)/MSSMRHN_high_scale_constraint.hpp \
		$(DIR)/MSSMRHN_info.hpp \
		$(DIR)/MSSMRHN_initial_guesser.hpp \
		$(DIR)/MSSMRHN_input_parameters.hpp \
		$(DIR)/MSSMRHN_low_scale_constraint.hpp \
		$(DIR)/MSSMRHN_mass_eigenstates.hpp \
		$(DIR)/MSSMRHN_model.hpp \
		$(DIR)/MSSMRHN_model_slha.hpp \
		$(DIR)/MSSMRHN_observables.hpp \
		$(DIR)/MSSMRHN_physical.hpp \
		$(DIR)/MSSMRHN_slha_io.hpp \
		$(DIR)/MSSMRHN_spectrum_generator.hpp \
		$(DIR)/MSSMRHN_spectrum_generator_interface.hpp \
		$(DIR)/MSSMRHN_soft_parameters.hpp \
		$(DIR)/MSSMRHN_susy_parameters.hpp \
		$(DIR)/MSSMRHN_susy_scale_constraint.hpp \
		$(DIR)/MSSMRHN_utilities.hpp \
		$(DIR)/MSSMRHN_weinberg_angle.hpp

ifneq ($(findstring two_scale,$(SOLVERS)),)
-include $(DIR)/two_scale.mk
endif
ifneq ($(findstring lattice,$(SOLVERS)),)
-include $(DIR)/lattice.mk
endif
ifneq ($(findstring semi_analytic,$(SOLVERS)),)
-include $(DIR)/semi_analytic.mk
endif

ifneq ($(MAKECMDGOALS),showbuild)
ifneq ($(MAKECMDGOALS),tag)
ifneq ($(MAKECMDGOALS),release)
ifneq ($(MAKECMDGOALS),doc)
-include $(MSSMRHN_SUSY_BETAS_MK)
-include $(MSSMRHN_SOFT_BETAS_MK)
-include $(MSSMRHN_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(MSSMRHN_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(MSSMRHN_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(MSSMRHN_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif

# remove duplicates in case all solvers are used
LIBMSSMRHN_SRC := $(sort $(LIBMSSMRHN_SRC))
EXEMSSMRHN_SRC := $(sort $(EXEMSSMRHN_SRC))

LIBMSSMRHN_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBMSSMRHN_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBMSSMRHN_SRC)))

EXEMSSMRHN_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEMSSMRHN_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEMSSMRHN_SRC)))

EXEMSSMRHN_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXEMSSMRHN_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXEMSSMRHN_SRC)))

LIBMSSMRHN_DEP := \
		$(LIBMSSMRHN_OBJ:.o=.d)

EXEMSSMRHN_DEP := \
		$(EXEMSSMRHN_OBJ:.o=.d)

LLMSSMRHN_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLMSSMRHN_SRC)))

LLMSSMRHN_OBJ  := $(LLMSSMRHN_SRC:.cpp=.o)
LLMSSMRHN_LIB  := $(LLMSSMRHN_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBMSSMRHN     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_MSSMRHN := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_MSSMRHN := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBMSSMRHN) $(EXEMSSMRHN_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(MSSMRHN_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBMSSMRHN_SRC) $(MSSMRHN_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBMSSMRHN_HDR) $(MSSMRHN_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXEMSSMRHN_SRC) $(MSSMRHN_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLMSSMRHN_SRC) $(MSSMRHN_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLMSSMRHN_MMA) $(MSSMRHN_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(MSSMRHN_MK) $(MSSMRHN_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(MSSMRHN_INCLUDE_MK) $(MSSMRHN_INSTALL_DIR)
ifneq ($(MSSMRHN_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(MSSMRHN_SLHA_INPUT) $(MSSMRHN_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(MSSMRHN_GNUPLOT) $(MSSMRHN_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBMSSMRHN_DEP)
		-rm -f $(EXEMSSMRHN_DEP)
		-rm -f $(LLMSSMRHN_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBMSSMRHN)
		-rm -f $(LLMSSMRHN_LIB)

clean-$(MODNAME)-obj:
		-rm -f $(LIBMSSMRHN_OBJ)
		-rm -f $(EXEMSSMRHN_OBJ)
		-rm -f $(LLMSSMRHN_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		-rm -f $(LIBMSSMRHN_SRC)
		-rm -f $(LIBMSSMRHN_HDR)
		-rm -f $(EXEMSSMRHN_SRC)
		-rm -f $(LLMSSMRHN_SRC)
		-rm -f $(LLMSSMRHN_MMA)
		-rm -f $(METACODE_STAMP_MSSMRHN)
		-rm -f $(MSSMRHN_INCLUDE_MK)
		-rm -f $(MSSMRHN_SLHA_INPUT)
		-rm -f $(MSSMRHN_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		-rm -f $(EXEMSSMRHN_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(MSSMRHN_TARBALL) \
		$(LIBMSSMRHN_SRC) $(LIBMSSMRHN_HDR) \
		$(EXEMSSMRHN_SRC) \
		$(LLMSSMRHN_SRC) $(LLMSSMRHN_MMA) \
		$(MSSMRHN_MK) $(MSSMRHN_INCLUDE_MK) \
		$(MSSMRHN_SLHA_INPUT) $(MSSMRHN_GNUPLOT)

$(LIBMSSMRHN_SRC) $(LIBMSSMRHN_HDR) $(EXEMSSMRHN_SRC) $(LLMSSMRHN_SRC) $(LLMSSMRHN_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_MSSMRHN)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_MSSMRHN): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_MSSMRHN)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_MSSMRHN)"
		@echo "Note: to regenerate MSSMRHN source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_MSSMRHN)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_MSSMRHN):
		@true
endif

$(LIBMSSMRHN_DEP) $(EXEMSSMRHN_DEP) $(LLMSSMRHN_DEP) $(LIBMSSMRHN_OBJ) $(EXEMSSMRHN_OBJ) $(LLMSSMRHN_OBJ) $(LLMSSMRHN_LIB): \
	CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBMSSMRHN_DEP) $(EXEMSSMRHN_DEP) $(LLMSSMRHN_DEP) $(LIBMSSMRHN_OBJ) $(EXEMSSMRHN_OBJ) $(LLMSSMRHN_OBJ) $(LLMSSMRHN_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLMSSMRHN_OBJ) $(LLMSSMRHN_LIB): \
	CPPFLAGS += $(shell $(MATH_INC_PATHS) --math-cmd="$(MATH)" -I --librarylink --mathlink)

$(LIBMSSMRHN): $(LIBMSSMRHN_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBMSSMRHN) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLMSSMRHN_LIB): $(LLMSSMRHN_OBJ) $(LIBMSSMRHN) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

ALLDEP += $(LIBMSSMRHN_DEP) $(EXEMSSMRHN_DEP)
ALLSRC += $(LIBMSSMRHN_SRC) $(EXEMSSMRHN_SRC)
ALLLIB += $(LIBMSSMRHN)
ALLEXE += $(EXEMSSMRHN_EXE)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLMSSMRHN_DEP)
ALLSRC += $(LLMSSMRHN_SRC)
ALLLL  += $(LLMSSMRHN_LIB)
endif
