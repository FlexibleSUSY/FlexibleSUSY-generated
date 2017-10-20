DIR          := models/SMSSM
MODNAME      := SMSSM
SARAH_MODEL  := SMSSM
WITH_$(MODNAME) := yes

SMSSM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

SMSSM_MK     := \
		$(DIR)/module.mk

SMSSM_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

SMSSM_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

SMSSM_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

SMSSM_INCLUDE_MK := \
		$(SMSSM_SUSY_BETAS_MK) \
		$(SMSSM_SOFT_BETAS_MK)

SMSSM_SLHA_INPUT := \
		$(DIR)/LesHouches.in.SMSSM_generated \
		$(DIR)/LesHouches.in.SMSSM

SMSSM_REFERENCES := \
		$(DIR)/SMSSM_references.tex

SMSSM_GNUPLOT := \
		$(DIR)/SMSSM_plot_rgflow.gnuplot \
		$(DIR)/SMSSM_plot_spectrum.gnuplot

SMSSM_TARBALL := \
		$(MODNAME).tar.gz

LIBSMSSM_SRC := \
		$(DIR)/SMSSM_a_muon.cpp \
		$(DIR)/SMSSM_edm.cpp \
		$(DIR)/SMSSM_effective_couplings.cpp \
		$(DIR)/SMSSM_info.cpp \
		$(DIR)/SMSSM_input_parameters.cpp \
		$(DIR)/SMSSM_mass_eigenstates.cpp \
		$(DIR)/SMSSM_observables.cpp \
		$(DIR)/SMSSM_physical.cpp \
		$(DIR)/SMSSM_slha_io.cpp \
		$(DIR)/SMSSM_soft_parameters.cpp \
		$(DIR)/SMSSM_susy_parameters.cpp \
		$(DIR)/SMSSM_utilities.cpp \
		$(DIR)/SMSSM_weinberg_angle.cpp

EXESMSSM_SRC := \
		$(DIR)/run_SMSSM.cpp \
		$(DIR)/run_cmd_line_SMSSM.cpp \
		$(DIR)/scan_SMSSM.cpp
LLSMSSM_LIB  :=
LLSMSSM_OBJ  :=
LLSMSSM_SRC  := \
		$(DIR)/SMSSM_librarylink.cpp

LLSMSSM_MMA  := \
		$(DIR)/SMSSM_librarylink.m \
		$(DIR)/run_SMSSM.m

LIBSMSSM_HDR := \
		$(DIR)/SMSSM_cxx_diagrams.hpp \
		$(DIR)/SMSSM_a_muon.hpp \
		$(DIR)/SMSSM_convergence_tester.hpp \
		$(DIR)/SMSSM_edm.hpp \
		$(DIR)/SMSSM_effective_couplings.hpp \
		$(DIR)/SMSSM_ewsb_solver.hpp \
		$(DIR)/SMSSM_ewsb_solver_interface.hpp \
		$(DIR)/SMSSM_high_scale_constraint.hpp \
		$(DIR)/SMSSM_info.hpp \
		$(DIR)/SMSSM_initial_guesser.hpp \
		$(DIR)/SMSSM_input_parameters.hpp \
		$(DIR)/SMSSM_low_scale_constraint.hpp \
		$(DIR)/SMSSM_mass_eigenstates.hpp \
		$(DIR)/SMSSM_model.hpp \
		$(DIR)/SMSSM_model_slha.hpp \
		$(DIR)/SMSSM_observables.hpp \
		$(DIR)/SMSSM_physical.hpp \
		$(DIR)/SMSSM_slha_io.hpp \
		$(DIR)/SMSSM_spectrum_generator.hpp \
		$(DIR)/SMSSM_spectrum_generator_interface.hpp \
		$(DIR)/SMSSM_soft_parameters.hpp \
		$(DIR)/SMSSM_susy_parameters.hpp \
		$(DIR)/SMSSM_susy_scale_constraint.hpp \
		$(DIR)/SMSSM_utilities.hpp \
		$(DIR)/SMSSM_weinberg_angle.hpp

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
-include $(SMSSM_SUSY_BETAS_MK)
-include $(SMSSM_SOFT_BETAS_MK)
-include $(SMSSM_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(SMSSM_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(SMSSM_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(SMSSM_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBSMSSM_SRC := $(sort $(LIBSMSSM_SRC))
EXESMSSM_SRC := $(sort $(EXESMSSM_SRC))

LIBSMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBSMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBSMSSM_SRC)))

EXESMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXESMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXESMSSM_SRC)))

EXESMSSM_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXESMSSM_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXESMSSM_SRC)))

LIBSMSSM_DEP := \
		$(LIBSMSSM_OBJ:.o=.d)

EXESMSSM_DEP := \
		$(EXESMSSM_OBJ:.o=.d)

LLSMSSM_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLSMSSM_SRC)))

LLSMSSM_OBJ  := $(LLSMSSM_SRC:.cpp=.o)
LLSMSSM_LIB  := $(LLSMSSM_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBSMSSM     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_SMSSM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_SMSSM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBSMSSM) $(EXESMSSM_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(SMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBSMSSM_SRC) $(SMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBSMSSM_HDR) $(SMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXESMSSM_SRC) $(SMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLSMSSM_SRC) $(SMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLSMSSM_MMA) $(SMSSM_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(SMSSM_MK) $(SMSSM_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(SMSSM_INCLUDE_MK) $(SMSSM_INSTALL_DIR)
ifneq ($(SMSSM_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(SMSSM_SLHA_INPUT) $(SMSSM_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(SMSSM_REFERENCES) $(SMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(SMSSM_GNUPLOT) $(SMSSM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBSMSSM_DEP)
		-rm -f $(EXESMSSM_DEP)
		-rm -f $(LLSMSSM_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBSMSSM)
		-rm -f $(LLSMSSM_LIB)

clean-$(MODNAME)-obj:
		-rm -f $(LIBSMSSM_OBJ)
		-rm -f $(EXESMSSM_OBJ)
		-rm -f $(LLSMSSM_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		-rm -f $(LIBSMSSM_SRC)
		-rm -f $(LIBSMSSM_HDR)
		-rm -f $(EXESMSSM_SRC)
		-rm -f $(LLSMSSM_SRC)
		-rm -f $(LLSMSSM_MMA)
		-rm -f $(METACODE_STAMP_SMSSM)
		-rm -f $(SMSSM_INCLUDE_MK)
		-rm -f $(SMSSM_SLHA_INPUT)
		-rm -f $(SMSSM_REFERENCES)
		-rm -f $(SMSSM_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		-rm -f $(EXESMSSM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(SMSSM_TARBALL) \
		$(LIBSMSSM_SRC) $(LIBSMSSM_HDR) \
		$(EXESMSSM_SRC) \
		$(LLSMSSM_SRC) $(LLSMSSM_MMA) \
		$(SMSSM_MK) $(SMSSM_INCLUDE_MK) \
		$(SMSSM_SLHA_INPUT) $(SMSSM_REFERENCES) \
		$(SMSSM_GNUPLOT)

$(LIBSMSSM_SRC) $(LIBSMSSM_HDR) $(EXESMSSM_SRC) $(LLSMSSM_SRC) $(LLSMSSM_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_SMSSM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_SMSSM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_SMSSM)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_SMSSM)"
		@echo "Note: to regenerate SMSSM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_SMSSM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_SMSSM):
		@true
endif

$(LIBSMSSM_DEP) $(EXESMSSM_DEP) $(LLSMSSM_DEP) $(LIBSMSSM_OBJ) $(EXESMSSM_OBJ) $(LLSMSSM_OBJ) $(LLSMSSM_LIB): \
	CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBSMSSM_DEP) $(EXESMSSM_DEP) $(LLSMSSM_DEP) $(LIBSMSSM_OBJ) $(EXESMSSM_OBJ) $(LLSMSSM_OBJ) $(LLSMSSM_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLSMSSM_OBJ) $(LLSMSSM_LIB): \
	CPPFLAGS += $(shell $(MATH_INC_PATHS) --math-cmd="$(MATH)" -I --librarylink --mathlink)

$(LIBSMSSM): $(LIBSMSSM_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBSMSSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLSMSSM_LIB): $(LLSMSSM_OBJ) $(LIBSMSSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

ALLDEP += $(LIBSMSSM_DEP) $(EXESMSSM_DEP)
ALLSRC += $(LIBSMSSM_SRC) $(EXESMSSM_SRC)
ALLLIB += $(LIBSMSSM)
ALLEXE += $(EXESMSSM_EXE)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLSMSSM_DEP)
ALLSRC += $(LLSMSSM_SRC)
ALLLL  += $(LLSMSSM_LIB)
endif
