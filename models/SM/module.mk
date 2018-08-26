DIR          := models/SM
MODNAME      := SM
SARAH_MODEL  := SM
WITH_$(MODNAME) := yes
MODSM_MOD := SM
MODSM_DEP := $(patsubst %,model_specific/%,$(MODSM_MOD))
MODSM_INC := $(patsubst %,-Imodel_specific/%,$(MODSM_MOD))
MODSM_LIB := $(foreach M,$(MODSM_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

SM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

SM_MK     := \
		$(DIR)/module.mk

SM_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

SM_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

SM_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

SM_INCLUDE_MK := \
		$(SM_SUSY_BETAS_MK) \
		$(SM_SOFT_BETAS_MK)

SM_SLHA_INPUT := \
		$(DIR)/LesHouches.in.SM_generated \
		$(DIR)/LesHouches.in.SM

SM_REFERENCES := \
		$(DIR)/SM_references.tex

SM_GNUPLOT := \
		$(DIR)/SM_plot_rgflow.gnuplot \
		$(DIR)/SM_plot_spectrum.gnuplot

SM_TARBALL := \
		$(MODNAME).tar.gz

LIBSM_SRC := \
		$(DIR)/SM_a_muon.cpp \
		$(DIR)/SM_edm.cpp \
		$(DIR)/SM_effective_couplings.cpp \
		$(DIR)/SM_info.cpp \
		$(DIR)/SM_input_parameters.cpp \
		$(DIR)/SM_mass_eigenstates.cpp \
		$(DIR)/SM_observables.cpp \
		$(DIR)/SM_physical.cpp \
		$(DIR)/SM_slha_io.cpp \
		$(DIR)/SM_soft_parameters.cpp \
		$(DIR)/SM_susy_parameters.cpp \
		$(DIR)/SM_utilities.cpp \
		$(DIR)/SM_weinberg_angle.cpp

EXESM_SRC := \
		$(DIR)/run_SM.cpp \
		$(DIR)/run_cmd_line_SM.cpp \
		$(DIR)/scan_SM.cpp
LLSM_LIB  :=
LLSM_OBJ  :=
LLSM_SRC  := \
		$(DIR)/SM_librarylink.cpp

LLSM_MMA  := \
		$(DIR)/SM_librarylink.m \
		$(DIR)/run_SM.m

LIBSM_HDR := \
		$(DIR)/SM_cxx_diagrams.hpp \
		$(DIR)/SM_a_muon.hpp \
		$(DIR)/SM_convergence_tester.hpp \
		$(DIR)/SM_edm.hpp \
		$(DIR)/SM_effective_couplings.hpp \
		$(DIR)/SM_ewsb_solver.hpp \
		$(DIR)/SM_ewsb_solver_interface.hpp \
		$(DIR)/SM_high_scale_constraint.hpp \
		$(DIR)/SM_info.hpp \
		$(DIR)/SM_initial_guesser.hpp \
		$(DIR)/SM_input_parameters.hpp \
		$(DIR)/SM_low_scale_constraint.hpp \
		$(DIR)/SM_mass_eigenstates.hpp \
		$(DIR)/SM_model.hpp \
		$(DIR)/SM_model_slha.hpp \
		$(DIR)/SM_observables.hpp \
		$(DIR)/SM_physical.hpp \
		$(DIR)/SM_slha_io.hpp \
		$(DIR)/SM_spectrum_generator.hpp \
		$(DIR)/SM_spectrum_generator_interface.hpp \
		$(DIR)/SM_soft_parameters.hpp \
		$(DIR)/SM_susy_parameters.hpp \
		$(DIR)/SM_susy_scale_constraint.hpp \
		$(DIR)/SM_utilities.hpp \
		$(DIR)/SM_weinberg_angle.hpp

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
-include $(SM_SUSY_BETAS_MK)
-include $(SM_SOFT_BETAS_MK)
-include $(SM_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(SM_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(SM_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(SM_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBSM_SRC := $(sort $(LIBSM_SRC))
EXESM_SRC := $(sort $(EXESM_SRC))

LIBSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBSM_SRC)))

EXESM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXESM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXESM_SRC)))

EXESM_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXESM_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXESM_SRC)))

LIBSM_DEP := \
		$(LIBSM_OBJ:.o=.d)

EXESM_DEP := \
		$(EXESM_OBJ:.o=.d)

LLSM_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLSM_SRC)))

LLSM_OBJ  := $(LLSM_SRC:.cpp=.o)
LLSM_LIB  := $(LLSM_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBSM     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_SM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_SM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL_DIR) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBSM) $(EXESM_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(SM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBSM_SRC) $(SM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBSM_HDR) $(SM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXESM_SRC) $(SM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLSM_SRC) $(SM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLSM_MMA) $(SM_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(SM_MK) $(SM_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(SM_INCLUDE_MK) $(SM_INSTALL_DIR)
ifneq ($(SM_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(SM_SLHA_INPUT) $(SM_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(SM_REFERENCES) $(SM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(SM_GNUPLOT) $(SM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBSM_DEP)
		-rm -f $(EXESM_DEP)
		-rm -f $(LLSM_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBSM)
		-rm -f $(LLSM_LIB)

clean-$(MODNAME)-obj:
		-rm -f $(LIBSM_OBJ)
		-rm -f $(EXESM_OBJ)
		-rm -f $(LLSM_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		-rm -f $(LIBSM_SRC)
		-rm -f $(LIBSM_HDR)
		-rm -f $(EXESM_SRC)
		-rm -f $(LLSM_SRC)
		-rm -f $(LLSM_MMA)
		-rm -f $(METACODE_STAMP_SM)
		-rm -f $(SM_INCLUDE_MK)
		-rm -f $(SM_SLHA_INPUT)
		-rm -f $(SM_REFERENCES)
		-rm -f $(SM_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		-rm -f $(EXESM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(SM_TARBALL) \
		$(LIBSM_SRC) $(LIBSM_HDR) \
		$(EXESM_SRC) \
		$(LLSM_SRC) $(LLSM_MMA) \
		$(SM_MK) $(SM_INCLUDE_MK) \
		$(SM_SLHA_INPUT) $(SM_REFERENCES) \
		$(SM_GNUPLOT)

$(LIBSM_SRC) $(LIBSM_HDR) $(EXESM_SRC) $(LLSM_SRC) $(LLSM_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_SM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_SM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_SM)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]" || (echo "Error: The code generation failed!"; exit 1)
		@touch "$(METACODE_STAMP_SM)"
		@echo "Note: to regenerate SM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_SM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_SM):
		@true
endif

$(LIBSM_DEP) $(EXESM_DEP) $(LLSM_DEP) $(LIBSM_OBJ) $(EXESM_OBJ) $(LLSM_OBJ) $(LLSM_LIB): \
	CPPFLAGS += $(MODSM_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBSM_DEP) $(EXESM_DEP) $(LLSM_DEP) $(LIBSM_OBJ) $(EXESM_OBJ) $(LLSM_OBJ) $(LLSM_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLSM_OBJ) $(LLSM_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBSM): $(LIBSM_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBSM) $(MODSM_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLSM_LIB): $(LLSM_OBJ) $(LIBSM) $(MODSM_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS)

ALLDEP += $(LIBSM_DEP) $(EXESM_DEP)
ALLSRC += $(LIBSM_SRC) $(EXESM_SRC)
ALLLIB += $(LIBSM)
ALLEXE += $(EXESM_EXE)
ALLMODDEP += $(MODSM_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLSM_DEP)
ALLSRC += $(LLSM_SRC)
ALLLL  += $(LLSM_LIB)
endif
