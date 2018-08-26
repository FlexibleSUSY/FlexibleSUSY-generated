DIR          := models/NUTSMSSM
MODNAME      := NUTSMSSM
SARAH_MODEL  := SMSSM
WITH_$(MODNAME) := yes
MODNUTSMSSM_MOD := SM MSSM_higgs NMSSM_higgs
MODNUTSMSSM_DEP := $(patsubst %,model_specific/%,$(MODNUTSMSSM_MOD))
MODNUTSMSSM_INC := $(patsubst %,-Imodel_specific/%,$(MODNUTSMSSM_MOD))
MODNUTSMSSM_LIB := $(foreach M,$(MODNUTSMSSM_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

NUTSMSSM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

NUTSMSSM_MK     := \
		$(DIR)/module.mk

NUTSMSSM_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

NUTSMSSM_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

NUTSMSSM_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

NUTSMSSM_INCLUDE_MK := \
		$(NUTSMSSM_SUSY_BETAS_MK) \
		$(NUTSMSSM_SOFT_BETAS_MK)

NUTSMSSM_SLHA_INPUT := \
		$(DIR)/LesHouches.in.NUTSMSSM_generated \
		$(DIR)/LesHouches.in.NUTSMSSM

NUTSMSSM_REFERENCES := \
		$(DIR)/NUTSMSSM_references.tex

NUTSMSSM_GNUPLOT := \
		$(DIR)/NUTSMSSM_plot_rgflow.gnuplot \
		$(DIR)/NUTSMSSM_plot_spectrum.gnuplot

NUTSMSSM_TARBALL := \
		$(MODNAME).tar.gz

LIBNUTSMSSM_SRC := \
		$(DIR)/NUTSMSSM_a_muon.cpp \
		$(DIR)/NUTSMSSM_edm.cpp \
		$(DIR)/NUTSMSSM_effective_couplings.cpp \
		$(DIR)/NUTSMSSM_info.cpp \
		$(DIR)/NUTSMSSM_input_parameters.cpp \
		$(DIR)/NUTSMSSM_mass_eigenstates.cpp \
		$(DIR)/NUTSMSSM_observables.cpp \
		$(DIR)/NUTSMSSM_physical.cpp \
		$(DIR)/NUTSMSSM_slha_io.cpp \
		$(DIR)/NUTSMSSM_soft_parameters.cpp \
		$(DIR)/NUTSMSSM_susy_parameters.cpp \
		$(DIR)/NUTSMSSM_utilities.cpp \
		$(DIR)/NUTSMSSM_weinberg_angle.cpp

EXENUTSMSSM_SRC := \
		$(DIR)/run_NUTSMSSM.cpp \
		$(DIR)/run_cmd_line_NUTSMSSM.cpp \
		$(DIR)/scan_NUTSMSSM.cpp
LLNUTSMSSM_LIB  :=
LLNUTSMSSM_OBJ  :=
LLNUTSMSSM_SRC  := \
		$(DIR)/NUTSMSSM_librarylink.cpp

LLNUTSMSSM_MMA  := \
		$(DIR)/NUTSMSSM_librarylink.m \
		$(DIR)/run_NUTSMSSM.m

LIBNUTSMSSM_HDR := \
		$(DIR)/NUTSMSSM_cxx_diagrams.hpp \
		$(DIR)/NUTSMSSM_a_muon.hpp \
		$(DIR)/NUTSMSSM_convergence_tester.hpp \
		$(DIR)/NUTSMSSM_edm.hpp \
		$(DIR)/NUTSMSSM_effective_couplings.hpp \
		$(DIR)/NUTSMSSM_ewsb_solver.hpp \
		$(DIR)/NUTSMSSM_ewsb_solver_interface.hpp \
		$(DIR)/NUTSMSSM_high_scale_constraint.hpp \
		$(DIR)/NUTSMSSM_info.hpp \
		$(DIR)/NUTSMSSM_initial_guesser.hpp \
		$(DIR)/NUTSMSSM_input_parameters.hpp \
		$(DIR)/NUTSMSSM_low_scale_constraint.hpp \
		$(DIR)/NUTSMSSM_mass_eigenstates.hpp \
		$(DIR)/NUTSMSSM_model.hpp \
		$(DIR)/NUTSMSSM_model_slha.hpp \
		$(DIR)/NUTSMSSM_observables.hpp \
		$(DIR)/NUTSMSSM_physical.hpp \
		$(DIR)/NUTSMSSM_slha_io.hpp \
		$(DIR)/NUTSMSSM_spectrum_generator.hpp \
		$(DIR)/NUTSMSSM_spectrum_generator_interface.hpp \
		$(DIR)/NUTSMSSM_soft_parameters.hpp \
		$(DIR)/NUTSMSSM_susy_parameters.hpp \
		$(DIR)/NUTSMSSM_susy_scale_constraint.hpp \
		$(DIR)/NUTSMSSM_utilities.hpp \
		$(DIR)/NUTSMSSM_weinberg_angle.hpp

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
-include $(NUTSMSSM_SUSY_BETAS_MK)
-include $(NUTSMSSM_SOFT_BETAS_MK)
-include $(NUTSMSSM_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(NUTSMSSM_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(NUTSMSSM_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(NUTSMSSM_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBNUTSMSSM_SRC := $(sort $(LIBNUTSMSSM_SRC))
EXENUTSMSSM_SRC := $(sort $(EXENUTSMSSM_SRC))

LIBNUTSMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBNUTSMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBNUTSMSSM_SRC)))

EXENUTSMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXENUTSMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXENUTSMSSM_SRC)))

EXENUTSMSSM_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXENUTSMSSM_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXENUTSMSSM_SRC)))

LIBNUTSMSSM_DEP := \
		$(LIBNUTSMSSM_OBJ:.o=.d)

EXENUTSMSSM_DEP := \
		$(EXENUTSMSSM_OBJ:.o=.d)

LLNUTSMSSM_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLNUTSMSSM_SRC)))

LLNUTSMSSM_OBJ  := $(LLNUTSMSSM_SRC:.cpp=.o)
LLNUTSMSSM_LIB  := $(LLNUTSMSSM_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBNUTSMSSM     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_NUTSMSSM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_NUTSMSSM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL_DIR) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBNUTSMSSM) $(EXENUTSMSSM_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(NUTSMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBNUTSMSSM_SRC) $(NUTSMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBNUTSMSSM_HDR) $(NUTSMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXENUTSMSSM_SRC) $(NUTSMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLNUTSMSSM_SRC) $(NUTSMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLNUTSMSSM_MMA) $(NUTSMSSM_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(NUTSMSSM_MK) $(NUTSMSSM_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(NUTSMSSM_INCLUDE_MK) $(NUTSMSSM_INSTALL_DIR)
ifneq ($(NUTSMSSM_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(NUTSMSSM_SLHA_INPUT) $(NUTSMSSM_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(NUTSMSSM_REFERENCES) $(NUTSMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(NUTSMSSM_GNUPLOT) $(NUTSMSSM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBNUTSMSSM_DEP)
		-rm -f $(EXENUTSMSSM_DEP)
		-rm -f $(LLNUTSMSSM_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBNUTSMSSM)
		-rm -f $(LLNUTSMSSM_LIB)

clean-$(MODNAME)-obj:
		-rm -f $(LIBNUTSMSSM_OBJ)
		-rm -f $(EXENUTSMSSM_OBJ)
		-rm -f $(LLNUTSMSSM_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		-rm -f $(LIBNUTSMSSM_SRC)
		-rm -f $(LIBNUTSMSSM_HDR)
		-rm -f $(EXENUTSMSSM_SRC)
		-rm -f $(LLNUTSMSSM_SRC)
		-rm -f $(LLNUTSMSSM_MMA)
		-rm -f $(METACODE_STAMP_NUTSMSSM)
		-rm -f $(NUTSMSSM_INCLUDE_MK)
		-rm -f $(NUTSMSSM_SLHA_INPUT)
		-rm -f $(NUTSMSSM_REFERENCES)
		-rm -f $(NUTSMSSM_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		-rm -f $(EXENUTSMSSM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(NUTSMSSM_TARBALL) \
		$(LIBNUTSMSSM_SRC) $(LIBNUTSMSSM_HDR) \
		$(EXENUTSMSSM_SRC) \
		$(LLNUTSMSSM_SRC) $(LLNUTSMSSM_MMA) \
		$(NUTSMSSM_MK) $(NUTSMSSM_INCLUDE_MK) \
		$(NUTSMSSM_SLHA_INPUT) $(NUTSMSSM_REFERENCES) \
		$(NUTSMSSM_GNUPLOT)

$(LIBNUTSMSSM_SRC) $(LIBNUTSMSSM_HDR) $(EXENUTSMSSM_SRC) $(LLNUTSMSSM_SRC) $(LLNUTSMSSM_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_NUTSMSSM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_NUTSMSSM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_NUTSMSSM)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]" || (echo "Error: The code generation failed!"; exit 1)
		@touch "$(METACODE_STAMP_NUTSMSSM)"
		@echo "Note: to regenerate NUTSMSSM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_NUTSMSSM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_NUTSMSSM):
		@true
endif

$(LIBNUTSMSSM_DEP) $(EXENUTSMSSM_DEP) $(LLNUTSMSSM_DEP) $(LIBNUTSMSSM_OBJ) $(EXENUTSMSSM_OBJ) $(LLNUTSMSSM_OBJ) $(LLNUTSMSSM_LIB): \
	CPPFLAGS += $(MODNUTSMSSM_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBNUTSMSSM_DEP) $(EXENUTSMSSM_DEP) $(LLNUTSMSSM_DEP) $(LIBNUTSMSSM_OBJ) $(EXENUTSMSSM_OBJ) $(LLNUTSMSSM_OBJ) $(LLNUTSMSSM_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLNUTSMSSM_OBJ) $(LLNUTSMSSM_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBNUTSMSSM): $(LIBNUTSMSSM_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBNUTSMSSM) $(MODNUTSMSSM_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLNUTSMSSM_LIB): $(LLNUTSMSSM_OBJ) $(LIBNUTSMSSM) $(MODNUTSMSSM_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS)

ALLDEP += $(LIBNUTSMSSM_DEP) $(EXENUTSMSSM_DEP)
ALLSRC += $(LIBNUTSMSSM_SRC) $(EXENUTSMSSM_SRC)
ALLLIB += $(LIBNUTSMSSM)
ALLEXE += $(EXENUTSMSSM_EXE)
ALLMODDEP += $(MODNUTSMSSM_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLNUTSMSSM_DEP)
ALLSRC += $(LLNUTSMSSM_SRC)
ALLLL  += $(LLNUTSMSSM_LIB)
endif
