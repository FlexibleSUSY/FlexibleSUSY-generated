DIR          := models/CMSSMSemiAnalytic
MODNAME      := CMSSMSemiAnalytic
SARAH_MODEL  := MSSM
WITH_$(MODNAME) := yes
MODCMSSMSemiAnalytic_MOD := SM MSSM_higgs
MODCMSSMSemiAnalytic_DEP := $(patsubst %,model_specific/%,$(MODCMSSMSemiAnalytic_MOD))
MODCMSSMSemiAnalytic_INC := $(patsubst %,-Imodel_specific/%,$(MODCMSSMSemiAnalytic_MOD))
MODCMSSMSemiAnalytic_LIB := $(foreach M,$(MODCMSSMSemiAnalytic_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

MODCMSSMSemiAnalytic_SUBMOD  := $(DIR)/cxx_qft
MODCMSSMSemiAnalytic_SUBMOD_INC := $(patsubst %,-I%,$(MODCMSSMSemiAnalytic_SUBMOD))

CMSSMSemiAnalytic_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)
CMSSMSemiAnalytic_INSTALL_CXXQFT_DIR := \
		$(CMSSMSemiAnalytic_INSTALL_DIR)/cxx_qft

CMSSMSemiAnalytic_MK     := \
		$(DIR)/module.mk

CMSSMSemiAnalytic_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

CMSSMSemiAnalytic_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

CMSSMSemiAnalytic_CXX_QFT_VERTICES_MK := \
		$(DIR)/cxx_qft/vertices.mk

CMSSMSemiAnalytic_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

CMSSMSemiAnalytic_INCLUDE_MK := \
		$(CMSSMSemiAnalytic_SUSY_BETAS_MK) \
		$(CMSSMSemiAnalytic_SOFT_BETAS_MK) \
		$(CMSSMSemiAnalytic_CXX_QFT_VERTICES_MK)

CMSSMSemiAnalytic_SLHA_INPUT := \
		$(DIR)/LesHouches.in.CMSSMSemiAnalytic_generated \
		$(DIR)/LesHouches.in.CMSSMSemiAnalytic

CMSSMSemiAnalytic_REFERENCES := \
		$(DIR)/CMSSMSemiAnalytic_references.tex

CMSSMSemiAnalytic_GNUPLOT := \
		$(DIR)/CMSSMSemiAnalytic_plot_rgflow.gnuplot \
		$(DIR)/CMSSMSemiAnalytic_plot_spectrum.gnuplot

CMSSMSemiAnalytic_TARBALL := \
		$(MODNAME).tar.gz

LIBCMSSMSemiAnalytic_SRC := \
		$(DIR)/CMSSMSemiAnalytic_a_muon.cpp \
		$(DIR)/CMSSMSemiAnalytic_edm.cpp \
		$(DIR)/CMSSMSemiAnalytic_FFV_form_factors.cpp \
		$(DIR)/CMSSMSemiAnalytic_l_to_lgamma.cpp \
		$(DIR)/CMSSMSemiAnalytic_effective_couplings.cpp \
		$(DIR)/CMSSMSemiAnalytic_info.cpp \
		$(DIR)/CMSSMSemiAnalytic_input_parameters.cpp \
		$(DIR)/CMSSMSemiAnalytic_mass_eigenstates.cpp \
		$(DIR)/CMSSMSemiAnalytic_observables.cpp \
		$(DIR)/CMSSMSemiAnalytic_physical.cpp \
		$(DIR)/CMSSMSemiAnalytic_slha_io.cpp \
		$(DIR)/CMSSMSemiAnalytic_soft_parameters.cpp \
		$(DIR)/CMSSMSemiAnalytic_susy_parameters.cpp \
		$(DIR)/CMSSMSemiAnalytic_utilities.cpp \
		$(DIR)/CMSSMSemiAnalytic_weinberg_angle.cpp

EXECMSSMSemiAnalytic_SRC := \
		$(DIR)/run_CMSSMSemiAnalytic.cpp \
		$(DIR)/run_cmd_line_CMSSMSemiAnalytic.cpp \
		$(DIR)/scan_CMSSMSemiAnalytic.cpp
LLCMSSMSemiAnalytic_LIB  :=
LLCMSSMSemiAnalytic_OBJ  :=
LLCMSSMSemiAnalytic_SRC  := \
		$(DIR)/CMSSMSemiAnalytic_librarylink.cpp

LLCMSSMSemiAnalytic_MMA  := \
		$(DIR)/CMSSMSemiAnalytic_librarylink.m \
		$(DIR)/run_CMSSMSemiAnalytic.m

LIBCMSSMSemiAnalytic_HDR := \
		$(DIR)/CMSSMSemiAnalytic_a_muon.hpp \
		$(DIR)/CMSSMSemiAnalytic_convergence_tester.hpp \
		$(DIR)/CMSSMSemiAnalytic_edm.hpp \
		$(DIR)/CMSSMSemiAnalytic_FFV_form_factors.hpp \
		$(DIR)/CMSSMSemiAnalytic_l_to_lgamma.hpp \
		$(DIR)/CMSSMSemiAnalytic_effective_couplings.hpp \
		$(DIR)/CMSSMSemiAnalytic_ewsb_solver.hpp \
		$(DIR)/CMSSMSemiAnalytic_ewsb_solver_interface.hpp \
		$(DIR)/CMSSMSemiAnalytic_high_scale_constraint.hpp \
		$(DIR)/CMSSMSemiAnalytic_info.hpp \
		$(DIR)/CMSSMSemiAnalytic_initial_guesser.hpp \
		$(DIR)/CMSSMSemiAnalytic_input_parameters.hpp \
		$(DIR)/CMSSMSemiAnalytic_low_scale_constraint.hpp \
		$(DIR)/CMSSMSemiAnalytic_mass_eigenstates.hpp \
		$(DIR)/CMSSMSemiAnalytic_model.hpp \
		$(DIR)/CMSSMSemiAnalytic_model_slha.hpp \
		$(DIR)/CMSSMSemiAnalytic_observables.hpp \
		$(DIR)/CMSSMSemiAnalytic_physical.hpp \
		$(DIR)/CMSSMSemiAnalytic_slha_io.hpp \
		$(DIR)/CMSSMSemiAnalytic_spectrum_generator.hpp \
		$(DIR)/CMSSMSemiAnalytic_spectrum_generator_interface.hpp \
		$(DIR)/CMSSMSemiAnalytic_soft_parameters.hpp \
		$(DIR)/CMSSMSemiAnalytic_susy_parameters.hpp \
		$(DIR)/CMSSMSemiAnalytic_susy_scale_constraint.hpp \
		$(DIR)/CMSSMSemiAnalytic_utilities.hpp \
		$(DIR)/CMSSMSemiAnalytic_weinberg_angle.hpp

LIBCMSSMSemiAnalytic_CXXQFT_HDR := \
		$(DIR)/cxx_qft/CMSSMSemiAnalytic_qft.hpp \
		$(DIR)/cxx_qft/CMSSMSemiAnalytic_fields.hpp \
		$(DIR)/cxx_qft/CMSSMSemiAnalytic_vertices.hpp \
		$(DIR)/cxx_qft/CMSSMSemiAnalytic_context_base.hpp \
		$(DIR)/cxx_qft/CMSSMSemiAnalytic_npointfunctions.hpp

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
-include $(CMSSMSemiAnalytic_SUSY_BETAS_MK)
-include $(CMSSMSemiAnalytic_SOFT_BETAS_MK)
-include $(CMSSMSemiAnalytic_CXX_QFT_VERTICES_MK)
-include $(CMSSMSemiAnalytic_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(CMSSMSemiAnalytic_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(CMSSMSemiAnalytic_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(CMSSMSemiAnalytic_CXX_QFT_VERTICES_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(CMSSMSemiAnalytic_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBCMSSMSemiAnalytic_SRC := $(sort $(LIBCMSSMSemiAnalytic_SRC))
EXECMSSMSemiAnalytic_SRC := $(sort $(EXECMSSMSemiAnalytic_SRC))

LIBCMSSMSemiAnalytic_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBCMSSMSemiAnalytic_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBCMSSMSemiAnalytic_SRC)))

EXECMSSMSemiAnalytic_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXECMSSMSemiAnalytic_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXECMSSMSemiAnalytic_SRC)))

EXECMSSMSemiAnalytic_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXECMSSMSemiAnalytic_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXECMSSMSemiAnalytic_SRC)))

LIBCMSSMSemiAnalytic_DEP := \
		$(LIBCMSSMSemiAnalytic_OBJ:.o=.d)

EXECMSSMSemiAnalytic_DEP := \
		$(EXECMSSMSemiAnalytic_OBJ:.o=.d)

LLCMSSMSemiAnalytic_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLCMSSMSemiAnalytic_SRC)))

LLCMSSMSemiAnalytic_OBJ  := $(LLCMSSMSemiAnalytic_SRC:.cpp=.o)
LLCMSSMSemiAnalytic_LIB  := $(LLCMSSMSemiAnalytic_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBCMSSMSemiAnalytic     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_CMSSMSemiAnalytic := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_CMSSMSemiAnalytic := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL_DIR) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBCMSSMSemiAnalytic) $(EXECMSSMSemiAnalytic_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(CMSSMSemiAnalytic_INSTALL_DIR)
		$(Q)install -d $(CMSSMSemiAnalytic_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBCMSSMSemiAnalytic_SRC) $(CMSSMSemiAnalytic_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBCMSSMSemiAnalytic_HDR) $(CMSSMSemiAnalytic_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBCMSSMSemiAnalytic_CXXQFT_HDR) $(CMSSMSemiAnalytic_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(EXECMSSMSemiAnalytic_SRC) $(CMSSMSemiAnalytic_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLCMSSMSemiAnalytic_SRC) $(CMSSMSemiAnalytic_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLCMSSMSemiAnalytic_MMA) $(CMSSMSemiAnalytic_INSTALL_DIR)
		$(Q)$(INSTALL_STRIPPED) $(CMSSMSemiAnalytic_MK) $(CMSSMSemiAnalytic_INSTALL_DIR) -m u=rw,g=r,o=r
		$(Q)install -m u=rw,g=r,o=r $(CMSSMSemiAnalytic_INCLUDE_MK) $(CMSSMSemiAnalytic_INSTALL_DIR)
ifneq ($(CMSSMSemiAnalytic_SLHA_INPUT),)
		$(Q)install -m u=rw,g=r,o=r $(CMSSMSemiAnalytic_SLHA_INPUT) $(CMSSMSemiAnalytic_INSTALL_DIR)
endif
		$(Q)install -m u=rw,g=r,o=r $(CMSSMSemiAnalytic_REFERENCES) $(CMSSMSemiAnalytic_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(CMSSMSemiAnalytic_GNUPLOT) $(CMSSMSemiAnalytic_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		$(Q)-rm -f $(LIBCMSSMSemiAnalytic_DEP)
		$(Q)-rm -f $(EXECMSSMSemiAnalytic_DEP)
		$(Q)-rm -f $(LLCMSSMSemiAnalytic_DEP)

clean-$(MODNAME)-lib:
		$(Q)-rm -f $(LIBCMSSMSemiAnalytic)
		$(Q)-rm -f $(LLCMSSMSemiAnalytic_LIB)

clean-$(MODNAME)-obj:
		$(Q)-rm -f $(LIBCMSSMSemiAnalytic_OBJ)
		$(Q)-rm -f $(EXECMSSMSemiAnalytic_OBJ)
		$(Q)-rm -f $(LLCMSSMSemiAnalytic_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		$(Q)-rm -f $(LIBCMSSMSemiAnalytic_SRC)
		$(Q)-rm -f $(LIBCMSSMSemiAnalytic_HDR)
		$(Q)-rm -f $(LIBCMSSMSemiAnalytic_CXXQFT_HDR)
		$(Q)-rm -f $(EXECMSSMSemiAnalytic_SRC)
		$(Q)-rm -f $(LLCMSSMSemiAnalytic_SRC)
		$(Q)-rm -f $(LLCMSSMSemiAnalytic_MMA)
		$(Q)-rm -f $(METACODE_STAMP_CMSSMSemiAnalytic)
		$(Q)-rm -f $(CMSSMSemiAnalytic_INCLUDE_MK)
		$(Q)-rm -f $(CMSSMSemiAnalytic_SLHA_INPUT)
		$(Q)-rm -f $(CMSSMSemiAnalytic_REFERENCES)
		$(Q)-rm -f $(CMSSMSemiAnalytic_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		$(Q)-rm -f $(EXECMSSMSemiAnalytic_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		$(Q)tar -czf $(CMSSMSemiAnalytic_TARBALL) \
		$(LIBCMSSMSemiAnalytic_SRC) $(LIBCMSSMSemiAnalytic_HDR) $(LIBCMSSMSemiAnalytic_CXXQFT_HDR) \
		$(EXECMSSMSemiAnalytic_SRC) \
		$(LLCMSSMSemiAnalytic_SRC) $(LLCMSSMSemiAnalytic_MMA) \
		$(CMSSMSemiAnalytic_MK) $(CMSSMSemiAnalytic_INCLUDE_MK) \
		$(CMSSMSemiAnalytic_SLHA_INPUT) $(CMSSMSemiAnalytic_REFERENCES) \
		$(CMSSMSemiAnalytic_GNUPLOT)

$(LIBCMSSMSemiAnalytic_SRC) $(LIBCMSSMSemiAnalytic_HDR) $(LIBCMSSMSemiAnalytic_CXXQFT_HDR) $(EXECMSSMSemiAnalytic_SRC) $(LLCMSSMSemiAnalytic_SRC) $(LLCMSSMSemiAnalytic_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_CMSSMSemiAnalytic)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_CMSSMSemiAnalytic): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_CMSSMSemiAnalytic)
		@$(MSG)
		$(Q)printf "%s" "Get[\"$<\"]; Quit[]" | "$(MATH)" || (echo "Error: The code generation failed!"; exit 1)
		@touch "$(METACODE_STAMP_CMSSMSemiAnalytic)"
		@echo "Note: to regenerate CMSSMSemiAnalytic source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_CMSSMSemiAnalytic)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_CMSSMSemiAnalytic):
		@true
endif

$(LIBCMSSMSemiAnalytic_DEP) $(EXECMSSMSemiAnalytic_DEP) $(LLCMSSMSemiAnalytic_DEP) $(LIBCMSSMSemiAnalytic_OBJ) $(EXECMSSMSemiAnalytic_OBJ) $(LLCMSSMSemiAnalytic_OBJ) $(LLCMSSMSemiAnalytic_LIB): \
	CPPFLAGS += $(MODCMSSMSemiAnalytic_SUBMOD_INC) $(MODCMSSMSemiAnalytic_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBCMSSMSemiAnalytic_DEP) $(EXECMSSMSemiAnalytic_DEP) $(LLCMSSMSemiAnalytic_DEP) $(LIBCMSSMSemiAnalytic_OBJ) $(EXECMSSMSemiAnalytic_OBJ) $(LLCMSSMSemiAnalytic_OBJ) $(LLCMSSMSemiAnalytic_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLCMSSMSemiAnalytic_OBJ) $(LLCMSSMSemiAnalytic_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBCMSSMSemiAnalytic): $(LIBCMSSMSemiAnalytic_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBCMSSMSemiAnalytic) $(MODCMSSMSemiAnalytic_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		@$(MSG)
		$(Q)$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLCMSSMSemiAnalytic_LIB): $(LLCMSSMSemiAnalytic_OBJ) $(LIBCMSSMSemiAnalytic) $(MODCMSSMSemiAnalytic_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		@$(MSG)
		$(Q)$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS)

ALLDEP += $(LIBCMSSMSemiAnalytic_DEP) $(EXECMSSMSemiAnalytic_DEP)
ALLSRC += $(LIBCMSSMSemiAnalytic_SRC) $(EXECMSSMSemiAnalytic_SRC)
ALLLIB += $(LIBCMSSMSemiAnalytic)
ALLEXE += $(EXECMSSMSemiAnalytic_EXE)
ALLMODDEP += $(MODCMSSMSemiAnalytic_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLCMSSMSemiAnalytic_DEP)
ALLSRC += $(LLCMSSMSemiAnalytic_SRC)
ALLLL  += $(LLCMSSMSemiAnalytic_LIB)
endif
