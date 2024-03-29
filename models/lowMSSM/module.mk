DIR          := models/lowMSSM
MODNAME      := lowMSSM
SARAH_MODEL  := MSSM
WITH_$(MODNAME) := yes
MODlowMSSM_MOD := SM MSSM_higgs
MODlowMSSM_DEP := $(patsubst %,model_specific/%,$(MODlowMSSM_MOD))
MODlowMSSM_INC := $(patsubst %,-Imodel_specific/%,$(MODlowMSSM_MOD))
MODlowMSSM_LIB := $(foreach M,$(MODlowMSSM_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

MODlowMSSM_SUBMOD  := $(DIR)/cxx_qft
MODlowMSSM_SUBMOD_INC := $(patsubst %,-I%,$(MODlowMSSM_SUBMOD))

lowMSSM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)
lowMSSM_INSTALL_CXXQFT_DIR := \
		$(lowMSSM_INSTALL_DIR)/cxx_qft

lowMSSM_MK     := \
		$(DIR)/module.mk

lowMSSM_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

lowMSSM_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

lowMSSM_CXXQFT_VERTICES_MK := \
		$(DIR)/cxx_qft/vertices.mk

-include $(lowMSSM_CXXQFT_VERTICES_MK)
LIBlowMSSM_CXXQFT_VERTICES_SRC ?= ''

lowMSSM_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

lowMSSM_FlexibleDecay_MK := \
		$(DIR)/decays/FlexibleDecay.mk

lowMSSM_INCLUDE_MK := \
		$(lowMSSM_SUSY_BETAS_MK) \
		$(lowMSSM_SOFT_BETAS_MK)

lowMSSM_SLHA_INPUT := \
		$(DIR)/LesHouches.in.lowMSSM_generated \
		$(DIR)/LesHouches.in.lowMSSM

lowMSSM_REFERENCES := \
		$(DIR)/lowMSSM_references.tex

lowMSSM_GNUPLOT := \
		$(DIR)/lowMSSM_plot_rgflow.gnuplot \
		$(DIR)/lowMSSM_plot_spectrum.gnuplot

lowMSSM_TARBALL := \
		$(MODNAME).tar.gz

LIBlowMSSM_SRC := \
		$(DIR)/lowMSSM_amm.cpp \
		$(DIR)/lowMSSM_edm.cpp \
		$(DIR)/lowMSSM_FFV_form_factors.cpp \
		$(wildcard $(DIR)/observables/lowMSSM*.cpp) \
		$(DIR)/lowMSSM_b_to_s_gamma.cpp \
		$(DIR)/lowMSSM_info.cpp \
		$(DIR)/lowMSSM_input_parameters.cpp \
		$(DIR)/lowMSSM_mass_eigenstates.cpp \
		$(DIR)/lowMSSM_mass_eigenstates_decoupling_scheme.cpp \
		$(DIR)/lowMSSM_model_slha.cpp \
		$(DIR)/lowMSSM_lepton_amm_wrapper.cpp \
		$(DIR)/lowMSSM_observables.cpp \
		$(DIR)/lowMSSM_physical.cpp \
		$(DIR)/lowMSSM_slha_io.cpp \
		$(DIR)/lowMSSM_soft_parameters.cpp \
		$(DIR)/lowMSSM_susy_parameters.cpp \
		$(DIR)/lowMSSM_unitarity.cpp \
		$(DIR)/lowMSSM_utilities.cpp \
		$(DIR)/lowMSSM_weinberg_angle.cpp

LIBlowMSSM_SRC += $(LIBlowMSSM_CXXQFT_VERTICES_SRC)

EXElowMSSM_SRC := \
		$(DIR)/run_lowMSSM.cpp \
		$(DIR)/run_cmd_line_lowMSSM.cpp \
		$(DIR)/scan_lowMSSM.cpp
LLlowMSSM_LIB  :=
LLlowMSSM_OBJ  :=
LLlowMSSM_SRC  := \
		$(DIR)/lowMSSM_librarylink.cpp

LLlowMSSM_MMA  := \
		$(DIR)/lowMSSM_librarylink.m \
		$(DIR)/run_lowMSSM.m

LIBlowMSSM_HDR := \
		$(DIR)/lowMSSM_amm.hpp \
		$(DIR)/lowMSSM_convergence_tester.hpp \
		$(DIR)/lowMSSM_edm.hpp \
		$(DIR)/lowMSSM_FFV_form_factors.hpp \
		$(wildcard $(DIR)/observables/lowMSSM*.hpp) \
		$(DIR)/lowMSSM_b_to_s_gamma.hpp \
		$(DIR)/lowMSSM_ewsb_solver.hpp \
		$(DIR)/lowMSSM_ewsb_solver_interface.hpp \
		$(DIR)/lowMSSM_high_scale_constraint.hpp \
		$(DIR)/lowMSSM_info.hpp \
		$(DIR)/lowMSSM_initial_guesser.hpp \
		$(DIR)/lowMSSM_input_parameters.hpp \
		$(DIR)/lowMSSM_low_scale_constraint.hpp \
		$(DIR)/lowMSSM_mass_eigenstates.hpp \
		$(DIR)/lowMSSM_mass_eigenstates_interface.hpp \
		$(DIR)/lowMSSM_mass_eigenstates_decoupling_scheme.hpp \
		$(DIR)/lowMSSM_model.hpp \
		$(DIR)/lowMSSM_model_slha.hpp \
		$(DIR)/lowMSSM_lepton_amm_wrapper.hpp \
		$(DIR)/lowMSSM_observables.hpp \
		$(DIR)/lowMSSM_physical.hpp \
		$(DIR)/lowMSSM_slha_io.hpp \
		$(DIR)/lowMSSM_spectrum_generator.hpp \
		$(DIR)/lowMSSM_spectrum_generator_interface.hpp \
		$(DIR)/lowMSSM_soft_parameters.hpp \
		$(DIR)/lowMSSM_susy_parameters.hpp \
		$(DIR)/lowMSSM_susy_scale_constraint.hpp \
		$(DIR)/lowMSSM_unitarity.hpp \
		$(DIR)/lowMSSM_utilities.hpp \
		$(DIR)/lowMSSM_weinberg_angle.hpp

LIBlowMSSM_CXXQFT_HDR := \
		$(DIR)/cxx_qft/lowMSSM_qft.hpp \
		$(DIR)/cxx_qft/lowMSSM_fields.hpp \
		$(DIR)/cxx_qft/lowMSSM_particle_aliases.hpp \
		$(DIR)/cxx_qft/lowMSSM_vertices.hpp \
		$(DIR)/cxx_qft/lowMSSM_context_base.hpp \
		$(DIR)/cxx_qft/lowMSSM_npointfunctions_wilsoncoeffs.hpp

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
-include $(lowMSSM_SUSY_BETAS_MK)
-include $(lowMSSM_SOFT_BETAS_MK)
-include $(lowMSSM_FlexibleDecay_MK)
-include $(lowMSSM_CXXQFT_VERTICES_MK)
-include $(lowMSSM_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(lowMSSM_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(lowMSSM_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@

$(lowMSSM_FlexibleDecay_MK): run-metacode-$(MODNAME)
$(lowMSSM_CXXQFT_VERTICES_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(lowMSSM_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBlowMSSM_SRC := $(sort $(LIBlowMSSM_SRC))
EXElowMSSM_SRC := $(sort $(EXElowMSSM_SRC))

LIBlowMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBlowMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBlowMSSM_SRC)))

EXElowMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXElowMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXElowMSSM_SRC)))

EXElowMSSM_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXElowMSSM_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXElowMSSM_SRC)))

LIBlowMSSM_DEP := \
		$(LIBlowMSSM_OBJ:.o=.d)

EXElowMSSM_DEP := \
		$(EXElowMSSM_OBJ:.o=.d)

LLlowMSSM_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLlowMSSM_SRC)))

LLlowMSSM_OBJ  := $(LLlowMSSM_SRC:.cpp=.o)
LLlowMSSM_LIB  := $(LLlowMSSM_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBlowMSSM     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_lowMSSM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_lowMSSM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL_DIR) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBlowMSSM) $(EXElowMSSM_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(lowMSSM_INSTALL_DIR)
		$(Q)install -d $(lowMSSM_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBlowMSSM_SRC) $(lowMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBlowMSSM_CXXQFT_VERTICES_SRC) $(lowMSSM_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBlowMSSM_HDR) $(lowMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBlowMSSM_CXXQFT_HDR) $(lowMSSM_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(EXElowMSSM_SRC) $(lowMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLlowMSSM_SRC) $(lowMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLlowMSSM_MMA) $(lowMSSM_INSTALL_DIR)
		$(Q)$(INSTALL_STRIPPED) $(lowMSSM_MK) $(lowMSSM_INSTALL_DIR) -m u=rw,g=r,o=r
		$(Q)install -m u=rw,g=r,o=r $(lowMSSM_INCLUDE_MK) $(lowMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(lowMSSM_CXXQFT_VERTICES_MK) $(lowMSSM_INSTALL_CXXQFT_DIR)

ifneq ($(lowMSSM_SLHA_INPUT),)
		$(Q)install -m u=rw,g=r,o=r $(lowMSSM_SLHA_INPUT) $(lowMSSM_INSTALL_DIR)
endif
		$(Q)install -m u=rw,g=r,o=r $(lowMSSM_REFERENCES) $(lowMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(lowMSSM_GNUPLOT) $(lowMSSM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		$(Q)-rm -f $(LIBlowMSSM_DEP)
		$(Q)-rm -f $(EXElowMSSM_DEP)
		$(Q)-rm -f $(LLlowMSSM_DEP)

clean-$(MODNAME)-lib:
		$(Q)-rm -f $(LIBlowMSSM)
		$(Q)-rm -f $(LLlowMSSM_LIB)

clean-$(MODNAME)-obj:
		$(Q)-rm -f $(LIBlowMSSM_OBJ)
		$(Q)-rm -f $(EXElowMSSM_OBJ)
		$(Q)-rm -f $(LLlowMSSM_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		$(Q)-rm -f $(LIBlowMSSM_SRC)
		$(Q)-rm -f $(LIBlowMSSM_HDR)
		$(Q)-rm -f $(LIBlowMSSM_CXXQFT_HDR)
		$(Q)-rm -f $(EXElowMSSM_SRC)
		$(Q)-rm -f $(LLlowMSSM_SRC)
		$(Q)-rm -f $(LLlowMSSM_MMA)
		$(Q)-rm -f $(METACODE_STAMP_lowMSSM)
		$(Q)-rm -f $(lowMSSM_INCLUDE_MK)
		$(Q)-rm -f $(lowMSSM_CXXQFT_VERTICES_MK)
		$(Q)-rm -f $(lowMSSM_SLHA_INPUT)
		$(Q)-rm -f $(lowMSSM_REFERENCES)
		$(Q)-rm -f $(lowMSSM_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		$(Q)-rm -f $(EXElowMSSM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		$(Q)tar -czf $(lowMSSM_TARBALL) \
		$(LIBlowMSSM_SRC) $(LIBlowMSSM_HDR) $(LIBlowMSSM_CXXQFT_HDR) \
		$(EXElowMSSM_SRC) \
		$(LLlowMSSM_SRC) $(LLlowMSSM_MMA) \
		$(lowMSSM_MK) $(lowMSSM_INCLUDE_MK) $(lowMSSM_CXXQFT_VERTICES_MK) \
		$(lowMSSM_SLHA_INPUT) $(lowMSSM_REFERENCES) \
		$(lowMSSM_GNUPLOT) \
		$(lowMSSM_FlexibleDecay_MK)

$(LIBlowMSSM_SRC) $(LIBlowMSSM_HDR) $(LIBlowMSSM_CXXQFT_HDR) $(EXElowMSSM_SRC) $(LLlowMSSM_SRC) $(LLlowMSSM_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_lowMSSM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_lowMSSM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_lowMSSM)
		@$(MSG)
		$(Q)printf "%s" "Get[\"$<\"]; Quit[]" | "$(MATH)" || (echo "Error: The code generation failed!"; exit 1)
		@touch "$(METACODE_STAMP_lowMSSM)"
		@echo "Note: to regenerate lowMSSM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_lowMSSM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_lowMSSM):
		@true
endif

$(LIBlowMSSM_DEP) $(EXElowMSSM_DEP) $(LLlowMSSM_DEP) $(LIBlowMSSM_OBJ) $(EXElowMSSM_OBJ) $(LLlowMSSM_OBJ) $(LLlowMSSM_LIB): \
	CPPFLAGS += $(MODlowMSSM_SUBMOD_INC) $(MODlowMSSM_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(GM2CALCFLAGS) $(HIGGSTOOLSFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBlowMSSM_DEP) $(EXElowMSSM_DEP) $(LLlowMSSM_DEP) $(LIBlowMSSM_OBJ) $(EXElowMSSM_OBJ) $(LLlowMSSM_OBJ) $(LLlowMSSM_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLlowMSSM_OBJ) $(LLlowMSSM_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBlowMSSM): $(LIBlowMSSM_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBlowMSSM) $(MODlowMSSM_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(FUTILIBS)
		@$(MSG)
		$(Q)$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^) $(filter -%,$(LOOPFUNCLIBS)) $(GM2CALCLIBS) $(HIGGSTOOLSLIBS) $(PYTHONLIBS) $(HIMALAYALIBS) $(GSLLIBS) $(SQLITELIBS) $(TSILLIBS) $(FLIBS) $(THREADLIBS) $(LDLIBS) $(FUTILIBS)

$(LLlowMSSM_LIB): $(LLlowMSSM_OBJ) $(LIBlowMSSM) $(MODlowMSSM_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(FUTILIBS)
		@$(MSG)
		$(Q)$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^) $(filter -%,$(LOOPFUNCLIBS)) $(GM2CALCLIBS) $(HIGGSTOOLSLIBS) $(PYTHONLIBS) $(HIMALAYALIBS) $(TSILLIBS) $(GSLLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS) $(FUTILIBS) $(FLIBS)

ALLDEP += $(LIBlowMSSM_DEP) $(EXElowMSSM_DEP)
ALLSRC += $(LIBlowMSSM_SRC) $(EXElowMSSM_SRC)
ALLLIB += $(LIBlowMSSM)
ALLEXE += $(EXElowMSSM_EXE)
ALLMODDEP += $(MODlowMSSM_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLlowMSSM_DEP)
ALLSRC += $(LLlowMSSM_SRC)
ALLLL  += $(LLlowMSSM_LIB)
endif
