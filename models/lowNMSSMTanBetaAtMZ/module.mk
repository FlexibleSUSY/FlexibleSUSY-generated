DIR          := models/lowNMSSMTanBetaAtMZ
MODNAME      := lowNMSSMTanBetaAtMZ
SARAH_MODEL  := NMSSM
WITH_$(MODNAME) := yes
MODlowNMSSMTanBetaAtMZ_MOD := SM MSSM_higgs NMSSM_higgs
MODlowNMSSMTanBetaAtMZ_DEP := $(patsubst %,model_specific/%,$(MODlowNMSSMTanBetaAtMZ_MOD))
MODlowNMSSMTanBetaAtMZ_INC := $(patsubst %,-Imodel_specific/%,$(MODlowNMSSMTanBetaAtMZ_MOD))
MODlowNMSSMTanBetaAtMZ_LIB := $(foreach M,$(MODlowNMSSMTanBetaAtMZ_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

MODlowNMSSMTanBetaAtMZ_SUBMOD  := $(DIR)/cxx_qft
MODlowNMSSMTanBetaAtMZ_SUBMOD_INC := $(patsubst %,-I%,$(MODlowNMSSMTanBetaAtMZ_SUBMOD))

lowNMSSMTanBetaAtMZ_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)
lowNMSSMTanBetaAtMZ_INSTALL_CXXQFT_DIR := \
		$(lowNMSSMTanBetaAtMZ_INSTALL_DIR)/cxx_qft

lowNMSSMTanBetaAtMZ_MK     := \
		$(DIR)/module.mk

lowNMSSMTanBetaAtMZ_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

lowNMSSMTanBetaAtMZ_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

lowNMSSMTanBetaAtMZ_CXXQFT_VERTICES_MK := \
		$(DIR)/cxx_qft/vertices.mk

-include $(lowNMSSMTanBetaAtMZ_CXXQFT_VERTICES_MK)
LIBlowNMSSMTanBetaAtMZ_CXXQFT_VERTICES_SRC ?= ''

lowNMSSMTanBetaAtMZ_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

lowNMSSMTanBetaAtMZ_FlexibleDecay_MK := \
		$(DIR)/decays/FlexibleDecay.mk

lowNMSSMTanBetaAtMZ_INCLUDE_MK := \
		$(lowNMSSMTanBetaAtMZ_SUSY_BETAS_MK) \
		$(lowNMSSMTanBetaAtMZ_SOFT_BETAS_MK)

lowNMSSMTanBetaAtMZ_SLHA_INPUT := \
		$(DIR)/LesHouches.in.lowNMSSMTanBetaAtMZ_generated \
		$(DIR)/LesHouches.in.lowNMSSMTanBetaAtMZ \
		$(DIR)/LesHouches.in.TP4 \
		$(DIR)/LesHouches.in.TP3 \
		$(DIR)/LesHouches.in.TP5 \
		$(DIR)/LesHouches.in.TP2 \
		$(DIR)/LesHouches.in.TP6 \
		$(DIR)/LesHouches.in.TP1

lowNMSSMTanBetaAtMZ_REFERENCES := \
		$(DIR)/lowNMSSMTanBetaAtMZ_references.tex

lowNMSSMTanBetaAtMZ_GNUPLOT := \
		$(DIR)/lowNMSSMTanBetaAtMZ_plot_rgflow.gnuplot \
		$(DIR)/lowNMSSMTanBetaAtMZ_plot_spectrum.gnuplot

lowNMSSMTanBetaAtMZ_TARBALL := \
		$(MODNAME).tar.gz

LIBlowNMSSMTanBetaAtMZ_SRC := \
		$(DIR)/lowNMSSMTanBetaAtMZ_amm.cpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_edm.cpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_FFV_form_factors.cpp \
		$(wildcard $(DIR)/observables/lowNMSSMTanBetaAtMZ*.cpp) \
		$(DIR)/lowNMSSMTanBetaAtMZ_b_to_s_gamma.cpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_info.cpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_input_parameters.cpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_mass_eigenstates.cpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_mass_eigenstates_decoupling_scheme.cpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_model_slha.cpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_lepton_amm_wrapper.cpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_observables.cpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_physical.cpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_slha_io.cpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_soft_parameters.cpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_susy_parameters.cpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_unitarity.cpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_utilities.cpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_weinberg_angle.cpp

LIBlowNMSSMTanBetaAtMZ_SRC += $(LIBlowNMSSMTanBetaAtMZ_CXXQFT_VERTICES_SRC)

EXElowNMSSMTanBetaAtMZ_SRC := \
		$(DIR)/run_lowNMSSMTanBetaAtMZ.cpp \
		$(DIR)/run_cmd_line_lowNMSSMTanBetaAtMZ.cpp \
		$(DIR)/scan_lowNMSSMTanBetaAtMZ.cpp
LLlowNMSSMTanBetaAtMZ_LIB  :=
LLlowNMSSMTanBetaAtMZ_OBJ  :=
LLlowNMSSMTanBetaAtMZ_SRC  := \
		$(DIR)/lowNMSSMTanBetaAtMZ_librarylink.cpp

LLlowNMSSMTanBetaAtMZ_MMA  := \
		$(DIR)/lowNMSSMTanBetaAtMZ_librarylink.m \
		$(DIR)/run_lowNMSSMTanBetaAtMZ.m

LIBlowNMSSMTanBetaAtMZ_HDR := \
		$(DIR)/lowNMSSMTanBetaAtMZ_amm.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_convergence_tester.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_edm.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_FFV_form_factors.hpp \
		$(wildcard $(DIR)/observables/lowNMSSMTanBetaAtMZ*.hpp) \
		$(DIR)/lowNMSSMTanBetaAtMZ_b_to_s_gamma.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_ewsb_solver.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_ewsb_solver_interface.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_high_scale_constraint.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_info.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_initial_guesser.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_input_parameters.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_low_scale_constraint.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_mass_eigenstates.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_mass_eigenstates_interface.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_mass_eigenstates_decoupling_scheme.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_model.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_model_slha.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_lepton_amm_wrapper.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_observables.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_physical.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_slha_io.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_spectrum_generator.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_spectrum_generator_interface.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_soft_parameters.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_susy_parameters.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_susy_scale_constraint.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_unitarity.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_utilities.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_weinberg_angle.hpp

LIBlowNMSSMTanBetaAtMZ_CXXQFT_HDR := \
		$(DIR)/cxx_qft/lowNMSSMTanBetaAtMZ_qft.hpp \
		$(DIR)/cxx_qft/lowNMSSMTanBetaAtMZ_fields.hpp \
		$(DIR)/cxx_qft/lowNMSSMTanBetaAtMZ_particle_aliases.hpp \
		$(DIR)/cxx_qft/lowNMSSMTanBetaAtMZ_vertices.hpp \
		$(DIR)/cxx_qft/lowNMSSMTanBetaAtMZ_context_base.hpp \
		$(DIR)/cxx_qft/lowNMSSMTanBetaAtMZ_npointfunctions_wilsoncoeffs.hpp

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
-include $(lowNMSSMTanBetaAtMZ_SUSY_BETAS_MK)
-include $(lowNMSSMTanBetaAtMZ_SOFT_BETAS_MK)
-include $(lowNMSSMTanBetaAtMZ_FlexibleDecay_MK)
-include $(lowNMSSMTanBetaAtMZ_CXXQFT_VERTICES_MK)
-include $(lowNMSSMTanBetaAtMZ_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(lowNMSSMTanBetaAtMZ_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(lowNMSSMTanBetaAtMZ_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@

$(lowNMSSMTanBetaAtMZ_FlexibleDecay_MK): run-metacode-$(MODNAME)
$(lowNMSSMTanBetaAtMZ_CXXQFT_VERTICES_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(lowNMSSMTanBetaAtMZ_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBlowNMSSMTanBetaAtMZ_SRC := $(sort $(LIBlowNMSSMTanBetaAtMZ_SRC))
EXElowNMSSMTanBetaAtMZ_SRC := $(sort $(EXElowNMSSMTanBetaAtMZ_SRC))

LIBlowNMSSMTanBetaAtMZ_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBlowNMSSMTanBetaAtMZ_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBlowNMSSMTanBetaAtMZ_SRC)))

EXElowNMSSMTanBetaAtMZ_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXElowNMSSMTanBetaAtMZ_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXElowNMSSMTanBetaAtMZ_SRC)))

EXElowNMSSMTanBetaAtMZ_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXElowNMSSMTanBetaAtMZ_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXElowNMSSMTanBetaAtMZ_SRC)))

LIBlowNMSSMTanBetaAtMZ_DEP := \
		$(LIBlowNMSSMTanBetaAtMZ_OBJ:.o=.d)

EXElowNMSSMTanBetaAtMZ_DEP := \
		$(EXElowNMSSMTanBetaAtMZ_OBJ:.o=.d)

LLlowNMSSMTanBetaAtMZ_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLlowNMSSMTanBetaAtMZ_SRC)))

LLlowNMSSMTanBetaAtMZ_OBJ  := $(LLlowNMSSMTanBetaAtMZ_SRC:.cpp=.o)
LLlowNMSSMTanBetaAtMZ_LIB  := $(LLlowNMSSMTanBetaAtMZ_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBlowNMSSMTanBetaAtMZ     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_lowNMSSMTanBetaAtMZ := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_lowNMSSMTanBetaAtMZ := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL_DIR) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBlowNMSSMTanBetaAtMZ) $(EXElowNMSSMTanBetaAtMZ_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(lowNMSSMTanBetaAtMZ_INSTALL_DIR)
		$(Q)install -d $(lowNMSSMTanBetaAtMZ_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBlowNMSSMTanBetaAtMZ_SRC) $(lowNMSSMTanBetaAtMZ_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBlowNMSSMTanBetaAtMZ_CXXQFT_VERTICES_SRC) $(lowNMSSMTanBetaAtMZ_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBlowNMSSMTanBetaAtMZ_HDR) $(lowNMSSMTanBetaAtMZ_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBlowNMSSMTanBetaAtMZ_CXXQFT_HDR) $(lowNMSSMTanBetaAtMZ_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(EXElowNMSSMTanBetaAtMZ_SRC) $(lowNMSSMTanBetaAtMZ_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLlowNMSSMTanBetaAtMZ_SRC) $(lowNMSSMTanBetaAtMZ_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLlowNMSSMTanBetaAtMZ_MMA) $(lowNMSSMTanBetaAtMZ_INSTALL_DIR)
		$(Q)$(INSTALL_STRIPPED) $(lowNMSSMTanBetaAtMZ_MK) $(lowNMSSMTanBetaAtMZ_INSTALL_DIR) -m u=rw,g=r,o=r
		$(Q)install -m u=rw,g=r,o=r $(lowNMSSMTanBetaAtMZ_INCLUDE_MK) $(lowNMSSMTanBetaAtMZ_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(lowNMSSMTanBetaAtMZ_CXXQFT_VERTICES_MK) $(lowNMSSMTanBetaAtMZ_INSTALL_CXXQFT_DIR)

ifneq ($(lowNMSSMTanBetaAtMZ_SLHA_INPUT),)
		$(Q)install -m u=rw,g=r,o=r $(lowNMSSMTanBetaAtMZ_SLHA_INPUT) $(lowNMSSMTanBetaAtMZ_INSTALL_DIR)
endif
		$(Q)install -m u=rw,g=r,o=r $(lowNMSSMTanBetaAtMZ_REFERENCES) $(lowNMSSMTanBetaAtMZ_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(lowNMSSMTanBetaAtMZ_GNUPLOT) $(lowNMSSMTanBetaAtMZ_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		$(Q)-rm -f $(LIBlowNMSSMTanBetaAtMZ_DEP)
		$(Q)-rm -f $(EXElowNMSSMTanBetaAtMZ_DEP)
		$(Q)-rm -f $(LLlowNMSSMTanBetaAtMZ_DEP)

clean-$(MODNAME)-lib:
		$(Q)-rm -f $(LIBlowNMSSMTanBetaAtMZ)
		$(Q)-rm -f $(LLlowNMSSMTanBetaAtMZ_LIB)

clean-$(MODNAME)-obj:
		$(Q)-rm -f $(LIBlowNMSSMTanBetaAtMZ_OBJ)
		$(Q)-rm -f $(EXElowNMSSMTanBetaAtMZ_OBJ)
		$(Q)-rm -f $(LLlowNMSSMTanBetaAtMZ_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		$(Q)-rm -f $(LIBlowNMSSMTanBetaAtMZ_SRC)
		$(Q)-rm -f $(LIBlowNMSSMTanBetaAtMZ_HDR)
		$(Q)-rm -f $(LIBlowNMSSMTanBetaAtMZ_CXXQFT_HDR)
		$(Q)-rm -f $(EXElowNMSSMTanBetaAtMZ_SRC)
		$(Q)-rm -f $(LLlowNMSSMTanBetaAtMZ_SRC)
		$(Q)-rm -f $(LLlowNMSSMTanBetaAtMZ_MMA)
		$(Q)-rm -f $(METACODE_STAMP_lowNMSSMTanBetaAtMZ)
		$(Q)-rm -f $(lowNMSSMTanBetaAtMZ_INCLUDE_MK)
		$(Q)-rm -f $(lowNMSSMTanBetaAtMZ_CXXQFT_VERTICES_MK)
		$(Q)-rm -f $(lowNMSSMTanBetaAtMZ_SLHA_INPUT)
		$(Q)-rm -f $(lowNMSSMTanBetaAtMZ_REFERENCES)
		$(Q)-rm -f $(lowNMSSMTanBetaAtMZ_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		$(Q)-rm -f $(EXElowNMSSMTanBetaAtMZ_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		$(Q)tar -czf $(lowNMSSMTanBetaAtMZ_TARBALL) \
		$(LIBlowNMSSMTanBetaAtMZ_SRC) $(LIBlowNMSSMTanBetaAtMZ_HDR) $(LIBlowNMSSMTanBetaAtMZ_CXXQFT_HDR) \
		$(EXElowNMSSMTanBetaAtMZ_SRC) \
		$(LLlowNMSSMTanBetaAtMZ_SRC) $(LLlowNMSSMTanBetaAtMZ_MMA) \
		$(lowNMSSMTanBetaAtMZ_MK) $(lowNMSSMTanBetaAtMZ_INCLUDE_MK) $(lowNMSSMTanBetaAtMZ_CXXQFT_VERTICES_MK) \
		$(lowNMSSMTanBetaAtMZ_SLHA_INPUT) $(lowNMSSMTanBetaAtMZ_REFERENCES) \
		$(lowNMSSMTanBetaAtMZ_GNUPLOT) \
		$(lowNMSSMTanBetaAtMZ_FlexibleDecay_MK)

$(LIBlowNMSSMTanBetaAtMZ_SRC) $(LIBlowNMSSMTanBetaAtMZ_HDR) $(LIBlowNMSSMTanBetaAtMZ_CXXQFT_HDR) $(EXElowNMSSMTanBetaAtMZ_SRC) $(LLlowNMSSMTanBetaAtMZ_SRC) $(LLlowNMSSMTanBetaAtMZ_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_lowNMSSMTanBetaAtMZ)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_lowNMSSMTanBetaAtMZ): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_lowNMSSMTanBetaAtMZ)
		@$(MSG)
		$(Q)printf "%s" "Get[\"$<\"]; Quit[]" | "$(MATH)" || (echo "Error: The code generation failed!"; exit 1)
		@touch "$(METACODE_STAMP_lowNMSSMTanBetaAtMZ)"
		@echo "Note: to regenerate lowNMSSMTanBetaAtMZ source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_lowNMSSMTanBetaAtMZ)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_lowNMSSMTanBetaAtMZ):
		@true
endif

$(LIBlowNMSSMTanBetaAtMZ_DEP) $(EXElowNMSSMTanBetaAtMZ_DEP) $(LLlowNMSSMTanBetaAtMZ_DEP) $(LIBlowNMSSMTanBetaAtMZ_OBJ) $(EXElowNMSSMTanBetaAtMZ_OBJ) $(LLlowNMSSMTanBetaAtMZ_OBJ) $(LLlowNMSSMTanBetaAtMZ_LIB): \
	CPPFLAGS += $(MODlowNMSSMTanBetaAtMZ_SUBMOD_INC) $(MODlowNMSSMTanBetaAtMZ_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(GM2CALCFLAGS) $(HIGGSTOOLSFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBlowNMSSMTanBetaAtMZ_DEP) $(EXElowNMSSMTanBetaAtMZ_DEP) $(LLlowNMSSMTanBetaAtMZ_DEP) $(LIBlowNMSSMTanBetaAtMZ_OBJ) $(EXElowNMSSMTanBetaAtMZ_OBJ) $(LLlowNMSSMTanBetaAtMZ_OBJ) $(LLlowNMSSMTanBetaAtMZ_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLlowNMSSMTanBetaAtMZ_OBJ) $(LLlowNMSSMTanBetaAtMZ_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBlowNMSSMTanBetaAtMZ): $(LIBlowNMSSMTanBetaAtMZ_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBlowNMSSMTanBetaAtMZ) $(MODlowNMSSMTanBetaAtMZ_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(FUTILIBS)
		@$(MSG)
		$(Q)$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^) $(filter -%,$(LOOPFUNCLIBS)) $(GM2CALCLIBS) $(HIGGSTOOLSLIBS) $(PYTHONLIBS) $(HIMALAYALIBS) $(GSLLIBS) $(SQLITELIBS) $(TSILLIBS) $(FLIBS) $(THREADLIBS) $(LDLIBS) $(FUTILIBS)

$(LLlowNMSSMTanBetaAtMZ_LIB): $(LLlowNMSSMTanBetaAtMZ_OBJ) $(LIBlowNMSSMTanBetaAtMZ) $(MODlowNMSSMTanBetaAtMZ_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(FUTILIBS)
		@$(MSG)
		$(Q)$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^) $(filter -%,$(LOOPFUNCLIBS)) $(GM2CALCLIBS) $(HIGGSTOOLSLIBS) $(PYTHONLIBS) $(HIMALAYALIBS) $(TSILLIBS) $(GSLLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS) $(FUTILIBS) $(FLIBS)

ALLDEP += $(LIBlowNMSSMTanBetaAtMZ_DEP) $(EXElowNMSSMTanBetaAtMZ_DEP)
ALLSRC += $(LIBlowNMSSMTanBetaAtMZ_SRC) $(EXElowNMSSMTanBetaAtMZ_SRC)
ALLLIB += $(LIBlowNMSSMTanBetaAtMZ)
ALLEXE += $(EXElowNMSSMTanBetaAtMZ_EXE)
ALLMODDEP += $(MODlowNMSSMTanBetaAtMZ_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLlowNMSSMTanBetaAtMZ_DEP)
ALLSRC += $(LLlowNMSSMTanBetaAtMZ_SRC)
ALLLL  += $(LLlowNMSSMTanBetaAtMZ_LIB)
endif
