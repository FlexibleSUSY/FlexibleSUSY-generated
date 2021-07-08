DIR          := models/lowNMSSM
MODNAME      := lowNMSSM
SARAH_MODEL  := NMSSM
WITH_$(MODNAME) := yes
MODlowNMSSM_MOD := SM MSSM_higgs NMSSM_higgs
MODlowNMSSM_DEP := $(patsubst %,model_specific/%,$(MODlowNMSSM_MOD))
MODlowNMSSM_INC := $(patsubst %,-Imodel_specific/%,$(MODlowNMSSM_MOD))
MODlowNMSSM_LIB := $(foreach M,$(MODlowNMSSM_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

MODlowNMSSM_SUBMOD  := $(DIR)/cxx_qft
MODlowNMSSM_SUBMOD_INC := $(patsubst %,-I%,$(MODlowNMSSM_SUBMOD))

lowNMSSM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)
lowNMSSM_INSTALL_CXXQFT_DIR := \
		$(lowNMSSM_INSTALL_DIR)/cxx_qft

lowNMSSM_MK     := \
		$(DIR)/module.mk

lowNMSSM_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

lowNMSSM_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

lowNMSSM_CXXQFT_VERTICES_MK := \
		$(DIR)/cxx_qft/vertices.mk

-include $(lowNMSSM_CXXQFT_VERTICES_MK)
LIBlowNMSSM_CXXQFT_VERTICES_SRC ?= ''

lowNMSSM_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

lowNMSSM_FlexibleDecay_MK := \
		$(DIR)/decays/FlexibleDecay.mk

lowNMSSM_INCLUDE_MK := \
		$(lowNMSSM_SUSY_BETAS_MK) \
		$(lowNMSSM_SOFT_BETAS_MK)

lowNMSSM_SLHA_INPUT := \
		$(DIR)/LesHouches.in.lowNMSSM_generated \
		$(DIR)/LesHouches.in.TP1 \
		$(DIR)/LesHouches.in.lowNMSSM \
		$(DIR)/LesHouches.in.TP3 \
		$(DIR)/LesHouches.in.TP4 \
		$(DIR)/LesHouches.in.TP2 \
		$(DIR)/LesHouches.in.TP5 \
		$(DIR)/LesHouches.in.TP6

lowNMSSM_REFERENCES := \
		$(DIR)/lowNMSSM_references.tex

lowNMSSM_GNUPLOT := \
		$(DIR)/lowNMSSM_plot_rgflow.gnuplot \
		$(DIR)/lowNMSSM_plot_spectrum.gnuplot

lowNMSSM_TARBALL := \
		$(MODNAME).tar.gz

LIBlowNMSSM_SRC := \
		$(DIR)/lowNMSSM_a_muon.cpp \
		$(DIR)/lowNMSSM_edm.cpp \
		$(DIR)/lowNMSSM_FFV_form_factors.cpp \
		$(DIR)/lowNMSSM_f_to_f_conversion.cpp \
		$(DIR)/lowNMSSM_l_to_lgamma.cpp \
		$(DIR)/lowNMSSM_b_to_s_gamma.cpp \
		$(DIR)/lowNMSSM_info.cpp \
		$(DIR)/lowNMSSM_input_parameters.cpp \
		$(DIR)/lowNMSSM_mass_eigenstates.cpp \
		$(DIR)/lowNMSSM_mass_eigenstates_decoupling_scheme.cpp \
		$(DIR)/lowNMSSM_model_slha.cpp \
		$(DIR)/lowNMSSM_observables.cpp \
		$(DIR)/lowNMSSM_physical.cpp \
		$(DIR)/lowNMSSM_slha_io.cpp \
		$(DIR)/lowNMSSM_soft_parameters.cpp \
		$(DIR)/lowNMSSM_susy_parameters.cpp \
		$(DIR)/lowNMSSM_utilities.cpp \
		$(DIR)/lowNMSSM_weinberg_angle.cpp

LIBlowNMSSM_SRC += $(LIBlowNMSSM_CXXQFT_VERTICES_SRC)

EXElowNMSSM_SRC := \
		$(DIR)/run_lowNMSSM.cpp \
		$(DIR)/run_cmd_line_lowNMSSM.cpp \
		$(DIR)/scan_lowNMSSM.cpp
LLlowNMSSM_LIB  :=
LLlowNMSSM_OBJ  :=
LLlowNMSSM_SRC  := \
		$(DIR)/lowNMSSM_librarylink.cpp

LLlowNMSSM_MMA  := \
		$(DIR)/lowNMSSM_librarylink.m \
		$(DIR)/run_lowNMSSM.m

LIBlowNMSSM_HDR := \
		$(DIR)/lowNMSSM_a_muon.hpp \
		$(DIR)/lowNMSSM_convergence_tester.hpp \
		$(DIR)/lowNMSSM_edm.hpp \
		$(DIR)/lowNMSSM_FFV_form_factors.hpp \
		$(DIR)/lowNMSSM_f_to_f_conversion.hpp \
		$(DIR)/lowNMSSM_l_to_lgamma.hpp \
		$(DIR)/lowNMSSM_b_to_s_gamma.hpp \
		$(DIR)/lowNMSSM_ewsb_solver.hpp \
		$(DIR)/lowNMSSM_ewsb_solver_interface.hpp \
		$(DIR)/lowNMSSM_high_scale_constraint.hpp \
		$(DIR)/lowNMSSM_info.hpp \
		$(DIR)/lowNMSSM_initial_guesser.hpp \
		$(DIR)/lowNMSSM_input_parameters.hpp \
		$(DIR)/lowNMSSM_low_scale_constraint.hpp \
		$(DIR)/lowNMSSM_mass_eigenstates.hpp \
		$(DIR)/lowNMSSM_mass_eigenstates_interface.hpp \
		$(DIR)/lowNMSSM_mass_eigenstates_decoupling_scheme.hpp \
		$(DIR)/lowNMSSM_model.hpp \
		$(DIR)/lowNMSSM_model_slha.hpp \
		$(DIR)/lowNMSSM_observables.hpp \
		$(DIR)/lowNMSSM_physical.hpp \
		$(DIR)/lowNMSSM_slha_io.hpp \
		$(DIR)/lowNMSSM_spectrum_generator.hpp \
		$(DIR)/lowNMSSM_spectrum_generator_interface.hpp \
		$(DIR)/lowNMSSM_soft_parameters.hpp \
		$(DIR)/lowNMSSM_susy_parameters.hpp \
		$(DIR)/lowNMSSM_susy_scale_constraint.hpp \
		$(DIR)/lowNMSSM_utilities.hpp \
		$(DIR)/lowNMSSM_weinberg_angle.hpp

LIBlowNMSSM_CXXQFT_HDR := \
		$(DIR)/cxx_qft/lowNMSSM_qft.hpp \
		$(DIR)/cxx_qft/lowNMSSM_fields.hpp \
		$(DIR)/cxx_qft/lowNMSSM_vertices.hpp \
		$(DIR)/cxx_qft/lowNMSSM_context_base.hpp \
		$(DIR)/cxx_qft/lowNMSSM_npointfunctions_wilsoncoeffs.hpp

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
-include $(lowNMSSM_SUSY_BETAS_MK)
-include $(lowNMSSM_SOFT_BETAS_MK)
-include $(lowNMSSM_FlexibleDecay_MK)
-include $(lowNMSSM_CXXQFT_VERTICES_MK)
-include $(lowNMSSM_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(lowNMSSM_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(lowNMSSM_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@

$(lowNMSSM_FlexibleDecay_MK): run-metacode-$(MODNAME)
$(lowNMSSM_CXXQFT_VERTICES_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(lowNMSSM_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBlowNMSSM_SRC := $(sort $(LIBlowNMSSM_SRC))
EXElowNMSSM_SRC := $(sort $(EXElowNMSSM_SRC))

LIBlowNMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBlowNMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBlowNMSSM_SRC)))

EXElowNMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXElowNMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXElowNMSSM_SRC)))

EXElowNMSSM_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXElowNMSSM_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXElowNMSSM_SRC)))

LIBlowNMSSM_DEP := \
		$(LIBlowNMSSM_OBJ:.o=.d)

EXElowNMSSM_DEP := \
		$(EXElowNMSSM_OBJ:.o=.d)

LLlowNMSSM_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLlowNMSSM_SRC)))

LLlowNMSSM_OBJ  := $(LLlowNMSSM_SRC:.cpp=.o)
LLlowNMSSM_LIB  := $(LLlowNMSSM_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBlowNMSSM     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_lowNMSSM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_lowNMSSM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL_DIR) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBlowNMSSM) $(EXElowNMSSM_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(lowNMSSM_INSTALL_DIR)
		$(Q)install -d $(lowNMSSM_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBlowNMSSM_SRC) $(lowNMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBlowNMSSM_CXXQFT_VERTICES_SRC) $(lowNMSSM_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBlowNMSSM_HDR) $(lowNMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBlowNMSSM_CXXQFT_HDR) $(lowNMSSM_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(EXElowNMSSM_SRC) $(lowNMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLlowNMSSM_SRC) $(lowNMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLlowNMSSM_MMA) $(lowNMSSM_INSTALL_DIR)
		$(Q)$(INSTALL_STRIPPED) $(lowNMSSM_MK) $(lowNMSSM_INSTALL_DIR) -m u=rw,g=r,o=r
		$(Q)install -m u=rw,g=r,o=r $(lowNMSSM_INCLUDE_MK) $(lowNMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(lowNMSSM_CXXQFT_VERTICES_MK) $(lowNMSSM_INSTALL_CXXQFT_DIR)

ifneq ($(lowNMSSM_SLHA_INPUT),)
		$(Q)install -m u=rw,g=r,o=r $(lowNMSSM_SLHA_INPUT) $(lowNMSSM_INSTALL_DIR)
endif
		$(Q)install -m u=rw,g=r,o=r $(lowNMSSM_REFERENCES) $(lowNMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(lowNMSSM_GNUPLOT) $(lowNMSSM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		$(Q)-rm -f $(LIBlowNMSSM_DEP)
		$(Q)-rm -f $(EXElowNMSSM_DEP)
		$(Q)-rm -f $(LLlowNMSSM_DEP)

clean-$(MODNAME)-lib:
		$(Q)-rm -f $(LIBlowNMSSM)
		$(Q)-rm -f $(LLlowNMSSM_LIB)

clean-$(MODNAME)-obj:
		$(Q)-rm -f $(LIBlowNMSSM_OBJ)
		$(Q)-rm -f $(EXElowNMSSM_OBJ)
		$(Q)-rm -f $(LLlowNMSSM_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		$(Q)-rm -f $(LIBlowNMSSM_SRC)
		$(Q)-rm -f $(LIBlowNMSSM_HDR)
		$(Q)-rm -f $(LIBlowNMSSM_CXXQFT_HDR)
		$(Q)-rm -f $(EXElowNMSSM_SRC)
		$(Q)-rm -f $(LLlowNMSSM_SRC)
		$(Q)-rm -f $(LLlowNMSSM_MMA)
		$(Q)-rm -f $(METACODE_STAMP_lowNMSSM)
		$(Q)-rm -f $(lowNMSSM_INCLUDE_MK)
		$(Q)-rm -f $(lowNMSSM_CXXQFT_VERTICES_MK)
		$(Q)-rm -f $(lowNMSSM_SLHA_INPUT)
		$(Q)-rm -f $(lowNMSSM_REFERENCES)
		$(Q)-rm -f $(lowNMSSM_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		$(Q)-rm -f $(EXElowNMSSM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		$(Q)tar -czf $(lowNMSSM_TARBALL) \
		$(LIBlowNMSSM_SRC) $(LIBlowNMSSM_HDR) $(LIBlowNMSSM_CXXQFT_HDR) \
		$(EXElowNMSSM_SRC) \
		$(LLlowNMSSM_SRC) $(LLlowNMSSM_MMA) \
		$(lowNMSSM_MK) $(lowNMSSM_INCLUDE_MK) $(lowNMSSM_CXXQFT_VERTICES_MK) \
		$(lowNMSSM_SLHA_INPUT) $(lowNMSSM_REFERENCES) \
		$(lowNMSSM_GNUPLOT) \
		$(lowNMSSM_FlexibleDecay_MK)

$(LIBlowNMSSM_SRC) $(LIBlowNMSSM_HDR) $(LIBlowNMSSM_CXXQFT_HDR) $(EXElowNMSSM_SRC) $(LLlowNMSSM_SRC) $(LLlowNMSSM_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_lowNMSSM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_lowNMSSM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_lowNMSSM)
		@$(MSG)
		$(Q)printf "%s" "Get[\"$<\"]; Quit[]" | "$(MATH)" || (echo "Error: The code generation failed!"; exit 1)
		@touch "$(METACODE_STAMP_lowNMSSM)"
		@echo "Note: to regenerate lowNMSSM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_lowNMSSM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_lowNMSSM):
		@true
endif

$(LIBlowNMSSM_DEP) $(EXElowNMSSM_DEP) $(LLlowNMSSM_DEP) $(LIBlowNMSSM_OBJ) $(EXElowNMSSM_OBJ) $(LLlowNMSSM_OBJ) $(LLlowNMSSM_LIB): \
	CPPFLAGS += $(MODlowNMSSM_SUBMOD_INC) $(MODlowNMSSM_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(GM2CALCFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBlowNMSSM_DEP) $(EXElowNMSSM_DEP) $(LLlowNMSSM_DEP) $(LIBlowNMSSM_OBJ) $(EXElowNMSSM_OBJ) $(LLlowNMSSM_OBJ) $(LLlowNMSSM_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLlowNMSSM_OBJ) $(LLlowNMSSM_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBlowNMSSM): $(LIBlowNMSSM_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBlowNMSSM) $(MODlowNMSSM_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(FUTILIBS)
		@$(MSG)
		$(Q)$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^) $(filter -%,$(LOOPFUNCLIBS)) $(GM2CALCLIBS) $(HIMALAYALIBS) $(GSLLIBS) $(SQLITELIBS) $(TSILLIBS) $(FLIBS) $(THREADLIBS) $(LDLIBS) $(FUTILIBS)

$(LLlowNMSSM_LIB): $(LLlowNMSSM_OBJ) $(LIBlowNMSSM) $(MODlowNMSSM_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(FUTILIBS)
		@$(MSG)
		$(Q)$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^) $(filter -%,$(LOOPFUNCLIBS)) $(GM2CALCLIBS) $(HIMALAYALIBS) $(TSILLIBS) $(GSLLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS) $(FUTILIBS) $(FLIBS)

ALLDEP += $(LIBlowNMSSM_DEP) $(EXElowNMSSM_DEP)
ALLSRC += $(LIBlowNMSSM_SRC) $(EXElowNMSSM_SRC)
ALLLIB += $(LIBlowNMSSM)
ALLEXE += $(EXElowNMSSM_EXE)
ALLMODDEP += $(MODlowNMSSM_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLlowNMSSM_DEP)
ALLSRC += $(LLlowNMSSM_SRC)
ALLLL  += $(LLlowNMSSM_LIB)
endif
