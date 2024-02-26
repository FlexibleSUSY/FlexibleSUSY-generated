DIR          := models/MSSMNoFVHimalaya
MODNAME      := MSSMNoFVHimalaya
SARAH_MODEL  := MSSMNoFV
WITH_$(MODNAME) := yes
MODMSSMNoFVHimalaya_MOD := SM MSSM_higgs MSSM_thresholds
MODMSSMNoFVHimalaya_DEP := $(patsubst %,model_specific/%,$(MODMSSMNoFVHimalaya_MOD))
MODMSSMNoFVHimalaya_INC := $(patsubst %,-Imodel_specific/%,$(MODMSSMNoFVHimalaya_MOD))
MODMSSMNoFVHimalaya_LIB := $(foreach M,$(MODMSSMNoFVHimalaya_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

MODMSSMNoFVHimalaya_SUBMOD  := $(DIR)/cxx_qft
MODMSSMNoFVHimalaya_SUBMOD_INC := $(patsubst %,-I%,$(MODMSSMNoFVHimalaya_SUBMOD))

MSSMNoFVHimalaya_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)
MSSMNoFVHimalaya_INSTALL_CXXQFT_DIR := \
		$(MSSMNoFVHimalaya_INSTALL_DIR)/cxx_qft

MSSMNoFVHimalaya_MK     := \
		$(DIR)/module.mk

MSSMNoFVHimalaya_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

MSSMNoFVHimalaya_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

MSSMNoFVHimalaya_CXXQFT_VERTICES_MK := \
		$(DIR)/cxx_qft/vertices.mk

-include $(MSSMNoFVHimalaya_CXXQFT_VERTICES_MK)
LIBMSSMNoFVHimalaya_CXXQFT_VERTICES_SRC ?= ''

MSSMNoFVHimalaya_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

MSSMNoFVHimalaya_FlexibleDecay_MK := \
		$(DIR)/decays/FlexibleDecay.mk

MSSMNoFVHimalaya_INCLUDE_MK := \
		$(MSSMNoFVHimalaya_SUSY_BETAS_MK) \
		$(MSSMNoFVHimalaya_SOFT_BETAS_MK)

MSSMNoFVHimalaya_SLHA_INPUT := \
		$(DIR)/LesHouches.in.MSSMNoFVHimalaya_generated \
		$(DIR)/LesHouches.in.MSSMNoFVHimalaya

MSSMNoFVHimalaya_REFERENCES := \
		$(DIR)/MSSMNoFVHimalaya_references.tex

MSSMNoFVHimalaya_GNUPLOT := \
		$(DIR)/MSSMNoFVHimalaya_plot_rgflow.gnuplot \
		$(DIR)/MSSMNoFVHimalaya_plot_spectrum.gnuplot

MSSMNoFVHimalaya_TARBALL := \
		$(MODNAME).tar.gz

LIBMSSMNoFVHimalaya_SRC := \
		$(DIR)/MSSMNoFVHimalaya_amm.cpp \
		$(DIR)/MSSMNoFVHimalaya_edm.cpp \
		$(DIR)/MSSMNoFVHimalaya_FFV_form_factors.cpp \
		$(wildcard $(DIR)/observables/MSSMNoFVHimalaya*.cpp) \
		$(DIR)/MSSMNoFVHimalaya_b_to_s_gamma.cpp \
		$(DIR)/MSSMNoFVHimalaya_info.cpp \
		$(DIR)/MSSMNoFVHimalaya_input_parameters.cpp \
		$(DIR)/MSSMNoFVHimalaya_mass_eigenstates.cpp \
		$(DIR)/MSSMNoFVHimalaya_mass_eigenstates_decoupling_scheme.cpp \
		$(DIR)/MSSMNoFVHimalaya_model_slha.cpp \
		$(DIR)/MSSMNoFVHimalaya_lepton_amm_wrapper.cpp \
		$(DIR)/MSSMNoFVHimalaya_observables.cpp \
		$(DIR)/MSSMNoFVHimalaya_physical.cpp \
		$(DIR)/MSSMNoFVHimalaya_slha_io.cpp \
		$(DIR)/MSSMNoFVHimalaya_soft_parameters.cpp \
		$(DIR)/MSSMNoFVHimalaya_susy_parameters.cpp \
		$(DIR)/MSSMNoFVHimalaya_unitarity.cpp \
		$(DIR)/MSSMNoFVHimalaya_utilities.cpp \
		$(DIR)/MSSMNoFVHimalaya_weinberg_angle.cpp

LIBMSSMNoFVHimalaya_SRC += $(LIBMSSMNoFVHimalaya_CXXQFT_VERTICES_SRC)

EXEMSSMNoFVHimalaya_SRC := \
		$(DIR)/run_MSSMNoFVHimalaya.cpp \
		$(DIR)/run_cmd_line_MSSMNoFVHimalaya.cpp \
		$(DIR)/scan_MSSMNoFVHimalaya.cpp
LLMSSMNoFVHimalaya_LIB  :=
LLMSSMNoFVHimalaya_OBJ  :=
LLMSSMNoFVHimalaya_SRC  := \
		$(DIR)/MSSMNoFVHimalaya_librarylink.cpp

LLMSSMNoFVHimalaya_MMA  := \
		$(DIR)/MSSMNoFVHimalaya_librarylink.m \
		$(DIR)/run_MSSMNoFVHimalaya.m

LIBMSSMNoFVHimalaya_HDR := \
		$(DIR)/MSSMNoFVHimalaya_amm.hpp \
		$(DIR)/MSSMNoFVHimalaya_convergence_tester.hpp \
		$(DIR)/MSSMNoFVHimalaya_edm.hpp \
		$(DIR)/MSSMNoFVHimalaya_FFV_form_factors.hpp \
		$(wildcard $(DIR)/observables/MSSMNoFVHimalaya*.hpp) \
		$(DIR)/MSSMNoFVHimalaya_b_to_s_gamma.hpp \
		$(DIR)/MSSMNoFVHimalaya_ewsb_solver.hpp \
		$(DIR)/MSSMNoFVHimalaya_ewsb_solver_interface.hpp \
		$(DIR)/MSSMNoFVHimalaya_high_scale_constraint.hpp \
		$(DIR)/MSSMNoFVHimalaya_info.hpp \
		$(DIR)/MSSMNoFVHimalaya_initial_guesser.hpp \
		$(DIR)/MSSMNoFVHimalaya_input_parameters.hpp \
		$(DIR)/MSSMNoFVHimalaya_low_scale_constraint.hpp \
		$(DIR)/MSSMNoFVHimalaya_mass_eigenstates.hpp \
		$(DIR)/MSSMNoFVHimalaya_mass_eigenstates_interface.hpp \
		$(DIR)/MSSMNoFVHimalaya_mass_eigenstates_decoupling_scheme.hpp \
		$(DIR)/MSSMNoFVHimalaya_model.hpp \
		$(DIR)/MSSMNoFVHimalaya_model_slha.hpp \
		$(DIR)/MSSMNoFVHimalaya_lepton_amm_wrapper.hpp \
		$(DIR)/MSSMNoFVHimalaya_observables.hpp \
		$(DIR)/MSSMNoFVHimalaya_physical.hpp \
		$(DIR)/MSSMNoFVHimalaya_slha_io.hpp \
		$(DIR)/MSSMNoFVHimalaya_spectrum_generator.hpp \
		$(DIR)/MSSMNoFVHimalaya_spectrum_generator_interface.hpp \
		$(DIR)/MSSMNoFVHimalaya_soft_parameters.hpp \
		$(DIR)/MSSMNoFVHimalaya_susy_parameters.hpp \
		$(DIR)/MSSMNoFVHimalaya_susy_scale_constraint.hpp \
		$(DIR)/MSSMNoFVHimalaya_unitarity.hpp \
		$(DIR)/MSSMNoFVHimalaya_utilities.hpp \
		$(DIR)/MSSMNoFVHimalaya_weinberg_angle.hpp

LIBMSSMNoFVHimalaya_CXXQFT_HDR := \
		$(DIR)/cxx_qft/MSSMNoFVHimalaya_qft.hpp \
		$(DIR)/cxx_qft/MSSMNoFVHimalaya_fields.hpp \
		$(DIR)/cxx_qft/MSSMNoFVHimalaya_particle_aliases.hpp \
		$(DIR)/cxx_qft/MSSMNoFVHimalaya_vertices.hpp \
		$(DIR)/cxx_qft/MSSMNoFVHimalaya_context_base.hpp \
		$(DIR)/cxx_qft/MSSMNoFVHimalaya_npointfunctions_wilsoncoeffs.hpp

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
-include $(MSSMNoFVHimalaya_SUSY_BETAS_MK)
-include $(MSSMNoFVHimalaya_SOFT_BETAS_MK)
-include $(MSSMNoFVHimalaya_FlexibleDecay_MK)
-include $(MSSMNoFVHimalaya_CXXQFT_VERTICES_MK)
-include $(MSSMNoFVHimalaya_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(MSSMNoFVHimalaya_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(MSSMNoFVHimalaya_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@

$(MSSMNoFVHimalaya_FlexibleDecay_MK): run-metacode-$(MODNAME)
$(MSSMNoFVHimalaya_CXXQFT_VERTICES_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(MSSMNoFVHimalaya_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBMSSMNoFVHimalaya_SRC := $(sort $(LIBMSSMNoFVHimalaya_SRC))
EXEMSSMNoFVHimalaya_SRC := $(sort $(EXEMSSMNoFVHimalaya_SRC))

LIBMSSMNoFVHimalaya_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBMSSMNoFVHimalaya_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBMSSMNoFVHimalaya_SRC)))

EXEMSSMNoFVHimalaya_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEMSSMNoFVHimalaya_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEMSSMNoFVHimalaya_SRC)))

EXEMSSMNoFVHimalaya_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXEMSSMNoFVHimalaya_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXEMSSMNoFVHimalaya_SRC)))

LIBMSSMNoFVHimalaya_DEP := \
		$(LIBMSSMNoFVHimalaya_OBJ:.o=.d)

EXEMSSMNoFVHimalaya_DEP := \
		$(EXEMSSMNoFVHimalaya_OBJ:.o=.d)

LLMSSMNoFVHimalaya_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLMSSMNoFVHimalaya_SRC)))

LLMSSMNoFVHimalaya_OBJ  := $(LLMSSMNoFVHimalaya_SRC:.cpp=.o)
LLMSSMNoFVHimalaya_LIB  := $(LLMSSMNoFVHimalaya_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBMSSMNoFVHimalaya     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_MSSMNoFVHimalaya := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_MSSMNoFVHimalaya := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL_DIR) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBMSSMNoFVHimalaya) $(EXEMSSMNoFVHimalaya_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(MSSMNoFVHimalaya_INSTALL_DIR)
		$(Q)install -d $(MSSMNoFVHimalaya_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBMSSMNoFVHimalaya_SRC) $(MSSMNoFVHimalaya_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBMSSMNoFVHimalaya_CXXQFT_VERTICES_SRC) $(MSSMNoFVHimalaya_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBMSSMNoFVHimalaya_HDR) $(MSSMNoFVHimalaya_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBMSSMNoFVHimalaya_CXXQFT_HDR) $(MSSMNoFVHimalaya_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(EXEMSSMNoFVHimalaya_SRC) $(MSSMNoFVHimalaya_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLMSSMNoFVHimalaya_SRC) $(MSSMNoFVHimalaya_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLMSSMNoFVHimalaya_MMA) $(MSSMNoFVHimalaya_INSTALL_DIR)
		$(Q)$(INSTALL_STRIPPED) $(MSSMNoFVHimalaya_MK) $(MSSMNoFVHimalaya_INSTALL_DIR) -m u=rw,g=r,o=r
		$(Q)install -m u=rw,g=r,o=r $(MSSMNoFVHimalaya_INCLUDE_MK) $(MSSMNoFVHimalaya_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(MSSMNoFVHimalaya_CXXQFT_VERTICES_MK) $(MSSMNoFVHimalaya_INSTALL_CXXQFT_DIR)

ifneq ($(MSSMNoFVHimalaya_SLHA_INPUT),)
		$(Q)install -m u=rw,g=r,o=r $(MSSMNoFVHimalaya_SLHA_INPUT) $(MSSMNoFVHimalaya_INSTALL_DIR)
endif
		$(Q)install -m u=rw,g=r,o=r $(MSSMNoFVHimalaya_REFERENCES) $(MSSMNoFVHimalaya_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(MSSMNoFVHimalaya_GNUPLOT) $(MSSMNoFVHimalaya_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		$(Q)-rm -f $(LIBMSSMNoFVHimalaya_DEP)
		$(Q)-rm -f $(EXEMSSMNoFVHimalaya_DEP)
		$(Q)-rm -f $(LLMSSMNoFVHimalaya_DEP)

clean-$(MODNAME)-lib:
		$(Q)-rm -f $(LIBMSSMNoFVHimalaya)
		$(Q)-rm -f $(LLMSSMNoFVHimalaya_LIB)

clean-$(MODNAME)-obj:
		$(Q)-rm -f $(LIBMSSMNoFVHimalaya_OBJ)
		$(Q)-rm -f $(EXEMSSMNoFVHimalaya_OBJ)
		$(Q)-rm -f $(LLMSSMNoFVHimalaya_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		$(Q)-rm -f $(LIBMSSMNoFVHimalaya_SRC)
		$(Q)-rm -f $(LIBMSSMNoFVHimalaya_HDR)
		$(Q)-rm -f $(LIBMSSMNoFVHimalaya_CXXQFT_HDR)
		$(Q)-rm -f $(EXEMSSMNoFVHimalaya_SRC)
		$(Q)-rm -f $(LLMSSMNoFVHimalaya_SRC)
		$(Q)-rm -f $(LLMSSMNoFVHimalaya_MMA)
		$(Q)-rm -f $(METACODE_STAMP_MSSMNoFVHimalaya)
		$(Q)-rm -f $(MSSMNoFVHimalaya_INCLUDE_MK)
		$(Q)-rm -f $(MSSMNoFVHimalaya_CXXQFT_VERTICES_MK)
		$(Q)-rm -f $(MSSMNoFVHimalaya_SLHA_INPUT)
		$(Q)-rm -f $(MSSMNoFVHimalaya_REFERENCES)
		$(Q)-rm -f $(MSSMNoFVHimalaya_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		$(Q)-rm -f $(EXEMSSMNoFVHimalaya_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		$(Q)tar -czf $(MSSMNoFVHimalaya_TARBALL) \
		$(LIBMSSMNoFVHimalaya_SRC) $(LIBMSSMNoFVHimalaya_HDR) $(LIBMSSMNoFVHimalaya_CXXQFT_HDR) \
		$(EXEMSSMNoFVHimalaya_SRC) \
		$(LLMSSMNoFVHimalaya_SRC) $(LLMSSMNoFVHimalaya_MMA) \
		$(MSSMNoFVHimalaya_MK) $(MSSMNoFVHimalaya_INCLUDE_MK) $(MSSMNoFVHimalaya_CXXQFT_VERTICES_MK) \
		$(MSSMNoFVHimalaya_SLHA_INPUT) $(MSSMNoFVHimalaya_REFERENCES) \
		$(MSSMNoFVHimalaya_GNUPLOT) \
		$(MSSMNoFVHimalaya_FlexibleDecay_MK)

$(LIBMSSMNoFVHimalaya_SRC) $(LIBMSSMNoFVHimalaya_HDR) $(LIBMSSMNoFVHimalaya_CXXQFT_HDR) $(EXEMSSMNoFVHimalaya_SRC) $(LLMSSMNoFVHimalaya_SRC) $(LLMSSMNoFVHimalaya_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_MSSMNoFVHimalaya)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_MSSMNoFVHimalaya): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_MSSMNoFVHimalaya)
		@$(MSG)
		$(Q)printf "%s" "Get[\"$<\"]; Quit[]" | "$(MATH)" || (echo "Error: The code generation failed!"; exit 1)
		@touch "$(METACODE_STAMP_MSSMNoFVHimalaya)"
		@echo "Note: to regenerate MSSMNoFVHimalaya source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_MSSMNoFVHimalaya)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_MSSMNoFVHimalaya):
		@true
endif

$(LIBMSSMNoFVHimalaya_DEP) $(EXEMSSMNoFVHimalaya_DEP) $(LLMSSMNoFVHimalaya_DEP) $(LIBMSSMNoFVHimalaya_OBJ) $(EXEMSSMNoFVHimalaya_OBJ) $(LLMSSMNoFVHimalaya_OBJ) $(LLMSSMNoFVHimalaya_LIB): \
	CPPFLAGS += $(MODMSSMNoFVHimalaya_SUBMOD_INC) $(MODMSSMNoFVHimalaya_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(GM2CALCFLAGS) $(HIGGSTOOLSFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBMSSMNoFVHimalaya_DEP) $(EXEMSSMNoFVHimalaya_DEP) $(LLMSSMNoFVHimalaya_DEP) $(LIBMSSMNoFVHimalaya_OBJ) $(EXEMSSMNoFVHimalaya_OBJ) $(LLMSSMNoFVHimalaya_OBJ) $(LLMSSMNoFVHimalaya_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLMSSMNoFVHimalaya_OBJ) $(LLMSSMNoFVHimalaya_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBMSSMNoFVHimalaya): $(LIBMSSMNoFVHimalaya_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBMSSMNoFVHimalaya) $(MODMSSMNoFVHimalaya_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(FUTILIBS)
		@$(MSG)
		$(Q)$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^) $(filter -%,$(LOOPFUNCLIBS)) $(GM2CALCLIBS) $(HIGGSTOOLSLIBS) $(PYTHONLIBS) $(HIMALAYALIBS) $(GSLLIBS) $(SQLITELIBS) $(TSILLIBS) $(FLIBS) $(THREADLIBS) $(LDLIBS) $(FUTILIBS)

$(LLMSSMNoFVHimalaya_LIB): $(LLMSSMNoFVHimalaya_OBJ) $(LIBMSSMNoFVHimalaya) $(MODMSSMNoFVHimalaya_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(FUTILIBS)
		@$(MSG)
		$(Q)$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^) $(filter -%,$(LOOPFUNCLIBS)) $(GM2CALCLIBS) $(HIGGSTOOLSLIBS) $(PYTHONLIBS) $(HIMALAYALIBS) $(TSILLIBS) $(GSLLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS) $(FUTILIBS) $(FLIBS)

ALLDEP += $(LIBMSSMNoFVHimalaya_DEP) $(EXEMSSMNoFVHimalaya_DEP)
ALLSRC += $(LIBMSSMNoFVHimalaya_SRC) $(EXEMSSMNoFVHimalaya_SRC)
ALLLIB += $(LIBMSSMNoFVHimalaya)
ALLEXE += $(EXEMSSMNoFVHimalaya_EXE)
ALLMODDEP += $(MODMSSMNoFVHimalaya_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLMSSMNoFVHimalaya_DEP)
ALLSRC += $(LLMSSMNoFVHimalaya_SRC)
ALLLL  += $(LLMSSMNoFVHimalaya_LIB)
endif
