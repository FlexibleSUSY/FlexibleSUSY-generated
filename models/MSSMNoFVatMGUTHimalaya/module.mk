DIR          := models/MSSMNoFVatMGUTHimalaya
MODNAME      := MSSMNoFVatMGUTHimalaya
SARAH_MODEL  := MSSMNoFV
WITH_$(MODNAME) := yes
MODMSSMNoFVatMGUTHimalaya_MOD := SM MSSM_higgs MSSM_thresholds
MODMSSMNoFVatMGUTHimalaya_DEP := $(patsubst %,model_specific/%,$(MODMSSMNoFVatMGUTHimalaya_MOD))
MODMSSMNoFVatMGUTHimalaya_INC := $(patsubst %,-Imodel_specific/%,$(MODMSSMNoFVatMGUTHimalaya_MOD))
MODMSSMNoFVatMGUTHimalaya_LIB := $(foreach M,$(MODMSSMNoFVatMGUTHimalaya_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

MODMSSMNoFVatMGUTHimalaya_SUBMOD  := $(DIR)/cxx_qft
MODMSSMNoFVatMGUTHimalaya_SUBMOD_INC := $(patsubst %,-I%,$(MODMSSMNoFVatMGUTHimalaya_SUBMOD))

MSSMNoFVatMGUTHimalaya_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)
MSSMNoFVatMGUTHimalaya_INSTALL_CXXQFT_DIR := \
		$(MSSMNoFVatMGUTHimalaya_INSTALL_DIR)/cxx_qft

MSSMNoFVatMGUTHimalaya_MK     := \
		$(DIR)/module.mk

MSSMNoFVatMGUTHimalaya_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

MSSMNoFVatMGUTHimalaya_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

MSSMNoFVatMGUTHimalaya_CXXQFT_VERTICES_MK := \
		$(DIR)/cxx_qft/vertices.mk

-include $(MSSMNoFVatMGUTHimalaya_CXXQFT_VERTICES_MK)
LIBMSSMNoFVatMGUTHimalaya_CXXQFT_VERTICES_SRC ?= ''

MSSMNoFVatMGUTHimalaya_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

MSSMNoFVatMGUTHimalaya_FlexibleDecay_MK := \
		$(DIR)/decays/FlexibleDecay.mk

MSSMNoFVatMGUTHimalaya_INCLUDE_MK := \
		$(MSSMNoFVatMGUTHimalaya_SUSY_BETAS_MK) \
		$(MSSMNoFVatMGUTHimalaya_SOFT_BETAS_MK)

MSSMNoFVatMGUTHimalaya_SLHA_INPUT := \
		$(DIR)/LesHouches.in.MSSMNoFVatMGUTHimalaya_generated \
		$(DIR)/LesHouches.in.MSSMNoFVatMGUTHimalaya

MSSMNoFVatMGUTHimalaya_REFERENCES := \
		$(DIR)/MSSMNoFVatMGUTHimalaya_references.tex

MSSMNoFVatMGUTHimalaya_GNUPLOT := \
		$(DIR)/MSSMNoFVatMGUTHimalaya_plot_rgflow.gnuplot \
		$(DIR)/MSSMNoFVatMGUTHimalaya_plot_spectrum.gnuplot

MSSMNoFVatMGUTHimalaya_TARBALL := \
		$(MODNAME).tar.gz

LIBMSSMNoFVatMGUTHimalaya_SRC := \
		$(DIR)/MSSMNoFVatMGUTHimalaya_a_muon.cpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_edm.cpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_FFV_form_factors.cpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_f_to_f_conversion.cpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_l_to_lgamma.cpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_b_to_s_gamma.cpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_info.cpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_input_parameters.cpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_mass_eigenstates.cpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_mass_eigenstates_decoupling_scheme.cpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_model_slha.cpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_observables.cpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_physical.cpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_slha_io.cpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_soft_parameters.cpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_susy_parameters.cpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_utilities.cpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_weinberg_angle.cpp

LIBMSSMNoFVatMGUTHimalaya_SRC += $(LIBMSSMNoFVatMGUTHimalaya_CXXQFT_VERTICES_SRC)

EXEMSSMNoFVatMGUTHimalaya_SRC := \
		$(DIR)/run_MSSMNoFVatMGUTHimalaya.cpp \
		$(DIR)/run_cmd_line_MSSMNoFVatMGUTHimalaya.cpp \
		$(DIR)/scan_MSSMNoFVatMGUTHimalaya.cpp
LLMSSMNoFVatMGUTHimalaya_LIB  :=
LLMSSMNoFVatMGUTHimalaya_OBJ  :=
LLMSSMNoFVatMGUTHimalaya_SRC  := \
		$(DIR)/MSSMNoFVatMGUTHimalaya_librarylink.cpp

LLMSSMNoFVatMGUTHimalaya_MMA  := \
		$(DIR)/MSSMNoFVatMGUTHimalaya_librarylink.m \
		$(DIR)/run_MSSMNoFVatMGUTHimalaya.m

LIBMSSMNoFVatMGUTHimalaya_HDR := \
		$(DIR)/MSSMNoFVatMGUTHimalaya_a_muon.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_convergence_tester.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_edm.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_FFV_form_factors.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_f_to_f_conversion.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_l_to_lgamma.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_b_to_s_gamma.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_ewsb_solver.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_ewsb_solver_interface.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_high_scale_constraint.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_info.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_initial_guesser.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_input_parameters.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_low_scale_constraint.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_mass_eigenstates.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_mass_eigenstates_interface.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_mass_eigenstates_decoupling_scheme.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_model.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_model_slha.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_observables.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_physical.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_slha_io.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_spectrum_generator.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_spectrum_generator_interface.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_soft_parameters.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_susy_parameters.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_susy_scale_constraint.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_utilities.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_weinberg_angle.hpp

LIBMSSMNoFVatMGUTHimalaya_CXXQFT_HDR := \
		$(DIR)/cxx_qft/MSSMNoFVatMGUTHimalaya_qft.hpp \
		$(DIR)/cxx_qft/MSSMNoFVatMGUTHimalaya_fields.hpp \
		$(DIR)/cxx_qft/MSSMNoFVatMGUTHimalaya_vertices.hpp \
		$(DIR)/cxx_qft/MSSMNoFVatMGUTHimalaya_context_base.hpp \
		$(DIR)/cxx_qft/MSSMNoFVatMGUTHimalaya_npointfunctions_wilsoncoeffs.hpp

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
-include $(MSSMNoFVatMGUTHimalaya_SUSY_BETAS_MK)
-include $(MSSMNoFVatMGUTHimalaya_SOFT_BETAS_MK)
-include $(MSSMNoFVatMGUTHimalaya_FlexibleDecay_MK)
-include $(MSSMNoFVatMGUTHimalaya_CXXQFT_VERTICES_MK)
-include $(MSSMNoFVatMGUTHimalaya_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(MSSMNoFVatMGUTHimalaya_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(MSSMNoFVatMGUTHimalaya_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@

$(MSSMNoFVatMGUTHimalaya_FlexibleDecay_MK): run-metacode-$(MODNAME)
$(MSSMNoFVatMGUTHimalaya_CXXQFT_VERTICES_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(MSSMNoFVatMGUTHimalaya_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBMSSMNoFVatMGUTHimalaya_SRC := $(sort $(LIBMSSMNoFVatMGUTHimalaya_SRC))
EXEMSSMNoFVatMGUTHimalaya_SRC := $(sort $(EXEMSSMNoFVatMGUTHimalaya_SRC))

LIBMSSMNoFVatMGUTHimalaya_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBMSSMNoFVatMGUTHimalaya_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBMSSMNoFVatMGUTHimalaya_SRC)))

EXEMSSMNoFVatMGUTHimalaya_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEMSSMNoFVatMGUTHimalaya_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEMSSMNoFVatMGUTHimalaya_SRC)))

EXEMSSMNoFVatMGUTHimalaya_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXEMSSMNoFVatMGUTHimalaya_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXEMSSMNoFVatMGUTHimalaya_SRC)))

LIBMSSMNoFVatMGUTHimalaya_DEP := \
		$(LIBMSSMNoFVatMGUTHimalaya_OBJ:.o=.d)

EXEMSSMNoFVatMGUTHimalaya_DEP := \
		$(EXEMSSMNoFVatMGUTHimalaya_OBJ:.o=.d)

LLMSSMNoFVatMGUTHimalaya_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLMSSMNoFVatMGUTHimalaya_SRC)))

LLMSSMNoFVatMGUTHimalaya_OBJ  := $(LLMSSMNoFVatMGUTHimalaya_SRC:.cpp=.o)
LLMSSMNoFVatMGUTHimalaya_LIB  := $(LLMSSMNoFVatMGUTHimalaya_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBMSSMNoFVatMGUTHimalaya     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_MSSMNoFVatMGUTHimalaya := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_MSSMNoFVatMGUTHimalaya := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL_DIR) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBMSSMNoFVatMGUTHimalaya) $(EXEMSSMNoFVatMGUTHimalaya_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(MSSMNoFVatMGUTHimalaya_INSTALL_DIR)
		$(Q)install -d $(MSSMNoFVatMGUTHimalaya_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBMSSMNoFVatMGUTHimalaya_SRC) $(MSSMNoFVatMGUTHimalaya_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBMSSMNoFVatMGUTHimalaya_CXXQFT_VERTICES_SRC) $(MSSMNoFVatMGUTHimalaya_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBMSSMNoFVatMGUTHimalaya_HDR) $(MSSMNoFVatMGUTHimalaya_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBMSSMNoFVatMGUTHimalaya_CXXQFT_HDR) $(MSSMNoFVatMGUTHimalaya_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(EXEMSSMNoFVatMGUTHimalaya_SRC) $(MSSMNoFVatMGUTHimalaya_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLMSSMNoFVatMGUTHimalaya_SRC) $(MSSMNoFVatMGUTHimalaya_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLMSSMNoFVatMGUTHimalaya_MMA) $(MSSMNoFVatMGUTHimalaya_INSTALL_DIR)
		$(Q)$(INSTALL_STRIPPED) $(MSSMNoFVatMGUTHimalaya_MK) $(MSSMNoFVatMGUTHimalaya_INSTALL_DIR) -m u=rw,g=r,o=r
		$(Q)install -m u=rw,g=r,o=r $(MSSMNoFVatMGUTHimalaya_INCLUDE_MK) $(MSSMNoFVatMGUTHimalaya_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(MSSMNoFVatMGUTHimalaya_CXXQFT_VERTICES_MK) $(MSSMNoFVatMGUTHimalaya_INSTALL_CXXQFT_DIR)

ifneq ($(MSSMNoFVatMGUTHimalaya_SLHA_INPUT),)
		$(Q)install -m u=rw,g=r,o=r $(MSSMNoFVatMGUTHimalaya_SLHA_INPUT) $(MSSMNoFVatMGUTHimalaya_INSTALL_DIR)
endif
		$(Q)install -m u=rw,g=r,o=r $(MSSMNoFVatMGUTHimalaya_REFERENCES) $(MSSMNoFVatMGUTHimalaya_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(MSSMNoFVatMGUTHimalaya_GNUPLOT) $(MSSMNoFVatMGUTHimalaya_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		$(Q)-rm -f $(LIBMSSMNoFVatMGUTHimalaya_DEP)
		$(Q)-rm -f $(EXEMSSMNoFVatMGUTHimalaya_DEP)
		$(Q)-rm -f $(LLMSSMNoFVatMGUTHimalaya_DEP)

clean-$(MODNAME)-lib:
		$(Q)-rm -f $(LIBMSSMNoFVatMGUTHimalaya)
		$(Q)-rm -f $(LLMSSMNoFVatMGUTHimalaya_LIB)

clean-$(MODNAME)-obj:
		$(Q)-rm -f $(LIBMSSMNoFVatMGUTHimalaya_OBJ)
		$(Q)-rm -f $(EXEMSSMNoFVatMGUTHimalaya_OBJ)
		$(Q)-rm -f $(LLMSSMNoFVatMGUTHimalaya_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		$(Q)-rm -f $(LIBMSSMNoFVatMGUTHimalaya_SRC)
		$(Q)-rm -f $(LIBMSSMNoFVatMGUTHimalaya_HDR)
		$(Q)-rm -f $(LIBMSSMNoFVatMGUTHimalaya_CXXQFT_HDR)
		$(Q)-rm -f $(EXEMSSMNoFVatMGUTHimalaya_SRC)
		$(Q)-rm -f $(LLMSSMNoFVatMGUTHimalaya_SRC)
		$(Q)-rm -f $(LLMSSMNoFVatMGUTHimalaya_MMA)
		$(Q)-rm -f $(METACODE_STAMP_MSSMNoFVatMGUTHimalaya)
		$(Q)-rm -f $(MSSMNoFVatMGUTHimalaya_INCLUDE_MK)
		$(Q)-rm -f $(MSSMNoFVatMGUTHimalaya_CXXQFT_VERTICES_MK)
		$(Q)-rm -f $(MSSMNoFVatMGUTHimalaya_SLHA_INPUT)
		$(Q)-rm -f $(MSSMNoFVatMGUTHimalaya_REFERENCES)
		$(Q)-rm -f $(MSSMNoFVatMGUTHimalaya_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		$(Q)-rm -f $(EXEMSSMNoFVatMGUTHimalaya_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		$(Q)tar -czf $(MSSMNoFVatMGUTHimalaya_TARBALL) \
		$(LIBMSSMNoFVatMGUTHimalaya_SRC) $(LIBMSSMNoFVatMGUTHimalaya_HDR) $(LIBMSSMNoFVatMGUTHimalaya_CXXQFT_HDR) \
		$(EXEMSSMNoFVatMGUTHimalaya_SRC) \
		$(LLMSSMNoFVatMGUTHimalaya_SRC) $(LLMSSMNoFVatMGUTHimalaya_MMA) \
		$(MSSMNoFVatMGUTHimalaya_MK) $(MSSMNoFVatMGUTHimalaya_INCLUDE_MK) $(MSSMNoFVatMGUTHimalaya_CXXQFT_VERTICES_MK) \
		$(MSSMNoFVatMGUTHimalaya_SLHA_INPUT) $(MSSMNoFVatMGUTHimalaya_REFERENCES) \
		$(MSSMNoFVatMGUTHimalaya_GNUPLOT) \
		$(MSSMNoFVatMGUTHimalaya_FlexibleDecay_MK)

$(LIBMSSMNoFVatMGUTHimalaya_SRC) $(LIBMSSMNoFVatMGUTHimalaya_HDR) $(LIBMSSMNoFVatMGUTHimalaya_CXXQFT_HDR) $(EXEMSSMNoFVatMGUTHimalaya_SRC) $(LLMSSMNoFVatMGUTHimalaya_SRC) $(LLMSSMNoFVatMGUTHimalaya_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_MSSMNoFVatMGUTHimalaya)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_MSSMNoFVatMGUTHimalaya): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_MSSMNoFVatMGUTHimalaya)
		@$(MSG)
		$(Q)printf "%s" "Get[\"$<\"]; Quit[]" | "$(MATH)" || (echo "Error: The code generation failed!"; exit 1)
		@touch "$(METACODE_STAMP_MSSMNoFVatMGUTHimalaya)"
		@echo "Note: to regenerate MSSMNoFVatMGUTHimalaya source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_MSSMNoFVatMGUTHimalaya)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_MSSMNoFVatMGUTHimalaya):
		@true
endif

$(LIBMSSMNoFVatMGUTHimalaya_DEP) $(EXEMSSMNoFVatMGUTHimalaya_DEP) $(LLMSSMNoFVatMGUTHimalaya_DEP) $(LIBMSSMNoFVatMGUTHimalaya_OBJ) $(EXEMSSMNoFVatMGUTHimalaya_OBJ) $(LLMSSMNoFVatMGUTHimalaya_OBJ) $(LLMSSMNoFVatMGUTHimalaya_LIB): \
	CPPFLAGS += $(MODMSSMNoFVatMGUTHimalaya_SUBMOD_INC) $(MODMSSMNoFVatMGUTHimalaya_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(GM2CALCFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBMSSMNoFVatMGUTHimalaya_DEP) $(EXEMSSMNoFVatMGUTHimalaya_DEP) $(LLMSSMNoFVatMGUTHimalaya_DEP) $(LIBMSSMNoFVatMGUTHimalaya_OBJ) $(EXEMSSMNoFVatMGUTHimalaya_OBJ) $(LLMSSMNoFVatMGUTHimalaya_OBJ) $(LLMSSMNoFVatMGUTHimalaya_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLMSSMNoFVatMGUTHimalaya_OBJ) $(LLMSSMNoFVatMGUTHimalaya_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBMSSMNoFVatMGUTHimalaya): $(LIBMSSMNoFVatMGUTHimalaya_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBMSSMNoFVatMGUTHimalaya) $(MODMSSMNoFVatMGUTHimalaya_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(FUTILIBS)
		@$(MSG)
		$(Q)$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^) $(filter -%,$(LOOPFUNCLIBS)) $(GM2CALCLIBS) $(HIMALAYALIBS) $(GSLLIBS) $(SQLITELIBS) $(TSILLIBS) $(FLIBS) $(THREADLIBS) $(LDLIBS) $(FUTILIBS)

$(LLMSSMNoFVatMGUTHimalaya_LIB): $(LLMSSMNoFVatMGUTHimalaya_OBJ) $(LIBMSSMNoFVatMGUTHimalaya) $(MODMSSMNoFVatMGUTHimalaya_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(FUTILIBS)
		@$(MSG)
		$(Q)$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^) $(filter -%,$(LOOPFUNCLIBS)) $(GM2CALCLIBS) $(HIMALAYALIBS) $(TSILLIBS) $(GSLLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS) $(FUTILIBS) $(FLIBS)

ALLDEP += $(LIBMSSMNoFVatMGUTHimalaya_DEP) $(EXEMSSMNoFVatMGUTHimalaya_DEP)
ALLSRC += $(LIBMSSMNoFVatMGUTHimalaya_SRC) $(EXEMSSMNoFVatMGUTHimalaya_SRC)
ALLLIB += $(LIBMSSMNoFVatMGUTHimalaya)
ALLEXE += $(EXEMSSMNoFVatMGUTHimalaya_EXE)
ALLMODDEP += $(MODMSSMNoFVatMGUTHimalaya_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLMSSMNoFVatMGUTHimalaya_DEP)
ALLSRC += $(LLMSSMNoFVatMGUTHimalaya_SRC)
ALLLL  += $(LLMSSMNoFVatMGUTHimalaya_LIB)
endif
