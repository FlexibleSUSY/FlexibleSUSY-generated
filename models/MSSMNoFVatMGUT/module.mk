DIR          := models/MSSMNoFVatMGUT
MODNAME      := MSSMNoFVatMGUT
SARAH_MODEL  := MSSMNoFV
WITH_$(MODNAME) := yes
MODMSSMNoFVatMGUT_MOD := SM MSSM_higgs
MODMSSMNoFVatMGUT_DEP := $(patsubst %,model_specific/%,$(MODMSSMNoFVatMGUT_MOD))
MODMSSMNoFVatMGUT_INC := $(patsubst %,-Imodel_specific/%,$(MODMSSMNoFVatMGUT_MOD))
MODMSSMNoFVatMGUT_LIB := $(foreach M,$(MODMSSMNoFVatMGUT_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

MODMSSMNoFVatMGUT_SUBMOD  := $(DIR)/cxx_qft
MODMSSMNoFVatMGUT_SUBMOD_INC := $(patsubst %,-I%,$(MODMSSMNoFVatMGUT_SUBMOD))

MSSMNoFVatMGUT_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)
MSSMNoFVatMGUT_INSTALL_CXXQFT_DIR := \
		$(MSSMNoFVatMGUT_INSTALL_DIR)/cxx_qft

MSSMNoFVatMGUT_MK     := \
		$(DIR)/module.mk

MSSMNoFVatMGUT_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

MSSMNoFVatMGUT_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

MSSMNoFVatMGUT_CXXQFT_VERTICES_MK := \
		$(DIR)/cxx_qft/vertices.mk

-include $(MSSMNoFVatMGUT_CXXQFT_VERTICES_MK)
LIBMSSMNoFVatMGUT_CXXQFT_VERTICES_SRC ?= ''

MSSMNoFVatMGUT_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

MSSMNoFVatMGUT_INCLUDE_MK := \
		$(MSSMNoFVatMGUT_SUSY_BETAS_MK) \
		$(MSSMNoFVatMGUT_SOFT_BETAS_MK)

MSSMNoFVatMGUT_SLHA_INPUT := \
		$(DIR)/LesHouches.in.MSSMNoFVatMGUT_generated \
		$(DIR)/LesHouches.in.MSSMNoFVatMGUT

MSSMNoFVatMGUT_REFERENCES := \
		$(DIR)/MSSMNoFVatMGUT_references.tex

MSSMNoFVatMGUT_GNUPLOT := \
		$(DIR)/MSSMNoFVatMGUT_plot_rgflow.gnuplot \
		$(DIR)/MSSMNoFVatMGUT_plot_spectrum.gnuplot

MSSMNoFVatMGUT_TARBALL := \
		$(MODNAME).tar.gz

LIBMSSMNoFVatMGUT_SRC := \
		$(DIR)/MSSMNoFVatMGUT_a_muon.cpp \
		$(DIR)/MSSMNoFVatMGUT_edm.cpp \
		$(DIR)/MSSMNoFVatMGUT_FFV_form_factors.cpp \
		$(DIR)/MSSMNoFVatMGUT_f_to_f_conversion.cpp \
		$(DIR)/MSSMNoFVatMGUT_l_to_lgamma.cpp \
		$(DIR)/MSSMNoFVatMGUT_b_to_s_gamma.cpp \
		$(DIR)/MSSMNoFVatMGUT_effective_couplings.cpp \
		$(DIR)/MSSMNoFVatMGUT_info.cpp \
		$(DIR)/MSSMNoFVatMGUT_input_parameters.cpp \
		$(DIR)/MSSMNoFVatMGUT_mass_eigenstates.cpp \
		$(DIR)/MSSMNoFVatMGUT_mass_eigenstates_decoupling_scheme.cpp \
		$(DIR)/MSSMNoFVatMGUT_model_slha.cpp \
		$(DIR)/MSSMNoFVatMGUT_observables.cpp \
		$(DIR)/MSSMNoFVatMGUT_physical.cpp \
		$(DIR)/MSSMNoFVatMGUT_slha_io.cpp \
		$(DIR)/MSSMNoFVatMGUT_soft_parameters.cpp \
		$(DIR)/MSSMNoFVatMGUT_susy_parameters.cpp \
		$(DIR)/MSSMNoFVatMGUT_utilities.cpp \
		$(DIR)/MSSMNoFVatMGUT_weinberg_angle.cpp

LIBMSSMNoFVatMGUT_SRC += $(LIBMSSMNoFVatMGUT_CXXQFT_VERTICES_SRC)

EXEMSSMNoFVatMGUT_SRC := \
		$(DIR)/run_MSSMNoFVatMGUT.cpp \
		$(DIR)/run_cmd_line_MSSMNoFVatMGUT.cpp \
		$(DIR)/scan_MSSMNoFVatMGUT.cpp
LLMSSMNoFVatMGUT_LIB  :=
LLMSSMNoFVatMGUT_OBJ  :=
LLMSSMNoFVatMGUT_SRC  := \
		$(DIR)/MSSMNoFVatMGUT_librarylink.cpp

LLMSSMNoFVatMGUT_MMA  := \
		$(DIR)/MSSMNoFVatMGUT_librarylink.m \
		$(DIR)/run_MSSMNoFVatMGUT.m

LIBMSSMNoFVatMGUT_HDR := \
		$(DIR)/MSSMNoFVatMGUT_a_muon.hpp \
		$(DIR)/MSSMNoFVatMGUT_convergence_tester.hpp \
		$(DIR)/MSSMNoFVatMGUT_edm.hpp \
		$(DIR)/MSSMNoFVatMGUT_FFV_form_factors.hpp \
		$(DIR)/MSSMNoFVatMGUT_f_to_f_conversion.hpp \
		$(DIR)/MSSMNoFVatMGUT_l_to_lgamma.hpp \
		$(DIR)/MSSMNoFVatMGUT_b_to_s_gamma.hpp \
		$(DIR)/MSSMNoFVatMGUT_effective_couplings.hpp \
		$(DIR)/MSSMNoFVatMGUT_ewsb_solver.hpp \
		$(DIR)/MSSMNoFVatMGUT_ewsb_solver_interface.hpp \
		$(DIR)/MSSMNoFVatMGUT_high_scale_constraint.hpp \
		$(DIR)/MSSMNoFVatMGUT_info.hpp \
		$(DIR)/MSSMNoFVatMGUT_initial_guesser.hpp \
		$(DIR)/MSSMNoFVatMGUT_input_parameters.hpp \
		$(DIR)/MSSMNoFVatMGUT_low_scale_constraint.hpp \
		$(DIR)/MSSMNoFVatMGUT_mass_eigenstates.hpp \
		$(DIR)/MSSMNoFVatMGUT_mass_eigenstates_interface.hpp \
		$(DIR)/MSSMNoFVatMGUT_mass_eigenstates_decoupling_scheme.hpp \
		$(DIR)/MSSMNoFVatMGUT_model.hpp \
		$(DIR)/MSSMNoFVatMGUT_model_slha.hpp \
		$(DIR)/MSSMNoFVatMGUT_observables.hpp \
		$(DIR)/MSSMNoFVatMGUT_physical.hpp \
		$(DIR)/MSSMNoFVatMGUT_slha_io.hpp \
		$(DIR)/MSSMNoFVatMGUT_spectrum_generator.hpp \
		$(DIR)/MSSMNoFVatMGUT_spectrum_generator_interface.hpp \
		$(DIR)/MSSMNoFVatMGUT_soft_parameters.hpp \
		$(DIR)/MSSMNoFVatMGUT_susy_parameters.hpp \
		$(DIR)/MSSMNoFVatMGUT_susy_scale_constraint.hpp \
		$(DIR)/MSSMNoFVatMGUT_utilities.hpp \
		$(DIR)/MSSMNoFVatMGUT_weinberg_angle.hpp

LIBMSSMNoFVatMGUT_CXXQFT_HDR := \
		$(DIR)/cxx_qft/MSSMNoFVatMGUT_qft.hpp \
		$(DIR)/cxx_qft/MSSMNoFVatMGUT_fields.hpp \
		$(DIR)/cxx_qft/MSSMNoFVatMGUT_vertices.hpp \
		$(DIR)/cxx_qft/MSSMNoFVatMGUT_context_base.hpp \
		$(DIR)/cxx_qft/MSSMNoFVatMGUT_npointfunctions_wilsoncoeffs.hpp

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
-include $(MSSMNoFVatMGUT_SUSY_BETAS_MK)
-include $(MSSMNoFVatMGUT_SOFT_BETAS_MK)
-include $(MSSMNoFVatMGUT_CXXQFT_VERTICES_MK)
-include $(MSSMNoFVatMGUT_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(MSSMNoFVatMGUT_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(MSSMNoFVatMGUT_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(MSSMNoFVatMGUT_CXXQFT_VERTICES_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(MSSMNoFVatMGUT_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBMSSMNoFVatMGUT_SRC := $(sort $(LIBMSSMNoFVatMGUT_SRC))
EXEMSSMNoFVatMGUT_SRC := $(sort $(EXEMSSMNoFVatMGUT_SRC))

LIBMSSMNoFVatMGUT_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBMSSMNoFVatMGUT_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBMSSMNoFVatMGUT_SRC)))

EXEMSSMNoFVatMGUT_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEMSSMNoFVatMGUT_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEMSSMNoFVatMGUT_SRC)))

EXEMSSMNoFVatMGUT_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXEMSSMNoFVatMGUT_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXEMSSMNoFVatMGUT_SRC)))

LIBMSSMNoFVatMGUT_DEP := \
		$(LIBMSSMNoFVatMGUT_OBJ:.o=.d)

EXEMSSMNoFVatMGUT_DEP := \
		$(EXEMSSMNoFVatMGUT_OBJ:.o=.d)

LLMSSMNoFVatMGUT_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLMSSMNoFVatMGUT_SRC)))

LLMSSMNoFVatMGUT_OBJ  := $(LLMSSMNoFVatMGUT_SRC:.cpp=.o)
LLMSSMNoFVatMGUT_LIB  := $(LLMSSMNoFVatMGUT_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBMSSMNoFVatMGUT     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_MSSMNoFVatMGUT := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_MSSMNoFVatMGUT := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL_DIR) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBMSSMNoFVatMGUT) $(EXEMSSMNoFVatMGUT_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(MSSMNoFVatMGUT_INSTALL_DIR)
		$(Q)install -d $(MSSMNoFVatMGUT_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBMSSMNoFVatMGUT_SRC) $(MSSMNoFVatMGUT_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBMSSMNoFVatMGUT_CXXQFT_VERTICES_SRC) $(MSSMNoFVatMGUT_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBMSSMNoFVatMGUT_HDR) $(MSSMNoFVatMGUT_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBMSSMNoFVatMGUT_CXXQFT_HDR) $(MSSMNoFVatMGUT_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(EXEMSSMNoFVatMGUT_SRC) $(MSSMNoFVatMGUT_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLMSSMNoFVatMGUT_SRC) $(MSSMNoFVatMGUT_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLMSSMNoFVatMGUT_MMA) $(MSSMNoFVatMGUT_INSTALL_DIR)
		$(Q)$(INSTALL_STRIPPED) $(MSSMNoFVatMGUT_MK) $(MSSMNoFVatMGUT_INSTALL_DIR) -m u=rw,g=r,o=r
		$(Q)install -m u=rw,g=r,o=r $(MSSMNoFVatMGUT_INCLUDE_MK) $(MSSMNoFVatMGUT_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(MSSMNoFVatMGUT_CXXQFT_VERTICES_MK) $(MSSMNoFVatMGUT_INSTALL_CXXQFT_DIR)

ifneq ($(MSSMNoFVatMGUT_SLHA_INPUT),)
		$(Q)install -m u=rw,g=r,o=r $(MSSMNoFVatMGUT_SLHA_INPUT) $(MSSMNoFVatMGUT_INSTALL_DIR)
endif
		$(Q)install -m u=rw,g=r,o=r $(MSSMNoFVatMGUT_REFERENCES) $(MSSMNoFVatMGUT_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(MSSMNoFVatMGUT_GNUPLOT) $(MSSMNoFVatMGUT_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		$(Q)-rm -f $(LIBMSSMNoFVatMGUT_DEP)
		$(Q)-rm -f $(EXEMSSMNoFVatMGUT_DEP)
		$(Q)-rm -f $(LLMSSMNoFVatMGUT_DEP)

clean-$(MODNAME)-lib:
		$(Q)-rm -f $(LIBMSSMNoFVatMGUT)
		$(Q)-rm -f $(LLMSSMNoFVatMGUT_LIB)

clean-$(MODNAME)-obj:
		$(Q)-rm -f $(LIBMSSMNoFVatMGUT_OBJ)
		$(Q)-rm -f $(EXEMSSMNoFVatMGUT_OBJ)
		$(Q)-rm -f $(LLMSSMNoFVatMGUT_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		$(Q)-rm -f $(LIBMSSMNoFVatMGUT_SRC)
		$(Q)-rm -f $(LIBMSSMNoFVatMGUT_HDR)
		$(Q)-rm -f $(LIBMSSMNoFVatMGUT_CXXQFT_HDR)
		$(Q)-rm -f $(EXEMSSMNoFVatMGUT_SRC)
		$(Q)-rm -f $(LLMSSMNoFVatMGUT_SRC)
		$(Q)-rm -f $(LLMSSMNoFVatMGUT_MMA)
		$(Q)-rm -f $(METACODE_STAMP_MSSMNoFVatMGUT)
		$(Q)-rm -f $(MSSMNoFVatMGUT_INCLUDE_MK)
		$(Q)-rm -f $(MSSMNoFVatMGUT_CXXQFT_VERTICES_MK)
		$(Q)-rm -f $(MSSMNoFVatMGUT_SLHA_INPUT)
		$(Q)-rm -f $(MSSMNoFVatMGUT_REFERENCES)
		$(Q)-rm -f $(MSSMNoFVatMGUT_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		$(Q)-rm -f $(EXEMSSMNoFVatMGUT_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		$(Q)tar -czf $(MSSMNoFVatMGUT_TARBALL) \
		$(LIBMSSMNoFVatMGUT_SRC) $(LIBMSSMNoFVatMGUT_HDR) $(LIBMSSMNoFVatMGUT_CXXQFT_HDR) \
		$(EXEMSSMNoFVatMGUT_SRC) \
		$(LLMSSMNoFVatMGUT_SRC) $(LLMSSMNoFVatMGUT_MMA) \
		$(MSSMNoFVatMGUT_MK) $(MSSMNoFVatMGUT_INCLUDE_MK) $(MSSMNoFVatMGUT_CXXQFT_VERTICES_MK) \
		$(MSSMNoFVatMGUT_SLHA_INPUT) $(MSSMNoFVatMGUT_REFERENCES) \
		$(MSSMNoFVatMGUT_GNUPLOT)

$(LIBMSSMNoFVatMGUT_SRC) $(LIBMSSMNoFVatMGUT_HDR) $(LIBMSSMNoFVatMGUT_CXXQFT_HDR) $(EXEMSSMNoFVatMGUT_SRC) $(LLMSSMNoFVatMGUT_SRC) $(LLMSSMNoFVatMGUT_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_MSSMNoFVatMGUT)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_MSSMNoFVatMGUT): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_MSSMNoFVatMGUT)
		@$(MSG)
		$(Q)printf "%s" "Get[\"$<\"]; Quit[]" | "$(MATH)" || (echo "Error: The code generation failed!"; exit 1)
		@touch "$(METACODE_STAMP_MSSMNoFVatMGUT)"
		@echo "Note: to regenerate MSSMNoFVatMGUT source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_MSSMNoFVatMGUT)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_MSSMNoFVatMGUT):
		@true
endif

$(LIBMSSMNoFVatMGUT_DEP) $(EXEMSSMNoFVatMGUT_DEP) $(LLMSSMNoFVatMGUT_DEP) $(LIBMSSMNoFVatMGUT_OBJ) $(EXEMSSMNoFVatMGUT_OBJ) $(LLMSSMNoFVatMGUT_OBJ) $(LLMSSMNoFVatMGUT_LIB): \
	CPPFLAGS += $(MODMSSMNoFVatMGUT_SUBMOD_INC) $(MODMSSMNoFVatMGUT_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(GM2CALCFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBMSSMNoFVatMGUT_DEP) $(EXEMSSMNoFVatMGUT_DEP) $(LLMSSMNoFVatMGUT_DEP) $(LIBMSSMNoFVatMGUT_OBJ) $(EXEMSSMNoFVatMGUT_OBJ) $(LLMSSMNoFVatMGUT_OBJ) $(LLMSSMNoFVatMGUT_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLMSSMNoFVatMGUT_OBJ) $(LLMSSMNoFVatMGUT_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBMSSMNoFVatMGUT): $(LIBMSSMNoFVatMGUT_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBMSSMNoFVatMGUT) $(MODMSSMNoFVatMGUT_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(FUTILIBS)
		@$(MSG)
		$(Q)$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^) $(filter -%,$(LOOPFUNCLIBS)) $(GM2CALCLIBS) $(HIMALAYALIBS) $(GSLLIBS) $(SQLITELIBS) $(TSILLIBS) $(FLIBS) $(THREADLIBS) $(LDLIBS) $(FUTILIBS)

$(LLMSSMNoFVatMGUT_LIB): $(LLMSSMNoFVatMGUT_OBJ) $(LIBMSSMNoFVatMGUT) $(MODMSSMNoFVatMGUT_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(FUTILIBS)
		@$(MSG)
		$(Q)$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^) $(filter -%,$(LOOPFUNCLIBS)) $(GM2CALCLIBS) $(HIMALAYALIBS) $(TSILLIBS) $(GSLLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS) $(FUTILIBS) $(FLIBS)

ALLDEP += $(LIBMSSMNoFVatMGUT_DEP) $(EXEMSSMNoFVatMGUT_DEP)
ALLSRC += $(LIBMSSMNoFVatMGUT_SRC) $(EXEMSSMNoFVatMGUT_SRC)
ALLLIB += $(LIBMSSMNoFVatMGUT)
ALLEXE += $(EXEMSSMNoFVatMGUT_EXE)
ALLMODDEP += $(MODMSSMNoFVatMGUT_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLMSSMNoFVatMGUT_DEP)
ALLSRC += $(LLMSSMNoFVatMGUT_SRC)
ALLLL  += $(LLMSSMNoFVatMGUT_LIB)
endif
