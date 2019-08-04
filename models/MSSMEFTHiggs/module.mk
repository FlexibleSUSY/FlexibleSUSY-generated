DIR          := models/MSSMEFTHiggs
MODNAME      := MSSMEFTHiggs
SARAH_MODEL  := MSSM
WITH_$(MODNAME) := yes
MODMSSMEFTHiggs_MOD := SM MSSM_higgs
MODMSSMEFTHiggs_DEP := $(patsubst %,model_specific/%,$(MODMSSMEFTHiggs_MOD))
MODMSSMEFTHiggs_INC := $(patsubst %,-Imodel_specific/%,$(MODMSSMEFTHiggs_MOD))
MODMSSMEFTHiggs_LIB := $(foreach M,$(MODMSSMEFTHiggs_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

MODMSSMEFTHiggs_SUBMOD  := $(DIR)/cxx_qft
MODMSSMEFTHiggs_SUBMOD_INC := $(patsubst %,-I%,$(MODMSSMEFTHiggs_SUBMOD))

MSSMEFTHiggs_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)
MSSMEFTHiggs_INSTALL_CXXQFT_DIR := \
		$(MSSMEFTHiggs_INSTALL_DIR)/cxx_qft

MSSMEFTHiggs_MK     := \
		$(DIR)/module.mk

MSSMEFTHiggs_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

MSSMEFTHiggs_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

MSSMEFTHiggs_CXX_QFT_VERTICES_MK := \
		$(DIR)/cxx_qft/vertices.mk

MSSMEFTHiggs_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

MSSMEFTHiggs_INCLUDE_MK := \
		$(MSSMEFTHiggs_SUSY_BETAS_MK) \
		$(MSSMEFTHiggs_SOFT_BETAS_MK) \
		$(MSSMEFTHiggs_CXX_QFT_VERTICES_MK)

MSSMEFTHiggs_SLHA_INPUT := \
		$(DIR)/LesHouches.in.MSSMEFTHiggs_generated \
		$(DIR)/LesHouches.in.MSSMEFTHiggs

MSSMEFTHiggs_REFERENCES := \
		$(DIR)/MSSMEFTHiggs_references.tex

MSSMEFTHiggs_GNUPLOT := \
		$(DIR)/MSSMEFTHiggs_plot_rgflow.gnuplot \
		$(DIR)/MSSMEFTHiggs_plot_spectrum.gnuplot

MSSMEFTHiggs_TARBALL := \
		$(MODNAME).tar.gz

LIBMSSMEFTHiggs_SRC := \
		$(DIR)/MSSMEFTHiggs_a_muon.cpp \
		$(DIR)/MSSMEFTHiggs_edm.cpp \
		$(DIR)/MSSMEFTHiggs_FFV_form_factors.cpp \
		$(DIR)/MSSMEFTHiggs_l_to_lgamma.cpp \
		$(DIR)/MSSMEFTHiggs_effective_couplings.cpp \
		$(DIR)/MSSMEFTHiggs_info.cpp \
		$(DIR)/MSSMEFTHiggs_input_parameters.cpp \
		$(DIR)/MSSMEFTHiggs_mass_eigenstates.cpp \
		$(DIR)/MSSMEFTHiggs_observables.cpp \
		$(DIR)/MSSMEFTHiggs_physical.cpp \
		$(DIR)/MSSMEFTHiggs_slha_io.cpp \
		$(DIR)/MSSMEFTHiggs_soft_parameters.cpp \
		$(DIR)/MSSMEFTHiggs_susy_parameters.cpp \
		$(DIR)/MSSMEFTHiggs_utilities.cpp \
		$(DIR)/MSSMEFTHiggs_weinberg_angle.cpp

EXEMSSMEFTHiggs_SRC := \
		$(DIR)/run_MSSMEFTHiggs.cpp \
		$(DIR)/run_cmd_line_MSSMEFTHiggs.cpp \
		$(DIR)/scan_MSSMEFTHiggs.cpp
LLMSSMEFTHiggs_LIB  :=
LLMSSMEFTHiggs_OBJ  :=
LLMSSMEFTHiggs_SRC  := \
		$(DIR)/MSSMEFTHiggs_librarylink.cpp

LLMSSMEFTHiggs_MMA  := \
		$(DIR)/MSSMEFTHiggs_librarylink.m \
		$(DIR)/run_MSSMEFTHiggs.m

LIBMSSMEFTHiggs_HDR := \
		$(DIR)/MSSMEFTHiggs_a_muon.hpp \
		$(DIR)/MSSMEFTHiggs_convergence_tester.hpp \
		$(DIR)/MSSMEFTHiggs_edm.hpp \
		$(DIR)/MSSMEFTHiggs_FFV_form_factors.hpp \
		$(DIR)/MSSMEFTHiggs_l_to_lgamma.hpp \
		$(DIR)/MSSMEFTHiggs_effective_couplings.hpp \
		$(DIR)/MSSMEFTHiggs_ewsb_solver.hpp \
		$(DIR)/MSSMEFTHiggs_ewsb_solver_interface.hpp \
		$(DIR)/MSSMEFTHiggs_high_scale_constraint.hpp \
		$(DIR)/MSSMEFTHiggs_info.hpp \
		$(DIR)/MSSMEFTHiggs_initial_guesser.hpp \
		$(DIR)/MSSMEFTHiggs_input_parameters.hpp \
		$(DIR)/MSSMEFTHiggs_low_scale_constraint.hpp \
		$(DIR)/MSSMEFTHiggs_mass_eigenstates.hpp \
		$(DIR)/MSSMEFTHiggs_model.hpp \
		$(DIR)/MSSMEFTHiggs_model_slha.hpp \
		$(DIR)/MSSMEFTHiggs_observables.hpp \
		$(DIR)/MSSMEFTHiggs_physical.hpp \
		$(DIR)/MSSMEFTHiggs_slha_io.hpp \
		$(DIR)/MSSMEFTHiggs_spectrum_generator.hpp \
		$(DIR)/MSSMEFTHiggs_spectrum_generator_interface.hpp \
		$(DIR)/MSSMEFTHiggs_soft_parameters.hpp \
		$(DIR)/MSSMEFTHiggs_susy_parameters.hpp \
		$(DIR)/MSSMEFTHiggs_susy_scale_constraint.hpp \
		$(DIR)/MSSMEFTHiggs_utilities.hpp \
		$(DIR)/MSSMEFTHiggs_weinberg_angle.hpp

LIBMSSMEFTHiggs_CXXQFT_HDR := \
		$(DIR)/cxx_qft/MSSMEFTHiggs_qft.hpp \
		$(DIR)/cxx_qft/MSSMEFTHiggs_fields.hpp \
		$(DIR)/cxx_qft/MSSMEFTHiggs_vertices.hpp \
		$(DIR)/cxx_qft/MSSMEFTHiggs_context_base.hpp \
		$(DIR)/cxx_qft/MSSMEFTHiggs_npointfunctions.hpp

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
-include $(MSSMEFTHiggs_SUSY_BETAS_MK)
-include $(MSSMEFTHiggs_SOFT_BETAS_MK)
-include $(MSSMEFTHiggs_CXX_QFT_VERTICES_MK)
-include $(MSSMEFTHiggs_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(MSSMEFTHiggs_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(MSSMEFTHiggs_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(MSSMEFTHiggs_CXX_QFT_VERTICES_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(MSSMEFTHiggs_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBMSSMEFTHiggs_SRC := $(sort $(LIBMSSMEFTHiggs_SRC))
EXEMSSMEFTHiggs_SRC := $(sort $(EXEMSSMEFTHiggs_SRC))

LIBMSSMEFTHiggs_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBMSSMEFTHiggs_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBMSSMEFTHiggs_SRC)))

EXEMSSMEFTHiggs_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEMSSMEFTHiggs_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEMSSMEFTHiggs_SRC)))

EXEMSSMEFTHiggs_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXEMSSMEFTHiggs_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXEMSSMEFTHiggs_SRC)))

LIBMSSMEFTHiggs_DEP := \
		$(LIBMSSMEFTHiggs_OBJ:.o=.d)

EXEMSSMEFTHiggs_DEP := \
		$(EXEMSSMEFTHiggs_OBJ:.o=.d)

LLMSSMEFTHiggs_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLMSSMEFTHiggs_SRC)))

LLMSSMEFTHiggs_OBJ  := $(LLMSSMEFTHiggs_SRC:.cpp=.o)
LLMSSMEFTHiggs_LIB  := $(LLMSSMEFTHiggs_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBMSSMEFTHiggs     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_MSSMEFTHiggs := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_MSSMEFTHiggs := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL_DIR) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBMSSMEFTHiggs) $(EXEMSSMEFTHiggs_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(MSSMEFTHiggs_INSTALL_DIR)
		$(Q)install -d $(MSSMEFTHiggs_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBMSSMEFTHiggs_SRC) $(MSSMEFTHiggs_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBMSSMEFTHiggs_HDR) $(MSSMEFTHiggs_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBMSSMEFTHiggs_CXXQFT_HDR) $(MSSMEFTHiggs_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(EXEMSSMEFTHiggs_SRC) $(MSSMEFTHiggs_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLMSSMEFTHiggs_SRC) $(MSSMEFTHiggs_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLMSSMEFTHiggs_MMA) $(MSSMEFTHiggs_INSTALL_DIR)
		$(Q)$(INSTALL_STRIPPED) $(MSSMEFTHiggs_MK) $(MSSMEFTHiggs_INSTALL_DIR) -m u=rw,g=r,o=r
		$(Q)install -m u=rw,g=r,o=r $(MSSMEFTHiggs_INCLUDE_MK) $(MSSMEFTHiggs_INSTALL_DIR)
ifneq ($(MSSMEFTHiggs_SLHA_INPUT),)
		$(Q)install -m u=rw,g=r,o=r $(MSSMEFTHiggs_SLHA_INPUT) $(MSSMEFTHiggs_INSTALL_DIR)
endif
		$(Q)install -m u=rw,g=r,o=r $(MSSMEFTHiggs_REFERENCES) $(MSSMEFTHiggs_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(MSSMEFTHiggs_GNUPLOT) $(MSSMEFTHiggs_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		$(Q)-rm -f $(LIBMSSMEFTHiggs_DEP)
		$(Q)-rm -f $(EXEMSSMEFTHiggs_DEP)
		$(Q)-rm -f $(LLMSSMEFTHiggs_DEP)

clean-$(MODNAME)-lib:
		$(Q)-rm -f $(LIBMSSMEFTHiggs)
		$(Q)-rm -f $(LLMSSMEFTHiggs_LIB)

clean-$(MODNAME)-obj:
		$(Q)-rm -f $(LIBMSSMEFTHiggs_OBJ)
		$(Q)-rm -f $(EXEMSSMEFTHiggs_OBJ)
		$(Q)-rm -f $(LLMSSMEFTHiggs_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		$(Q)-rm -f $(LIBMSSMEFTHiggs_SRC)
		$(Q)-rm -f $(LIBMSSMEFTHiggs_HDR)
		$(Q)-rm -f $(LIBMSSMEFTHiggs_CXXQFT_HDR)
		$(Q)-rm -f $(EXEMSSMEFTHiggs_SRC)
		$(Q)-rm -f $(LLMSSMEFTHiggs_SRC)
		$(Q)-rm -f $(LLMSSMEFTHiggs_MMA)
		$(Q)-rm -f $(METACODE_STAMP_MSSMEFTHiggs)
		$(Q)-rm -f $(MSSMEFTHiggs_INCLUDE_MK)
		$(Q)-rm -f $(MSSMEFTHiggs_SLHA_INPUT)
		$(Q)-rm -f $(MSSMEFTHiggs_REFERENCES)
		$(Q)-rm -f $(MSSMEFTHiggs_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		$(Q)-rm -f $(EXEMSSMEFTHiggs_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		$(Q)tar -czf $(MSSMEFTHiggs_TARBALL) \
		$(LIBMSSMEFTHiggs_SRC) $(LIBMSSMEFTHiggs_HDR) $(LIBMSSMEFTHiggs_CXXQFT_HDR) \
		$(EXEMSSMEFTHiggs_SRC) \
		$(LLMSSMEFTHiggs_SRC) $(LLMSSMEFTHiggs_MMA) \
		$(MSSMEFTHiggs_MK) $(MSSMEFTHiggs_INCLUDE_MK) \
		$(MSSMEFTHiggs_SLHA_INPUT) $(MSSMEFTHiggs_REFERENCES) \
		$(MSSMEFTHiggs_GNUPLOT)

$(LIBMSSMEFTHiggs_SRC) $(LIBMSSMEFTHiggs_HDR) $(LIBMSSMEFTHiggs_CXXQFT_HDR) $(EXEMSSMEFTHiggs_SRC) $(LLMSSMEFTHiggs_SRC) $(LLMSSMEFTHiggs_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_MSSMEFTHiggs)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_MSSMEFTHiggs): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_MSSMEFTHiggs)
		@$(MSG)
		$(Q)"$(MATH)" -run "Get[\"$<\"]; Quit[]" || (echo "Error: The code generation failed!"; exit 1)
		@touch "$(METACODE_STAMP_MSSMEFTHiggs)"
		@echo "Note: to regenerate MSSMEFTHiggs source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_MSSMEFTHiggs)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_MSSMEFTHiggs):
		@true
endif

$(LIBMSSMEFTHiggs_DEP) $(EXEMSSMEFTHiggs_DEP) $(LLMSSMEFTHiggs_DEP) $(LIBMSSMEFTHiggs_OBJ) $(EXEMSSMEFTHiggs_OBJ) $(LLMSSMEFTHiggs_OBJ) $(LLMSSMEFTHiggs_LIB): \
	CPPFLAGS += $(MODMSSMEFTHiggs_SUBMOD_INC) $(MODMSSMEFTHiggs_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBMSSMEFTHiggs_DEP) $(EXEMSSMEFTHiggs_DEP) $(LLMSSMEFTHiggs_DEP) $(LIBMSSMEFTHiggs_OBJ) $(EXEMSSMEFTHiggs_OBJ) $(LLMSSMEFTHiggs_OBJ) $(LLMSSMEFTHiggs_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLMSSMEFTHiggs_OBJ) $(LLMSSMEFTHiggs_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBMSSMEFTHiggs): $(LIBMSSMEFTHiggs_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBMSSMEFTHiggs) $(MODMSSMEFTHiggs_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		@$(MSG)
		$(Q)$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLMSSMEFTHiggs_LIB): $(LLMSSMEFTHiggs_OBJ) $(LIBMSSMEFTHiggs) $(MODMSSMEFTHiggs_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		@$(MSG)
		$(Q)$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS)

ALLDEP += $(LIBMSSMEFTHiggs_DEP) $(EXEMSSMEFTHiggs_DEP)
ALLSRC += $(LIBMSSMEFTHiggs_SRC) $(EXEMSSMEFTHiggs_SRC)
ALLLIB += $(LIBMSSMEFTHiggs)
ALLEXE += $(EXEMSSMEFTHiggs_EXE)
ALLMODDEP += $(MODMSSMEFTHiggs_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLMSSMEFTHiggs_DEP)
ALLSRC += $(LLMSSMEFTHiggs_SRC)
ALLLL  += $(LLMSSMEFTHiggs_LIB)
endif
