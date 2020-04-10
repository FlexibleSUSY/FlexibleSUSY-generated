DIR          := models/MSSMatMGUT
MODNAME      := MSSMatMGUT
SARAH_MODEL  := MSSM
WITH_$(MODNAME) := yes
MODMSSMatMGUT_MOD := SM MSSM_higgs
MODMSSMatMGUT_DEP := $(patsubst %,model_specific/%,$(MODMSSMatMGUT_MOD))
MODMSSMatMGUT_INC := $(patsubst %,-Imodel_specific/%,$(MODMSSMatMGUT_MOD))
MODMSSMatMGUT_LIB := $(foreach M,$(MODMSSMatMGUT_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

MODMSSMatMGUT_SUBMOD  := $(DIR)/cxx_qft
MODMSSMatMGUT_SUBMOD_INC := $(patsubst %,-I%,$(MODMSSMatMGUT_SUBMOD))

MSSMatMGUT_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)
MSSMatMGUT_INSTALL_CXXQFT_DIR := \
		$(MSSMatMGUT_INSTALL_DIR)/cxx_qft

MSSMatMGUT_MK     := \
		$(DIR)/module.mk

MSSMatMGUT_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

MSSMatMGUT_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

MSSMatMGUT_CXX_QFT_VERTICES_MK := \
		$(DIR)/cxx_qft/vertices.mk

MSSMatMGUT_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

MSSMatMGUT_INCLUDE_MK := \
		$(MSSMatMGUT_SUSY_BETAS_MK) \
		$(MSSMatMGUT_SOFT_BETAS_MK) \
		$(MSSMatMGUT_CXX_QFT_VERTICES_MK)

MSSMatMGUT_SLHA_INPUT := \
		$(DIR)/LesHouches.in.MSSMatMGUT_generated \
		$(DIR)/LesHouches.in.MSSMatMGUT

MSSMatMGUT_REFERENCES := \
		$(DIR)/MSSMatMGUT_references.tex

MSSMatMGUT_GNUPLOT := \
		$(DIR)/MSSMatMGUT_plot_rgflow.gnuplot \
		$(DIR)/MSSMatMGUT_plot_spectrum.gnuplot

MSSMatMGUT_TARBALL := \
		$(MODNAME).tar.gz

LIBMSSMatMGUT_SRC := \
		$(DIR)/MSSMatMGUT_a_muon.cpp \
		$(DIR)/MSSMatMGUT_edm.cpp \
		$(DIR)/MSSMatMGUT_FFV_form_factors.cpp \
		$(DIR)/MSSMatMGUT_l_to_lgamma.cpp \
		$(DIR)/MSSMatMGUT_effective_couplings.cpp \
		$(DIR)/MSSMatMGUT_info.cpp \
		$(DIR)/MSSMatMGUT_input_parameters.cpp \
		$(DIR)/MSSMatMGUT_mass_eigenstates.cpp \
		$(DIR)/MSSMatMGUT_observables.cpp \
		$(DIR)/MSSMatMGUT_physical.cpp \
		$(DIR)/MSSMatMGUT_slha_io.cpp \
		$(DIR)/MSSMatMGUT_soft_parameters.cpp \
		$(DIR)/MSSMatMGUT_susy_parameters.cpp \
		$(DIR)/MSSMatMGUT_utilities.cpp \
		$(DIR)/MSSMatMGUT_weinberg_angle.cpp

EXEMSSMatMGUT_SRC := \
		$(DIR)/run_MSSMatMGUT.cpp \
		$(DIR)/run_cmd_line_MSSMatMGUT.cpp \
		$(DIR)/scan_MSSMatMGUT.cpp
LLMSSMatMGUT_LIB  :=
LLMSSMatMGUT_OBJ  :=
LLMSSMatMGUT_SRC  := \
		$(DIR)/MSSMatMGUT_librarylink.cpp

LLMSSMatMGUT_MMA  := \
		$(DIR)/MSSMatMGUT_librarylink.m \
		$(DIR)/run_MSSMatMGUT.m

LIBMSSMatMGUT_HDR := \
		$(DIR)/MSSMatMGUT_a_muon.hpp \
		$(DIR)/MSSMatMGUT_convergence_tester.hpp \
		$(DIR)/MSSMatMGUT_edm.hpp \
		$(DIR)/MSSMatMGUT_FFV_form_factors.hpp \
		$(DIR)/MSSMatMGUT_l_to_lgamma.hpp \
		$(DIR)/MSSMatMGUT_effective_couplings.hpp \
		$(DIR)/MSSMatMGUT_ewsb_solver.hpp \
		$(DIR)/MSSMatMGUT_ewsb_solver_interface.hpp \
		$(DIR)/MSSMatMGUT_high_scale_constraint.hpp \
		$(DIR)/MSSMatMGUT_info.hpp \
		$(DIR)/MSSMatMGUT_initial_guesser.hpp \
		$(DIR)/MSSMatMGUT_input_parameters.hpp \
		$(DIR)/MSSMatMGUT_low_scale_constraint.hpp \
		$(DIR)/MSSMatMGUT_mass_eigenstates.hpp \
		$(DIR)/MSSMatMGUT_model.hpp \
		$(DIR)/MSSMatMGUT_model_slha.hpp \
		$(DIR)/MSSMatMGUT_observables.hpp \
		$(DIR)/MSSMatMGUT_physical.hpp \
		$(DIR)/MSSMatMGUT_slha_io.hpp \
		$(DIR)/MSSMatMGUT_spectrum_generator.hpp \
		$(DIR)/MSSMatMGUT_spectrum_generator_interface.hpp \
		$(DIR)/MSSMatMGUT_soft_parameters.hpp \
		$(DIR)/MSSMatMGUT_susy_parameters.hpp \
		$(DIR)/MSSMatMGUT_susy_scale_constraint.hpp \
		$(DIR)/MSSMatMGUT_utilities.hpp \
		$(DIR)/MSSMatMGUT_weinberg_angle.hpp

LIBMSSMatMGUT_CXXQFT_HDR := \
		$(DIR)/cxx_qft/MSSMatMGUT_qft.hpp \
		$(DIR)/cxx_qft/MSSMatMGUT_fields.hpp \
		$(DIR)/cxx_qft/MSSMatMGUT_vertices.hpp \
		$(DIR)/cxx_qft/MSSMatMGUT_context_base.hpp \
		$(DIR)/cxx_qft/MSSMatMGUT_npointfunctions.hpp

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
-include $(MSSMatMGUT_SUSY_BETAS_MK)
-include $(MSSMatMGUT_SOFT_BETAS_MK)
-include $(MSSMatMGUT_CXX_QFT_VERTICES_MK)
-include $(MSSMatMGUT_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(MSSMatMGUT_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(MSSMatMGUT_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(MSSMatMGUT_CXX_QFT_VERTICES_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(MSSMatMGUT_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBMSSMatMGUT_SRC := $(sort $(LIBMSSMatMGUT_SRC))
EXEMSSMatMGUT_SRC := $(sort $(EXEMSSMatMGUT_SRC))

LIBMSSMatMGUT_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBMSSMatMGUT_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBMSSMatMGUT_SRC)))

EXEMSSMatMGUT_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEMSSMatMGUT_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEMSSMatMGUT_SRC)))

EXEMSSMatMGUT_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXEMSSMatMGUT_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXEMSSMatMGUT_SRC)))

LIBMSSMatMGUT_DEP := \
		$(LIBMSSMatMGUT_OBJ:.o=.d)

EXEMSSMatMGUT_DEP := \
		$(EXEMSSMatMGUT_OBJ:.o=.d)

LLMSSMatMGUT_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLMSSMatMGUT_SRC)))

LLMSSMatMGUT_OBJ  := $(LLMSSMatMGUT_SRC:.cpp=.o)
LLMSSMatMGUT_LIB  := $(LLMSSMatMGUT_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBMSSMatMGUT     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_MSSMatMGUT := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_MSSMatMGUT := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL_DIR) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBMSSMatMGUT) $(EXEMSSMatMGUT_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(MSSMatMGUT_INSTALL_DIR)
		$(Q)install -d $(MSSMatMGUT_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBMSSMatMGUT_SRC) $(MSSMatMGUT_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBMSSMatMGUT_HDR) $(MSSMatMGUT_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBMSSMatMGUT_CXXQFT_HDR) $(MSSMatMGUT_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(EXEMSSMatMGUT_SRC) $(MSSMatMGUT_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLMSSMatMGUT_SRC) $(MSSMatMGUT_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLMSSMatMGUT_MMA) $(MSSMatMGUT_INSTALL_DIR)
		$(Q)$(INSTALL_STRIPPED) $(MSSMatMGUT_MK) $(MSSMatMGUT_INSTALL_DIR) -m u=rw,g=r,o=r
		$(Q)install -m u=rw,g=r,o=r $(MSSMatMGUT_INCLUDE_MK) $(MSSMatMGUT_INSTALL_DIR)
ifneq ($(MSSMatMGUT_SLHA_INPUT),)
		$(Q)install -m u=rw,g=r,o=r $(MSSMatMGUT_SLHA_INPUT) $(MSSMatMGUT_INSTALL_DIR)
endif
		$(Q)install -m u=rw,g=r,o=r $(MSSMatMGUT_REFERENCES) $(MSSMatMGUT_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(MSSMatMGUT_GNUPLOT) $(MSSMatMGUT_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		$(Q)-rm -f $(LIBMSSMatMGUT_DEP)
		$(Q)-rm -f $(EXEMSSMatMGUT_DEP)
		$(Q)-rm -f $(LLMSSMatMGUT_DEP)

clean-$(MODNAME)-lib:
		$(Q)-rm -f $(LIBMSSMatMGUT)
		$(Q)-rm -f $(LLMSSMatMGUT_LIB)

clean-$(MODNAME)-obj:
		$(Q)-rm -f $(LIBMSSMatMGUT_OBJ)
		$(Q)-rm -f $(EXEMSSMatMGUT_OBJ)
		$(Q)-rm -f $(LLMSSMatMGUT_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		$(Q)-rm -f $(LIBMSSMatMGUT_SRC)
		$(Q)-rm -f $(LIBMSSMatMGUT_HDR)
		$(Q)-rm -f $(LIBMSSMatMGUT_CXXQFT_HDR)
		$(Q)-rm -f $(EXEMSSMatMGUT_SRC)
		$(Q)-rm -f $(LLMSSMatMGUT_SRC)
		$(Q)-rm -f $(LLMSSMatMGUT_MMA)
		$(Q)-rm -f $(METACODE_STAMP_MSSMatMGUT)
		$(Q)-rm -f $(MSSMatMGUT_INCLUDE_MK)
		$(Q)-rm -f $(MSSMatMGUT_SLHA_INPUT)
		$(Q)-rm -f $(MSSMatMGUT_REFERENCES)
		$(Q)-rm -f $(MSSMatMGUT_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		$(Q)-rm -f $(EXEMSSMatMGUT_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		$(Q)tar -czf $(MSSMatMGUT_TARBALL) \
		$(LIBMSSMatMGUT_SRC) $(LIBMSSMatMGUT_HDR) $(LIBMSSMatMGUT_CXXQFT_HDR) \
		$(EXEMSSMatMGUT_SRC) \
		$(LLMSSMatMGUT_SRC) $(LLMSSMatMGUT_MMA) \
		$(MSSMatMGUT_MK) $(MSSMatMGUT_INCLUDE_MK) \
		$(MSSMatMGUT_SLHA_INPUT) $(MSSMatMGUT_REFERENCES) \
		$(MSSMatMGUT_GNUPLOT)

$(LIBMSSMatMGUT_SRC) $(LIBMSSMatMGUT_HDR) $(LIBMSSMatMGUT_CXXQFT_HDR) $(EXEMSSMatMGUT_SRC) $(LLMSSMatMGUT_SRC) $(LLMSSMatMGUT_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_MSSMatMGUT)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_MSSMatMGUT): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_MSSMatMGUT)
		@$(MSG)
		$(Q)printf "%s" "Get[\"$<\"]; Quit[]" | "$(MATH)" || (echo "Error: The code generation failed!"; exit 1)
		@touch "$(METACODE_STAMP_MSSMatMGUT)"
		@echo "Note: to regenerate MSSMatMGUT source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_MSSMatMGUT)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_MSSMatMGUT):
		@true
endif

$(LIBMSSMatMGUT_DEP) $(EXEMSSMatMGUT_DEP) $(LLMSSMatMGUT_DEP) $(LIBMSSMatMGUT_OBJ) $(EXEMSSMatMGUT_OBJ) $(LLMSSMatMGUT_OBJ) $(LLMSSMatMGUT_LIB): \
	CPPFLAGS += $(MODMSSMatMGUT_SUBMOD_INC) $(MODMSSMatMGUT_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBMSSMatMGUT_DEP) $(EXEMSSMatMGUT_DEP) $(LLMSSMatMGUT_DEP) $(LIBMSSMatMGUT_OBJ) $(EXEMSSMatMGUT_OBJ) $(LLMSSMatMGUT_OBJ) $(LLMSSMatMGUT_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLMSSMatMGUT_OBJ) $(LLMSSMatMGUT_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBMSSMatMGUT): $(LIBMSSMatMGUT_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBMSSMatMGUT) $(MODMSSMatMGUT_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		@$(MSG)
		$(Q)$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLMSSMatMGUT_LIB): $(LLMSSMatMGUT_OBJ) $(LIBMSSMatMGUT) $(MODMSSMatMGUT_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		@$(MSG)
		$(Q)$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS)

ALLDEP += $(LIBMSSMatMGUT_DEP) $(EXEMSSMatMGUT_DEP)
ALLSRC += $(LIBMSSMatMGUT_SRC) $(EXEMSSMatMGUT_SRC)
ALLLIB += $(LIBMSSMatMGUT)
ALLEXE += $(EXEMSSMatMGUT_EXE)
ALLMODDEP += $(MODMSSMatMGUT_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLMSSMatMGUT_DEP)
ALLSRC += $(LLMSSMatMGUT_SRC)
ALLLL  += $(LLMSSMatMGUT_LIB)
endif
