DIR          := models/CMSSMNoFV
MODNAME      := CMSSMNoFV
SARAH_MODEL  := MSSMNoFV
WITH_$(MODNAME) := yes
MODCMSSMNoFV_MOD := SM MSSM_higgs
MODCMSSMNoFV_DEP := $(patsubst %,model_specific/%,$(MODCMSSMNoFV_MOD))
MODCMSSMNoFV_INC := $(patsubst %,-Imodel_specific/%,$(MODCMSSMNoFV_MOD))
MODCMSSMNoFV_LIB := $(foreach M,$(MODCMSSMNoFV_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

MODCMSSMNoFV_SUBMOD  := $(DIR)/cxx_qft
MODCMSSMNoFV_SUBMOD_INC := $(patsubst %,-I%,$(MODCMSSMNoFV_SUBMOD))

CMSSMNoFV_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)
CMSSMNoFV_INSTALL_CXXQFT_DIR := \
		$(CMSSMNoFV_INSTALL_DIR)/cxx_qft

CMSSMNoFV_MK     := \
		$(DIR)/module.mk

CMSSMNoFV_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

CMSSMNoFV_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

CMSSMNoFV_CXXQFT_VERTICES_MK := \
		$(DIR)/cxx_qft/vertices.mk

-include $(CMSSMNoFV_CXXQFT_VERTICES_MK)
LIBCMSSMNoFV_CXXQFT_VERTICES_SRC ?= ''

CMSSMNoFV_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

CMSSMNoFV_INCLUDE_MK := \
		$(CMSSMNoFV_SUSY_BETAS_MK) \
		$(CMSSMNoFV_SOFT_BETAS_MK)

CMSSMNoFV_SLHA_INPUT := \
		$(DIR)/LesHouches.in.CMSSMNoFV_generated \
		$(DIR)/LesHouches.in.CMSSMNoFV

CMSSMNoFV_REFERENCES := \
		$(DIR)/CMSSMNoFV_references.tex

CMSSMNoFV_GNUPLOT := \
		$(DIR)/CMSSMNoFV_plot_rgflow.gnuplot \
		$(DIR)/CMSSMNoFV_plot_spectrum.gnuplot

CMSSMNoFV_TARBALL := \
		$(MODNAME).tar.gz

LIBCMSSMNoFV_SRC := \
		$(DIR)/CMSSMNoFV_a_muon.cpp \
		$(DIR)/CMSSMNoFV_edm.cpp \
		$(DIR)/CMSSMNoFV_FFV_form_factors.cpp \
		$(DIR)/CMSSMNoFV_f_to_f_conversion.cpp \
		$(DIR)/CMSSMNoFV_l_to_lgamma.cpp \
		$(DIR)/CMSSMNoFV_b_to_s_gamma.cpp \
		$(DIR)/CMSSMNoFV_effective_couplings.cpp \
		$(DIR)/CMSSMNoFV_info.cpp \
		$(DIR)/CMSSMNoFV_input_parameters.cpp \
		$(DIR)/CMSSMNoFV_mass_eigenstates.cpp \
		$(DIR)/CMSSMNoFV_mass_eigenstates_decoupling_scheme.cpp \
		$(DIR)/CMSSMNoFV_model_slha.cpp \
		$(DIR)/CMSSMNoFV_observables.cpp \
		$(DIR)/CMSSMNoFV_physical.cpp \
		$(DIR)/CMSSMNoFV_slha_io.cpp \
		$(DIR)/CMSSMNoFV_soft_parameters.cpp \
		$(DIR)/CMSSMNoFV_susy_parameters.cpp \
		$(DIR)/CMSSMNoFV_utilities.cpp \
		$(DIR)/CMSSMNoFV_weinberg_angle.cpp

LIBCMSSMNoFV_SRC += $(LIBCMSSMNoFV_CXXQFT_VERTICES_SRC)

EXECMSSMNoFV_SRC := \
		$(DIR)/run_CMSSMNoFV.cpp \
		$(DIR)/run_cmd_line_CMSSMNoFV.cpp \
		$(DIR)/scan_CMSSMNoFV.cpp
LLCMSSMNoFV_LIB  :=
LLCMSSMNoFV_OBJ  :=
LLCMSSMNoFV_SRC  := \
		$(DIR)/CMSSMNoFV_librarylink.cpp

LLCMSSMNoFV_MMA  := \
		$(DIR)/CMSSMNoFV_librarylink.m \
		$(DIR)/run_CMSSMNoFV.m

LIBCMSSMNoFV_HDR := \
		$(DIR)/CMSSMNoFV_a_muon.hpp \
		$(DIR)/CMSSMNoFV_convergence_tester.hpp \
		$(DIR)/CMSSMNoFV_edm.hpp \
		$(DIR)/CMSSMNoFV_FFV_form_factors.hpp \
		$(DIR)/CMSSMNoFV_f_to_f_conversion.hpp \
		$(DIR)/CMSSMNoFV_l_to_lgamma.hpp \
		$(DIR)/CMSSMNoFV_b_to_s_gamma.hpp \
		$(DIR)/CMSSMNoFV_effective_couplings.hpp \
		$(DIR)/CMSSMNoFV_ewsb_solver.hpp \
		$(DIR)/CMSSMNoFV_ewsb_solver_interface.hpp \
		$(DIR)/CMSSMNoFV_high_scale_constraint.hpp \
		$(DIR)/CMSSMNoFV_info.hpp \
		$(DIR)/CMSSMNoFV_initial_guesser.hpp \
		$(DIR)/CMSSMNoFV_input_parameters.hpp \
		$(DIR)/CMSSMNoFV_low_scale_constraint.hpp \
		$(DIR)/CMSSMNoFV_mass_eigenstates.hpp \
		$(DIR)/CMSSMNoFV_mass_eigenstates_interface.hpp \
		$(DIR)/CMSSMNoFV_mass_eigenstates_decoupling_scheme.hpp \
		$(DIR)/CMSSMNoFV_model.hpp \
		$(DIR)/CMSSMNoFV_model_slha.hpp \
		$(DIR)/CMSSMNoFV_observables.hpp \
		$(DIR)/CMSSMNoFV_physical.hpp \
		$(DIR)/CMSSMNoFV_slha_io.hpp \
		$(DIR)/CMSSMNoFV_spectrum_generator.hpp \
		$(DIR)/CMSSMNoFV_spectrum_generator_interface.hpp \
		$(DIR)/CMSSMNoFV_soft_parameters.hpp \
		$(DIR)/CMSSMNoFV_susy_parameters.hpp \
		$(DIR)/CMSSMNoFV_susy_scale_constraint.hpp \
		$(DIR)/CMSSMNoFV_utilities.hpp \
		$(DIR)/CMSSMNoFV_weinberg_angle.hpp

LIBCMSSMNoFV_CXXQFT_HDR := \
		$(DIR)/cxx_qft/CMSSMNoFV_qft.hpp \
		$(DIR)/cxx_qft/CMSSMNoFV_fields.hpp \
		$(DIR)/cxx_qft/CMSSMNoFV_vertices.hpp \
		$(DIR)/cxx_qft/CMSSMNoFV_context_base.hpp \
		$(DIR)/cxx_qft/CMSSMNoFV_npointfunctions_wilsoncoeffs.hpp

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
-include $(CMSSMNoFV_SUSY_BETAS_MK)
-include $(CMSSMNoFV_SOFT_BETAS_MK)
-include $(CMSSMNoFV_CXXQFT_VERTICES_MK)
-include $(CMSSMNoFV_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(CMSSMNoFV_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(CMSSMNoFV_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(CMSSMNoFV_CXXQFT_VERTICES_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(CMSSMNoFV_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBCMSSMNoFV_SRC := $(sort $(LIBCMSSMNoFV_SRC))
EXECMSSMNoFV_SRC := $(sort $(EXECMSSMNoFV_SRC))

LIBCMSSMNoFV_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBCMSSMNoFV_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBCMSSMNoFV_SRC)))

EXECMSSMNoFV_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXECMSSMNoFV_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXECMSSMNoFV_SRC)))

EXECMSSMNoFV_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXECMSSMNoFV_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXECMSSMNoFV_SRC)))

LIBCMSSMNoFV_DEP := \
		$(LIBCMSSMNoFV_OBJ:.o=.d)

EXECMSSMNoFV_DEP := \
		$(EXECMSSMNoFV_OBJ:.o=.d)

LLCMSSMNoFV_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLCMSSMNoFV_SRC)))

LLCMSSMNoFV_OBJ  := $(LLCMSSMNoFV_SRC:.cpp=.o)
LLCMSSMNoFV_LIB  := $(LLCMSSMNoFV_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBCMSSMNoFV     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_CMSSMNoFV := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_CMSSMNoFV := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL_DIR) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBCMSSMNoFV) $(EXECMSSMNoFV_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(CMSSMNoFV_INSTALL_DIR)
		$(Q)install -d $(CMSSMNoFV_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBCMSSMNoFV_SRC) $(CMSSMNoFV_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBCMSSMNoFV_CXXQFT_VERTICES_SRC) $(CMSSMNoFV_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBCMSSMNoFV_HDR) $(CMSSMNoFV_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBCMSSMNoFV_CXXQFT_HDR) $(CMSSMNoFV_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(EXECMSSMNoFV_SRC) $(CMSSMNoFV_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLCMSSMNoFV_SRC) $(CMSSMNoFV_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLCMSSMNoFV_MMA) $(CMSSMNoFV_INSTALL_DIR)
		$(Q)$(INSTALL_STRIPPED) $(CMSSMNoFV_MK) $(CMSSMNoFV_INSTALL_DIR) -m u=rw,g=r,o=r
		$(Q)install -m u=rw,g=r,o=r $(CMSSMNoFV_INCLUDE_MK) $(CMSSMNoFV_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(CMSSMNoFV_CXXQFT_VERTICES_MK) $(CMSSMNoFV_INSTALL_CXXQFT_DIR)

ifneq ($(CMSSMNoFV_SLHA_INPUT),)
		$(Q)install -m u=rw,g=r,o=r $(CMSSMNoFV_SLHA_INPUT) $(CMSSMNoFV_INSTALL_DIR)
endif
		$(Q)install -m u=rw,g=r,o=r $(CMSSMNoFV_REFERENCES) $(CMSSMNoFV_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(CMSSMNoFV_GNUPLOT) $(CMSSMNoFV_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		$(Q)-rm -f $(LIBCMSSMNoFV_DEP)
		$(Q)-rm -f $(EXECMSSMNoFV_DEP)
		$(Q)-rm -f $(LLCMSSMNoFV_DEP)

clean-$(MODNAME)-lib:
		$(Q)-rm -f $(LIBCMSSMNoFV)
		$(Q)-rm -f $(LLCMSSMNoFV_LIB)

clean-$(MODNAME)-obj:
		$(Q)-rm -f $(LIBCMSSMNoFV_OBJ)
		$(Q)-rm -f $(EXECMSSMNoFV_OBJ)
		$(Q)-rm -f $(LLCMSSMNoFV_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		$(Q)-rm -f $(LIBCMSSMNoFV_SRC)
		$(Q)-rm -f $(LIBCMSSMNoFV_HDR)
		$(Q)-rm -f $(LIBCMSSMNoFV_CXXQFT_HDR)
		$(Q)-rm -f $(EXECMSSMNoFV_SRC)
		$(Q)-rm -f $(LLCMSSMNoFV_SRC)
		$(Q)-rm -f $(LLCMSSMNoFV_MMA)
		$(Q)-rm -f $(METACODE_STAMP_CMSSMNoFV)
		$(Q)-rm -f $(CMSSMNoFV_INCLUDE_MK)
		$(Q)-rm -f $(CMSSMNoFV_CXXQFT_VERTICES_MK)
		$(Q)-rm -f $(CMSSMNoFV_SLHA_INPUT)
		$(Q)-rm -f $(CMSSMNoFV_REFERENCES)
		$(Q)-rm -f $(CMSSMNoFV_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		$(Q)-rm -f $(EXECMSSMNoFV_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		$(Q)tar -czf $(CMSSMNoFV_TARBALL) \
		$(LIBCMSSMNoFV_SRC) $(LIBCMSSMNoFV_HDR) $(LIBCMSSMNoFV_CXXQFT_HDR) \
		$(EXECMSSMNoFV_SRC) \
		$(LLCMSSMNoFV_SRC) $(LLCMSSMNoFV_MMA) \
		$(CMSSMNoFV_MK) $(CMSSMNoFV_INCLUDE_MK) $(CMSSMNoFV_CXXQFT_VERTICES_MK) \
		$(CMSSMNoFV_SLHA_INPUT) $(CMSSMNoFV_REFERENCES) \
		$(CMSSMNoFV_GNUPLOT)

$(LIBCMSSMNoFV_SRC) $(LIBCMSSMNoFV_HDR) $(LIBCMSSMNoFV_CXXQFT_HDR) $(EXECMSSMNoFV_SRC) $(LLCMSSMNoFV_SRC) $(LLCMSSMNoFV_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_CMSSMNoFV)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_CMSSMNoFV): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_CMSSMNoFV)
		@$(MSG)
		$(Q)printf "%s" "Get[\"$<\"]; Quit[]" | "$(MATH)" || (echo "Error: The code generation failed!"; exit 1)
		@touch "$(METACODE_STAMP_CMSSMNoFV)"
		@echo "Note: to regenerate CMSSMNoFV source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_CMSSMNoFV)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_CMSSMNoFV):
		@true
endif

$(LIBCMSSMNoFV_DEP) $(EXECMSSMNoFV_DEP) $(LLCMSSMNoFV_DEP) $(LIBCMSSMNoFV_OBJ) $(EXECMSSMNoFV_OBJ) $(LLCMSSMNoFV_OBJ) $(LLCMSSMNoFV_LIB): \
	CPPFLAGS += $(MODCMSSMNoFV_SUBMOD_INC) $(MODCMSSMNoFV_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(GM2CALCFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBCMSSMNoFV_DEP) $(EXECMSSMNoFV_DEP) $(LLCMSSMNoFV_DEP) $(LIBCMSSMNoFV_OBJ) $(EXECMSSMNoFV_OBJ) $(LLCMSSMNoFV_OBJ) $(LLCMSSMNoFV_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLCMSSMNoFV_OBJ) $(LLCMSSMNoFV_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBCMSSMNoFV): $(LIBCMSSMNoFV_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBCMSSMNoFV) $(MODCMSSMNoFV_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(FUTILIBS)
		@$(MSG)
		$(Q)$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^) $(filter -%,$(LOOPFUNCLIBS)) $(GM2CALCLIBS) $(HIMALAYALIBS) $(GSLLIBS) $(SQLITELIBS) $(TSILLIBS) $(FLIBS) $(THREADLIBS) $(LDLIBS) $(FUTILIBS)

$(LLCMSSMNoFV_LIB): $(LLCMSSMNoFV_OBJ) $(LIBCMSSMNoFV) $(MODCMSSMNoFV_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(FUTILIBS)
		@$(MSG)
		$(Q)$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^) $(filter -%,$(LOOPFUNCLIBS)) $(GM2CALCLIBS) $(HIMALAYALIBS) $(TSILLIBS) $(GSLLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS) $(FUTILIBS) $(FLIBS)

ALLDEP += $(LIBCMSSMNoFV_DEP) $(EXECMSSMNoFV_DEP)
ALLSRC += $(LIBCMSSMNoFV_SRC) $(EXECMSSMNoFV_SRC)
ALLLIB += $(LIBCMSSMNoFV)
ALLEXE += $(EXECMSSMNoFV_EXE)
ALLMODDEP += $(MODCMSSMNoFV_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLCMSSMNoFV_DEP)
ALLSRC += $(LLCMSSMNoFV_SRC)
ALLLL  += $(LLCMSSMNoFV_LIB)
endif
