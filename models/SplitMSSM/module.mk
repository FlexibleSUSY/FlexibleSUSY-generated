DIR          := models/SplitMSSM
MODNAME      := SplitMSSM
SARAH_MODEL  := SplitMSSM
WITH_$(MODNAME) := yes
MODSplitMSSM_MOD := SM SplitMSSM
MODSplitMSSM_DEP := $(patsubst %,model_specific/%,$(MODSplitMSSM_MOD))
MODSplitMSSM_INC := $(patsubst %,-Imodel_specific/%,$(MODSplitMSSM_MOD))
MODSplitMSSM_LIB := $(foreach M,$(MODSplitMSSM_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

MODSplitMSSM_SUBMOD  := $(DIR)/cxx_qft
MODSplitMSSM_SUBMOD_INC := $(patsubst %,-I%,$(MODSplitMSSM_SUBMOD))

SplitMSSM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)
SplitMSSM_INSTALL_CXXQFT_DIR := \
		$(SplitMSSM_INSTALL_DIR)/cxx_qft

SplitMSSM_MK     := \
		$(DIR)/module.mk

SplitMSSM_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

SplitMSSM_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

SplitMSSM_CXXQFT_VERTICES_MK := \
		$(DIR)/cxx_qft/vertices.mk

-include $(SplitMSSM_CXXQFT_VERTICES_MK)
LIBSplitMSSM_CXXQFT_VERTICES_SRC ?= ''

SplitMSSM_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

SplitMSSM_FlexibleDecay_MK := \
		$(DIR)/decays/FlexibleDecay.mk

SplitMSSM_INCLUDE_MK := \
		$(SplitMSSM_SUSY_BETAS_MK) \
		$(SplitMSSM_SOFT_BETAS_MK)

SplitMSSM_SLHA_INPUT := \
		$(DIR)/LesHouches.in.SplitMSSM_generated \
		$(DIR)/LesHouches.in.SplitMSSM

SplitMSSM_REFERENCES := \
		$(DIR)/SplitMSSM_references.tex

SplitMSSM_GNUPLOT := \
		$(DIR)/SplitMSSM_plot_rgflow.gnuplot \
		$(DIR)/SplitMSSM_plot_spectrum.gnuplot

SplitMSSM_TARBALL := \
		$(MODNAME).tar.gz

LIBSplitMSSM_SRC := \
		$(DIR)/SplitMSSM_amm.cpp \
		$(DIR)/SplitMSSM_edm.cpp \
		$(DIR)/SplitMSSM_FFV_form_factors.cpp \
		$(wildcard $(DIR)/observables/SplitMSSM*.cpp) \
		$(DIR)/SplitMSSM_b_to_s_gamma.cpp \
		$(DIR)/SplitMSSM_info.cpp \
		$(DIR)/SplitMSSM_input_parameters.cpp \
		$(DIR)/SplitMSSM_mass_eigenstates.cpp \
		$(DIR)/SplitMSSM_mass_eigenstates_decoupling_scheme.cpp \
		$(DIR)/SplitMSSM_model_slha.cpp \
		$(DIR)/SplitMSSM_lepton_amm_wrapper.cpp \
		$(DIR)/SplitMSSM_observables.cpp \
		$(DIR)/SplitMSSM_physical.cpp \
		$(DIR)/SplitMSSM_slha_io.cpp \
		$(DIR)/SplitMSSM_soft_parameters.cpp \
		$(DIR)/SplitMSSM_susy_parameters.cpp \
		$(DIR)/SplitMSSM_unitarity.cpp \
		$(DIR)/SplitMSSM_utilities.cpp \
		$(DIR)/SplitMSSM_weinberg_angle.cpp

LIBSplitMSSM_SRC += $(LIBSplitMSSM_CXXQFT_VERTICES_SRC)

EXESplitMSSM_SRC := \
		$(DIR)/run_SplitMSSM.cpp \
		$(DIR)/run_cmd_line_SplitMSSM.cpp \
		$(DIR)/scan_SplitMSSM.cpp
LLSplitMSSM_LIB  :=
LLSplitMSSM_OBJ  :=
LLSplitMSSM_SRC  := \
		$(DIR)/SplitMSSM_librarylink.cpp

LLSplitMSSM_MMA  := \
		$(DIR)/SplitMSSM_librarylink.m \
		$(DIR)/run_SplitMSSM.m

LIBSplitMSSM_HDR := \
		$(DIR)/SplitMSSM_amm.hpp \
		$(DIR)/SplitMSSM_convergence_tester.hpp \
		$(DIR)/SplitMSSM_edm.hpp \
		$(DIR)/SplitMSSM_FFV_form_factors.hpp \
		$(wildcard $(DIR)/observables/SplitMSSM*.hpp) \
		$(DIR)/SplitMSSM_b_to_s_gamma.hpp \
		$(DIR)/SplitMSSM_ewsb_solver.hpp \
		$(DIR)/SplitMSSM_ewsb_solver_interface.hpp \
		$(DIR)/SplitMSSM_high_scale_constraint.hpp \
		$(DIR)/SplitMSSM_info.hpp \
		$(DIR)/SplitMSSM_initial_guesser.hpp \
		$(DIR)/SplitMSSM_input_parameters.hpp \
		$(DIR)/SplitMSSM_low_scale_constraint.hpp \
		$(DIR)/SplitMSSM_mass_eigenstates.hpp \
		$(DIR)/SplitMSSM_mass_eigenstates_interface.hpp \
		$(DIR)/SplitMSSM_mass_eigenstates_decoupling_scheme.hpp \
		$(DIR)/SplitMSSM_model.hpp \
		$(DIR)/SplitMSSM_model_slha.hpp \
		$(DIR)/SplitMSSM_lepton_amm_wrapper.hpp \
		$(DIR)/SplitMSSM_observables.hpp \
		$(DIR)/SplitMSSM_physical.hpp \
		$(DIR)/SplitMSSM_slha_io.hpp \
		$(DIR)/SplitMSSM_spectrum_generator.hpp \
		$(DIR)/SplitMSSM_spectrum_generator_interface.hpp \
		$(DIR)/SplitMSSM_soft_parameters.hpp \
		$(DIR)/SplitMSSM_susy_parameters.hpp \
		$(DIR)/SplitMSSM_susy_scale_constraint.hpp \
		$(DIR)/SplitMSSM_unitarity.hpp \
		$(DIR)/SplitMSSM_utilities.hpp \
		$(DIR)/SplitMSSM_weinberg_angle.hpp

LIBSplitMSSM_CXXQFT_HDR := \
		$(DIR)/cxx_qft/SplitMSSM_qft.hpp \
		$(DIR)/cxx_qft/SplitMSSM_fields.hpp \
		$(DIR)/cxx_qft/SplitMSSM_particle_aliases.hpp \
		$(DIR)/cxx_qft/SplitMSSM_vertices.hpp \
		$(DIR)/cxx_qft/SplitMSSM_context_base.hpp \
		$(DIR)/cxx_qft/SplitMSSM_npointfunctions_wilsoncoeffs.hpp

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
-include $(SplitMSSM_SUSY_BETAS_MK)
-include $(SplitMSSM_SOFT_BETAS_MK)
-include $(SplitMSSM_FlexibleDecay_MK)
-include $(SplitMSSM_CXXQFT_VERTICES_MK)
-include $(SplitMSSM_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(SplitMSSM_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(SplitMSSM_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@

$(SplitMSSM_FlexibleDecay_MK): run-metacode-$(MODNAME)
$(SplitMSSM_CXXQFT_VERTICES_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(SplitMSSM_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBSplitMSSM_SRC := $(sort $(LIBSplitMSSM_SRC))
EXESplitMSSM_SRC := $(sort $(EXESplitMSSM_SRC))

LIBSplitMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBSplitMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBSplitMSSM_SRC)))

EXESplitMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXESplitMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXESplitMSSM_SRC)))

EXESplitMSSM_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXESplitMSSM_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXESplitMSSM_SRC)))

LIBSplitMSSM_DEP := \
		$(LIBSplitMSSM_OBJ:.o=.d)

EXESplitMSSM_DEP := \
		$(EXESplitMSSM_OBJ:.o=.d)

LLSplitMSSM_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLSplitMSSM_SRC)))

LLSplitMSSM_OBJ  := $(LLSplitMSSM_SRC:.cpp=.o)
LLSplitMSSM_LIB  := $(LLSplitMSSM_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBSplitMSSM     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_SplitMSSM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_SplitMSSM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL_DIR) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBSplitMSSM) $(EXESplitMSSM_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(SplitMSSM_INSTALL_DIR)
		$(Q)install -d $(SplitMSSM_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBSplitMSSM_SRC) $(SplitMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBSplitMSSM_CXXQFT_VERTICES_SRC) $(SplitMSSM_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBSplitMSSM_HDR) $(SplitMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBSplitMSSM_CXXQFT_HDR) $(SplitMSSM_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(EXESplitMSSM_SRC) $(SplitMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLSplitMSSM_SRC) $(SplitMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLSplitMSSM_MMA) $(SplitMSSM_INSTALL_DIR)
		$(Q)$(INSTALL_STRIPPED) $(SplitMSSM_MK) $(SplitMSSM_INSTALL_DIR) -m u=rw,g=r,o=r
		$(Q)install -m u=rw,g=r,o=r $(SplitMSSM_INCLUDE_MK) $(SplitMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(SplitMSSM_CXXQFT_VERTICES_MK) $(SplitMSSM_INSTALL_CXXQFT_DIR)

ifneq ($(SplitMSSM_SLHA_INPUT),)
		$(Q)install -m u=rw,g=r,o=r $(SplitMSSM_SLHA_INPUT) $(SplitMSSM_INSTALL_DIR)
endif
		$(Q)install -m u=rw,g=r,o=r $(SplitMSSM_REFERENCES) $(SplitMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(SplitMSSM_GNUPLOT) $(SplitMSSM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		$(Q)-rm -f $(LIBSplitMSSM_DEP)
		$(Q)-rm -f $(EXESplitMSSM_DEP)
		$(Q)-rm -f $(LLSplitMSSM_DEP)

clean-$(MODNAME)-lib:
		$(Q)-rm -f $(LIBSplitMSSM)
		$(Q)-rm -f $(LLSplitMSSM_LIB)

clean-$(MODNAME)-obj:
		$(Q)-rm -f $(LIBSplitMSSM_OBJ)
		$(Q)-rm -f $(EXESplitMSSM_OBJ)
		$(Q)-rm -f $(LLSplitMSSM_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		$(Q)-rm -f $(LIBSplitMSSM_SRC)
		$(Q)-rm -f $(LIBSplitMSSM_HDR)
		$(Q)-rm -f $(LIBSplitMSSM_CXXQFT_HDR)
		$(Q)-rm -f $(EXESplitMSSM_SRC)
		$(Q)-rm -f $(LLSplitMSSM_SRC)
		$(Q)-rm -f $(LLSplitMSSM_MMA)
		$(Q)-rm -f $(METACODE_STAMP_SplitMSSM)
		$(Q)-rm -f $(SplitMSSM_INCLUDE_MK)
		$(Q)-rm -f $(SplitMSSM_CXXQFT_VERTICES_MK)
		$(Q)-rm -f $(SplitMSSM_SLHA_INPUT)
		$(Q)-rm -f $(SplitMSSM_REFERENCES)
		$(Q)-rm -f $(SplitMSSM_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		$(Q)-rm -f $(EXESplitMSSM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		$(Q)tar -czf $(SplitMSSM_TARBALL) \
		$(LIBSplitMSSM_SRC) $(LIBSplitMSSM_HDR) $(LIBSplitMSSM_CXXQFT_HDR) \
		$(EXESplitMSSM_SRC) \
		$(LLSplitMSSM_SRC) $(LLSplitMSSM_MMA) \
		$(SplitMSSM_MK) $(SplitMSSM_INCLUDE_MK) $(SplitMSSM_CXXQFT_VERTICES_MK) \
		$(SplitMSSM_SLHA_INPUT) $(SplitMSSM_REFERENCES) \
		$(SplitMSSM_GNUPLOT) \
		$(SplitMSSM_FlexibleDecay_MK)

$(LIBSplitMSSM_SRC) $(LIBSplitMSSM_HDR) $(LIBSplitMSSM_CXXQFT_HDR) $(EXESplitMSSM_SRC) $(LLSplitMSSM_SRC) $(LLSplitMSSM_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_SplitMSSM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_SplitMSSM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_SplitMSSM)
		@$(MSG)
		$(Q)printf "%s" "Get[\"$<\"]; Quit[]" | "$(MATH)" || (echo "Error: The code generation failed!"; exit 1)
		@touch "$(METACODE_STAMP_SplitMSSM)"
		@echo "Note: to regenerate SplitMSSM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_SplitMSSM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_SplitMSSM):
		@true
endif

$(LIBSplitMSSM_DEP) $(EXESplitMSSM_DEP) $(LLSplitMSSM_DEP) $(LIBSplitMSSM_OBJ) $(EXESplitMSSM_OBJ) $(LLSplitMSSM_OBJ) $(LLSplitMSSM_LIB): \
	CPPFLAGS += $(MODSplitMSSM_SUBMOD_INC) $(MODSplitMSSM_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(GM2CALCFLAGS) $(HIGGSTOOLSFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBSplitMSSM_DEP) $(EXESplitMSSM_DEP) $(LLSplitMSSM_DEP) $(LIBSplitMSSM_OBJ) $(EXESplitMSSM_OBJ) $(LLSplitMSSM_OBJ) $(LLSplitMSSM_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLSplitMSSM_OBJ) $(LLSplitMSSM_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBSplitMSSM): $(LIBSplitMSSM_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBSplitMSSM) $(MODSplitMSSM_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(FUTILIBS)
		@$(MSG)
		$(Q)$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^) $(filter -%,$(LOOPFUNCLIBS)) $(GM2CALCLIBS) $(HIGGSTOOLSLIBS) $(PYTHONLIBS) $(HIMALAYALIBS) $(GSLLIBS) $(SQLITELIBS) $(TSILLIBS) $(FLIBS) $(THREADLIBS) $(LDLIBS) $(FUTILIBS)

$(LLSplitMSSM_LIB): $(LLSplitMSSM_OBJ) $(LIBSplitMSSM) $(MODSplitMSSM_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(FUTILIBS)
		@$(MSG)
		$(Q)$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^) $(filter -%,$(LOOPFUNCLIBS)) $(GM2CALCLIBS) $(HIGGSTOOLSLIBS) $(PYTHONLIBS) $(HIMALAYALIBS) $(TSILLIBS) $(GSLLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS) $(FUTILIBS) $(FLIBS)

ALLDEP += $(LIBSplitMSSM_DEP) $(EXESplitMSSM_DEP)
ALLSRC += $(LIBSplitMSSM_SRC) $(EXESplitMSSM_SRC)
ALLLIB += $(LIBSplitMSSM)
ALLEXE += $(EXESplitMSSM_EXE)
ALLMODDEP += $(MODSplitMSSM_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLSplitMSSM_DEP)
ALLSRC += $(LLSplitMSSM_SRC)
ALLLL  += $(LLSplitMSSM_LIB)
endif
