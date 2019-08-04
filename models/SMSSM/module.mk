DIR          := models/SMSSM
MODNAME      := SMSSM
SARAH_MODEL  := SMSSM
WITH_$(MODNAME) := yes
MODSMSSM_MOD := SM MSSM_higgs NMSSM_higgs
MODSMSSM_DEP := $(patsubst %,model_specific/%,$(MODSMSSM_MOD))
MODSMSSM_INC := $(patsubst %,-Imodel_specific/%,$(MODSMSSM_MOD))
MODSMSSM_LIB := $(foreach M,$(MODSMSSM_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

MODSMSSM_SUBMOD  := $(DIR)/cxx_qft
MODSMSSM_SUBMOD_INC := $(patsubst %,-I%,$(MODSMSSM_SUBMOD))

SMSSM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)
SMSSM_INSTALL_CXXQFT_DIR := \
		$(SMSSM_INSTALL_DIR)/cxx_qft

SMSSM_MK     := \
		$(DIR)/module.mk

SMSSM_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

SMSSM_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

SMSSM_CXX_QFT_VERTICES_MK := \
		$(DIR)/cxx_qft/vertices.mk

SMSSM_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

SMSSM_INCLUDE_MK := \
		$(SMSSM_SUSY_BETAS_MK) \
		$(SMSSM_SOFT_BETAS_MK) \
		$(SMSSM_CXX_QFT_VERTICES_MK)

SMSSM_SLHA_INPUT := \
		$(DIR)/LesHouches.in.SMSSM_generated \
		$(DIR)/LesHouches.in.SMSSM

SMSSM_REFERENCES := \
		$(DIR)/SMSSM_references.tex

SMSSM_GNUPLOT := \
		$(DIR)/SMSSM_plot_rgflow.gnuplot \
		$(DIR)/SMSSM_plot_spectrum.gnuplot

SMSSM_TARBALL := \
		$(MODNAME).tar.gz

LIBSMSSM_SRC := \
		$(DIR)/SMSSM_a_muon.cpp \
		$(DIR)/SMSSM_edm.cpp \
		$(DIR)/SMSSM_FFV_form_factors.cpp \
		$(DIR)/SMSSM_l_to_lgamma.cpp \
		$(DIR)/SMSSM_effective_couplings.cpp \
		$(DIR)/SMSSM_info.cpp \
		$(DIR)/SMSSM_input_parameters.cpp \
		$(DIR)/SMSSM_mass_eigenstates.cpp \
		$(DIR)/SMSSM_observables.cpp \
		$(DIR)/SMSSM_physical.cpp \
		$(DIR)/SMSSM_slha_io.cpp \
		$(DIR)/SMSSM_soft_parameters.cpp \
		$(DIR)/SMSSM_susy_parameters.cpp \
		$(DIR)/SMSSM_utilities.cpp \
		$(DIR)/SMSSM_weinberg_angle.cpp

EXESMSSM_SRC := \
		$(DIR)/run_SMSSM.cpp \
		$(DIR)/run_cmd_line_SMSSM.cpp \
		$(DIR)/scan_SMSSM.cpp
LLSMSSM_LIB  :=
LLSMSSM_OBJ  :=
LLSMSSM_SRC  := \
		$(DIR)/SMSSM_librarylink.cpp

LLSMSSM_MMA  := \
		$(DIR)/SMSSM_librarylink.m \
		$(DIR)/run_SMSSM.m

LIBSMSSM_HDR := \
		$(DIR)/SMSSM_a_muon.hpp \
		$(DIR)/SMSSM_convergence_tester.hpp \
		$(DIR)/SMSSM_edm.hpp \
		$(DIR)/SMSSM_FFV_form_factors.hpp \
		$(DIR)/SMSSM_l_to_lgamma.hpp \
		$(DIR)/SMSSM_effective_couplings.hpp \
		$(DIR)/SMSSM_ewsb_solver.hpp \
		$(DIR)/SMSSM_ewsb_solver_interface.hpp \
		$(DIR)/SMSSM_high_scale_constraint.hpp \
		$(DIR)/SMSSM_info.hpp \
		$(DIR)/SMSSM_initial_guesser.hpp \
		$(DIR)/SMSSM_input_parameters.hpp \
		$(DIR)/SMSSM_low_scale_constraint.hpp \
		$(DIR)/SMSSM_mass_eigenstates.hpp \
		$(DIR)/SMSSM_model.hpp \
		$(DIR)/SMSSM_model_slha.hpp \
		$(DIR)/SMSSM_observables.hpp \
		$(DIR)/SMSSM_physical.hpp \
		$(DIR)/SMSSM_slha_io.hpp \
		$(DIR)/SMSSM_spectrum_generator.hpp \
		$(DIR)/SMSSM_spectrum_generator_interface.hpp \
		$(DIR)/SMSSM_soft_parameters.hpp \
		$(DIR)/SMSSM_susy_parameters.hpp \
		$(DIR)/SMSSM_susy_scale_constraint.hpp \
		$(DIR)/SMSSM_utilities.hpp \
		$(DIR)/SMSSM_weinberg_angle.hpp

LIBSMSSM_CXXQFT_HDR := \
		$(DIR)/cxx_qft/SMSSM_qft.hpp \
		$(DIR)/cxx_qft/SMSSM_fields.hpp \
		$(DIR)/cxx_qft/SMSSM_vertices.hpp \
		$(DIR)/cxx_qft/SMSSM_context_base.hpp \
		$(DIR)/cxx_qft/SMSSM_npointfunctions.hpp

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
-include $(SMSSM_SUSY_BETAS_MK)
-include $(SMSSM_SOFT_BETAS_MK)
-include $(SMSSM_CXX_QFT_VERTICES_MK)
-include $(SMSSM_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(SMSSM_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(SMSSM_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(SMSSM_CXX_QFT_VERTICES_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(SMSSM_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBSMSSM_SRC := $(sort $(LIBSMSSM_SRC))
EXESMSSM_SRC := $(sort $(EXESMSSM_SRC))

LIBSMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBSMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBSMSSM_SRC)))

EXESMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXESMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXESMSSM_SRC)))

EXESMSSM_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXESMSSM_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXESMSSM_SRC)))

LIBSMSSM_DEP := \
		$(LIBSMSSM_OBJ:.o=.d)

EXESMSSM_DEP := \
		$(EXESMSSM_OBJ:.o=.d)

LLSMSSM_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLSMSSM_SRC)))

LLSMSSM_OBJ  := $(LLSMSSM_SRC:.cpp=.o)
LLSMSSM_LIB  := $(LLSMSSM_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBSMSSM     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_SMSSM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_SMSSM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL_DIR) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBSMSSM) $(EXESMSSM_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(SMSSM_INSTALL_DIR)
		$(Q)install -d $(SMSSM_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBSMSSM_SRC) $(SMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBSMSSM_HDR) $(SMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBSMSSM_CXXQFT_HDR) $(SMSSM_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(EXESMSSM_SRC) $(SMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLSMSSM_SRC) $(SMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLSMSSM_MMA) $(SMSSM_INSTALL_DIR)
		$(Q)$(INSTALL_STRIPPED) $(SMSSM_MK) $(SMSSM_INSTALL_DIR) -m u=rw,g=r,o=r
		$(Q)install -m u=rw,g=r,o=r $(SMSSM_INCLUDE_MK) $(SMSSM_INSTALL_DIR)
ifneq ($(SMSSM_SLHA_INPUT),)
		$(Q)install -m u=rw,g=r,o=r $(SMSSM_SLHA_INPUT) $(SMSSM_INSTALL_DIR)
endif
		$(Q)install -m u=rw,g=r,o=r $(SMSSM_REFERENCES) $(SMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(SMSSM_GNUPLOT) $(SMSSM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		$(Q)-rm -f $(LIBSMSSM_DEP)
		$(Q)-rm -f $(EXESMSSM_DEP)
		$(Q)-rm -f $(LLSMSSM_DEP)

clean-$(MODNAME)-lib:
		$(Q)-rm -f $(LIBSMSSM)
		$(Q)-rm -f $(LLSMSSM_LIB)

clean-$(MODNAME)-obj:
		$(Q)-rm -f $(LIBSMSSM_OBJ)
		$(Q)-rm -f $(EXESMSSM_OBJ)
		$(Q)-rm -f $(LLSMSSM_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		$(Q)-rm -f $(LIBSMSSM_SRC)
		$(Q)-rm -f $(LIBSMSSM_HDR)
		$(Q)-rm -f $(LIBSMSSM_CXXQFT_HDR)
		$(Q)-rm -f $(EXESMSSM_SRC)
		$(Q)-rm -f $(LLSMSSM_SRC)
		$(Q)-rm -f $(LLSMSSM_MMA)
		$(Q)-rm -f $(METACODE_STAMP_SMSSM)
		$(Q)-rm -f $(SMSSM_INCLUDE_MK)
		$(Q)-rm -f $(SMSSM_SLHA_INPUT)
		$(Q)-rm -f $(SMSSM_REFERENCES)
		$(Q)-rm -f $(SMSSM_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		$(Q)-rm -f $(EXESMSSM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		$(Q)tar -czf $(SMSSM_TARBALL) \
		$(LIBSMSSM_SRC) $(LIBSMSSM_HDR) $(LIBSMSSM_CXXQFT_HDR) \
		$(EXESMSSM_SRC) \
		$(LLSMSSM_SRC) $(LLSMSSM_MMA) \
		$(SMSSM_MK) $(SMSSM_INCLUDE_MK) \
		$(SMSSM_SLHA_INPUT) $(SMSSM_REFERENCES) \
		$(SMSSM_GNUPLOT)

$(LIBSMSSM_SRC) $(LIBSMSSM_HDR) $(LIBSMSSM_CXXQFT_HDR) $(EXESMSSM_SRC) $(LLSMSSM_SRC) $(LLSMSSM_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_SMSSM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_SMSSM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_SMSSM)
		@$(MSG)
		$(Q)"$(MATH)" -run "Get[\"$<\"]; Quit[]" || (echo "Error: The code generation failed!"; exit 1)
		@touch "$(METACODE_STAMP_SMSSM)"
		@echo "Note: to regenerate SMSSM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_SMSSM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_SMSSM):
		@true
endif

$(LIBSMSSM_DEP) $(EXESMSSM_DEP) $(LLSMSSM_DEP) $(LIBSMSSM_OBJ) $(EXESMSSM_OBJ) $(LLSMSSM_OBJ) $(LLSMSSM_LIB): \
	CPPFLAGS += $(MODSMSSM_SUBMOD_INC) $(MODSMSSM_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBSMSSM_DEP) $(EXESMSSM_DEP) $(LLSMSSM_DEP) $(LIBSMSSM_OBJ) $(EXESMSSM_OBJ) $(LLSMSSM_OBJ) $(LLSMSSM_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLSMSSM_OBJ) $(LLSMSSM_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBSMSSM): $(LIBSMSSM_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBSMSSM) $(MODSMSSM_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		@$(MSG)
		$(Q)$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLSMSSM_LIB): $(LLSMSSM_OBJ) $(LIBSMSSM) $(MODSMSSM_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		@$(MSG)
		$(Q)$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS)

ALLDEP += $(LIBSMSSM_DEP) $(EXESMSSM_DEP)
ALLSRC += $(LIBSMSSM_SRC) $(EXESMSSM_SRC)
ALLLIB += $(LIBSMSSM)
ALLEXE += $(EXESMSSM_EXE)
ALLMODDEP += $(MODSMSSM_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLSMSSM_DEP)
ALLSRC += $(LLSMSSM_SRC)
ALLLL  += $(LLSMSSM_LIB)
endif
