DIR          := models/MRSSM
MODNAME      := MRSSM
SARAH_MODEL  := MRSSM
WITH_$(MODNAME) := yes
MODMRSSM_MOD := SM
MODMRSSM_DEP := $(patsubst %,model_specific/%,$(MODMRSSM_MOD))
MODMRSSM_INC := $(patsubst %,-Imodel_specific/%,$(MODMRSSM_MOD))
MODMRSSM_LIB := $(foreach M,$(MODMRSSM_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

MODMRSSM_SUBMOD  := $(DIR)/cxx_qft
MODMRSSM_SUBMOD_INC := $(patsubst %,-I%,$(MODMRSSM_SUBMOD))

MRSSM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)
MRSSM_INSTALL_CXXQFT_DIR := \
		$(MRSSM_INSTALL_DIR)/cxx_qft

MRSSM_MK     := \
		$(DIR)/module.mk

MRSSM_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

MRSSM_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

MRSSM_CXX_QFT_VERTICES_MK := \
		$(DIR)/cxx_qft/vertices.mk

MRSSM_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

MRSSM_INCLUDE_MK := \
		$(MRSSM_SUSY_BETAS_MK) \
		$(MRSSM_SOFT_BETAS_MK) \
		$(MRSSM_CXX_QFT_VERTICES_MK)

MRSSM_SLHA_INPUT := \
		$(DIR)/LesHouches.in.MRSSM_generated \
		$(DIR)/LesHouches.in.MRSSM

MRSSM_REFERENCES := \
		$(DIR)/MRSSM_references.tex

MRSSM_GNUPLOT := \
		$(DIR)/MRSSM_plot_rgflow.gnuplot \
		$(DIR)/MRSSM_plot_spectrum.gnuplot

MRSSM_TARBALL := \
		$(MODNAME).tar.gz

LIBMRSSM_SRC := \
		$(DIR)/MRSSM_a_muon.cpp \
		$(DIR)/MRSSM_edm.cpp \
		$(DIR)/MRSSM_FFV_form_factors.cpp \
		$(DIR)/MRSSM_l_to_lgamma.cpp \
		$(DIR)/MRSSM_effective_couplings.cpp \
		$(DIR)/MRSSM_info.cpp \
		$(DIR)/MRSSM_input_parameters.cpp \
		$(DIR)/MRSSM_mass_eigenstates.cpp \
		$(DIR)/MRSSM_observables.cpp \
		$(DIR)/MRSSM_physical.cpp \
		$(DIR)/MRSSM_slha_io.cpp \
		$(DIR)/MRSSM_soft_parameters.cpp \
		$(DIR)/MRSSM_susy_parameters.cpp \
		$(DIR)/MRSSM_utilities.cpp \
		$(DIR)/MRSSM_weinberg_angle.cpp

EXEMRSSM_SRC := \
		$(DIR)/run_MRSSM.cpp \
		$(DIR)/run_cmd_line_MRSSM.cpp \
		$(DIR)/scan_MRSSM.cpp
LLMRSSM_LIB  :=
LLMRSSM_OBJ  :=
LLMRSSM_SRC  := \
		$(DIR)/MRSSM_librarylink.cpp

LLMRSSM_MMA  := \
		$(DIR)/MRSSM_librarylink.m \
		$(DIR)/run_MRSSM.m

LIBMRSSM_HDR := \
		$(DIR)/MRSSM_a_muon.hpp \
		$(DIR)/MRSSM_convergence_tester.hpp \
		$(DIR)/MRSSM_edm.hpp \
		$(DIR)/MRSSM_FFV_form_factors.hpp \
		$(DIR)/MRSSM_l_to_lgamma.hpp \
		$(DIR)/MRSSM_effective_couplings.hpp \
		$(DIR)/MRSSM_ewsb_solver.hpp \
		$(DIR)/MRSSM_ewsb_solver_interface.hpp \
		$(DIR)/MRSSM_high_scale_constraint.hpp \
		$(DIR)/MRSSM_info.hpp \
		$(DIR)/MRSSM_initial_guesser.hpp \
		$(DIR)/MRSSM_input_parameters.hpp \
		$(DIR)/MRSSM_low_scale_constraint.hpp \
		$(DIR)/MRSSM_mass_eigenstates.hpp \
		$(DIR)/MRSSM_model.hpp \
		$(DIR)/MRSSM_model_slha.hpp \
		$(DIR)/MRSSM_observables.hpp \
		$(DIR)/MRSSM_physical.hpp \
		$(DIR)/MRSSM_slha_io.hpp \
		$(DIR)/MRSSM_spectrum_generator.hpp \
		$(DIR)/MRSSM_spectrum_generator_interface.hpp \
		$(DIR)/MRSSM_soft_parameters.hpp \
		$(DIR)/MRSSM_susy_parameters.hpp \
		$(DIR)/MRSSM_susy_scale_constraint.hpp \
		$(DIR)/MRSSM_utilities.hpp \
		$(DIR)/MRSSM_weinberg_angle.hpp

LIBMRSSM_CXXQFT_HDR := \
		$(DIR)/cxx_qft/MRSSM_qft.hpp \
		$(DIR)/cxx_qft/MRSSM_fields.hpp \
		$(DIR)/cxx_qft/MRSSM_vertices.hpp \
		$(DIR)/cxx_qft/MRSSM_context_base.hpp \
		$(DIR)/cxx_qft/MRSSM_npointfunctions.hpp

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
-include $(MRSSM_SUSY_BETAS_MK)
-include $(MRSSM_SOFT_BETAS_MK)
-include $(MRSSM_CXX_QFT_VERTICES_MK)
-include $(MRSSM_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(MRSSM_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(MRSSM_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(MRSSM_CXX_QFT_VERTICES_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(MRSSM_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBMRSSM_SRC := $(sort $(LIBMRSSM_SRC))
EXEMRSSM_SRC := $(sort $(EXEMRSSM_SRC))

LIBMRSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBMRSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBMRSSM_SRC)))

EXEMRSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEMRSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEMRSSM_SRC)))

EXEMRSSM_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXEMRSSM_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXEMRSSM_SRC)))

LIBMRSSM_DEP := \
		$(LIBMRSSM_OBJ:.o=.d)

EXEMRSSM_DEP := \
		$(EXEMRSSM_OBJ:.o=.d)

LLMRSSM_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLMRSSM_SRC)))

LLMRSSM_OBJ  := $(LLMRSSM_SRC:.cpp=.o)
LLMRSSM_LIB  := $(LLMRSSM_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBMRSSM     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_MRSSM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_MRSSM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL_DIR) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBMRSSM) $(EXEMRSSM_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(MRSSM_INSTALL_DIR)
		$(Q)install -d $(MRSSM_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBMRSSM_SRC) $(MRSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBMRSSM_HDR) $(MRSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBMRSSM_CXXQFT_HDR) $(MRSSM_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(EXEMRSSM_SRC) $(MRSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLMRSSM_SRC) $(MRSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLMRSSM_MMA) $(MRSSM_INSTALL_DIR)
		$(Q)$(INSTALL_STRIPPED) $(MRSSM_MK) $(MRSSM_INSTALL_DIR) -m u=rw,g=r,o=r
		$(Q)install -m u=rw,g=r,o=r $(MRSSM_INCLUDE_MK) $(MRSSM_INSTALL_DIR)
ifneq ($(MRSSM_SLHA_INPUT),)
		$(Q)install -m u=rw,g=r,o=r $(MRSSM_SLHA_INPUT) $(MRSSM_INSTALL_DIR)
endif
		$(Q)install -m u=rw,g=r,o=r $(MRSSM_REFERENCES) $(MRSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(MRSSM_GNUPLOT) $(MRSSM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		$(Q)-rm -f $(LIBMRSSM_DEP)
		$(Q)-rm -f $(EXEMRSSM_DEP)
		$(Q)-rm -f $(LLMRSSM_DEP)

clean-$(MODNAME)-lib:
		$(Q)-rm -f $(LIBMRSSM)
		$(Q)-rm -f $(LLMRSSM_LIB)

clean-$(MODNAME)-obj:
		$(Q)-rm -f $(LIBMRSSM_OBJ)
		$(Q)-rm -f $(EXEMRSSM_OBJ)
		$(Q)-rm -f $(LLMRSSM_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		$(Q)-rm -f $(LIBMRSSM_SRC)
		$(Q)-rm -f $(LIBMRSSM_HDR)
		$(Q)-rm -f $(LIBMRSSM_CXXQFT_HDR)
		$(Q)-rm -f $(EXEMRSSM_SRC)
		$(Q)-rm -f $(LLMRSSM_SRC)
		$(Q)-rm -f $(LLMRSSM_MMA)
		$(Q)-rm -f $(METACODE_STAMP_MRSSM)
		$(Q)-rm -f $(MRSSM_INCLUDE_MK)
		$(Q)-rm -f $(MRSSM_SLHA_INPUT)
		$(Q)-rm -f $(MRSSM_REFERENCES)
		$(Q)-rm -f $(MRSSM_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		$(Q)-rm -f $(EXEMRSSM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		$(Q)tar -czf $(MRSSM_TARBALL) \
		$(LIBMRSSM_SRC) $(LIBMRSSM_HDR) $(LIBMRSSM_CXXQFT_HDR) \
		$(EXEMRSSM_SRC) \
		$(LLMRSSM_SRC) $(LLMRSSM_MMA) \
		$(MRSSM_MK) $(MRSSM_INCLUDE_MK) \
		$(MRSSM_SLHA_INPUT) $(MRSSM_REFERENCES) \
		$(MRSSM_GNUPLOT)

$(LIBMRSSM_SRC) $(LIBMRSSM_HDR) $(LIBMRSSM_CXXQFT_HDR) $(EXEMRSSM_SRC) $(LLMRSSM_SRC) $(LLMRSSM_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_MRSSM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_MRSSM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_MRSSM)
		@$(MSG)
		$(Q)"$(MATH)" -run "Get[\"$<\"]; Quit[]" || (echo "Error: The code generation failed!"; exit 1)
		@touch "$(METACODE_STAMP_MRSSM)"
		@echo "Note: to regenerate MRSSM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_MRSSM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_MRSSM):
		@true
endif

$(LIBMRSSM_DEP) $(EXEMRSSM_DEP) $(LLMRSSM_DEP) $(LIBMRSSM_OBJ) $(EXEMRSSM_OBJ) $(LLMRSSM_OBJ) $(LLMRSSM_LIB): \
	CPPFLAGS += $(MODMRSSM_SUBMOD_INC) $(MODMRSSM_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBMRSSM_DEP) $(EXEMRSSM_DEP) $(LLMRSSM_DEP) $(LIBMRSSM_OBJ) $(EXEMRSSM_OBJ) $(LLMRSSM_OBJ) $(LLMRSSM_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLMRSSM_OBJ) $(LLMRSSM_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBMRSSM): $(LIBMRSSM_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBMRSSM) $(MODMRSSM_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		@$(MSG)
		$(Q)$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLMRSSM_LIB): $(LLMRSSM_OBJ) $(LIBMRSSM) $(MODMRSSM_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		@$(MSG)
		$(Q)$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS)

ALLDEP += $(LIBMRSSM_DEP) $(EXEMRSSM_DEP)
ALLSRC += $(LIBMRSSM_SRC) $(EXEMRSSM_SRC)
ALLLIB += $(LIBMRSSM)
ALLEXE += $(EXEMRSSM_EXE)
ALLMODDEP += $(MODMRSSM_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLMRSSM_DEP)
ALLSRC += $(LLMRSSM_SRC)
ALLLL  += $(LLMRSSM_LIB)
endif
