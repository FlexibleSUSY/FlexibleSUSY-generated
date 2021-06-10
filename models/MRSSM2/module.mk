DIR          := models/MRSSM2
MODNAME      := MRSSM2
SARAH_MODEL  := MRSSM
WITH_$(MODNAME) := yes
MODMRSSM2_MOD := SM
MODMRSSM2_DEP := $(patsubst %,model_specific/%,$(MODMRSSM2_MOD))
MODMRSSM2_INC := $(patsubst %,-Imodel_specific/%,$(MODMRSSM2_MOD))
MODMRSSM2_LIB := $(foreach M,$(MODMRSSM2_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

MODMRSSM2_SUBMOD  := $(DIR)/cxx_qft
MODMRSSM2_SUBMOD_INC := $(patsubst %,-I%,$(MODMRSSM2_SUBMOD))

MRSSM2_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)
MRSSM2_INSTALL_CXXQFT_DIR := \
		$(MRSSM2_INSTALL_DIR)/cxx_qft

MRSSM2_MK     := \
		$(DIR)/module.mk

MRSSM2_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

MRSSM2_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

MRSSM2_CXXQFT_VERTICES_MK := \
		$(DIR)/cxx_qft/vertices.mk

-include $(MRSSM2_CXXQFT_VERTICES_MK)
LIBMRSSM2_CXXQFT_VERTICES_SRC ?= ''

MRSSM2_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

MRSSM2_FlexibleDecay_MK := \
		$(DIR)/decays/FlexibleDecay.mk

MRSSM2_INCLUDE_MK := \
		$(MRSSM2_SUSY_BETAS_MK) \
		$(MRSSM2_SOFT_BETAS_MK)

MRSSM2_SLHA_INPUT := \
		$(DIR)/LesHouches.in.MRSSM2_generated \
		$(DIR)/LesHouches.in.MRSSM2

MRSSM2_REFERENCES := \
		$(DIR)/MRSSM2_references.tex

MRSSM2_GNUPLOT := \
		$(DIR)/MRSSM2_plot_rgflow.gnuplot \
		$(DIR)/MRSSM2_plot_spectrum.gnuplot

MRSSM2_TARBALL := \
		$(MODNAME).tar.gz

LIBMRSSM2_SRC := \
		$(DIR)/MRSSM2_a_muon.cpp \
		$(DIR)/MRSSM2_edm.cpp \
		$(DIR)/MRSSM2_FFV_form_factors.cpp \
		$(DIR)/MRSSM2_f_to_f_conversion.cpp \
		$(DIR)/MRSSM2_l_to_lgamma.cpp \
		$(DIR)/MRSSM2_b_to_s_gamma.cpp \
		$(DIR)/MRSSM2_effective_couplings.cpp \
		$(DIR)/MRSSM2_info.cpp \
		$(DIR)/MRSSM2_input_parameters.cpp \
		$(DIR)/MRSSM2_mass_eigenstates.cpp \
		$(DIR)/MRSSM2_mass_eigenstates_decoupling_scheme.cpp \
		$(DIR)/MRSSM2_model_slha.cpp \
		$(DIR)/MRSSM2_observables.cpp \
		$(DIR)/MRSSM2_physical.cpp \
		$(DIR)/MRSSM2_slha_io.cpp \
		$(DIR)/MRSSM2_soft_parameters.cpp \
		$(DIR)/MRSSM2_susy_parameters.cpp \
		$(DIR)/MRSSM2_utilities.cpp \
		$(DIR)/MRSSM2_weinberg_angle.cpp

LIBMRSSM2_SRC += $(LIBMRSSM2_CXXQFT_VERTICES_SRC)

EXEMRSSM2_SRC := \
		$(DIR)/run_MRSSM2.cpp \
		$(DIR)/run_cmd_line_MRSSM2.cpp \
		$(DIR)/scan_MRSSM2.cpp
LLMRSSM2_LIB  :=
LLMRSSM2_OBJ  :=
LLMRSSM2_SRC  := \
		$(DIR)/MRSSM2_librarylink.cpp

LLMRSSM2_MMA  := \
		$(DIR)/MRSSM2_librarylink.m \
		$(DIR)/run_MRSSM2.m

LIBMRSSM2_HDR := \
		$(DIR)/MRSSM2_a_muon.hpp \
		$(DIR)/MRSSM2_convergence_tester.hpp \
		$(DIR)/MRSSM2_edm.hpp \
		$(DIR)/MRSSM2_FFV_form_factors.hpp \
		$(DIR)/MRSSM2_f_to_f_conversion.hpp \
		$(DIR)/MRSSM2_l_to_lgamma.hpp \
		$(DIR)/MRSSM2_b_to_s_gamma.hpp \
		$(DIR)/MRSSM2_effective_couplings.hpp \
		$(DIR)/MRSSM2_ewsb_solver.hpp \
		$(DIR)/MRSSM2_ewsb_solver_interface.hpp \
		$(DIR)/MRSSM2_high_scale_constraint.hpp \
		$(DIR)/MRSSM2_info.hpp \
		$(DIR)/MRSSM2_initial_guesser.hpp \
		$(DIR)/MRSSM2_input_parameters.hpp \
		$(DIR)/MRSSM2_low_scale_constraint.hpp \
		$(DIR)/MRSSM2_mass_eigenstates.hpp \
		$(DIR)/MRSSM2_mass_eigenstates_interface.hpp \
		$(DIR)/MRSSM2_mass_eigenstates_decoupling_scheme.hpp \
		$(DIR)/MRSSM2_model.hpp \
		$(DIR)/MRSSM2_model_slha.hpp \
		$(DIR)/MRSSM2_observables.hpp \
		$(DIR)/MRSSM2_physical.hpp \
		$(DIR)/MRSSM2_slha_io.hpp \
		$(DIR)/MRSSM2_spectrum_generator.hpp \
		$(DIR)/MRSSM2_spectrum_generator_interface.hpp \
		$(DIR)/MRSSM2_soft_parameters.hpp \
		$(DIR)/MRSSM2_susy_parameters.hpp \
		$(DIR)/MRSSM2_susy_scale_constraint.hpp \
		$(DIR)/MRSSM2_utilities.hpp \
		$(DIR)/MRSSM2_weinberg_angle.hpp

LIBMRSSM2_CXXQFT_HDR := \
		$(DIR)/cxx_qft/MRSSM2_qft.hpp \
		$(DIR)/cxx_qft/MRSSM2_fields.hpp \
		$(DIR)/cxx_qft/MRSSM2_vertices.hpp \
		$(DIR)/cxx_qft/MRSSM2_context_base.hpp \
		$(DIR)/cxx_qft/MRSSM2_npointfunctions_wilsoncoeffs.hpp

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
-include $(MRSSM2_SUSY_BETAS_MK)
-include $(MRSSM2_SOFT_BETAS_MK)
-include $(MRSSM2_FlexibleDecay_MK)
-include $(MRSSM2_CXXQFT_VERTICES_MK)
-include $(MRSSM2_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(MRSSM2_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(MRSSM2_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@

$(MRSSM2_FlexibleDecay_MK): run-metacode-$(MODNAME)
$(MRSSM2_CXXQFT_VERTICES_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(MRSSM2_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBMRSSM2_SRC := $(sort $(LIBMRSSM2_SRC))
EXEMRSSM2_SRC := $(sort $(EXEMRSSM2_SRC))

LIBMRSSM2_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBMRSSM2_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBMRSSM2_SRC)))

EXEMRSSM2_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEMRSSM2_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEMRSSM2_SRC)))

EXEMRSSM2_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXEMRSSM2_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXEMRSSM2_SRC)))

LIBMRSSM2_DEP := \
		$(LIBMRSSM2_OBJ:.o=.d)

EXEMRSSM2_DEP := \
		$(EXEMRSSM2_OBJ:.o=.d)

LLMRSSM2_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLMRSSM2_SRC)))

LLMRSSM2_OBJ  := $(LLMRSSM2_SRC:.cpp=.o)
LLMRSSM2_LIB  := $(LLMRSSM2_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBMRSSM2     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_MRSSM2 := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_MRSSM2 := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL_DIR) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBMRSSM2) $(EXEMRSSM2_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(MRSSM2_INSTALL_DIR)
		$(Q)install -d $(MRSSM2_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBMRSSM2_SRC) $(MRSSM2_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBMRSSM2_CXXQFT_VERTICES_SRC) $(MRSSM2_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBMRSSM2_HDR) $(MRSSM2_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBMRSSM2_CXXQFT_HDR) $(MRSSM2_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(EXEMRSSM2_SRC) $(MRSSM2_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLMRSSM2_SRC) $(MRSSM2_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLMRSSM2_MMA) $(MRSSM2_INSTALL_DIR)
		$(Q)$(INSTALL_STRIPPED) $(MRSSM2_MK) $(MRSSM2_INSTALL_DIR) -m u=rw,g=r,o=r
		$(Q)install -m u=rw,g=r,o=r $(MRSSM2_INCLUDE_MK) $(MRSSM2_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(MRSSM2_CXXQFT_VERTICES_MK) $(MRSSM2_INSTALL_CXXQFT_DIR)

ifneq ($(MRSSM2_SLHA_INPUT),)
		$(Q)install -m u=rw,g=r,o=r $(MRSSM2_SLHA_INPUT) $(MRSSM2_INSTALL_DIR)
endif
		$(Q)install -m u=rw,g=r,o=r $(MRSSM2_REFERENCES) $(MRSSM2_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(MRSSM2_GNUPLOT) $(MRSSM2_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		$(Q)-rm -f $(LIBMRSSM2_DEP)
		$(Q)-rm -f $(EXEMRSSM2_DEP)
		$(Q)-rm -f $(LLMRSSM2_DEP)

clean-$(MODNAME)-lib:
		$(Q)-rm -f $(LIBMRSSM2)
		$(Q)-rm -f $(LLMRSSM2_LIB)

clean-$(MODNAME)-obj:
		$(Q)-rm -f $(LIBMRSSM2_OBJ)
		$(Q)-rm -f $(EXEMRSSM2_OBJ)
		$(Q)-rm -f $(LLMRSSM2_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		$(Q)-rm -f $(LIBMRSSM2_SRC)
		$(Q)-rm -f $(LIBMRSSM2_HDR)
		$(Q)-rm -f $(LIBMRSSM2_CXXQFT_HDR)
		$(Q)-rm -f $(EXEMRSSM2_SRC)
		$(Q)-rm -f $(LLMRSSM2_SRC)
		$(Q)-rm -f $(LLMRSSM2_MMA)
		$(Q)-rm -f $(METACODE_STAMP_MRSSM2)
		$(Q)-rm -f $(MRSSM2_INCLUDE_MK)
		$(Q)-rm -f $(MRSSM2_CXXQFT_VERTICES_MK)
		$(Q)-rm -f $(MRSSM2_SLHA_INPUT)
		$(Q)-rm -f $(MRSSM2_REFERENCES)
		$(Q)-rm -f $(MRSSM2_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		$(Q)-rm -f $(EXEMRSSM2_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		$(Q)tar -czf $(MRSSM2_TARBALL) \
		$(LIBMRSSM2_SRC) $(LIBMRSSM2_HDR) $(LIBMRSSM2_CXXQFT_HDR) \
		$(EXEMRSSM2_SRC) \
		$(LLMRSSM2_SRC) $(LLMRSSM2_MMA) \
		$(MRSSM2_MK) $(MRSSM2_INCLUDE_MK) $(MRSSM2_CXXQFT_VERTICES_MK) \
		$(MRSSM2_SLHA_INPUT) $(MRSSM2_REFERENCES) \
		$(MRSSM2_GNUPLOT) \
		$(MRSSM2_FlexibleDecay_MK)

$(LIBMRSSM2_SRC) $(LIBMRSSM2_HDR) $(LIBMRSSM2_CXXQFT_HDR) $(EXEMRSSM2_SRC) $(LLMRSSM2_SRC) $(LLMRSSM2_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_MRSSM2)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_MRSSM2): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_MRSSM2)
		@$(MSG)
		$(Q)printf "%s" "Get[\"$<\"]; Quit[]" | "$(MATH)" || (echo "Error: The code generation failed!"; exit 1)
		@touch "$(METACODE_STAMP_MRSSM2)"
		@echo "Note: to regenerate MRSSM2 source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_MRSSM2)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_MRSSM2):
		@true
endif

$(LIBMRSSM2_DEP) $(EXEMRSSM2_DEP) $(LLMRSSM2_DEP) $(LIBMRSSM2_OBJ) $(EXEMRSSM2_OBJ) $(LLMRSSM2_OBJ) $(LLMRSSM2_LIB): \
	CPPFLAGS += $(MODMRSSM2_SUBMOD_INC) $(MODMRSSM2_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(GM2CALCFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBMRSSM2_DEP) $(EXEMRSSM2_DEP) $(LLMRSSM2_DEP) $(LIBMRSSM2_OBJ) $(EXEMRSSM2_OBJ) $(LLMRSSM2_OBJ) $(LLMRSSM2_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLMRSSM2_OBJ) $(LLMRSSM2_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBMRSSM2): $(LIBMRSSM2_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBMRSSM2) $(MODMRSSM2_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(FUTILIBS)
		@$(MSG)
		$(Q)$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^) $(filter -%,$(LOOPFUNCLIBS)) $(GM2CALCLIBS) $(HIMALAYALIBS) $(GSLLIBS) $(SQLITELIBS) $(TSILLIBS) $(FLIBS) $(THREADLIBS) $(LDLIBS) $(FUTILIBS)

$(LLMRSSM2_LIB): $(LLMRSSM2_OBJ) $(LIBMRSSM2) $(MODMRSSM2_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(FUTILIBS)
		@$(MSG)
		$(Q)$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^) $(filter -%,$(LOOPFUNCLIBS)) $(GM2CALCLIBS) $(HIMALAYALIBS) $(TSILLIBS) $(GSLLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS) $(FUTILIBS) $(FLIBS)

ALLDEP += $(LIBMRSSM2_DEP) $(EXEMRSSM2_DEP)
ALLSRC += $(LIBMRSSM2_SRC) $(EXEMRSSM2_SRC)
ALLLIB += $(LIBMRSSM2)
ALLEXE += $(EXEMRSSM2_EXE)
ALLMODDEP += $(MODMRSSM2_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLMRSSM2_DEP)
ALLSRC += $(LLMRSSM2_SRC)
ALLLL  += $(LLMRSSM2_LIB)
endif
