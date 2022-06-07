DIR          := models/NUTNMSSM
MODNAME      := NUTNMSSM
SARAH_MODEL  := NMSSM
WITH_$(MODNAME) := yes
MODNUTNMSSM_MOD := SM MSSM_higgs NMSSM_higgs
MODNUTNMSSM_DEP := $(patsubst %,model_specific/%,$(MODNUTNMSSM_MOD))
MODNUTNMSSM_INC := $(patsubst %,-Imodel_specific/%,$(MODNUTNMSSM_MOD))
MODNUTNMSSM_LIB := $(foreach M,$(MODNUTNMSSM_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

MODNUTNMSSM_SUBMOD  := $(DIR)/cxx_qft
MODNUTNMSSM_SUBMOD_INC := $(patsubst %,-I%,$(MODNUTNMSSM_SUBMOD))

NUTNMSSM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)
NUTNMSSM_INSTALL_CXXQFT_DIR := \
		$(NUTNMSSM_INSTALL_DIR)/cxx_qft

NUTNMSSM_MK     := \
		$(DIR)/module.mk

NUTNMSSM_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

NUTNMSSM_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

NUTNMSSM_CXXQFT_VERTICES_MK := \
		$(DIR)/cxx_qft/vertices.mk

-include $(NUTNMSSM_CXXQFT_VERTICES_MK)
LIBNUTNMSSM_CXXQFT_VERTICES_SRC ?= ''

NUTNMSSM_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

NUTNMSSM_FlexibleDecay_MK := \
		$(DIR)/decays/FlexibleDecay.mk

NUTNMSSM_INCLUDE_MK := \
		$(NUTNMSSM_SUSY_BETAS_MK) \
		$(NUTNMSSM_SOFT_BETAS_MK)

NUTNMSSM_SLHA_INPUT := \
		$(DIR)/LesHouches.in.NUTNMSSM_generated \
		$(DIR)/LesHouches.in.NUTNMSSM \
		$(DIR)/LesHouches.in.NUTNMSSM_GTP1 \
		$(DIR)/LesHouches.in.NUTNMSSM_1308.1333_BP3 \
		$(DIR)/LesHouches.in.NUTNMSSM_1308.1333_BP2 \
		$(DIR)/LesHouches.in.NUTNMSSM_1308.1333_BP1 \
		$(DIR)/LesHouches.in.NUTNMSSM_GTP2

NUTNMSSM_REFERENCES := \
		$(DIR)/NUTNMSSM_references.tex

NUTNMSSM_GNUPLOT := \
		$(DIR)/NUTNMSSM_plot_rgflow.gnuplot \
		$(DIR)/NUTNMSSM_plot_spectrum.gnuplot

NUTNMSSM_TARBALL := \
		$(MODNAME).tar.gz

LIBNUTNMSSM_SRC := \
		$(DIR)/NUTNMSSM_a_muon.cpp \
		$(DIR)/NUTNMSSM_edm.cpp \
		$(DIR)/NUTNMSSM_FFV_form_factors.cpp \
		$(DIR)/NUTNMSSM_f_to_f_conversion.cpp \
		$(DIR)/NUTNMSSM_l_to_lgamma.cpp \
		$(DIR)/NUTNMSSM_b_to_s_gamma.cpp \
		$(DIR)/NUTNMSSM_info.cpp \
		$(DIR)/NUTNMSSM_input_parameters.cpp \
		$(DIR)/NUTNMSSM_mass_eigenstates.cpp \
		$(DIR)/NUTNMSSM_mass_eigenstates_decoupling_scheme.cpp \
		$(DIR)/NUTNMSSM_model_slha.cpp \
		$(DIR)/NUTNMSSM_observables.cpp \
		$(DIR)/NUTNMSSM_physical.cpp \
		$(DIR)/NUTNMSSM_slha_io.cpp \
		$(DIR)/NUTNMSSM_soft_parameters.cpp \
		$(DIR)/NUTNMSSM_susy_parameters.cpp \
		$(DIR)/NUTNMSSM_utilities.cpp \
		$(DIR)/NUTNMSSM_weinberg_angle.cpp

LIBNUTNMSSM_SRC += $(LIBNUTNMSSM_CXXQFT_VERTICES_SRC)

EXENUTNMSSM_SRC := \
		$(DIR)/run_NUTNMSSM.cpp \
		$(DIR)/run_cmd_line_NUTNMSSM.cpp \
		$(DIR)/scan_NUTNMSSM.cpp
LLNUTNMSSM_LIB  :=
LLNUTNMSSM_OBJ  :=
LLNUTNMSSM_SRC  := \
		$(DIR)/NUTNMSSM_librarylink.cpp

LLNUTNMSSM_MMA  := \
		$(DIR)/NUTNMSSM_librarylink.m \
		$(DIR)/run_NUTNMSSM.m

LIBNUTNMSSM_HDR := \
		$(DIR)/NUTNMSSM_a_muon.hpp \
		$(DIR)/NUTNMSSM_convergence_tester.hpp \
		$(DIR)/NUTNMSSM_edm.hpp \
		$(DIR)/NUTNMSSM_FFV_form_factors.hpp \
		$(DIR)/NUTNMSSM_f_to_f_conversion.hpp \
		$(DIR)/NUTNMSSM_l_to_lgamma.hpp \
		$(DIR)/NUTNMSSM_b_to_s_gamma.hpp \
		$(DIR)/NUTNMSSM_ewsb_solver.hpp \
		$(DIR)/NUTNMSSM_ewsb_solver_interface.hpp \
		$(DIR)/NUTNMSSM_high_scale_constraint.hpp \
		$(DIR)/NUTNMSSM_info.hpp \
		$(DIR)/NUTNMSSM_initial_guesser.hpp \
		$(DIR)/NUTNMSSM_input_parameters.hpp \
		$(DIR)/NUTNMSSM_low_scale_constraint.hpp \
		$(DIR)/NUTNMSSM_mass_eigenstates.hpp \
		$(DIR)/NUTNMSSM_mass_eigenstates_interface.hpp \
		$(DIR)/NUTNMSSM_mass_eigenstates_decoupling_scheme.hpp \
		$(DIR)/NUTNMSSM_model.hpp \
		$(DIR)/NUTNMSSM_model_slha.hpp \
		$(DIR)/NUTNMSSM_observables.hpp \
		$(DIR)/NUTNMSSM_physical.hpp \
		$(DIR)/NUTNMSSM_slha_io.hpp \
		$(DIR)/NUTNMSSM_spectrum_generator.hpp \
		$(DIR)/NUTNMSSM_spectrum_generator_interface.hpp \
		$(DIR)/NUTNMSSM_soft_parameters.hpp \
		$(DIR)/NUTNMSSM_susy_parameters.hpp \
		$(DIR)/NUTNMSSM_susy_scale_constraint.hpp \
		$(DIR)/NUTNMSSM_utilities.hpp \
		$(DIR)/NUTNMSSM_weinberg_angle.hpp

LIBNUTNMSSM_CXXQFT_HDR := \
		$(DIR)/cxx_qft/NUTNMSSM_qft.hpp \
		$(DIR)/cxx_qft/NUTNMSSM_fields.hpp \
		$(DIR)/cxx_qft/NUTNMSSM_vertices.hpp \
		$(DIR)/cxx_qft/NUTNMSSM_context_base.hpp \
		$(DIR)/cxx_qft/NUTNMSSM_npointfunctions_wilsoncoeffs.hpp

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
-include $(NUTNMSSM_SUSY_BETAS_MK)
-include $(NUTNMSSM_SOFT_BETAS_MK)
-include $(NUTNMSSM_FlexibleDecay_MK)
-include $(NUTNMSSM_CXXQFT_VERTICES_MK)
-include $(NUTNMSSM_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(NUTNMSSM_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(NUTNMSSM_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@

$(NUTNMSSM_FlexibleDecay_MK): run-metacode-$(MODNAME)
$(NUTNMSSM_CXXQFT_VERTICES_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(NUTNMSSM_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBNUTNMSSM_SRC := $(sort $(LIBNUTNMSSM_SRC))
EXENUTNMSSM_SRC := $(sort $(EXENUTNMSSM_SRC))

LIBNUTNMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBNUTNMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBNUTNMSSM_SRC)))

EXENUTNMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXENUTNMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXENUTNMSSM_SRC)))

EXENUTNMSSM_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXENUTNMSSM_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXENUTNMSSM_SRC)))

LIBNUTNMSSM_DEP := \
		$(LIBNUTNMSSM_OBJ:.o=.d)

EXENUTNMSSM_DEP := \
		$(EXENUTNMSSM_OBJ:.o=.d)

LLNUTNMSSM_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLNUTNMSSM_SRC)))

LLNUTNMSSM_OBJ  := $(LLNUTNMSSM_SRC:.cpp=.o)
LLNUTNMSSM_LIB  := $(LLNUTNMSSM_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBNUTNMSSM     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_NUTNMSSM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_NUTNMSSM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL_DIR) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBNUTNMSSM) $(EXENUTNMSSM_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(NUTNMSSM_INSTALL_DIR)
		$(Q)install -d $(NUTNMSSM_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBNUTNMSSM_SRC) $(NUTNMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBNUTNMSSM_CXXQFT_VERTICES_SRC) $(NUTNMSSM_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBNUTNMSSM_HDR) $(NUTNMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBNUTNMSSM_CXXQFT_HDR) $(NUTNMSSM_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(EXENUTNMSSM_SRC) $(NUTNMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLNUTNMSSM_SRC) $(NUTNMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLNUTNMSSM_MMA) $(NUTNMSSM_INSTALL_DIR)
		$(Q)$(INSTALL_STRIPPED) $(NUTNMSSM_MK) $(NUTNMSSM_INSTALL_DIR) -m u=rw,g=r,o=r
		$(Q)install -m u=rw,g=r,o=r $(NUTNMSSM_INCLUDE_MK) $(NUTNMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(NUTNMSSM_CXXQFT_VERTICES_MK) $(NUTNMSSM_INSTALL_CXXQFT_DIR)

ifneq ($(NUTNMSSM_SLHA_INPUT),)
		$(Q)install -m u=rw,g=r,o=r $(NUTNMSSM_SLHA_INPUT) $(NUTNMSSM_INSTALL_DIR)
endif
		$(Q)install -m u=rw,g=r,o=r $(NUTNMSSM_REFERENCES) $(NUTNMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(NUTNMSSM_GNUPLOT) $(NUTNMSSM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		$(Q)-rm -f $(LIBNUTNMSSM_DEP)
		$(Q)-rm -f $(EXENUTNMSSM_DEP)
		$(Q)-rm -f $(LLNUTNMSSM_DEP)

clean-$(MODNAME)-lib:
		$(Q)-rm -f $(LIBNUTNMSSM)
		$(Q)-rm -f $(LLNUTNMSSM_LIB)

clean-$(MODNAME)-obj:
		$(Q)-rm -f $(LIBNUTNMSSM_OBJ)
		$(Q)-rm -f $(EXENUTNMSSM_OBJ)
		$(Q)-rm -f $(LLNUTNMSSM_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		$(Q)-rm -f $(LIBNUTNMSSM_SRC)
		$(Q)-rm -f $(LIBNUTNMSSM_HDR)
		$(Q)-rm -f $(LIBNUTNMSSM_CXXQFT_HDR)
		$(Q)-rm -f $(EXENUTNMSSM_SRC)
		$(Q)-rm -f $(LLNUTNMSSM_SRC)
		$(Q)-rm -f $(LLNUTNMSSM_MMA)
		$(Q)-rm -f $(METACODE_STAMP_NUTNMSSM)
		$(Q)-rm -f $(NUTNMSSM_INCLUDE_MK)
		$(Q)-rm -f $(NUTNMSSM_CXXQFT_VERTICES_MK)
		$(Q)-rm -f $(NUTNMSSM_SLHA_INPUT)
		$(Q)-rm -f $(NUTNMSSM_REFERENCES)
		$(Q)-rm -f $(NUTNMSSM_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		$(Q)-rm -f $(EXENUTNMSSM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		$(Q)tar -czf $(NUTNMSSM_TARBALL) \
		$(LIBNUTNMSSM_SRC) $(LIBNUTNMSSM_HDR) $(LIBNUTNMSSM_CXXQFT_HDR) \
		$(EXENUTNMSSM_SRC) \
		$(LLNUTNMSSM_SRC) $(LLNUTNMSSM_MMA) \
		$(NUTNMSSM_MK) $(NUTNMSSM_INCLUDE_MK) $(NUTNMSSM_CXXQFT_VERTICES_MK) \
		$(NUTNMSSM_SLHA_INPUT) $(NUTNMSSM_REFERENCES) \
		$(NUTNMSSM_GNUPLOT) \
		$(NUTNMSSM_FlexibleDecay_MK)

$(LIBNUTNMSSM_SRC) $(LIBNUTNMSSM_HDR) $(LIBNUTNMSSM_CXXQFT_HDR) $(EXENUTNMSSM_SRC) $(LLNUTNMSSM_SRC) $(LLNUTNMSSM_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_NUTNMSSM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_NUTNMSSM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_NUTNMSSM)
		@$(MSG)
		$(Q)printf "%s" "Get[\"$<\"]; Quit[]" | "$(MATH)" || (echo "Error: The code generation failed!"; exit 1)
		@touch "$(METACODE_STAMP_NUTNMSSM)"
		@echo "Note: to regenerate NUTNMSSM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_NUTNMSSM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_NUTNMSSM):
		@true
endif

$(LIBNUTNMSSM_DEP) $(EXENUTNMSSM_DEP) $(LLNUTNMSSM_DEP) $(LIBNUTNMSSM_OBJ) $(EXENUTNMSSM_OBJ) $(LLNUTNMSSM_OBJ) $(LLNUTNMSSM_LIB): \
	CPPFLAGS += $(MODNUTNMSSM_SUBMOD_INC) $(MODNUTNMSSM_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(GM2CALCFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBNUTNMSSM_DEP) $(EXENUTNMSSM_DEP) $(LLNUTNMSSM_DEP) $(LIBNUTNMSSM_OBJ) $(EXENUTNMSSM_OBJ) $(LLNUTNMSSM_OBJ) $(LLNUTNMSSM_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLNUTNMSSM_OBJ) $(LLNUTNMSSM_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBNUTNMSSM): $(LIBNUTNMSSM_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBNUTNMSSM) $(MODNUTNMSSM_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(FUTILIBS)
		@$(MSG)
		$(Q)$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^) $(filter -%,$(LOOPFUNCLIBS)) $(GM2CALCLIBS) $(HIMALAYALIBS) $(GSLLIBS) $(SQLITELIBS) $(TSILLIBS) $(FLIBS) $(THREADLIBS) $(LDLIBS) $(FUTILIBS)

$(LLNUTNMSSM_LIB): $(LLNUTNMSSM_OBJ) $(LIBNUTNMSSM) $(MODNUTNMSSM_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(FUTILIBS)
		@$(MSG)
		$(Q)$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^) $(filter -%,$(LOOPFUNCLIBS)) $(GM2CALCLIBS) $(HIMALAYALIBS) $(TSILLIBS) $(GSLLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS) $(FUTILIBS) $(FLIBS)

ALLDEP += $(LIBNUTNMSSM_DEP) $(EXENUTNMSSM_DEP)
ALLSRC += $(LIBNUTNMSSM_SRC) $(EXENUTNMSSM_SRC)
ALLLIB += $(LIBNUTNMSSM)
ALLEXE += $(EXENUTNMSSM_EXE)
ALLMODDEP += $(MODNUTNMSSM_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLNUTNMSSM_DEP)
ALLSRC += $(LLNUTNMSSM_SRC)
ALLLL  += $(LLNUTNMSSM_LIB)
endif
