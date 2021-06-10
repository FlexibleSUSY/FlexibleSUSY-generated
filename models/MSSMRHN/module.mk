DIR          := models/MSSMRHN
MODNAME      := MSSMRHN
SARAH_MODEL  := MSSMRHN
WITH_$(MODNAME) := yes
MODMSSMRHN_MOD := SM
MODMSSMRHN_DEP := $(patsubst %,model_specific/%,$(MODMSSMRHN_MOD))
MODMSSMRHN_INC := $(patsubst %,-Imodel_specific/%,$(MODMSSMRHN_MOD))
MODMSSMRHN_LIB := $(foreach M,$(MODMSSMRHN_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

MODMSSMRHN_SUBMOD  := $(DIR)/cxx_qft
MODMSSMRHN_SUBMOD_INC := $(patsubst %,-I%,$(MODMSSMRHN_SUBMOD))

MSSMRHN_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)
MSSMRHN_INSTALL_CXXQFT_DIR := \
		$(MSSMRHN_INSTALL_DIR)/cxx_qft

MSSMRHN_MK     := \
		$(DIR)/module.mk

MSSMRHN_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

MSSMRHN_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

MSSMRHN_CXXQFT_VERTICES_MK := \
		$(DIR)/cxx_qft/vertices.mk

-include $(MSSMRHN_CXXQFT_VERTICES_MK)
LIBMSSMRHN_CXXQFT_VERTICES_SRC ?= ''

MSSMRHN_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

MSSMRHN_FlexibleDecay_MK := \
		$(DIR)/decays/FlexibleDecay.mk

MSSMRHN_INCLUDE_MK := \
		$(MSSMRHN_SUSY_BETAS_MK) \
		$(MSSMRHN_SOFT_BETAS_MK)

MSSMRHN_SLHA_INPUT := \
		$(DIR)/LesHouches.in.MSSMRHN_generated \
		$(DIR)/LesHouches.in.MSSMRHN

MSSMRHN_REFERENCES := \
		$(DIR)/MSSMRHN_references.tex

MSSMRHN_GNUPLOT := \
		$(DIR)/MSSMRHN_plot_rgflow.gnuplot \
		$(DIR)/MSSMRHN_plot_spectrum.gnuplot

MSSMRHN_TARBALL := \
		$(MODNAME).tar.gz

LIBMSSMRHN_SRC := \
		$(DIR)/MSSMRHN_a_muon.cpp \
		$(DIR)/MSSMRHN_edm.cpp \
		$(DIR)/MSSMRHN_FFV_form_factors.cpp \
		$(DIR)/MSSMRHN_f_to_f_conversion.cpp \
		$(DIR)/MSSMRHN_l_to_lgamma.cpp \
		$(DIR)/MSSMRHN_b_to_s_gamma.cpp \
		$(DIR)/MSSMRHN_effective_couplings.cpp \
		$(DIR)/MSSMRHN_info.cpp \
		$(DIR)/MSSMRHN_input_parameters.cpp \
		$(DIR)/MSSMRHN_mass_eigenstates.cpp \
		$(DIR)/MSSMRHN_mass_eigenstates_decoupling_scheme.cpp \
		$(DIR)/MSSMRHN_model_slha.cpp \
		$(DIR)/MSSMRHN_observables.cpp \
		$(DIR)/MSSMRHN_physical.cpp \
		$(DIR)/MSSMRHN_slha_io.cpp \
		$(DIR)/MSSMRHN_soft_parameters.cpp \
		$(DIR)/MSSMRHN_susy_parameters.cpp \
		$(DIR)/MSSMRHN_utilities.cpp \
		$(DIR)/MSSMRHN_weinberg_angle.cpp

LIBMSSMRHN_SRC += $(LIBMSSMRHN_CXXQFT_VERTICES_SRC)

EXEMSSMRHN_SRC := \
		$(DIR)/run_MSSMRHN.cpp \
		$(DIR)/run_cmd_line_MSSMRHN.cpp \
		$(DIR)/scan_MSSMRHN.cpp
LLMSSMRHN_LIB  :=
LLMSSMRHN_OBJ  :=
LLMSSMRHN_SRC  := \
		$(DIR)/MSSMRHN_librarylink.cpp

LLMSSMRHN_MMA  := \
		$(DIR)/MSSMRHN_librarylink.m \
		$(DIR)/run_MSSMRHN.m

LIBMSSMRHN_HDR := \
		$(DIR)/MSSMRHN_a_muon.hpp \
		$(DIR)/MSSMRHN_convergence_tester.hpp \
		$(DIR)/MSSMRHN_edm.hpp \
		$(DIR)/MSSMRHN_FFV_form_factors.hpp \
		$(DIR)/MSSMRHN_f_to_f_conversion.hpp \
		$(DIR)/MSSMRHN_l_to_lgamma.hpp \
		$(DIR)/MSSMRHN_b_to_s_gamma.hpp \
		$(DIR)/MSSMRHN_effective_couplings.hpp \
		$(DIR)/MSSMRHN_ewsb_solver.hpp \
		$(DIR)/MSSMRHN_ewsb_solver_interface.hpp \
		$(DIR)/MSSMRHN_high_scale_constraint.hpp \
		$(DIR)/MSSMRHN_info.hpp \
		$(DIR)/MSSMRHN_initial_guesser.hpp \
		$(DIR)/MSSMRHN_input_parameters.hpp \
		$(DIR)/MSSMRHN_low_scale_constraint.hpp \
		$(DIR)/MSSMRHN_mass_eigenstates.hpp \
		$(DIR)/MSSMRHN_mass_eigenstates_interface.hpp \
		$(DIR)/MSSMRHN_mass_eigenstates_decoupling_scheme.hpp \
		$(DIR)/MSSMRHN_model.hpp \
		$(DIR)/MSSMRHN_model_slha.hpp \
		$(DIR)/MSSMRHN_observables.hpp \
		$(DIR)/MSSMRHN_physical.hpp \
		$(DIR)/MSSMRHN_slha_io.hpp \
		$(DIR)/MSSMRHN_spectrum_generator.hpp \
		$(DIR)/MSSMRHN_spectrum_generator_interface.hpp \
		$(DIR)/MSSMRHN_soft_parameters.hpp \
		$(DIR)/MSSMRHN_susy_parameters.hpp \
		$(DIR)/MSSMRHN_susy_scale_constraint.hpp \
		$(DIR)/MSSMRHN_utilities.hpp \
		$(DIR)/MSSMRHN_weinberg_angle.hpp

LIBMSSMRHN_CXXQFT_HDR := \
		$(DIR)/cxx_qft/MSSMRHN_qft.hpp \
		$(DIR)/cxx_qft/MSSMRHN_fields.hpp \
		$(DIR)/cxx_qft/MSSMRHN_vertices.hpp \
		$(DIR)/cxx_qft/MSSMRHN_context_base.hpp \
		$(DIR)/cxx_qft/MSSMRHN_npointfunctions_wilsoncoeffs.hpp

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
-include $(MSSMRHN_SUSY_BETAS_MK)
-include $(MSSMRHN_SOFT_BETAS_MK)
-include $(MSSMRHN_FlexibleDecay_MK)
-include $(MSSMRHN_CXXQFT_VERTICES_MK)
-include $(MSSMRHN_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(MSSMRHN_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(MSSMRHN_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@

$(MSSMRHN_FlexibleDecay_MK): run-metacode-$(MODNAME)
$(MSSMRHN_CXXQFT_VERTICES_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(MSSMRHN_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBMSSMRHN_SRC := $(sort $(LIBMSSMRHN_SRC))
EXEMSSMRHN_SRC := $(sort $(EXEMSSMRHN_SRC))

LIBMSSMRHN_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBMSSMRHN_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBMSSMRHN_SRC)))

EXEMSSMRHN_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEMSSMRHN_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEMSSMRHN_SRC)))

EXEMSSMRHN_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXEMSSMRHN_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXEMSSMRHN_SRC)))

LIBMSSMRHN_DEP := \
		$(LIBMSSMRHN_OBJ:.o=.d)

EXEMSSMRHN_DEP := \
		$(EXEMSSMRHN_OBJ:.o=.d)

LLMSSMRHN_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLMSSMRHN_SRC)))

LLMSSMRHN_OBJ  := $(LLMSSMRHN_SRC:.cpp=.o)
LLMSSMRHN_LIB  := $(LLMSSMRHN_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBMSSMRHN     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_MSSMRHN := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_MSSMRHN := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL_DIR) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBMSSMRHN) $(EXEMSSMRHN_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(MSSMRHN_INSTALL_DIR)
		$(Q)install -d $(MSSMRHN_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBMSSMRHN_SRC) $(MSSMRHN_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBMSSMRHN_CXXQFT_VERTICES_SRC) $(MSSMRHN_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBMSSMRHN_HDR) $(MSSMRHN_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBMSSMRHN_CXXQFT_HDR) $(MSSMRHN_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(EXEMSSMRHN_SRC) $(MSSMRHN_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLMSSMRHN_SRC) $(MSSMRHN_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLMSSMRHN_MMA) $(MSSMRHN_INSTALL_DIR)
		$(Q)$(INSTALL_STRIPPED) $(MSSMRHN_MK) $(MSSMRHN_INSTALL_DIR) -m u=rw,g=r,o=r
		$(Q)install -m u=rw,g=r,o=r $(MSSMRHN_INCLUDE_MK) $(MSSMRHN_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(MSSMRHN_CXXQFT_VERTICES_MK) $(MSSMRHN_INSTALL_CXXQFT_DIR)

ifneq ($(MSSMRHN_SLHA_INPUT),)
		$(Q)install -m u=rw,g=r,o=r $(MSSMRHN_SLHA_INPUT) $(MSSMRHN_INSTALL_DIR)
endif
		$(Q)install -m u=rw,g=r,o=r $(MSSMRHN_REFERENCES) $(MSSMRHN_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(MSSMRHN_GNUPLOT) $(MSSMRHN_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		$(Q)-rm -f $(LIBMSSMRHN_DEP)
		$(Q)-rm -f $(EXEMSSMRHN_DEP)
		$(Q)-rm -f $(LLMSSMRHN_DEP)

clean-$(MODNAME)-lib:
		$(Q)-rm -f $(LIBMSSMRHN)
		$(Q)-rm -f $(LLMSSMRHN_LIB)

clean-$(MODNAME)-obj:
		$(Q)-rm -f $(LIBMSSMRHN_OBJ)
		$(Q)-rm -f $(EXEMSSMRHN_OBJ)
		$(Q)-rm -f $(LLMSSMRHN_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		$(Q)-rm -f $(LIBMSSMRHN_SRC)
		$(Q)-rm -f $(LIBMSSMRHN_HDR)
		$(Q)-rm -f $(LIBMSSMRHN_CXXQFT_HDR)
		$(Q)-rm -f $(EXEMSSMRHN_SRC)
		$(Q)-rm -f $(LLMSSMRHN_SRC)
		$(Q)-rm -f $(LLMSSMRHN_MMA)
		$(Q)-rm -f $(METACODE_STAMP_MSSMRHN)
		$(Q)-rm -f $(MSSMRHN_INCLUDE_MK)
		$(Q)-rm -f $(MSSMRHN_CXXQFT_VERTICES_MK)
		$(Q)-rm -f $(MSSMRHN_SLHA_INPUT)
		$(Q)-rm -f $(MSSMRHN_REFERENCES)
		$(Q)-rm -f $(MSSMRHN_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		$(Q)-rm -f $(EXEMSSMRHN_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		$(Q)tar -czf $(MSSMRHN_TARBALL) \
		$(LIBMSSMRHN_SRC) $(LIBMSSMRHN_HDR) $(LIBMSSMRHN_CXXQFT_HDR) \
		$(EXEMSSMRHN_SRC) \
		$(LLMSSMRHN_SRC) $(LLMSSMRHN_MMA) \
		$(MSSMRHN_MK) $(MSSMRHN_INCLUDE_MK) $(MSSMRHN_CXXQFT_VERTICES_MK) \
		$(MSSMRHN_SLHA_INPUT) $(MSSMRHN_REFERENCES) \
		$(MSSMRHN_GNUPLOT) \
		$(MSSMRHN_FlexibleDecay_MK)

$(LIBMSSMRHN_SRC) $(LIBMSSMRHN_HDR) $(LIBMSSMRHN_CXXQFT_HDR) $(EXEMSSMRHN_SRC) $(LLMSSMRHN_SRC) $(LLMSSMRHN_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_MSSMRHN)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_MSSMRHN): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_MSSMRHN)
		@$(MSG)
		$(Q)printf "%s" "Get[\"$<\"]; Quit[]" | "$(MATH)" || (echo "Error: The code generation failed!"; exit 1)
		@touch "$(METACODE_STAMP_MSSMRHN)"
		@echo "Note: to regenerate MSSMRHN source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_MSSMRHN)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_MSSMRHN):
		@true
endif

$(LIBMSSMRHN_DEP) $(EXEMSSMRHN_DEP) $(LLMSSMRHN_DEP) $(LIBMSSMRHN_OBJ) $(EXEMSSMRHN_OBJ) $(LLMSSMRHN_OBJ) $(LLMSSMRHN_LIB): \
	CPPFLAGS += $(MODMSSMRHN_SUBMOD_INC) $(MODMSSMRHN_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(GM2CALCFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBMSSMRHN_DEP) $(EXEMSSMRHN_DEP) $(LLMSSMRHN_DEP) $(LIBMSSMRHN_OBJ) $(EXEMSSMRHN_OBJ) $(LLMSSMRHN_OBJ) $(LLMSSMRHN_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLMSSMRHN_OBJ) $(LLMSSMRHN_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBMSSMRHN): $(LIBMSSMRHN_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBMSSMRHN) $(MODMSSMRHN_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(FUTILIBS)
		@$(MSG)
		$(Q)$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^) $(filter -%,$(LOOPFUNCLIBS)) $(GM2CALCLIBS) $(HIMALAYALIBS) $(GSLLIBS) $(SQLITELIBS) $(TSILLIBS) $(FLIBS) $(THREADLIBS) $(LDLIBS) $(FUTILIBS)

$(LLMSSMRHN_LIB): $(LLMSSMRHN_OBJ) $(LIBMSSMRHN) $(MODMSSMRHN_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(FUTILIBS)
		@$(MSG)
		$(Q)$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^) $(filter -%,$(LOOPFUNCLIBS)) $(GM2CALCLIBS) $(HIMALAYALIBS) $(TSILLIBS) $(GSLLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS) $(FUTILIBS) $(FLIBS)

ALLDEP += $(LIBMSSMRHN_DEP) $(EXEMSSMRHN_DEP)
ALLSRC += $(LIBMSSMRHN_SRC) $(EXEMSSMRHN_SRC)
ALLLIB += $(LIBMSSMRHN)
ALLEXE += $(EXEMSSMRHN_EXE)
ALLMODDEP += $(MODMSSMRHN_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLMSSMRHN_DEP)
ALLSRC += $(LLMSSMRHN_SRC)
ALLLL  += $(LLMSSMRHN_LIB)
endif
