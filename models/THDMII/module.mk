DIR          := models/THDMII
MODNAME      := THDMII
SARAH_MODEL  := THDM-II
WITH_$(MODNAME) := yes
MODTHDMII_MOD := SM
MODTHDMII_DEP := $(patsubst %,model_specific/%,$(MODTHDMII_MOD))
MODTHDMII_INC := $(patsubst %,-Imodel_specific/%,$(MODTHDMII_MOD))
MODTHDMII_LIB := $(foreach M,$(MODTHDMII_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

MODTHDMII_SUBMOD  := $(DIR)/cxx_qft
MODTHDMII_SUBMOD_INC := $(patsubst %,-I%,$(MODTHDMII_SUBMOD))

THDMII_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)
THDMII_INSTALL_CXXQFT_DIR := \
		$(THDMII_INSTALL_DIR)/cxx_qft

THDMII_MK     := \
		$(DIR)/module.mk

THDMII_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

THDMII_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

THDMII_CXXQFT_VERTICES_MK := \
		$(DIR)/cxx_qft/vertices.mk

-include $(THDMII_CXXQFT_VERTICES_MK)
LIBTHDMII_CXXQFT_VERTICES_SRC ?= ''

THDMII_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

THDMII_FlexibleDecay_MK := \
		$(DIR)/decays/FlexibleDecay.mk

THDMII_INCLUDE_MK := \
		$(THDMII_SUSY_BETAS_MK) \
		$(THDMII_SOFT_BETAS_MK)

THDMII_SLHA_INPUT := \
		$(DIR)/LesHouches.in.THDMII_generated \
		$(DIR)/LesHouches.in.THDMII

THDMII_REFERENCES := \
		$(DIR)/THDMII_references.tex

THDMII_GNUPLOT := \
		$(DIR)/THDMII_plot_rgflow.gnuplot \
		$(DIR)/THDMII_plot_spectrum.gnuplot

THDMII_TARBALL := \
		$(MODNAME).tar.gz

LIBTHDMII_SRC := \
		$(DIR)/THDMII_amm.cpp \
		$(DIR)/THDMII_edm.cpp \
		$(DIR)/THDMII_FFV_form_factors.cpp \
		$(wildcard $(DIR)/observables/THDMII*.cpp) \
		$(DIR)/THDMII_b_to_s_gamma.cpp \
		$(DIR)/THDMII_info.cpp \
		$(DIR)/THDMII_input_parameters.cpp \
		$(DIR)/THDMII_mass_eigenstates.cpp \
		$(DIR)/THDMII_mass_eigenstates_decoupling_scheme.cpp \
		$(DIR)/THDMII_model_slha.cpp \
		$(DIR)/THDMII_lepton_amm_wrapper.cpp \
		$(DIR)/THDMII_observables.cpp \
		$(DIR)/THDMII_physical.cpp \
		$(DIR)/THDMII_slha_io.cpp \
		$(DIR)/THDMII_soft_parameters.cpp \
		$(DIR)/THDMII_susy_parameters.cpp \
		$(DIR)/THDMII_unitarity.cpp \
		$(DIR)/THDMII_utilities.cpp \
		$(DIR)/THDMII_weinberg_angle.cpp

LIBTHDMII_SRC += $(LIBTHDMII_CXXQFT_VERTICES_SRC)

EXETHDMII_SRC := \
		$(DIR)/run_THDMII.cpp \
		$(DIR)/run_cmd_line_THDMII.cpp \
		$(DIR)/scan_THDMII.cpp
LLTHDMII_LIB  :=
LLTHDMII_OBJ  :=
LLTHDMII_SRC  := \
		$(DIR)/THDMII_librarylink.cpp

LLTHDMII_MMA  := \
		$(DIR)/THDMII_librarylink.m \
		$(DIR)/run_THDMII.m

LIBTHDMII_HDR := \
		$(DIR)/THDMII_amm.hpp \
		$(DIR)/THDMII_convergence_tester.hpp \
		$(DIR)/THDMII_edm.hpp \
		$(DIR)/THDMII_FFV_form_factors.hpp \
		$(wildcard $(DIR)/observables/THDMII*.hpp) \
		$(DIR)/THDMII_b_to_s_gamma.hpp \
		$(DIR)/THDMII_ewsb_solver.hpp \
		$(DIR)/THDMII_ewsb_solver_interface.hpp \
		$(DIR)/THDMII_high_scale_constraint.hpp \
		$(DIR)/THDMII_info.hpp \
		$(DIR)/THDMII_initial_guesser.hpp \
		$(DIR)/THDMII_input_parameters.hpp \
		$(DIR)/THDMII_low_scale_constraint.hpp \
		$(DIR)/THDMII_mass_eigenstates.hpp \
		$(DIR)/THDMII_mass_eigenstates_interface.hpp \
		$(DIR)/THDMII_mass_eigenstates_decoupling_scheme.hpp \
		$(DIR)/THDMII_model.hpp \
		$(DIR)/THDMII_model_slha.hpp \
		$(DIR)/THDMII_lepton_amm_wrapper.hpp \
		$(DIR)/THDMII_observables.hpp \
		$(DIR)/THDMII_physical.hpp \
		$(DIR)/THDMII_slha_io.hpp \
		$(DIR)/THDMII_spectrum_generator.hpp \
		$(DIR)/THDMII_spectrum_generator_interface.hpp \
		$(DIR)/THDMII_soft_parameters.hpp \
		$(DIR)/THDMII_susy_parameters.hpp \
		$(DIR)/THDMII_susy_scale_constraint.hpp \
		$(DIR)/THDMII_unitarity.hpp \
		$(DIR)/THDMII_utilities.hpp \
		$(DIR)/THDMII_weinberg_angle.hpp

LIBTHDMII_CXXQFT_HDR := \
		$(DIR)/cxx_qft/THDMII_qft.hpp \
		$(DIR)/cxx_qft/THDMII_fields.hpp \
		$(DIR)/cxx_qft/THDMII_particle_aliases.hpp \
		$(DIR)/cxx_qft/THDMII_vertices.hpp \
		$(DIR)/cxx_qft/THDMII_context_base.hpp \
		$(DIR)/cxx_qft/THDMII_npointfunctions_wilsoncoeffs.hpp

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
-include $(THDMII_SUSY_BETAS_MK)
-include $(THDMII_SOFT_BETAS_MK)
-include $(THDMII_FlexibleDecay_MK)
-include $(THDMII_CXXQFT_VERTICES_MK)
-include $(THDMII_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(THDMII_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(THDMII_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@

$(THDMII_FlexibleDecay_MK): run-metacode-$(MODNAME)
$(THDMII_CXXQFT_VERTICES_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(THDMII_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBTHDMII_SRC := $(sort $(LIBTHDMII_SRC))
EXETHDMII_SRC := $(sort $(EXETHDMII_SRC))

LIBTHDMII_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBTHDMII_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBTHDMII_SRC)))

EXETHDMII_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXETHDMII_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXETHDMII_SRC)))

EXETHDMII_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXETHDMII_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXETHDMII_SRC)))

LIBTHDMII_DEP := \
		$(LIBTHDMII_OBJ:.o=.d)

EXETHDMII_DEP := \
		$(EXETHDMII_OBJ:.o=.d)

LLTHDMII_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLTHDMII_SRC)))

LLTHDMII_OBJ  := $(LLTHDMII_SRC:.cpp=.o)
LLTHDMII_LIB  := $(LLTHDMII_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBTHDMII     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_THDMII := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_THDMII := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL_DIR) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBTHDMII) $(EXETHDMII_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(THDMII_INSTALL_DIR)
		$(Q)install -d $(THDMII_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBTHDMII_SRC) $(THDMII_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBTHDMII_CXXQFT_VERTICES_SRC) $(THDMII_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBTHDMII_HDR) $(THDMII_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBTHDMII_CXXQFT_HDR) $(THDMII_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(EXETHDMII_SRC) $(THDMII_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLTHDMII_SRC) $(THDMII_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLTHDMII_MMA) $(THDMII_INSTALL_DIR)
		$(Q)$(INSTALL_STRIPPED) $(THDMII_MK) $(THDMII_INSTALL_DIR) -m u=rw,g=r,o=r
		$(Q)install -m u=rw,g=r,o=r $(THDMII_INCLUDE_MK) $(THDMII_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(THDMII_CXXQFT_VERTICES_MK) $(THDMII_INSTALL_CXXQFT_DIR)

ifneq ($(THDMII_SLHA_INPUT),)
		$(Q)install -m u=rw,g=r,o=r $(THDMII_SLHA_INPUT) $(THDMII_INSTALL_DIR)
endif
		$(Q)install -m u=rw,g=r,o=r $(THDMII_REFERENCES) $(THDMII_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(THDMII_GNUPLOT) $(THDMII_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		$(Q)-rm -f $(LIBTHDMII_DEP)
		$(Q)-rm -f $(EXETHDMII_DEP)
		$(Q)-rm -f $(LLTHDMII_DEP)

clean-$(MODNAME)-lib:
		$(Q)-rm -f $(LIBTHDMII)
		$(Q)-rm -f $(LLTHDMII_LIB)

clean-$(MODNAME)-obj:
		$(Q)-rm -f $(LIBTHDMII_OBJ)
		$(Q)-rm -f $(EXETHDMII_OBJ)
		$(Q)-rm -f $(LLTHDMII_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		$(Q)-rm -f $(LIBTHDMII_SRC)
		$(Q)-rm -f $(LIBTHDMII_HDR)
		$(Q)-rm -f $(LIBTHDMII_CXXQFT_HDR)
		$(Q)-rm -f $(EXETHDMII_SRC)
		$(Q)-rm -f $(LLTHDMII_SRC)
		$(Q)-rm -f $(LLTHDMII_MMA)
		$(Q)-rm -f $(METACODE_STAMP_THDMII)
		$(Q)-rm -f $(THDMII_INCLUDE_MK)
		$(Q)-rm -f $(THDMII_CXXQFT_VERTICES_MK)
		$(Q)-rm -f $(THDMII_SLHA_INPUT)
		$(Q)-rm -f $(THDMII_REFERENCES)
		$(Q)-rm -f $(THDMII_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		$(Q)-rm -f $(EXETHDMII_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		$(Q)tar -czf $(THDMII_TARBALL) \
		$(LIBTHDMII_SRC) $(LIBTHDMII_HDR) $(LIBTHDMII_CXXQFT_HDR) \
		$(EXETHDMII_SRC) \
		$(LLTHDMII_SRC) $(LLTHDMII_MMA) \
		$(THDMII_MK) $(THDMII_INCLUDE_MK) $(THDMII_CXXQFT_VERTICES_MK) \
		$(THDMII_SLHA_INPUT) $(THDMII_REFERENCES) \
		$(THDMII_GNUPLOT) \
		$(THDMII_FlexibleDecay_MK)

$(LIBTHDMII_SRC) $(LIBTHDMII_HDR) $(LIBTHDMII_CXXQFT_HDR) $(EXETHDMII_SRC) $(LLTHDMII_SRC) $(LLTHDMII_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_THDMII)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_THDMII): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_THDMII)
		@$(MSG)
		$(Q)printf "%s" "Get[\"$<\"]; Quit[]" | "$(MATH)" || (echo "Error: The code generation failed!"; exit 1)
		@touch "$(METACODE_STAMP_THDMII)"
		@echo "Note: to regenerate THDMII source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_THDMII)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_THDMII):
		@true
endif

$(LIBTHDMII_DEP) $(EXETHDMII_DEP) $(LLTHDMII_DEP) $(LIBTHDMII_OBJ) $(EXETHDMII_OBJ) $(LLTHDMII_OBJ) $(LLTHDMII_LIB): \
	CPPFLAGS += $(MODTHDMII_SUBMOD_INC) $(MODTHDMII_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(GM2CALCFLAGS) $(HIGGSTOOLSFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBTHDMII_DEP) $(EXETHDMII_DEP) $(LLTHDMII_DEP) $(LIBTHDMII_OBJ) $(EXETHDMII_OBJ) $(LLTHDMII_OBJ) $(LLTHDMII_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLTHDMII_OBJ) $(LLTHDMII_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBTHDMII): $(LIBTHDMII_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBTHDMII) $(MODTHDMII_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(FUTILIBS)
		@$(MSG)
		$(Q)$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^) $(filter -%,$(LOOPFUNCLIBS)) $(GM2CALCLIBS) $(HIGGSTOOLSLIBS) $(PYTHONLIBS) $(HIMALAYALIBS) $(GSLLIBS) $(SQLITELIBS) $(TSILLIBS) $(FLIBS) $(THREADLIBS) $(LDLIBS) $(FUTILIBS)

$(LLTHDMII_LIB): $(LLTHDMII_OBJ) $(LIBTHDMII) $(MODTHDMII_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(FUTILIBS)
		@$(MSG)
		$(Q)$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^) $(filter -%,$(LOOPFUNCLIBS)) $(GM2CALCLIBS) $(HIGGSTOOLSLIBS) $(PYTHONLIBS) $(HIMALAYALIBS) $(TSILLIBS) $(GSLLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS) $(FUTILIBS) $(FLIBS)

ALLDEP += $(LIBTHDMII_DEP) $(EXETHDMII_DEP)
ALLSRC += $(LIBTHDMII_SRC) $(EXETHDMII_SRC)
ALLLIB += $(LIBTHDMII)
ALLEXE += $(EXETHDMII_EXE)
ALLMODDEP += $(MODTHDMII_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLTHDMII_DEP)
ALLSRC += $(LLTHDMII_SRC)
ALLLL  += $(LLTHDMII_LIB)
endif
