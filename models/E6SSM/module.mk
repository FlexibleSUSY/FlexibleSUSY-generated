DIR          := models/E6SSM
MODNAME      := E6SSM
SARAH_MODEL  := E6SSM
WITH_$(MODNAME) := yes
MODE6SSM_MOD := SM MSSM_higgs NMSSM_higgs
MODE6SSM_DEP := $(patsubst %,model_specific/%,$(MODE6SSM_MOD))
MODE6SSM_INC := $(patsubst %,-Imodel_specific/%,$(MODE6SSM_MOD))
MODE6SSM_LIB := $(foreach M,$(MODE6SSM_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

MODE6SSM_SUBMOD  := $(DIR)/cxx_qft
MODE6SSM_SUBMOD_INC := $(patsubst %,-I%,$(MODE6SSM_SUBMOD))

E6SSM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)
E6SSM_INSTALL_CXXQFT_DIR := \
		$(E6SSM_INSTALL_DIR)/cxx_qft

E6SSM_MK     := \
		$(DIR)/module.mk

E6SSM_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

E6SSM_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

E6SSM_CXXQFT_VERTICES_MK := \
		$(DIR)/cxx_qft/vertices.mk

-include $(E6SSM_CXXQFT_VERTICES_MK)
LIBE6SSM_CXXQFT_VERTICES_SRC ?= ''

E6SSM_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

E6SSM_FlexibleDecay_MK := \
		$(DIR)/decays/FlexibleDecay.mk

E6SSM_INCLUDE_MK := \
		$(E6SSM_SUSY_BETAS_MK) \
		$(E6SSM_SOFT_BETAS_MK)

E6SSM_SLHA_INPUT := \
		$(DIR)/LesHouches.in.E6SSM_generated \
		$(DIR)/LesHouches.in.E6SSM

E6SSM_REFERENCES := \
		$(DIR)/E6SSM_references.tex

E6SSM_GNUPLOT := \
		$(DIR)/E6SSM_plot_rgflow.gnuplot \
		$(DIR)/E6SSM_plot_spectrum.gnuplot

E6SSM_TARBALL := \
		$(MODNAME).tar.gz

LIBE6SSM_SRC := \
		$(DIR)/E6SSM_amm.cpp \
		$(DIR)/E6SSM_edm.cpp \
		$(DIR)/E6SSM_FFV_form_factors.cpp \
		$(wildcard $(DIR)/observables/E6SSM*.cpp) \
		$(DIR)/E6SSM_b_to_s_gamma.cpp \
		$(DIR)/E6SSM_info.cpp \
		$(DIR)/E6SSM_input_parameters.cpp \
		$(DIR)/E6SSM_mass_eigenstates.cpp \
		$(DIR)/E6SSM_mass_eigenstates_decoupling_scheme.cpp \
		$(DIR)/E6SSM_model_slha.cpp \
		$(DIR)/E6SSM_lepton_amm_wrapper.cpp \
		$(DIR)/E6SSM_observables.cpp \
		$(DIR)/E6SSM_physical.cpp \
		$(DIR)/E6SSM_slha_io.cpp \
		$(DIR)/E6SSM_soft_parameters.cpp \
		$(DIR)/E6SSM_susy_parameters.cpp \
		$(DIR)/E6SSM_unitarity.cpp \
		$(DIR)/E6SSM_utilities.cpp \
		$(DIR)/E6SSM_weinberg_angle.cpp

LIBE6SSM_SRC += $(LIBE6SSM_CXXQFT_VERTICES_SRC)

EXEE6SSM_SRC := \
		$(DIR)/run_E6SSM.cpp \
		$(DIR)/run_cmd_line_E6SSM.cpp \
		$(DIR)/scan_E6SSM.cpp
LLE6SSM_LIB  :=
LLE6SSM_OBJ  :=
LLE6SSM_SRC  := \
		$(DIR)/E6SSM_librarylink.cpp

LLE6SSM_MMA  := \
		$(DIR)/E6SSM_librarylink.m \
		$(DIR)/run_E6SSM.m

LIBE6SSM_HDR := \
		$(DIR)/E6SSM_amm.hpp \
		$(DIR)/E6SSM_convergence_tester.hpp \
		$(DIR)/E6SSM_edm.hpp \
		$(DIR)/E6SSM_FFV_form_factors.hpp \
		$(wildcard $(DIR)/observables/E6SSM*.hpp) \
		$(DIR)/E6SSM_b_to_s_gamma.hpp \
		$(DIR)/E6SSM_ewsb_solver.hpp \
		$(DIR)/E6SSM_ewsb_solver_interface.hpp \
		$(DIR)/E6SSM_high_scale_constraint.hpp \
		$(DIR)/E6SSM_info.hpp \
		$(DIR)/E6SSM_initial_guesser.hpp \
		$(DIR)/E6SSM_input_parameters.hpp \
		$(DIR)/E6SSM_low_scale_constraint.hpp \
		$(DIR)/E6SSM_mass_eigenstates.hpp \
		$(DIR)/E6SSM_mass_eigenstates_interface.hpp \
		$(DIR)/E6SSM_mass_eigenstates_decoupling_scheme.hpp \
		$(DIR)/E6SSM_model.hpp \
		$(DIR)/E6SSM_model_slha.hpp \
		$(DIR)/E6SSM_lepton_amm_wrapper.hpp \
		$(DIR)/E6SSM_observables.hpp \
		$(DIR)/E6SSM_physical.hpp \
		$(DIR)/E6SSM_slha_io.hpp \
		$(DIR)/E6SSM_spectrum_generator.hpp \
		$(DIR)/E6SSM_spectrum_generator_interface.hpp \
		$(DIR)/E6SSM_soft_parameters.hpp \
		$(DIR)/E6SSM_susy_parameters.hpp \
		$(DIR)/E6SSM_susy_scale_constraint.hpp \
		$(DIR)/E6SSM_unitarity.hpp \
		$(DIR)/E6SSM_utilities.hpp \
		$(DIR)/E6SSM_weinberg_angle.hpp

LIBE6SSM_CXXQFT_HDR := \
		$(DIR)/cxx_qft/E6SSM_qft.hpp \
		$(DIR)/cxx_qft/E6SSM_fields.hpp \
		$(DIR)/cxx_qft/E6SSM_particle_aliases.hpp \
		$(DIR)/cxx_qft/E6SSM_vertices.hpp \
		$(DIR)/cxx_qft/E6SSM_context_base.hpp \
		$(DIR)/cxx_qft/E6SSM_npointfunctions_wilsoncoeffs.hpp

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
-include $(E6SSM_SUSY_BETAS_MK)
-include $(E6SSM_SOFT_BETAS_MK)
-include $(E6SSM_FlexibleDecay_MK)
-include $(E6SSM_CXXQFT_VERTICES_MK)
-include $(E6SSM_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(E6SSM_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(E6SSM_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@

$(E6SSM_FlexibleDecay_MK): run-metacode-$(MODNAME)
$(E6SSM_CXXQFT_VERTICES_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(E6SSM_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBE6SSM_SRC := $(sort $(LIBE6SSM_SRC))
EXEE6SSM_SRC := $(sort $(EXEE6SSM_SRC))

LIBE6SSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBE6SSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBE6SSM_SRC)))

EXEE6SSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEE6SSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEE6SSM_SRC)))

EXEE6SSM_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXEE6SSM_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXEE6SSM_SRC)))

LIBE6SSM_DEP := \
		$(LIBE6SSM_OBJ:.o=.d)

EXEE6SSM_DEP := \
		$(EXEE6SSM_OBJ:.o=.d)

LLE6SSM_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLE6SSM_SRC)))

LLE6SSM_OBJ  := $(LLE6SSM_SRC:.cpp=.o)
LLE6SSM_LIB  := $(LLE6SSM_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBE6SSM     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_E6SSM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_E6SSM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL_DIR) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBE6SSM) $(EXEE6SSM_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(E6SSM_INSTALL_DIR)
		$(Q)install -d $(E6SSM_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBE6SSM_SRC) $(E6SSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBE6SSM_CXXQFT_VERTICES_SRC) $(E6SSM_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBE6SSM_HDR) $(E6SSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBE6SSM_CXXQFT_HDR) $(E6SSM_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(EXEE6SSM_SRC) $(E6SSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLE6SSM_SRC) $(E6SSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLE6SSM_MMA) $(E6SSM_INSTALL_DIR)
		$(Q)$(INSTALL_STRIPPED) $(E6SSM_MK) $(E6SSM_INSTALL_DIR) -m u=rw,g=r,o=r
		$(Q)install -m u=rw,g=r,o=r $(E6SSM_INCLUDE_MK) $(E6SSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(E6SSM_CXXQFT_VERTICES_MK) $(E6SSM_INSTALL_CXXQFT_DIR)

ifneq ($(E6SSM_SLHA_INPUT),)
		$(Q)install -m u=rw,g=r,o=r $(E6SSM_SLHA_INPUT) $(E6SSM_INSTALL_DIR)
endif
		$(Q)install -m u=rw,g=r,o=r $(E6SSM_REFERENCES) $(E6SSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(E6SSM_GNUPLOT) $(E6SSM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		$(Q)-rm -f $(LIBE6SSM_DEP)
		$(Q)-rm -f $(EXEE6SSM_DEP)
		$(Q)-rm -f $(LLE6SSM_DEP)

clean-$(MODNAME)-lib:
		$(Q)-rm -f $(LIBE6SSM)
		$(Q)-rm -f $(LLE6SSM_LIB)

clean-$(MODNAME)-obj:
		$(Q)-rm -f $(LIBE6SSM_OBJ)
		$(Q)-rm -f $(EXEE6SSM_OBJ)
		$(Q)-rm -f $(LLE6SSM_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		$(Q)-rm -f $(LIBE6SSM_SRC)
		$(Q)-rm -f $(LIBE6SSM_HDR)
		$(Q)-rm -f $(LIBE6SSM_CXXQFT_HDR)
		$(Q)-rm -f $(EXEE6SSM_SRC)
		$(Q)-rm -f $(LLE6SSM_SRC)
		$(Q)-rm -f $(LLE6SSM_MMA)
		$(Q)-rm -f $(METACODE_STAMP_E6SSM)
		$(Q)-rm -f $(E6SSM_INCLUDE_MK)
		$(Q)-rm -f $(E6SSM_CXXQFT_VERTICES_MK)
		$(Q)-rm -f $(E6SSM_SLHA_INPUT)
		$(Q)-rm -f $(E6SSM_REFERENCES)
		$(Q)-rm -f $(E6SSM_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		$(Q)-rm -f $(EXEE6SSM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		$(Q)tar -czf $(E6SSM_TARBALL) \
		$(LIBE6SSM_SRC) $(LIBE6SSM_HDR) $(LIBE6SSM_CXXQFT_HDR) \
		$(EXEE6SSM_SRC) \
		$(LLE6SSM_SRC) $(LLE6SSM_MMA) \
		$(E6SSM_MK) $(E6SSM_INCLUDE_MK) $(E6SSM_CXXQFT_VERTICES_MK) \
		$(E6SSM_SLHA_INPUT) $(E6SSM_REFERENCES) \
		$(E6SSM_GNUPLOT) \
		$(E6SSM_FlexibleDecay_MK)

$(LIBE6SSM_SRC) $(LIBE6SSM_HDR) $(LIBE6SSM_CXXQFT_HDR) $(EXEE6SSM_SRC) $(LLE6SSM_SRC) $(LLE6SSM_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_E6SSM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_E6SSM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_E6SSM)
		@$(MSG)
		$(Q)printf "%s" "Get[\"$<\"]; Quit[]" | "$(MATH)" || (echo "Error: The code generation failed!"; exit 1)
		@touch "$(METACODE_STAMP_E6SSM)"
		@echo "Note: to regenerate E6SSM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_E6SSM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_E6SSM):
		@true
endif

$(LIBE6SSM_DEP) $(EXEE6SSM_DEP) $(LLE6SSM_DEP) $(LIBE6SSM_OBJ) $(EXEE6SSM_OBJ) $(LLE6SSM_OBJ) $(LLE6SSM_LIB): \
	CPPFLAGS += $(MODE6SSM_SUBMOD_INC) $(MODE6SSM_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(GM2CALCFLAGS) $(HIGGSTOOLSFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBE6SSM_DEP) $(EXEE6SSM_DEP) $(LLE6SSM_DEP) $(LIBE6SSM_OBJ) $(EXEE6SSM_OBJ) $(LLE6SSM_OBJ) $(LLE6SSM_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLE6SSM_OBJ) $(LLE6SSM_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBE6SSM): $(LIBE6SSM_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBE6SSM) $(MODE6SSM_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(FUTILIBS)
		@$(MSG)
		$(Q)$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^) $(filter -%,$(LOOPFUNCLIBS)) $(GM2CALCLIBS) $(HIGGSTOOLSLIBS) $(PYTHONLIBS) $(HIMALAYALIBS) $(GSLLIBS) $(SQLITELIBS) $(TSILLIBS) $(FLIBS) $(THREADLIBS) $(LDLIBS) $(FUTILIBS)

$(LLE6SSM_LIB): $(LLE6SSM_OBJ) $(LIBE6SSM) $(MODE6SSM_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(FUTILIBS)
		@$(MSG)
		$(Q)$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^) $(filter -%,$(LOOPFUNCLIBS)) $(GM2CALCLIBS) $(HIGGSTOOLSLIBS) $(PYTHONLIBS) $(HIMALAYALIBS) $(TSILLIBS) $(GSLLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS) $(FUTILIBS) $(FLIBS)

ALLDEP += $(LIBE6SSM_DEP) $(EXEE6SSM_DEP)
ALLSRC += $(LIBE6SSM_SRC) $(EXEE6SSM_SRC)
ALLLIB += $(LIBE6SSM)
ALLEXE += $(EXEE6SSM_EXE)
ALLMODDEP += $(MODE6SSM_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLE6SSM_DEP)
ALLSRC += $(LLE6SSM_SRC)
ALLLL  += $(LLE6SSM_LIB)
endif
