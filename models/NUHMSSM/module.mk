DIR          := models/NUHMSSM
MODNAME      := NUHMSSM
SARAH_MODEL  := MSSM
WITH_$(MODNAME) := yes
MODNUHMSSM_MOD := SM MSSM_higgs
MODNUHMSSM_DEP := $(patsubst %,model_specific/%,$(MODNUHMSSM_MOD))
MODNUHMSSM_INC := $(patsubst %,-Imodel_specific/%,$(MODNUHMSSM_MOD))
MODNUHMSSM_LIB := $(foreach M,$(MODNUHMSSM_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

MODNUHMSSM_SUBMOD  := $(DIR)/cxx_qft
MODNUHMSSM_SUBMOD_INC := $(patsubst %,-I%,$(MODNUHMSSM_SUBMOD))

NUHMSSM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)
NUHMSSM_INSTALL_CXXQFT_DIR := \
		$(NUHMSSM_INSTALL_DIR)/cxx_qft

NUHMSSM_MK     := \
		$(DIR)/module.mk

NUHMSSM_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

NUHMSSM_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

NUHMSSM_CXXQFT_VERTICES_MK := \
		$(DIR)/cxx_qft/vertices.mk

-include $(NUHMSSM_CXXQFT_VERTICES_MK)
LIBNUHMSSM_CXXQFT_VERTICES_SRC ?= ''

NUHMSSM_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

NUHMSSM_INCLUDE_MK := \
		$(NUHMSSM_SUSY_BETAS_MK) \
		$(NUHMSSM_SOFT_BETAS_MK)

NUHMSSM_SLHA_INPUT := \
		$(DIR)/LesHouches.in.NUHMSSM_generated \
		$(DIR)/LesHouches.in.NUHMSSM

NUHMSSM_REFERENCES := \
		$(DIR)/NUHMSSM_references.tex

NUHMSSM_GNUPLOT := \
		$(DIR)/NUHMSSM_plot_rgflow.gnuplot \
		$(DIR)/NUHMSSM_plot_spectrum.gnuplot

NUHMSSM_TARBALL := \
		$(MODNAME).tar.gz

LIBNUHMSSM_SRC := \
		$(DIR)/NUHMSSM_a_muon.cpp \
		$(DIR)/NUHMSSM_edm.cpp \
		$(DIR)/NUHMSSM_FFV_form_factors.cpp \
		$(DIR)/NUHMSSM_f_to_f_conversion.cpp \
		$(DIR)/NUHMSSM_l_to_lgamma.cpp \
		$(DIR)/NUHMSSM_b_to_s_gamma.cpp \
		$(DIR)/NUHMSSM_effective_couplings.cpp \
		$(DIR)/NUHMSSM_info.cpp \
		$(DIR)/NUHMSSM_input_parameters.cpp \
		$(DIR)/NUHMSSM_mass_eigenstates.cpp \
		$(DIR)/NUHMSSM_mass_eigenstates_decoupling_scheme.cpp \
		$(DIR)/NUHMSSM_model_slha.cpp \
		$(DIR)/NUHMSSM_observables.cpp \
		$(DIR)/NUHMSSM_physical.cpp \
		$(DIR)/NUHMSSM_slha_io.cpp \
		$(DIR)/NUHMSSM_soft_parameters.cpp \
		$(DIR)/NUHMSSM_susy_parameters.cpp \
		$(DIR)/NUHMSSM_utilities.cpp \
		$(DIR)/NUHMSSM_weinberg_angle.cpp

LIBNUHMSSM_SRC += $(LIBNUHMSSM_CXXQFT_VERTICES_SRC)

EXENUHMSSM_SRC := \
		$(DIR)/run_NUHMSSM.cpp \
		$(DIR)/run_cmd_line_NUHMSSM.cpp \
		$(DIR)/scan_NUHMSSM.cpp
LLNUHMSSM_LIB  :=
LLNUHMSSM_OBJ  :=
LLNUHMSSM_SRC  := \
		$(DIR)/NUHMSSM_librarylink.cpp

LLNUHMSSM_MMA  := \
		$(DIR)/NUHMSSM_librarylink.m \
		$(DIR)/run_NUHMSSM.m

LIBNUHMSSM_HDR := \
		$(DIR)/NUHMSSM_a_muon.hpp \
		$(DIR)/NUHMSSM_convergence_tester.hpp \
		$(DIR)/NUHMSSM_edm.hpp \
		$(DIR)/NUHMSSM_FFV_form_factors.hpp \
		$(DIR)/NUHMSSM_f_to_f_conversion.hpp \
		$(DIR)/NUHMSSM_l_to_lgamma.hpp \
		$(DIR)/NUHMSSM_b_to_s_gamma.hpp \
		$(DIR)/NUHMSSM_effective_couplings.hpp \
		$(DIR)/NUHMSSM_ewsb_solver.hpp \
		$(DIR)/NUHMSSM_ewsb_solver_interface.hpp \
		$(DIR)/NUHMSSM_high_scale_constraint.hpp \
		$(DIR)/NUHMSSM_info.hpp \
		$(DIR)/NUHMSSM_initial_guesser.hpp \
		$(DIR)/NUHMSSM_input_parameters.hpp \
		$(DIR)/NUHMSSM_low_scale_constraint.hpp \
		$(DIR)/NUHMSSM_mass_eigenstates.hpp \
		$(DIR)/NUHMSSM_mass_eigenstates_interface.hpp \
		$(DIR)/NUHMSSM_mass_eigenstates_decoupling_scheme.hpp \
		$(DIR)/NUHMSSM_model.hpp \
		$(DIR)/NUHMSSM_model_slha.hpp \
		$(DIR)/NUHMSSM_observables.hpp \
		$(DIR)/NUHMSSM_physical.hpp \
		$(DIR)/NUHMSSM_slha_io.hpp \
		$(DIR)/NUHMSSM_spectrum_generator.hpp \
		$(DIR)/NUHMSSM_spectrum_generator_interface.hpp \
		$(DIR)/NUHMSSM_soft_parameters.hpp \
		$(DIR)/NUHMSSM_susy_parameters.hpp \
		$(DIR)/NUHMSSM_susy_scale_constraint.hpp \
		$(DIR)/NUHMSSM_utilities.hpp \
		$(DIR)/NUHMSSM_weinberg_angle.hpp

LIBNUHMSSM_CXXQFT_HDR := \
		$(DIR)/cxx_qft/NUHMSSM_qft.hpp \
		$(DIR)/cxx_qft/NUHMSSM_fields.hpp \
		$(DIR)/cxx_qft/NUHMSSM_vertices.hpp \
		$(DIR)/cxx_qft/NUHMSSM_context_base.hpp \
		$(DIR)/cxx_qft/NUHMSSM_npointfunctions_wilsoncoeffs.hpp

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
-include $(NUHMSSM_SUSY_BETAS_MK)
-include $(NUHMSSM_SOFT_BETAS_MK)
-include $(NUHMSSM_CXXQFT_VERTICES_MK)
-include $(NUHMSSM_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(NUHMSSM_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(NUHMSSM_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(NUHMSSM_CXXQFT_VERTICES_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(NUHMSSM_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBNUHMSSM_SRC := $(sort $(LIBNUHMSSM_SRC))
EXENUHMSSM_SRC := $(sort $(EXENUHMSSM_SRC))

LIBNUHMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBNUHMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBNUHMSSM_SRC)))

EXENUHMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXENUHMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXENUHMSSM_SRC)))

EXENUHMSSM_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXENUHMSSM_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXENUHMSSM_SRC)))

LIBNUHMSSM_DEP := \
		$(LIBNUHMSSM_OBJ:.o=.d)

EXENUHMSSM_DEP := \
		$(EXENUHMSSM_OBJ:.o=.d)

LLNUHMSSM_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLNUHMSSM_SRC)))

LLNUHMSSM_OBJ  := $(LLNUHMSSM_SRC:.cpp=.o)
LLNUHMSSM_LIB  := $(LLNUHMSSM_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBNUHMSSM     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_NUHMSSM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_NUHMSSM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL_DIR) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBNUHMSSM) $(EXENUHMSSM_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(NUHMSSM_INSTALL_DIR)
		$(Q)install -d $(NUHMSSM_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBNUHMSSM_SRC) $(NUHMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBNUHMSSM_CXXQFT_VERTICES_SRC) $(NUHMSSM_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBNUHMSSM_HDR) $(NUHMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBNUHMSSM_CXXQFT_HDR) $(NUHMSSM_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(EXENUHMSSM_SRC) $(NUHMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLNUHMSSM_SRC) $(NUHMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLNUHMSSM_MMA) $(NUHMSSM_INSTALL_DIR)
		$(Q)$(INSTALL_STRIPPED) $(NUHMSSM_MK) $(NUHMSSM_INSTALL_DIR) -m u=rw,g=r,o=r
		$(Q)install -m u=rw,g=r,o=r $(NUHMSSM_INCLUDE_MK) $(NUHMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(NUHMSSM_CXXQFT_VERTICES_MK) $(NUHMSSM_INSTALL_CXXQFT_DIR)

ifneq ($(NUHMSSM_SLHA_INPUT),)
		$(Q)install -m u=rw,g=r,o=r $(NUHMSSM_SLHA_INPUT) $(NUHMSSM_INSTALL_DIR)
endif
		$(Q)install -m u=rw,g=r,o=r $(NUHMSSM_REFERENCES) $(NUHMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(NUHMSSM_GNUPLOT) $(NUHMSSM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		$(Q)-rm -f $(LIBNUHMSSM_DEP)
		$(Q)-rm -f $(EXENUHMSSM_DEP)
		$(Q)-rm -f $(LLNUHMSSM_DEP)

clean-$(MODNAME)-lib:
		$(Q)-rm -f $(LIBNUHMSSM)
		$(Q)-rm -f $(LLNUHMSSM_LIB)

clean-$(MODNAME)-obj:
		$(Q)-rm -f $(LIBNUHMSSM_OBJ)
		$(Q)-rm -f $(EXENUHMSSM_OBJ)
		$(Q)-rm -f $(LLNUHMSSM_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		$(Q)-rm -f $(LIBNUHMSSM_SRC)
		$(Q)-rm -f $(LIBNUHMSSM_HDR)
		$(Q)-rm -f $(LIBNUHMSSM_CXXQFT_HDR)
		$(Q)-rm -f $(EXENUHMSSM_SRC)
		$(Q)-rm -f $(LLNUHMSSM_SRC)
		$(Q)-rm -f $(LLNUHMSSM_MMA)
		$(Q)-rm -f $(METACODE_STAMP_NUHMSSM)
		$(Q)-rm -f $(NUHMSSM_INCLUDE_MK)
		$(Q)-rm -f $(NUHMSSM_CXXQFT_VERTICES_MK)
		$(Q)-rm -f $(NUHMSSM_SLHA_INPUT)
		$(Q)-rm -f $(NUHMSSM_REFERENCES)
		$(Q)-rm -f $(NUHMSSM_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		$(Q)-rm -f $(EXENUHMSSM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		$(Q)tar -czf $(NUHMSSM_TARBALL) \
		$(LIBNUHMSSM_SRC) $(LIBNUHMSSM_HDR) $(LIBNUHMSSM_CXXQFT_HDR) \
		$(EXENUHMSSM_SRC) \
		$(LLNUHMSSM_SRC) $(LLNUHMSSM_MMA) \
		$(NUHMSSM_MK) $(NUHMSSM_INCLUDE_MK) $(NUHMSSM_CXXQFT_VERTICES_MK) \
		$(NUHMSSM_SLHA_INPUT) $(NUHMSSM_REFERENCES) \
		$(NUHMSSM_GNUPLOT)

$(LIBNUHMSSM_SRC) $(LIBNUHMSSM_HDR) $(LIBNUHMSSM_CXXQFT_HDR) $(EXENUHMSSM_SRC) $(LLNUHMSSM_SRC) $(LLNUHMSSM_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_NUHMSSM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_NUHMSSM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_NUHMSSM)
		@$(MSG)
		$(Q)printf "%s" "Get[\"$<\"]; Quit[]" | "$(MATH)" || (echo "Error: The code generation failed!"; exit 1)
		@touch "$(METACODE_STAMP_NUHMSSM)"
		@echo "Note: to regenerate NUHMSSM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_NUHMSSM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_NUHMSSM):
		@true
endif

$(LIBNUHMSSM_DEP) $(EXENUHMSSM_DEP) $(LLNUHMSSM_DEP) $(LIBNUHMSSM_OBJ) $(EXENUHMSSM_OBJ) $(LLNUHMSSM_OBJ) $(LLNUHMSSM_LIB): \
	CPPFLAGS += $(MODNUHMSSM_SUBMOD_INC) $(MODNUHMSSM_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS)  $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBNUHMSSM_DEP) $(EXENUHMSSM_DEP) $(LLNUHMSSM_DEP) $(LIBNUHMSSM_OBJ) $(EXENUHMSSM_OBJ) $(LLNUHMSSM_OBJ) $(LLNUHMSSM_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLNUHMSSM_OBJ) $(LLNUHMSSM_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBNUHMSSM): $(LIBNUHMSSM_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBNUHMSSM) $(MODNUHMSSM_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(FUTILIBS)
		@$(MSG)
		$(Q)$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(SQLITELIBS) $(TSILLIBS) $(FLIBS) $(THREADLIBS) $(LDLIBS) $(FUTILIBS)

$(LLNUHMSSM_LIB): $(LLNUHMSSM_OBJ) $(LIBNUHMSSM) $(MODNUHMSSM_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS)) $(FUTILIBS)
		@$(MSG)
		$(Q)$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS) $(FUTILIBS)

ALLDEP += $(LIBNUHMSSM_DEP) $(EXENUHMSSM_DEP)
ALLSRC += $(LIBNUHMSSM_SRC) $(EXENUHMSSM_SRC)
ALLLIB += $(LIBNUHMSSM)
ALLEXE += $(EXENUHMSSM_EXE)
ALLMODDEP += $(MODNUHMSSM_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLNUHMSSM_DEP)
ALLSRC += $(LLNUHMSSM_SRC)
ALLLL  += $(LLNUHMSSM_LIB)
endif
