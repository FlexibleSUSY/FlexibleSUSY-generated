DIR          := models/UMSSM
MODNAME      := UMSSM
SARAH_MODEL  := UMSSM
WITH_$(MODNAME) := yes
MODUMSSM_MOD := SM MSSM_higgs NMSSM_higgs
MODUMSSM_DEP := $(patsubst %,model_specific/%,$(MODUMSSM_MOD))
MODUMSSM_INC := $(patsubst %,-Imodel_specific/%,$(MODUMSSM_MOD))
MODUMSSM_LIB := $(foreach M,$(MODUMSSM_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

MODUMSSM_SUBMOD  := $(DIR)/cxx_qft
MODUMSSM_SUBMOD_INC := $(patsubst %,-I%,$(MODUMSSM_SUBMOD))

UMSSM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)
UMSSM_INSTALL_CXXQFT_DIR := \
		$(UMSSM_INSTALL_DIR)/cxx_qft

UMSSM_MK     := \
		$(DIR)/module.mk

UMSSM_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

UMSSM_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

UMSSM_CXX_QFT_VERTICES_MK := \
		$(DIR)/cxx_qft/vertices.mk

UMSSM_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

UMSSM_INCLUDE_MK := \
		$(UMSSM_SUSY_BETAS_MK) \
		$(UMSSM_SOFT_BETAS_MK) \
		$(UMSSM_CXX_QFT_VERTICES_MK)

UMSSM_SLHA_INPUT := \
		$(DIR)/LesHouches.in.UMSSM_generated \
		$(DIR)/LesHouches.in.UMSSM

UMSSM_REFERENCES := \
		$(DIR)/UMSSM_references.tex

UMSSM_GNUPLOT := \
		$(DIR)/UMSSM_plot_rgflow.gnuplot \
		$(DIR)/UMSSM_plot_spectrum.gnuplot

UMSSM_TARBALL := \
		$(MODNAME).tar.gz

LIBUMSSM_SRC := \
		$(DIR)/UMSSM_a_muon.cpp \
		$(DIR)/UMSSM_edm.cpp \
		$(DIR)/UMSSM_FFV_form_factors.cpp \
		$(DIR)/UMSSM_l_to_lgamma.cpp \
		$(DIR)/UMSSM_effective_couplings.cpp \
		$(DIR)/UMSSM_info.cpp \
		$(DIR)/UMSSM_input_parameters.cpp \
		$(DIR)/UMSSM_mass_eigenstates.cpp \
		$(DIR)/UMSSM_observables.cpp \
		$(DIR)/UMSSM_physical.cpp \
		$(DIR)/UMSSM_slha_io.cpp \
		$(DIR)/UMSSM_soft_parameters.cpp \
		$(DIR)/UMSSM_susy_parameters.cpp \
		$(DIR)/UMSSM_utilities.cpp \
		$(DIR)/UMSSM_weinberg_angle.cpp

EXEUMSSM_SRC := \
		$(DIR)/run_UMSSM.cpp \
		$(DIR)/run_cmd_line_UMSSM.cpp \
		$(DIR)/scan_UMSSM.cpp
LLUMSSM_LIB  :=
LLUMSSM_OBJ  :=
LLUMSSM_SRC  := \
		$(DIR)/UMSSM_librarylink.cpp

LLUMSSM_MMA  := \
		$(DIR)/UMSSM_librarylink.m \
		$(DIR)/run_UMSSM.m

LIBUMSSM_HDR := \
		$(DIR)/UMSSM_a_muon.hpp \
		$(DIR)/UMSSM_convergence_tester.hpp \
		$(DIR)/UMSSM_edm.hpp \
		$(DIR)/UMSSM_FFV_form_factors.hpp \
		$(DIR)/UMSSM_l_to_lgamma.hpp \
		$(DIR)/UMSSM_effective_couplings.hpp \
		$(DIR)/UMSSM_ewsb_solver.hpp \
		$(DIR)/UMSSM_ewsb_solver_interface.hpp \
		$(DIR)/UMSSM_high_scale_constraint.hpp \
		$(DIR)/UMSSM_info.hpp \
		$(DIR)/UMSSM_initial_guesser.hpp \
		$(DIR)/UMSSM_input_parameters.hpp \
		$(DIR)/UMSSM_low_scale_constraint.hpp \
		$(DIR)/UMSSM_mass_eigenstates.hpp \
		$(DIR)/UMSSM_model.hpp \
		$(DIR)/UMSSM_model_slha.hpp \
		$(DIR)/UMSSM_observables.hpp \
		$(DIR)/UMSSM_physical.hpp \
		$(DIR)/UMSSM_slha_io.hpp \
		$(DIR)/UMSSM_spectrum_generator.hpp \
		$(DIR)/UMSSM_spectrum_generator_interface.hpp \
		$(DIR)/UMSSM_soft_parameters.hpp \
		$(DIR)/UMSSM_susy_parameters.hpp \
		$(DIR)/UMSSM_susy_scale_constraint.hpp \
		$(DIR)/UMSSM_utilities.hpp \
		$(DIR)/UMSSM_weinberg_angle.hpp

LIBUMSSM_CXXQFT_HDR := \
		$(DIR)/cxx_qft/UMSSM_qft.hpp \
		$(DIR)/cxx_qft/UMSSM_fields.hpp \
		$(DIR)/cxx_qft/UMSSM_vertices.hpp \
		$(DIR)/cxx_qft/UMSSM_context_base.hpp \
		$(DIR)/cxx_qft/UMSSM_npointfunctions.hpp

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
-include $(UMSSM_SUSY_BETAS_MK)
-include $(UMSSM_SOFT_BETAS_MK)
-include $(UMSSM_CXX_QFT_VERTICES_MK)
-include $(UMSSM_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(UMSSM_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(UMSSM_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(UMSSM_CXX_QFT_VERTICES_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(UMSSM_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBUMSSM_SRC := $(sort $(LIBUMSSM_SRC))
EXEUMSSM_SRC := $(sort $(EXEUMSSM_SRC))

LIBUMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBUMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBUMSSM_SRC)))

EXEUMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEUMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEUMSSM_SRC)))

EXEUMSSM_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXEUMSSM_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXEUMSSM_SRC)))

LIBUMSSM_DEP := \
		$(LIBUMSSM_OBJ:.o=.d)

EXEUMSSM_DEP := \
		$(EXEUMSSM_OBJ:.o=.d)

LLUMSSM_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLUMSSM_SRC)))

LLUMSSM_OBJ  := $(LLUMSSM_SRC:.cpp=.o)
LLUMSSM_LIB  := $(LLUMSSM_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBUMSSM     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_UMSSM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_UMSSM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL_DIR) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBUMSSM) $(EXEUMSSM_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(UMSSM_INSTALL_DIR)
		$(Q)install -d $(UMSSM_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBUMSSM_SRC) $(UMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBUMSSM_HDR) $(UMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBUMSSM_CXXQFT_HDR) $(UMSSM_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(EXEUMSSM_SRC) $(UMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLUMSSM_SRC) $(UMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLUMSSM_MMA) $(UMSSM_INSTALL_DIR)
		$(Q)$(INSTALL_STRIPPED) $(UMSSM_MK) $(UMSSM_INSTALL_DIR) -m u=rw,g=r,o=r
		$(Q)install -m u=rw,g=r,o=r $(UMSSM_INCLUDE_MK) $(UMSSM_INSTALL_DIR)
ifneq ($(UMSSM_SLHA_INPUT),)
		$(Q)install -m u=rw,g=r,o=r $(UMSSM_SLHA_INPUT) $(UMSSM_INSTALL_DIR)
endif
		$(Q)install -m u=rw,g=r,o=r $(UMSSM_REFERENCES) $(UMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(UMSSM_GNUPLOT) $(UMSSM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		$(Q)-rm -f $(LIBUMSSM_DEP)
		$(Q)-rm -f $(EXEUMSSM_DEP)
		$(Q)-rm -f $(LLUMSSM_DEP)

clean-$(MODNAME)-lib:
		$(Q)-rm -f $(LIBUMSSM)
		$(Q)-rm -f $(LLUMSSM_LIB)

clean-$(MODNAME)-obj:
		$(Q)-rm -f $(LIBUMSSM_OBJ)
		$(Q)-rm -f $(EXEUMSSM_OBJ)
		$(Q)-rm -f $(LLUMSSM_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		$(Q)-rm -f $(LIBUMSSM_SRC)
		$(Q)-rm -f $(LIBUMSSM_HDR)
		$(Q)-rm -f $(LIBUMSSM_CXXQFT_HDR)
		$(Q)-rm -f $(EXEUMSSM_SRC)
		$(Q)-rm -f $(LLUMSSM_SRC)
		$(Q)-rm -f $(LLUMSSM_MMA)
		$(Q)-rm -f $(METACODE_STAMP_UMSSM)
		$(Q)-rm -f $(UMSSM_INCLUDE_MK)
		$(Q)-rm -f $(UMSSM_SLHA_INPUT)
		$(Q)-rm -f $(UMSSM_REFERENCES)
		$(Q)-rm -f $(UMSSM_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		$(Q)-rm -f $(EXEUMSSM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		$(Q)tar -czf $(UMSSM_TARBALL) \
		$(LIBUMSSM_SRC) $(LIBUMSSM_HDR) $(LIBUMSSM_CXXQFT_HDR) \
		$(EXEUMSSM_SRC) \
		$(LLUMSSM_SRC) $(LLUMSSM_MMA) \
		$(UMSSM_MK) $(UMSSM_INCLUDE_MK) \
		$(UMSSM_SLHA_INPUT) $(UMSSM_REFERENCES) \
		$(UMSSM_GNUPLOT)

$(LIBUMSSM_SRC) $(LIBUMSSM_HDR) $(LIBUMSSM_CXXQFT_HDR) $(EXEUMSSM_SRC) $(LLUMSSM_SRC) $(LLUMSSM_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_UMSSM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_UMSSM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_UMSSM)
		@$(MSG)
		$(Q)printf "%s" "Get[\"$<\"]; Quit[]" | "$(MATH)" || (echo "Error: The code generation failed!"; exit 1)
		@touch "$(METACODE_STAMP_UMSSM)"
		@echo "Note: to regenerate UMSSM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_UMSSM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_UMSSM):
		@true
endif

$(LIBUMSSM_DEP) $(EXEUMSSM_DEP) $(LLUMSSM_DEP) $(LIBUMSSM_OBJ) $(EXEUMSSM_OBJ) $(LLUMSSM_OBJ) $(LLUMSSM_LIB): \
	CPPFLAGS += $(MODUMSSM_SUBMOD_INC) $(MODUMSSM_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBUMSSM_DEP) $(EXEUMSSM_DEP) $(LLUMSSM_DEP) $(LIBUMSSM_OBJ) $(EXEUMSSM_OBJ) $(LLUMSSM_OBJ) $(LLUMSSM_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLUMSSM_OBJ) $(LLUMSSM_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBUMSSM): $(LIBUMSSM_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBUMSSM) $(MODUMSSM_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		@$(MSG)
		$(Q)$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLUMSSM_LIB): $(LLUMSSM_OBJ) $(LIBUMSSM) $(MODUMSSM_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		@$(MSG)
		$(Q)$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS)

ALLDEP += $(LIBUMSSM_DEP) $(EXEUMSSM_DEP)
ALLSRC += $(LIBUMSSM_SRC) $(EXEUMSSM_SRC)
ALLLIB += $(LIBUMSSM)
ALLEXE += $(EXEUMSSM_EXE)
ALLMODDEP += $(MODUMSSM_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLUMSSM_DEP)
ALLSRC += $(LLUMSSM_SRC)
ALLLL  += $(LLUMSSM_LIB)
endif
