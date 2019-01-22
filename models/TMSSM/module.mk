DIR          := models/TMSSM
MODNAME      := TMSSM
SARAH_MODEL  := TMSSM
WITH_$(MODNAME) := yes
MODTMSSM_MOD := SM
MODTMSSM_DEP := $(patsubst %,model_specific/%,$(MODTMSSM_MOD))
MODTMSSM_INC := $(patsubst %,-Imodel_specific/%,$(MODTMSSM_MOD))
MODTMSSM_LIB := $(foreach M,$(MODTMSSM_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

MODTMSSM_SUBMOD  := $(DIR)/cxx_qft
MODTMSSM_SUBMOD_INC := $(patsubst %,-I%,$(MODTMSSM_SUBMOD))

TMSSM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)
TMSSM_INSTALL_CXXQFT_DIR := \
		$(TMSSM_INSTALL_DIR)/cxx_qft

TMSSM_MK     := \
		$(DIR)/module.mk

TMSSM_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

TMSSM_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

TMSSM_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

TMSSM_INCLUDE_MK := \
		$(TMSSM_SUSY_BETAS_MK) \
		$(TMSSM_SOFT_BETAS_MK)

TMSSM_SLHA_INPUT := \
		$(DIR)/LesHouches.in.TMSSM_generated \
		$(DIR)/LesHouches.in.TMSSM

TMSSM_REFERENCES := \
		$(DIR)/TMSSM_references.tex

TMSSM_GNUPLOT := \
		$(DIR)/TMSSM_plot_rgflow.gnuplot \
		$(DIR)/TMSSM_plot_spectrum.gnuplot

TMSSM_TARBALL := \
		$(MODNAME).tar.gz

LIBTMSSM_SRC := \
		$(DIR)/TMSSM_a_muon.cpp \
		$(DIR)/TMSSM_edm.cpp \
		$(DIR)/TMSSM_effective_couplings.cpp \
		$(DIR)/TMSSM_info.cpp \
		$(DIR)/TMSSM_input_parameters.cpp \
		$(DIR)/TMSSM_mass_eigenstates.cpp \
		$(DIR)/TMSSM_observables.cpp \
		$(DIR)/TMSSM_physical.cpp \
		$(DIR)/TMSSM_slha_io.cpp \
		$(DIR)/TMSSM_soft_parameters.cpp \
		$(DIR)/TMSSM_susy_parameters.cpp \
		$(DIR)/TMSSM_utilities.cpp \
		$(DIR)/TMSSM_weinberg_angle.cpp

EXETMSSM_SRC := \
		$(DIR)/run_TMSSM.cpp \
		$(DIR)/run_cmd_line_TMSSM.cpp \
		$(DIR)/scan_TMSSM.cpp
LLTMSSM_LIB  :=
LLTMSSM_OBJ  :=
LLTMSSM_SRC  := \
		$(DIR)/TMSSM_librarylink.cpp

LLTMSSM_MMA  := \
		$(DIR)/TMSSM_librarylink.m \
		$(DIR)/run_TMSSM.m

LIBTMSSM_HDR := \
		$(DIR)/TMSSM_a_muon.hpp \
		$(DIR)/TMSSM_convergence_tester.hpp \
		$(DIR)/TMSSM_edm.hpp \
		$(DIR)/TMSSM_effective_couplings.hpp \
		$(DIR)/TMSSM_ewsb_solver.hpp \
		$(DIR)/TMSSM_ewsb_solver_interface.hpp \
		$(DIR)/TMSSM_high_scale_constraint.hpp \
		$(DIR)/TMSSM_info.hpp \
		$(DIR)/TMSSM_initial_guesser.hpp \
		$(DIR)/TMSSM_input_parameters.hpp \
		$(DIR)/TMSSM_low_scale_constraint.hpp \
		$(DIR)/TMSSM_mass_eigenstates.hpp \
		$(DIR)/TMSSM_model.hpp \
		$(DIR)/TMSSM_model_slha.hpp \
		$(DIR)/TMSSM_observables.hpp \
		$(DIR)/TMSSM_physical.hpp \
		$(DIR)/TMSSM_slha_io.hpp \
		$(DIR)/TMSSM_spectrum_generator.hpp \
		$(DIR)/TMSSM_spectrum_generator_interface.hpp \
		$(DIR)/TMSSM_soft_parameters.hpp \
		$(DIR)/TMSSM_susy_parameters.hpp \
		$(DIR)/TMSSM_susy_scale_constraint.hpp \
		$(DIR)/TMSSM_utilities.hpp \
		$(DIR)/TMSSM_weinberg_angle.hpp

LIBTMSSM_CXXQFT_HDR := \
		$(DIR)/cxx_qft/TMSSM_qft.hpp \
		$(DIR)/cxx_qft/TMSSM_fields.hpp \
		$(DIR)/cxx_qft/TMSSM_vertices.hpp \
		$(DIR)/cxx_qft/TMSSM_context_base.hpp \
		$(DIR)/cxx_qft/TMSSM_npointfunctions.hpp

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
-include $(TMSSM_SUSY_BETAS_MK)
-include $(TMSSM_SOFT_BETAS_MK)
-include $(TMSSM_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(TMSSM_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(TMSSM_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(TMSSM_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBTMSSM_SRC := $(sort $(LIBTMSSM_SRC))
EXETMSSM_SRC := $(sort $(EXETMSSM_SRC))

LIBTMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBTMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBTMSSM_SRC)))

EXETMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXETMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXETMSSM_SRC)))

EXETMSSM_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXETMSSM_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXETMSSM_SRC)))

LIBTMSSM_DEP := \
		$(LIBTMSSM_OBJ:.o=.d)

EXETMSSM_DEP := \
		$(EXETMSSM_OBJ:.o=.d)

LLTMSSM_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLTMSSM_SRC)))

LLTMSSM_OBJ  := $(LLTMSSM_SRC:.cpp=.o)
LLTMSSM_LIB  := $(LLTMSSM_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBTMSSM     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_TMSSM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_TMSSM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL_DIR) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBTMSSM) $(EXETMSSM_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(TMSSM_INSTALL_DIR)
		install -d $(TMSSM_INSTALL_CXXQFT_DIR)
		install -m u=rw,g=r,o=r $(LIBTMSSM_SRC) $(TMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBTMSSM_HDR) $(TMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBTMSSM_CXXQFT_HDR) $(TMSSM_INSTALL_CXXQFT_DIR)
		install -m u=rw,g=r,o=r $(EXETMSSM_SRC) $(TMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLTMSSM_SRC) $(TMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLTMSSM_MMA) $(TMSSM_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(TMSSM_MK) $(TMSSM_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(TMSSM_INCLUDE_MK) $(TMSSM_INSTALL_DIR)
ifneq ($(TMSSM_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(TMSSM_SLHA_INPUT) $(TMSSM_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(TMSSM_REFERENCES) $(TMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(TMSSM_GNUPLOT) $(TMSSM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBTMSSM_DEP)
		-rm -f $(EXETMSSM_DEP)
		-rm -f $(LLTMSSM_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBTMSSM)
		-rm -f $(LLTMSSM_LIB)

clean-$(MODNAME)-obj:
		-rm -f $(LIBTMSSM_OBJ)
		-rm -f $(EXETMSSM_OBJ)
		-rm -f $(LLTMSSM_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		-rm -f $(LIBTMSSM_SRC)
		-rm -f $(LIBTMSSM_HDR)
		-rm -f $(LIBTMSSM_CXXQFT_HDR)
		-rm -f $(EXETMSSM_SRC)
		-rm -f $(LLTMSSM_SRC)
		-rm -f $(LLTMSSM_MMA)
		-rm -f $(METACODE_STAMP_TMSSM)
		-rm -f $(TMSSM_INCLUDE_MK)
		-rm -f $(TMSSM_SLHA_INPUT)
		-rm -f $(TMSSM_REFERENCES)
		-rm -f $(TMSSM_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		-rm -f $(EXETMSSM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(TMSSM_TARBALL) \
		$(LIBTMSSM_SRC) $(LIBTMSSM_HDR) $(LIBTMSSM_CXXQFT_HDR) \
		$(EXETMSSM_SRC) \
		$(LLTMSSM_SRC) $(LLTMSSM_MMA) \
		$(TMSSM_MK) $(TMSSM_INCLUDE_MK) \
		$(TMSSM_SLHA_INPUT) $(TMSSM_REFERENCES) \
		$(TMSSM_GNUPLOT)

$(LIBTMSSM_SRC) $(LIBTMSSM_HDR) $(LIBTMSSM_CXXQFT_HDR) $(EXETMSSM_SRC) $(LLTMSSM_SRC) $(LLTMSSM_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_TMSSM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_TMSSM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_TMSSM)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]" || (echo "Error: The code generation failed!"; exit 1)
		@touch "$(METACODE_STAMP_TMSSM)"
		@echo "Note: to regenerate TMSSM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_TMSSM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_TMSSM):
		@true
endif

$(LIBTMSSM_DEP) $(EXETMSSM_DEP) $(LLTMSSM_DEP) $(LIBTMSSM_OBJ) $(EXETMSSM_OBJ) $(LLTMSSM_OBJ) $(LLTMSSM_LIB): \
	CPPFLAGS += $(MODTMSSM_SUBMOD_INC) $(MODTMSSM_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBTMSSM_DEP) $(EXETMSSM_DEP) $(LLTMSSM_DEP) $(LIBTMSSM_OBJ) $(EXETMSSM_OBJ) $(LLTMSSM_OBJ) $(LLTMSSM_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLTMSSM_OBJ) $(LLTMSSM_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBTMSSM): $(LIBTMSSM_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBTMSSM) $(MODTMSSM_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLTMSSM_LIB): $(LLTMSSM_OBJ) $(LIBTMSSM) $(MODTMSSM_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS)

ALLDEP += $(LIBTMSSM_DEP) $(EXETMSSM_DEP)
ALLSRC += $(LIBTMSSM_SRC) $(EXETMSSM_SRC)
ALLLIB += $(LIBTMSSM)
ALLEXE += $(EXETMSSM_EXE)
ALLMODDEP += $(MODTMSSM_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLTMSSM_DEP)
ALLSRC += $(LLTMSSM_SRC)
ALLLL  += $(LLTMSSM_LIB)
endif
