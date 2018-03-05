DIR          := models/CE6SSM
MODNAME      := CE6SSM
SARAH_MODEL  := E6SSM
WITH_$(MODNAME) := yes
MODCE6SSM_MOD := SM MSSM_higgs NMSSM_higgs
MODCE6SSM_DEP := $(patsubst %,model_specific/%,$(MODCE6SSM_MOD))
MODCE6SSM_INC := $(patsubst %,-Imodel_specific/%,$(MODCE6SSM_MOD))
MODCE6SSM_LIB := $(foreach M,$(MODCE6SSM_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

CE6SSM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

CE6SSM_MK     := \
		$(DIR)/module.mk

CE6SSM_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

CE6SSM_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

CE6SSM_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

CE6SSM_INCLUDE_MK := \
		$(CE6SSM_SUSY_BETAS_MK) \
		$(CE6SSM_SOFT_BETAS_MK)

CE6SSM_SLHA_INPUT := \
		$(DIR)/LesHouches.in.CE6SSM_generated \
		$(DIR)/LesHouches.in.CE6SSM

CE6SSM_REFERENCES := \
		$(DIR)/CE6SSM_references.tex

CE6SSM_GNUPLOT := \
		$(DIR)/CE6SSM_plot_rgflow.gnuplot \
		$(DIR)/CE6SSM_plot_spectrum.gnuplot

CE6SSM_TARBALL := \
		$(MODNAME).tar.gz

LIBCE6SSM_SRC := \
		$(DIR)/CE6SSM_a_muon.cpp \
		$(DIR)/CE6SSM_edm.cpp \
		$(DIR)/CE6SSM_effective_couplings.cpp \
		$(DIR)/CE6SSM_info.cpp \
		$(DIR)/CE6SSM_input_parameters.cpp \
		$(DIR)/CE6SSM_mass_eigenstates.cpp \
		$(DIR)/CE6SSM_observables.cpp \
		$(DIR)/CE6SSM_physical.cpp \
		$(DIR)/CE6SSM_slha_io.cpp \
		$(DIR)/CE6SSM_soft_parameters.cpp \
		$(DIR)/CE6SSM_susy_parameters.cpp \
		$(DIR)/CE6SSM_utilities.cpp \
		$(DIR)/CE6SSM_weinberg_angle.cpp

EXECE6SSM_SRC := \
		$(DIR)/run_CE6SSM.cpp \
		$(DIR)/run_cmd_line_CE6SSM.cpp \
		$(DIR)/scan_CE6SSM.cpp
LLCE6SSM_LIB  :=
LLCE6SSM_OBJ  :=
LLCE6SSM_SRC  := \
		$(DIR)/CE6SSM_librarylink.cpp

LLCE6SSM_MMA  := \
		$(DIR)/CE6SSM_librarylink.m \
		$(DIR)/run_CE6SSM.m

LIBCE6SSM_HDR := \
		$(DIR)/CE6SSM_cxx_diagrams.hpp \
		$(DIR)/CE6SSM_a_muon.hpp \
		$(DIR)/CE6SSM_convergence_tester.hpp \
		$(DIR)/CE6SSM_edm.hpp \
		$(DIR)/CE6SSM_effective_couplings.hpp \
		$(DIR)/CE6SSM_ewsb_solver.hpp \
		$(DIR)/CE6SSM_ewsb_solver_interface.hpp \
		$(DIR)/CE6SSM_high_scale_constraint.hpp \
		$(DIR)/CE6SSM_info.hpp \
		$(DIR)/CE6SSM_initial_guesser.hpp \
		$(DIR)/CE6SSM_input_parameters.hpp \
		$(DIR)/CE6SSM_low_scale_constraint.hpp \
		$(DIR)/CE6SSM_mass_eigenstates.hpp \
		$(DIR)/CE6SSM_model.hpp \
		$(DIR)/CE6SSM_model_slha.hpp \
		$(DIR)/CE6SSM_observables.hpp \
		$(DIR)/CE6SSM_physical.hpp \
		$(DIR)/CE6SSM_slha_io.hpp \
		$(DIR)/CE6SSM_spectrum_generator.hpp \
		$(DIR)/CE6SSM_spectrum_generator_interface.hpp \
		$(DIR)/CE6SSM_soft_parameters.hpp \
		$(DIR)/CE6SSM_susy_parameters.hpp \
		$(DIR)/CE6SSM_susy_scale_constraint.hpp \
		$(DIR)/CE6SSM_utilities.hpp \
		$(DIR)/CE6SSM_weinberg_angle.hpp

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
-include $(CE6SSM_SUSY_BETAS_MK)
-include $(CE6SSM_SOFT_BETAS_MK)
-include $(CE6SSM_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(CE6SSM_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(CE6SSM_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(CE6SSM_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBCE6SSM_SRC := $(sort $(LIBCE6SSM_SRC))
EXECE6SSM_SRC := $(sort $(EXECE6SSM_SRC))

LIBCE6SSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBCE6SSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBCE6SSM_SRC)))

EXECE6SSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXECE6SSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXECE6SSM_SRC)))

EXECE6SSM_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXECE6SSM_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXECE6SSM_SRC)))

LIBCE6SSM_DEP := \
		$(LIBCE6SSM_OBJ:.o=.d)

EXECE6SSM_DEP := \
		$(EXECE6SSM_OBJ:.o=.d)

LLCE6SSM_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLCE6SSM_SRC)))

LLCE6SSM_OBJ  := $(LLCE6SSM_SRC:.cpp=.o)
LLCE6SSM_LIB  := $(LLCE6SSM_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBCE6SSM     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_CE6SSM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_CE6SSM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBCE6SSM) $(EXECE6SSM_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(CE6SSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBCE6SSM_SRC) $(CE6SSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBCE6SSM_HDR) $(CE6SSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXECE6SSM_SRC) $(CE6SSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLCE6SSM_SRC) $(CE6SSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLCE6SSM_MMA) $(CE6SSM_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(CE6SSM_MK) $(CE6SSM_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(CE6SSM_INCLUDE_MK) $(CE6SSM_INSTALL_DIR)
ifneq ($(CE6SSM_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(CE6SSM_SLHA_INPUT) $(CE6SSM_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(CE6SSM_REFERENCES) $(CE6SSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(CE6SSM_GNUPLOT) $(CE6SSM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBCE6SSM_DEP)
		-rm -f $(EXECE6SSM_DEP)
		-rm -f $(LLCE6SSM_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBCE6SSM)
		-rm -f $(LLCE6SSM_LIB)

clean-$(MODNAME)-obj:
		-rm -f $(LIBCE6SSM_OBJ)
		-rm -f $(EXECE6SSM_OBJ)
		-rm -f $(LLCE6SSM_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		-rm -f $(LIBCE6SSM_SRC)
		-rm -f $(LIBCE6SSM_HDR)
		-rm -f $(EXECE6SSM_SRC)
		-rm -f $(LLCE6SSM_SRC)
		-rm -f $(LLCE6SSM_MMA)
		-rm -f $(METACODE_STAMP_CE6SSM)
		-rm -f $(CE6SSM_INCLUDE_MK)
		-rm -f $(CE6SSM_SLHA_INPUT)
		-rm -f $(CE6SSM_REFERENCES)
		-rm -f $(CE6SSM_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		-rm -f $(EXECE6SSM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(CE6SSM_TARBALL) \
		$(LIBCE6SSM_SRC) $(LIBCE6SSM_HDR) \
		$(EXECE6SSM_SRC) \
		$(LLCE6SSM_SRC) $(LLCE6SSM_MMA) \
		$(CE6SSM_MK) $(CE6SSM_INCLUDE_MK) \
		$(CE6SSM_SLHA_INPUT) $(CE6SSM_REFERENCES) \
		$(CE6SSM_GNUPLOT)

$(LIBCE6SSM_SRC) $(LIBCE6SSM_HDR) $(EXECE6SSM_SRC) $(LLCE6SSM_SRC) $(LLCE6SSM_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_CE6SSM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_CE6SSM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_CE6SSM)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_CE6SSM)"
		@echo "Note: to regenerate CE6SSM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_CE6SSM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_CE6SSM):
		@true
endif

$(LIBCE6SSM_DEP) $(EXECE6SSM_DEP) $(LLCE6SSM_DEP) $(LIBCE6SSM_OBJ) $(EXECE6SSM_OBJ) $(LLCE6SSM_OBJ) $(LLCE6SSM_LIB): \
	CPPFLAGS += $(MODCE6SSM_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBCE6SSM_DEP) $(EXECE6SSM_DEP) $(LLCE6SSM_DEP) $(LIBCE6SSM_OBJ) $(EXECE6SSM_OBJ) $(LLCE6SSM_OBJ) $(LLCE6SSM_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLCE6SSM_OBJ) $(LLCE6SSM_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBCE6SSM): $(LIBCE6SSM_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBCE6SSM) $(MODCE6SSM_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLCE6SSM_LIB): $(LLCE6SSM_OBJ) $(LIBCE6SSM) $(MODCE6SSM_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS)

ALLDEP += $(LIBCE6SSM_DEP) $(EXECE6SSM_DEP)
ALLSRC += $(LIBCE6SSM_SRC) $(EXECE6SSM_SRC)
ALLLIB += $(LIBCE6SSM)
ALLEXE += $(EXECE6SSM_EXE)
ALLMODDEP += $(MODCE6SSM_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLCE6SSM_DEP)
ALLSRC += $(LLCE6SSM_SRC)
ALLLL  += $(LLCE6SSM_LIB)
endif
