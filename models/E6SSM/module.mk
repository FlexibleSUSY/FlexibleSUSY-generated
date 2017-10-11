DIR          := models/E6SSM
MODNAME      := E6SSM
SARAH_MODEL  := E6SSM
WITH_$(MODNAME) := yes

E6SSM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

E6SSM_MK     := \
		$(DIR)/module.mk

E6SSM_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

E6SSM_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

E6SSM_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

E6SSM_INCLUDE_MK := \
		$(E6SSM_SUSY_BETAS_MK) \
		$(E6SSM_SOFT_BETAS_MK)

E6SSM_SLHA_INPUT := \
		$(DIR)/LesHouches.in.E6SSM_generated \
		$(DIR)/LesHouches.in.E6SSM

E6SSM_GNUPLOT := \
		$(DIR)/E6SSM_plot_rgflow.gnuplot \
		$(DIR)/E6SSM_plot_spectrum.gnuplot

E6SSM_TARBALL := \
		$(MODNAME).tar.gz

LIBE6SSM_SRC := \
		$(DIR)/E6SSM_a_muon.cpp \
		$(DIR)/E6SSM_edm.cpp \
		$(DIR)/E6SSM_effective_couplings.cpp \
		$(DIR)/E6SSM_info.cpp \
		$(DIR)/E6SSM_input_parameters.cpp \
		$(DIR)/E6SSM_mass_eigenstates.cpp \
		$(DIR)/E6SSM_observables.cpp \
		$(DIR)/E6SSM_physical.cpp \
		$(DIR)/E6SSM_slha_io.cpp \
		$(DIR)/E6SSM_soft_parameters.cpp \
		$(DIR)/E6SSM_susy_parameters.cpp \
		$(DIR)/E6SSM_utilities.cpp \
		$(DIR)/E6SSM_weinberg_angle.cpp

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
		$(DIR)/E6SSM_cxx_diagrams.hpp \
		$(DIR)/E6SSM_a_muon.hpp \
		$(DIR)/E6SSM_convergence_tester.hpp \
		$(DIR)/E6SSM_edm.hpp \
		$(DIR)/E6SSM_effective_couplings.hpp \
		$(DIR)/E6SSM_ewsb_solver.hpp \
		$(DIR)/E6SSM_ewsb_solver_interface.hpp \
		$(DIR)/E6SSM_high_scale_constraint.hpp \
		$(DIR)/E6SSM_info.hpp \
		$(DIR)/E6SSM_initial_guesser.hpp \
		$(DIR)/E6SSM_input_parameters.hpp \
		$(DIR)/E6SSM_low_scale_constraint.hpp \
		$(DIR)/E6SSM_mass_eigenstates.hpp \
		$(DIR)/E6SSM_model.hpp \
		$(DIR)/E6SSM_model_slha.hpp \
		$(DIR)/E6SSM_observables.hpp \
		$(DIR)/E6SSM_physical.hpp \
		$(DIR)/E6SSM_slha_io.hpp \
		$(DIR)/E6SSM_spectrum_generator.hpp \
		$(DIR)/E6SSM_spectrum_generator_interface.hpp \
		$(DIR)/E6SSM_soft_parameters.hpp \
		$(DIR)/E6SSM_susy_parameters.hpp \
		$(DIR)/E6SSM_susy_scale_constraint.hpp \
		$(DIR)/E6SSM_utilities.hpp \
		$(DIR)/E6SSM_weinberg_angle.hpp

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
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBE6SSM) $(EXEE6SSM_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(E6SSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBE6SSM_SRC) $(E6SSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBE6SSM_HDR) $(E6SSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXEE6SSM_SRC) $(E6SSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLE6SSM_SRC) $(E6SSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLE6SSM_MMA) $(E6SSM_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(E6SSM_MK) $(E6SSM_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(E6SSM_INCLUDE_MK) $(E6SSM_INSTALL_DIR)
ifneq ($(E6SSM_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(E6SSM_SLHA_INPUT) $(E6SSM_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(E6SSM_GNUPLOT) $(E6SSM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBE6SSM_DEP)
		-rm -f $(EXEE6SSM_DEP)
		-rm -f $(LLE6SSM_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBE6SSM)
		-rm -f $(LLE6SSM_LIB)

clean-$(MODNAME)-obj:
		-rm -f $(LIBE6SSM_OBJ)
		-rm -f $(EXEE6SSM_OBJ)
		-rm -f $(LLE6SSM_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		-rm -f $(LIBE6SSM_SRC)
		-rm -f $(LIBE6SSM_HDR)
		-rm -f $(EXEE6SSM_SRC)
		-rm -f $(LLE6SSM_SRC)
		-rm -f $(LLE6SSM_MMA)
		-rm -f $(METACODE_STAMP_E6SSM)
		-rm -f $(E6SSM_INCLUDE_MK)
		-rm -f $(E6SSM_SLHA_INPUT)
		-rm -f $(E6SSM_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		-rm -f $(EXEE6SSM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(E6SSM_TARBALL) \
		$(LIBE6SSM_SRC) $(LIBE6SSM_HDR) \
		$(EXEE6SSM_SRC) \
		$(LLE6SSM_SRC) $(LLE6SSM_MMA) \
		$(E6SSM_MK) $(E6SSM_INCLUDE_MK) \
		$(E6SSM_SLHA_INPUT) $(E6SSM_GNUPLOT)

$(LIBE6SSM_SRC) $(LIBE6SSM_HDR) $(EXEE6SSM_SRC) $(LLE6SSM_SRC) $(LLE6SSM_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_E6SSM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_E6SSM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_E6SSM)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
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
	CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBE6SSM_DEP) $(EXEE6SSM_DEP) $(LLE6SSM_DEP) $(LIBE6SSM_OBJ) $(EXEE6SSM_OBJ) $(LLE6SSM_OBJ) $(LLE6SSM_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLE6SSM_OBJ) $(LLE6SSM_LIB): \
	CPPFLAGS += $(shell $(MATH_INC_PATHS) --math-cmd="$(MATH)" -I --librarylink --mathlink)

$(LIBE6SSM): $(LIBE6SSM_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBE6SSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLE6SSM_LIB): $(LLE6SSM_OBJ) $(LIBE6SSM) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

ALLDEP += $(LIBE6SSM_DEP) $(EXEE6SSM_DEP)
ALLSRC += $(LIBE6SSM_SRC) $(EXEE6SSM_SRC)
ALLLIB += $(LIBE6SSM)
ALLEXE += $(EXEE6SSM_EXE)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLE6SSM_DEP)
ALLSRC += $(LLE6SSM_SRC)
ALLLL  += $(LLE6SSM_LIB)
endif
