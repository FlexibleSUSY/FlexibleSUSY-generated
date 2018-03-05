DIR          := models/THDMII
MODNAME      := THDMII
SARAH_MODEL  := THDM-II
WITH_$(MODNAME) := yes
MODTHDMII_MOD := SM
MODTHDMII_DEP := $(patsubst %,model_specific/%,$(MODTHDMII_MOD))
MODTHDMII_INC := $(patsubst %,-Imodel_specific/%,$(MODTHDMII_MOD))
MODTHDMII_LIB := $(foreach M,$(MODTHDMII_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

THDMII_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

THDMII_MK     := \
		$(DIR)/module.mk

THDMII_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

THDMII_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

THDMII_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

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
		$(DIR)/THDMII_a_muon.cpp \
		$(DIR)/THDMII_edm.cpp \
		$(DIR)/THDMII_effective_couplings.cpp \
		$(DIR)/THDMII_info.cpp \
		$(DIR)/THDMII_input_parameters.cpp \
		$(DIR)/THDMII_mass_eigenstates.cpp \
		$(DIR)/THDMII_observables.cpp \
		$(DIR)/THDMII_physical.cpp \
		$(DIR)/THDMII_slha_io.cpp \
		$(DIR)/THDMII_soft_parameters.cpp \
		$(DIR)/THDMII_susy_parameters.cpp \
		$(DIR)/THDMII_utilities.cpp \
		$(DIR)/THDMII_weinberg_angle.cpp

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
		$(DIR)/THDMII_cxx_diagrams.hpp \
		$(DIR)/THDMII_a_muon.hpp \
		$(DIR)/THDMII_convergence_tester.hpp \
		$(DIR)/THDMII_edm.hpp \
		$(DIR)/THDMII_effective_couplings.hpp \
		$(DIR)/THDMII_ewsb_solver.hpp \
		$(DIR)/THDMII_ewsb_solver_interface.hpp \
		$(DIR)/THDMII_high_scale_constraint.hpp \
		$(DIR)/THDMII_info.hpp \
		$(DIR)/THDMII_initial_guesser.hpp \
		$(DIR)/THDMII_input_parameters.hpp \
		$(DIR)/THDMII_low_scale_constraint.hpp \
		$(DIR)/THDMII_mass_eigenstates.hpp \
		$(DIR)/THDMII_model.hpp \
		$(DIR)/THDMII_model_slha.hpp \
		$(DIR)/THDMII_observables.hpp \
		$(DIR)/THDMII_physical.hpp \
		$(DIR)/THDMII_slha_io.hpp \
		$(DIR)/THDMII_spectrum_generator.hpp \
		$(DIR)/THDMII_spectrum_generator_interface.hpp \
		$(DIR)/THDMII_soft_parameters.hpp \
		$(DIR)/THDMII_susy_parameters.hpp \
		$(DIR)/THDMII_susy_scale_constraint.hpp \
		$(DIR)/THDMII_utilities.hpp \
		$(DIR)/THDMII_weinberg_angle.hpp

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
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBTHDMII) $(EXETHDMII_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(THDMII_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBTHDMII_SRC) $(THDMII_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBTHDMII_HDR) $(THDMII_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXETHDMII_SRC) $(THDMII_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLTHDMII_SRC) $(THDMII_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLTHDMII_MMA) $(THDMII_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(THDMII_MK) $(THDMII_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(THDMII_INCLUDE_MK) $(THDMII_INSTALL_DIR)
ifneq ($(THDMII_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(THDMII_SLHA_INPUT) $(THDMII_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(THDMII_REFERENCES) $(THDMII_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(THDMII_GNUPLOT) $(THDMII_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBTHDMII_DEP)
		-rm -f $(EXETHDMII_DEP)
		-rm -f $(LLTHDMII_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBTHDMII)
		-rm -f $(LLTHDMII_LIB)

clean-$(MODNAME)-obj:
		-rm -f $(LIBTHDMII_OBJ)
		-rm -f $(EXETHDMII_OBJ)
		-rm -f $(LLTHDMII_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		-rm -f $(LIBTHDMII_SRC)
		-rm -f $(LIBTHDMII_HDR)
		-rm -f $(EXETHDMII_SRC)
		-rm -f $(LLTHDMII_SRC)
		-rm -f $(LLTHDMII_MMA)
		-rm -f $(METACODE_STAMP_THDMII)
		-rm -f $(THDMII_INCLUDE_MK)
		-rm -f $(THDMII_SLHA_INPUT)
		-rm -f $(THDMII_REFERENCES)
		-rm -f $(THDMII_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		-rm -f $(EXETHDMII_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(THDMII_TARBALL) \
		$(LIBTHDMII_SRC) $(LIBTHDMII_HDR) \
		$(EXETHDMII_SRC) \
		$(LLTHDMII_SRC) $(LLTHDMII_MMA) \
		$(THDMII_MK) $(THDMII_INCLUDE_MK) \
		$(THDMII_SLHA_INPUT) $(THDMII_REFERENCES) \
		$(THDMII_GNUPLOT)

$(LIBTHDMII_SRC) $(LIBTHDMII_HDR) $(EXETHDMII_SRC) $(LLTHDMII_SRC) $(LLTHDMII_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_THDMII)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_THDMII): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_THDMII)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
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
	CPPFLAGS += $(MODTHDMII_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBTHDMII_DEP) $(EXETHDMII_DEP) $(LLTHDMII_DEP) $(LIBTHDMII_OBJ) $(EXETHDMII_OBJ) $(LLTHDMII_OBJ) $(LLTHDMII_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLTHDMII_OBJ) $(LLTHDMII_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBTHDMII): $(LIBTHDMII_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBTHDMII) $(MODTHDMII_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLTHDMII_LIB): $(LLTHDMII_OBJ) $(LIBTHDMII) $(MODTHDMII_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS)

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
