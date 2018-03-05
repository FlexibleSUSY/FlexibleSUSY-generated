DIR          := models/THDMIIMSSMBC
MODNAME      := THDMIIMSSMBC
SARAH_MODEL  := THDM-II
WITH_$(MODNAME) := yes
MODTHDMIIMSSMBC_MOD := SM
MODTHDMIIMSSMBC_DEP := $(patsubst %,model_specific/%,$(MODTHDMIIMSSMBC_MOD))
MODTHDMIIMSSMBC_INC := $(patsubst %,-Imodel_specific/%,$(MODTHDMIIMSSMBC_MOD))
MODTHDMIIMSSMBC_LIB := $(foreach M,$(MODTHDMIIMSSMBC_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

THDMIIMSSMBC_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

THDMIIMSSMBC_MK     := \
		$(DIR)/module.mk

THDMIIMSSMBC_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

THDMIIMSSMBC_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

THDMIIMSSMBC_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

THDMIIMSSMBC_INCLUDE_MK := \
		$(THDMIIMSSMBC_SUSY_BETAS_MK) \
		$(THDMIIMSSMBC_SOFT_BETAS_MK)

THDMIIMSSMBC_SLHA_INPUT := \
		$(DIR)/LesHouches.in.THDMIIMSSMBC_generated \
		$(DIR)/LesHouches.in.THDMIIMSSMBC

THDMIIMSSMBC_REFERENCES := \
		$(DIR)/THDMIIMSSMBC_references.tex

THDMIIMSSMBC_GNUPLOT := \
		$(DIR)/THDMIIMSSMBC_plot_rgflow.gnuplot \
		$(DIR)/THDMIIMSSMBC_plot_spectrum.gnuplot

THDMIIMSSMBC_TARBALL := \
		$(MODNAME).tar.gz

LIBTHDMIIMSSMBC_SRC := \
		$(DIR)/THDMIIMSSMBC_a_muon.cpp \
		$(DIR)/THDMIIMSSMBC_edm.cpp \
		$(DIR)/THDMIIMSSMBC_effective_couplings.cpp \
		$(DIR)/THDMIIMSSMBC_info.cpp \
		$(DIR)/THDMIIMSSMBC_input_parameters.cpp \
		$(DIR)/THDMIIMSSMBC_mass_eigenstates.cpp \
		$(DIR)/THDMIIMSSMBC_observables.cpp \
		$(DIR)/THDMIIMSSMBC_physical.cpp \
		$(DIR)/THDMIIMSSMBC_slha_io.cpp \
		$(DIR)/THDMIIMSSMBC_soft_parameters.cpp \
		$(DIR)/THDMIIMSSMBC_susy_parameters.cpp \
		$(DIR)/THDMIIMSSMBC_utilities.cpp \
		$(DIR)/THDMIIMSSMBC_weinberg_angle.cpp

EXETHDMIIMSSMBC_SRC := \
		$(DIR)/run_THDMIIMSSMBC.cpp \
		$(DIR)/run_cmd_line_THDMIIMSSMBC.cpp \
		$(DIR)/scan_THDMIIMSSMBC.cpp
LLTHDMIIMSSMBC_LIB  :=
LLTHDMIIMSSMBC_OBJ  :=
LLTHDMIIMSSMBC_SRC  := \
		$(DIR)/THDMIIMSSMBC_librarylink.cpp

LLTHDMIIMSSMBC_MMA  := \
		$(DIR)/THDMIIMSSMBC_librarylink.m \
		$(DIR)/run_THDMIIMSSMBC.m

LIBTHDMIIMSSMBC_HDR := \
		$(DIR)/THDMIIMSSMBC_cxx_diagrams.hpp \
		$(DIR)/THDMIIMSSMBC_a_muon.hpp \
		$(DIR)/THDMIIMSSMBC_convergence_tester.hpp \
		$(DIR)/THDMIIMSSMBC_edm.hpp \
		$(DIR)/THDMIIMSSMBC_effective_couplings.hpp \
		$(DIR)/THDMIIMSSMBC_ewsb_solver.hpp \
		$(DIR)/THDMIIMSSMBC_ewsb_solver_interface.hpp \
		$(DIR)/THDMIIMSSMBC_high_scale_constraint.hpp \
		$(DIR)/THDMIIMSSMBC_info.hpp \
		$(DIR)/THDMIIMSSMBC_initial_guesser.hpp \
		$(DIR)/THDMIIMSSMBC_input_parameters.hpp \
		$(DIR)/THDMIIMSSMBC_low_scale_constraint.hpp \
		$(DIR)/THDMIIMSSMBC_mass_eigenstates.hpp \
		$(DIR)/THDMIIMSSMBC_model.hpp \
		$(DIR)/THDMIIMSSMBC_model_slha.hpp \
		$(DIR)/THDMIIMSSMBC_observables.hpp \
		$(DIR)/THDMIIMSSMBC_physical.hpp \
		$(DIR)/THDMIIMSSMBC_slha_io.hpp \
		$(DIR)/THDMIIMSSMBC_spectrum_generator.hpp \
		$(DIR)/THDMIIMSSMBC_spectrum_generator_interface.hpp \
		$(DIR)/THDMIIMSSMBC_soft_parameters.hpp \
		$(DIR)/THDMIIMSSMBC_susy_parameters.hpp \
		$(DIR)/THDMIIMSSMBC_susy_scale_constraint.hpp \
		$(DIR)/THDMIIMSSMBC_utilities.hpp \
		$(DIR)/THDMIIMSSMBC_weinberg_angle.hpp

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
-include $(THDMIIMSSMBC_SUSY_BETAS_MK)
-include $(THDMIIMSSMBC_SOFT_BETAS_MK)
-include $(THDMIIMSSMBC_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(THDMIIMSSMBC_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(THDMIIMSSMBC_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(THDMIIMSSMBC_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBTHDMIIMSSMBC_SRC := $(sort $(LIBTHDMIIMSSMBC_SRC))
EXETHDMIIMSSMBC_SRC := $(sort $(EXETHDMIIMSSMBC_SRC))

LIBTHDMIIMSSMBC_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBTHDMIIMSSMBC_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBTHDMIIMSSMBC_SRC)))

EXETHDMIIMSSMBC_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXETHDMIIMSSMBC_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXETHDMIIMSSMBC_SRC)))

EXETHDMIIMSSMBC_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXETHDMIIMSSMBC_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXETHDMIIMSSMBC_SRC)))

LIBTHDMIIMSSMBC_DEP := \
		$(LIBTHDMIIMSSMBC_OBJ:.o=.d)

EXETHDMIIMSSMBC_DEP := \
		$(EXETHDMIIMSSMBC_OBJ:.o=.d)

LLTHDMIIMSSMBC_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLTHDMIIMSSMBC_SRC)))

LLTHDMIIMSSMBC_OBJ  := $(LLTHDMIIMSSMBC_SRC:.cpp=.o)
LLTHDMIIMSSMBC_LIB  := $(LLTHDMIIMSSMBC_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBTHDMIIMSSMBC     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_THDMIIMSSMBC := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_THDMIIMSSMBC := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBTHDMIIMSSMBC) $(EXETHDMIIMSSMBC_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(THDMIIMSSMBC_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBTHDMIIMSSMBC_SRC) $(THDMIIMSSMBC_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBTHDMIIMSSMBC_HDR) $(THDMIIMSSMBC_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXETHDMIIMSSMBC_SRC) $(THDMIIMSSMBC_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLTHDMIIMSSMBC_SRC) $(THDMIIMSSMBC_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLTHDMIIMSSMBC_MMA) $(THDMIIMSSMBC_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(THDMIIMSSMBC_MK) $(THDMIIMSSMBC_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(THDMIIMSSMBC_INCLUDE_MK) $(THDMIIMSSMBC_INSTALL_DIR)
ifneq ($(THDMIIMSSMBC_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(THDMIIMSSMBC_SLHA_INPUT) $(THDMIIMSSMBC_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(THDMIIMSSMBC_REFERENCES) $(THDMIIMSSMBC_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(THDMIIMSSMBC_GNUPLOT) $(THDMIIMSSMBC_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBTHDMIIMSSMBC_DEP)
		-rm -f $(EXETHDMIIMSSMBC_DEP)
		-rm -f $(LLTHDMIIMSSMBC_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBTHDMIIMSSMBC)
		-rm -f $(LLTHDMIIMSSMBC_LIB)

clean-$(MODNAME)-obj:
		-rm -f $(LIBTHDMIIMSSMBC_OBJ)
		-rm -f $(EXETHDMIIMSSMBC_OBJ)
		-rm -f $(LLTHDMIIMSSMBC_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		-rm -f $(LIBTHDMIIMSSMBC_SRC)
		-rm -f $(LIBTHDMIIMSSMBC_HDR)
		-rm -f $(EXETHDMIIMSSMBC_SRC)
		-rm -f $(LLTHDMIIMSSMBC_SRC)
		-rm -f $(LLTHDMIIMSSMBC_MMA)
		-rm -f $(METACODE_STAMP_THDMIIMSSMBC)
		-rm -f $(THDMIIMSSMBC_INCLUDE_MK)
		-rm -f $(THDMIIMSSMBC_SLHA_INPUT)
		-rm -f $(THDMIIMSSMBC_REFERENCES)
		-rm -f $(THDMIIMSSMBC_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		-rm -f $(EXETHDMIIMSSMBC_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(THDMIIMSSMBC_TARBALL) \
		$(LIBTHDMIIMSSMBC_SRC) $(LIBTHDMIIMSSMBC_HDR) \
		$(EXETHDMIIMSSMBC_SRC) \
		$(LLTHDMIIMSSMBC_SRC) $(LLTHDMIIMSSMBC_MMA) \
		$(THDMIIMSSMBC_MK) $(THDMIIMSSMBC_INCLUDE_MK) \
		$(THDMIIMSSMBC_SLHA_INPUT) $(THDMIIMSSMBC_REFERENCES) \
		$(THDMIIMSSMBC_GNUPLOT)

$(LIBTHDMIIMSSMBC_SRC) $(LIBTHDMIIMSSMBC_HDR) $(EXETHDMIIMSSMBC_SRC) $(LLTHDMIIMSSMBC_SRC) $(LLTHDMIIMSSMBC_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_THDMIIMSSMBC)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_THDMIIMSSMBC): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_THDMIIMSSMBC)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_THDMIIMSSMBC)"
		@echo "Note: to regenerate THDMIIMSSMBC source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_THDMIIMSSMBC)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_THDMIIMSSMBC):
		@true
endif

$(LIBTHDMIIMSSMBC_DEP) $(EXETHDMIIMSSMBC_DEP) $(LLTHDMIIMSSMBC_DEP) $(LIBTHDMIIMSSMBC_OBJ) $(EXETHDMIIMSSMBC_OBJ) $(LLTHDMIIMSSMBC_OBJ) $(LLTHDMIIMSSMBC_LIB): \
	CPPFLAGS += $(MODTHDMIIMSSMBC_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBTHDMIIMSSMBC_DEP) $(EXETHDMIIMSSMBC_DEP) $(LLTHDMIIMSSMBC_DEP) $(LIBTHDMIIMSSMBC_OBJ) $(EXETHDMIIMSSMBC_OBJ) $(LLTHDMIIMSSMBC_OBJ) $(LLTHDMIIMSSMBC_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLTHDMIIMSSMBC_OBJ) $(LLTHDMIIMSSMBC_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBTHDMIIMSSMBC): $(LIBTHDMIIMSSMBC_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBTHDMIIMSSMBC) $(MODTHDMIIMSSMBC_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLTHDMIIMSSMBC_LIB): $(LLTHDMIIMSSMBC_OBJ) $(LIBTHDMIIMSSMBC) $(MODTHDMIIMSSMBC_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS)

ALLDEP += $(LIBTHDMIIMSSMBC_DEP) $(EXETHDMIIMSSMBC_DEP)
ALLSRC += $(LIBTHDMIIMSSMBC_SRC) $(EXETHDMIIMSSMBC_SRC)
ALLLIB += $(LIBTHDMIIMSSMBC)
ALLEXE += $(EXETHDMIIMSSMBC_EXE)
ALLMODDEP += $(MODTHDMIIMSSMBC_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLTHDMIIMSSMBC_DEP)
ALLSRC += $(LLTHDMIIMSSMBC_SRC)
ALLLL  += $(LLTHDMIIMSSMBC_LIB)
endif
