DIR          := models/HGTHDMIIMSSMBC
MODNAME      := HGTHDMIIMSSMBC
SARAH_MODEL  := HGTHDM-II
WITH_$(MODNAME) := yes
MODHGTHDMIIMSSMBC_MOD := SM
MODHGTHDMIIMSSMBC_DEP := $(patsubst %,model_specific/%,$(MODHGTHDMIIMSSMBC_MOD))
MODHGTHDMIIMSSMBC_INC := $(patsubst %,-Imodel_specific/%,$(MODHGTHDMIIMSSMBC_MOD))
MODHGTHDMIIMSSMBC_LIB := $(foreach M,$(MODHGTHDMIIMSSMBC_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

HGTHDMIIMSSMBC_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

HGTHDMIIMSSMBC_MK     := \
		$(DIR)/module.mk

HGTHDMIIMSSMBC_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

HGTHDMIIMSSMBC_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

HGTHDMIIMSSMBC_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

HGTHDMIIMSSMBC_INCLUDE_MK := \
		$(HGTHDMIIMSSMBC_SUSY_BETAS_MK) \
		$(HGTHDMIIMSSMBC_SOFT_BETAS_MK)

HGTHDMIIMSSMBC_SLHA_INPUT := \
		$(DIR)/LesHouches.in.HGTHDMIIMSSMBC_generated \
		$(DIR)/LesHouches.in.HGTHDMIIMSSMBC

HGTHDMIIMSSMBC_REFERENCES := \
		$(DIR)/HGTHDMIIMSSMBC_references.tex

HGTHDMIIMSSMBC_GNUPLOT := \
		$(DIR)/HGTHDMIIMSSMBC_plot_rgflow.gnuplot \
		$(DIR)/HGTHDMIIMSSMBC_plot_spectrum.gnuplot

HGTHDMIIMSSMBC_TARBALL := \
		$(MODNAME).tar.gz

LIBHGTHDMIIMSSMBC_SRC := \
		$(DIR)/HGTHDMIIMSSMBC_a_muon.cpp \
		$(DIR)/HGTHDMIIMSSMBC_edm.cpp \
		$(DIR)/HGTHDMIIMSSMBC_effective_couplings.cpp \
		$(DIR)/HGTHDMIIMSSMBC_info.cpp \
		$(DIR)/HGTHDMIIMSSMBC_input_parameters.cpp \
		$(DIR)/HGTHDMIIMSSMBC_mass_eigenstates.cpp \
		$(DIR)/HGTHDMIIMSSMBC_observables.cpp \
		$(DIR)/HGTHDMIIMSSMBC_physical.cpp \
		$(DIR)/HGTHDMIIMSSMBC_slha_io.cpp \
		$(DIR)/HGTHDMIIMSSMBC_soft_parameters.cpp \
		$(DIR)/HGTHDMIIMSSMBC_susy_parameters.cpp \
		$(DIR)/HGTHDMIIMSSMBC_utilities.cpp \
		$(DIR)/HGTHDMIIMSSMBC_weinberg_angle.cpp

EXEHGTHDMIIMSSMBC_SRC := \
		$(DIR)/run_HGTHDMIIMSSMBC.cpp \
		$(DIR)/run_cmd_line_HGTHDMIIMSSMBC.cpp \
		$(DIR)/scan_HGTHDMIIMSSMBC.cpp
LLHGTHDMIIMSSMBC_LIB  :=
LLHGTHDMIIMSSMBC_OBJ  :=
LLHGTHDMIIMSSMBC_SRC  := \
		$(DIR)/HGTHDMIIMSSMBC_librarylink.cpp

LLHGTHDMIIMSSMBC_MMA  := \
		$(DIR)/HGTHDMIIMSSMBC_librarylink.m \
		$(DIR)/run_HGTHDMIIMSSMBC.m

LIBHGTHDMIIMSSMBC_HDR := \
		$(DIR)/HGTHDMIIMSSMBC_cxx_diagrams.hpp \
		$(DIR)/HGTHDMIIMSSMBC_a_muon.hpp \
		$(DIR)/HGTHDMIIMSSMBC_convergence_tester.hpp \
		$(DIR)/HGTHDMIIMSSMBC_edm.hpp \
		$(DIR)/HGTHDMIIMSSMBC_effective_couplings.hpp \
		$(DIR)/HGTHDMIIMSSMBC_ewsb_solver.hpp \
		$(DIR)/HGTHDMIIMSSMBC_ewsb_solver_interface.hpp \
		$(DIR)/HGTHDMIIMSSMBC_high_scale_constraint.hpp \
		$(DIR)/HGTHDMIIMSSMBC_info.hpp \
		$(DIR)/HGTHDMIIMSSMBC_initial_guesser.hpp \
		$(DIR)/HGTHDMIIMSSMBC_input_parameters.hpp \
		$(DIR)/HGTHDMIIMSSMBC_low_scale_constraint.hpp \
		$(DIR)/HGTHDMIIMSSMBC_mass_eigenstates.hpp \
		$(DIR)/HGTHDMIIMSSMBC_model.hpp \
		$(DIR)/HGTHDMIIMSSMBC_model_slha.hpp \
		$(DIR)/HGTHDMIIMSSMBC_observables.hpp \
		$(DIR)/HGTHDMIIMSSMBC_physical.hpp \
		$(DIR)/HGTHDMIIMSSMBC_slha_io.hpp \
		$(DIR)/HGTHDMIIMSSMBC_spectrum_generator.hpp \
		$(DIR)/HGTHDMIIMSSMBC_spectrum_generator_interface.hpp \
		$(DIR)/HGTHDMIIMSSMBC_soft_parameters.hpp \
		$(DIR)/HGTHDMIIMSSMBC_susy_parameters.hpp \
		$(DIR)/HGTHDMIIMSSMBC_susy_scale_constraint.hpp \
		$(DIR)/HGTHDMIIMSSMBC_utilities.hpp \
		$(DIR)/HGTHDMIIMSSMBC_weinberg_angle.hpp

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
-include $(HGTHDMIIMSSMBC_SUSY_BETAS_MK)
-include $(HGTHDMIIMSSMBC_SOFT_BETAS_MK)
-include $(HGTHDMIIMSSMBC_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(HGTHDMIIMSSMBC_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(HGTHDMIIMSSMBC_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(HGTHDMIIMSSMBC_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBHGTHDMIIMSSMBC_SRC := $(sort $(LIBHGTHDMIIMSSMBC_SRC))
EXEHGTHDMIIMSSMBC_SRC := $(sort $(EXEHGTHDMIIMSSMBC_SRC))

LIBHGTHDMIIMSSMBC_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBHGTHDMIIMSSMBC_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBHGTHDMIIMSSMBC_SRC)))

EXEHGTHDMIIMSSMBC_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEHGTHDMIIMSSMBC_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEHGTHDMIIMSSMBC_SRC)))

EXEHGTHDMIIMSSMBC_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXEHGTHDMIIMSSMBC_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXEHGTHDMIIMSSMBC_SRC)))

LIBHGTHDMIIMSSMBC_DEP := \
		$(LIBHGTHDMIIMSSMBC_OBJ:.o=.d)

EXEHGTHDMIIMSSMBC_DEP := \
		$(EXEHGTHDMIIMSSMBC_OBJ:.o=.d)

LLHGTHDMIIMSSMBC_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLHGTHDMIIMSSMBC_SRC)))

LLHGTHDMIIMSSMBC_OBJ  := $(LLHGTHDMIIMSSMBC_SRC:.cpp=.o)
LLHGTHDMIIMSSMBC_LIB  := $(LLHGTHDMIIMSSMBC_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBHGTHDMIIMSSMBC     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_HGTHDMIIMSSMBC := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_HGTHDMIIMSSMBC := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL_DIR) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBHGTHDMIIMSSMBC) $(EXEHGTHDMIIMSSMBC_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(HGTHDMIIMSSMBC_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBHGTHDMIIMSSMBC_SRC) $(HGTHDMIIMSSMBC_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBHGTHDMIIMSSMBC_HDR) $(HGTHDMIIMSSMBC_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXEHGTHDMIIMSSMBC_SRC) $(HGTHDMIIMSSMBC_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLHGTHDMIIMSSMBC_SRC) $(HGTHDMIIMSSMBC_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLHGTHDMIIMSSMBC_MMA) $(HGTHDMIIMSSMBC_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(HGTHDMIIMSSMBC_MK) $(HGTHDMIIMSSMBC_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(HGTHDMIIMSSMBC_INCLUDE_MK) $(HGTHDMIIMSSMBC_INSTALL_DIR)
ifneq ($(HGTHDMIIMSSMBC_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(HGTHDMIIMSSMBC_SLHA_INPUT) $(HGTHDMIIMSSMBC_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(HGTHDMIIMSSMBC_REFERENCES) $(HGTHDMIIMSSMBC_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(HGTHDMIIMSSMBC_GNUPLOT) $(HGTHDMIIMSSMBC_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBHGTHDMIIMSSMBC_DEP)
		-rm -f $(EXEHGTHDMIIMSSMBC_DEP)
		-rm -f $(LLHGTHDMIIMSSMBC_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBHGTHDMIIMSSMBC)
		-rm -f $(LLHGTHDMIIMSSMBC_LIB)

clean-$(MODNAME)-obj:
		-rm -f $(LIBHGTHDMIIMSSMBC_OBJ)
		-rm -f $(EXEHGTHDMIIMSSMBC_OBJ)
		-rm -f $(LLHGTHDMIIMSSMBC_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		-rm -f $(LIBHGTHDMIIMSSMBC_SRC)
		-rm -f $(LIBHGTHDMIIMSSMBC_HDR)
		-rm -f $(EXEHGTHDMIIMSSMBC_SRC)
		-rm -f $(LLHGTHDMIIMSSMBC_SRC)
		-rm -f $(LLHGTHDMIIMSSMBC_MMA)
		-rm -f $(METACODE_STAMP_HGTHDMIIMSSMBC)
		-rm -f $(HGTHDMIIMSSMBC_INCLUDE_MK)
		-rm -f $(HGTHDMIIMSSMBC_SLHA_INPUT)
		-rm -f $(HGTHDMIIMSSMBC_REFERENCES)
		-rm -f $(HGTHDMIIMSSMBC_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		-rm -f $(EXEHGTHDMIIMSSMBC_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(HGTHDMIIMSSMBC_TARBALL) \
		$(LIBHGTHDMIIMSSMBC_SRC) $(LIBHGTHDMIIMSSMBC_HDR) \
		$(EXEHGTHDMIIMSSMBC_SRC) \
		$(LLHGTHDMIIMSSMBC_SRC) $(LLHGTHDMIIMSSMBC_MMA) \
		$(HGTHDMIIMSSMBC_MK) $(HGTHDMIIMSSMBC_INCLUDE_MK) \
		$(HGTHDMIIMSSMBC_SLHA_INPUT) $(HGTHDMIIMSSMBC_REFERENCES) \
		$(HGTHDMIIMSSMBC_GNUPLOT)

$(LIBHGTHDMIIMSSMBC_SRC) $(LIBHGTHDMIIMSSMBC_HDR) $(EXEHGTHDMIIMSSMBC_SRC) $(LLHGTHDMIIMSSMBC_SRC) $(LLHGTHDMIIMSSMBC_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_HGTHDMIIMSSMBC)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_HGTHDMIIMSSMBC): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_HGTHDMIIMSSMBC)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]" || (echo "Error: The code generation failed!"; exit 1)
		@touch "$(METACODE_STAMP_HGTHDMIIMSSMBC)"
		@echo "Note: to regenerate HGTHDMIIMSSMBC source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_HGTHDMIIMSSMBC)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_HGTHDMIIMSSMBC):
		@true
endif

$(LIBHGTHDMIIMSSMBC_DEP) $(EXEHGTHDMIIMSSMBC_DEP) $(LLHGTHDMIIMSSMBC_DEP) $(LIBHGTHDMIIMSSMBC_OBJ) $(EXEHGTHDMIIMSSMBC_OBJ) $(LLHGTHDMIIMSSMBC_OBJ) $(LLHGTHDMIIMSSMBC_LIB): \
	CPPFLAGS += $(MODHGTHDMIIMSSMBC_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBHGTHDMIIMSSMBC_DEP) $(EXEHGTHDMIIMSSMBC_DEP) $(LLHGTHDMIIMSSMBC_DEP) $(LIBHGTHDMIIMSSMBC_OBJ) $(EXEHGTHDMIIMSSMBC_OBJ) $(LLHGTHDMIIMSSMBC_OBJ) $(LLHGTHDMIIMSSMBC_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLHGTHDMIIMSSMBC_OBJ) $(LLHGTHDMIIMSSMBC_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBHGTHDMIIMSSMBC): $(LIBHGTHDMIIMSSMBC_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBHGTHDMIIMSSMBC) $(MODHGTHDMIIMSSMBC_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLHGTHDMIIMSSMBC_LIB): $(LLHGTHDMIIMSSMBC_OBJ) $(LIBHGTHDMIIMSSMBC) $(MODHGTHDMIIMSSMBC_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS)

ALLDEP += $(LIBHGTHDMIIMSSMBC_DEP) $(EXEHGTHDMIIMSSMBC_DEP)
ALLSRC += $(LIBHGTHDMIIMSSMBC_SRC) $(EXEHGTHDMIIMSSMBC_SRC)
ALLLIB += $(LIBHGTHDMIIMSSMBC)
ALLEXE += $(EXEHGTHDMIIMSSMBC_EXE)
ALLMODDEP += $(MODHGTHDMIIMSSMBC_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLHGTHDMIIMSSMBC_DEP)
ALLSRC += $(LLHGTHDMIIMSSMBC_SRC)
ALLLL  += $(LLHGTHDMIIMSSMBC_LIB)
endif
