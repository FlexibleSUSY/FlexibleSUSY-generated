DIR          := models/E6SSMEFTHiggs
MODNAME      := E6SSMEFTHiggs
SARAH_MODEL  := E6SSM
WITH_$(MODNAME) := yes
MODE6SSMEFTHiggs_MOD := SM MSSM_higgs NMSSM_higgs
MODE6SSMEFTHiggs_DEP := $(patsubst %,model_specific/%,$(MODE6SSMEFTHiggs_MOD))
MODE6SSMEFTHiggs_INC := $(patsubst %,-Imodel_specific/%,$(MODE6SSMEFTHiggs_MOD))
MODE6SSMEFTHiggs_LIB := $(foreach M,$(MODE6SSMEFTHiggs_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

E6SSMEFTHiggs_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

E6SSMEFTHiggs_MK     := \
		$(DIR)/module.mk

E6SSMEFTHiggs_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

E6SSMEFTHiggs_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

E6SSMEFTHiggs_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

E6SSMEFTHiggs_INCLUDE_MK := \
		$(E6SSMEFTHiggs_SUSY_BETAS_MK) \
		$(E6SSMEFTHiggs_SOFT_BETAS_MK)

E6SSMEFTHiggs_SLHA_INPUT := \
		$(DIR)/LesHouches.in.E6SSMEFTHiggs_generated \
		$(DIR)/LesHouches.in.E6SSMEFTHiggs

E6SSMEFTHiggs_REFERENCES := \
		$(DIR)/E6SSMEFTHiggs_references.tex

E6SSMEFTHiggs_GNUPLOT := \
		$(DIR)/E6SSMEFTHiggs_plot_rgflow.gnuplot \
		$(DIR)/E6SSMEFTHiggs_plot_spectrum.gnuplot

E6SSMEFTHiggs_TARBALL := \
		$(MODNAME).tar.gz

LIBE6SSMEFTHiggs_SRC := \
		$(DIR)/E6SSMEFTHiggs_a_muon.cpp \
		$(DIR)/E6SSMEFTHiggs_edm.cpp \
		$(DIR)/E6SSMEFTHiggs_effective_couplings.cpp \
		$(DIR)/E6SSMEFTHiggs_info.cpp \
		$(DIR)/E6SSMEFTHiggs_input_parameters.cpp \
		$(DIR)/E6SSMEFTHiggs_mass_eigenstates.cpp \
		$(DIR)/E6SSMEFTHiggs_observables.cpp \
		$(DIR)/E6SSMEFTHiggs_physical.cpp \
		$(DIR)/E6SSMEFTHiggs_slha_io.cpp \
		$(DIR)/E6SSMEFTHiggs_soft_parameters.cpp \
		$(DIR)/E6SSMEFTHiggs_susy_parameters.cpp \
		$(DIR)/E6SSMEFTHiggs_utilities.cpp \
		$(DIR)/E6SSMEFTHiggs_weinberg_angle.cpp

EXEE6SSMEFTHiggs_SRC := \
		$(DIR)/run_E6SSMEFTHiggs.cpp \
		$(DIR)/run_cmd_line_E6SSMEFTHiggs.cpp \
		$(DIR)/scan_E6SSMEFTHiggs.cpp
LLE6SSMEFTHiggs_LIB  :=
LLE6SSMEFTHiggs_OBJ  :=
LLE6SSMEFTHiggs_SRC  := \
		$(DIR)/E6SSMEFTHiggs_librarylink.cpp

LLE6SSMEFTHiggs_MMA  := \
		$(DIR)/E6SSMEFTHiggs_librarylink.m \
		$(DIR)/run_E6SSMEFTHiggs.m

LIBE6SSMEFTHiggs_HDR := \
		$(DIR)/E6SSMEFTHiggs_cxx_diagrams.hpp \
		$(DIR)/E6SSMEFTHiggs_a_muon.hpp \
		$(DIR)/E6SSMEFTHiggs_convergence_tester.hpp \
		$(DIR)/E6SSMEFTHiggs_edm.hpp \
		$(DIR)/E6SSMEFTHiggs_effective_couplings.hpp \
		$(DIR)/E6SSMEFTHiggs_ewsb_solver.hpp \
		$(DIR)/E6SSMEFTHiggs_ewsb_solver_interface.hpp \
		$(DIR)/E6SSMEFTHiggs_high_scale_constraint.hpp \
		$(DIR)/E6SSMEFTHiggs_info.hpp \
		$(DIR)/E6SSMEFTHiggs_initial_guesser.hpp \
		$(DIR)/E6SSMEFTHiggs_input_parameters.hpp \
		$(DIR)/E6SSMEFTHiggs_low_scale_constraint.hpp \
		$(DIR)/E6SSMEFTHiggs_mass_eigenstates.hpp \
		$(DIR)/E6SSMEFTHiggs_model.hpp \
		$(DIR)/E6SSMEFTHiggs_model_slha.hpp \
		$(DIR)/E6SSMEFTHiggs_observables.hpp \
		$(DIR)/E6SSMEFTHiggs_physical.hpp \
		$(DIR)/E6SSMEFTHiggs_slha_io.hpp \
		$(DIR)/E6SSMEFTHiggs_spectrum_generator.hpp \
		$(DIR)/E6SSMEFTHiggs_spectrum_generator_interface.hpp \
		$(DIR)/E6SSMEFTHiggs_soft_parameters.hpp \
		$(DIR)/E6SSMEFTHiggs_susy_parameters.hpp \
		$(DIR)/E6SSMEFTHiggs_susy_scale_constraint.hpp \
		$(DIR)/E6SSMEFTHiggs_utilities.hpp \
		$(DIR)/E6SSMEFTHiggs_weinberg_angle.hpp

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
-include $(E6SSMEFTHiggs_SUSY_BETAS_MK)
-include $(E6SSMEFTHiggs_SOFT_BETAS_MK)
-include $(E6SSMEFTHiggs_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(E6SSMEFTHiggs_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(E6SSMEFTHiggs_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(E6SSMEFTHiggs_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBE6SSMEFTHiggs_SRC := $(sort $(LIBE6SSMEFTHiggs_SRC))
EXEE6SSMEFTHiggs_SRC := $(sort $(EXEE6SSMEFTHiggs_SRC))

LIBE6SSMEFTHiggs_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBE6SSMEFTHiggs_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBE6SSMEFTHiggs_SRC)))

EXEE6SSMEFTHiggs_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEE6SSMEFTHiggs_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEE6SSMEFTHiggs_SRC)))

EXEE6SSMEFTHiggs_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXEE6SSMEFTHiggs_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXEE6SSMEFTHiggs_SRC)))

LIBE6SSMEFTHiggs_DEP := \
		$(LIBE6SSMEFTHiggs_OBJ:.o=.d)

EXEE6SSMEFTHiggs_DEP := \
		$(EXEE6SSMEFTHiggs_OBJ:.o=.d)

LLE6SSMEFTHiggs_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLE6SSMEFTHiggs_SRC)))

LLE6SSMEFTHiggs_OBJ  := $(LLE6SSMEFTHiggs_SRC:.cpp=.o)
LLE6SSMEFTHiggs_LIB  := $(LLE6SSMEFTHiggs_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBE6SSMEFTHiggs     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_E6SSMEFTHiggs := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_E6SSMEFTHiggs := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBE6SSMEFTHiggs) $(EXEE6SSMEFTHiggs_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(E6SSMEFTHiggs_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBE6SSMEFTHiggs_SRC) $(E6SSMEFTHiggs_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBE6SSMEFTHiggs_HDR) $(E6SSMEFTHiggs_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXEE6SSMEFTHiggs_SRC) $(E6SSMEFTHiggs_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLE6SSMEFTHiggs_SRC) $(E6SSMEFTHiggs_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLE6SSMEFTHiggs_MMA) $(E6SSMEFTHiggs_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(E6SSMEFTHiggs_MK) $(E6SSMEFTHiggs_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(E6SSMEFTHiggs_INCLUDE_MK) $(E6SSMEFTHiggs_INSTALL_DIR)
ifneq ($(E6SSMEFTHiggs_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(E6SSMEFTHiggs_SLHA_INPUT) $(E6SSMEFTHiggs_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(E6SSMEFTHiggs_REFERENCES) $(E6SSMEFTHiggs_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(E6SSMEFTHiggs_GNUPLOT) $(E6SSMEFTHiggs_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBE6SSMEFTHiggs_DEP)
		-rm -f $(EXEE6SSMEFTHiggs_DEP)
		-rm -f $(LLE6SSMEFTHiggs_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBE6SSMEFTHiggs)
		-rm -f $(LLE6SSMEFTHiggs_LIB)

clean-$(MODNAME)-obj:
		-rm -f $(LIBE6SSMEFTHiggs_OBJ)
		-rm -f $(EXEE6SSMEFTHiggs_OBJ)
		-rm -f $(LLE6SSMEFTHiggs_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		-rm -f $(LIBE6SSMEFTHiggs_SRC)
		-rm -f $(LIBE6SSMEFTHiggs_HDR)
		-rm -f $(EXEE6SSMEFTHiggs_SRC)
		-rm -f $(LLE6SSMEFTHiggs_SRC)
		-rm -f $(LLE6SSMEFTHiggs_MMA)
		-rm -f $(METACODE_STAMP_E6SSMEFTHiggs)
		-rm -f $(E6SSMEFTHiggs_INCLUDE_MK)
		-rm -f $(E6SSMEFTHiggs_SLHA_INPUT)
		-rm -f $(E6SSMEFTHiggs_REFERENCES)
		-rm -f $(E6SSMEFTHiggs_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		-rm -f $(EXEE6SSMEFTHiggs_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(E6SSMEFTHiggs_TARBALL) \
		$(LIBE6SSMEFTHiggs_SRC) $(LIBE6SSMEFTHiggs_HDR) \
		$(EXEE6SSMEFTHiggs_SRC) \
		$(LLE6SSMEFTHiggs_SRC) $(LLE6SSMEFTHiggs_MMA) \
		$(E6SSMEFTHiggs_MK) $(E6SSMEFTHiggs_INCLUDE_MK) \
		$(E6SSMEFTHiggs_SLHA_INPUT) $(E6SSMEFTHiggs_REFERENCES) \
		$(E6SSMEFTHiggs_GNUPLOT)

$(LIBE6SSMEFTHiggs_SRC) $(LIBE6SSMEFTHiggs_HDR) $(EXEE6SSMEFTHiggs_SRC) $(LLE6SSMEFTHiggs_SRC) $(LLE6SSMEFTHiggs_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_E6SSMEFTHiggs)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_E6SSMEFTHiggs): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_E6SSMEFTHiggs)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_E6SSMEFTHiggs)"
		@echo "Note: to regenerate E6SSMEFTHiggs source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_E6SSMEFTHiggs)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_E6SSMEFTHiggs):
		@true
endif

$(LIBE6SSMEFTHiggs_DEP) $(EXEE6SSMEFTHiggs_DEP) $(LLE6SSMEFTHiggs_DEP) $(LIBE6SSMEFTHiggs_OBJ) $(EXEE6SSMEFTHiggs_OBJ) $(LLE6SSMEFTHiggs_OBJ) $(LLE6SSMEFTHiggs_LIB): \
	CPPFLAGS += $(MODE6SSMEFTHiggs_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBE6SSMEFTHiggs_DEP) $(EXEE6SSMEFTHiggs_DEP) $(LLE6SSMEFTHiggs_DEP) $(LIBE6SSMEFTHiggs_OBJ) $(EXEE6SSMEFTHiggs_OBJ) $(LLE6SSMEFTHiggs_OBJ) $(LLE6SSMEFTHiggs_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLE6SSMEFTHiggs_OBJ) $(LLE6SSMEFTHiggs_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBE6SSMEFTHiggs): $(LIBE6SSMEFTHiggs_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBE6SSMEFTHiggs) $(MODE6SSMEFTHiggs_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLE6SSMEFTHiggs_LIB): $(LLE6SSMEFTHiggs_OBJ) $(LIBE6SSMEFTHiggs) $(MODE6SSMEFTHiggs_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS)

ALLDEP += $(LIBE6SSMEFTHiggs_DEP) $(EXEE6SSMEFTHiggs_DEP)
ALLSRC += $(LIBE6SSMEFTHiggs_SRC) $(EXEE6SSMEFTHiggs_SRC)
ALLLIB += $(LIBE6SSMEFTHiggs)
ALLEXE += $(EXEE6SSMEFTHiggs_EXE)
ALLMODDEP += $(MODE6SSMEFTHiggs_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLE6SSMEFTHiggs_DEP)
ALLSRC += $(LLE6SSMEFTHiggs_SRC)
ALLLL  += $(LLE6SSMEFTHiggs_LIB)
endif
