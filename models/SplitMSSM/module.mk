DIR          := models/SplitMSSM
MODNAME      := SplitMSSM
SARAH_MODEL  := SplitMSSM
WITH_$(MODNAME) := yes
MODSplitMSSM_MOD := SM SplitMSSM
MODSplitMSSM_DEP := $(patsubst %,model_specific/%,$(MODSplitMSSM_MOD))
MODSplitMSSM_INC := $(patsubst %,-Imodel_specific/%,$(MODSplitMSSM_MOD))
MODSplitMSSM_LIB := $(foreach M,$(MODSplitMSSM_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

SplitMSSM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

SplitMSSM_MK     := \
		$(DIR)/module.mk

SplitMSSM_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

SplitMSSM_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

SplitMSSM_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

SplitMSSM_INCLUDE_MK := \
		$(SplitMSSM_SUSY_BETAS_MK) \
		$(SplitMSSM_SOFT_BETAS_MK)

SplitMSSM_SLHA_INPUT := \
		$(DIR)/LesHouches.in.SplitMSSM_generated \
		$(DIR)/LesHouches.in.SplitMSSM

SplitMSSM_REFERENCES := \
		$(DIR)/SplitMSSM_references.tex

SplitMSSM_GNUPLOT := \
		$(DIR)/SplitMSSM_plot_rgflow.gnuplot \
		$(DIR)/SplitMSSM_plot_spectrum.gnuplot

SplitMSSM_TARBALL := \
		$(MODNAME).tar.gz

LIBSplitMSSM_SRC := \
		$(DIR)/SplitMSSM_a_muon.cpp \
		$(DIR)/SplitMSSM_edm.cpp \
		$(DIR)/SplitMSSM_effective_couplings.cpp \
		$(DIR)/SplitMSSM_info.cpp \
		$(DIR)/SplitMSSM_input_parameters.cpp \
		$(DIR)/SplitMSSM_mass_eigenstates.cpp \
		$(DIR)/SplitMSSM_observables.cpp \
		$(DIR)/SplitMSSM_physical.cpp \
		$(DIR)/SplitMSSM_slha_io.cpp \
		$(DIR)/SplitMSSM_soft_parameters.cpp \
		$(DIR)/SplitMSSM_susy_parameters.cpp \
		$(DIR)/SplitMSSM_utilities.cpp \
		$(DIR)/SplitMSSM_weinberg_angle.cpp

EXESplitMSSM_SRC := \
		$(DIR)/run_SplitMSSM.cpp \
		$(DIR)/run_cmd_line_SplitMSSM.cpp \
		$(DIR)/scan_SplitMSSM.cpp
LLSplitMSSM_LIB  :=
LLSplitMSSM_OBJ  :=
LLSplitMSSM_SRC  := \
		$(DIR)/SplitMSSM_librarylink.cpp

LLSplitMSSM_MMA  := \
		$(DIR)/SplitMSSM_librarylink.m \
		$(DIR)/run_SplitMSSM.m

LIBSplitMSSM_HDR := \
		$(DIR)/SplitMSSM_cxx_diagrams.hpp \
		$(DIR)/SplitMSSM_a_muon.hpp \
		$(DIR)/SplitMSSM_convergence_tester.hpp \
		$(DIR)/SplitMSSM_edm.hpp \
		$(DIR)/SplitMSSM_effective_couplings.hpp \
		$(DIR)/SplitMSSM_ewsb_solver.hpp \
		$(DIR)/SplitMSSM_ewsb_solver_interface.hpp \
		$(DIR)/SplitMSSM_high_scale_constraint.hpp \
		$(DIR)/SplitMSSM_info.hpp \
		$(DIR)/SplitMSSM_initial_guesser.hpp \
		$(DIR)/SplitMSSM_input_parameters.hpp \
		$(DIR)/SplitMSSM_low_scale_constraint.hpp \
		$(DIR)/SplitMSSM_mass_eigenstates.hpp \
		$(DIR)/SplitMSSM_model.hpp \
		$(DIR)/SplitMSSM_model_slha.hpp \
		$(DIR)/SplitMSSM_observables.hpp \
		$(DIR)/SplitMSSM_physical.hpp \
		$(DIR)/SplitMSSM_slha_io.hpp \
		$(DIR)/SplitMSSM_spectrum_generator.hpp \
		$(DIR)/SplitMSSM_spectrum_generator_interface.hpp \
		$(DIR)/SplitMSSM_soft_parameters.hpp \
		$(DIR)/SplitMSSM_susy_parameters.hpp \
		$(DIR)/SplitMSSM_susy_scale_constraint.hpp \
		$(DIR)/SplitMSSM_utilities.hpp \
		$(DIR)/SplitMSSM_weinberg_angle.hpp

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
-include $(SplitMSSM_SUSY_BETAS_MK)
-include $(SplitMSSM_SOFT_BETAS_MK)
-include $(SplitMSSM_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(SplitMSSM_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(SplitMSSM_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(SplitMSSM_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBSplitMSSM_SRC := $(sort $(LIBSplitMSSM_SRC))
EXESplitMSSM_SRC := $(sort $(EXESplitMSSM_SRC))

LIBSplitMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBSplitMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBSplitMSSM_SRC)))

EXESplitMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXESplitMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXESplitMSSM_SRC)))

EXESplitMSSM_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXESplitMSSM_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXESplitMSSM_SRC)))

LIBSplitMSSM_DEP := \
		$(LIBSplitMSSM_OBJ:.o=.d)

EXESplitMSSM_DEP := \
		$(EXESplitMSSM_OBJ:.o=.d)

LLSplitMSSM_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLSplitMSSM_SRC)))

LLSplitMSSM_OBJ  := $(LLSplitMSSM_SRC:.cpp=.o)
LLSplitMSSM_LIB  := $(LLSplitMSSM_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBSplitMSSM     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_SplitMSSM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_SplitMSSM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBSplitMSSM) $(EXESplitMSSM_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(SplitMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBSplitMSSM_SRC) $(SplitMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBSplitMSSM_HDR) $(SplitMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXESplitMSSM_SRC) $(SplitMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLSplitMSSM_SRC) $(SplitMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLSplitMSSM_MMA) $(SplitMSSM_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(SplitMSSM_MK) $(SplitMSSM_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(SplitMSSM_INCLUDE_MK) $(SplitMSSM_INSTALL_DIR)
ifneq ($(SplitMSSM_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(SplitMSSM_SLHA_INPUT) $(SplitMSSM_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(SplitMSSM_REFERENCES) $(SplitMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(SplitMSSM_GNUPLOT) $(SplitMSSM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBSplitMSSM_DEP)
		-rm -f $(EXESplitMSSM_DEP)
		-rm -f $(LLSplitMSSM_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBSplitMSSM)
		-rm -f $(LLSplitMSSM_LIB)

clean-$(MODNAME)-obj:
		-rm -f $(LIBSplitMSSM_OBJ)
		-rm -f $(EXESplitMSSM_OBJ)
		-rm -f $(LLSplitMSSM_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		-rm -f $(LIBSplitMSSM_SRC)
		-rm -f $(LIBSplitMSSM_HDR)
		-rm -f $(EXESplitMSSM_SRC)
		-rm -f $(LLSplitMSSM_SRC)
		-rm -f $(LLSplitMSSM_MMA)
		-rm -f $(METACODE_STAMP_SplitMSSM)
		-rm -f $(SplitMSSM_INCLUDE_MK)
		-rm -f $(SplitMSSM_SLHA_INPUT)
		-rm -f $(SplitMSSM_REFERENCES)
		-rm -f $(SplitMSSM_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		-rm -f $(EXESplitMSSM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(SplitMSSM_TARBALL) \
		$(LIBSplitMSSM_SRC) $(LIBSplitMSSM_HDR) \
		$(EXESplitMSSM_SRC) \
		$(LLSplitMSSM_SRC) $(LLSplitMSSM_MMA) \
		$(SplitMSSM_MK) $(SplitMSSM_INCLUDE_MK) \
		$(SplitMSSM_SLHA_INPUT) $(SplitMSSM_REFERENCES) \
		$(SplitMSSM_GNUPLOT)

$(LIBSplitMSSM_SRC) $(LIBSplitMSSM_HDR) $(EXESplitMSSM_SRC) $(LLSplitMSSM_SRC) $(LLSplitMSSM_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_SplitMSSM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_SplitMSSM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_SplitMSSM)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_SplitMSSM)"
		@echo "Note: to regenerate SplitMSSM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_SplitMSSM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_SplitMSSM):
		@true
endif

$(LIBSplitMSSM_DEP) $(EXESplitMSSM_DEP) $(LLSplitMSSM_DEP) $(LIBSplitMSSM_OBJ) $(EXESplitMSSM_OBJ) $(LLSplitMSSM_OBJ) $(LLSplitMSSM_LIB): \
	CPPFLAGS += $(MODSplitMSSM_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBSplitMSSM_DEP) $(EXESplitMSSM_DEP) $(LLSplitMSSM_DEP) $(LIBSplitMSSM_OBJ) $(EXESplitMSSM_OBJ) $(LLSplitMSSM_OBJ) $(LLSplitMSSM_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLSplitMSSM_OBJ) $(LLSplitMSSM_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBSplitMSSM): $(LIBSplitMSSM_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBSplitMSSM) $(MODSplitMSSM_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLSplitMSSM_LIB): $(LLSplitMSSM_OBJ) $(LIBSplitMSSM) $(MODSplitMSSM_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS)

ALLDEP += $(LIBSplitMSSM_DEP) $(EXESplitMSSM_DEP)
ALLSRC += $(LIBSplitMSSM_SRC) $(EXESplitMSSM_SRC)
ALLLIB += $(LIBSplitMSSM)
ALLEXE += $(EXESplitMSSM_EXE)
ALLMODDEP += $(MODSplitMSSM_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLSplitMSSM_DEP)
ALLSRC += $(LLSplitMSSM_SRC)
ALLLL  += $(LLSplitMSSM_LIB)
endif
