DIR          := models/MRSSMEFTHiggs
MODNAME      := MRSSMEFTHiggs
SARAH_MODEL  := MRSSM
WITH_$(MODNAME) := yes

MRSSMEFTHiggs_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

MRSSMEFTHiggs_MK     := \
		$(DIR)/module.mk

MRSSMEFTHiggs_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

MRSSMEFTHiggs_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

MRSSMEFTHiggs_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

MRSSMEFTHiggs_INCLUDE_MK := \
		$(MRSSMEFTHiggs_SUSY_BETAS_MK) \
		$(MRSSMEFTHiggs_SOFT_BETAS_MK)

MRSSMEFTHiggs_SLHA_INPUT := \
		$(DIR)/LesHouches.in.MRSSMEFTHiggs_generated \
		$(DIR)/LesHouches.in.MRSSMEFTHiggs

MRSSMEFTHiggs_GNUPLOT := \
		$(DIR)/MRSSMEFTHiggs_plot_rgflow.gnuplot \
		$(DIR)/MRSSMEFTHiggs_plot_spectrum.gnuplot

MRSSMEFTHiggs_TARBALL := \
		$(MODNAME).tar.gz

LIBMRSSMEFTHiggs_SRC := \
		$(DIR)/MRSSMEFTHiggs_a_muon.cpp \
		$(DIR)/MRSSMEFTHiggs_edm.cpp \
		$(DIR)/MRSSMEFTHiggs_effective_couplings.cpp \
		$(DIR)/MRSSMEFTHiggs_info.cpp \
		$(DIR)/MRSSMEFTHiggs_input_parameters.cpp \
		$(DIR)/MRSSMEFTHiggs_mass_eigenstates.cpp \
		$(DIR)/MRSSMEFTHiggs_observables.cpp \
		$(DIR)/MRSSMEFTHiggs_physical.cpp \
		$(DIR)/MRSSMEFTHiggs_slha_io.cpp \
		$(DIR)/MRSSMEFTHiggs_soft_parameters.cpp \
		$(DIR)/MRSSMEFTHiggs_susy_parameters.cpp \
		$(DIR)/MRSSMEFTHiggs_utilities.cpp \
		$(DIR)/MRSSMEFTHiggs_weinberg_angle.cpp

EXEMRSSMEFTHiggs_SRC := \
		$(DIR)/run_MRSSMEFTHiggs.cpp \
		$(DIR)/run_cmd_line_MRSSMEFTHiggs.cpp \
		$(DIR)/scan_MRSSMEFTHiggs.cpp
LLMRSSMEFTHiggs_LIB  :=
LLMRSSMEFTHiggs_OBJ  :=
LLMRSSMEFTHiggs_SRC  := \
		$(DIR)/MRSSMEFTHiggs_librarylink.cpp

LLMRSSMEFTHiggs_MMA  := \
		$(DIR)/MRSSMEFTHiggs_librarylink.m \
		$(DIR)/run_MRSSMEFTHiggs.m

LIBMRSSMEFTHiggs_HDR := \
		$(DIR)/MRSSMEFTHiggs_cxx_diagrams.hpp \
		$(DIR)/MRSSMEFTHiggs_a_muon.hpp \
		$(DIR)/MRSSMEFTHiggs_convergence_tester.hpp \
		$(DIR)/MRSSMEFTHiggs_edm.hpp \
		$(DIR)/MRSSMEFTHiggs_effective_couplings.hpp \
		$(DIR)/MRSSMEFTHiggs_ewsb_solver.hpp \
		$(DIR)/MRSSMEFTHiggs_ewsb_solver_interface.hpp \
		$(DIR)/MRSSMEFTHiggs_high_scale_constraint.hpp \
		$(DIR)/MRSSMEFTHiggs_info.hpp \
		$(DIR)/MRSSMEFTHiggs_initial_guesser.hpp \
		$(DIR)/MRSSMEFTHiggs_input_parameters.hpp \
		$(DIR)/MRSSMEFTHiggs_low_scale_constraint.hpp \
		$(DIR)/MRSSMEFTHiggs_mass_eigenstates.hpp \
		$(DIR)/MRSSMEFTHiggs_model.hpp \
		$(DIR)/MRSSMEFTHiggs_model_slha.hpp \
		$(DIR)/MRSSMEFTHiggs_observables.hpp \
		$(DIR)/MRSSMEFTHiggs_physical.hpp \
		$(DIR)/MRSSMEFTHiggs_slha_io.hpp \
		$(DIR)/MRSSMEFTHiggs_spectrum_generator.hpp \
		$(DIR)/MRSSMEFTHiggs_spectrum_generator_interface.hpp \
		$(DIR)/MRSSMEFTHiggs_soft_parameters.hpp \
		$(DIR)/MRSSMEFTHiggs_susy_parameters.hpp \
		$(DIR)/MRSSMEFTHiggs_susy_scale_constraint.hpp \
		$(DIR)/MRSSMEFTHiggs_utilities.hpp \
		$(DIR)/MRSSMEFTHiggs_weinberg_angle.hpp

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
-include $(MRSSMEFTHiggs_SUSY_BETAS_MK)
-include $(MRSSMEFTHiggs_SOFT_BETAS_MK)
-include $(MRSSMEFTHiggs_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(MRSSMEFTHiggs_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(MRSSMEFTHiggs_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(MRSSMEFTHiggs_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBMRSSMEFTHiggs_SRC := $(sort $(LIBMRSSMEFTHiggs_SRC))
EXEMRSSMEFTHiggs_SRC := $(sort $(EXEMRSSMEFTHiggs_SRC))

LIBMRSSMEFTHiggs_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBMRSSMEFTHiggs_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBMRSSMEFTHiggs_SRC)))

EXEMRSSMEFTHiggs_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEMRSSMEFTHiggs_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEMRSSMEFTHiggs_SRC)))

EXEMRSSMEFTHiggs_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXEMRSSMEFTHiggs_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXEMRSSMEFTHiggs_SRC)))

LIBMRSSMEFTHiggs_DEP := \
		$(LIBMRSSMEFTHiggs_OBJ:.o=.d)

EXEMRSSMEFTHiggs_DEP := \
		$(EXEMRSSMEFTHiggs_OBJ:.o=.d)

LLMRSSMEFTHiggs_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLMRSSMEFTHiggs_SRC)))

LLMRSSMEFTHiggs_OBJ  := $(LLMRSSMEFTHiggs_SRC:.cpp=.o)
LLMRSSMEFTHiggs_LIB  := $(LLMRSSMEFTHiggs_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBMRSSMEFTHiggs     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_MRSSMEFTHiggs := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_MRSSMEFTHiggs := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBMRSSMEFTHiggs) $(EXEMRSSMEFTHiggs_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(MRSSMEFTHiggs_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBMRSSMEFTHiggs_SRC) $(MRSSMEFTHiggs_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBMRSSMEFTHiggs_HDR) $(MRSSMEFTHiggs_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXEMRSSMEFTHiggs_SRC) $(MRSSMEFTHiggs_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLMRSSMEFTHiggs_SRC) $(MRSSMEFTHiggs_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLMRSSMEFTHiggs_MMA) $(MRSSMEFTHiggs_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(MRSSMEFTHiggs_MK) $(MRSSMEFTHiggs_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(MRSSMEFTHiggs_INCLUDE_MK) $(MRSSMEFTHiggs_INSTALL_DIR)
ifneq ($(MRSSMEFTHiggs_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(MRSSMEFTHiggs_SLHA_INPUT) $(MRSSMEFTHiggs_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(MRSSMEFTHiggs_GNUPLOT) $(MRSSMEFTHiggs_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBMRSSMEFTHiggs_DEP)
		-rm -f $(EXEMRSSMEFTHiggs_DEP)
		-rm -f $(LLMRSSMEFTHiggs_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBMRSSMEFTHiggs)
		-rm -f $(LLMRSSMEFTHiggs_LIB)

clean-$(MODNAME)-obj:
		-rm -f $(LIBMRSSMEFTHiggs_OBJ)
		-rm -f $(EXEMRSSMEFTHiggs_OBJ)
		-rm -f $(LLMRSSMEFTHiggs_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		-rm -f $(LIBMRSSMEFTHiggs_SRC)
		-rm -f $(LIBMRSSMEFTHiggs_HDR)
		-rm -f $(EXEMRSSMEFTHiggs_SRC)
		-rm -f $(LLMRSSMEFTHiggs_SRC)
		-rm -f $(LLMRSSMEFTHiggs_MMA)
		-rm -f $(METACODE_STAMP_MRSSMEFTHiggs)
		-rm -f $(MRSSMEFTHiggs_INCLUDE_MK)
		-rm -f $(MRSSMEFTHiggs_SLHA_INPUT)
		-rm -f $(MRSSMEFTHiggs_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		-rm -f $(EXEMRSSMEFTHiggs_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(MRSSMEFTHiggs_TARBALL) \
		$(LIBMRSSMEFTHiggs_SRC) $(LIBMRSSMEFTHiggs_HDR) \
		$(EXEMRSSMEFTHiggs_SRC) \
		$(LLMRSSMEFTHiggs_SRC) $(LLMRSSMEFTHiggs_MMA) \
		$(MRSSMEFTHiggs_MK) $(MRSSMEFTHiggs_INCLUDE_MK) \
		$(MRSSMEFTHiggs_SLHA_INPUT) $(MRSSMEFTHiggs_GNUPLOT)

$(LIBMRSSMEFTHiggs_SRC) $(LIBMRSSMEFTHiggs_HDR) $(EXEMRSSMEFTHiggs_SRC) $(LLMRSSMEFTHiggs_SRC) $(LLMRSSMEFTHiggs_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_MRSSMEFTHiggs)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_MRSSMEFTHiggs): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_MRSSMEFTHiggs)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_MRSSMEFTHiggs)"
		@echo "Note: to regenerate MRSSMEFTHiggs source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_MRSSMEFTHiggs)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_MRSSMEFTHiggs):
		@true
endif

$(LIBMRSSMEFTHiggs_DEP) $(EXEMRSSMEFTHiggs_DEP) $(LLMRSSMEFTHiggs_DEP) $(LIBMRSSMEFTHiggs_OBJ) $(EXEMRSSMEFTHiggs_OBJ) $(LLMRSSMEFTHiggs_OBJ) $(LLMRSSMEFTHiggs_LIB): \
	CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBMRSSMEFTHiggs_DEP) $(EXEMRSSMEFTHiggs_DEP) $(LLMRSSMEFTHiggs_DEP) $(LIBMRSSMEFTHiggs_OBJ) $(EXEMRSSMEFTHiggs_OBJ) $(LLMRSSMEFTHiggs_OBJ) $(LLMRSSMEFTHiggs_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLMRSSMEFTHiggs_OBJ) $(LLMRSSMEFTHiggs_LIB): \
	CPPFLAGS += $(shell $(MATH_INC_PATHS) --math-cmd="$(MATH)" -I --librarylink --mathlink)

$(LIBMRSSMEFTHiggs): $(LIBMRSSMEFTHiggs_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBMRSSMEFTHiggs) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLMRSSMEFTHiggs_LIB): $(LLMRSSMEFTHiggs_OBJ) $(LIBMRSSMEFTHiggs) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

ALLDEP += $(LIBMRSSMEFTHiggs_DEP) $(EXEMRSSMEFTHiggs_DEP)
ALLSRC += $(LIBMRSSMEFTHiggs_SRC) $(EXEMRSSMEFTHiggs_SRC)
ALLLIB += $(LIBMRSSMEFTHiggs)
ALLEXE += $(EXEMRSSMEFTHiggs_EXE)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLMRSSMEFTHiggs_DEP)
ALLSRC += $(LLMRSSMEFTHiggs_SRC)
ALLLL  += $(LLMRSSMEFTHiggs_LIB)
endif
