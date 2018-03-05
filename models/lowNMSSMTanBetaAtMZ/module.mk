DIR          := models/lowNMSSMTanBetaAtMZ
MODNAME      := lowNMSSMTanBetaAtMZ
SARAH_MODEL  := NMSSM
WITH_$(MODNAME) := yes
MODlowNMSSMTanBetaAtMZ_MOD := SM MSSM_higgs NMSSM_higgs
MODlowNMSSMTanBetaAtMZ_DEP := $(patsubst %,model_specific/%,$(MODlowNMSSMTanBetaAtMZ_MOD))
MODlowNMSSMTanBetaAtMZ_INC := $(patsubst %,-Imodel_specific/%,$(MODlowNMSSMTanBetaAtMZ_MOD))
MODlowNMSSMTanBetaAtMZ_LIB := $(foreach M,$(MODlowNMSSMTanBetaAtMZ_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

lowNMSSMTanBetaAtMZ_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

lowNMSSMTanBetaAtMZ_MK     := \
		$(DIR)/module.mk

lowNMSSMTanBetaAtMZ_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

lowNMSSMTanBetaAtMZ_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

lowNMSSMTanBetaAtMZ_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

lowNMSSMTanBetaAtMZ_INCLUDE_MK := \
		$(lowNMSSMTanBetaAtMZ_SUSY_BETAS_MK) \
		$(lowNMSSMTanBetaAtMZ_SOFT_BETAS_MK)

lowNMSSMTanBetaAtMZ_SLHA_INPUT := \
		$(DIR)/LesHouches.in.lowNMSSMTanBetaAtMZ_generated \
		$(DIR)/LesHouches.in.lowNMSSMTanBetaAtMZ \
		$(DIR)/LesHouches.in.TP4 \
		$(DIR)/LesHouches.in.TP3 \
		$(DIR)/LesHouches.in.TP5 \
		$(DIR)/LesHouches.in.TP2 \
		$(DIR)/LesHouches.in.TP6 \
		$(DIR)/LesHouches.in.TP1

lowNMSSMTanBetaAtMZ_REFERENCES := \
		$(DIR)/lowNMSSMTanBetaAtMZ_references.tex

lowNMSSMTanBetaAtMZ_GNUPLOT := \
		$(DIR)/lowNMSSMTanBetaAtMZ_plot_rgflow.gnuplot \
		$(DIR)/lowNMSSMTanBetaAtMZ_plot_spectrum.gnuplot

lowNMSSMTanBetaAtMZ_TARBALL := \
		$(MODNAME).tar.gz

LIBlowNMSSMTanBetaAtMZ_SRC := \
		$(DIR)/lowNMSSMTanBetaAtMZ_a_muon.cpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_edm.cpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_effective_couplings.cpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_info.cpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_input_parameters.cpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_mass_eigenstates.cpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_observables.cpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_physical.cpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_slha_io.cpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_soft_parameters.cpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_susy_parameters.cpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_utilities.cpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_weinberg_angle.cpp

EXElowNMSSMTanBetaAtMZ_SRC := \
		$(DIR)/run_lowNMSSMTanBetaAtMZ.cpp \
		$(DIR)/run_cmd_line_lowNMSSMTanBetaAtMZ.cpp \
		$(DIR)/scan_lowNMSSMTanBetaAtMZ.cpp
LLlowNMSSMTanBetaAtMZ_LIB  :=
LLlowNMSSMTanBetaAtMZ_OBJ  :=
LLlowNMSSMTanBetaAtMZ_SRC  := \
		$(DIR)/lowNMSSMTanBetaAtMZ_librarylink.cpp

LLlowNMSSMTanBetaAtMZ_MMA  := \
		$(DIR)/lowNMSSMTanBetaAtMZ_librarylink.m \
		$(DIR)/run_lowNMSSMTanBetaAtMZ.m

LIBlowNMSSMTanBetaAtMZ_HDR := \
		$(DIR)/lowNMSSMTanBetaAtMZ_cxx_diagrams.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_a_muon.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_convergence_tester.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_edm.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_effective_couplings.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_ewsb_solver.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_ewsb_solver_interface.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_high_scale_constraint.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_info.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_initial_guesser.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_input_parameters.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_low_scale_constraint.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_mass_eigenstates.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_model.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_model_slha.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_observables.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_physical.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_slha_io.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_spectrum_generator.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_spectrum_generator_interface.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_soft_parameters.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_susy_parameters.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_susy_scale_constraint.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_utilities.hpp \
		$(DIR)/lowNMSSMTanBetaAtMZ_weinberg_angle.hpp

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
-include $(lowNMSSMTanBetaAtMZ_SUSY_BETAS_MK)
-include $(lowNMSSMTanBetaAtMZ_SOFT_BETAS_MK)
-include $(lowNMSSMTanBetaAtMZ_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(lowNMSSMTanBetaAtMZ_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(lowNMSSMTanBetaAtMZ_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(lowNMSSMTanBetaAtMZ_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBlowNMSSMTanBetaAtMZ_SRC := $(sort $(LIBlowNMSSMTanBetaAtMZ_SRC))
EXElowNMSSMTanBetaAtMZ_SRC := $(sort $(EXElowNMSSMTanBetaAtMZ_SRC))

LIBlowNMSSMTanBetaAtMZ_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBlowNMSSMTanBetaAtMZ_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBlowNMSSMTanBetaAtMZ_SRC)))

EXElowNMSSMTanBetaAtMZ_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXElowNMSSMTanBetaAtMZ_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXElowNMSSMTanBetaAtMZ_SRC)))

EXElowNMSSMTanBetaAtMZ_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXElowNMSSMTanBetaAtMZ_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXElowNMSSMTanBetaAtMZ_SRC)))

LIBlowNMSSMTanBetaAtMZ_DEP := \
		$(LIBlowNMSSMTanBetaAtMZ_OBJ:.o=.d)

EXElowNMSSMTanBetaAtMZ_DEP := \
		$(EXElowNMSSMTanBetaAtMZ_OBJ:.o=.d)

LLlowNMSSMTanBetaAtMZ_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLlowNMSSMTanBetaAtMZ_SRC)))

LLlowNMSSMTanBetaAtMZ_OBJ  := $(LLlowNMSSMTanBetaAtMZ_SRC:.cpp=.o)
LLlowNMSSMTanBetaAtMZ_LIB  := $(LLlowNMSSMTanBetaAtMZ_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBlowNMSSMTanBetaAtMZ     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_lowNMSSMTanBetaAtMZ := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_lowNMSSMTanBetaAtMZ := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBlowNMSSMTanBetaAtMZ) $(EXElowNMSSMTanBetaAtMZ_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(lowNMSSMTanBetaAtMZ_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBlowNMSSMTanBetaAtMZ_SRC) $(lowNMSSMTanBetaAtMZ_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBlowNMSSMTanBetaAtMZ_HDR) $(lowNMSSMTanBetaAtMZ_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXElowNMSSMTanBetaAtMZ_SRC) $(lowNMSSMTanBetaAtMZ_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLlowNMSSMTanBetaAtMZ_SRC) $(lowNMSSMTanBetaAtMZ_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLlowNMSSMTanBetaAtMZ_MMA) $(lowNMSSMTanBetaAtMZ_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(lowNMSSMTanBetaAtMZ_MK) $(lowNMSSMTanBetaAtMZ_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(lowNMSSMTanBetaAtMZ_INCLUDE_MK) $(lowNMSSMTanBetaAtMZ_INSTALL_DIR)
ifneq ($(lowNMSSMTanBetaAtMZ_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(lowNMSSMTanBetaAtMZ_SLHA_INPUT) $(lowNMSSMTanBetaAtMZ_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(lowNMSSMTanBetaAtMZ_REFERENCES) $(lowNMSSMTanBetaAtMZ_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(lowNMSSMTanBetaAtMZ_GNUPLOT) $(lowNMSSMTanBetaAtMZ_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBlowNMSSMTanBetaAtMZ_DEP)
		-rm -f $(EXElowNMSSMTanBetaAtMZ_DEP)
		-rm -f $(LLlowNMSSMTanBetaAtMZ_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBlowNMSSMTanBetaAtMZ)
		-rm -f $(LLlowNMSSMTanBetaAtMZ_LIB)

clean-$(MODNAME)-obj:
		-rm -f $(LIBlowNMSSMTanBetaAtMZ_OBJ)
		-rm -f $(EXElowNMSSMTanBetaAtMZ_OBJ)
		-rm -f $(LLlowNMSSMTanBetaAtMZ_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		-rm -f $(LIBlowNMSSMTanBetaAtMZ_SRC)
		-rm -f $(LIBlowNMSSMTanBetaAtMZ_HDR)
		-rm -f $(EXElowNMSSMTanBetaAtMZ_SRC)
		-rm -f $(LLlowNMSSMTanBetaAtMZ_SRC)
		-rm -f $(LLlowNMSSMTanBetaAtMZ_MMA)
		-rm -f $(METACODE_STAMP_lowNMSSMTanBetaAtMZ)
		-rm -f $(lowNMSSMTanBetaAtMZ_INCLUDE_MK)
		-rm -f $(lowNMSSMTanBetaAtMZ_SLHA_INPUT)
		-rm -f $(lowNMSSMTanBetaAtMZ_REFERENCES)
		-rm -f $(lowNMSSMTanBetaAtMZ_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		-rm -f $(EXElowNMSSMTanBetaAtMZ_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(lowNMSSMTanBetaAtMZ_TARBALL) \
		$(LIBlowNMSSMTanBetaAtMZ_SRC) $(LIBlowNMSSMTanBetaAtMZ_HDR) \
		$(EXElowNMSSMTanBetaAtMZ_SRC) \
		$(LLlowNMSSMTanBetaAtMZ_SRC) $(LLlowNMSSMTanBetaAtMZ_MMA) \
		$(lowNMSSMTanBetaAtMZ_MK) $(lowNMSSMTanBetaAtMZ_INCLUDE_MK) \
		$(lowNMSSMTanBetaAtMZ_SLHA_INPUT) $(lowNMSSMTanBetaAtMZ_REFERENCES) \
		$(lowNMSSMTanBetaAtMZ_GNUPLOT)

$(LIBlowNMSSMTanBetaAtMZ_SRC) $(LIBlowNMSSMTanBetaAtMZ_HDR) $(EXElowNMSSMTanBetaAtMZ_SRC) $(LLlowNMSSMTanBetaAtMZ_SRC) $(LLlowNMSSMTanBetaAtMZ_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_lowNMSSMTanBetaAtMZ)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_lowNMSSMTanBetaAtMZ): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_lowNMSSMTanBetaAtMZ)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_lowNMSSMTanBetaAtMZ)"
		@echo "Note: to regenerate lowNMSSMTanBetaAtMZ source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_lowNMSSMTanBetaAtMZ)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_lowNMSSMTanBetaAtMZ):
		@true
endif

$(LIBlowNMSSMTanBetaAtMZ_DEP) $(EXElowNMSSMTanBetaAtMZ_DEP) $(LLlowNMSSMTanBetaAtMZ_DEP) $(LIBlowNMSSMTanBetaAtMZ_OBJ) $(EXElowNMSSMTanBetaAtMZ_OBJ) $(LLlowNMSSMTanBetaAtMZ_OBJ) $(LLlowNMSSMTanBetaAtMZ_LIB): \
	CPPFLAGS += $(MODlowNMSSMTanBetaAtMZ_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBlowNMSSMTanBetaAtMZ_DEP) $(EXElowNMSSMTanBetaAtMZ_DEP) $(LLlowNMSSMTanBetaAtMZ_DEP) $(LIBlowNMSSMTanBetaAtMZ_OBJ) $(EXElowNMSSMTanBetaAtMZ_OBJ) $(LLlowNMSSMTanBetaAtMZ_OBJ) $(LLlowNMSSMTanBetaAtMZ_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLlowNMSSMTanBetaAtMZ_OBJ) $(LLlowNMSSMTanBetaAtMZ_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBlowNMSSMTanBetaAtMZ): $(LIBlowNMSSMTanBetaAtMZ_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBlowNMSSMTanBetaAtMZ) $(MODlowNMSSMTanBetaAtMZ_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLlowNMSSMTanBetaAtMZ_LIB): $(LLlowNMSSMTanBetaAtMZ_OBJ) $(LIBlowNMSSMTanBetaAtMZ) $(MODlowNMSSMTanBetaAtMZ_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS)

ALLDEP += $(LIBlowNMSSMTanBetaAtMZ_DEP) $(EXElowNMSSMTanBetaAtMZ_DEP)
ALLSRC += $(LIBlowNMSSMTanBetaAtMZ_SRC) $(EXElowNMSSMTanBetaAtMZ_SRC)
ALLLIB += $(LIBlowNMSSMTanBetaAtMZ)
ALLEXE += $(EXElowNMSSMTanBetaAtMZ_EXE)
ALLMODDEP += $(MODlowNMSSMTanBetaAtMZ_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLlowNMSSMTanBetaAtMZ_DEP)
ALLSRC += $(LLlowNMSSMTanBetaAtMZ_SRC)
ALLLL  += $(LLlowNMSSMTanBetaAtMZ_LIB)
endif
