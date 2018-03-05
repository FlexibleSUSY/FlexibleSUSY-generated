DIR          := models/MSSMNoFVatMGUTHimalaya
MODNAME      := MSSMNoFVatMGUTHimalaya
SARAH_MODEL  := MSSMNoFV
WITH_$(MODNAME) := yes
MODMSSMNoFVatMGUTHimalaya_MOD := SM MSSM_higgs MSSM_thresholds
MODMSSMNoFVatMGUTHimalaya_DEP := $(patsubst %,model_specific/%,$(MODMSSMNoFVatMGUTHimalaya_MOD))
MODMSSMNoFVatMGUTHimalaya_INC := $(patsubst %,-Imodel_specific/%,$(MODMSSMNoFVatMGUTHimalaya_MOD))
MODMSSMNoFVatMGUTHimalaya_LIB := $(foreach M,$(MODMSSMNoFVatMGUTHimalaya_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

MSSMNoFVatMGUTHimalaya_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

MSSMNoFVatMGUTHimalaya_MK     := \
		$(DIR)/module.mk

MSSMNoFVatMGUTHimalaya_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

MSSMNoFVatMGUTHimalaya_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

MSSMNoFVatMGUTHimalaya_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

MSSMNoFVatMGUTHimalaya_INCLUDE_MK := \
		$(MSSMNoFVatMGUTHimalaya_SUSY_BETAS_MK) \
		$(MSSMNoFVatMGUTHimalaya_SOFT_BETAS_MK)

MSSMNoFVatMGUTHimalaya_SLHA_INPUT := \
		$(DIR)/LesHouches.in.MSSMNoFVatMGUTHimalaya_generated \
		$(DIR)/LesHouches.in.MSSMNoFVatMGUTHimalaya

MSSMNoFVatMGUTHimalaya_REFERENCES := \
		$(DIR)/MSSMNoFVatMGUTHimalaya_references.tex

MSSMNoFVatMGUTHimalaya_GNUPLOT := \
		$(DIR)/MSSMNoFVatMGUTHimalaya_plot_rgflow.gnuplot \
		$(DIR)/MSSMNoFVatMGUTHimalaya_plot_spectrum.gnuplot

MSSMNoFVatMGUTHimalaya_TARBALL := \
		$(MODNAME).tar.gz

LIBMSSMNoFVatMGUTHimalaya_SRC := \
		$(DIR)/MSSMNoFVatMGUTHimalaya_a_muon.cpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_edm.cpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_effective_couplings.cpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_info.cpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_input_parameters.cpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_mass_eigenstates.cpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_observables.cpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_physical.cpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_slha_io.cpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_soft_parameters.cpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_susy_parameters.cpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_utilities.cpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_weinberg_angle.cpp

EXEMSSMNoFVatMGUTHimalaya_SRC := \
		$(DIR)/run_MSSMNoFVatMGUTHimalaya.cpp \
		$(DIR)/run_cmd_line_MSSMNoFVatMGUTHimalaya.cpp \
		$(DIR)/scan_MSSMNoFVatMGUTHimalaya.cpp
LLMSSMNoFVatMGUTHimalaya_LIB  :=
LLMSSMNoFVatMGUTHimalaya_OBJ  :=
LLMSSMNoFVatMGUTHimalaya_SRC  := \
		$(DIR)/MSSMNoFVatMGUTHimalaya_librarylink.cpp

LLMSSMNoFVatMGUTHimalaya_MMA  := \
		$(DIR)/MSSMNoFVatMGUTHimalaya_librarylink.m \
		$(DIR)/run_MSSMNoFVatMGUTHimalaya.m

LIBMSSMNoFVatMGUTHimalaya_HDR := \
		$(DIR)/MSSMNoFVatMGUTHimalaya_cxx_diagrams.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_a_muon.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_convergence_tester.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_edm.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_effective_couplings.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_ewsb_solver.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_ewsb_solver_interface.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_high_scale_constraint.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_info.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_initial_guesser.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_input_parameters.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_low_scale_constraint.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_mass_eigenstates.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_model.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_model_slha.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_observables.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_physical.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_slha_io.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_spectrum_generator.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_spectrum_generator_interface.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_soft_parameters.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_susy_parameters.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_susy_scale_constraint.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_utilities.hpp \
		$(DIR)/MSSMNoFVatMGUTHimalaya_weinberg_angle.hpp

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
-include $(MSSMNoFVatMGUTHimalaya_SUSY_BETAS_MK)
-include $(MSSMNoFVatMGUTHimalaya_SOFT_BETAS_MK)
-include $(MSSMNoFVatMGUTHimalaya_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(MSSMNoFVatMGUTHimalaya_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(MSSMNoFVatMGUTHimalaya_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(MSSMNoFVatMGUTHimalaya_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBMSSMNoFVatMGUTHimalaya_SRC := $(sort $(LIBMSSMNoFVatMGUTHimalaya_SRC))
EXEMSSMNoFVatMGUTHimalaya_SRC := $(sort $(EXEMSSMNoFVatMGUTHimalaya_SRC))

LIBMSSMNoFVatMGUTHimalaya_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBMSSMNoFVatMGUTHimalaya_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBMSSMNoFVatMGUTHimalaya_SRC)))

EXEMSSMNoFVatMGUTHimalaya_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEMSSMNoFVatMGUTHimalaya_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEMSSMNoFVatMGUTHimalaya_SRC)))

EXEMSSMNoFVatMGUTHimalaya_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXEMSSMNoFVatMGUTHimalaya_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXEMSSMNoFVatMGUTHimalaya_SRC)))

LIBMSSMNoFVatMGUTHimalaya_DEP := \
		$(LIBMSSMNoFVatMGUTHimalaya_OBJ:.o=.d)

EXEMSSMNoFVatMGUTHimalaya_DEP := \
		$(EXEMSSMNoFVatMGUTHimalaya_OBJ:.o=.d)

LLMSSMNoFVatMGUTHimalaya_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLMSSMNoFVatMGUTHimalaya_SRC)))

LLMSSMNoFVatMGUTHimalaya_OBJ  := $(LLMSSMNoFVatMGUTHimalaya_SRC:.cpp=.o)
LLMSSMNoFVatMGUTHimalaya_LIB  := $(LLMSSMNoFVatMGUTHimalaya_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBMSSMNoFVatMGUTHimalaya     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_MSSMNoFVatMGUTHimalaya := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_MSSMNoFVatMGUTHimalaya := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBMSSMNoFVatMGUTHimalaya) $(EXEMSSMNoFVatMGUTHimalaya_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(MSSMNoFVatMGUTHimalaya_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBMSSMNoFVatMGUTHimalaya_SRC) $(MSSMNoFVatMGUTHimalaya_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBMSSMNoFVatMGUTHimalaya_HDR) $(MSSMNoFVatMGUTHimalaya_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXEMSSMNoFVatMGUTHimalaya_SRC) $(MSSMNoFVatMGUTHimalaya_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLMSSMNoFVatMGUTHimalaya_SRC) $(MSSMNoFVatMGUTHimalaya_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLMSSMNoFVatMGUTHimalaya_MMA) $(MSSMNoFVatMGUTHimalaya_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(MSSMNoFVatMGUTHimalaya_MK) $(MSSMNoFVatMGUTHimalaya_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(MSSMNoFVatMGUTHimalaya_INCLUDE_MK) $(MSSMNoFVatMGUTHimalaya_INSTALL_DIR)
ifneq ($(MSSMNoFVatMGUTHimalaya_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(MSSMNoFVatMGUTHimalaya_SLHA_INPUT) $(MSSMNoFVatMGUTHimalaya_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(MSSMNoFVatMGUTHimalaya_REFERENCES) $(MSSMNoFVatMGUTHimalaya_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(MSSMNoFVatMGUTHimalaya_GNUPLOT) $(MSSMNoFVatMGUTHimalaya_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBMSSMNoFVatMGUTHimalaya_DEP)
		-rm -f $(EXEMSSMNoFVatMGUTHimalaya_DEP)
		-rm -f $(LLMSSMNoFVatMGUTHimalaya_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBMSSMNoFVatMGUTHimalaya)
		-rm -f $(LLMSSMNoFVatMGUTHimalaya_LIB)

clean-$(MODNAME)-obj:
		-rm -f $(LIBMSSMNoFVatMGUTHimalaya_OBJ)
		-rm -f $(EXEMSSMNoFVatMGUTHimalaya_OBJ)
		-rm -f $(LLMSSMNoFVatMGUTHimalaya_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		-rm -f $(LIBMSSMNoFVatMGUTHimalaya_SRC)
		-rm -f $(LIBMSSMNoFVatMGUTHimalaya_HDR)
		-rm -f $(EXEMSSMNoFVatMGUTHimalaya_SRC)
		-rm -f $(LLMSSMNoFVatMGUTHimalaya_SRC)
		-rm -f $(LLMSSMNoFVatMGUTHimalaya_MMA)
		-rm -f $(METACODE_STAMP_MSSMNoFVatMGUTHimalaya)
		-rm -f $(MSSMNoFVatMGUTHimalaya_INCLUDE_MK)
		-rm -f $(MSSMNoFVatMGUTHimalaya_SLHA_INPUT)
		-rm -f $(MSSMNoFVatMGUTHimalaya_REFERENCES)
		-rm -f $(MSSMNoFVatMGUTHimalaya_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		-rm -f $(EXEMSSMNoFVatMGUTHimalaya_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(MSSMNoFVatMGUTHimalaya_TARBALL) \
		$(LIBMSSMNoFVatMGUTHimalaya_SRC) $(LIBMSSMNoFVatMGUTHimalaya_HDR) \
		$(EXEMSSMNoFVatMGUTHimalaya_SRC) \
		$(LLMSSMNoFVatMGUTHimalaya_SRC) $(LLMSSMNoFVatMGUTHimalaya_MMA) \
		$(MSSMNoFVatMGUTHimalaya_MK) $(MSSMNoFVatMGUTHimalaya_INCLUDE_MK) \
		$(MSSMNoFVatMGUTHimalaya_SLHA_INPUT) $(MSSMNoFVatMGUTHimalaya_REFERENCES) \
		$(MSSMNoFVatMGUTHimalaya_GNUPLOT)

$(LIBMSSMNoFVatMGUTHimalaya_SRC) $(LIBMSSMNoFVatMGUTHimalaya_HDR) $(EXEMSSMNoFVatMGUTHimalaya_SRC) $(LLMSSMNoFVatMGUTHimalaya_SRC) $(LLMSSMNoFVatMGUTHimalaya_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_MSSMNoFVatMGUTHimalaya)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_MSSMNoFVatMGUTHimalaya): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_MSSMNoFVatMGUTHimalaya)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_MSSMNoFVatMGUTHimalaya)"
		@echo "Note: to regenerate MSSMNoFVatMGUTHimalaya source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_MSSMNoFVatMGUTHimalaya)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_MSSMNoFVatMGUTHimalaya):
		@true
endif

$(LIBMSSMNoFVatMGUTHimalaya_DEP) $(EXEMSSMNoFVatMGUTHimalaya_DEP) $(LLMSSMNoFVatMGUTHimalaya_DEP) $(LIBMSSMNoFVatMGUTHimalaya_OBJ) $(EXEMSSMNoFVatMGUTHimalaya_OBJ) $(LLMSSMNoFVatMGUTHimalaya_OBJ) $(LLMSSMNoFVatMGUTHimalaya_LIB): \
	CPPFLAGS += $(MODMSSMNoFVatMGUTHimalaya_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBMSSMNoFVatMGUTHimalaya_DEP) $(EXEMSSMNoFVatMGUTHimalaya_DEP) $(LLMSSMNoFVatMGUTHimalaya_DEP) $(LIBMSSMNoFVatMGUTHimalaya_OBJ) $(EXEMSSMNoFVatMGUTHimalaya_OBJ) $(LLMSSMNoFVatMGUTHimalaya_OBJ) $(LLMSSMNoFVatMGUTHimalaya_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLMSSMNoFVatMGUTHimalaya_OBJ) $(LLMSSMNoFVatMGUTHimalaya_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBMSSMNoFVatMGUTHimalaya): $(LIBMSSMNoFVatMGUTHimalaya_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBMSSMNoFVatMGUTHimalaya) $(MODMSSMNoFVatMGUTHimalaya_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLMSSMNoFVatMGUTHimalaya_LIB): $(LLMSSMNoFVatMGUTHimalaya_OBJ) $(LIBMSSMNoFVatMGUTHimalaya) $(MODMSSMNoFVatMGUTHimalaya_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS)

ALLDEP += $(LIBMSSMNoFVatMGUTHimalaya_DEP) $(EXEMSSMNoFVatMGUTHimalaya_DEP)
ALLSRC += $(LIBMSSMNoFVatMGUTHimalaya_SRC) $(EXEMSSMNoFVatMGUTHimalaya_SRC)
ALLLIB += $(LIBMSSMNoFVatMGUTHimalaya)
ALLEXE += $(EXEMSSMNoFVatMGUTHimalaya_EXE)
ALLMODDEP += $(MODMSSMNoFVatMGUTHimalaya_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLMSSMNoFVatMGUTHimalaya_DEP)
ALLSRC += $(LLMSSMNoFVatMGUTHimalaya_SRC)
ALLLL  += $(LLMSSMNoFVatMGUTHimalaya_LIB)
endif