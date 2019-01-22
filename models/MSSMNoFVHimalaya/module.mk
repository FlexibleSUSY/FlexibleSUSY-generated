DIR          := models/MSSMNoFVHimalaya
MODNAME      := MSSMNoFVHimalaya
SARAH_MODEL  := MSSMNoFV
WITH_$(MODNAME) := yes
MODMSSMNoFVHimalaya_MOD := SM MSSM_higgs MSSM_thresholds
MODMSSMNoFVHimalaya_DEP := $(patsubst %,model_specific/%,$(MODMSSMNoFVHimalaya_MOD))
MODMSSMNoFVHimalaya_INC := $(patsubst %,-Imodel_specific/%,$(MODMSSMNoFVHimalaya_MOD))
MODMSSMNoFVHimalaya_LIB := $(foreach M,$(MODMSSMNoFVHimalaya_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

MODMSSMNoFVHimalaya_SUBMOD  := $(DIR)/cxx_qft
MODMSSMNoFVHimalaya_SUBMOD_INC := $(patsubst %,-I%,$(MODMSSMNoFVHimalaya_SUBMOD))

MSSMNoFVHimalaya_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)
MSSMNoFVHimalaya_INSTALL_CXXQFT_DIR := \
		$(MSSMNoFVHimalaya_INSTALL_DIR)/cxx_qft

MSSMNoFVHimalaya_MK     := \
		$(DIR)/module.mk

MSSMNoFVHimalaya_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

MSSMNoFVHimalaya_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

MSSMNoFVHimalaya_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

MSSMNoFVHimalaya_INCLUDE_MK := \
		$(MSSMNoFVHimalaya_SUSY_BETAS_MK) \
		$(MSSMNoFVHimalaya_SOFT_BETAS_MK)

MSSMNoFVHimalaya_SLHA_INPUT := \
		$(DIR)/LesHouches.in.MSSMNoFVHimalaya_generated \
		$(DIR)/LesHouches.in.MSSMNoFVHimalaya

MSSMNoFVHimalaya_REFERENCES := \
		$(DIR)/MSSMNoFVHimalaya_references.tex

MSSMNoFVHimalaya_GNUPLOT := \
		$(DIR)/MSSMNoFVHimalaya_plot_rgflow.gnuplot \
		$(DIR)/MSSMNoFVHimalaya_plot_spectrum.gnuplot

MSSMNoFVHimalaya_TARBALL := \
		$(MODNAME).tar.gz

LIBMSSMNoFVHimalaya_SRC := \
		$(DIR)/MSSMNoFVHimalaya_a_muon.cpp \
		$(DIR)/MSSMNoFVHimalaya_edm.cpp \
		$(DIR)/MSSMNoFVHimalaya_effective_couplings.cpp \
		$(DIR)/MSSMNoFVHimalaya_info.cpp \
		$(DIR)/MSSMNoFVHimalaya_input_parameters.cpp \
		$(DIR)/MSSMNoFVHimalaya_mass_eigenstates.cpp \
		$(DIR)/MSSMNoFVHimalaya_observables.cpp \
		$(DIR)/MSSMNoFVHimalaya_physical.cpp \
		$(DIR)/MSSMNoFVHimalaya_slha_io.cpp \
		$(DIR)/MSSMNoFVHimalaya_soft_parameters.cpp \
		$(DIR)/MSSMNoFVHimalaya_susy_parameters.cpp \
		$(DIR)/MSSMNoFVHimalaya_utilities.cpp \
		$(DIR)/MSSMNoFVHimalaya_weinberg_angle.cpp

EXEMSSMNoFVHimalaya_SRC := \
		$(DIR)/run_MSSMNoFVHimalaya.cpp \
		$(DIR)/run_cmd_line_MSSMNoFVHimalaya.cpp \
		$(DIR)/scan_MSSMNoFVHimalaya.cpp
LLMSSMNoFVHimalaya_LIB  :=
LLMSSMNoFVHimalaya_OBJ  :=
LLMSSMNoFVHimalaya_SRC  := \
		$(DIR)/MSSMNoFVHimalaya_librarylink.cpp

LLMSSMNoFVHimalaya_MMA  := \
		$(DIR)/MSSMNoFVHimalaya_librarylink.m \
		$(DIR)/run_MSSMNoFVHimalaya.m

LIBMSSMNoFVHimalaya_HDR := \
		$(DIR)/MSSMNoFVHimalaya_a_muon.hpp \
		$(DIR)/MSSMNoFVHimalaya_convergence_tester.hpp \
		$(DIR)/MSSMNoFVHimalaya_edm.hpp \
		$(DIR)/MSSMNoFVHimalaya_effective_couplings.hpp \
		$(DIR)/MSSMNoFVHimalaya_ewsb_solver.hpp \
		$(DIR)/MSSMNoFVHimalaya_ewsb_solver_interface.hpp \
		$(DIR)/MSSMNoFVHimalaya_high_scale_constraint.hpp \
		$(DIR)/MSSMNoFVHimalaya_info.hpp \
		$(DIR)/MSSMNoFVHimalaya_initial_guesser.hpp \
		$(DIR)/MSSMNoFVHimalaya_input_parameters.hpp \
		$(DIR)/MSSMNoFVHimalaya_low_scale_constraint.hpp \
		$(DIR)/MSSMNoFVHimalaya_mass_eigenstates.hpp \
		$(DIR)/MSSMNoFVHimalaya_model.hpp \
		$(DIR)/MSSMNoFVHimalaya_model_slha.hpp \
		$(DIR)/MSSMNoFVHimalaya_observables.hpp \
		$(DIR)/MSSMNoFVHimalaya_physical.hpp \
		$(DIR)/MSSMNoFVHimalaya_slha_io.hpp \
		$(DIR)/MSSMNoFVHimalaya_spectrum_generator.hpp \
		$(DIR)/MSSMNoFVHimalaya_spectrum_generator_interface.hpp \
		$(DIR)/MSSMNoFVHimalaya_soft_parameters.hpp \
		$(DIR)/MSSMNoFVHimalaya_susy_parameters.hpp \
		$(DIR)/MSSMNoFVHimalaya_susy_scale_constraint.hpp \
		$(DIR)/MSSMNoFVHimalaya_utilities.hpp \
		$(DIR)/MSSMNoFVHimalaya_weinberg_angle.hpp

LIBMSSMNoFVHimalaya_CXXQFT_HDR := \
		$(DIR)/cxx_qft/MSSMNoFVHimalaya_qft.hpp \
		$(DIR)/cxx_qft/MSSMNoFVHimalaya_fields.hpp \
		$(DIR)/cxx_qft/MSSMNoFVHimalaya_vertices.hpp \
		$(DIR)/cxx_qft/MSSMNoFVHimalaya_context_base.hpp \
		$(DIR)/cxx_qft/MSSMNoFVHimalaya_npointfunctions.hpp

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
-include $(MSSMNoFVHimalaya_SUSY_BETAS_MK)
-include $(MSSMNoFVHimalaya_SOFT_BETAS_MK)
-include $(MSSMNoFVHimalaya_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(MSSMNoFVHimalaya_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(MSSMNoFVHimalaya_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(MSSMNoFVHimalaya_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBMSSMNoFVHimalaya_SRC := $(sort $(LIBMSSMNoFVHimalaya_SRC))
EXEMSSMNoFVHimalaya_SRC := $(sort $(EXEMSSMNoFVHimalaya_SRC))

LIBMSSMNoFVHimalaya_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBMSSMNoFVHimalaya_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBMSSMNoFVHimalaya_SRC)))

EXEMSSMNoFVHimalaya_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEMSSMNoFVHimalaya_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEMSSMNoFVHimalaya_SRC)))

EXEMSSMNoFVHimalaya_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXEMSSMNoFVHimalaya_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXEMSSMNoFVHimalaya_SRC)))

LIBMSSMNoFVHimalaya_DEP := \
		$(LIBMSSMNoFVHimalaya_OBJ:.o=.d)

EXEMSSMNoFVHimalaya_DEP := \
		$(EXEMSSMNoFVHimalaya_OBJ:.o=.d)

LLMSSMNoFVHimalaya_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLMSSMNoFVHimalaya_SRC)))

LLMSSMNoFVHimalaya_OBJ  := $(LLMSSMNoFVHimalaya_SRC:.cpp=.o)
LLMSSMNoFVHimalaya_LIB  := $(LLMSSMNoFVHimalaya_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBMSSMNoFVHimalaya     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_MSSMNoFVHimalaya := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_MSSMNoFVHimalaya := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL_DIR) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBMSSMNoFVHimalaya) $(EXEMSSMNoFVHimalaya_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(MSSMNoFVHimalaya_INSTALL_DIR)
		install -d $(MSSMNoFVHimalaya_INSTALL_CXXQFT_DIR)
		install -m u=rw,g=r,o=r $(LIBMSSMNoFVHimalaya_SRC) $(MSSMNoFVHimalaya_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBMSSMNoFVHimalaya_HDR) $(MSSMNoFVHimalaya_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBMSSMNoFVHimalaya_CXXQFT_HDR) $(MSSMNoFVHimalaya_INSTALL_CXXQFT_DIR)
		install -m u=rw,g=r,o=r $(EXEMSSMNoFVHimalaya_SRC) $(MSSMNoFVHimalaya_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLMSSMNoFVHimalaya_SRC) $(MSSMNoFVHimalaya_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLMSSMNoFVHimalaya_MMA) $(MSSMNoFVHimalaya_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(MSSMNoFVHimalaya_MK) $(MSSMNoFVHimalaya_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(MSSMNoFVHimalaya_INCLUDE_MK) $(MSSMNoFVHimalaya_INSTALL_DIR)
ifneq ($(MSSMNoFVHimalaya_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(MSSMNoFVHimalaya_SLHA_INPUT) $(MSSMNoFVHimalaya_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(MSSMNoFVHimalaya_REFERENCES) $(MSSMNoFVHimalaya_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(MSSMNoFVHimalaya_GNUPLOT) $(MSSMNoFVHimalaya_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBMSSMNoFVHimalaya_DEP)
		-rm -f $(EXEMSSMNoFVHimalaya_DEP)
		-rm -f $(LLMSSMNoFVHimalaya_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBMSSMNoFVHimalaya)
		-rm -f $(LLMSSMNoFVHimalaya_LIB)

clean-$(MODNAME)-obj:
		-rm -f $(LIBMSSMNoFVHimalaya_OBJ)
		-rm -f $(EXEMSSMNoFVHimalaya_OBJ)
		-rm -f $(LLMSSMNoFVHimalaya_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		-rm -f $(LIBMSSMNoFVHimalaya_SRC)
		-rm -f $(LIBMSSMNoFVHimalaya_HDR)
		-rm -f $(LIBMSSMNoFVHimalaya_CXXQFT_HDR)
		-rm -f $(EXEMSSMNoFVHimalaya_SRC)
		-rm -f $(LLMSSMNoFVHimalaya_SRC)
		-rm -f $(LLMSSMNoFVHimalaya_MMA)
		-rm -f $(METACODE_STAMP_MSSMNoFVHimalaya)
		-rm -f $(MSSMNoFVHimalaya_INCLUDE_MK)
		-rm -f $(MSSMNoFVHimalaya_SLHA_INPUT)
		-rm -f $(MSSMNoFVHimalaya_REFERENCES)
		-rm -f $(MSSMNoFVHimalaya_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		-rm -f $(EXEMSSMNoFVHimalaya_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(MSSMNoFVHimalaya_TARBALL) \
		$(LIBMSSMNoFVHimalaya_SRC) $(LIBMSSMNoFVHimalaya_HDR) $(LIBMSSMNoFVHimalaya_CXXQFT_HDR) \
		$(EXEMSSMNoFVHimalaya_SRC) \
		$(LLMSSMNoFVHimalaya_SRC) $(LLMSSMNoFVHimalaya_MMA) \
		$(MSSMNoFVHimalaya_MK) $(MSSMNoFVHimalaya_INCLUDE_MK) \
		$(MSSMNoFVHimalaya_SLHA_INPUT) $(MSSMNoFVHimalaya_REFERENCES) \
		$(MSSMNoFVHimalaya_GNUPLOT)

$(LIBMSSMNoFVHimalaya_SRC) $(LIBMSSMNoFVHimalaya_HDR) $(LIBMSSMNoFVHimalaya_CXXQFT_HDR) $(EXEMSSMNoFVHimalaya_SRC) $(LLMSSMNoFVHimalaya_SRC) $(LLMSSMNoFVHimalaya_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_MSSMNoFVHimalaya)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_MSSMNoFVHimalaya): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_MSSMNoFVHimalaya)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]" || (echo "Error: The code generation failed!"; exit 1)
		@touch "$(METACODE_STAMP_MSSMNoFVHimalaya)"
		@echo "Note: to regenerate MSSMNoFVHimalaya source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_MSSMNoFVHimalaya)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_MSSMNoFVHimalaya):
		@true
endif

$(LIBMSSMNoFVHimalaya_DEP) $(EXEMSSMNoFVHimalaya_DEP) $(LLMSSMNoFVHimalaya_DEP) $(LIBMSSMNoFVHimalaya_OBJ) $(EXEMSSMNoFVHimalaya_OBJ) $(LLMSSMNoFVHimalaya_OBJ) $(LLMSSMNoFVHimalaya_LIB): \
	CPPFLAGS += $(MODMSSMNoFVHimalaya_SUBMOD_INC) $(MODMSSMNoFVHimalaya_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBMSSMNoFVHimalaya_DEP) $(EXEMSSMNoFVHimalaya_DEP) $(LLMSSMNoFVHimalaya_DEP) $(LIBMSSMNoFVHimalaya_OBJ) $(EXEMSSMNoFVHimalaya_OBJ) $(LLMSSMNoFVHimalaya_OBJ) $(LLMSSMNoFVHimalaya_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLMSSMNoFVHimalaya_OBJ) $(LLMSSMNoFVHimalaya_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBMSSMNoFVHimalaya): $(LIBMSSMNoFVHimalaya_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBMSSMNoFVHimalaya) $(MODMSSMNoFVHimalaya_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLMSSMNoFVHimalaya_LIB): $(LLMSSMNoFVHimalaya_OBJ) $(LIBMSSMNoFVHimalaya) $(MODMSSMNoFVHimalaya_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS)

ALLDEP += $(LIBMSSMNoFVHimalaya_DEP) $(EXEMSSMNoFVHimalaya_DEP)
ALLSRC += $(LIBMSSMNoFVHimalaya_SRC) $(EXEMSSMNoFVHimalaya_SRC)
ALLLIB += $(LIBMSSMNoFVHimalaya)
ALLEXE += $(EXEMSSMNoFVHimalaya_EXE)
ALLMODDEP += $(MODMSSMNoFVHimalaya_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLMSSMNoFVHimalaya_DEP)
ALLSRC += $(LLMSSMNoFVHimalaya_SRC)
ALLLL  += $(LLMSSMNoFVHimalaya_LIB)
endif
