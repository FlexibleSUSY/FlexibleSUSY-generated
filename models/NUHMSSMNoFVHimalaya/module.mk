DIR          := models/NUHMSSMNoFVHimalaya
MODNAME      := NUHMSSMNoFVHimalaya
SARAH_MODEL  := MSSMNoFV
WITH_$(MODNAME) := yes
MODNUHMSSMNoFVHimalaya_MOD := SM MSSM_higgs MSSM_thresholds
MODNUHMSSMNoFVHimalaya_DEP := $(patsubst %,model_specific/%,$(MODNUHMSSMNoFVHimalaya_MOD))
MODNUHMSSMNoFVHimalaya_INC := $(patsubst %,-Imodel_specific/%,$(MODNUHMSSMNoFVHimalaya_MOD))
MODNUHMSSMNoFVHimalaya_LIB := $(foreach M,$(MODNUHMSSMNoFVHimalaya_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

NUHMSSMNoFVHimalaya_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

NUHMSSMNoFVHimalaya_MK     := \
		$(DIR)/module.mk

NUHMSSMNoFVHimalaya_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

NUHMSSMNoFVHimalaya_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

NUHMSSMNoFVHimalaya_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

NUHMSSMNoFVHimalaya_INCLUDE_MK := \
		$(NUHMSSMNoFVHimalaya_SUSY_BETAS_MK) \
		$(NUHMSSMNoFVHimalaya_SOFT_BETAS_MK)

NUHMSSMNoFVHimalaya_SLHA_INPUT := \
		$(DIR)/LesHouches.in.NUHMSSMNoFVHimalaya_generated \
		$(DIR)/LesHouches.in.NUHMSSMNoFVHimalaya

NUHMSSMNoFVHimalaya_REFERENCES := \
		$(DIR)/NUHMSSMNoFVHimalaya_references.tex

NUHMSSMNoFVHimalaya_GNUPLOT := \
		$(DIR)/NUHMSSMNoFVHimalaya_plot_rgflow.gnuplot \
		$(DIR)/NUHMSSMNoFVHimalaya_plot_spectrum.gnuplot

NUHMSSMNoFVHimalaya_TARBALL := \
		$(MODNAME).tar.gz

LIBNUHMSSMNoFVHimalaya_SRC := \
		$(DIR)/NUHMSSMNoFVHimalaya_a_muon.cpp \
		$(DIR)/NUHMSSMNoFVHimalaya_edm.cpp \
		$(DIR)/NUHMSSMNoFVHimalaya_effective_couplings.cpp \
		$(DIR)/NUHMSSMNoFVHimalaya_info.cpp \
		$(DIR)/NUHMSSMNoFVHimalaya_input_parameters.cpp \
		$(DIR)/NUHMSSMNoFVHimalaya_mass_eigenstates.cpp \
		$(DIR)/NUHMSSMNoFVHimalaya_observables.cpp \
		$(DIR)/NUHMSSMNoFVHimalaya_physical.cpp \
		$(DIR)/NUHMSSMNoFVHimalaya_slha_io.cpp \
		$(DIR)/NUHMSSMNoFVHimalaya_soft_parameters.cpp \
		$(DIR)/NUHMSSMNoFVHimalaya_susy_parameters.cpp \
		$(DIR)/NUHMSSMNoFVHimalaya_utilities.cpp \
		$(DIR)/NUHMSSMNoFVHimalaya_weinberg_angle.cpp

EXENUHMSSMNoFVHimalaya_SRC := \
		$(DIR)/run_NUHMSSMNoFVHimalaya.cpp \
		$(DIR)/run_cmd_line_NUHMSSMNoFVHimalaya.cpp \
		$(DIR)/scan_NUHMSSMNoFVHimalaya.cpp
LLNUHMSSMNoFVHimalaya_LIB  :=
LLNUHMSSMNoFVHimalaya_OBJ  :=
LLNUHMSSMNoFVHimalaya_SRC  := \
		$(DIR)/NUHMSSMNoFVHimalaya_librarylink.cpp

LLNUHMSSMNoFVHimalaya_MMA  := \
		$(DIR)/NUHMSSMNoFVHimalaya_librarylink.m \
		$(DIR)/run_NUHMSSMNoFVHimalaya.m

LIBNUHMSSMNoFVHimalaya_HDR := \
		$(DIR)/NUHMSSMNoFVHimalaya_cxx_diagrams.hpp \
		$(DIR)/NUHMSSMNoFVHimalaya_a_muon.hpp \
		$(DIR)/NUHMSSMNoFVHimalaya_convergence_tester.hpp \
		$(DIR)/NUHMSSMNoFVHimalaya_edm.hpp \
		$(DIR)/NUHMSSMNoFVHimalaya_effective_couplings.hpp \
		$(DIR)/NUHMSSMNoFVHimalaya_ewsb_solver.hpp \
		$(DIR)/NUHMSSMNoFVHimalaya_ewsb_solver_interface.hpp \
		$(DIR)/NUHMSSMNoFVHimalaya_high_scale_constraint.hpp \
		$(DIR)/NUHMSSMNoFVHimalaya_info.hpp \
		$(DIR)/NUHMSSMNoFVHimalaya_initial_guesser.hpp \
		$(DIR)/NUHMSSMNoFVHimalaya_input_parameters.hpp \
		$(DIR)/NUHMSSMNoFVHimalaya_low_scale_constraint.hpp \
		$(DIR)/NUHMSSMNoFVHimalaya_mass_eigenstates.hpp \
		$(DIR)/NUHMSSMNoFVHimalaya_model.hpp \
		$(DIR)/NUHMSSMNoFVHimalaya_model_slha.hpp \
		$(DIR)/NUHMSSMNoFVHimalaya_observables.hpp \
		$(DIR)/NUHMSSMNoFVHimalaya_physical.hpp \
		$(DIR)/NUHMSSMNoFVHimalaya_slha_io.hpp \
		$(DIR)/NUHMSSMNoFVHimalaya_spectrum_generator.hpp \
		$(DIR)/NUHMSSMNoFVHimalaya_spectrum_generator_interface.hpp \
		$(DIR)/NUHMSSMNoFVHimalaya_soft_parameters.hpp \
		$(DIR)/NUHMSSMNoFVHimalaya_susy_parameters.hpp \
		$(DIR)/NUHMSSMNoFVHimalaya_susy_scale_constraint.hpp \
		$(DIR)/NUHMSSMNoFVHimalaya_utilities.hpp \
		$(DIR)/NUHMSSMNoFVHimalaya_weinberg_angle.hpp

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
-include $(NUHMSSMNoFVHimalaya_SUSY_BETAS_MK)
-include $(NUHMSSMNoFVHimalaya_SOFT_BETAS_MK)
-include $(NUHMSSMNoFVHimalaya_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(NUHMSSMNoFVHimalaya_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(NUHMSSMNoFVHimalaya_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(NUHMSSMNoFVHimalaya_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBNUHMSSMNoFVHimalaya_SRC := $(sort $(LIBNUHMSSMNoFVHimalaya_SRC))
EXENUHMSSMNoFVHimalaya_SRC := $(sort $(EXENUHMSSMNoFVHimalaya_SRC))

LIBNUHMSSMNoFVHimalaya_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBNUHMSSMNoFVHimalaya_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBNUHMSSMNoFVHimalaya_SRC)))

EXENUHMSSMNoFVHimalaya_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXENUHMSSMNoFVHimalaya_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXENUHMSSMNoFVHimalaya_SRC)))

EXENUHMSSMNoFVHimalaya_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXENUHMSSMNoFVHimalaya_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXENUHMSSMNoFVHimalaya_SRC)))

LIBNUHMSSMNoFVHimalaya_DEP := \
		$(LIBNUHMSSMNoFVHimalaya_OBJ:.o=.d)

EXENUHMSSMNoFVHimalaya_DEP := \
		$(EXENUHMSSMNoFVHimalaya_OBJ:.o=.d)

LLNUHMSSMNoFVHimalaya_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLNUHMSSMNoFVHimalaya_SRC)))

LLNUHMSSMNoFVHimalaya_OBJ  := $(LLNUHMSSMNoFVHimalaya_SRC:.cpp=.o)
LLNUHMSSMNoFVHimalaya_LIB  := $(LLNUHMSSMNoFVHimalaya_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBNUHMSSMNoFVHimalaya     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_NUHMSSMNoFVHimalaya := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_NUHMSSMNoFVHimalaya := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL_DIR) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBNUHMSSMNoFVHimalaya) $(EXENUHMSSMNoFVHimalaya_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(NUHMSSMNoFVHimalaya_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBNUHMSSMNoFVHimalaya_SRC) $(NUHMSSMNoFVHimalaya_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBNUHMSSMNoFVHimalaya_HDR) $(NUHMSSMNoFVHimalaya_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXENUHMSSMNoFVHimalaya_SRC) $(NUHMSSMNoFVHimalaya_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLNUHMSSMNoFVHimalaya_SRC) $(NUHMSSMNoFVHimalaya_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLNUHMSSMNoFVHimalaya_MMA) $(NUHMSSMNoFVHimalaya_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(NUHMSSMNoFVHimalaya_MK) $(NUHMSSMNoFVHimalaya_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(NUHMSSMNoFVHimalaya_INCLUDE_MK) $(NUHMSSMNoFVHimalaya_INSTALL_DIR)
ifneq ($(NUHMSSMNoFVHimalaya_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(NUHMSSMNoFVHimalaya_SLHA_INPUT) $(NUHMSSMNoFVHimalaya_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(NUHMSSMNoFVHimalaya_REFERENCES) $(NUHMSSMNoFVHimalaya_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(NUHMSSMNoFVHimalaya_GNUPLOT) $(NUHMSSMNoFVHimalaya_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBNUHMSSMNoFVHimalaya_DEP)
		-rm -f $(EXENUHMSSMNoFVHimalaya_DEP)
		-rm -f $(LLNUHMSSMNoFVHimalaya_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBNUHMSSMNoFVHimalaya)
		-rm -f $(LLNUHMSSMNoFVHimalaya_LIB)

clean-$(MODNAME)-obj:
		-rm -f $(LIBNUHMSSMNoFVHimalaya_OBJ)
		-rm -f $(EXENUHMSSMNoFVHimalaya_OBJ)
		-rm -f $(LLNUHMSSMNoFVHimalaya_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		-rm -f $(LIBNUHMSSMNoFVHimalaya_SRC)
		-rm -f $(LIBNUHMSSMNoFVHimalaya_HDR)
		-rm -f $(EXENUHMSSMNoFVHimalaya_SRC)
		-rm -f $(LLNUHMSSMNoFVHimalaya_SRC)
		-rm -f $(LLNUHMSSMNoFVHimalaya_MMA)
		-rm -f $(METACODE_STAMP_NUHMSSMNoFVHimalaya)
		-rm -f $(NUHMSSMNoFVHimalaya_INCLUDE_MK)
		-rm -f $(NUHMSSMNoFVHimalaya_SLHA_INPUT)
		-rm -f $(NUHMSSMNoFVHimalaya_REFERENCES)
		-rm -f $(NUHMSSMNoFVHimalaya_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		-rm -f $(EXENUHMSSMNoFVHimalaya_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(NUHMSSMNoFVHimalaya_TARBALL) \
		$(LIBNUHMSSMNoFVHimalaya_SRC) $(LIBNUHMSSMNoFVHimalaya_HDR) \
		$(EXENUHMSSMNoFVHimalaya_SRC) \
		$(LLNUHMSSMNoFVHimalaya_SRC) $(LLNUHMSSMNoFVHimalaya_MMA) \
		$(NUHMSSMNoFVHimalaya_MK) $(NUHMSSMNoFVHimalaya_INCLUDE_MK) \
		$(NUHMSSMNoFVHimalaya_SLHA_INPUT) $(NUHMSSMNoFVHimalaya_REFERENCES) \
		$(NUHMSSMNoFVHimalaya_GNUPLOT)

$(LIBNUHMSSMNoFVHimalaya_SRC) $(LIBNUHMSSMNoFVHimalaya_HDR) $(EXENUHMSSMNoFVHimalaya_SRC) $(LLNUHMSSMNoFVHimalaya_SRC) $(LLNUHMSSMNoFVHimalaya_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_NUHMSSMNoFVHimalaya)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_NUHMSSMNoFVHimalaya): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_NUHMSSMNoFVHimalaya)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]" || (echo "Error: The code generation failed!"; exit 1)
		@touch "$(METACODE_STAMP_NUHMSSMNoFVHimalaya)"
		@echo "Note: to regenerate NUHMSSMNoFVHimalaya source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_NUHMSSMNoFVHimalaya)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_NUHMSSMNoFVHimalaya):
		@true
endif

$(LIBNUHMSSMNoFVHimalaya_DEP) $(EXENUHMSSMNoFVHimalaya_DEP) $(LLNUHMSSMNoFVHimalaya_DEP) $(LIBNUHMSSMNoFVHimalaya_OBJ) $(EXENUHMSSMNoFVHimalaya_OBJ) $(LLNUHMSSMNoFVHimalaya_OBJ) $(LLNUHMSSMNoFVHimalaya_LIB): \
	CPPFLAGS += $(MODNUHMSSMNoFVHimalaya_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBNUHMSSMNoFVHimalaya_DEP) $(EXENUHMSSMNoFVHimalaya_DEP) $(LLNUHMSSMNoFVHimalaya_DEP) $(LIBNUHMSSMNoFVHimalaya_OBJ) $(EXENUHMSSMNoFVHimalaya_OBJ) $(LLNUHMSSMNoFVHimalaya_OBJ) $(LLNUHMSSMNoFVHimalaya_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLNUHMSSMNoFVHimalaya_OBJ) $(LLNUHMSSMNoFVHimalaya_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBNUHMSSMNoFVHimalaya): $(LIBNUHMSSMNoFVHimalaya_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBNUHMSSMNoFVHimalaya) $(MODNUHMSSMNoFVHimalaya_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLNUHMSSMNoFVHimalaya_LIB): $(LLNUHMSSMNoFVHimalaya_OBJ) $(LIBNUHMSSMNoFVHimalaya) $(MODNUHMSSMNoFVHimalaya_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS)

ALLDEP += $(LIBNUHMSSMNoFVHimalaya_DEP) $(EXENUHMSSMNoFVHimalaya_DEP)
ALLSRC += $(LIBNUHMSSMNoFVHimalaya_SRC) $(EXENUHMSSMNoFVHimalaya_SRC)
ALLLIB += $(LIBNUHMSSMNoFVHimalaya)
ALLEXE += $(EXENUHMSSMNoFVHimalaya_EXE)
ALLMODDEP += $(MODNUHMSSMNoFVHimalaya_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLNUHMSSMNoFVHimalaya_DEP)
ALLSRC += $(LLNUHMSSMNoFVHimalaya_SRC)
ALLLL  += $(LLNUHMSSMNoFVHimalaya_LIB)
endif
