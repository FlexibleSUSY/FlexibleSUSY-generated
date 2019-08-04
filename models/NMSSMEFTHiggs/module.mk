DIR          := models/NMSSMEFTHiggs
MODNAME      := NMSSMEFTHiggs
SARAH_MODEL  := NMSSM
WITH_$(MODNAME) := yes
MODNMSSMEFTHiggs_MOD := SM MSSM_higgs NMSSM_higgs
MODNMSSMEFTHiggs_DEP := $(patsubst %,model_specific/%,$(MODNMSSMEFTHiggs_MOD))
MODNMSSMEFTHiggs_INC := $(patsubst %,-Imodel_specific/%,$(MODNMSSMEFTHiggs_MOD))
MODNMSSMEFTHiggs_LIB := $(foreach M,$(MODNMSSMEFTHiggs_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

MODNMSSMEFTHiggs_SUBMOD  := $(DIR)/cxx_qft
MODNMSSMEFTHiggs_SUBMOD_INC := $(patsubst %,-I%,$(MODNMSSMEFTHiggs_SUBMOD))

NMSSMEFTHiggs_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)
NMSSMEFTHiggs_INSTALL_CXXQFT_DIR := \
		$(NMSSMEFTHiggs_INSTALL_DIR)/cxx_qft

NMSSMEFTHiggs_MK     := \
		$(DIR)/module.mk

NMSSMEFTHiggs_SUSY_BETAS_MK := \
		$(DIR)/susy_betas.mk

NMSSMEFTHiggs_SOFT_BETAS_MK := \
		$(DIR)/soft_betas.mk

NMSSMEFTHiggs_CXX_QFT_VERTICES_MK := \
		$(DIR)/cxx_qft/vertices.mk

NMSSMEFTHiggs_FlexibleEFTHiggs_MK := \
		$(DIR)/FlexibleEFTHiggs.mk

NMSSMEFTHiggs_INCLUDE_MK := \
		$(NMSSMEFTHiggs_SUSY_BETAS_MK) \
		$(NMSSMEFTHiggs_SOFT_BETAS_MK) \
		$(NMSSMEFTHiggs_CXX_QFT_VERTICES_MK)

NMSSMEFTHiggs_SLHA_INPUT := \
		$(DIR)/LesHouches.in.NMSSMEFTHiggs_generated \
		$(DIR)/LesHouches.in.NMSSMEFTHiggs \
		$(DIR)/LesHouches.in.NMSSMEFTHiggs_1507.05093_TP3

NMSSMEFTHiggs_REFERENCES := \
		$(DIR)/NMSSMEFTHiggs_references.tex

NMSSMEFTHiggs_GNUPLOT := \
		$(DIR)/NMSSMEFTHiggs_plot_rgflow.gnuplot \
		$(DIR)/NMSSMEFTHiggs_plot_spectrum.gnuplot

NMSSMEFTHiggs_TARBALL := \
		$(MODNAME).tar.gz

LIBNMSSMEFTHiggs_SRC := \
		$(DIR)/NMSSMEFTHiggs_a_muon.cpp \
		$(DIR)/NMSSMEFTHiggs_edm.cpp \
		$(DIR)/NMSSMEFTHiggs_FFV_form_factors.cpp \
		$(DIR)/NMSSMEFTHiggs_l_to_lgamma.cpp \
		$(DIR)/NMSSMEFTHiggs_effective_couplings.cpp \
		$(DIR)/NMSSMEFTHiggs_info.cpp \
		$(DIR)/NMSSMEFTHiggs_input_parameters.cpp \
		$(DIR)/NMSSMEFTHiggs_mass_eigenstates.cpp \
		$(DIR)/NMSSMEFTHiggs_observables.cpp \
		$(DIR)/NMSSMEFTHiggs_physical.cpp \
		$(DIR)/NMSSMEFTHiggs_slha_io.cpp \
		$(DIR)/NMSSMEFTHiggs_soft_parameters.cpp \
		$(DIR)/NMSSMEFTHiggs_susy_parameters.cpp \
		$(DIR)/NMSSMEFTHiggs_utilities.cpp \
		$(DIR)/NMSSMEFTHiggs_weinberg_angle.cpp

EXENMSSMEFTHiggs_SRC := \
		$(DIR)/run_NMSSMEFTHiggs.cpp \
		$(DIR)/run_cmd_line_NMSSMEFTHiggs.cpp \
		$(DIR)/scan_NMSSMEFTHiggs.cpp
LLNMSSMEFTHiggs_LIB  :=
LLNMSSMEFTHiggs_OBJ  :=
LLNMSSMEFTHiggs_SRC  := \
		$(DIR)/NMSSMEFTHiggs_librarylink.cpp

LLNMSSMEFTHiggs_MMA  := \
		$(DIR)/NMSSMEFTHiggs_librarylink.m \
		$(DIR)/run_NMSSMEFTHiggs.m

LIBNMSSMEFTHiggs_HDR := \
		$(DIR)/NMSSMEFTHiggs_a_muon.hpp \
		$(DIR)/NMSSMEFTHiggs_convergence_tester.hpp \
		$(DIR)/NMSSMEFTHiggs_edm.hpp \
		$(DIR)/NMSSMEFTHiggs_FFV_form_factors.hpp \
		$(DIR)/NMSSMEFTHiggs_l_to_lgamma.hpp \
		$(DIR)/NMSSMEFTHiggs_effective_couplings.hpp \
		$(DIR)/NMSSMEFTHiggs_ewsb_solver.hpp \
		$(DIR)/NMSSMEFTHiggs_ewsb_solver_interface.hpp \
		$(DIR)/NMSSMEFTHiggs_high_scale_constraint.hpp \
		$(DIR)/NMSSMEFTHiggs_info.hpp \
		$(DIR)/NMSSMEFTHiggs_initial_guesser.hpp \
		$(DIR)/NMSSMEFTHiggs_input_parameters.hpp \
		$(DIR)/NMSSMEFTHiggs_low_scale_constraint.hpp \
		$(DIR)/NMSSMEFTHiggs_mass_eigenstates.hpp \
		$(DIR)/NMSSMEFTHiggs_model.hpp \
		$(DIR)/NMSSMEFTHiggs_model_slha.hpp \
		$(DIR)/NMSSMEFTHiggs_observables.hpp \
		$(DIR)/NMSSMEFTHiggs_physical.hpp \
		$(DIR)/NMSSMEFTHiggs_slha_io.hpp \
		$(DIR)/NMSSMEFTHiggs_spectrum_generator.hpp \
		$(DIR)/NMSSMEFTHiggs_spectrum_generator_interface.hpp \
		$(DIR)/NMSSMEFTHiggs_soft_parameters.hpp \
		$(DIR)/NMSSMEFTHiggs_susy_parameters.hpp \
		$(DIR)/NMSSMEFTHiggs_susy_scale_constraint.hpp \
		$(DIR)/NMSSMEFTHiggs_utilities.hpp \
		$(DIR)/NMSSMEFTHiggs_weinberg_angle.hpp

LIBNMSSMEFTHiggs_CXXQFT_HDR := \
		$(DIR)/cxx_qft/NMSSMEFTHiggs_qft.hpp \
		$(DIR)/cxx_qft/NMSSMEFTHiggs_fields.hpp \
		$(DIR)/cxx_qft/NMSSMEFTHiggs_vertices.hpp \
		$(DIR)/cxx_qft/NMSSMEFTHiggs_context_base.hpp \
		$(DIR)/cxx_qft/NMSSMEFTHiggs_npointfunctions.hpp

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
-include $(NMSSMEFTHiggs_SUSY_BETAS_MK)
-include $(NMSSMEFTHiggs_SOFT_BETAS_MK)
-include $(NMSSMEFTHiggs_CXX_QFT_VERTICES_MK)
-include $(NMSSMEFTHiggs_FlexibleEFTHiggs_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(NMSSMEFTHiggs_SUSY_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(NMSSMEFTHiggs_SOFT_BETAS_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(NMSSMEFTHiggs_CXX_QFT_VERTICES_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(NMSSMEFTHiggs_FlexibleEFTHiggs_MK): run-metacode-$(MODNAME)
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
LIBNMSSMEFTHiggs_SRC := $(sort $(LIBNMSSMEFTHiggs_SRC))
EXENMSSMEFTHiggs_SRC := $(sort $(EXENMSSMEFTHiggs_SRC))

LIBNMSSMEFTHiggs_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBNMSSMEFTHiggs_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBNMSSMEFTHiggs_SRC)))

EXENMSSMEFTHiggs_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXENMSSMEFTHiggs_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXENMSSMEFTHiggs_SRC)))

EXENMSSMEFTHiggs_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXENMSSMEFTHiggs_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXENMSSMEFTHiggs_SRC)))

LIBNMSSMEFTHiggs_DEP := \
		$(LIBNMSSMEFTHiggs_OBJ:.o=.d)

EXENMSSMEFTHiggs_DEP := \
		$(EXENMSSMEFTHiggs_OBJ:.o=.d)

LLNMSSMEFTHiggs_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLNMSSMEFTHiggs_SRC)))

LLNMSSMEFTHiggs_OBJ  := $(LLNMSSMEFTHiggs_SRC:.cpp=.o)
LLNMSSMEFTHiggs_LIB  := $(LLNMSSMEFTHiggs_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBNMSSMEFTHiggs     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_NMSSMEFTHiggs := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_NMSSMEFTHiggs := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL_DIR) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBNMSSMEFTHiggs) $(EXENMSSMEFTHiggs_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(NMSSMEFTHiggs_INSTALL_DIR)
		$(Q)install -d $(NMSSMEFTHiggs_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBNMSSMEFTHiggs_SRC) $(NMSSMEFTHiggs_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBNMSSMEFTHiggs_HDR) $(NMSSMEFTHiggs_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBNMSSMEFTHiggs_CXXQFT_HDR) $(NMSSMEFTHiggs_INSTALL_CXXQFT_DIR)
		$(Q)install -m u=rw,g=r,o=r $(EXENMSSMEFTHiggs_SRC) $(NMSSMEFTHiggs_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLNMSSMEFTHiggs_SRC) $(NMSSMEFTHiggs_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LLNMSSMEFTHiggs_MMA) $(NMSSMEFTHiggs_INSTALL_DIR)
		$(Q)$(INSTALL_STRIPPED) $(NMSSMEFTHiggs_MK) $(NMSSMEFTHiggs_INSTALL_DIR) -m u=rw,g=r,o=r
		$(Q)install -m u=rw,g=r,o=r $(NMSSMEFTHiggs_INCLUDE_MK) $(NMSSMEFTHiggs_INSTALL_DIR)
ifneq ($(NMSSMEFTHiggs_SLHA_INPUT),)
		$(Q)install -m u=rw,g=r,o=r $(NMSSMEFTHiggs_SLHA_INPUT) $(NMSSMEFTHiggs_INSTALL_DIR)
endif
		$(Q)install -m u=rw,g=r,o=r $(NMSSMEFTHiggs_REFERENCES) $(NMSSMEFTHiggs_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(NMSSMEFTHiggs_GNUPLOT) $(NMSSMEFTHiggs_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		$(Q)-rm -f $(LIBNMSSMEFTHiggs_DEP)
		$(Q)-rm -f $(EXENMSSMEFTHiggs_DEP)
		$(Q)-rm -f $(LLNMSSMEFTHiggs_DEP)

clean-$(MODNAME)-lib:
		$(Q)-rm -f $(LIBNMSSMEFTHiggs)
		$(Q)-rm -f $(LLNMSSMEFTHiggs_LIB)

clean-$(MODNAME)-obj:
		$(Q)-rm -f $(LIBNMSSMEFTHiggs_OBJ)
		$(Q)-rm -f $(EXENMSSMEFTHiggs_OBJ)
		$(Q)-rm -f $(LLNMSSMEFTHiggs_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		$(Q)-rm -f $(LIBNMSSMEFTHiggs_SRC)
		$(Q)-rm -f $(LIBNMSSMEFTHiggs_HDR)
		$(Q)-rm -f $(LIBNMSSMEFTHiggs_CXXQFT_HDR)
		$(Q)-rm -f $(EXENMSSMEFTHiggs_SRC)
		$(Q)-rm -f $(LLNMSSMEFTHiggs_SRC)
		$(Q)-rm -f $(LLNMSSMEFTHiggs_MMA)
		$(Q)-rm -f $(METACODE_STAMP_NMSSMEFTHiggs)
		$(Q)-rm -f $(NMSSMEFTHiggs_INCLUDE_MK)
		$(Q)-rm -f $(NMSSMEFTHiggs_SLHA_INPUT)
		$(Q)-rm -f $(NMSSMEFTHiggs_REFERENCES)
		$(Q)-rm -f $(NMSSMEFTHiggs_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		$(Q)-rm -f $(EXENMSSMEFTHiggs_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		$(Q)tar -czf $(NMSSMEFTHiggs_TARBALL) \
		$(LIBNMSSMEFTHiggs_SRC) $(LIBNMSSMEFTHiggs_HDR) $(LIBNMSSMEFTHiggs_CXXQFT_HDR) \
		$(EXENMSSMEFTHiggs_SRC) \
		$(LLNMSSMEFTHiggs_SRC) $(LLNMSSMEFTHiggs_MMA) \
		$(NMSSMEFTHiggs_MK) $(NMSSMEFTHiggs_INCLUDE_MK) \
		$(NMSSMEFTHiggs_SLHA_INPUT) $(NMSSMEFTHiggs_REFERENCES) \
		$(NMSSMEFTHiggs_GNUPLOT)

$(LIBNMSSMEFTHiggs_SRC) $(LIBNMSSMEFTHiggs_HDR) $(LIBNMSSMEFTHiggs_CXXQFT_HDR) $(EXENMSSMEFTHiggs_SRC) $(LLNMSSMEFTHiggs_SRC) $(LLNMSSMEFTHiggs_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_NMSSMEFTHiggs)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_NMSSMEFTHiggs): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_NMSSMEFTHiggs)
		@$(MSG)
		$(Q)"$(MATH)" -run "Get[\"$<\"]; Quit[]" || (echo "Error: The code generation failed!"; exit 1)
		@touch "$(METACODE_STAMP_NMSSMEFTHiggs)"
		@echo "Note: to regenerate NMSSMEFTHiggs source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_NMSSMEFTHiggs)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_NMSSMEFTHiggs):
		@true
endif

$(LIBNMSSMEFTHiggs_DEP) $(EXENMSSMEFTHiggs_DEP) $(LLNMSSMEFTHiggs_DEP) $(LIBNMSSMEFTHiggs_OBJ) $(EXENMSSMEFTHiggs_OBJ) $(LLNMSSMEFTHiggs_OBJ) $(LLNMSSMEFTHiggs_LIB): \
	CPPFLAGS += $(MODNMSSMEFTHiggs_SUBMOD_INC) $(MODNMSSMEFTHiggs_INC) $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS) $(HIMALAYAFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBNMSSMEFTHiggs_DEP) $(EXENMSSMEFTHiggs_DEP) $(LLNMSSMEFTHiggs_DEP) $(LIBNMSSMEFTHiggs_OBJ) $(EXENMSSMEFTHiggs_OBJ) $(LLNMSSMEFTHiggs_OBJ) $(LLNMSSMEFTHiggs_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLNMSSMEFTHiggs_OBJ) $(LLNMSSMEFTHiggs_LIB): \
	CPPFLAGS += $(LLFLAGS)

$(LIBNMSSMEFTHiggs): $(LIBNMSSMEFTHiggs_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBNMSSMEFTHiggs) $(MODNMSSMEFTHiggs_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		@$(MSG)
		$(Q)$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLNMSSMEFTHiggs_LIB): $(LLNMSSMEFTHiggs_OBJ) $(LIBNMSSMEFTHiggs) $(MODNMSSMEFTHiggs_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		@$(MSG)
		$(Q)$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$(ADDONLIBS) $^ $(LIBGM2Calc)) $(filter -%,$(LOOPFUNCLIBS)) $(HIMALAYALIBS) $(GSLLIBS) $(BOOSTTHREADLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS) $(LLLIBS)

ALLDEP += $(LIBNMSSMEFTHiggs_DEP) $(EXENMSSMEFTHiggs_DEP)
ALLSRC += $(LIBNMSSMEFTHiggs_SRC) $(EXENMSSMEFTHiggs_SRC)
ALLLIB += $(LIBNMSSMEFTHiggs)
ALLEXE += $(EXENMSSMEFTHiggs_EXE)
ALLMODDEP += $(MODNMSSMEFTHiggs_DEP)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLNMSSMEFTHiggs_DEP)
ALLSRC += $(LLNMSSMEFTHiggs_SRC)
ALLLL  += $(LLNMSSMEFTHiggs_LIB)
endif
