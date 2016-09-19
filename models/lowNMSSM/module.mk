DIR          := models/lowNMSSM
MODNAME      := lowNMSSM
SARAH_MODEL  := NMSSM
WITH_$(MODNAME) := yes

lowNMSSM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

lowNMSSM_MK     := \
		$(DIR)/module.mk

lowNMSSM_TWO_SCALE_SUSY_MK := \
		$(DIR)/two_scale_susy.mk

lowNMSSM_TWO_SCALE_SOFT_MK := \
		$(DIR)/two_scale_soft.mk

lowNMSSM_TWO_SCALE_MK := \
		$(lowNMSSM_TWO_SCALE_SUSY_MK) \
		$(lowNMSSM_TWO_SCALE_SOFT_MK)

lowNMSSM_SLHA_INPUT := \
		$(DIR)/LesHouches.in.lowNMSSM_generated \
		$(DIR)/LesHouches.in.TP4 \
		$(DIR)/LesHouches.in.TP3 \
		$(DIR)/LesHouches.in.TP5 \
		$(DIR)/LesHouches.in.TP2 \
		$(DIR)/LesHouches.in.TP6 \
		$(DIR)/LesHouches.in.lowNMSSM \
		$(DIR)/LesHouches.in.TP1

lowNMSSM_GNUPLOT := \
		$(DIR)/lowNMSSM_plot_rgflow.gnuplot \
		$(DIR)/lowNMSSM_plot_spectrum.gnuplot

lowNMSSM_TARBALL := \
		$(MODNAME).tar.gz

LIBlowNMSSM_SRC :=
EXElowNMSSM_SRC :=
LLlowNMSSM_LIB  :=
LLlowNMSSM_OBJ  :=
LLlowNMSSM_SRC  :=
LLlowNMSSM_MMA  :=

LIBlowNMSSM_HDR :=

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBlowNMSSM_SRC += \
		$(DIR)/lowNMSSM_effective_couplings.cpp \
		$(DIR)/lowNMSSM_mass_eigenstates.cpp \
		$(DIR)/lowNMSSM_info.cpp \
		$(DIR)/lowNMSSM_input_parameters.cpp \
		$(DIR)/lowNMSSM_observables.cpp \
		$(DIR)/lowNMSSM_slha_io.cpp \
		$(DIR)/lowNMSSM_physical.cpp \
		$(DIR)/lowNMSSM_utilities.cpp \
		$(DIR)/lowNMSSM_standard_model_matching.cpp \
		$(DIR)/lowNMSSM_standard_model_two_scale_matching.cpp \
		$(DIR)/lowNMSSM_two_scale_convergence_tester.cpp \
		$(DIR)/lowNMSSM_two_scale_high_scale_constraint.cpp \
		$(DIR)/lowNMSSM_two_scale_initial_guesser.cpp \
		$(DIR)/lowNMSSM_two_scale_low_scale_constraint.cpp \
		$(DIR)/lowNMSSM_two_scale_model.cpp \
		$(DIR)/lowNMSSM_two_scale_model_slha.cpp \
		$(DIR)/lowNMSSM_two_scale_susy_parameters.cpp \
		$(DIR)/lowNMSSM_two_scale_soft_parameters.cpp \
		$(DIR)/lowNMSSM_two_scale_susy_scale_constraint.cpp
EXElowNMSSM_SRC += \
		$(DIR)/run_lowNMSSM.cpp \
		$(DIR)/run_cmd_line_lowNMSSM.cpp \
		$(DIR)/scan_lowNMSSM.cpp
LIBlowNMSSM_HDR += \
		$(DIR)/lowNMSSM_convergence_tester.hpp \
		$(DIR)/lowNMSSM_effective_couplings.hpp \
		$(DIR)/lowNMSSM_high_scale_constraint.hpp \
		$(DIR)/lowNMSSM_mass_eigenstates.hpp \
		$(DIR)/lowNMSSM_info.hpp \
		$(DIR)/lowNMSSM_initial_guesser.hpp \
		$(DIR)/lowNMSSM_input_parameters.hpp \
		$(DIR)/lowNMSSM_low_scale_constraint.hpp \
		$(DIR)/lowNMSSM_model.hpp \
		$(DIR)/lowNMSSM_model_slha.hpp \
		$(DIR)/lowNMSSM_observables.hpp \
		$(DIR)/lowNMSSM_physical.hpp \
		$(DIR)/lowNMSSM_slha_io.hpp \
		$(DIR)/lowNMSSM_spectrum_generator_interface.hpp \
		$(DIR)/lowNMSSM_spectrum_generator.hpp \
		$(DIR)/lowNMSSM_standard_model_matching.hpp \
		$(DIR)/lowNMSSM_standard_model_two_scale_matching.hpp \
		$(DIR)/lowNMSSM_susy_scale_constraint.hpp \
		$(DIR)/lowNMSSM_utilities.hpp \
		$(DIR)/lowNMSSM_two_scale_convergence_tester.hpp \
		$(DIR)/lowNMSSM_two_scale_high_scale_constraint.hpp \
		$(DIR)/lowNMSSM_two_scale_initial_guesser.hpp \
		$(DIR)/lowNMSSM_two_scale_low_scale_constraint.hpp \
		$(DIR)/lowNMSSM_two_scale_model.hpp \
		$(DIR)/lowNMSSM_two_scale_model_slha.hpp \
		$(DIR)/lowNMSSM_two_scale_soft_parameters.hpp \
		$(DIR)/lowNMSSM_two_scale_susy_parameters.hpp \
		$(DIR)/lowNMSSM_two_scale_susy_scale_constraint.hpp
LLlowNMSSM_SRC  += \
		$(DIR)/lowNMSSM_librarylink.cpp

LLlowNMSSM_MMA  += \
		$(DIR)/lowNMSSM_librarylink.m \
		$(DIR)/run_lowNMSSM.m

ifneq ($(MAKECMDGOALS),showbuild)
ifneq ($(MAKECMDGOALS),tag)
ifneq ($(MAKECMDGOALS),release)
ifneq ($(MAKECMDGOALS),doc)
-include $(lowNMSSM_TWO_SCALE_SUSY_MK)
-include $(lowNMSSM_TWO_SCALE_SOFT_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(lowNMSSM_TWO_SCALE_SUSY_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(lowNMSSM_TWO_SCALE_SOFT_MK): run-metacode-$(MODNAME)
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

endif

# remove duplicates in case all algorithms are used
LIBlowNMSSM_SRC := $(sort $(LIBlowNMSSM_SRC))
EXElowNMSSM_SRC := $(sort $(EXElowNMSSM_SRC))

LIBlowNMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBlowNMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBlowNMSSM_SRC)))

EXElowNMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXElowNMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXElowNMSSM_SRC)))

EXElowNMSSM_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXElowNMSSM_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXElowNMSSM_SRC)))

LIBlowNMSSM_DEP := \
		$(LIBlowNMSSM_OBJ:.o=.d)

EXElowNMSSM_DEP := \
		$(EXElowNMSSM_OBJ:.o=.d)

LLlowNMSSM_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLlowNMSSM_SRC)))

LLlowNMSSM_OBJ  := $(LLlowNMSSM_SRC:.cpp=.o)
LLlowNMSSM_LIB  := $(LLlowNMSSM_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBlowNMSSM     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_lowNMSSM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_lowNMSSM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBlowNMSSM) $(EXElowNMSSM_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(lowNMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBlowNMSSM_SRC) $(lowNMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBlowNMSSM_HDR) $(lowNMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXElowNMSSM_SRC) $(lowNMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLlowNMSSM_SRC) $(lowNMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLlowNMSSM_MMA) $(lowNMSSM_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(lowNMSSM_MK) $(lowNMSSM_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(lowNMSSM_TWO_SCALE_MK) $(lowNMSSM_INSTALL_DIR)
ifneq ($(lowNMSSM_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(lowNMSSM_SLHA_INPUT) $(lowNMSSM_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(lowNMSSM_GNUPLOT) $(lowNMSSM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBlowNMSSM_DEP)
		-rm -f $(EXElowNMSSM_DEP)
		-rm -f $(LLlowNMSSM_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBlowNMSSM)
		-rm -f $(LLlowNMSSM_LIB)

clean-$(MODNAME)-obj:
		-rm -f $(LIBlowNMSSM_OBJ)
		-rm -f $(EXElowNMSSM_OBJ)
		-rm -f $(LLlowNMSSM_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		-rm -f $(LIBlowNMSSM_SRC)
		-rm -f $(LIBlowNMSSM_HDR)
		-rm -f $(EXElowNMSSM_SRC)
		-rm -f $(LLlowNMSSM_SRC)
		-rm -f $(LLlowNMSSM_MMA)
		-rm -f $(METACODE_STAMP_lowNMSSM)
		-rm -f $(lowNMSSM_TWO_SCALE_MK)
		-rm -f $(lowNMSSM_SLHA_INPUT)
		-rm -f $(lowNMSSM_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		-rm -f $(EXElowNMSSM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(lowNMSSM_TARBALL) \
		$(LIBlowNMSSM_SRC) $(LIBlowNMSSM_HDR) \
		$(EXElowNMSSM_SRC) \
		$(lowNMSSM_MK) $(lowNMSSM_TWO_SCALE_MK) \
		$(lowNMSSM_SLHA_INPUT) $(lowNMSSM_GNUPLOT)

$(LIBlowNMSSM_SRC) $(LIBlowNMSSM_HDR) $(EXElowNMSSM_SRC) $(LLlowNMSSM_SRC) $(LLlowNMSSM_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_lowNMSSM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_lowNMSSM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_lowNMSSM)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_lowNMSSM)"
		@echo "Note: to regenerate lowNMSSM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_lowNMSSM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_lowNMSSM):
		@true
endif

$(LIBlowNMSSM_DEP) $(EXElowNMSSM_DEP) $(LLlowNMSSM_DEP) $(LIBlowNMSSM_OBJ) $(EXElowNMSSM_OBJ) $(LLlowNMSSM_OBJ) $(LLlowNMSSM_LIB): \
	CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBlowNMSSM_DEP) $(EXElowNMSSM_DEP) $(LLlowNMSSM_DEP) $(LIBlowNMSSM_OBJ) $(EXElowNMSSM_OBJ) $(LLlowNMSSM_OBJ) $(LLlowNMSSM_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLlowNMSSM_OBJ) $(LLlowNMSSM_LIB): \
	CPPFLAGS += $(shell $(MATH_INC_PATHS) --math-cmd="$(MATH)" -I --librarylink --mathlink)

$(LIBlowNMSSM): $(LIBlowNMSSM_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBlowNMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$^ $(ADDONLIBS)) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLlowNMSSM_LIB): $(LLlowNMSSM_OBJ) $(LIBlowNMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$^) $(ADDONLIBS) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

ALLDEP += $(LIBlowNMSSM_DEP) $(EXElowNMSSM_DEP)
ALLSRC += $(LIBlowNMSSM_SRC) $(EXElowNMSSM_SRC)
ALLLIB += $(LIBlowNMSSM)
ALLEXE += $(EXElowNMSSM_EXE)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLlowNMSSM_DEP)
ALLSRC += $(LLlowNMSSM_SRC)
ALLLL  += $(LLlowNMSSM_LIB)
endif
