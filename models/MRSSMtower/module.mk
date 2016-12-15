DIR          := models/MRSSMtower
MODNAME      := MRSSMtower
SARAH_MODEL  := MRSSM
WITH_$(MODNAME) := yes

MRSSMtower_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

MRSSMtower_MK     := \
		$(DIR)/module.mk

MRSSMtower_TWO_SCALE_SUSY_MK := \
		$(DIR)/two_scale_susy.mk

MRSSMtower_TWO_SCALE_SOFT_MK := \
		$(DIR)/two_scale_soft.mk

MRSSMtower_TWO_SCALE_MK := \
		$(MRSSMtower_TWO_SCALE_SUSY_MK) \
		$(MRSSMtower_TWO_SCALE_SOFT_MK)

MRSSMtower_SLHA_INPUT := \
		$(DIR)/LesHouches.in.MRSSMtower_generated \
		$(DIR)/LesHouches.in.MRSSMtower

MRSSMtower_GNUPLOT := \
		$(DIR)/MRSSMtower_plot_rgflow.gnuplot \
		$(DIR)/MRSSMtower_plot_spectrum.gnuplot

MRSSMtower_TARBALL := \
		$(MODNAME).tar.gz

LIBMRSSMtower_SRC :=
EXEMRSSMtower_SRC :=
LLMRSSMtower_LIB  :=
LLMRSSMtower_OBJ  :=
LLMRSSMtower_SRC  :=
LLMRSSMtower_MMA  :=

LIBMRSSMtower_HDR :=

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBMRSSMtower_SRC += \
		$(DIR)/MRSSMtower_effective_couplings.cpp \
		$(DIR)/MRSSMtower_mass_eigenstates.cpp \
		$(DIR)/MRSSMtower_info.cpp \
		$(DIR)/MRSSMtower_input_parameters.cpp \
		$(DIR)/MRSSMtower_observables.cpp \
		$(DIR)/MRSSMtower_slha_io.cpp \
		$(DIR)/MRSSMtower_physical.cpp \
		$(DIR)/MRSSMtower_utilities.cpp \
		$(DIR)/MRSSMtower_standard_model_matching.cpp \
		$(DIR)/MRSSMtower_standard_model_two_scale_matching.cpp \
		$(DIR)/MRSSMtower_two_scale_convergence_tester.cpp \
		$(DIR)/MRSSMtower_two_scale_high_scale_constraint.cpp \
		$(DIR)/MRSSMtower_two_scale_initial_guesser.cpp \
		$(DIR)/MRSSMtower_two_scale_low_scale_constraint.cpp \
		$(DIR)/MRSSMtower_two_scale_model.cpp \
		$(DIR)/MRSSMtower_two_scale_model_slha.cpp \
		$(DIR)/MRSSMtower_two_scale_susy_parameters.cpp \
		$(DIR)/MRSSMtower_two_scale_soft_parameters.cpp \
		$(DIR)/MRSSMtower_two_scale_susy_scale_constraint.cpp
EXEMRSSMtower_SRC += \
		$(DIR)/run_MRSSMtower.cpp \
		$(DIR)/run_cmd_line_MRSSMtower.cpp \
		$(DIR)/scan_MRSSMtower.cpp
LIBMRSSMtower_HDR += \
		$(DIR)/MRSSMtower_convergence_tester.hpp \
		$(DIR)/MRSSMtower_effective_couplings.hpp \
		$(DIR)/MRSSMtower_high_scale_constraint.hpp \
		$(DIR)/MRSSMtower_mass_eigenstates.hpp \
		$(DIR)/MRSSMtower_info.hpp \
		$(DIR)/MRSSMtower_initial_guesser.hpp \
		$(DIR)/MRSSMtower_input_parameters.hpp \
		$(DIR)/MRSSMtower_low_scale_constraint.hpp \
		$(DIR)/MRSSMtower_model.hpp \
		$(DIR)/MRSSMtower_model_slha.hpp \
		$(DIR)/MRSSMtower_observables.hpp \
		$(DIR)/MRSSMtower_physical.hpp \
		$(DIR)/MRSSMtower_slha_io.hpp \
		$(DIR)/MRSSMtower_spectrum_generator_interface.hpp \
		$(DIR)/MRSSMtower_spectrum_generator.hpp \
		$(DIR)/MRSSMtower_standard_model_matching.hpp \
		$(DIR)/MRSSMtower_standard_model_two_scale_matching.hpp \
		$(DIR)/MRSSMtower_susy_scale_constraint.hpp \
		$(DIR)/MRSSMtower_utilities.hpp \
		$(DIR)/MRSSMtower_two_scale_convergence_tester.hpp \
		$(DIR)/MRSSMtower_two_scale_high_scale_constraint.hpp \
		$(DIR)/MRSSMtower_two_scale_initial_guesser.hpp \
		$(DIR)/MRSSMtower_two_scale_low_scale_constraint.hpp \
		$(DIR)/MRSSMtower_two_scale_model.hpp \
		$(DIR)/MRSSMtower_two_scale_model_slha.hpp \
		$(DIR)/MRSSMtower_two_scale_soft_parameters.hpp \
		$(DIR)/MRSSMtower_two_scale_susy_parameters.hpp \
		$(DIR)/MRSSMtower_two_scale_susy_scale_constraint.hpp
LLMRSSMtower_SRC  += \
		$(DIR)/MRSSMtower_librarylink.cpp

LLMRSSMtower_MMA  += \
		$(DIR)/MRSSMtower_librarylink.m \
		$(DIR)/run_MRSSMtower.m

ifneq ($(MAKECMDGOALS),showbuild)
ifneq ($(MAKECMDGOALS),tag)
ifneq ($(MAKECMDGOALS),release)
ifneq ($(MAKECMDGOALS),doc)
-include $(MRSSMtower_TWO_SCALE_SUSY_MK)
-include $(MRSSMtower_TWO_SCALE_SOFT_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(MRSSMtower_TWO_SCALE_SUSY_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(MRSSMtower_TWO_SCALE_SOFT_MK): run-metacode-$(MODNAME)
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
LIBMRSSMtower_SRC := $(sort $(LIBMRSSMtower_SRC))
EXEMRSSMtower_SRC := $(sort $(EXEMRSSMtower_SRC))

LIBMRSSMtower_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBMRSSMtower_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBMRSSMtower_SRC)))

EXEMRSSMtower_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEMRSSMtower_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEMRSSMtower_SRC)))

EXEMRSSMtower_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXEMRSSMtower_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXEMRSSMtower_SRC)))

LIBMRSSMtower_DEP := \
		$(LIBMRSSMtower_OBJ:.o=.d)

EXEMRSSMtower_DEP := \
		$(EXEMRSSMtower_OBJ:.o=.d)

LLMRSSMtower_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLMRSSMtower_SRC)))

LLMRSSMtower_OBJ  := $(LLMRSSMtower_SRC:.cpp=.o)
LLMRSSMtower_LIB  := $(LLMRSSMtower_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBMRSSMtower     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_MRSSMtower := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_MRSSMtower := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBMRSSMtower) $(EXEMRSSMtower_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(MRSSMtower_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBMRSSMtower_SRC) $(MRSSMtower_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBMRSSMtower_HDR) $(MRSSMtower_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXEMRSSMtower_SRC) $(MRSSMtower_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLMRSSMtower_SRC) $(MRSSMtower_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLMRSSMtower_MMA) $(MRSSMtower_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(MRSSMtower_MK) $(MRSSMtower_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(MRSSMtower_TWO_SCALE_MK) $(MRSSMtower_INSTALL_DIR)
ifneq ($(MRSSMtower_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(MRSSMtower_SLHA_INPUT) $(MRSSMtower_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(MRSSMtower_GNUPLOT) $(MRSSMtower_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBMRSSMtower_DEP)
		-rm -f $(EXEMRSSMtower_DEP)
		-rm -f $(LLMRSSMtower_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBMRSSMtower)
		-rm -f $(LLMRSSMtower_LIB)

clean-$(MODNAME)-obj:
		-rm -f $(LIBMRSSMtower_OBJ)
		-rm -f $(EXEMRSSMtower_OBJ)
		-rm -f $(LLMRSSMtower_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		-rm -f $(LIBMRSSMtower_SRC)
		-rm -f $(LIBMRSSMtower_HDR)
		-rm -f $(EXEMRSSMtower_SRC)
		-rm -f $(LLMRSSMtower_SRC)
		-rm -f $(LLMRSSMtower_MMA)
		-rm -f $(METACODE_STAMP_MRSSMtower)
		-rm -f $(MRSSMtower_TWO_SCALE_MK)
		-rm -f $(MRSSMtower_SLHA_INPUT)
		-rm -f $(MRSSMtower_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		-rm -f $(EXEMRSSMtower_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(MRSSMtower_TARBALL) \
		$(LIBMRSSMtower_SRC) $(LIBMRSSMtower_HDR) \
		$(EXEMRSSMtower_SRC) \
		$(LLMRSSMtower_SRC) $(LLMRSSMtower_MMA) \
		$(MRSSMtower_MK) $(MRSSMtower_TWO_SCALE_MK) \
		$(MRSSMtower_SLHA_INPUT) $(MRSSMtower_GNUPLOT)

$(LIBMRSSMtower_SRC) $(LIBMRSSMtower_HDR) $(EXEMRSSMtower_SRC) $(LLMRSSMtower_SRC) $(LLMRSSMtower_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_MRSSMtower)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_MRSSMtower): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_MRSSMtower)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_MRSSMtower)"
		@echo "Note: to regenerate MRSSMtower source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_MRSSMtower)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_MRSSMtower):
		@true
endif

$(LIBMRSSMtower_DEP) $(EXEMRSSMtower_DEP) $(LLMRSSMtower_DEP) $(LIBMRSSMtower_OBJ) $(EXEMRSSMtower_OBJ) $(LLMRSSMtower_OBJ) $(LLMRSSMtower_LIB): \
	CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBMRSSMtower_DEP) $(EXEMRSSMtower_DEP) $(LLMRSSMtower_DEP) $(LIBMRSSMtower_OBJ) $(EXEMRSSMtower_OBJ) $(LLMRSSMtower_OBJ) $(LLMRSSMtower_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLMRSSMtower_OBJ) $(LLMRSSMtower_LIB): \
	CPPFLAGS += $(shell $(MATH_INC_PATHS) --math-cmd="$(MATH)" -I --librarylink --mathlink)

$(LIBMRSSMtower): $(LIBMRSSMtower_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBMRSSMtower) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$^ $(ADDONLIBS)) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLMRSSMtower_LIB): $(LLMRSSMtower_OBJ) $(LIBMRSSMtower) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$^) $(ADDONLIBS) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

ALLDEP += $(LIBMRSSMtower_DEP) $(EXEMRSSMtower_DEP)
ALLSRC += $(LIBMRSSMtower_SRC) $(EXEMRSSMtower_SRC)
ALLLIB += $(LIBMRSSMtower)
ALLEXE += $(EXEMRSSMtower_EXE)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLMRSSMtower_DEP)
ALLSRC += $(LLMRSSMtower_SRC)
ALLLL  += $(LLMRSSMtower_LIB)
endif
