DIR          := models/MSSMtower
MODNAME      := MSSMtower
SARAH_MODEL  := MSSM
WITH_$(MODNAME) := yes

MSSMtower_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

MSSMtower_MK     := \
		$(DIR)/module.mk

MSSMtower_TWO_SCALE_SUSY_MK := \
		$(DIR)/two_scale_susy.mk

MSSMtower_TWO_SCALE_SOFT_MK := \
		$(DIR)/two_scale_soft.mk

MSSMtower_TWO_SCALE_MK := \
		$(MSSMtower_TWO_SCALE_SUSY_MK) \
		$(MSSMtower_TWO_SCALE_SOFT_MK)

MSSMtower_SLHA_INPUT := \
		$(DIR)/LesHouches.in.MSSMtower_generated \
		$(DIR)/LesHouches.in.MSSMtower

MSSMtower_GNUPLOT := \
		$(DIR)/MSSMtower_plot_rgflow.gnuplot \
		$(DIR)/MSSMtower_plot_spectrum.gnuplot

MSSMtower_TARBALL := \
		$(MODNAME).tar.gz

LIBMSSMtower_SRC :=
EXEMSSMtower_SRC :=
LLMSSMtower_LIB  :=
LLMSSMtower_OBJ  :=
LLMSSMtower_SRC  :=
LLMSSMtower_MMA  :=

LIBMSSMtower_HDR :=

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBMSSMtower_SRC += \
		$(DIR)/MSSMtower_effective_couplings.cpp \
		$(DIR)/MSSMtower_mass_eigenstates.cpp \
		$(DIR)/MSSMtower_info.cpp \
		$(DIR)/MSSMtower_input_parameters.cpp \
		$(DIR)/MSSMtower_observables.cpp \
		$(DIR)/MSSMtower_slha_io.cpp \
		$(DIR)/MSSMtower_physical.cpp \
		$(DIR)/MSSMtower_utilities.cpp \
		$(DIR)/MSSMtower_standard_model_matching.cpp \
		$(DIR)/MSSMtower_standard_model_two_scale_matching.cpp \
		$(DIR)/MSSMtower_two_scale_convergence_tester.cpp \
		$(DIR)/MSSMtower_two_scale_high_scale_constraint.cpp \
		$(DIR)/MSSMtower_two_scale_initial_guesser.cpp \
		$(DIR)/MSSMtower_two_scale_low_scale_constraint.cpp \
		$(DIR)/MSSMtower_two_scale_model.cpp \
		$(DIR)/MSSMtower_two_scale_model_slha.cpp \
		$(DIR)/MSSMtower_two_scale_susy_parameters.cpp \
		$(DIR)/MSSMtower_two_scale_soft_parameters.cpp \
		$(DIR)/MSSMtower_two_scale_susy_scale_constraint.cpp
EXEMSSMtower_SRC += \
		$(DIR)/run_MSSMtower.cpp \
		$(DIR)/run_cmd_line_MSSMtower.cpp \
		$(DIR)/scan_MSSMtower.cpp
LIBMSSMtower_HDR += \
		$(DIR)/MSSMtower_convergence_tester.hpp \
		$(DIR)/MSSMtower_effective_couplings.hpp \
		$(DIR)/MSSMtower_high_scale_constraint.hpp \
		$(DIR)/MSSMtower_mass_eigenstates.hpp \
		$(DIR)/MSSMtower_info.hpp \
		$(DIR)/MSSMtower_initial_guesser.hpp \
		$(DIR)/MSSMtower_input_parameters.hpp \
		$(DIR)/MSSMtower_low_scale_constraint.hpp \
		$(DIR)/MSSMtower_model.hpp \
		$(DIR)/MSSMtower_model_slha.hpp \
		$(DIR)/MSSMtower_observables.hpp \
		$(DIR)/MSSMtower_physical.hpp \
		$(DIR)/MSSMtower_slha_io.hpp \
		$(DIR)/MSSMtower_spectrum_generator_interface.hpp \
		$(DIR)/MSSMtower_spectrum_generator.hpp \
		$(DIR)/MSSMtower_standard_model_matching.hpp \
		$(DIR)/MSSMtower_standard_model_two_scale_matching.hpp \
		$(DIR)/MSSMtower_susy_scale_constraint.hpp \
		$(DIR)/MSSMtower_utilities.hpp \
		$(DIR)/MSSMtower_two_scale_convergence_tester.hpp \
		$(DIR)/MSSMtower_two_scale_high_scale_constraint.hpp \
		$(DIR)/MSSMtower_two_scale_initial_guesser.hpp \
		$(DIR)/MSSMtower_two_scale_low_scale_constraint.hpp \
		$(DIR)/MSSMtower_two_scale_model.hpp \
		$(DIR)/MSSMtower_two_scale_model_slha.hpp \
		$(DIR)/MSSMtower_two_scale_soft_parameters.hpp \
		$(DIR)/MSSMtower_two_scale_susy_parameters.hpp \
		$(DIR)/MSSMtower_two_scale_susy_scale_constraint.hpp
LLMSSMtower_SRC  += \
		$(DIR)/MSSMtower_librarylink.cpp

LLMSSMtower_MMA  += \
		$(DIR)/MSSMtower_librarylink.m \
		$(DIR)/run_MSSMtower.m

ifneq ($(MAKECMDGOALS),showbuild)
ifneq ($(MAKECMDGOALS),tag)
ifneq ($(MAKECMDGOALS),release)
ifneq ($(MAKECMDGOALS),doc)
-include $(MSSMtower_TWO_SCALE_SUSY_MK)
-include $(MSSMtower_TWO_SCALE_SOFT_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(MSSMtower_TWO_SCALE_SUSY_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(MSSMtower_TWO_SCALE_SOFT_MK): run-metacode-$(MODNAME)
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
LIBMSSMtower_SRC := $(sort $(LIBMSSMtower_SRC))
EXEMSSMtower_SRC := $(sort $(EXEMSSMtower_SRC))

LIBMSSMtower_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBMSSMtower_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBMSSMtower_SRC)))

EXEMSSMtower_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEMSSMtower_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEMSSMtower_SRC)))

EXEMSSMtower_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXEMSSMtower_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXEMSSMtower_SRC)))

LIBMSSMtower_DEP := \
		$(LIBMSSMtower_OBJ:.o=.d)

EXEMSSMtower_DEP := \
		$(EXEMSSMtower_OBJ:.o=.d)

LLMSSMtower_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLMSSMtower_SRC)))

LLMSSMtower_OBJ  := $(LLMSSMtower_SRC:.cpp=.o)
LLMSSMtower_LIB  := $(LLMSSMtower_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBMSSMtower     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_MSSMtower := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_MSSMtower := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBMSSMtower) $(EXEMSSMtower_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(MSSMtower_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBMSSMtower_SRC) $(MSSMtower_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBMSSMtower_HDR) $(MSSMtower_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXEMSSMtower_SRC) $(MSSMtower_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLMSSMtower_SRC) $(MSSMtower_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLMSSMtower_MMA) $(MSSMtower_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(MSSMtower_MK) $(MSSMtower_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(MSSMtower_TWO_SCALE_MK) $(MSSMtower_INSTALL_DIR)
ifneq ($(MSSMtower_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(MSSMtower_SLHA_INPUT) $(MSSMtower_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(MSSMtower_GNUPLOT) $(MSSMtower_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBMSSMtower_DEP)
		-rm -f $(EXEMSSMtower_DEP)
		-rm -f $(LLMSSMtower_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBMSSMtower)
		-rm -f $(LLMSSMtower_LIB)

clean-$(MODNAME)-obj:
		-rm -f $(LIBMSSMtower_OBJ)
		-rm -f $(EXEMSSMtower_OBJ)
		-rm -f $(LLMSSMtower_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		-rm -f $(LIBMSSMtower_SRC)
		-rm -f $(LIBMSSMtower_HDR)
		-rm -f $(EXEMSSMtower_SRC)
		-rm -f $(LLMSSMtower_SRC)
		-rm -f $(LLMSSMtower_MMA)
		-rm -f $(METACODE_STAMP_MSSMtower)
		-rm -f $(MSSMtower_TWO_SCALE_MK)
		-rm -f $(MSSMtower_SLHA_INPUT)
		-rm -f $(MSSMtower_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		-rm -f $(EXEMSSMtower_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(MSSMtower_TARBALL) \
		$(LIBMSSMtower_SRC) $(LIBMSSMtower_HDR) \
		$(EXEMSSMtower_SRC) \
		$(LLMSSMtower_SRC) $(LLMSSMtower_MMA) \
		$(MSSMtower_MK) $(MSSMtower_TWO_SCALE_MK) \
		$(MSSMtower_SLHA_INPUT) $(MSSMtower_GNUPLOT)

$(LIBMSSMtower_SRC) $(LIBMSSMtower_HDR) $(EXEMSSMtower_SRC) $(LLMSSMtower_SRC) $(LLMSSMtower_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_MSSMtower)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_MSSMtower): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_MSSMtower)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_MSSMtower)"
		@echo "Note: to regenerate MSSMtower source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_MSSMtower)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_MSSMtower):
		@true
endif

$(LIBMSSMtower_DEP) $(EXEMSSMtower_DEP) $(LLMSSMtower_DEP) $(LIBMSSMtower_OBJ) $(EXEMSSMtower_OBJ) $(LLMSSMtower_OBJ) $(LLMSSMtower_LIB): \
	CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBMSSMtower_DEP) $(EXEMSSMtower_DEP) $(LLMSSMtower_DEP) $(LIBMSSMtower_OBJ) $(EXEMSSMtower_OBJ) $(LLMSSMtower_OBJ) $(LLMSSMtower_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLMSSMtower_OBJ) $(LLMSSMtower_LIB): \
	CPPFLAGS += $(shell $(MATH_INC_PATHS) --math-cmd="$(MATH)" -I --librarylink --mathlink)

$(LIBMSSMtower): $(LIBMSSMtower_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBMSSMtower) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$^ $(ADDONLIBS)) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLMSSMtower_LIB): $(LLMSSMtower_OBJ) $(LIBMSSMtower) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$^) $(ADDONLIBS) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

ALLDEP += $(LIBMSSMtower_DEP) $(EXEMSSMtower_DEP)
ALLSRC += $(LIBMSSMtower_SRC) $(EXEMSSMtower_SRC)
ALLLIB += $(LIBMSSMtower)
ALLEXE += $(EXEMSSMtower_EXE)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLMSSMtower_DEP)
ALLSRC += $(LLMSSMtower_SRC)
ALLLL  += $(LLMSSMtower_LIB)
endif
