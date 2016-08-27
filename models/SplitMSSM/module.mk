DIR          := models/SplitMSSM
MODNAME      := SplitMSSM
SARAH_MODEL  := SplitMSSM
WITH_$(MODNAME) := yes

SplitMSSM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

SplitMSSM_MK     := \
		$(DIR)/module.mk

SplitMSSM_TWO_SCALE_SUSY_MK := \
		$(DIR)/two_scale_susy.mk

SplitMSSM_TWO_SCALE_SOFT_MK := \
		$(DIR)/two_scale_soft.mk

SplitMSSM_TWO_SCALE_MK := \
		$(SplitMSSM_TWO_SCALE_SUSY_MK) \
		$(SplitMSSM_TWO_SCALE_SOFT_MK)

SplitMSSM_SLHA_INPUT := \
		$(DIR)/LesHouches.in.SplitMSSM_generated \
		$(DIR)/LesHouches.in.SplitMSSM

SplitMSSM_GNUPLOT := \
		$(DIR)/SplitMSSM_plot_rgflow.gnuplot \
		$(DIR)/SplitMSSM_plot_spectrum.gnuplot

SplitMSSM_TARBALL := \
		$(MODNAME).tar.gz

LIBSplitMSSM_SRC :=
EXESplitMSSM_SRC :=

LIBSplitMSSM_HDR :=

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBSplitMSSM_SRC += \
		$(DIR)/SplitMSSM_effective_couplings.cpp \
		$(DIR)/SplitMSSM_mass_eigenstates.cpp \
		$(DIR)/SplitMSSM_info.cpp \
		$(DIR)/SplitMSSM_input_parameters.cpp \
		$(DIR)/SplitMSSM_observables.cpp \
		$(DIR)/SplitMSSM_slha_io.cpp \
		$(DIR)/SplitMSSM_physical.cpp \
		$(DIR)/SplitMSSM_utilities.cpp \
		$(DIR)/SplitMSSM_two_scale_convergence_tester.cpp \
		$(DIR)/SplitMSSM_two_scale_high_scale_constraint.cpp \
		$(DIR)/SplitMSSM_two_scale_initial_guesser.cpp \
		$(DIR)/SplitMSSM_two_scale_low_scale_constraint.cpp \
		$(DIR)/SplitMSSM_two_scale_model.cpp \
		$(DIR)/SplitMSSM_two_scale_model_slha.cpp \
		$(DIR)/SplitMSSM_two_scale_susy_parameters.cpp \
		$(DIR)/SplitMSSM_two_scale_soft_parameters.cpp \
		$(DIR)/SplitMSSM_two_scale_susy_scale_constraint.cpp
EXESplitMSSM_SRC += \
		$(DIR)/run_SplitMSSM.cpp \
		$(DIR)/run_cmd_line_SplitMSSM.cpp \
		$(DIR)/scan_SplitMSSM.cpp
LIBSplitMSSM_HDR += \
		$(DIR)/SplitMSSM_convergence_tester.hpp \
		$(DIR)/SplitMSSM_effective_couplings.hpp \
		$(DIR)/SplitMSSM_high_scale_constraint.hpp \
		$(DIR)/SplitMSSM_mass_eigenstates.hpp \
		$(DIR)/SplitMSSM_info.hpp \
		$(DIR)/SplitMSSM_initial_guesser.hpp \
		$(DIR)/SplitMSSM_input_parameters.hpp \
		$(DIR)/SplitMSSM_low_scale_constraint.hpp \
		$(DIR)/SplitMSSM_model.hpp \
		$(DIR)/SplitMSSM_model_slha.hpp \
		$(DIR)/SplitMSSM_observables.hpp \
		$(DIR)/SplitMSSM_physical.hpp \
		$(DIR)/SplitMSSM_slha_io.hpp \
		$(DIR)/SplitMSSM_spectrum_generator_interface.hpp \
		$(DIR)/SplitMSSM_spectrum_generator.hpp \
		$(DIR)/SplitMSSM_susy_scale_constraint.hpp \
		$(DIR)/SplitMSSM_utilities.hpp \
		$(DIR)/SplitMSSM_two_scale_convergence_tester.hpp \
		$(DIR)/SplitMSSM_two_scale_high_scale_constraint.hpp \
		$(DIR)/SplitMSSM_two_scale_initial_guesser.hpp \
		$(DIR)/SplitMSSM_two_scale_low_scale_constraint.hpp \
		$(DIR)/SplitMSSM_two_scale_model.hpp \
		$(DIR)/SplitMSSM_two_scale_model_slha.hpp \
		$(DIR)/SplitMSSM_two_scale_soft_parameters.hpp \
		$(DIR)/SplitMSSM_two_scale_susy_parameters.hpp \
		$(DIR)/SplitMSSM_two_scale_susy_scale_constraint.hpp

ifneq ($(MAKECMDGOALS),showbuild)
ifneq ($(MAKECMDGOALS),tag)
ifneq ($(MAKECMDGOALS),release)
ifneq ($(MAKECMDGOALS),doc)
-include $(SplitMSSM_TWO_SCALE_SUSY_MK)
-include $(SplitMSSM_TWO_SCALE_SOFT_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(SplitMSSM_TWO_SCALE_SUSY_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(SplitMSSM_TWO_SCALE_SOFT_MK): run-metacode-$(MODNAME)
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

LIBSplitMSSM     := $(DIR)/lib$(MODNAME)$(LIBEXT)

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
		$(INSTALL_STRIPPED) $(SplitMSSM_MK) $(SplitMSSM_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(SplitMSSM_TWO_SCALE_MK) $(SplitMSSM_INSTALL_DIR)
ifneq ($(SplitMSSM_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(SplitMSSM_SLHA_INPUT) $(SplitMSSM_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(SplitMSSM_GNUPLOT) $(SplitMSSM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBSplitMSSM_DEP)
		-rm -f $(EXESplitMSSM_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBSplitMSSM)

clean-$(MODNAME)-obj:
		-rm -f $(LIBSplitMSSM_OBJ)
		-rm -f $(EXESplitMSSM_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		-rm -f $(LIBSplitMSSM_SRC)
		-rm -f $(LIBSplitMSSM_HDR)
		-rm -f $(EXESplitMSSM_SRC)
		-rm -f $(METACODE_STAMP_SplitMSSM)
		-rm -f $(SplitMSSM_TWO_SCALE_MK)
		-rm -f $(SplitMSSM_SLHA_INPUT)
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
		$(SplitMSSM_MK) $(SplitMSSM_TWO_SCALE_MK) \
		$(SplitMSSM_SLHA_INPUT) $(SplitMSSM_GNUPLOT)

$(LIBSplitMSSM_SRC) $(LIBSplitMSSM_HDR) $(EXESplitMSSM_SRC) \
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

$(LIBSplitMSSM_DEP) $(EXESplitMSSM_DEP) $(LIBSplitMSSM_OBJ) $(EXESplitMSSM_OBJ): CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBSplitMSSM_DEP) $(EXESplitMSSM_DEP) $(LIBSplitMSSM_OBJ) $(EXESplitMSSM_OBJ): CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LIBSplitMSSM): $(LIBSplitMSSM_OBJ)
		$(MAKELIB) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBSplitMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$^ $(ADDONLIBS)) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(LDLIBS)

ALLDEP += $(LIBSplitMSSM_DEP) $(EXESplitMSSM_DEP)
ALLSRC += $(LIBSplitMSSM_SRC) $(EXESplitMSSM_SRC)
ALLLIB += $(LIBSplitMSSM)
ALLEXE += $(EXESplitMSSM_EXE)
