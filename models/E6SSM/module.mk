DIR          := models/E6SSM
MODNAME      := E6SSM
SARAH_MODEL  := E6SSM

E6SSM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

E6SSM_MK     := \
		$(DIR)/module.mk

E6SSM_TWO_SCALE_SUSY_MK := \
		$(DIR)/two_scale_susy.mk

E6SSM_TWO_SCALE_SOFT_MK := \
		$(DIR)/two_scale_soft.mk

E6SSM_TWO_SCALE_MK := \
		$(E6SSM_TWO_SCALE_SUSY_MK) \
		$(E6SSM_TWO_SCALE_SOFT_MK)

E6SSM_SLHA_INPUT := \
		$(DIR)/LesHouches.in.E6SSM \
		$(DIR)/LesHouches.in.E6SSM~

E6SSM_GNUPLOT := \
		$(DIR)/E6SSM_plot_rgflow.gnuplot \
		$(DIR)/E6SSM_plot_spectrum.gnuplot

E6SSM_TARBALL := \
		$(MODNAME).tar.gz

LIBE6SSM_SRC :=
EXEE6SSM_SRC :=

LIBE6SSM_HDR :=

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBE6SSM_SRC += \
		$(DIR)/E6SSM_mass_eigenstates.cpp \
		$(DIR)/E6SSM_info.cpp \
		$(DIR)/E6SSM_input_parameters.cpp \
		$(DIR)/E6SSM_slha_io.cpp \
		$(DIR)/E6SSM_physical.cpp \
		$(DIR)/E6SSM_utilities.cpp \
		$(DIR)/E6SSM_two_scale_convergence_tester.cpp \
		$(DIR)/E6SSM_two_scale_high_scale_constraint.cpp \
		$(DIR)/E6SSM_two_scale_initial_guesser.cpp \
		$(DIR)/E6SSM_two_scale_low_scale_constraint.cpp \
		$(DIR)/E6SSM_two_scale_model.cpp \
		$(DIR)/E6SSM_two_scale_model_slha.cpp \
		$(DIR)/E6SSM_two_scale_susy_parameters.cpp \
		$(DIR)/E6SSM_two_scale_soft_parameters.cpp \
		$(DIR)/E6SSM_two_scale_susy_scale_constraint.cpp
EXEE6SSM_SRC += \
		$(DIR)/run_E6SSM.cpp \
		$(DIR)/run_cmd_line_E6SSM.cpp \
		$(DIR)/scan_E6SSM.cpp
LIBE6SSM_HDR += \
		$(DIR)/E6SSM_convergence_tester.hpp \
		$(DIR)/E6SSM_high_scale_constraint.hpp \
		$(DIR)/E6SSM_mass_eigenstates.hpp \
		$(DIR)/E6SSM_info.hpp \
		$(DIR)/E6SSM_initial_guesser.hpp \
		$(DIR)/E6SSM_input_parameters.hpp \
		$(DIR)/E6SSM_low_scale_constraint.hpp \
		$(DIR)/E6SSM_model.hpp \
		$(DIR)/E6SSM_model_slha.hpp \
		$(DIR)/E6SSM_physical.hpp \
		$(DIR)/E6SSM_slha_io.hpp \
		$(DIR)/E6SSM_spectrum_generator.hpp \
		$(DIR)/E6SSM_susy_scale_constraint.hpp \
		$(DIR)/E6SSM_utilities.hpp \
		$(DIR)/E6SSM_two_scale_convergence_tester.hpp \
		$(DIR)/E6SSM_two_scale_high_scale_constraint.hpp \
		$(DIR)/E6SSM_two_scale_initial_guesser.hpp \
		$(DIR)/E6SSM_two_scale_low_scale_constraint.hpp \
		$(DIR)/E6SSM_two_scale_model.hpp \
		$(DIR)/E6SSM_two_scale_model_slha.hpp \
		$(DIR)/E6SSM_two_scale_soft_parameters.hpp \
		$(DIR)/E6SSM_two_scale_susy_parameters.hpp \
		$(DIR)/E6SSM_two_scale_susy_scale_constraint.hpp

ifneq ($(MAKECMDGOALS),showbuild)
ifneq ($(MAKECMDGOALS),tag)
ifneq ($(MAKECMDGOALS),release)
ifneq ($(MAKECMDGOALS),doc)
-include $(E6SSM_TWO_SCALE_SUSY_MK)
-include $(E6SSM_TWO_SCALE_SOFT_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(E6SSM_TWO_SCALE_SUSY_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(E6SSM_TWO_SCALE_SOFT_MK): run-metacode-$(MODNAME)
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
LIBE6SSM_SRC := $(sort $(LIBE6SSM_SRC))
EXEE6SSM_SRC := $(sort $(EXEE6SSM_SRC))

LIBE6SSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBE6SSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBE6SSM_SRC)))

EXEE6SSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEE6SSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEE6SSM_SRC)))

LIBE6SSM_DEP := \
		$(LIBE6SSM_OBJ:.o=.d)

EXEE6SSM_DEP := \
		$(EXEE6SSM_OBJ:.o=.d)

LIBE6SSM     := $(DIR)/lib$(MODNAME)$(LIBEXT)

RUN_E6SSM_OBJ := $(DIR)/run_E6SSM.o
RUN_E6SSM_EXE := $(DIR)/run_E6SSM.x

RUN_CMD_LINE_E6SSM_OBJ := $(DIR)/run_cmd_line_E6SSM.o
RUN_CMD_LINE_E6SSM_EXE := $(DIR)/run_cmd_line_E6SSM.x

SCAN_E6SSM_OBJ := $(DIR)/scan_E6SSM.o
SCAN_E6SSM_EXE := $(DIR)/scan_E6SSM.x

METACODE_STAMP_E6SSM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_E6SSM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-obj \
		distclean-$(MODNAME) run-metacode-$(MODNAME) \
		pack-$(MODNAME)-src

all-$(MODNAME): $(LIBE6SSM)

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(E6SSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBE6SSM_SRC) $(E6SSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBE6SSM_HDR) $(E6SSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXEE6SSM_SRC) $(E6SSM_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(E6SSM_MK) $(E6SSM_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(E6SSM_TWO_SCALE_MK) $(E6SSM_INSTALL_DIR)
ifneq ($(E6SSM_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(E6SSM_SLHA_INPUT) $(E6SSM_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(E6SSM_GNUPLOT) $(E6SSM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBE6SSM_DEP)
		-rm -f $(EXEE6SSM_DEP)

clean-$(MODNAME)-obj:
		-rm -f $(LIBE6SSM_OBJ)
		-rm -f $(EXEE6SSM_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		-rm -f $(LIBE6SSM_SRC)
		-rm -f $(LIBE6SSM_HDR)
		-rm -f $(EXEE6SSM_SRC)
		-rm -f $(METACODE_STAMP_E6SSM)
		-rm -f $(E6SSM_TWO_SCALE_MK)
		-rm -f $(E6SSM_GNUPLOT)

clean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj
		-rm -f $(LIBE6SSM)
		-rm -f $(RUN_E6SSM_EXE)
		-rm -f $(RUN_CMD_LINE_E6SSM_EXE)
		-rm -f $(SCAN_E6SSM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(E6SSM_TARBALL) \
		$(LIBE6SSM_SRC) $(LIBE6SSM_HDR) \
		$(EXEE6SSM_SRC) \
		$(E6SSM_MK) $(E6SSM_TWO_SCALE_MK) \
		$(E6SSM_SLHA_INPUT) $(E6SSM_GNUPLOT)

$(LIBE6SSM_SRC) $(LIBE6SSM_HDR) $(EXEE6SSM_SRC) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_E6SSM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_E6SSM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_E6SSM)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_E6SSM)"
		@echo "Note: to regenerate E6SSM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_E6SSM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_E6SSM):
		@true
endif

$(LIBE6SSM_DEP) $(EXEE6SSM_DEP) $(LIBE6SSM_OBJ) $(EXEE6SSM_OBJ): CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBE6SSM_DEP) $(EXEE6SSM_DEP) $(LIBE6SSM_OBJ) $(EXEE6SSM_OBJ): CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LIBE6SSM): $(LIBE6SSM_OBJ)
		$(MAKELIB) $@ $^

$(RUN_E6SSM_EXE): $(RUN_E6SSM_OBJ) $(LIBE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS)

$(RUN_CMD_LINE_E6SSM_EXE): $(RUN_CMD_LINE_E6SSM_OBJ) $(LIBE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS)

$(SCAN_E6SSM_EXE): $(SCAN_E6SSM_OBJ) $(LIBE6SSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS)

ALLDEP += $(LIBE6SSM_DEP) $(EXEE6SSM_DEP)
ALLSRC += $(LIBE6SSM_SRC) $(EXEE6SSM_SRC)
ALLLIB += $(LIBE6SSM)
ALLEXE += $(RUN_E6SSM_EXE) $(RUN_CMD_LINE_E6SSM_EXE) $(SCAN_E6SSM_EXE)
