DIR          := models/SMSSM
MODNAME      := SMSSM
SARAH_MODEL  := SMSSM

SMSSM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

SMSSM_MK     := \
		$(DIR)/module.mk

SMSSM_TWO_SCALE_SUSY_MK := \
		$(DIR)/two_scale_susy.mk

SMSSM_TWO_SCALE_SOFT_MK := \
		$(DIR)/two_scale_soft.mk

SMSSM_TWO_SCALE_MK := \
		$(SMSSM_TWO_SCALE_SUSY_MK) \
		$(SMSSM_TWO_SCALE_SOFT_MK)

SMSSM_SLHA_INPUT := \
		$(DIR)/LesHouches.in.SMSSM

SMSSM_GNUPLOT := \
		$(DIR)/SMSSM_plot_rgflow.gnuplot \
		$(DIR)/SMSSM_plot_spectrum.gnuplot

SMSSM_TARBALL := \
		$(MODNAME).tar.gz

LIBSMSSM_SRC :=
EXESMSSM_SRC :=

LIBSMSSM_HDR :=

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBSMSSM_SRC += \
		$(DIR)/SMSSM_info.cpp \
		$(DIR)/SMSSM_input_parameters.cpp \
		$(DIR)/SMSSM_slha_io.cpp \
		$(DIR)/SMSSM_physical.cpp \
		$(DIR)/SMSSM_utilities.cpp \
		$(DIR)/SMSSM_two_scale_convergence_tester.cpp \
		$(DIR)/SMSSM_two_scale_high_scale_constraint.cpp \
		$(DIR)/SMSSM_two_scale_initial_guesser.cpp \
		$(DIR)/SMSSM_two_scale_low_scale_constraint.cpp \
		$(DIR)/SMSSM_two_scale_model.cpp \
		$(DIR)/SMSSM_two_scale_model_slha.cpp \
		$(DIR)/SMSSM_two_scale_susy_parameters.cpp \
		$(DIR)/SMSSM_two_scale_soft_parameters.cpp \
		$(DIR)/SMSSM_two_scale_susy_scale_constraint.cpp
EXESMSSM_SRC += \
		$(DIR)/run_SMSSM.cpp \
		$(DIR)/run_cmd_line_SMSSM.cpp \
		$(DIR)/scan_SMSSM.cpp
LIBSMSSM_HDR += \
		$(DIR)/SMSSM_convergence_tester.hpp \
		$(DIR)/SMSSM_high_scale_constraint.hpp \
		$(DIR)/SMSSM_info.hpp \
		$(DIR)/SMSSM_initial_guesser.hpp \
		$(DIR)/SMSSM_input_parameters.hpp \
		$(DIR)/SMSSM_low_scale_constraint.hpp \
		$(DIR)/SMSSM_model.hpp \
		$(DIR)/SMSSM_model_slha.hpp \
		$(DIR)/SMSSM_physical.hpp \
		$(DIR)/SMSSM_slha_io.hpp \
		$(DIR)/SMSSM_spectrum_generator.hpp \
		$(DIR)/SMSSM_susy_scale_constraint.hpp \
		$(DIR)/SMSSM_utilities.hpp \
		$(DIR)/SMSSM_two_scale_convergence_tester.hpp \
		$(DIR)/SMSSM_two_scale_high_scale_constraint.hpp \
		$(DIR)/SMSSM_two_scale_initial_guesser.hpp \
		$(DIR)/SMSSM_two_scale_low_scale_constraint.hpp \
		$(DIR)/SMSSM_two_scale_model.hpp \
		$(DIR)/SMSSM_two_scale_model_slha.hpp \
		$(DIR)/SMSSM_two_scale_soft_parameters.hpp \
		$(DIR)/SMSSM_two_scale_susy_parameters.hpp \
		$(DIR)/SMSSM_two_scale_susy_scale_constraint.hpp

ifneq ($(MAKECMDGOALS),showbuild)
ifneq ($(MAKECMDGOALS),tag)
ifneq ($(MAKECMDGOALS),release)
ifneq ($(MAKECMDGOALS),doc)
-include $(SMSSM_TWO_SCALE_SUSY_MK)
-include $(SMSSM_TWO_SCALE_SOFT_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(SMSSM_TWO_SCALE_SUSY_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(SMSSM_TWO_SCALE_SOFT_MK): run-metacode-$(MODNAME)
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
LIBSMSSM_SRC := $(sort $(LIBSMSSM_SRC))
EXESMSSM_SRC := $(sort $(EXESMSSM_SRC))

LIBSMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBSMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBSMSSM_SRC)))

EXESMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXESMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXESMSSM_SRC)))

LIBSMSSM_DEP := \
		$(LIBSMSSM_OBJ:.o=.d)

EXESMSSM_DEP := \
		$(EXESMSSM_OBJ:.o=.d)

LIBSMSSM     := $(DIR)/lib$(MODNAME)$(LIBEXT)

RUN_SMSSM_OBJ := $(DIR)/run_SMSSM.o
RUN_SMSSM_EXE := $(DIR)/run_SMSSM.x

RUN_CMD_LINE_SMSSM_OBJ := $(DIR)/run_cmd_line_SMSSM.o
RUN_CMD_LINE_SMSSM_EXE := $(DIR)/run_cmd_line_SMSSM.x

SCAN_SMSSM_OBJ := $(DIR)/scan_SMSSM.o
SCAN_SMSSM_EXE := $(DIR)/scan_SMSSM.x

METACODE_STAMP_SMSSM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_SMSSM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-obj \
		distclean-$(MODNAME) run-metacode-$(MODNAME) \
		pack-$(MODNAME)-src

all-$(MODNAME): $(LIBSMSSM)

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(SMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBSMSSM_SRC) $(SMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBSMSSM_HDR) $(SMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXESMSSM_SRC) $(SMSSM_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(SMSSM_MK) $(SMSSM_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(SMSSM_TWO_SCALE_MK) $(SMSSM_INSTALL_DIR)
ifneq ($(SMSSM_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(SMSSM_SLHA_INPUT) $(SMSSM_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(SMSSM_GNUPLOT) $(SMSSM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBSMSSM_DEP)
		-rm -f $(EXESMSSM_DEP)

clean-$(MODNAME)-obj:
		-rm -f $(LIBSMSSM_OBJ)
		-rm -f $(EXESMSSM_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		-rm -f $(LIBSMSSM_SRC)
		-rm -f $(LIBSMSSM_HDR)
		-rm -f $(EXESMSSM_SRC)
		-rm -f $(METACODE_STAMP_SMSSM)
		-rm -f $(SMSSM_TWO_SCALE_MK)
		-rm -f $(SMSSM_GNUPLOT)

clean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj
		-rm -f $(LIBSMSSM)
		-rm -f $(RUN_SMSSM_EXE)
		-rm -f $(RUN_CMD_LINE_SMSSM_EXE)
		-rm -f $(SCAN_SMSSM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(SMSSM_TARBALL) \
		$(LIBSMSSM_SRC) $(LIBSMSSM_HDR) \
		$(EXESMSSM_SRC) \
		$(SMSSM_MK) $(SMSSM_TWO_SCALE_MK) \
		$(SMSSM_SLHA_INPUT) $(SMSSM_GNUPLOT)

$(LIBSMSSM_SRC) $(LIBSMSSM_HDR) $(EXESMSSM_SRC) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_SMSSM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_SMSSM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_SMSSM)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_SMSSM)"
		@echo "Note: to regenerate SMSSM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_SMSSM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_SMSSM):
		@true
endif

$(LIBSMSSM_DEP) $(EXESMSSM_DEP) $(LIBSMSSM_OBJ) $(EXESMSSM_OBJ): CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBSMSSM_DEP) $(EXESMSSM_DEP) $(LIBSMSSM_OBJ) $(EXESMSSM_OBJ): CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LIBSMSSM): $(LIBSMSSM_OBJ)
		$(MAKELIB) $@ $^

$(RUN_SMSSM_EXE): $(RUN_SMSSM_OBJ) $(LIBSMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS)

$(RUN_CMD_LINE_SMSSM_EXE): $(RUN_CMD_LINE_SMSSM_OBJ) $(LIBSMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS)

$(SCAN_SMSSM_EXE): $(SCAN_SMSSM_OBJ) $(LIBSMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS)

ALLDEP += $(LIBSMSSM_DEP) $(EXESMSSM_DEP)
ALLSRC += $(LIBSMSSM_SRC) $(EXESMSSM_SRC)
ALLLIB += $(LIBSMSSM)
ALLEXE += $(RUN_SMSSM_EXE) $(RUN_CMD_LINE_SMSSM_EXE) $(SCAN_SMSSM_EXE)
