DIR          := models/MSSMRHN
MODNAME      := MSSMRHN
SARAH_MODEL  := MSSMRHN

MSSMRHN_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

MSSMRHN_MK     := \
		$(DIR)/module.mk

MSSMRHN_TWO_SCALE_SUSY_MK := \
		$(DIR)/two_scale_susy.mk

MSSMRHN_TWO_SCALE_SOFT_MK := \
		$(DIR)/two_scale_soft.mk

MSSMRHN_TWO_SCALE_MK := \
		$(MSSMRHN_TWO_SCALE_SUSY_MK) \
		$(MSSMRHN_TWO_SCALE_SOFT_MK)

MSSMRHN_SLHA_INPUT := \


MSSMRHN_GNUPLOT := \
		$(DIR)/MSSMRHN_plot_rgflow.gnuplot \
		$(DIR)/MSSMRHN_plot_spectrum.gnuplot

MSSMRHN_TARBALL := \
		$(MODNAME).tar.gz

LIBMSSMRHN_SRC :=
EXEMSSMRHN_SRC :=

LIBMSSMRHN_HDR :=

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBMSSMRHN_SRC += \
		$(DIR)/MSSMRHN_info.cpp \
		$(DIR)/MSSMRHN_input_parameters.cpp \
		$(DIR)/MSSMRHN_slha_io.cpp \
		$(DIR)/MSSMRHN_physical.cpp \
		$(DIR)/MSSMRHN_utilities.cpp \
		$(DIR)/MSSMRHN_two_scale_convergence_tester.cpp \
		$(DIR)/MSSMRHN_two_scale_high_scale_constraint.cpp \
		$(DIR)/MSSMRHN_two_scale_initial_guesser.cpp \
		$(DIR)/MSSMRHN_two_scale_low_scale_constraint.cpp \
		$(DIR)/MSSMRHN_two_scale_model.cpp \
		$(DIR)/MSSMRHN_two_scale_model_slha.cpp \
		$(DIR)/MSSMRHN_two_scale_susy_parameters.cpp \
		$(DIR)/MSSMRHN_two_scale_soft_parameters.cpp \
		$(DIR)/MSSMRHN_two_scale_susy_scale_constraint.cpp
EXEMSSMRHN_SRC += \
		$(DIR)/run_MSSMRHN.cpp \
		$(DIR)/run_cmd_line_MSSMRHN.cpp \
		$(DIR)/scan_MSSMRHN.cpp
LIBMSSMRHN_HDR += \
		$(DIR)/MSSMRHN_convergence_tester.hpp \
		$(DIR)/MSSMRHN_high_scale_constraint.hpp \
		$(DIR)/MSSMRHN_info.hpp \
		$(DIR)/MSSMRHN_initial_guesser.hpp \
		$(DIR)/MSSMRHN_input_parameters.hpp \
		$(DIR)/MSSMRHN_low_scale_constraint.hpp \
		$(DIR)/MSSMRHN_model.hpp \
		$(DIR)/MSSMRHN_model_slha.hpp \
		$(DIR)/MSSMRHN_physical.hpp \
		$(DIR)/MSSMRHN_slha_io.hpp \
		$(DIR)/MSSMRHN_spectrum_generator.hpp \
		$(DIR)/MSSMRHN_susy_scale_constraint.hpp \
		$(DIR)/MSSMRHN_utilities.hpp \
		$(DIR)/MSSMRHN_two_scale_convergence_tester.hpp \
		$(DIR)/MSSMRHN_two_scale_high_scale_constraint.hpp \
		$(DIR)/MSSMRHN_two_scale_initial_guesser.hpp \
		$(DIR)/MSSMRHN_two_scale_low_scale_constraint.hpp \
		$(DIR)/MSSMRHN_two_scale_model.hpp \
		$(DIR)/MSSMRHN_two_scale_model_slha.hpp \
		$(DIR)/MSSMRHN_two_scale_soft_parameters.hpp \
		$(DIR)/MSSMRHN_two_scale_susy_parameters.hpp \
		$(DIR)/MSSMRHN_two_scale_susy_scale_constraint.hpp

ifneq ($(MAKECMDGOALS),showbuild)
ifneq ($(MAKECMDGOALS),tag)
ifneq ($(MAKECMDGOALS),release)
ifneq ($(MAKECMDGOALS),doc)
-include $(MSSMRHN_TWO_SCALE_SUSY_MK)
-include $(MSSMRHN_TWO_SCALE_SOFT_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(MSSMRHN_TWO_SCALE_SUSY_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(MSSMRHN_TWO_SCALE_SOFT_MK): run-metacode-$(MODNAME)
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
LIBMSSMRHN_SRC := $(sort $(LIBMSSMRHN_SRC))
EXEMSSMRHN_SRC := $(sort $(EXEMSSMRHN_SRC))

LIBMSSMRHN_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBMSSMRHN_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBMSSMRHN_SRC)))

EXEMSSMRHN_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEMSSMRHN_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEMSSMRHN_SRC)))

LIBMSSMRHN_DEP := \
		$(LIBMSSMRHN_OBJ:.o=.d)

EXEMSSMRHN_DEP := \
		$(EXEMSSMRHN_OBJ:.o=.d)

LIBMSSMRHN     := $(DIR)/lib$(MODNAME)$(LIBEXT)

RUN_MSSMRHN_OBJ := $(DIR)/run_MSSMRHN.o
RUN_MSSMRHN_EXE := $(DIR)/run_MSSMRHN.x

RUN_CMD_LINE_MSSMRHN_OBJ := $(DIR)/run_cmd_line_MSSMRHN.o
RUN_CMD_LINE_MSSMRHN_EXE := $(DIR)/run_cmd_line_MSSMRHN.x

SCAN_MSSMRHN_OBJ := $(DIR)/scan_MSSMRHN.o
SCAN_MSSMRHN_EXE := $(DIR)/scan_MSSMRHN.x

METACODE_STAMP_MSSMRHN := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_MSSMRHN := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-obj \
		distclean-$(MODNAME) run-metacode-$(MODNAME) \
		pack-$(MODNAME)-src

all-$(MODNAME): $(LIBMSSMRHN)

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(MSSMRHN_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBMSSMRHN_SRC) $(MSSMRHN_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBMSSMRHN_HDR) $(MSSMRHN_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXEMSSMRHN_SRC) $(MSSMRHN_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(MSSMRHN_MK) $(MSSMRHN_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(MSSMRHN_TWO_SCALE_MK) $(MSSMRHN_INSTALL_DIR)
ifneq ($(MSSMRHN_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(MSSMRHN_SLHA_INPUT) $(MSSMRHN_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(MSSMRHN_GNUPLOT) $(MSSMRHN_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBMSSMRHN_DEP)
		-rm -f $(EXEMSSMRHN_DEP)

clean-$(MODNAME)-obj:
		-rm -f $(LIBMSSMRHN_OBJ)
		-rm -f $(EXEMSSMRHN_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		-rm -f $(LIBMSSMRHN_SRC)
		-rm -f $(LIBMSSMRHN_HDR)
		-rm -f $(EXEMSSMRHN_SRC)
		-rm -f $(METACODE_STAMP_MSSMRHN)
		-rm -f $(MSSMRHN_TWO_SCALE_MK)
		-rm -f $(MSSMRHN_GNUPLOT)

clean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj
		-rm -f $(LIBMSSMRHN)
		-rm -f $(RUN_MSSMRHN_EXE)
		-rm -f $(RUN_CMD_LINE_MSSMRHN_EXE)
		-rm -f $(SCAN_MSSMRHN_EXE)

distclean-$(MODNAME): clean-$(MODNAME)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(MSSMRHN_TARBALL) \
		$(LIBMSSMRHN_SRC) $(LIBMSSMRHN_HDR) \
		$(EXEMSSMRHN_SRC) \
		$(MSSMRHN_MK) $(MSSMRHN_TWO_SCALE_MK) \
		$(MSSMRHN_SLHA_INPUT) $(MSSMRHN_GNUPLOT)

$(LIBMSSMRHN_SRC) $(LIBMSSMRHN_HDR) $(EXEMSSMRHN_SRC) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_MSSMRHN)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_MSSMRHN): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_MSSMRHN)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_MSSMRHN)"
		@echo "Note: to regenerate MSSMRHN source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_MSSMRHN)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_MSSMRHN):
		@true
endif

$(LIBMSSMRHN_DEP) $(EXEMSSMRHN_DEP) $(LIBMSSMRHN_OBJ) $(EXEMSSMRHN_OBJ): CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBMSSMRHN_DEP) $(EXEMSSMRHN_DEP) $(LIBMSSMRHN_OBJ) $(EXEMSSMRHN_OBJ): CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LIBMSSMRHN): $(LIBMSSMRHN_OBJ)
		$(MAKELIB) $@ $^

$(RUN_MSSMRHN_EXE): $(RUN_MSSMRHN_OBJ) $(LIBMSSMRHN) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS)

$(RUN_CMD_LINE_MSSMRHN_EXE): $(RUN_CMD_LINE_MSSMRHN_OBJ) $(LIBMSSMRHN) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS)

$(SCAN_MSSMRHN_EXE): $(SCAN_MSSMRHN_OBJ) $(LIBMSSMRHN) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS)

ALLDEP += $(LIBMSSMRHN_DEP) $(EXEMSSMRHN_DEP)
ALLSRC += $(LIBMSSMRHN_SRC) $(EXEMSSMRHN_SRC)
ALLLIB += $(LIBMSSMRHN)
ALLEXE += $(RUN_MSSMRHN_EXE) $(RUN_CMD_LINE_MSSMRHN_EXE) $(SCAN_MSSMRHN_EXE)
