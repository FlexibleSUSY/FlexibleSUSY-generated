DIR          := models/NUTNMSSM
MODNAME      := NUTNMSSM
SARAH_MODEL  := NMSSM

NUTNMSSM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

NUTNMSSM_MK     := \
		$(DIR)/module.mk

NUTNMSSM_TWO_SCALE_SUSY_MK := \
		$(DIR)/two_scale_susy.mk

NUTNMSSM_TWO_SCALE_SOFT_MK := \
		$(DIR)/two_scale_soft.mk

NUTNMSSM_TWO_SCALE_MK := \
		$(NUTNMSSM_TWO_SCALE_SUSY_MK) \
		$(NUTNMSSM_TWO_SCALE_SOFT_MK)

NUTNMSSM_SLHA_INPUT := \
		$(DIR)/LesHouches.in.NUTNMSSM_1308.1333_BP2 \
		$(DIR)/LesHouches.in.NUTNMSSM_1308.1333_BP1 \
		$(DIR)/LesHouches.in.NUTNMSSM_1308.1333_BP3

NUTNMSSM_GNUPLOT := \
		$(DIR)/NUTNMSSM_plot_rgflow.gnuplot \
		$(DIR)/NUTNMSSM_plot_spectrum.gnuplot

NUTNMSSM_TARBALL := \
		$(MODNAME).tar.gz

LIBNUTNMSSM_SRC :=
EXENUTNMSSM_SRC :=

LIBNUTNMSSM_HDR :=

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBNUTNMSSM_SRC += \
		$(DIR)/NUTNMSSM_mass_eigenstates.cpp \
		$(DIR)/NUTNMSSM_info.cpp \
		$(DIR)/NUTNMSSM_input_parameters.cpp \
		$(DIR)/NUTNMSSM_slha_io.cpp \
		$(DIR)/NUTNMSSM_physical.cpp \
		$(DIR)/NUTNMSSM_utilities.cpp \
		$(DIR)/NUTNMSSM_two_scale_convergence_tester.cpp \
		$(DIR)/NUTNMSSM_two_scale_high_scale_constraint.cpp \
		$(DIR)/NUTNMSSM_two_scale_initial_guesser.cpp \
		$(DIR)/NUTNMSSM_two_scale_low_scale_constraint.cpp \
		$(DIR)/NUTNMSSM_two_scale_model.cpp \
		$(DIR)/NUTNMSSM_two_scale_model_slha.cpp \
		$(DIR)/NUTNMSSM_two_scale_susy_parameters.cpp \
		$(DIR)/NUTNMSSM_two_scale_soft_parameters.cpp \
		$(DIR)/NUTNMSSM_two_scale_susy_scale_constraint.cpp
EXENUTNMSSM_SRC += \
		$(DIR)/run_NUTNMSSM.cpp \
		$(DIR)/run_cmd_line_NUTNMSSM.cpp \
		$(DIR)/scan_NUTNMSSM.cpp
LIBNUTNMSSM_HDR += \
		$(DIR)/NUTNMSSM_convergence_tester.hpp \
		$(DIR)/NUTNMSSM_high_scale_constraint.hpp \
		$(DIR)/NUTNMSSM_mass_eigenstates.hpp \
		$(DIR)/NUTNMSSM_info.hpp \
		$(DIR)/NUTNMSSM_initial_guesser.hpp \
		$(DIR)/NUTNMSSM_input_parameters.hpp \
		$(DIR)/NUTNMSSM_low_scale_constraint.hpp \
		$(DIR)/NUTNMSSM_model.hpp \
		$(DIR)/NUTNMSSM_model_slha.hpp \
		$(DIR)/NUTNMSSM_physical.hpp \
		$(DIR)/NUTNMSSM_slha_io.hpp \
		$(DIR)/NUTNMSSM_spectrum_generator.hpp \
		$(DIR)/NUTNMSSM_susy_scale_constraint.hpp \
		$(DIR)/NUTNMSSM_utilities.hpp \
		$(DIR)/NUTNMSSM_two_scale_convergence_tester.hpp \
		$(DIR)/NUTNMSSM_two_scale_high_scale_constraint.hpp \
		$(DIR)/NUTNMSSM_two_scale_initial_guesser.hpp \
		$(DIR)/NUTNMSSM_two_scale_low_scale_constraint.hpp \
		$(DIR)/NUTNMSSM_two_scale_model.hpp \
		$(DIR)/NUTNMSSM_two_scale_model_slha.hpp \
		$(DIR)/NUTNMSSM_two_scale_soft_parameters.hpp \
		$(DIR)/NUTNMSSM_two_scale_susy_parameters.hpp \
		$(DIR)/NUTNMSSM_two_scale_susy_scale_constraint.hpp

ifneq ($(MAKECMDGOALS),showbuild)
ifneq ($(MAKECMDGOALS),tag)
ifneq ($(MAKECMDGOALS),release)
ifneq ($(MAKECMDGOALS),doc)
-include $(NUTNMSSM_TWO_SCALE_SUSY_MK)
-include $(NUTNMSSM_TWO_SCALE_SOFT_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(NUTNMSSM_TWO_SCALE_SUSY_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(NUTNMSSM_TWO_SCALE_SOFT_MK): run-metacode-$(MODNAME)
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
LIBNUTNMSSM_SRC := $(sort $(LIBNUTNMSSM_SRC))
EXENUTNMSSM_SRC := $(sort $(EXENUTNMSSM_SRC))

LIBNUTNMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBNUTNMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBNUTNMSSM_SRC)))

EXENUTNMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXENUTNMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXENUTNMSSM_SRC)))

LIBNUTNMSSM_DEP := \
		$(LIBNUTNMSSM_OBJ:.o=.d)

EXENUTNMSSM_DEP := \
		$(EXENUTNMSSM_OBJ:.o=.d)

LIBNUTNMSSM     := $(DIR)/lib$(MODNAME)$(LIBEXT)

RUN_NUTNMSSM_OBJ := $(DIR)/run_NUTNMSSM.o
RUN_NUTNMSSM_EXE := $(DIR)/run_NUTNMSSM.x

RUN_CMD_LINE_NUTNMSSM_OBJ := $(DIR)/run_cmd_line_NUTNMSSM.o
RUN_CMD_LINE_NUTNMSSM_EXE := $(DIR)/run_cmd_line_NUTNMSSM.x

SCAN_NUTNMSSM_OBJ := $(DIR)/scan_NUTNMSSM.o
SCAN_NUTNMSSM_EXE := $(DIR)/scan_NUTNMSSM.x

METACODE_STAMP_NUTNMSSM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_NUTNMSSM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-obj \
		distclean-$(MODNAME) run-metacode-$(MODNAME) \
		pack-$(MODNAME)-src

all-$(MODNAME): $(LIBNUTNMSSM)

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(NUTNMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBNUTNMSSM_SRC) $(NUTNMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBNUTNMSSM_HDR) $(NUTNMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXENUTNMSSM_SRC) $(NUTNMSSM_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(NUTNMSSM_MK) $(NUTNMSSM_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(NUTNMSSM_TWO_SCALE_MK) $(NUTNMSSM_INSTALL_DIR)
ifneq ($(NUTNMSSM_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(NUTNMSSM_SLHA_INPUT) $(NUTNMSSM_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(NUTNMSSM_GNUPLOT) $(NUTNMSSM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBNUTNMSSM_DEP)
		-rm -f $(EXENUTNMSSM_DEP)

clean-$(MODNAME)-obj:
		-rm -f $(LIBNUTNMSSM_OBJ)
		-rm -f $(EXENUTNMSSM_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		-rm -f $(LIBNUTNMSSM_SRC)
		-rm -f $(LIBNUTNMSSM_HDR)
		-rm -f $(EXENUTNMSSM_SRC)
		-rm -f $(METACODE_STAMP_NUTNMSSM)
		-rm -f $(NUTNMSSM_TWO_SCALE_MK)
		-rm -f $(NUTNMSSM_GNUPLOT)

clean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj
		-rm -f $(LIBNUTNMSSM)
		-rm -f $(RUN_NUTNMSSM_EXE)
		-rm -f $(RUN_CMD_LINE_NUTNMSSM_EXE)
		-rm -f $(SCAN_NUTNMSSM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(NUTNMSSM_TARBALL) \
		$(LIBNUTNMSSM_SRC) $(LIBNUTNMSSM_HDR) \
		$(EXENUTNMSSM_SRC) \
		$(NUTNMSSM_MK) $(NUTNMSSM_TWO_SCALE_MK) \
		$(NUTNMSSM_SLHA_INPUT) $(NUTNMSSM_GNUPLOT)

$(LIBNUTNMSSM_SRC) $(LIBNUTNMSSM_HDR) $(EXENUTNMSSM_SRC) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_NUTNMSSM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_NUTNMSSM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_NUTNMSSM)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_NUTNMSSM)"
		@echo "Note: to regenerate NUTNMSSM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_NUTNMSSM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_NUTNMSSM):
		@true
endif

$(LIBNUTNMSSM_DEP) $(EXENUTNMSSM_DEP) $(LIBNUTNMSSM_OBJ) $(EXENUTNMSSM_OBJ): CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBNUTNMSSM_DEP) $(EXENUTNMSSM_DEP) $(LIBNUTNMSSM_OBJ) $(EXENUTNMSSM_OBJ): CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LIBNUTNMSSM): $(LIBNUTNMSSM_OBJ)
		$(MAKELIB) $@ $^

$(RUN_NUTNMSSM_EXE): $(RUN_NUTNMSSM_OBJ) $(LIBNUTNMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS)

$(RUN_CMD_LINE_NUTNMSSM_EXE): $(RUN_CMD_LINE_NUTNMSSM_OBJ) $(LIBNUTNMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS)

$(SCAN_NUTNMSSM_EXE): $(SCAN_NUTNMSSM_OBJ) $(LIBNUTNMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS)

ALLDEP += $(LIBNUTNMSSM_DEP) $(EXENUTNMSSM_DEP)
ALLSRC += $(LIBNUTNMSSM_SRC) $(EXENUTNMSSM_SRC)
ALLLIB += $(LIBNUTNMSSM)
ALLEXE += $(RUN_NUTNMSSM_EXE) $(RUN_CMD_LINE_NUTNMSSM_EXE) $(SCAN_NUTNMSSM_EXE)
