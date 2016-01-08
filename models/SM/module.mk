DIR          := models/SM
MODNAME      := SM
SARAH_MODEL  := SM

SM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

SM_MK     := \
		$(DIR)/module.mk

SM_TWO_SCALE_SUSY_MK := \
		$(DIR)/two_scale_susy.mk

SM_TWO_SCALE_SOFT_MK := \
		$(DIR)/two_scale_soft.mk

SM_TWO_SCALE_MK := \
		$(SM_TWO_SCALE_SUSY_MK) \
		$(SM_TWO_SCALE_SOFT_MK)

SM_SLHA_INPUT := \
		$(DIR)/LesHouches.in.SM

SM_GNUPLOT := \
		$(DIR)/SM_plot_rgflow.gnuplot \
		$(DIR)/SM_plot_spectrum.gnuplot

SM_TARBALL := \
		$(MODNAME).tar.gz

LIBSM_SRC :=
EXESM_SRC :=

LIBSM_HDR :=

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBSM_SRC += \
		$(DIR)/SM_mass_eigenstates.cpp \
		$(DIR)/SM_info.cpp \
		$(DIR)/SM_input_parameters.cpp \
		$(DIR)/SM_observables.cpp \
		$(DIR)/SM_slha_io.cpp \
		$(DIR)/SM_physical.cpp \
		$(DIR)/SM_utilities.cpp \
		$(DIR)/SM_two_scale_convergence_tester.cpp \
		$(DIR)/SM_two_scale_high_scale_constraint.cpp \
		$(DIR)/SM_two_scale_initial_guesser.cpp \
		$(DIR)/SM_two_scale_low_scale_constraint.cpp \
		$(DIR)/SM_two_scale_model.cpp \
		$(DIR)/SM_two_scale_model_slha.cpp \
		$(DIR)/SM_two_scale_susy_parameters.cpp \
		$(DIR)/SM_two_scale_soft_parameters.cpp \
		$(DIR)/SM_two_scale_susy_scale_constraint.cpp
EXESM_SRC += \
		$(DIR)/run_SM.cpp \
		$(DIR)/run_cmd_line_SM.cpp \
		$(DIR)/scan_SM.cpp
LIBSM_HDR += \
		$(DIR)/SM_convergence_tester.hpp \
		$(DIR)/SM_high_scale_constraint.hpp \
		$(DIR)/SM_mass_eigenstates.hpp \
		$(DIR)/SM_info.hpp \
		$(DIR)/SM_initial_guesser.hpp \
		$(DIR)/SM_input_parameters.hpp \
		$(DIR)/SM_low_scale_constraint.hpp \
		$(DIR)/SM_model.hpp \
		$(DIR)/SM_model_slha.hpp \
		$(DIR)/SM_observables.hpp \
		$(DIR)/SM_physical.hpp \
		$(DIR)/SM_slha_io.hpp \
		$(DIR)/SM_spectrum_generator_interface.hpp \
		$(DIR)/SM_spectrum_generator.hpp \
		$(DIR)/SM_susy_scale_constraint.hpp \
		$(DIR)/SM_utilities.hpp \
		$(DIR)/SM_two_scale_convergence_tester.hpp \
		$(DIR)/SM_two_scale_high_scale_constraint.hpp \
		$(DIR)/SM_two_scale_initial_guesser.hpp \
		$(DIR)/SM_two_scale_low_scale_constraint.hpp \
		$(DIR)/SM_two_scale_model.hpp \
		$(DIR)/SM_two_scale_model_slha.hpp \
		$(DIR)/SM_two_scale_soft_parameters.hpp \
		$(DIR)/SM_two_scale_susy_parameters.hpp \
		$(DIR)/SM_two_scale_susy_scale_constraint.hpp

ifneq ($(MAKECMDGOALS),showbuild)
ifneq ($(MAKECMDGOALS),tag)
ifneq ($(MAKECMDGOALS),release)
ifneq ($(MAKECMDGOALS),doc)
-include $(SM_TWO_SCALE_SUSY_MK)
-include $(SM_TWO_SCALE_SOFT_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(SM_TWO_SCALE_SUSY_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(SM_TWO_SCALE_SOFT_MK): run-metacode-$(MODNAME)
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
LIBSM_SRC := $(sort $(LIBSM_SRC))
EXESM_SRC := $(sort $(EXESM_SRC))

LIBSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBSM_SRC)))

EXESM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXESM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXESM_SRC)))

EXESM_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXESM_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXESM_SRC)))

LIBSM_DEP := \
		$(LIBSM_OBJ:.o=.d)

EXESM_DEP := \
		$(EXESM_OBJ:.o=.d)

LIBSM     := $(DIR)/lib$(MODNAME)$(LIBEXT)

METACODE_STAMP_SM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_SM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-obj \
		distclean-$(MODNAME) run-metacode-$(MODNAME) \
		pack-$(MODNAME)-src

all-$(MODNAME): $(LIBSM)

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(SM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBSM_SRC) $(SM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBSM_HDR) $(SM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXESM_SRC) $(SM_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(SM_MK) $(SM_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(SM_TWO_SCALE_MK) $(SM_INSTALL_DIR)
ifneq ($(SM_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(SM_SLHA_INPUT) $(SM_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(SM_GNUPLOT) $(SM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBSM_DEP)
		-rm -f $(EXESM_DEP)

clean-$(MODNAME)-obj:
		-rm -f $(LIBSM_OBJ)
		-rm -f $(EXESM_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		-rm -f $(LIBSM_SRC)
		-rm -f $(LIBSM_HDR)
		-rm -f $(EXESM_SRC)
		-rm -f $(METACODE_STAMP_SM)
		-rm -f $(SM_TWO_SCALE_MK)
		-rm -f $(SM_GNUPLOT)

clean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj
		-rm -f $(LIBSM)
		-rm -f $(EXESM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(SM_TARBALL) \
		$(LIBSM_SRC) $(LIBSM_HDR) \
		$(EXESM_SRC) \
		$(SM_MK) $(SM_TWO_SCALE_MK) \
		$(SM_SLHA_INPUT) $(SM_GNUPLOT)

$(LIBSM_SRC) $(LIBSM_HDR) $(EXESM_SRC) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_SM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_SM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_SM)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_SM)"
		@echo "Note: to regenerate SM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_SM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_SM):
		@true
endif

$(LIBSM_DEP) $(EXESM_DEP) $(LIBSM_OBJ) $(EXESM_OBJ): CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBSM_DEP) $(EXESM_DEP) $(LIBSM_OBJ) $(EXESM_OBJ): CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LIBSM): $(LIBSM_OBJ)
		$(MAKELIB) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$^) $(ADDONLIBS) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(LDLIBS)

ALLDEP += $(LIBSM_DEP) $(EXESM_DEP)
ALLSRC += $(LIBSM_SRC) $(EXESM_SRC)
ALLLIB += $(LIBSM)
ALLEXE += $(EXESM_EXE)
