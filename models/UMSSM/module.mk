DIR          := models/UMSSM
MODNAME      := UMSSM
SARAH_MODEL  := UMSSM

UMSSM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

UMSSM_MK     := \
		$(DIR)/module.mk

UMSSM_TWO_SCALE_SUSY_MK := \
		$(DIR)/two_scale_susy.mk

UMSSM_TWO_SCALE_SOFT_MK := \
		$(DIR)/two_scale_soft.mk

UMSSM_TWO_SCALE_MK := \
		$(UMSSM_TWO_SCALE_SUSY_MK) \
		$(UMSSM_TWO_SCALE_SOFT_MK)

UMSSM_SLHA_INPUT := \
		$(DIR)/LesHouches.in.UMSSM

UMSSM_GNUPLOT := \
		$(DIR)/UMSSM_plot_rgflow.gnuplot \
		$(DIR)/UMSSM_plot_spectrum.gnuplot

UMSSM_TARBALL := \
		$(MODNAME).tar.gz

LIBUMSSM_SRC :=
EXEUMSSM_SRC :=

LIBUMSSM_HDR :=

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBUMSSM_SRC += \
		$(DIR)/UMSSM_mass_eigenstates.cpp \
		$(DIR)/UMSSM_info.cpp \
		$(DIR)/UMSSM_input_parameters.cpp \
		$(DIR)/UMSSM_slha_io.cpp \
		$(DIR)/UMSSM_physical.cpp \
		$(DIR)/UMSSM_utilities.cpp \
		$(DIR)/UMSSM_two_scale_convergence_tester.cpp \
		$(DIR)/UMSSM_two_scale_high_scale_constraint.cpp \
		$(DIR)/UMSSM_two_scale_initial_guesser.cpp \
		$(DIR)/UMSSM_two_scale_low_scale_constraint.cpp \
		$(DIR)/UMSSM_two_scale_model.cpp \
		$(DIR)/UMSSM_two_scale_model_slha.cpp \
		$(DIR)/UMSSM_two_scale_susy_parameters.cpp \
		$(DIR)/UMSSM_two_scale_soft_parameters.cpp \
		$(DIR)/UMSSM_two_scale_susy_scale_constraint.cpp
EXEUMSSM_SRC += \
		$(DIR)/run_UMSSM.cpp \
		$(DIR)/run_cmd_line_UMSSM.cpp \
		$(DIR)/scan_UMSSM.cpp
LIBUMSSM_HDR += \
		$(DIR)/UMSSM_convergence_tester.hpp \
		$(DIR)/UMSSM_high_scale_constraint.hpp \
		$(DIR)/UMSSM_mass_eigenstates.hpp \
		$(DIR)/UMSSM_info.hpp \
		$(DIR)/UMSSM_initial_guesser.hpp \
		$(DIR)/UMSSM_input_parameters.hpp \
		$(DIR)/UMSSM_low_scale_constraint.hpp \
		$(DIR)/UMSSM_model.hpp \
		$(DIR)/UMSSM_model_slha.hpp \
		$(DIR)/UMSSM_physical.hpp \
		$(DIR)/UMSSM_slha_io.hpp \
		$(DIR)/UMSSM_spectrum_generator_interface.hpp \
		$(DIR)/UMSSM_spectrum_generator.hpp \
		$(DIR)/UMSSM_susy_scale_constraint.hpp \
		$(DIR)/UMSSM_utilities.hpp \
		$(DIR)/UMSSM_two_scale_convergence_tester.hpp \
		$(DIR)/UMSSM_two_scale_high_scale_constraint.hpp \
		$(DIR)/UMSSM_two_scale_initial_guesser.hpp \
		$(DIR)/UMSSM_two_scale_low_scale_constraint.hpp \
		$(DIR)/UMSSM_two_scale_model.hpp \
		$(DIR)/UMSSM_two_scale_model_slha.hpp \
		$(DIR)/UMSSM_two_scale_soft_parameters.hpp \
		$(DIR)/UMSSM_two_scale_susy_parameters.hpp \
		$(DIR)/UMSSM_two_scale_susy_scale_constraint.hpp

ifneq ($(MAKECMDGOALS),showbuild)
ifneq ($(MAKECMDGOALS),tag)
ifneq ($(MAKECMDGOALS),release)
ifneq ($(MAKECMDGOALS),doc)
-include $(UMSSM_TWO_SCALE_SUSY_MK)
-include $(UMSSM_TWO_SCALE_SOFT_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(UMSSM_TWO_SCALE_SUSY_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(UMSSM_TWO_SCALE_SOFT_MK): run-metacode-$(MODNAME)
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
LIBUMSSM_SRC := $(sort $(LIBUMSSM_SRC))
EXEUMSSM_SRC := $(sort $(EXEUMSSM_SRC))

LIBUMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBUMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBUMSSM_SRC)))

EXEUMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEUMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEUMSSM_SRC)))

EXEUMSSM_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXEUMSSM_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXEUMSSM_SRC)))

LIBUMSSM_DEP := \
		$(LIBUMSSM_OBJ:.o=.d)

EXEUMSSM_DEP := \
		$(EXEUMSSM_OBJ:.o=.d)

LIBUMSSM     := $(DIR)/lib$(MODNAME)$(LIBEXT)

METACODE_STAMP_UMSSM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_UMSSM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-obj \
		distclean-$(MODNAME) run-metacode-$(MODNAME) \
		pack-$(MODNAME)-src

all-$(MODNAME): $(LIBUMSSM)

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(UMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBUMSSM_SRC) $(UMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBUMSSM_HDR) $(UMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXEUMSSM_SRC) $(UMSSM_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(UMSSM_MK) $(UMSSM_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(UMSSM_TWO_SCALE_MK) $(UMSSM_INSTALL_DIR)
ifneq ($(UMSSM_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(UMSSM_SLHA_INPUT) $(UMSSM_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(UMSSM_GNUPLOT) $(UMSSM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBUMSSM_DEP)
		-rm -f $(EXEUMSSM_DEP)

clean-$(MODNAME)-obj:
		-rm -f $(LIBUMSSM_OBJ)
		-rm -f $(EXEUMSSM_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		-rm -f $(LIBUMSSM_SRC)
		-rm -f $(LIBUMSSM_HDR)
		-rm -f $(EXEUMSSM_SRC)
		-rm -f $(METACODE_STAMP_UMSSM)
		-rm -f $(UMSSM_TWO_SCALE_MK)
		-rm -f $(UMSSM_GNUPLOT)

clean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj
		-rm -f $(LIBUMSSM)
		-rm -f $(EXEUMSSM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(UMSSM_TARBALL) \
		$(LIBUMSSM_SRC) $(LIBUMSSM_HDR) \
		$(EXEUMSSM_SRC) \
		$(UMSSM_MK) $(UMSSM_TWO_SCALE_MK) \
		$(UMSSM_SLHA_INPUT) $(UMSSM_GNUPLOT)

$(LIBUMSSM_SRC) $(LIBUMSSM_HDR) $(EXEUMSSM_SRC) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_UMSSM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_UMSSM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_UMSSM)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_UMSSM)"
		@echo "Note: to regenerate UMSSM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_UMSSM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_UMSSM):
		@true
endif

$(LIBUMSSM_DEP) $(EXEUMSSM_DEP) $(LIBUMSSM_OBJ) $(EXEUMSSM_OBJ): CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBUMSSM_DEP) $(EXEUMSSM_DEP) $(LIBUMSSM_OBJ) $(EXEUMSSM_OBJ): CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LIBUMSSM): $(LIBUMSSM_OBJ)
		$(MAKELIB) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBUMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(LDLIBS)

ALLDEP += $(LIBUMSSM_DEP) $(EXEUMSSM_DEP)
ALLSRC += $(LIBUMSSM_SRC) $(EXEUMSSM_SRC)
ALLLIB += $(LIBUMSSM)
ALLEXE += $(EXEUMSSM_EXE)
