DIR          := models/NUTSMSSM
MODNAME      := NUTSMSSM
SARAH_MODEL  := SMSSM

NUTSMSSM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

NUTSMSSM_MK     := \
		$(DIR)/module.mk

NUTSMSSM_TWO_SCALE_SUSY_MK := \
		$(DIR)/two_scale_susy.mk

NUTSMSSM_TWO_SCALE_SOFT_MK := \
		$(DIR)/two_scale_soft.mk

NUTSMSSM_TWO_SCALE_MK := \
		$(NUTSMSSM_TWO_SCALE_SUSY_MK) \
		$(NUTSMSSM_TWO_SCALE_SOFT_MK)

NUTSMSSM_SLHA_INPUT := \
		$(DIR)/LesHouches.in.NUTSMSSM

NUTSMSSM_GNUPLOT := \
		$(DIR)/NUTSMSSM_plot_rgflow.gnuplot \
		$(DIR)/NUTSMSSM_plot_spectrum.gnuplot

NUTSMSSM_TARBALL := \
		$(MODNAME).tar.gz

LIBNUTSMSSM_SRC :=
EXENUTSMSSM_SRC :=

LIBNUTSMSSM_HDR :=

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBNUTSMSSM_SRC += \
		$(DIR)/NUTSMSSM_effective_couplings.cpp \
		$(DIR)/NUTSMSSM_mass_eigenstates.cpp \
		$(DIR)/NUTSMSSM_info.cpp \
		$(DIR)/NUTSMSSM_input_parameters.cpp \
		$(DIR)/NUTSMSSM_observables.cpp \
		$(DIR)/NUTSMSSM_slha_io.cpp \
		$(DIR)/NUTSMSSM_physical.cpp \
		$(DIR)/NUTSMSSM_utilities.cpp \
		$(DIR)/NUTSMSSM_two_scale_convergence_tester.cpp \
		$(DIR)/NUTSMSSM_two_scale_high_scale_constraint.cpp \
		$(DIR)/NUTSMSSM_two_scale_initial_guesser.cpp \
		$(DIR)/NUTSMSSM_two_scale_low_scale_constraint.cpp \
		$(DIR)/NUTSMSSM_two_scale_model.cpp \
		$(DIR)/NUTSMSSM_two_scale_model_slha.cpp \
		$(DIR)/NUTSMSSM_two_scale_susy_parameters.cpp \
		$(DIR)/NUTSMSSM_two_scale_soft_parameters.cpp \
		$(DIR)/NUTSMSSM_two_scale_susy_scale_constraint.cpp
EXENUTSMSSM_SRC += \
		$(DIR)/run_NUTSMSSM.cpp \
		$(DIR)/run_cmd_line_NUTSMSSM.cpp \
		$(DIR)/scan_NUTSMSSM.cpp
LIBNUTSMSSM_HDR += \
		$(DIR)/NUTSMSSM_convergence_tester.hpp \
		$(DIR)/NUTSMSSM_effective_couplings.hpp \
		$(DIR)/NUTSMSSM_high_scale_constraint.hpp \
		$(DIR)/NUTSMSSM_mass_eigenstates.hpp \
		$(DIR)/NUTSMSSM_info.hpp \
		$(DIR)/NUTSMSSM_initial_guesser.hpp \
		$(DIR)/NUTSMSSM_input_parameters.hpp \
		$(DIR)/NUTSMSSM_low_scale_constraint.hpp \
		$(DIR)/NUTSMSSM_model.hpp \
		$(DIR)/NUTSMSSM_model_slha.hpp \
		$(DIR)/NUTSMSSM_observables.hpp \
		$(DIR)/NUTSMSSM_physical.hpp \
		$(DIR)/NUTSMSSM_slha_io.hpp \
		$(DIR)/NUTSMSSM_spectrum_generator_interface.hpp \
		$(DIR)/NUTSMSSM_spectrum_generator.hpp \
		$(DIR)/NUTSMSSM_susy_scale_constraint.hpp \
		$(DIR)/NUTSMSSM_utilities.hpp \
		$(DIR)/NUTSMSSM_two_scale_convergence_tester.hpp \
		$(DIR)/NUTSMSSM_two_scale_high_scale_constraint.hpp \
		$(DIR)/NUTSMSSM_two_scale_initial_guesser.hpp \
		$(DIR)/NUTSMSSM_two_scale_low_scale_constraint.hpp \
		$(DIR)/NUTSMSSM_two_scale_model.hpp \
		$(DIR)/NUTSMSSM_two_scale_model_slha.hpp \
		$(DIR)/NUTSMSSM_two_scale_soft_parameters.hpp \
		$(DIR)/NUTSMSSM_two_scale_susy_parameters.hpp \
		$(DIR)/NUTSMSSM_two_scale_susy_scale_constraint.hpp

ifneq ($(MAKECMDGOALS),showbuild)
ifneq ($(MAKECMDGOALS),tag)
ifneq ($(MAKECMDGOALS),release)
ifneq ($(MAKECMDGOALS),doc)
-include $(NUTSMSSM_TWO_SCALE_SUSY_MK)
-include $(NUTSMSSM_TWO_SCALE_SOFT_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(NUTSMSSM_TWO_SCALE_SUSY_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(NUTSMSSM_TWO_SCALE_SOFT_MK): run-metacode-$(MODNAME)
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
LIBNUTSMSSM_SRC := $(sort $(LIBNUTSMSSM_SRC))
EXENUTSMSSM_SRC := $(sort $(EXENUTSMSSM_SRC))

LIBNUTSMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBNUTSMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBNUTSMSSM_SRC)))

EXENUTSMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXENUTSMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXENUTSMSSM_SRC)))

EXENUTSMSSM_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXENUTSMSSM_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXENUTSMSSM_SRC)))

LIBNUTSMSSM_DEP := \
		$(LIBNUTSMSSM_OBJ:.o=.d)

EXENUTSMSSM_DEP := \
		$(EXENUTSMSSM_OBJ:.o=.d)

LIBNUTSMSSM     := $(DIR)/lib$(MODNAME)$(LIBEXT)

METACODE_STAMP_NUTSMSSM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_NUTSMSSM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-obj \
		distclean-$(MODNAME) run-metacode-$(MODNAME) \
		pack-$(MODNAME)-src

all-$(MODNAME): $(LIBNUTSMSSM)

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(NUTSMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBNUTSMSSM_SRC) $(NUTSMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBNUTSMSSM_HDR) $(NUTSMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXENUTSMSSM_SRC) $(NUTSMSSM_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(NUTSMSSM_MK) $(NUTSMSSM_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(NUTSMSSM_TWO_SCALE_MK) $(NUTSMSSM_INSTALL_DIR)
ifneq ($(NUTSMSSM_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(NUTSMSSM_SLHA_INPUT) $(NUTSMSSM_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(NUTSMSSM_GNUPLOT) $(NUTSMSSM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBNUTSMSSM_DEP)
		-rm -f $(EXENUTSMSSM_DEP)

clean-$(MODNAME)-obj:
		-rm -f $(LIBNUTSMSSM_OBJ)
		-rm -f $(EXENUTSMSSM_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		-rm -f $(LIBNUTSMSSM_SRC)
		-rm -f $(LIBNUTSMSSM_HDR)
		-rm -f $(EXENUTSMSSM_SRC)
		-rm -f $(METACODE_STAMP_NUTSMSSM)
		-rm -f $(NUTSMSSM_TWO_SCALE_MK)
		-rm -f $(NUTSMSSM_GNUPLOT)

clean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj
		-rm -f $(LIBNUTSMSSM)
		-rm -f $(EXENUTSMSSM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(NUTSMSSM_TARBALL) \
		$(LIBNUTSMSSM_SRC) $(LIBNUTSMSSM_HDR) \
		$(EXENUTSMSSM_SRC) \
		$(NUTSMSSM_MK) $(NUTSMSSM_TWO_SCALE_MK) \
		$(NUTSMSSM_SLHA_INPUT) $(NUTSMSSM_GNUPLOT)

$(LIBNUTSMSSM_SRC) $(LIBNUTSMSSM_HDR) $(EXENUTSMSSM_SRC) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_NUTSMSSM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_NUTSMSSM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_NUTSMSSM)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_NUTSMSSM)"
		@echo "Note: to regenerate NUTSMSSM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_NUTSMSSM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_NUTSMSSM):
		@true
endif

$(LIBNUTSMSSM_DEP) $(EXENUTSMSSM_DEP) $(LIBNUTSMSSM_OBJ) $(EXENUTSMSSM_OBJ): CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBNUTSMSSM_DEP) $(EXENUTSMSSM_DEP) $(LIBNUTSMSSM_OBJ) $(EXENUTSMSSM_OBJ): CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LIBNUTSMSSM): $(LIBNUTSMSSM_OBJ)
		$(MAKELIB) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBNUTSMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$^) $(ADDONLIBS) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(LDLIBS)

ALLDEP += $(LIBNUTSMSSM_DEP) $(EXENUTSMSSM_DEP)
ALLSRC += $(LIBNUTSMSSM_SRC) $(EXENUTSMSSM_SRC)
ALLLIB += $(LIBNUTSMSSM)
ALLEXE += $(EXENUTSMSSM_EXE)
