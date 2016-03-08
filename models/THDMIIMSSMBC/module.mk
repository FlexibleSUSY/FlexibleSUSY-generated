DIR          := models/THDMIIMSSMBC
MODNAME      := THDMIIMSSMBC
SARAH_MODEL  := THDM-II

THDMIIMSSMBC_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

THDMIIMSSMBC_MK     := \
		$(DIR)/module.mk

THDMIIMSSMBC_TWO_SCALE_SUSY_MK := \
		$(DIR)/two_scale_susy.mk

THDMIIMSSMBC_TWO_SCALE_SOFT_MK := \
		$(DIR)/two_scale_soft.mk

THDMIIMSSMBC_TWO_SCALE_MK := \
		$(THDMIIMSSMBC_TWO_SCALE_SUSY_MK) \
		$(THDMIIMSSMBC_TWO_SCALE_SOFT_MK)

THDMIIMSSMBC_SLHA_INPUT := \
		$(DIR)/LesHouches.in.THDMIIMSSMBC

THDMIIMSSMBC_GNUPLOT := \
		$(DIR)/THDMIIMSSMBC_plot_rgflow.gnuplot \
		$(DIR)/THDMIIMSSMBC_plot_spectrum.gnuplot

THDMIIMSSMBC_TARBALL := \
		$(MODNAME).tar.gz

LIBTHDMIIMSSMBC_SRC :=
EXETHDMIIMSSMBC_SRC :=

LIBTHDMIIMSSMBC_HDR :=

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBTHDMIIMSSMBC_SRC += \
		$(DIR)/THDMIIMSSMBC_effective_couplings.cpp \
		$(DIR)/THDMIIMSSMBC_mass_eigenstates.cpp \
		$(DIR)/THDMIIMSSMBC_info.cpp \
		$(DIR)/THDMIIMSSMBC_input_parameters.cpp \
		$(DIR)/THDMIIMSSMBC_observables.cpp \
		$(DIR)/THDMIIMSSMBC_slha_io.cpp \
		$(DIR)/THDMIIMSSMBC_physical.cpp \
		$(DIR)/THDMIIMSSMBC_utilities.cpp \
		$(DIR)/THDMIIMSSMBC_two_scale_convergence_tester.cpp \
		$(DIR)/THDMIIMSSMBC_two_scale_high_scale_constraint.cpp \
		$(DIR)/THDMIIMSSMBC_two_scale_initial_guesser.cpp \
		$(DIR)/THDMIIMSSMBC_two_scale_low_scale_constraint.cpp \
		$(DIR)/THDMIIMSSMBC_two_scale_model.cpp \
		$(DIR)/THDMIIMSSMBC_two_scale_model_slha.cpp \
		$(DIR)/THDMIIMSSMBC_two_scale_susy_parameters.cpp \
		$(DIR)/THDMIIMSSMBC_two_scale_soft_parameters.cpp \
		$(DIR)/THDMIIMSSMBC_two_scale_susy_scale_constraint.cpp
EXETHDMIIMSSMBC_SRC += \
		$(DIR)/run_THDMIIMSSMBC.cpp \
		$(DIR)/run_cmd_line_THDMIIMSSMBC.cpp \
		$(DIR)/scan_THDMIIMSSMBC.cpp
LIBTHDMIIMSSMBC_HDR += \
		$(DIR)/THDMIIMSSMBC_convergence_tester.hpp \
		$(DIR)/THDMIIMSSMBC_effective_couplings.hpp \
		$(DIR)/THDMIIMSSMBC_high_scale_constraint.hpp \
		$(DIR)/THDMIIMSSMBC_mass_eigenstates.hpp \
		$(DIR)/THDMIIMSSMBC_info.hpp \
		$(DIR)/THDMIIMSSMBC_initial_guesser.hpp \
		$(DIR)/THDMIIMSSMBC_input_parameters.hpp \
		$(DIR)/THDMIIMSSMBC_low_scale_constraint.hpp \
		$(DIR)/THDMIIMSSMBC_model.hpp \
		$(DIR)/THDMIIMSSMBC_model_slha.hpp \
		$(DIR)/THDMIIMSSMBC_observables.hpp \
		$(DIR)/THDMIIMSSMBC_physical.hpp \
		$(DIR)/THDMIIMSSMBC_slha_io.hpp \
		$(DIR)/THDMIIMSSMBC_spectrum_generator_interface.hpp \
		$(DIR)/THDMIIMSSMBC_spectrum_generator.hpp \
		$(DIR)/THDMIIMSSMBC_susy_scale_constraint.hpp \
		$(DIR)/THDMIIMSSMBC_utilities.hpp \
		$(DIR)/THDMIIMSSMBC_two_scale_convergence_tester.hpp \
		$(DIR)/THDMIIMSSMBC_two_scale_high_scale_constraint.hpp \
		$(DIR)/THDMIIMSSMBC_two_scale_initial_guesser.hpp \
		$(DIR)/THDMIIMSSMBC_two_scale_low_scale_constraint.hpp \
		$(DIR)/THDMIIMSSMBC_two_scale_model.hpp \
		$(DIR)/THDMIIMSSMBC_two_scale_model_slha.hpp \
		$(DIR)/THDMIIMSSMBC_two_scale_soft_parameters.hpp \
		$(DIR)/THDMIIMSSMBC_two_scale_susy_parameters.hpp \
		$(DIR)/THDMIIMSSMBC_two_scale_susy_scale_constraint.hpp

ifneq ($(MAKECMDGOALS),showbuild)
ifneq ($(MAKECMDGOALS),tag)
ifneq ($(MAKECMDGOALS),release)
ifneq ($(MAKECMDGOALS),doc)
-include $(THDMIIMSSMBC_TWO_SCALE_SUSY_MK)
-include $(THDMIIMSSMBC_TWO_SCALE_SOFT_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(THDMIIMSSMBC_TWO_SCALE_SUSY_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(THDMIIMSSMBC_TWO_SCALE_SOFT_MK): run-metacode-$(MODNAME)
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
LIBTHDMIIMSSMBC_SRC := $(sort $(LIBTHDMIIMSSMBC_SRC))
EXETHDMIIMSSMBC_SRC := $(sort $(EXETHDMIIMSSMBC_SRC))

LIBTHDMIIMSSMBC_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBTHDMIIMSSMBC_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBTHDMIIMSSMBC_SRC)))

EXETHDMIIMSSMBC_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXETHDMIIMSSMBC_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXETHDMIIMSSMBC_SRC)))

EXETHDMIIMSSMBC_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXETHDMIIMSSMBC_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXETHDMIIMSSMBC_SRC)))

LIBTHDMIIMSSMBC_DEP := \
		$(LIBTHDMIIMSSMBC_OBJ:.o=.d)

EXETHDMIIMSSMBC_DEP := \
		$(EXETHDMIIMSSMBC_OBJ:.o=.d)

LIBTHDMIIMSSMBC     := $(DIR)/lib$(MODNAME)$(LIBEXT)

METACODE_STAMP_THDMIIMSSMBC := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_THDMIIMSSMBC := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-obj \
		distclean-$(MODNAME) run-metacode-$(MODNAME) \
		pack-$(MODNAME)-src

all-$(MODNAME): $(LIBTHDMIIMSSMBC)

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(THDMIIMSSMBC_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBTHDMIIMSSMBC_SRC) $(THDMIIMSSMBC_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBTHDMIIMSSMBC_HDR) $(THDMIIMSSMBC_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXETHDMIIMSSMBC_SRC) $(THDMIIMSSMBC_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(THDMIIMSSMBC_MK) $(THDMIIMSSMBC_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(THDMIIMSSMBC_TWO_SCALE_MK) $(THDMIIMSSMBC_INSTALL_DIR)
ifneq ($(THDMIIMSSMBC_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(THDMIIMSSMBC_SLHA_INPUT) $(THDMIIMSSMBC_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(THDMIIMSSMBC_GNUPLOT) $(THDMIIMSSMBC_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBTHDMIIMSSMBC_DEP)
		-rm -f $(EXETHDMIIMSSMBC_DEP)

clean-$(MODNAME)-obj:
		-rm -f $(LIBTHDMIIMSSMBC_OBJ)
		-rm -f $(EXETHDMIIMSSMBC_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		-rm -f $(LIBTHDMIIMSSMBC_SRC)
		-rm -f $(LIBTHDMIIMSSMBC_HDR)
		-rm -f $(EXETHDMIIMSSMBC_SRC)
		-rm -f $(METACODE_STAMP_THDMIIMSSMBC)
		-rm -f $(THDMIIMSSMBC_TWO_SCALE_MK)
		-rm -f $(THDMIIMSSMBC_GNUPLOT)

clean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj
		-rm -f $(LIBTHDMIIMSSMBC)
		-rm -f $(EXETHDMIIMSSMBC_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(THDMIIMSSMBC_TARBALL) \
		$(LIBTHDMIIMSSMBC_SRC) $(LIBTHDMIIMSSMBC_HDR) \
		$(EXETHDMIIMSSMBC_SRC) \
		$(THDMIIMSSMBC_MK) $(THDMIIMSSMBC_TWO_SCALE_MK) \
		$(THDMIIMSSMBC_SLHA_INPUT) $(THDMIIMSSMBC_GNUPLOT)

$(LIBTHDMIIMSSMBC_SRC) $(LIBTHDMIIMSSMBC_HDR) $(EXETHDMIIMSSMBC_SRC) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_THDMIIMSSMBC)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_THDMIIMSSMBC): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_THDMIIMSSMBC)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_THDMIIMSSMBC)"
		@echo "Note: to regenerate THDMIIMSSMBC source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_THDMIIMSSMBC)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_THDMIIMSSMBC):
		@true
endif

$(LIBTHDMIIMSSMBC_DEP) $(EXETHDMIIMSSMBC_DEP) $(LIBTHDMIIMSSMBC_OBJ) $(EXETHDMIIMSSMBC_OBJ): CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBTHDMIIMSSMBC_DEP) $(EXETHDMIIMSSMBC_DEP) $(LIBTHDMIIMSSMBC_OBJ) $(EXETHDMIIMSSMBC_OBJ): CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LIBTHDMIIMSSMBC): $(LIBTHDMIIMSSMBC_OBJ)
		$(MAKELIB) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBTHDMIIMSSMBC) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$^) $(ADDONLIBS) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(LDLIBS)

ALLDEP += $(LIBTHDMIIMSSMBC_DEP) $(EXETHDMIIMSSMBC_DEP)
ALLSRC += $(LIBTHDMIIMSSMBC_SRC) $(EXETHDMIIMSSMBC_SRC)
ALLLIB += $(LIBTHDMIIMSSMBC)
ALLEXE += $(EXETHDMIIMSSMBC_EXE)
