DIR          := models/TMSSM
MODNAME      := TMSSM
SARAH_MODEL  := TMSSM
WITH_$(MODNAME) := yes

TMSSM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

TMSSM_MK     := \
		$(DIR)/module.mk

TMSSM_TWO_SCALE_SUSY_MK := \
		$(DIR)/two_scale_susy.mk

TMSSM_TWO_SCALE_SOFT_MK := \
		$(DIR)/two_scale_soft.mk

TMSSM_TWO_SCALE_MK := \
		$(TMSSM_TWO_SCALE_SUSY_MK) \
		$(TMSSM_TWO_SCALE_SOFT_MK)

TMSSM_SLHA_INPUT := \
		$(DIR)/LesHouches.in.TMSSM_generated \
		$(DIR)/LesHouches.in.TMSSM

TMSSM_GNUPLOT := \
		$(DIR)/TMSSM_plot_rgflow.gnuplot \
		$(DIR)/TMSSM_plot_spectrum.gnuplot

TMSSM_TARBALL := \
		$(MODNAME).tar.gz

LIBTMSSM_SRC :=
EXETMSSM_SRC :=

LIBTMSSM_HDR :=

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBTMSSM_SRC += \
		$(DIR)/TMSSM_effective_couplings.cpp \
		$(DIR)/TMSSM_mass_eigenstates.cpp \
		$(DIR)/TMSSM_info.cpp \
		$(DIR)/TMSSM_input_parameters.cpp \
		$(DIR)/TMSSM_observables.cpp \
		$(DIR)/TMSSM_slha_io.cpp \
		$(DIR)/TMSSM_physical.cpp \
		$(DIR)/TMSSM_utilities.cpp \
		$(DIR)/TMSSM_two_scale_convergence_tester.cpp \
		$(DIR)/TMSSM_two_scale_high_scale_constraint.cpp \
		$(DIR)/TMSSM_two_scale_initial_guesser.cpp \
		$(DIR)/TMSSM_two_scale_low_scale_constraint.cpp \
		$(DIR)/TMSSM_two_scale_model.cpp \
		$(DIR)/TMSSM_two_scale_model_slha.cpp \
		$(DIR)/TMSSM_two_scale_susy_parameters.cpp \
		$(DIR)/TMSSM_two_scale_soft_parameters.cpp \
		$(DIR)/TMSSM_two_scale_susy_scale_constraint.cpp
EXETMSSM_SRC += \
		$(DIR)/run_TMSSM.cpp \
		$(DIR)/run_cmd_line_TMSSM.cpp \
		$(DIR)/scan_TMSSM.cpp
LIBTMSSM_HDR += \
		$(DIR)/TMSSM_convergence_tester.hpp \
		$(DIR)/TMSSM_effective_couplings.hpp \
		$(DIR)/TMSSM_high_scale_constraint.hpp \
		$(DIR)/TMSSM_mass_eigenstates.hpp \
		$(DIR)/TMSSM_info.hpp \
		$(DIR)/TMSSM_initial_guesser.hpp \
		$(DIR)/TMSSM_input_parameters.hpp \
		$(DIR)/TMSSM_low_scale_constraint.hpp \
		$(DIR)/TMSSM_model.hpp \
		$(DIR)/TMSSM_model_slha.hpp \
		$(DIR)/TMSSM_observables.hpp \
		$(DIR)/TMSSM_physical.hpp \
		$(DIR)/TMSSM_slha_io.hpp \
		$(DIR)/TMSSM_spectrum_generator_interface.hpp \
		$(DIR)/TMSSM_spectrum_generator.hpp \
		$(DIR)/TMSSM_susy_scale_constraint.hpp \
		$(DIR)/TMSSM_utilities.hpp \
		$(DIR)/TMSSM_two_scale_convergence_tester.hpp \
		$(DIR)/TMSSM_two_scale_high_scale_constraint.hpp \
		$(DIR)/TMSSM_two_scale_initial_guesser.hpp \
		$(DIR)/TMSSM_two_scale_low_scale_constraint.hpp \
		$(DIR)/TMSSM_two_scale_model.hpp \
		$(DIR)/TMSSM_two_scale_model_slha.hpp \
		$(DIR)/TMSSM_two_scale_soft_parameters.hpp \
		$(DIR)/TMSSM_two_scale_susy_parameters.hpp \
		$(DIR)/TMSSM_two_scale_susy_scale_constraint.hpp

ifneq ($(MAKECMDGOALS),showbuild)
ifneq ($(MAKECMDGOALS),tag)
ifneq ($(MAKECMDGOALS),release)
ifneq ($(MAKECMDGOALS),doc)
-include $(TMSSM_TWO_SCALE_SUSY_MK)
-include $(TMSSM_TWO_SCALE_SOFT_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(TMSSM_TWO_SCALE_SUSY_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(TMSSM_TWO_SCALE_SOFT_MK): run-metacode-$(MODNAME)
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
LIBTMSSM_SRC := $(sort $(LIBTMSSM_SRC))
EXETMSSM_SRC := $(sort $(EXETMSSM_SRC))

LIBTMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBTMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBTMSSM_SRC)))

EXETMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXETMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXETMSSM_SRC)))

EXETMSSM_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXETMSSM_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXETMSSM_SRC)))

LIBTMSSM_DEP := \
		$(LIBTMSSM_OBJ:.o=.d)

EXETMSSM_DEP := \
		$(EXETMSSM_OBJ:.o=.d)

LIBTMSSM     := $(DIR)/lib$(MODNAME)$(LIBEXT)

METACODE_STAMP_TMSSM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_TMSSM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBTMSSM) $(EXETMSSM_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(TMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBTMSSM_SRC) $(TMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBTMSSM_HDR) $(TMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXETMSSM_SRC) $(TMSSM_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(TMSSM_MK) $(TMSSM_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(TMSSM_TWO_SCALE_MK) $(TMSSM_INSTALL_DIR)
ifneq ($(TMSSM_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(TMSSM_SLHA_INPUT) $(TMSSM_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(TMSSM_GNUPLOT) $(TMSSM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBTMSSM_DEP)
		-rm -f $(EXETMSSM_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBTMSSM)

clean-$(MODNAME)-obj:
		-rm -f $(LIBTMSSM_OBJ)
		-rm -f $(EXETMSSM_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		-rm -f $(LIBTMSSM_SRC)
		-rm -f $(LIBTMSSM_HDR)
		-rm -f $(EXETMSSM_SRC)
		-rm -f $(METACODE_STAMP_TMSSM)
		-rm -f $(TMSSM_TWO_SCALE_MK)
		-rm -f $(TMSSM_SLHA_INPUT)
		-rm -f $(TMSSM_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		-rm -f $(EXETMSSM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(TMSSM_TARBALL) \
		$(LIBTMSSM_SRC) $(LIBTMSSM_HDR) \
		$(EXETMSSM_SRC) \
		$(TMSSM_MK) $(TMSSM_TWO_SCALE_MK) \
		$(TMSSM_SLHA_INPUT) $(TMSSM_GNUPLOT)

$(LIBTMSSM_SRC) $(LIBTMSSM_HDR) $(EXETMSSM_SRC) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_TMSSM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_TMSSM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_TMSSM)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_TMSSM)"
		@echo "Note: to regenerate TMSSM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_TMSSM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_TMSSM):
		@true
endif

$(LIBTMSSM_DEP) $(EXETMSSM_DEP) $(LIBTMSSM_OBJ) $(EXETMSSM_OBJ): CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBTMSSM_DEP) $(EXETMSSM_DEP) $(LIBTMSSM_OBJ) $(EXETMSSM_OBJ): CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LIBTMSSM): $(LIBTMSSM_OBJ)
		$(MAKELIB) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBTMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$^ $(ADDONLIBS)) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(LDLIBS)

ALLDEP += $(LIBTMSSM_DEP) $(EXETMSSM_DEP)
ALLSRC += $(LIBTMSSM_SRC) $(EXETMSSM_SRC)
ALLLIB += $(LIBTMSSM)
ALLEXE += $(EXETMSSM_EXE)
