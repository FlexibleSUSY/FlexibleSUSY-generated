DIR          := models/NMSSM
MODNAME      := NMSSM
SARAH_MODEL  := NMSSM
WITH_$(MODNAME) := yes

NMSSM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

NMSSM_MK     := \
		$(DIR)/module.mk

NMSSM_TWO_SCALE_SUSY_MK := \
		$(DIR)/two_scale_susy.mk

NMSSM_TWO_SCALE_SOFT_MK := \
		$(DIR)/two_scale_soft.mk

NMSSM_TWO_SCALE_MK := \
		$(NMSSM_TWO_SCALE_SUSY_MK) \
		$(NMSSM_TWO_SCALE_SOFT_MK)

NMSSM_SLHA_INPUT := \
		$(DIR)/LesHouches.in.NMSSM_generated \
		$(DIR)/LesHouches.in.NMSSM

NMSSM_GNUPLOT := \
		$(DIR)/NMSSM_plot_rgflow.gnuplot \
		$(DIR)/NMSSM_plot_spectrum.gnuplot

NMSSM_TARBALL := \
		$(MODNAME).tar.gz

LIBNMSSM_SRC :=
EXENMSSM_SRC :=

LIBNMSSM_HDR :=

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBNMSSM_SRC += \
		$(DIR)/NMSSM_effective_couplings.cpp \
		$(DIR)/NMSSM_mass_eigenstates.cpp \
		$(DIR)/NMSSM_info.cpp \
		$(DIR)/NMSSM_input_parameters.cpp \
		$(DIR)/NMSSM_observables.cpp \
		$(DIR)/NMSSM_slha_io.cpp \
		$(DIR)/NMSSM_physical.cpp \
		$(DIR)/NMSSM_utilities.cpp \
		$(DIR)/NMSSM_two_scale_convergence_tester.cpp \
		$(DIR)/NMSSM_two_scale_high_scale_constraint.cpp \
		$(DIR)/NMSSM_two_scale_initial_guesser.cpp \
		$(DIR)/NMSSM_two_scale_low_scale_constraint.cpp \
		$(DIR)/NMSSM_two_scale_model.cpp \
		$(DIR)/NMSSM_two_scale_model_slha.cpp \
		$(DIR)/NMSSM_two_scale_susy_parameters.cpp \
		$(DIR)/NMSSM_two_scale_soft_parameters.cpp \
		$(DIR)/NMSSM_two_scale_susy_scale_constraint.cpp
EXENMSSM_SRC += \
		$(DIR)/run_NMSSM.cpp \
		$(DIR)/run_cmd_line_NMSSM.cpp \
		$(DIR)/scan_NMSSM.cpp
LIBNMSSM_HDR += \
		$(DIR)/NMSSM_convergence_tester.hpp \
		$(DIR)/NMSSM_effective_couplings.hpp \
		$(DIR)/NMSSM_high_scale_constraint.hpp \
		$(DIR)/NMSSM_mass_eigenstates.hpp \
		$(DIR)/NMSSM_info.hpp \
		$(DIR)/NMSSM_initial_guesser.hpp \
		$(DIR)/NMSSM_input_parameters.hpp \
		$(DIR)/NMSSM_low_scale_constraint.hpp \
		$(DIR)/NMSSM_model.hpp \
		$(DIR)/NMSSM_model_slha.hpp \
		$(DIR)/NMSSM_observables.hpp \
		$(DIR)/NMSSM_physical.hpp \
		$(DIR)/NMSSM_slha_io.hpp \
		$(DIR)/NMSSM_spectrum_generator_interface.hpp \
		$(DIR)/NMSSM_spectrum_generator.hpp \
		$(DIR)/NMSSM_susy_scale_constraint.hpp \
		$(DIR)/NMSSM_utilities.hpp \
		$(DIR)/NMSSM_two_scale_convergence_tester.hpp \
		$(DIR)/NMSSM_two_scale_high_scale_constraint.hpp \
		$(DIR)/NMSSM_two_scale_initial_guesser.hpp \
		$(DIR)/NMSSM_two_scale_low_scale_constraint.hpp \
		$(DIR)/NMSSM_two_scale_model.hpp \
		$(DIR)/NMSSM_two_scale_model_slha.hpp \
		$(DIR)/NMSSM_two_scale_soft_parameters.hpp \
		$(DIR)/NMSSM_two_scale_susy_parameters.hpp \
		$(DIR)/NMSSM_two_scale_susy_scale_constraint.hpp

ifneq ($(MAKECMDGOALS),showbuild)
ifneq ($(MAKECMDGOALS),tag)
ifneq ($(MAKECMDGOALS),release)
ifneq ($(MAKECMDGOALS),doc)
-include $(NMSSM_TWO_SCALE_SUSY_MK)
-include $(NMSSM_TWO_SCALE_SOFT_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(NMSSM_TWO_SCALE_SUSY_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(NMSSM_TWO_SCALE_SOFT_MK): run-metacode-$(MODNAME)
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
LIBNMSSM_SRC := $(sort $(LIBNMSSM_SRC))
EXENMSSM_SRC := $(sort $(EXENMSSM_SRC))

LIBNMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBNMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBNMSSM_SRC)))

EXENMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXENMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXENMSSM_SRC)))

EXENMSSM_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXENMSSM_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXENMSSM_SRC)))

LIBNMSSM_DEP := \
		$(LIBNMSSM_OBJ:.o=.d)

EXENMSSM_DEP := \
		$(EXENMSSM_OBJ:.o=.d)

LIBNMSSM     := $(DIR)/lib$(MODNAME)$(LIBEXT)

METACODE_STAMP_NMSSM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_NMSSM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBNMSSM) $(EXENMSSM_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(NMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBNMSSM_SRC) $(NMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBNMSSM_HDR) $(NMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXENMSSM_SRC) $(NMSSM_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(NMSSM_MK) $(NMSSM_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(NMSSM_TWO_SCALE_MK) $(NMSSM_INSTALL_DIR)
ifneq ($(NMSSM_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(NMSSM_SLHA_INPUT) $(NMSSM_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(NMSSM_GNUPLOT) $(NMSSM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBNMSSM_DEP)
		-rm -f $(EXENMSSM_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBNMSSM)

clean-$(MODNAME)-obj:
		-rm -f $(LIBNMSSM_OBJ)
		-rm -f $(EXENMSSM_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		-rm -f $(LIBNMSSM_SRC)
		-rm -f $(LIBNMSSM_HDR)
		-rm -f $(EXENMSSM_SRC)
		-rm -f $(METACODE_STAMP_NMSSM)
		-rm -f $(NMSSM_TWO_SCALE_MK)
		-rm -f $(NMSSM_SLHA_INPUT)
		-rm -f $(NMSSM_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		-rm -f $(EXENMSSM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(NMSSM_TARBALL) \
		$(LIBNMSSM_SRC) $(LIBNMSSM_HDR) \
		$(EXENMSSM_SRC) \
		$(NMSSM_MK) $(NMSSM_TWO_SCALE_MK) \
		$(NMSSM_SLHA_INPUT) $(NMSSM_GNUPLOT)

$(LIBNMSSM_SRC) $(LIBNMSSM_HDR) $(EXENMSSM_SRC) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_NMSSM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_NMSSM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_NMSSM)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_NMSSM)"
		@echo "Note: to regenerate NMSSM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_NMSSM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_NMSSM):
		@true
endif

$(LIBNMSSM_DEP) $(EXENMSSM_DEP) $(LIBNMSSM_OBJ) $(EXENMSSM_OBJ): CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBNMSSM_DEP) $(EXENMSSM_DEP) $(LIBNMSSM_OBJ) $(EXENMSSM_OBJ): CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LIBNMSSM): $(LIBNMSSM_OBJ)
		$(MAKELIB) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBNMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$^ $(ADDONLIBS)) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(LDLIBS)

ALLDEP += $(LIBNMSSM_DEP) $(EXENMSSM_DEP)
ALLSRC += $(LIBNMSSM_SRC) $(EXENMSSM_SRC)
ALLLIB += $(LIBNMSSM)
ALLEXE += $(EXENMSSM_EXE)
