DIR          := models/HSSUSY
MODNAME      := HSSUSY
SARAH_MODEL  := SM
WITH_$(MODNAME) := yes

HSSUSY_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

HSSUSY_MK     := \
		$(DIR)/module.mk

HSSUSY_TWO_SCALE_SUSY_MK := \
		$(DIR)/two_scale_susy.mk

HSSUSY_TWO_SCALE_SOFT_MK := \
		$(DIR)/two_scale_soft.mk

HSSUSY_TWO_SCALE_MK := \
		$(HSSUSY_TWO_SCALE_SUSY_MK) \
		$(HSSUSY_TWO_SCALE_SOFT_MK)

HSSUSY_SLHA_INPUT := \
		$(DIR)/LesHouches.in.HSSUSY_generated \
		$(DIR)/LesHouches.in.HSSUSY

HSSUSY_GNUPLOT := \
		$(DIR)/HSSUSY_plot_rgflow.gnuplot \
		$(DIR)/HSSUSY_plot_spectrum.gnuplot

HSSUSY_TARBALL := \
		$(MODNAME).tar.gz

LIBHSSUSY_SRC :=
EXEHSSUSY_SRC :=
LLHSSUSY_LIB  :=
LLHSSUSY_OBJ  :=
LLHSSUSY_SRC  :=
LLHSSUSY_MMA  :=

LIBHSSUSY_HDR :=

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBHSSUSY_SRC += \
		$(DIR)/HSSUSY_effective_couplings.cpp \
		$(DIR)/HSSUSY_mass_eigenstates.cpp \
		$(DIR)/HSSUSY_info.cpp \
		$(DIR)/HSSUSY_input_parameters.cpp \
		$(DIR)/HSSUSY_observables.cpp \
		$(DIR)/HSSUSY_slha_io.cpp \
		$(DIR)/HSSUSY_physical.cpp \
		$(DIR)/HSSUSY_utilities.cpp \
		$(DIR)/HSSUSY_standard_model_matching.cpp \
		$(DIR)/HSSUSY_standard_model_two_scale_matching.cpp \
		$(DIR)/HSSUSY_two_scale_convergence_tester.cpp \
		$(DIR)/HSSUSY_two_scale_high_scale_constraint.cpp \
		$(DIR)/HSSUSY_two_scale_initial_guesser.cpp \
		$(DIR)/HSSUSY_two_scale_low_scale_constraint.cpp \
		$(DIR)/HSSUSY_two_scale_model.cpp \
		$(DIR)/HSSUSY_two_scale_model_slha.cpp \
		$(DIR)/HSSUSY_two_scale_susy_parameters.cpp \
		$(DIR)/HSSUSY_two_scale_soft_parameters.cpp \
		$(DIR)/HSSUSY_two_scale_susy_scale_constraint.cpp
EXEHSSUSY_SRC += \
		$(DIR)/run_HSSUSY.cpp \
		$(DIR)/run_cmd_line_HSSUSY.cpp \
		$(DIR)/scan_HSSUSY.cpp
LIBHSSUSY_HDR += \
		$(DIR)/HSSUSY_convergence_tester.hpp \
		$(DIR)/HSSUSY_effective_couplings.hpp \
		$(DIR)/HSSUSY_high_scale_constraint.hpp \
		$(DIR)/HSSUSY_mass_eigenstates.hpp \
		$(DIR)/HSSUSY_info.hpp \
		$(DIR)/HSSUSY_initial_guesser.hpp \
		$(DIR)/HSSUSY_input_parameters.hpp \
		$(DIR)/HSSUSY_low_scale_constraint.hpp \
		$(DIR)/HSSUSY_model.hpp \
		$(DIR)/HSSUSY_model_slha.hpp \
		$(DIR)/HSSUSY_observables.hpp \
		$(DIR)/HSSUSY_physical.hpp \
		$(DIR)/HSSUSY_slha_io.hpp \
		$(DIR)/HSSUSY_spectrum_generator_interface.hpp \
		$(DIR)/HSSUSY_spectrum_generator.hpp \
		$(DIR)/HSSUSY_standard_model_matching.hpp \
		$(DIR)/HSSUSY_standard_model_two_scale_matching.hpp \
		$(DIR)/HSSUSY_susy_scale_constraint.hpp \
		$(DIR)/HSSUSY_utilities.hpp \
		$(DIR)/HSSUSY_two_scale_convergence_tester.hpp \
		$(DIR)/HSSUSY_two_scale_high_scale_constraint.hpp \
		$(DIR)/HSSUSY_two_scale_initial_guesser.hpp \
		$(DIR)/HSSUSY_two_scale_low_scale_constraint.hpp \
		$(DIR)/HSSUSY_two_scale_model.hpp \
		$(DIR)/HSSUSY_two_scale_model_slha.hpp \
		$(DIR)/HSSUSY_two_scale_soft_parameters.hpp \
		$(DIR)/HSSUSY_two_scale_susy_parameters.hpp \
		$(DIR)/HSSUSY_two_scale_susy_scale_constraint.hpp
LLHSSUSY_SRC  += \
		$(DIR)/HSSUSY_librarylink.cpp

LLHSSUSY_MMA  += \
		$(DIR)/HSSUSY_librarylink.m \
		$(DIR)/run_HSSUSY.m

ifneq ($(MAKECMDGOALS),showbuild)
ifneq ($(MAKECMDGOALS),tag)
ifneq ($(MAKECMDGOALS),release)
ifneq ($(MAKECMDGOALS),doc)
-include $(HSSUSY_TWO_SCALE_SUSY_MK)
-include $(HSSUSY_TWO_SCALE_SOFT_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(HSSUSY_TWO_SCALE_SUSY_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(HSSUSY_TWO_SCALE_SOFT_MK): run-metacode-$(MODNAME)
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
LIBHSSUSY_SRC := $(sort $(LIBHSSUSY_SRC))
EXEHSSUSY_SRC := $(sort $(EXEHSSUSY_SRC))

LIBHSSUSY_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBHSSUSY_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBHSSUSY_SRC)))

EXEHSSUSY_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEHSSUSY_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEHSSUSY_SRC)))

EXEHSSUSY_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXEHSSUSY_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXEHSSUSY_SRC)))

LIBHSSUSY_DEP := \
		$(LIBHSSUSY_OBJ:.o=.d)

EXEHSSUSY_DEP := \
		$(EXEHSSUSY_OBJ:.o=.d)

LLHSSUSY_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLHSSUSY_SRC)))

LLHSSUSY_OBJ  := $(LLHSSUSY_SRC:.cpp=.o)
LLHSSUSY_LIB  := $(LLHSSUSY_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBHSSUSY     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_HSSUSY := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_HSSUSY := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBHSSUSY) $(EXEHSSUSY_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(HSSUSY_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBHSSUSY_SRC) $(HSSUSY_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBHSSUSY_HDR) $(HSSUSY_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXEHSSUSY_SRC) $(HSSUSY_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLHSSUSY_SRC) $(HSSUSY_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLHSSUSY_MMA) $(HSSUSY_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(HSSUSY_MK) $(HSSUSY_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(HSSUSY_TWO_SCALE_MK) $(HSSUSY_INSTALL_DIR)
ifneq ($(HSSUSY_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(HSSUSY_SLHA_INPUT) $(HSSUSY_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(HSSUSY_GNUPLOT) $(HSSUSY_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBHSSUSY_DEP)
		-rm -f $(EXEHSSUSY_DEP)
		-rm -f $(LLHSSUSY_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBHSSUSY)
		-rm -f $(LLHSSUSY_LIB)

clean-$(MODNAME)-obj:
		-rm -f $(LIBHSSUSY_OBJ)
		-rm -f $(EXEHSSUSY_OBJ)
		-rm -f $(LLHSSUSY_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		-rm -f $(LIBHSSUSY_SRC)
		-rm -f $(LIBHSSUSY_HDR)
		-rm -f $(EXEHSSUSY_SRC)
		-rm -f $(LLHSSUSY_SRC)
		-rm -f $(LLHSSUSY_MMA)
		-rm -f $(METACODE_STAMP_HSSUSY)
		-rm -f $(HSSUSY_TWO_SCALE_MK)
		-rm -f $(HSSUSY_SLHA_INPUT)
		-rm -f $(HSSUSY_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		-rm -f $(EXEHSSUSY_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(HSSUSY_TARBALL) \
		$(LIBHSSUSY_SRC) $(LIBHSSUSY_HDR) \
		$(EXEHSSUSY_SRC) \
		$(HSSUSY_MK) $(HSSUSY_TWO_SCALE_MK) \
		$(HSSUSY_SLHA_INPUT) $(HSSUSY_GNUPLOT)

$(LIBHSSUSY_SRC) $(LIBHSSUSY_HDR) $(EXEHSSUSY_SRC) $(LLHSSUSY_SRC) $(LLHSSUSY_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_HSSUSY)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_HSSUSY): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_HSSUSY)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_HSSUSY)"
		@echo "Note: to regenerate HSSUSY source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_HSSUSY)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_HSSUSY):
		@true
endif

$(LIBHSSUSY_DEP) $(EXEHSSUSY_DEP) $(LLHSSUSY_DEP) $(LIBHSSUSY_OBJ) $(EXEHSSUSY_OBJ) $(LLHSSUSY_OBJ) $(LLHSSUSY_LIB): \
	CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBHSSUSY_DEP) $(EXEHSSUSY_DEP) $(LLHSSUSY_DEP) $(LIBHSSUSY_OBJ) $(EXEHSSUSY_OBJ) $(LLHSSUSY_OBJ) $(LLHSSUSY_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLHSSUSY_OBJ) $(LLHSSUSY_LIB): \
	CPPFLAGS += $(shell $(MATH_INC_PATHS) --math-cmd="$(MATH)" -I --librarylink --mathlink)

$(LIBHSSUSY): $(LIBHSSUSY_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBHSSUSY) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$^ $(ADDONLIBS)) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLHSSUSY_LIB): $(LLHSSUSY_OBJ) $(LIBHSSUSY) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$^) $(ADDONLIBS) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

ALLDEP += $(LIBHSSUSY_DEP) $(EXEHSSUSY_DEP)
ALLSRC += $(LIBHSSUSY_SRC) $(EXEHSSUSY_SRC)
ALLLIB += $(LIBHSSUSY)
ALLEXE += $(EXEHSSUSY_EXE)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLHSSUSY_DEP)
ALLSRC += $(LLHSSUSY_SRC)
ALLLL  += $(LLHSSUSY_LIB)
endif
