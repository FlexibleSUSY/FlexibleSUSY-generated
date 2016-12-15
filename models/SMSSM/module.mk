DIR          := models/SMSSM
MODNAME      := SMSSM
SARAH_MODEL  := SMSSM
WITH_$(MODNAME) := yes

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
		$(DIR)/LesHouches.in.SMSSM_generated \
		$(DIR)/LesHouches.in.SMSSM

SMSSM_GNUPLOT := \
		$(DIR)/SMSSM_plot_rgflow.gnuplot \
		$(DIR)/SMSSM_plot_spectrum.gnuplot

SMSSM_TARBALL := \
		$(MODNAME).tar.gz

LIBSMSSM_SRC :=
EXESMSSM_SRC :=
LLSMSSM_LIB  :=
LLSMSSM_OBJ  :=
LLSMSSM_SRC  :=
LLSMSSM_MMA  :=

LIBSMSSM_HDR :=

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBSMSSM_SRC += \
		$(DIR)/SMSSM_effective_couplings.cpp \
		$(DIR)/SMSSM_mass_eigenstates.cpp \
		$(DIR)/SMSSM_info.cpp \
		$(DIR)/SMSSM_input_parameters.cpp \
		$(DIR)/SMSSM_observables.cpp \
		$(DIR)/SMSSM_slha_io.cpp \
		$(DIR)/SMSSM_physical.cpp \
		$(DIR)/SMSSM_utilities.cpp \
		$(DIR)/SMSSM_standard_model_matching.cpp \
		$(DIR)/SMSSM_standard_model_two_scale_matching.cpp \
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
		$(DIR)/SMSSM_effective_couplings.hpp \
		$(DIR)/SMSSM_high_scale_constraint.hpp \
		$(DIR)/SMSSM_mass_eigenstates.hpp \
		$(DIR)/SMSSM_info.hpp \
		$(DIR)/SMSSM_initial_guesser.hpp \
		$(DIR)/SMSSM_input_parameters.hpp \
		$(DIR)/SMSSM_low_scale_constraint.hpp \
		$(DIR)/SMSSM_model.hpp \
		$(DIR)/SMSSM_model_slha.hpp \
		$(DIR)/SMSSM_observables.hpp \
		$(DIR)/SMSSM_physical.hpp \
		$(DIR)/SMSSM_slha_io.hpp \
		$(DIR)/SMSSM_spectrum_generator_interface.hpp \
		$(DIR)/SMSSM_spectrum_generator.hpp \
		$(DIR)/SMSSM_standard_model_matching.hpp \
		$(DIR)/SMSSM_standard_model_two_scale_matching.hpp \
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
LLSMSSM_SRC  += \
		$(DIR)/SMSSM_librarylink.cpp

LLSMSSM_MMA  += \
		$(DIR)/SMSSM_librarylink.m \
		$(DIR)/run_SMSSM.m

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

EXESMSSM_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXESMSSM_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXESMSSM_SRC)))

LIBSMSSM_DEP := \
		$(LIBSMSSM_OBJ:.o=.d)

EXESMSSM_DEP := \
		$(EXESMSSM_OBJ:.o=.d)

LLSMSSM_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLSMSSM_SRC)))

LLSMSSM_OBJ  := $(LLSMSSM_SRC:.cpp=.o)
LLSMSSM_LIB  := $(LLSMSSM_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBSMSSM     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_SMSSM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_SMSSM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBSMSSM) $(EXESMSSM_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(SMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBSMSSM_SRC) $(SMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBSMSSM_HDR) $(SMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXESMSSM_SRC) $(SMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLSMSSM_SRC) $(SMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLSMSSM_MMA) $(SMSSM_INSTALL_DIR)
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
		-rm -f $(LLSMSSM_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBSMSSM)
		-rm -f $(LLSMSSM_LIB)

clean-$(MODNAME)-obj:
		-rm -f $(LIBSMSSM_OBJ)
		-rm -f $(EXESMSSM_OBJ)
		-rm -f $(LLSMSSM_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		-rm -f $(LIBSMSSM_SRC)
		-rm -f $(LIBSMSSM_HDR)
		-rm -f $(EXESMSSM_SRC)
		-rm -f $(LLSMSSM_SRC)
		-rm -f $(LLSMSSM_MMA)
		-rm -f $(METACODE_STAMP_SMSSM)
		-rm -f $(SMSSM_TWO_SCALE_MK)
		-rm -f $(SMSSM_SLHA_INPUT)
		-rm -f $(SMSSM_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		-rm -f $(EXESMSSM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(SMSSM_TARBALL) \
		$(LIBSMSSM_SRC) $(LIBSMSSM_HDR) \
		$(EXESMSSM_SRC) \
		$(LLSMSSM_SRC) $(LLSMSSM_MMA) \
		$(SMSSM_MK) $(SMSSM_TWO_SCALE_MK) \
		$(SMSSM_SLHA_INPUT) $(SMSSM_GNUPLOT)

$(LIBSMSSM_SRC) $(LIBSMSSM_HDR) $(EXESMSSM_SRC) $(LLSMSSM_SRC) $(LLSMSSM_MMA) \
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

$(LIBSMSSM_DEP) $(EXESMSSM_DEP) $(LLSMSSM_DEP) $(LIBSMSSM_OBJ) $(EXESMSSM_OBJ) $(LLSMSSM_OBJ) $(LLSMSSM_LIB): \
	CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBSMSSM_DEP) $(EXESMSSM_DEP) $(LLSMSSM_DEP) $(LIBSMSSM_OBJ) $(EXESMSSM_OBJ) $(LLSMSSM_OBJ) $(LLSMSSM_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLSMSSM_OBJ) $(LLSMSSM_LIB): \
	CPPFLAGS += $(shell $(MATH_INC_PATHS) --math-cmd="$(MATH)" -I --librarylink --mathlink)

$(LIBSMSSM): $(LIBSMSSM_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBSMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$^ $(ADDONLIBS)) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLSMSSM_LIB): $(LLSMSSM_OBJ) $(LIBSMSSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$^) $(ADDONLIBS) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

ALLDEP += $(LIBSMSSM_DEP) $(EXESMSSM_DEP)
ALLSRC += $(LIBSMSSM_SRC) $(EXESMSSM_SRC)
ALLLIB += $(LIBSMSSM)
ALLEXE += $(EXESMSSM_EXE)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLSMSSM_DEP)
ALLSRC += $(LLSMSSM_SRC)
ALLLL  += $(LLSMSSM_LIB)
endif
