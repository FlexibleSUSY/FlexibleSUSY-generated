DIR          := models/E6SSMtower
MODNAME      := E6SSMtower
SARAH_MODEL  := E6SSM
WITH_$(MODNAME) := yes

E6SSMtower_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

E6SSMtower_MK     := \
		$(DIR)/module.mk

E6SSMtower_TWO_SCALE_SUSY_MK := \
		$(DIR)/two_scale_susy.mk

E6SSMtower_TWO_SCALE_SOFT_MK := \
		$(DIR)/two_scale_soft.mk

E6SSMtower_TWO_SCALE_MK := \
		$(E6SSMtower_TWO_SCALE_SUSY_MK) \
		$(E6SSMtower_TWO_SCALE_SOFT_MK)

E6SSMtower_SLHA_INPUT := \
		$(DIR)/LesHouches.in.E6SSMtower_generated \
		$(DIR)/LesHouches.in.E6SSMtower

E6SSMtower_GNUPLOT := \
		$(DIR)/E6SSMtower_plot_rgflow.gnuplot \
		$(DIR)/E6SSMtower_plot_spectrum.gnuplot

E6SSMtower_TARBALL := \
		$(MODNAME).tar.gz

LIBE6SSMtower_SRC :=
EXEE6SSMtower_SRC :=
LLE6SSMtower_LIB  :=
LLE6SSMtower_OBJ  :=
LLE6SSMtower_SRC  :=
LLE6SSMtower_MMA  :=

LIBE6SSMtower_HDR :=

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBE6SSMtower_SRC += \
		$(DIR)/E6SSMtower_effective_couplings.cpp \
		$(DIR)/E6SSMtower_mass_eigenstates.cpp \
		$(DIR)/E6SSMtower_info.cpp \
		$(DIR)/E6SSMtower_input_parameters.cpp \
		$(DIR)/E6SSMtower_observables.cpp \
		$(DIR)/E6SSMtower_slha_io.cpp \
		$(DIR)/E6SSMtower_physical.cpp \
		$(DIR)/E6SSMtower_utilities.cpp \
		$(DIR)/E6SSMtower_standard_model_matching.cpp \
		$(DIR)/E6SSMtower_standard_model_two_scale_matching.cpp \
		$(DIR)/E6SSMtower_two_scale_convergence_tester.cpp \
		$(DIR)/E6SSMtower_two_scale_high_scale_constraint.cpp \
		$(DIR)/E6SSMtower_two_scale_initial_guesser.cpp \
		$(DIR)/E6SSMtower_two_scale_low_scale_constraint.cpp \
		$(DIR)/E6SSMtower_two_scale_model.cpp \
		$(DIR)/E6SSMtower_two_scale_model_slha.cpp \
		$(DIR)/E6SSMtower_two_scale_susy_parameters.cpp \
		$(DIR)/E6SSMtower_two_scale_soft_parameters.cpp \
		$(DIR)/E6SSMtower_two_scale_susy_scale_constraint.cpp
EXEE6SSMtower_SRC += \
		$(DIR)/run_E6SSMtower.cpp \
		$(DIR)/run_cmd_line_E6SSMtower.cpp \
		$(DIR)/scan_E6SSMtower.cpp
LIBE6SSMtower_HDR += \
		$(DIR)/E6SSMtower_convergence_tester.hpp \
		$(DIR)/E6SSMtower_effective_couplings.hpp \
		$(DIR)/E6SSMtower_high_scale_constraint.hpp \
		$(DIR)/E6SSMtower_mass_eigenstates.hpp \
		$(DIR)/E6SSMtower_info.hpp \
		$(DIR)/E6SSMtower_initial_guesser.hpp \
		$(DIR)/E6SSMtower_input_parameters.hpp \
		$(DIR)/E6SSMtower_low_scale_constraint.hpp \
		$(DIR)/E6SSMtower_model.hpp \
		$(DIR)/E6SSMtower_model_slha.hpp \
		$(DIR)/E6SSMtower_observables.hpp \
		$(DIR)/E6SSMtower_physical.hpp \
		$(DIR)/E6SSMtower_slha_io.hpp \
		$(DIR)/E6SSMtower_spectrum_generator_interface.hpp \
		$(DIR)/E6SSMtower_spectrum_generator.hpp \
		$(DIR)/E6SSMtower_standard_model_matching.hpp \
		$(DIR)/E6SSMtower_standard_model_two_scale_matching.hpp \
		$(DIR)/E6SSMtower_susy_scale_constraint.hpp \
		$(DIR)/E6SSMtower_utilities.hpp \
		$(DIR)/E6SSMtower_two_scale_convergence_tester.hpp \
		$(DIR)/E6SSMtower_two_scale_high_scale_constraint.hpp \
		$(DIR)/E6SSMtower_two_scale_initial_guesser.hpp \
		$(DIR)/E6SSMtower_two_scale_low_scale_constraint.hpp \
		$(DIR)/E6SSMtower_two_scale_model.hpp \
		$(DIR)/E6SSMtower_two_scale_model_slha.hpp \
		$(DIR)/E6SSMtower_two_scale_soft_parameters.hpp \
		$(DIR)/E6SSMtower_two_scale_susy_parameters.hpp \
		$(DIR)/E6SSMtower_two_scale_susy_scale_constraint.hpp
LLE6SSMtower_SRC  += \
		$(DIR)/E6SSMtower_librarylink.cpp

LLE6SSMtower_MMA  += \
		$(DIR)/E6SSMtower_librarylink.m \
		$(DIR)/run_E6SSMtower.m

ifneq ($(MAKECMDGOALS),showbuild)
ifneq ($(MAKECMDGOALS),tag)
ifneq ($(MAKECMDGOALS),release)
ifneq ($(MAKECMDGOALS),doc)
-include $(E6SSMtower_TWO_SCALE_SUSY_MK)
-include $(E6SSMtower_TWO_SCALE_SOFT_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(E6SSMtower_TWO_SCALE_SUSY_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(E6SSMtower_TWO_SCALE_SOFT_MK): run-metacode-$(MODNAME)
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
LIBE6SSMtower_SRC := $(sort $(LIBE6SSMtower_SRC))
EXEE6SSMtower_SRC := $(sort $(EXEE6SSMtower_SRC))

LIBE6SSMtower_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBE6SSMtower_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBE6SSMtower_SRC)))

EXEE6SSMtower_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEE6SSMtower_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEE6SSMtower_SRC)))

EXEE6SSMtower_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXEE6SSMtower_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXEE6SSMtower_SRC)))

LIBE6SSMtower_DEP := \
		$(LIBE6SSMtower_OBJ:.o=.d)

EXEE6SSMtower_DEP := \
		$(EXEE6SSMtower_OBJ:.o=.d)

LLE6SSMtower_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLE6SSMtower_SRC)))

LLE6SSMtower_OBJ  := $(LLE6SSMtower_SRC:.cpp=.o)
LLE6SSMtower_LIB  := $(LLE6SSMtower_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBE6SSMtower     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_E6SSMtower := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_E6SSMtower := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBE6SSMtower) $(EXEE6SSMtower_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(E6SSMtower_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBE6SSMtower_SRC) $(E6SSMtower_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBE6SSMtower_HDR) $(E6SSMtower_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXEE6SSMtower_SRC) $(E6SSMtower_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLE6SSMtower_SRC) $(E6SSMtower_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLE6SSMtower_MMA) $(E6SSMtower_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(E6SSMtower_MK) $(E6SSMtower_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(E6SSMtower_TWO_SCALE_MK) $(E6SSMtower_INSTALL_DIR)
ifneq ($(E6SSMtower_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(E6SSMtower_SLHA_INPUT) $(E6SSMtower_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(E6SSMtower_GNUPLOT) $(E6SSMtower_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBE6SSMtower_DEP)
		-rm -f $(EXEE6SSMtower_DEP)
		-rm -f $(LLE6SSMtower_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBE6SSMtower)
		-rm -f $(LLE6SSMtower_LIB)

clean-$(MODNAME)-obj:
		-rm -f $(LIBE6SSMtower_OBJ)
		-rm -f $(EXEE6SSMtower_OBJ)
		-rm -f $(LLE6SSMtower_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		-rm -f $(LIBE6SSMtower_SRC)
		-rm -f $(LIBE6SSMtower_HDR)
		-rm -f $(EXEE6SSMtower_SRC)
		-rm -f $(LLE6SSMtower_SRC)
		-rm -f $(LLE6SSMtower_MMA)
		-rm -f $(METACODE_STAMP_E6SSMtower)
		-rm -f $(E6SSMtower_TWO_SCALE_MK)
		-rm -f $(E6SSMtower_SLHA_INPUT)
		-rm -f $(E6SSMtower_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		-rm -f $(EXEE6SSMtower_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(E6SSMtower_TARBALL) \
		$(LIBE6SSMtower_SRC) $(LIBE6SSMtower_HDR) \
		$(EXEE6SSMtower_SRC) \
		$(E6SSMtower_MK) $(E6SSMtower_TWO_SCALE_MK) \
		$(E6SSMtower_SLHA_INPUT) $(E6SSMtower_GNUPLOT)

$(LIBE6SSMtower_SRC) $(LIBE6SSMtower_HDR) $(EXEE6SSMtower_SRC) $(LLE6SSMtower_SRC) $(LLE6SSMtower_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_E6SSMtower)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_E6SSMtower): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_E6SSMtower)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_E6SSMtower)"
		@echo "Note: to regenerate E6SSMtower source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_E6SSMtower)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_E6SSMtower):
		@true
endif

$(LIBE6SSMtower_DEP) $(EXEE6SSMtower_DEP) $(LLE6SSMtower_DEP) $(LIBE6SSMtower_OBJ) $(EXEE6SSMtower_OBJ) $(LLE6SSMtower_OBJ) $(LLE6SSMtower_LIB): \
	CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBE6SSMtower_DEP) $(EXEE6SSMtower_DEP) $(LLE6SSMtower_DEP) $(LIBE6SSMtower_OBJ) $(EXEE6SSMtower_OBJ) $(LLE6SSMtower_OBJ) $(LLE6SSMtower_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLE6SSMtower_OBJ) $(LLE6SSMtower_LIB): \
	CPPFLAGS += $(shell $(MATH_INC_PATHS) --math-cmd="$(MATH)" -I --librarylink --mathlink)

$(LIBE6SSMtower): $(LIBE6SSMtower_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBE6SSMtower) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$^ $(ADDONLIBS)) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLE6SSMtower_LIB): $(LLE6SSMtower_OBJ) $(LIBE6SSMtower) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$^) $(ADDONLIBS) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

ALLDEP += $(LIBE6SSMtower_DEP) $(EXEE6SSMtower_DEP)
ALLSRC += $(LIBE6SSMtower_SRC) $(EXEE6SSMtower_SRC)
ALLLIB += $(LIBE6SSMtower)
ALLEXE += $(EXEE6SSMtower_EXE)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLE6SSMtower_DEP)
ALLSRC += $(LLE6SSMtower_SRC)
ALLLL  += $(LLE6SSMtower_LIB)
endif
