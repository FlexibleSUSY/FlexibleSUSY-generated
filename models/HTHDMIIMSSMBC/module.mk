DIR          := models/HTHDMIIMSSMBC
MODNAME      := HTHDMIIMSSMBC
SARAH_MODEL  := HTHDM-II
WITH_$(MODNAME) := yes

HTHDMIIMSSMBC_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

HTHDMIIMSSMBC_MK     := \
		$(DIR)/module.mk

HTHDMIIMSSMBC_TWO_SCALE_SUSY_MK := \
		$(DIR)/two_scale_susy.mk

HTHDMIIMSSMBC_TWO_SCALE_SOFT_MK := \
		$(DIR)/two_scale_soft.mk

HTHDMIIMSSMBC_TWO_SCALE_MK := \
		$(HTHDMIIMSSMBC_TWO_SCALE_SUSY_MK) \
		$(HTHDMIIMSSMBC_TWO_SCALE_SOFT_MK)

HTHDMIIMSSMBC_SLHA_INPUT := \
		$(DIR)/LesHouches.in.HTHDMIIMSSMBC_generated \
		$(DIR)/LesHouches.in.HTHDMIIMSSMBC

HTHDMIIMSSMBC_GNUPLOT := \
		$(DIR)/HTHDMIIMSSMBC_plot_rgflow.gnuplot \
		$(DIR)/HTHDMIIMSSMBC_plot_spectrum.gnuplot

HTHDMIIMSSMBC_TARBALL := \
		$(MODNAME).tar.gz

LIBHTHDMIIMSSMBC_SRC :=
EXEHTHDMIIMSSMBC_SRC :=
LLHTHDMIIMSSMBC_LIB  :=
LLHTHDMIIMSSMBC_OBJ  :=
LLHTHDMIIMSSMBC_SRC  :=
LLHTHDMIIMSSMBC_MMA  :=

LIBHTHDMIIMSSMBC_HDR :=

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBHTHDMIIMSSMBC_SRC += \
		$(DIR)/HTHDMIIMSSMBC_effective_couplings.cpp \
		$(DIR)/HTHDMIIMSSMBC_mass_eigenstates.cpp \
		$(DIR)/HTHDMIIMSSMBC_info.cpp \
		$(DIR)/HTHDMIIMSSMBC_input_parameters.cpp \
		$(DIR)/HTHDMIIMSSMBC_observables.cpp \
		$(DIR)/HTHDMIIMSSMBC_slha_io.cpp \
		$(DIR)/HTHDMIIMSSMBC_physical.cpp \
		$(DIR)/HTHDMIIMSSMBC_utilities.cpp \
		$(DIR)/HTHDMIIMSSMBC_standard_model_matching.cpp \
		$(DIR)/HTHDMIIMSSMBC_standard_model_two_scale_matching.cpp \
		$(DIR)/HTHDMIIMSSMBC_two_scale_convergence_tester.cpp \
		$(DIR)/HTHDMIIMSSMBC_two_scale_high_scale_constraint.cpp \
		$(DIR)/HTHDMIIMSSMBC_two_scale_initial_guesser.cpp \
		$(DIR)/HTHDMIIMSSMBC_two_scale_low_scale_constraint.cpp \
		$(DIR)/HTHDMIIMSSMBC_two_scale_model.cpp \
		$(DIR)/HTHDMIIMSSMBC_two_scale_model_slha.cpp \
		$(DIR)/HTHDMIIMSSMBC_two_scale_susy_parameters.cpp \
		$(DIR)/HTHDMIIMSSMBC_two_scale_soft_parameters.cpp \
		$(DIR)/HTHDMIIMSSMBC_two_scale_susy_scale_constraint.cpp
EXEHTHDMIIMSSMBC_SRC += \
		$(DIR)/run_HTHDMIIMSSMBC.cpp \
		$(DIR)/run_cmd_line_HTHDMIIMSSMBC.cpp \
		$(DIR)/scan_HTHDMIIMSSMBC.cpp
LIBHTHDMIIMSSMBC_HDR += \
		$(DIR)/HTHDMIIMSSMBC_convergence_tester.hpp \
		$(DIR)/HTHDMIIMSSMBC_effective_couplings.hpp \
		$(DIR)/HTHDMIIMSSMBC_high_scale_constraint.hpp \
		$(DIR)/HTHDMIIMSSMBC_mass_eigenstates.hpp \
		$(DIR)/HTHDMIIMSSMBC_info.hpp \
		$(DIR)/HTHDMIIMSSMBC_initial_guesser.hpp \
		$(DIR)/HTHDMIIMSSMBC_input_parameters.hpp \
		$(DIR)/HTHDMIIMSSMBC_low_scale_constraint.hpp \
		$(DIR)/HTHDMIIMSSMBC_model.hpp \
		$(DIR)/HTHDMIIMSSMBC_model_slha.hpp \
		$(DIR)/HTHDMIIMSSMBC_observables.hpp \
		$(DIR)/HTHDMIIMSSMBC_physical.hpp \
		$(DIR)/HTHDMIIMSSMBC_slha_io.hpp \
		$(DIR)/HTHDMIIMSSMBC_spectrum_generator_interface.hpp \
		$(DIR)/HTHDMIIMSSMBC_spectrum_generator.hpp \
		$(DIR)/HTHDMIIMSSMBC_standard_model_matching.hpp \
		$(DIR)/HTHDMIIMSSMBC_standard_model_two_scale_matching.hpp \
		$(DIR)/HTHDMIIMSSMBC_susy_scale_constraint.hpp \
		$(DIR)/HTHDMIIMSSMBC_utilities.hpp \
		$(DIR)/HTHDMIIMSSMBC_two_scale_convergence_tester.hpp \
		$(DIR)/HTHDMIIMSSMBC_two_scale_high_scale_constraint.hpp \
		$(DIR)/HTHDMIIMSSMBC_two_scale_initial_guesser.hpp \
		$(DIR)/HTHDMIIMSSMBC_two_scale_low_scale_constraint.hpp \
		$(DIR)/HTHDMIIMSSMBC_two_scale_model.hpp \
		$(DIR)/HTHDMIIMSSMBC_two_scale_model_slha.hpp \
		$(DIR)/HTHDMIIMSSMBC_two_scale_soft_parameters.hpp \
		$(DIR)/HTHDMIIMSSMBC_two_scale_susy_parameters.hpp \
		$(DIR)/HTHDMIIMSSMBC_two_scale_susy_scale_constraint.hpp
LLHTHDMIIMSSMBC_SRC  += \
		$(DIR)/HTHDMIIMSSMBC_librarylink.cpp

LLHTHDMIIMSSMBC_MMA  += \
		$(DIR)/HTHDMIIMSSMBC_librarylink.m \
		$(DIR)/run_HTHDMIIMSSMBC.m

ifneq ($(MAKECMDGOALS),showbuild)
ifneq ($(MAKECMDGOALS),tag)
ifneq ($(MAKECMDGOALS),release)
ifneq ($(MAKECMDGOALS),doc)
-include $(HTHDMIIMSSMBC_TWO_SCALE_SUSY_MK)
-include $(HTHDMIIMSSMBC_TWO_SCALE_SOFT_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(HTHDMIIMSSMBC_TWO_SCALE_SUSY_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(HTHDMIIMSSMBC_TWO_SCALE_SOFT_MK): run-metacode-$(MODNAME)
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
LIBHTHDMIIMSSMBC_SRC := $(sort $(LIBHTHDMIIMSSMBC_SRC))
EXEHTHDMIIMSSMBC_SRC := $(sort $(EXEHTHDMIIMSSMBC_SRC))

LIBHTHDMIIMSSMBC_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBHTHDMIIMSSMBC_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBHTHDMIIMSSMBC_SRC)))

EXEHTHDMIIMSSMBC_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEHTHDMIIMSSMBC_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEHTHDMIIMSSMBC_SRC)))

EXEHTHDMIIMSSMBC_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXEHTHDMIIMSSMBC_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXEHTHDMIIMSSMBC_SRC)))

LIBHTHDMIIMSSMBC_DEP := \
		$(LIBHTHDMIIMSSMBC_OBJ:.o=.d)

EXEHTHDMIIMSSMBC_DEP := \
		$(EXEHTHDMIIMSSMBC_OBJ:.o=.d)

LLHTHDMIIMSSMBC_DEP  := \
		$(patsubst %.cpp, %.d, $(filter %.cpp, $(LLHTHDMIIMSSMBC_SRC)))

LLHTHDMIIMSSMBC_OBJ  := $(LLHTHDMIIMSSMBC_SRC:.cpp=.o)
LLHTHDMIIMSSMBC_LIB  := $(LLHTHDMIIMSSMBC_SRC:.cpp=$(LIBLNK_LIBEXT))

LIBHTHDMIIMSSMBC     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

METACODE_STAMP_HTHDMIIMSSMBC := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_HTHDMIIMSSMBC := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBHTHDMIIMSSMBC) $(EXEHTHDMIIMSSMBC_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(HTHDMIIMSSMBC_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBHTHDMIIMSSMBC_SRC) $(HTHDMIIMSSMBC_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBHTHDMIIMSSMBC_HDR) $(HTHDMIIMSSMBC_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXEHTHDMIIMSSMBC_SRC) $(HTHDMIIMSSMBC_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLHTHDMIIMSSMBC_SRC) $(HTHDMIIMSSMBC_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LLHTHDMIIMSSMBC_MMA) $(HTHDMIIMSSMBC_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(HTHDMIIMSSMBC_MK) $(HTHDMIIMSSMBC_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(HTHDMIIMSSMBC_TWO_SCALE_MK) $(HTHDMIIMSSMBC_INSTALL_DIR)
ifneq ($(HTHDMIIMSSMBC_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(HTHDMIIMSSMBC_SLHA_INPUT) $(HTHDMIIMSSMBC_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(HTHDMIIMSSMBC_GNUPLOT) $(HTHDMIIMSSMBC_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBHTHDMIIMSSMBC_DEP)
		-rm -f $(EXEHTHDMIIMSSMBC_DEP)
		-rm -f $(LLHTHDMIIMSSMBC_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBHTHDMIIMSSMBC)
		-rm -f $(LLHTHDMIIMSSMBC_LIB)

clean-$(MODNAME)-obj:
		-rm -f $(LIBHTHDMIIMSSMBC_OBJ)
		-rm -f $(EXEHTHDMIIMSSMBC_OBJ)
		-rm -f $(LLHTHDMIIMSSMBC_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		-rm -f $(LIBHTHDMIIMSSMBC_SRC)
		-rm -f $(LIBHTHDMIIMSSMBC_HDR)
		-rm -f $(EXEHTHDMIIMSSMBC_SRC)
		-rm -f $(LLHTHDMIIMSSMBC_SRC)
		-rm -f $(LLHTHDMIIMSSMBC_MMA)
		-rm -f $(METACODE_STAMP_HTHDMIIMSSMBC)
		-rm -f $(HTHDMIIMSSMBC_TWO_SCALE_MK)
		-rm -f $(HTHDMIIMSSMBC_SLHA_INPUT)
		-rm -f $(HTHDMIIMSSMBC_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		-rm -f $(EXEHTHDMIIMSSMBC_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(HTHDMIIMSSMBC_TARBALL) \
		$(LIBHTHDMIIMSSMBC_SRC) $(LIBHTHDMIIMSSMBC_HDR) \
		$(EXEHTHDMIIMSSMBC_SRC) \
		$(LLHTHDMIIMSSMBC_SRC) $(LLHTHDMIIMSSMBC_MMA) \
		$(HTHDMIIMSSMBC_MK) $(HTHDMIIMSSMBC_TWO_SCALE_MK) \
		$(HTHDMIIMSSMBC_SLHA_INPUT) $(HTHDMIIMSSMBC_GNUPLOT)

$(LIBHTHDMIIMSSMBC_SRC) $(LIBHTHDMIIMSSMBC_HDR) $(EXEHTHDMIIMSSMBC_SRC) $(LLHTHDMIIMSSMBC_SRC) $(LLHTHDMIIMSSMBC_MMA) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_HTHDMIIMSSMBC)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_HTHDMIIMSSMBC): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_HTHDMIIMSSMBC)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_HTHDMIIMSSMBC)"
		@echo "Note: to regenerate HTHDMIIMSSMBC source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_HTHDMIIMSSMBC)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_HTHDMIIMSSMBC):
		@true
endif

$(LIBHTHDMIIMSSMBC_DEP) $(EXEHTHDMIIMSSMBC_DEP) $(LLHTHDMIIMSSMBC_DEP) $(LIBHTHDMIIMSSMBC_OBJ) $(EXEHTHDMIIMSSMBC_OBJ) $(LLHTHDMIIMSSMBC_OBJ) $(LLHTHDMIIMSSMBC_LIB): \
	CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBHTHDMIIMSSMBC_DEP) $(EXEHTHDMIIMSSMBC_DEP) $(LLHTHDMIIMSSMBC_DEP) $(LIBHTHDMIIMSSMBC_OBJ) $(EXEHTHDMIIMSSMBC_OBJ) $(LLHTHDMIIMSSMBC_OBJ) $(LLHTHDMIIMSSMBC_LIB): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LLHTHDMIIMSSMBC_OBJ) $(LLHTHDMIIMSSMBC_LIB): \
	CPPFLAGS += $(shell $(MATH_INC_PATHS) --math-cmd="$(MATH)" -I --librarylink --mathlink)

$(LIBHTHDMIIMSSMBC): $(LIBHTHDMIIMSSMBC_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBHTHDMIIMSSMBC) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$^ $(ADDONLIBS)) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

$(LLHTHDMIIMSSMBC_LIB): $(LLHTHDMIIMSSMBC_OBJ) $(LIBHTHDMIIMSSMBC) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(LIBLNK_MAKE_LIB_CMD) $@ $(CPPFLAGS) $(CFLAGS) $(call abspathx,$^) $(ADDONLIBS) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS) $(LDLIBS)

ALLDEP += $(LIBHTHDMIIMSSMBC_DEP) $(EXEHTHDMIIMSSMBC_DEP)
ALLSRC += $(LIBHTHDMIIMSSMBC_SRC) $(EXEHTHDMIIMSSMBC_SRC)
ALLLIB += $(LIBHTHDMIIMSSMBC)
ALLEXE += $(EXEHTHDMIIMSSMBC_EXE)

ifeq ($(ENABLE_LIBRARYLINK),yes)
ALLDEP += $(LLHTHDMIIMSSMBC_DEP)
ALLSRC += $(LLHTHDMIIMSSMBC_SRC)
ALLLL  += $(LLHTHDMIIMSSMBC_LIB)
endif
