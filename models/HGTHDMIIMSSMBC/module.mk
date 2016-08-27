DIR          := models/HGTHDMIIMSSMBC
MODNAME      := HGTHDMIIMSSMBC
SARAH_MODEL  := HGTHDM-II
WITH_$(MODNAME) := yes

HGTHDMIIMSSMBC_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

HGTHDMIIMSSMBC_MK     := \
		$(DIR)/module.mk

HGTHDMIIMSSMBC_TWO_SCALE_SUSY_MK := \
		$(DIR)/two_scale_susy.mk

HGTHDMIIMSSMBC_TWO_SCALE_SOFT_MK := \
		$(DIR)/two_scale_soft.mk

HGTHDMIIMSSMBC_TWO_SCALE_MK := \
		$(HGTHDMIIMSSMBC_TWO_SCALE_SUSY_MK) \
		$(HGTHDMIIMSSMBC_TWO_SCALE_SOFT_MK)

HGTHDMIIMSSMBC_SLHA_INPUT := \
		$(DIR)/LesHouches.in.HGTHDMIIMSSMBC_generated \
		$(DIR)/LesHouches.in.HGTHDMIIMSSMBC

HGTHDMIIMSSMBC_GNUPLOT := \
		$(DIR)/HGTHDMIIMSSMBC_plot_rgflow.gnuplot \
		$(DIR)/HGTHDMIIMSSMBC_plot_spectrum.gnuplot

HGTHDMIIMSSMBC_TARBALL := \
		$(MODNAME).tar.gz

LIBHGTHDMIIMSSMBC_SRC :=
EXEHGTHDMIIMSSMBC_SRC :=

LIBHGTHDMIIMSSMBC_HDR :=

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBHGTHDMIIMSSMBC_SRC += \
		$(DIR)/HGTHDMIIMSSMBC_effective_couplings.cpp \
		$(DIR)/HGTHDMIIMSSMBC_mass_eigenstates.cpp \
		$(DIR)/HGTHDMIIMSSMBC_info.cpp \
		$(DIR)/HGTHDMIIMSSMBC_input_parameters.cpp \
		$(DIR)/HGTHDMIIMSSMBC_observables.cpp \
		$(DIR)/HGTHDMIIMSSMBC_slha_io.cpp \
		$(DIR)/HGTHDMIIMSSMBC_physical.cpp \
		$(DIR)/HGTHDMIIMSSMBC_utilities.cpp \
		$(DIR)/HGTHDMIIMSSMBC_two_scale_convergence_tester.cpp \
		$(DIR)/HGTHDMIIMSSMBC_two_scale_high_scale_constraint.cpp \
		$(DIR)/HGTHDMIIMSSMBC_two_scale_initial_guesser.cpp \
		$(DIR)/HGTHDMIIMSSMBC_two_scale_low_scale_constraint.cpp \
		$(DIR)/HGTHDMIIMSSMBC_two_scale_model.cpp \
		$(DIR)/HGTHDMIIMSSMBC_two_scale_model_slha.cpp \
		$(DIR)/HGTHDMIIMSSMBC_two_scale_susy_parameters.cpp \
		$(DIR)/HGTHDMIIMSSMBC_two_scale_soft_parameters.cpp \
		$(DIR)/HGTHDMIIMSSMBC_two_scale_susy_scale_constraint.cpp
EXEHGTHDMIIMSSMBC_SRC += \
		$(DIR)/run_HGTHDMIIMSSMBC.cpp \
		$(DIR)/run_cmd_line_HGTHDMIIMSSMBC.cpp \
		$(DIR)/scan_HGTHDMIIMSSMBC.cpp
LIBHGTHDMIIMSSMBC_HDR += \
		$(DIR)/HGTHDMIIMSSMBC_convergence_tester.hpp \
		$(DIR)/HGTHDMIIMSSMBC_effective_couplings.hpp \
		$(DIR)/HGTHDMIIMSSMBC_high_scale_constraint.hpp \
		$(DIR)/HGTHDMIIMSSMBC_mass_eigenstates.hpp \
		$(DIR)/HGTHDMIIMSSMBC_info.hpp \
		$(DIR)/HGTHDMIIMSSMBC_initial_guesser.hpp \
		$(DIR)/HGTHDMIIMSSMBC_input_parameters.hpp \
		$(DIR)/HGTHDMIIMSSMBC_low_scale_constraint.hpp \
		$(DIR)/HGTHDMIIMSSMBC_model.hpp \
		$(DIR)/HGTHDMIIMSSMBC_model_slha.hpp \
		$(DIR)/HGTHDMIIMSSMBC_observables.hpp \
		$(DIR)/HGTHDMIIMSSMBC_physical.hpp \
		$(DIR)/HGTHDMIIMSSMBC_slha_io.hpp \
		$(DIR)/HGTHDMIIMSSMBC_spectrum_generator_interface.hpp \
		$(DIR)/HGTHDMIIMSSMBC_spectrum_generator.hpp \
		$(DIR)/HGTHDMIIMSSMBC_susy_scale_constraint.hpp \
		$(DIR)/HGTHDMIIMSSMBC_utilities.hpp \
		$(DIR)/HGTHDMIIMSSMBC_two_scale_convergence_tester.hpp \
		$(DIR)/HGTHDMIIMSSMBC_two_scale_high_scale_constraint.hpp \
		$(DIR)/HGTHDMIIMSSMBC_two_scale_initial_guesser.hpp \
		$(DIR)/HGTHDMIIMSSMBC_two_scale_low_scale_constraint.hpp \
		$(DIR)/HGTHDMIIMSSMBC_two_scale_model.hpp \
		$(DIR)/HGTHDMIIMSSMBC_two_scale_model_slha.hpp \
		$(DIR)/HGTHDMIIMSSMBC_two_scale_soft_parameters.hpp \
		$(DIR)/HGTHDMIIMSSMBC_two_scale_susy_parameters.hpp \
		$(DIR)/HGTHDMIIMSSMBC_two_scale_susy_scale_constraint.hpp

ifneq ($(MAKECMDGOALS),showbuild)
ifneq ($(MAKECMDGOALS),tag)
ifneq ($(MAKECMDGOALS),release)
ifneq ($(MAKECMDGOALS),doc)
-include $(HGTHDMIIMSSMBC_TWO_SCALE_SUSY_MK)
-include $(HGTHDMIIMSSMBC_TWO_SCALE_SOFT_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(HGTHDMIIMSSMBC_TWO_SCALE_SUSY_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(HGTHDMIIMSSMBC_TWO_SCALE_SOFT_MK): run-metacode-$(MODNAME)
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
LIBHGTHDMIIMSSMBC_SRC := $(sort $(LIBHGTHDMIIMSSMBC_SRC))
EXEHGTHDMIIMSSMBC_SRC := $(sort $(EXEHGTHDMIIMSSMBC_SRC))

LIBHGTHDMIIMSSMBC_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBHGTHDMIIMSSMBC_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBHGTHDMIIMSSMBC_SRC)))

EXEHGTHDMIIMSSMBC_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEHGTHDMIIMSSMBC_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEHGTHDMIIMSSMBC_SRC)))

EXEHGTHDMIIMSSMBC_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXEHGTHDMIIMSSMBC_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXEHGTHDMIIMSSMBC_SRC)))

LIBHGTHDMIIMSSMBC_DEP := \
		$(LIBHGTHDMIIMSSMBC_OBJ:.o=.d)

EXEHGTHDMIIMSSMBC_DEP := \
		$(EXEHGTHDMIIMSSMBC_OBJ:.o=.d)

LIBHGTHDMIIMSSMBC     := $(DIR)/lib$(MODNAME)$(LIBEXT)

METACODE_STAMP_HGTHDMIIMSSMBC := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_HGTHDMIIMSSMBC := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib \
		clean-$(MODNAME)-obj distclean-$(MODNAME) \
		run-metacode-$(MODNAME) pack-$(MODNAME)-src

all-$(MODNAME): $(LIBHGTHDMIIMSSMBC) $(EXEHGTHDMIIMSSMBC_EXE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(HGTHDMIIMSSMBC_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBHGTHDMIIMSSMBC_SRC) $(HGTHDMIIMSSMBC_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBHGTHDMIIMSSMBC_HDR) $(HGTHDMIIMSSMBC_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXEHGTHDMIIMSSMBC_SRC) $(HGTHDMIIMSSMBC_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(HGTHDMIIMSSMBC_MK) $(HGTHDMIIMSSMBC_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(HGTHDMIIMSSMBC_TWO_SCALE_MK) $(HGTHDMIIMSSMBC_INSTALL_DIR)
ifneq ($(HGTHDMIIMSSMBC_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(HGTHDMIIMSSMBC_SLHA_INPUT) $(HGTHDMIIMSSMBC_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(HGTHDMIIMSSMBC_GNUPLOT) $(HGTHDMIIMSSMBC_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBHGTHDMIIMSSMBC_DEP)
		-rm -f $(EXEHGTHDMIIMSSMBC_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBHGTHDMIIMSSMBC)

clean-$(MODNAME)-obj:
		-rm -f $(LIBHGTHDMIIMSSMBC_OBJ)
		-rm -f $(EXEHGTHDMIIMSSMBC_OBJ)

# BEGIN: NOT EXPORTED ##########################################
clean-$(MODNAME)-src:
		-rm -f $(LIBHGTHDMIIMSSMBC_SRC)
		-rm -f $(LIBHGTHDMIIMSSMBC_HDR)
		-rm -f $(EXEHGTHDMIIMSSMBC_SRC)
		-rm -f $(METACODE_STAMP_HGTHDMIIMSSMBC)
		-rm -f $(HGTHDMIIMSSMBC_TWO_SCALE_MK)
		-rm -f $(HGTHDMIIMSSMBC_SLHA_INPUT)
		-rm -f $(HGTHDMIIMSSMBC_GNUPLOT)

distclean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		-rm -f $(EXEHGTHDMIIMSSMBC_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		@true

clean-generated:: clean-$(MODNAME)-src

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(HGTHDMIIMSSMBC_TARBALL) \
		$(LIBHGTHDMIIMSSMBC_SRC) $(LIBHGTHDMIIMSSMBC_HDR) \
		$(EXEHGTHDMIIMSSMBC_SRC) \
		$(HGTHDMIIMSSMBC_MK) $(HGTHDMIIMSSMBC_TWO_SCALE_MK) \
		$(HGTHDMIIMSSMBC_SLHA_INPUT) $(HGTHDMIIMSSMBC_GNUPLOT)

$(LIBHGTHDMIIMSSMBC_SRC) $(LIBHGTHDMIIMSSMBC_HDR) $(EXEHGTHDMIIMSSMBC_SRC) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_HGTHDMIIMSSMBC)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_HGTHDMIIMSSMBC): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_HGTHDMIIMSSMBC)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_HGTHDMIIMSSMBC)"
		@echo "Note: to regenerate HGTHDMIIMSSMBC source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_HGTHDMIIMSSMBC)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_HGTHDMIIMSSMBC):
		@true
endif

$(LIBHGTHDMIIMSSMBC_DEP) $(EXEHGTHDMIIMSSMBC_DEP) $(LIBHGTHDMIIMSSMBC_OBJ) $(EXEHGTHDMIIMSSMBC_OBJ): CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBHGTHDMIIMSSMBC_DEP) $(EXEHGTHDMIIMSSMBC_DEP) $(LIBHGTHDMIIMSSMBC_OBJ) $(EXEHGTHDMIIMSSMBC_OBJ): CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LIBHGTHDMIIMSSMBC): $(LIBHGTHDMIIMSSMBC_OBJ)
		$(MAKELIB) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBHGTHDMIIMSSMBC) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) $(LDFLAGS) -o $@ $(call abspathx,$^ $(ADDONLIBS)) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(LDLIBS)

ALLDEP += $(LIBHGTHDMIIMSSMBC_DEP) $(EXEHGTHDMIIMSSMBC_DEP)
ALLSRC += $(LIBHGTHDMIIMSSMBC_SRC) $(EXEHGTHDMIIMSSMBC_SRC)
ALLLIB += $(LIBHGTHDMIIMSSMBC)
ALLEXE += $(EXEHGTHDMIIMSSMBC_EXE)
