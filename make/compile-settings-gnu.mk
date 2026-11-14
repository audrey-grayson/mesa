FLAGS_CODE_SANITY := \
  -Wall \
  -Wextra \
  -Wno-unused-label \
  -Wno-unused-parameter \
  -D_FORTIFY_SOURCE=2 \

# Position independent code. By default everything is built with -fPIC.
# With PIC_ONLY_IF_DYNAMIC=yes, -fPIC is only used when DYNAMIC=yes; otherwise
# -fPIE (cheaper, still linkable into PIE executables) is used instead.
ifeq ($(PIC_ONLY_IF_DYNAMIC),yes)
  ifeq ($(DYNAMIC),yes)
    FLAGS_PIC := -fPIC
  else
    FLAGS_PIC := -fPIE
  endif
else
  FLAGS_PIC := -fPIC
endif
ifeq ($(NO_SEMANTIC_INTERPOSITION),yes)
  FLAGS_PIC += -fno-semantic-interposition
endif

FFLAGS_FP_SANITY :=
ifeq ($(WITH_FPE_CHECKS),yes)
FFLAGS_FP_SANITY += -ffpe-trap=invalid,overflow,zero -finit-real=snan
endif
FFLAGS_FORTRAN_SANITY := -std=f2008 -ffree-line-length-none -ffixed-line-length-none -Wno-unused-dummy-argument -Wno-compare-reals -Wno-do-subscript
FFLAGS_TRAMPOLINE := -ftrampoline-impl=heap
FLAGS_REPRO := -ffp-contract=off
FFLAGS_PREPROCESSOR := -cpp

ifeq ($(PROFILE),release)
  FLAGS_OPT := -O2
  FLAGS_DEBUG :=
else ifeq ($(PROFILE),release-with-dbg-info)
  FLAGS_OPT := -O2
  FLAGS_DEBUG := -ggdb
else ifeq ($(PROFILE),debug)
  FLAGS_OPT := -Og
  FLAGS_DEBUG := -ggdb
else ifeq ($(PROFILE),release-native)
  # Not portable: the binaries only run on CPUs similar to the build machine.
  # -ffp-contract=off (FLAGS_REPRO) is kept, but -O3/-march=native can still
  # change floating-point results through vectorization and libm selection.
  FLAGS_OPT := -O3
  FLAGS_DEBUG := -ggdb
  FLAGS_ARCH := -march=native
else
  $(error Unknown or unset PROFILE)
endif

# Link-time optimization (opt-in). Flags that influence code generation must be
# repeated at link time, since that is when LTO generates the code.
ifeq ($(WITH_LTO),yes)
  FLAGS_LTO := -flto=auto
  FLAGS_LTO_LINK := -flto=auto $(FLAGS_OPT) $(FLAGS_DEBUG) $(FLAGS_ARCH) $(FLAGS_REPRO)
  AR_TOOL := gcc-ar
else
  FLAGS_LTO :=
  FLAGS_LTO_LINK :=
  AR_TOOL := ar
endif

ifeq ($(WITH_OPENMP),yes)
  FLAGS_OPENMP := -fopenmp
else
  FLAGS_OPENMP :=
endif

# no-sign-zero only affects output formatting
FFLAGS_COMPAT := -fno-sign-zero
FFLAGS_FREE := -ffree-form -fimplicit-none
FFLAGS_FIXED := -ffixed-form
FLAGS_DEPS := $(call pkg-config,--cflags,$(EXTERNAL_DEPENDS_ON) $(INTERNAL_DEPENDS_ON)) $(INCLUDE_DIRS)

ifeq ($(WITH_COVERAGE),yes)
  FLAGS_COVERAGE := --coverage
else
  FLAGS_COVERAGE :=
endif

FFLAGS_SHARED  := \
  $(FLAGS_CODE_SANITY) \
  $(FFLAGS_FP_SANITY) \
  $(FFLAGS_FORTRAN_SANITY) \
  $(FFLAGS_TRAMPOLINE) \
  $(FLAGS_REPRO) \
  $(FFLAGS_PREPROCESSOR) \
  $(FLAGS_OPT) \
  $(FLAGS_DEBUG) \
  $(FLAGS_OPENMP) \
  $(FFLAGS_COMPAT) \
  $(FLAGS_DEPS) \
  $(FLAGS_COVERAGE)
FFLAGS_FIXED := $(FFLAGS_SHARED) $(FFLAGS_FIXED) $(FFLAGS)
_FFLAGS := $(FFLAGS_SHARED) $(FFLAGS_FREE) $(FFLAGS)
_CFLAGS := $(FLAGS_CODE_SANITY) $(FLAGS_REPRO) $(FLAGS_OPT) $(FLAGS_ARCH) $(FLAGS_PIC) $(FLAGS_LTO) $(FLAGS_DEBUG) $(FLAGS_DEPS) $(FLAGS_OPENMP) $(FLAGS_COVERAGE) $(CFLAGS)

PREPROCESS := gfortran -cpp -E
FCOMPILE := gfortran $(_FFLAGS) -c
FCOMPILE_MODULE := gfortran $(_FFLAGS) -w -c -fsyntax-only
FCOMPILE_FIXED := gfortran $(FFLAGS_FIXED) -c
FCOMPILE_MODULE_FIXED:= gfortran $(FFLAGS_FIXED) -w -c -fsyntax-only
CCOMPILE := gcc $(_CFLAGS) -c
LIB_DEP_ARGS := $(call pkg-config, --libs,$(EXTERNAL_DEPENDS_ON)) $(call pkg-config, --libs --static,$(INTERNAL_DEPENDS_ON))
LIB_TOOL_STATIC := $(AR_TOOL) rcs
LIB_TOOL_DYNAMIC := gfortran -shared $(FLAGS_LTO_LINK) $(FLAGS_OPENMP) $(FLAGS_COVERAGE) -o
EXECUTABLE := gfortran $(FLAGS_LTO_LINK) $(FLAGS_OPENMP) $(FLAGS_COVERAGE)
