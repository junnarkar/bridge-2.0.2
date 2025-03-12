# For GNU MAKE

### User modification part ########################### START ###

#####  OPTIONS  ################################################

##--------------------------------------------------------------
## 1. Machine environment and build conditions

## 1-a. Target platform
target = PC_GNU       # PC using GNU C++
#target = PC_CLANG     # PC using Clang C++
#target = PC_INTEL     # PC using intel compiler
#target = PC_NVIDIA    # PC using NVIDIA HPC SDK compiler
#target = Fugaku_CLANG # Fugaku using clang compiler
#target = SX_AURORA    # NEC SX-Aurora at KEK

## 1-b. Communication
## use_mpi = {yes, no, fjmpi}
use_mpi = no

## 1-c. Multi-threading
## use_thread = {yes, omp, no}, NB. yes == omp
use_thread = no

## 1-d. Optimized code
## use_opt_code = {yes, no}
use_opt_code = yes

## 1-e. Compiler optimization level
## use_opt_level = {low, std, high, trace}
##                 ('trace' is available for SX-Aurora)
use_opt_level = std

## 1-f. Static linkage
## use_static = {yes, no}
##   "use_static=yes" needs static system libraries.
use_static = no

## 1-g. Debug option
## use_debug = {yes, no}
use_debug = no

## 1-h. Complex type
## use_complex = {std, c99}
##   c99 option is compatible with gcc version 4.4.x
##   CAUTION: ver.2.0 requires std::complex and some codes fails for c99
use_complex = std

## use_cpp11 = {yes, no}
##   set "use_cpp11=yes", unless your compiler does not support cpp11.
##   CAUTION: ver.2.0 assumes cpp11 is available.
use_cpp11 = yes

## 1-i. python interface
## use_python = {yes, no}
##   set "yes" if you use python interface in sub-package.
##   CAUTION: in ver.2.0 python interface is not available.
use_python = no

## 1-j. Alternative code
## use_alternative = {yes, no}
##   alternative implementation that enables multi-precision and accelerator.
use_alternative = no

## use_alt_qxs = {yes, no}
use_alt_qxs = no

## use_qxs_arch = {acle, general}
use_qxs_arch = acle

## 1-k. Build directory
build_base_dir = ./build

#ifeq ($(use_mpi),yes)
#  build_base_dir:=$(build_base_dir)_mpi
#endif

build_postfix =
ifeq ($(use_opt_level),low)
  build_postfix:=$(build_postfix)_low
endif
ifeq ($(use_opt_level),high)
  build_postfix:=$(build_postfix)_high
endif
ifeq ($(use_opt_level),trace)
  build_postfix:=$(build_postfix)_trace
endif
ifeq ($(use_debug),yes)
  build_postfix:=$(build_postfix)_debug
endif

#build_dir = $(build_base_dir)$(build_postfix)
build_dir = $(build_base_dir)

test_dir = $(build_dir)/tests


##--------------------------------------------------------------
## 2. Simulation settings

## 2-a. Gauge Group
## use_gauge_group = {su3, su2, general}
use_gauge_group = su3

## 2-b. Target test name
test = all
#test = none
#test = Eigensolver
#test = FFT
#test = Gauge
#test = GradientFlow
#test = HMC
#test = HotStart
#test = IO
#test = Mult
#test = PolyakovLoop
#test = QuarkNumberSusceptibility
#test = RandomNumbers
#test = Rational
#test = SF_fAfP
#test = Solver
#test = Spectrum
#test = TopologicalCharge
#test = WilsonLoop

## 2-c. Use Test Manager
## use_testmanager = {yes, no}
use_testmanager = yes

##--------------------------------------------------------------
## 3. External libraries

extra_root_dir = extra

## 3-a. LIME
## use_lime_library = {yes, no}
use_lime_library = no

lime_library_path = $(extra_root_dir)/lime-1.3.2
lime_library_path_includes = $(lime_library_path)/include
lime_library_path_libs = $(lime_library_path)/lib
lime_library_libs = -llime

## 3-b. SFMT (Mersenne Twister)
## use_sfmt_library = {yes, no}
use_sfmt_library = no

sfmt_library_path = $(extra_root_dir)/SFMT-1.5.1_bridge
sfmt_library_path_includes = $(sfmt_library_path)
sfmt_library_path_libs = $(sfmt_library_path)
sfmt_library_libs = -lsfmt

## 3-c. NTL (used with SFMT)
## use_ntl_library = {yes, no}
use_ntl_library = no

ntl_library_path = $(extra_root_dir)/ntl
ntl_library_path_includes = $(ntl_library_path)/include
ntl_library_path_libs = $(ntl_library_path)/lib
ntl_library_libs = -lntl -lgmp

## 3-d. FFTW
## use_fftw_library = {yes, no}
use_fftw_library = no
ifeq ($(test),FFT)
  use_fftw_library = yes
endif

fftw_library_path = $(extra_root_dir)/fftw
fftw_library_path_includes = $(fftw_library_path)/include
fftw_library_path_libs = $(fftw_library_path)/lib
ifeq ($(use_mpi),yes)
  fftw_library_libs += -lfftw3_mpi
endif
ifeq ($(use_mpi),fjmpi)
  fftw_library_libs += -lfftw3_mpi
endif
ifeq ($(use_thread),omp)
  fftw_library_libs += -lfftw3_omp
endif
fftw_library_libs += -lfftw3

## 3-e. yaml-cpp
## use_yamlcpp_library = {yes, no}
use_yamlcpp_library = no

yamlcpp_library_path = $(extra_root_dir)/yaml-cpp
yamlcpp_library_path_includes = ${yamlcpp_library_path}/include
yamlcpp_library_path_libs = ${yamlcpp_library_path}/lib
yamlcpp_library_libs += -lyaml-cpp

## 3-f. tinyxml2
## use_tinyxml2_library = {yes, no}
use_tinyxml2_library = no

tinyxml2_library_path = $(extra_root_dir)/tinyxml2
tinyxml2_library_path_includes = ${tinyxml2_library_path}/include
tinyxml2_library_path_libs = ${tinyxml2_library_path}/lib
tinyxml2_library_libs += -ltinyxml2

## 3-g. QWS
## use_qws_library = {yes, no}
use_qws_library = no

qws_library_path = $(extra_root_dir)/qws-modified
qws_library_path_includes = $(qws_library_path)
qws_library_path_libs = $(qws_library_path)
qws_library_libs = -lqws


################################################################

## check compiler option conflicts
#ifeq ($(use_thread),omp)
#  ifeq ($(use_opt_code),no)
#    $(error Incompatible options [use_thread=omp,use_opt_code=no])
#  endif
#endif

ifeq ($(test),none)
  ifeq ($(use_testmanager),yes)
    $(warning Conflicting options [test=none,use_testmanager=yes]: No test is registered on TestManager)
  endif
endif

## replace "use_thread = yes" with "use_thread = omp"
ifeq ($(use_thread),yes)
  use_thread = omp
endif

ifeq ($(use_complex),c99)
  $(error Unavailable use_complex=c99: in ver.2.0, some classes do not work correctly.)
endif

include Makefile_target.inc

### User modification part ############################# END ###


### System part ###################################### START ###

## flag setup (default).
SRCCOMMUN = Single
SRCOPT  = Imp
SRCOPT2 = Imp
EXTRA_INCLUDES=
EXTRA_LIBS=
CHECK_EXTRA_LIBS=

CPP_DEFS += -D$(strip $(target))

ifeq ($(use_lime_library),yes)
  EXTRA_INCLUDES   += -I$(lime_library_path_includes)
  EXTRA_LIBS       += -L$(lime_library_path_libs) $(lime_library_libs)
  CHECK_EXTRA_LIBS += check-liblime
  CPP_DEFS         += -DUSE_LIMELIB
else ifeq ($(use_lime_library),no)
else
  $(error use_lime_library=$(use_lime_library), which must be {yes, no})
endif

ifeq ($(use_sfmt_library),yes)
  EXTRA_INCLUDES   += -I$(sfmt_library_path_includes)
  EXTRA_LIBS       += -L$(sfmt_library_path_libs) $(sfmt_library_libs)
  CHECK_EXTRA_LIBS += check-libsfmt
  CPP_DEFS         += -DUSE_SFMTLIB
else ifeq ($(use_sfmt_library),no)
else
  $(error use_sfmt_library=$(use_sfmt_library), which must be {yes, no})
endif

ifeq ($(use_ntl_library),yes)
  EXTRA_INCLUDES   += -I$(ntl_library_path_includes)
  EXTRA_LIBS       += -L$(ntl_library_path_libs) $(ntl_library_libs)
  CHECK_EXTRA_LIBS += check-libntl
  CPP_DEFS         += -DUSE_NTLLIB
else ifeq ($(use_ntl_library),no)
else
  $(error use_ntl_library=$(use_ntl_library), which must be {yes, no})
endif

ifeq ($(use_fftw_library),yes)
  EXTRA_INCLUDES   += -I$(fftw_library_path_includes)
  EXTRA_LIBS       += -L$(fftw_library_path_libs) $(fftw_library_libs)
  CHECK_EXTRA_LIBS += check-libfftw
  CPP_DEFS         += -DUSE_FFTWLIB
else ifeq ($(use_fftw_library),no)
else
  $(error use_fftw_library=$(use_fftw_library), which must be {yes, no})
endif

ifeq ($(use_yamlcpp_library),yes)
  EXTRA_INCLUDES   += -I$(yamlcpp_library_path_includes)
  EXTRA_LIBS       += -L$(yamlcpp_library_path_libs) $(yamlcpp_library_libs)
  CHECK_EXTRA_LIBS += check-libyamlcpp
  CPP_DEFS         += -DUSE_YAMLCPPLIB
else ifeq ($(use_yamlcpp_library),no)
else
  $(error use_yamlcpp_library=$(use_yamlcpp_library), which must be {yes, no})
endif

ifeq ($(use_tinyxml2_library),yes)
  EXTRA_INCLUDES   += -I$(tinyxml2_library_path_includes)
  EXTRA_LIBS       += -L$(tinyxml2_library_path_libs) $(tinyxml2_library_libs)
  CHECK_EXTRA_LIBS += check-libtinyxml2
  CPP_DEFS         += -DUSE_XML -DUSE_TINYXML2LIB
else ifeq ($(use_tinyxml2_library),no)
else
  $(error use_tinyxml2_library=$(use_tinyxml2_library), which must be {yes, no})
endif

ifeq ($(use_qws_library),yes)
  EXTRA_INCLUDES   += -I$(qws_library_path_includes)
  EXTRA_LIBS       += -L$(qws_library_path_libs) $(qws_library_libs)
  CHECK_EXTRA_LIBS += check-libqws
  CPP_DEFS         += -DUSE_QWSLIB
else ifeq ($(use_qws_library),no)
else
  $(error use_qws_library=$(use_qws_library), which must be {yes, no})
endif

ifeq ($(use_mpi),yes)
  SRCCOMMUN = MPI
  CPP_DEFS += -DUSE_MPI
else ifeq ($(use_mpi),fjmpi)
  SRCCOMMUN = FJMPI
  CPP_DEFS += -DUSE_MPI -DUSE_FJMPI
else ifeq ($(use_mpi),no)
  SRCCOMMUN = Single
else
  $(error use_mpi=$(use_mpi), which must be {yes, no, fjmpi})
endif

ifeq ($(use_opt_code),yes)
  SRCOPT  = Imp
  SRCOPT2 = Imp
  CPP_DEFS += -DUSE_IMP
else ifeq ($(use_opt_code),no)
#  ifeq ($(use_thread),no)
    SRCOPT  = Org
    SRCOPT2 = Org
    CPP_DEFS += -DUSE_ORG
#  else
#    $(error Incompatible options: Org version does not support multithreading)
#  endif
else
  $(error use_opt_code=$(use_opt_code), which must be {yes, no})
endif

ifeq ($(use_gauge_group),su2)
  CPP_DEFS += -DUSE_GROUP_SU2
else ifeq ($(use_gauge_group),su3)
  CPP_DEFS += -DUSE_GROUP_SU3
else ifeq ($(use_gauge_group),general)
  CPP_DEFS += -DUSE_GROUP_SU_N
else
  $(error use_gauge_group=$(use_gauge_group), which must be {su3, su2, general})
endif

ifeq ($(use_thread),omp)
  CPP_DEFS += -DUSE_OPENMP
  SRCTHREAD = OpenMP
else ifeq ($(use_thread),no)
  SRCTHREAD = OpenMP_stub
else
  $(error use_thread=$(use_thread), which must be {omp, no})
endif

ifeq ($(use_complex),std)
  CPP_DEFS += -DUSE_STD_COMPLEX
else ifeq ($(use_complex),c99)
  CPP_DEFS += -DUSE_C99_COMPLEX
else
  $(error use_complex=$(use_complex), which must be {std, c99})
endif

ifeq ($(use_cpp11),yes)
  CPP_DEFS += -DLIB_CPP11
else ifeq ($(use_cpp11),no)
else
  $(error use_cpp11=$(use_cpp11), which must be {yes, no})
endif

CPP_DEFS += -DUSE_FACTORY # -DUSE_FACTORY_AUTOREGISTER


## apply only to tests
ifeq ($(use_testmanager),yes)
  CPP_DEFS_TEST += -DUSE_TESTMANAGER -DUSE_TESTMANAGER_AUTOREGISTER
else ifeq ($(use_testmanager),no)
else
  $(error use_testmanager=$(use_testmanager), which must be {yes, no})
endif


# for alt-code
alt_inc_dirs =
ifeq ($(use_alternative),yes)
 alt_inc_dirs += -Iinclude/bridge/lib_alt
 CPP_DEFS += -DUSE_ALT_CODE
ifeq ($(use_alt_qxs),yes)
 alt_inc_dirs += -Iinclude/bridge/lib_alt_QXS
 CPP_DEFS += -DUSE_ALT_QXS
endif
endif


##--

CPPFLAGS =
CPPFLAGS += $(EXTRA_INCLUDES)
CPPFLAGS += $(CPP_DEFS)

CPPFLAGS_TEST = $(CPP_DEFS_TEST)

##--


build_library_base_dir = $(build_dir)
build_test_base_dir = $(test_dir)

msg_file = $(build_library_base_dir)/message.txt
make_inc_file = $(build_library_base_dir)/make.inc

##--

.PHONY: all msg check-extra_libs all-done
.PHONY: program program-with-lib
.PHONY: lib make-inc
.PHONY: make-build_dir make-test_dir

## export variables to subsequent make
export

## make all

#all: msg check-extra_libs lib program-with-lib all-done
#all: msg check-extra_libs link-inline lib program-with-lib all-done
all: msg lib program-with-lib all-done

## build library
lib: msg check-extra_libs link-inline make-inc
	$(MAKE) -f src/Makefile build_dir=$(build_dir) lib

## build test application
program: msg make-test_dir
	build_dir_path=`cd $(build_dir) && pwd`; \
	test_dir_path=`cd $(test_dir) && pwd`; \
	cd tests && $(MAKE) bridge_install_path=$$build_dir_path test_dir=$$test_dir_path test="""$(test)""" test_flags="""$(CPPFLAGS_TEST)"""

program-with-lib: msg lib make-test_dir
	build_dir_path=`cd $(build_dir) && pwd`; \
	test_dir_path=`cd $(test_dir) && pwd`; \
	cd tests && $(MAKE) bridge_install_path=$$build_dir_path test_dir=$$test_dir_path test="""$(test)""" test_flags="""$(CPPFLAGS_TEST)"""

## make work directories
make-build_dir:
	@if [ ! -d $(build_dir) ]; then \
		mkdir -p $(build_dir); \
		echo "create directory: $(build_dir)"; \
	fi

make-test_dir:
	@if [ ! -d $(test_dir) ]; then \
		mkdir -p $(test_dir); \
		echo "create directory: $(test_dir)"; \
	fi

## show summary of build options
msg: make-build_dir
	@( \
	echo "----------------------------------------------------------------"; \
	echo "Compilation Environment Summary" ; \
	echo ; \
	echo "Options:" ; \
	echo "  target              =" $(target) "(" $(TARGET_MACHINE) ")" ; \
	echo "  use_mpi             =" $(use_mpi) ; \
	echo "  use_thread          =" $(use_thread) ; \
	echo "  use_opt_code        =" $(use_opt_code) ; \
	echo "  use_opt_level       =" $(use_opt_level) ; \
	echo "  use_static          =" $(use_static) ; \
	echo "  use_debug           =" $(use_debug) ; \
	echo "  use_gauge_group     =" $(use_gauge_group) ; \
	echo "  use_complex         =" $(use_complex) ; \
	echo "  use_cpp11           =" $(use_cpp11) ; \
	echo "  use_testmanager     =" $(use_testmanager) ; \
	echo "  use_lime_library    =" $(use_lime_library) ; \
	echo "  use_sfmt_library    =" $(use_sfmt_library) ; \
	echo "  use_ntl_library     =" $(use_ntl_library) ; \
	echo "  use_fftw_library    =" $(use_fftw_library) ; \
	echo "  use_yamlcpp_library =" $(use_yamlcpp_library) ; \
	echo "  use_tinyxml2_library=" $(use_tinyxml2_library) ; \
	echo "  use_qws_library     =" $(use_qws_library) ; \
	echo ; \
	echo "Compiler and options:" ; \
	echo "  CXX         =" $(CXX) ; \
	echo "  CXXFLAGS    =" $(CXXFLAGS) ; \
	echo "  LD          =" $(LD) ; \
	echo "  LDFLAGS     =" $(LDFLAGS) ; \
	echo "  LDLIBS      =" $(LDLIBS) ; \
	echo "  CXXDEP      =" $(CXXDEP) ; \
	echo "  CXXDEPFLAGS =" $(CXXDEPFLAGS) ; \
	echo "  CPPFLAGS    =" $(CPPFLAGS) ; \
	echo "----------------------------------------------------------------" ; \
	) | tee $(msg_file)


make-inc: make-build_dir
	@echo generating make.inc file ...
	@( \
	echo "# make.inc" ; \
	echo "# See $(msg_file) for Compilation Environment Summary." ; \
	echo "# This file is automatically created. Do not edit this file manually."; \
	echo ; \
	echo "BRIDGE_ROOT =" `cd $(build_library_base_dir) && pwd` ; \
	echo "BRIDGE_EXTRA_PATH =" `cd $(extra_root_dir) && pwd` ; \
	echo ; \
	echo "CXX         =" $(CXX) ; \
	echo "CXXFLAGS    =" $(CXXFLAGS) ; \
	echo "LD          =" $(LD) ; \
	echo "LDFLAGS     =" $(LDFLAGS) ; \
	echo "LDLIBS      =" $(LDLIBS) ; \
	echo "CXXDEP      =" $(CXXDEP) ; \
	echo "CXXDEPFLAGS =" $(CXXDEPFLAGS) ; \
	echo ; \
	echo "CPPFLAGS    =" $(CPP_DEFS) ; \
	echo "CPPFLAGS_TEST  =" $(CPP_DEFS_TEST) ; \
	echo ; \
	echo 'BRIDGE_INCLUDE = -I$$(BRIDGE_ROOT)/include/bridge -I$$(BRIDGE_ROOT)/include/bridge/lib' ; \
        echo "BRIDGE_ALT_INCLUDE =" `echo $(alt_inc_dirs) | sed -e 's,I,I$$(BRIDGE_ROOT)/,g'` ; \
	echo 'BRIDGE_LIB     = -L$$(BRIDGE_ROOT) -lbridge' ; \
	echo ; \
	echo "EXTRA_INCLUDES =" `echo $(EXTRA_INCLUDES) | sed -e 's,$(extra_root_dir),$$(BRIDGE_EXTRA_PATH),g'` ; \
	echo "EXTRA_LIBS     =" `echo $(EXTRA_LIBS) | sed -e 's,$(extra_root_dir),$$(BRIDGE_EXTRA_PATH),g'` ; \
	echo ; \
	echo 'CXXFLAGS   += $$(CPPFLAGS) $$(BRIDGE_INCLUDE) $$(BRIDGE_ALT_INCLUDE) $$(EXTRA_INCLUDES)' ; \
	echo 'LDLIBS     += $$(BRIDGE_LIB) $$(EXTRA_LIBS)' ; \
	) > $(make_inc_file)


## check extra/
check-extra_libs: $(CHECK_EXTRA_LIBS)

check-liblime:
	@if [ ! -e $(lime_library_path_libs)/liblime.a ]; then \
		echo "Error: liblime not found. Maybe extra/lime is not installed."; \
		exit 1; \
	fi

check-libsfmt:
	@if [ ! -e $(sfmt_library_path_libs)/libsfmt.a ]; then \
		echo "Error: libsfmt not found. Maybe extra/sfmt is not installed."; \
		exit 1; \
	fi

check-libntl:
	@if [ ! -e $(ntl_library_path_libs)/libntl.a ]; then \
		echo "Error: libntl not found. Maybe extra/ntl is not installed."; \
		exit 1; \
	fi

check-libfftw:
	@if [ ! -e $(fftw_library_path_libs)/libfftw3.a ]; then \
		echo "Error: libfftw not found. Maybe extra/fftw is not installed."; \
		exit 1; \
	fi

check-libyamlcpp:
	@if [ ! -e $(yamlcpp_library_path_libs)/libyaml-cpp.a ]; then \
		echo "Error: libyaml-cpp not found. Maybe extra/yaml-cpp is not installed."; \
		exit 1; \
	fi

check-libtinyxml2:
	@if [ ! -e $(tinyxml2_library_path_libs)/libtinyxml2.a ]; then \
		echo "Error: tinyxml2 not found. Maybe extra/tinyxml2 is not installed."; \
		exit 1; \
	fi

check-libqws:
	@if [ ! -e $(qws_library_path_libs)/libqws.a ]; then \
		echo "Error: libqws not found. Maybe extra/qws-modified is not installed."; \
		exit 1; \
	fi


## final messages

all-done: lib program-with-lib
	@echo ""
	@echo "Bridge code compilation is successful"
	@echo ""


# symbolic link of inline directory
.PHONY: link-inline
link-inline:
ifeq ($(use_alt_qxs),yes)
#ifeq ($(use_qws_library),yes)
#	- cp --no-clobber ./extra/qws-master/qws.h ./src/lib_alt_QXS/extra/
#	- cp --no-clobber ./extra/qws-master/qws_xbound.h ./src/lib_alt_QXS/extra/
#else
#	- cd ./src/lib_alt_QXS/extra ; ln -s qws_dummy.h qws.h ; cd ../../../
#endif
ifeq ($(use_qxs_arch),acle)
	cd ./src/lib_alt_QXS/; ln -s inline_ACLE inline ; cd ../../
endif
ifeq ($(use_qxs_arch),general)
	cd ./src/lib_alt_QXS/; ln -s inline_General inline ; cd ../../
endif
endif


# symbolic link of inline directory
.PHONY: delete_links
delete_links:
ifeq ($(use_alt_qxs),yes)
	rm ./src/lib_alt_QXS/inline
endif


## generate doxygen html files to ./docs/html
.PHONY: doxygen
doxygen:
	@if [ -e ./docs/doxygen.conf ]; then \
	  echo "generating doxygen documents ..."; \
	  ( \
	    cd ./docs; \
	    doxygen doxygen.conf \
	  ) \
	fi

## make package
VERSION=$(shell cat VERSION)
NAME=bridge
PACK_DIR=$(NAME)-$(VERSION)

package_contents=\
  ChangeLog.txt \
  INSTALL.txt \
  LICENSE.txt \
  Makefile \
  Makefile.in \
  Makefile_target.inc \
  README.txt \
  VERSION \
  docs \
  extra \
  makeconfig.sh \
  sample_hmc \
  sample_spectrum \
  test_alt_spectrum \
  test_alt_multigrid \
  src \
  tests
#  python \

.PHONY: package
package:
	cd .. && mkdir $(PACK_DIR)
	cp -rf $(package_contents) ../$(PACK_DIR)
	find ../$(PACK_DIR) -name ".svn" | xargs rm -rf
	cd .. && tar -zcvf $(addsuffix .tar.gz, $(PACK_DIR)) $(PACK_DIR) && rm -rf $(PACK_DIR)

## file clean up
.PHONY: clean

clean:
	if [ -d $(build_dir) ]; then \
		build_dir_path=`cd $(build_dir) && pwd`; \
		test_dir_path=`cd $(test_dir) && pwd`; \
		cd tests && $(MAKE) bridge_install_path=$$build_dir_path test_dir=$$test_dir_path test=$(test) test_flags="""$(CPPFLAGS_TEST)""" clean; \
	fi
	$(MAKE) -f src/Makefile clean
	$(MAKE) -f src/Makefile real-clean
	$(MAKE) delete_links
#	
#	

################################################################################
#        @author  Shinji Motoki (smotoki) 
#                 $LastChangedBy: namekawa $
#        @date    $LastChangedDate:: 2023-04-09 20:20:34 #$
#        @version $LastChangedRevision: 2509 $
################################################################################
