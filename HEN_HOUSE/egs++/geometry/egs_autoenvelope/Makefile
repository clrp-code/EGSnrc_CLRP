#*****************************************************************************
#
#  $Id: Makefile,v 1.1 2005/04/29 20:50:19 iwan Exp $
#
#*****************************************************************************

include $(EGS_CONFIG)
include $(SPEC_DIR)egspp.spec
include $(SPEC_DIR)egspp_$(my_machine).conf

SOBOL =
SOBOL_DEF =
SOBOL_H =
ifneq ("$(wildcard egs_sobol.cpp)","")
SOBOL = sobol egs_sobol
SOBOL_DEF = -DBUILD_SOBOL_DLL -DBUILD_AE_SOBOL -DHAS_SOBOL
SOBOL_H = sobol.h egs_sobol.h
endif

GZSTREAM =
GZSTREAM_DEF =
GZSTREAM_H =
GZSTREAM_LIB =
ifneq ("$(wildcard gzstream.cpp)","")
GZSTREAM = gzstream
GZSTREAM_DEF = -DBUILD_GZSTREAM -DHAS_GZSTREAM
GZSTREAM_H = gzstream.h
GZSTREAM_LIB = z
endif

DEFS = $(DEF1) -DBUILD_AENVELOPE_DLL $(SOBOL_DEF) $(GZSTREAM_DEF)


library = egs_autoenvelope
lib_files = egs_autoenvelope $(SOBOL) $(GZSTREAM)

link2_libs = egs_gtransformed egspp $(GZSTREAM_LIB)


include $(SPEC_DIR)egspp_libs.spec

$(make_depend)


$(DSO2)autoenvelope.$(obje): volcor.h $(GZSTREAM_H) \
                          ..$(DSEP)egs_gtransformed$(DSEP)egs_gtransformed.h

$(ABS_DSO)$(libpre)egs_autoenvelope$(libext): volcor.h $(GZSTREAM_H) \
						..$(DSEP)egs_gtransformed$(DSEP)egs_gtransformed.h

clean:
	rm -f $(ABS_DSO)$(libpre)egs_autoenvelope$(libext) $(ABS_DSO)$(libpre)egs_sobol$(libext) $(ABS_DSO)$(libpre)sobol$(libext) \
		$(ABS_DSO)$(libpre)gzstream$(libext) $(ABS_DSO)$(libpre)gzstream$(libext)
	rm -f $(DSO2)autoenvelope.$(obje) $(DSO2)sobol.$(obje) $(DSO2)egs_sobol.$(obje) \
		$(DSO2)gzstream.$(obje) $(DSO2)gzstream.$(obje) \

touchup:
	touch egs_autoenvelope.cpp

fresh: clean touchup $(lib_objects)


.PHONY: clean touchup fresh
