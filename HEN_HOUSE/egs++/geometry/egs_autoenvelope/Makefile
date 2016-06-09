#*****************************************************************************
#
#  $Id: Makefile,v 1.1 2005/04/29 20:50:19 iwan Exp $
#
#*****************************************************************************

include $(EGS_CONFIG)
include $(SPEC_DIR)egspp.spec
include $(SPEC_DIR)egspp_$(my_machine).conf

DEFS = $(DEF1) -DBUILD_AENVELOPE_DLL -DBUILD_SOBOL_DLL -DBUILD_AE_SOBOL

library = egs_autoenvelope
lib_files = egs_autoenvelope sobol egs_sobol gzstream

link2_libs = egs_gtransformed egspp


include $(SPEC_DIR)egspp_libs.spec

$(make_depend)


$(DSO2)autoenvelope.$(obje): volcor.h gzstream.h \
                          ..$(DSEP)egs_gtransformed$(DSEP)egs_gtransformed.h

$(ABS_DSO)$(libpre)egs_autoenvelope$(libext): volcor.h gzstream.h \
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
