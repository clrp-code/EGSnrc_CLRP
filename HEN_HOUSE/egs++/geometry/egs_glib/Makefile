#*****************************************************************************
#
#  $Id: Makefile,v 1.1 2005/04/29 20:50:20 iwan Exp $
#
#*****************************************************************************

include $(EGS_CONFIG)
include $(SPEC_DIR)egspp.spec
include $(SPEC_DIR)egspp_$(my_machine).conf

DEFS = $(DEF1)

library = egs_glib
lib_files = egs_glib
link2_libs = egspp egs_ndgeometry egs_planes

.PHONY: ndgeom clean

include $(SPEC_DIR)egspp_libs.spec

$(make_depend)

$(DSO2)$(library).$(obje): ndgeom

$(ABS_DSO)$(libpre)$(library)$(libext): ndgeom

ndgeom:
	cd ..$(DSEP)egs_nd_geometry && make

clean:
	rm -f $(ABS_DSO)$(libpre)$(library)$(libext)
	rm -f $(DSO2)$(library).$(obje)
