#*****************************************************************************
#
#  $Id: Makefile,v 1.1 2005/04/29 20:50:20 iwan Exp $
#
#*****************************************************************************

include $(EGS_CONFIG)
include $(SPEC_DIR)egspp.spec
include $(SPEC_DIR)egspp_$(my_machine).conf

DEFS = $(DEF1) -Wall

library = egs_rz
lib_files = egs_rz
link2_libs = egspp egs_planes egs_cylinders egs_ndgeometry

.PHONY: ndgeom cyls planes clean

include $(SPEC_DIR)egspp_libs.spec

$(make_depend)

$(DSO2)$(library).$(obje): ndgeom cyls planes

$(ABS_DSO)$(libpre)$(library)$(libext): ndgeom cyls planes

ndgeom:
	cd ..$(DSEP)egs_nd_geometry && make

cyls:
	cd ..$(DSEP)egs_cylinders && make

planes:
	cd ..$(DSEP)egs_planes && make

clean:
	rm -f $(ABS_DSO)$(libpre)$(library)$(libext)
	rm -f $(DSO2)$(library).$(obje)
