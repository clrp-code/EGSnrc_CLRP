#*****************************************************************************
#
#  $Id: Makefile,v 1.1 2005/04/29 20:50:20 iwan Exp $
#
#*****************************************************************************

include $(EGS_CONFIG)
include $(SPEC_DIR)egspp.spec
include $(SPEC_DIR)egspp_$(my_machine).conf

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

DEFS = $(DEF1) $(GZSTREAM_DEF)

library = egs_glib
lib_files = egs_glib $(GZSTREAM)
link2_libs = egspp egs_ndgeometry egs_planes $(GZSTREAM_LIB)

.PHONY: ndgeom clean

include $(SPEC_DIR)egspp_libs.spec

$(make_depend)

$(DSO2)$(library).$(obje): ndgeom $(GZSTREAM_H)

$(ABS_DSO)$(libpre)$(library)$(libext): ndgeom $(GZSTREAM_H)

ndgeom:
	cd ..$(DSEP)egs_nd_geometry && make

clean:
	rm -f $(ABS_DSO)$(libpre)$(library)$(libext)
	rm -f $(DSO2)$(library).$(obje)
