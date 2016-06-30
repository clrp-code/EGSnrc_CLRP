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
ifneq ("$(wildcard gzstream.cpp)","")
GZSTREAM = gzstream
GZSTREAM_DEF = -DBUILD_GZSTREAM -DHAS_GZSTREAM
GZSTREAM_H = gzstream.h
endif

DEFS = $(DEF1) $(GZSTREAM_DEF)

library = egs_glib
lib_files = egs_glib $(GZSTREAM)

link2_libs = egspp egs_ndgeometry egs_planes

include $(SPEC_DIR)egspp_libs.spec

$(DSO2)glib.$(obje): $(GZSTREAM_H)

$(ABS_DSO)$(libpre)egs_glib$(libext): $(GZSTREAM_H)

$(make_depend)


clean:
	rm -f $(ABS_DSO)$(libpre)egs_glib$(libext)
	rm -f $(DSO2)egs_glib.$(obje)

