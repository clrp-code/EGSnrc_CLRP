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

extra_dep = ndgeom

include $(SPEC_DIR)egspp_libs.spec

$(make_depend)

ndgeom:
	cd ..$(DSEP)egs_nd_geometry$(DSEP) && make

$(DSO1)egs_spectra.$(obje): $(config1h) \
    egs_rndm.h egs_object_factory.h egs_alias_table.h egs_vector.h \
    egs_input.h egs_math.h

$(DSO2)brachy_seeds.$(obje): ..$(DSEP)egs_planes$(DSEP)egs_planes.h \
                         ..$(DSEP)egs_cylinders$(DSEP)egs_cylinders.h \
                         ..$(DSEP)egs_nd_geometry$(DSEP)egs_nd_geometry.h \
