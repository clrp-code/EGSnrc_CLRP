################################################################################
#
#  egs_mird Makefile
#  Copyright (C) 2021 Ting Lee, Rowan Thomson, Martin Martinov
#
#  This file is made for use in egs_mird
#
#  egs_mird is free software: you can redistribute it and/or modify it
#  under the terms of the GNU Affero General Public License as published
#  by the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  egs_mird is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Affero General Public License for more details:
#  <http://www.gnu.org/licenses/>.
#
################################################################################
#
#  When egs_mird is used for publications, please cite our paper:
#  
#  			TO BE ANNOUNCED
#  
#
################################################################################
#
#  Author:        Martin Martinov, 2021
#
#  Contributors:  
#
################################################################################

include $(EGS_CONFIG)
include $(SPEC_DIR)egspp.spec
include $(SPEC_DIR)egspp_$(my_machine).conf

DEFS = $(DEF1) -DBUILD_INTERNAL_SOURCE_DLL

library = egs_internal_source
lib_files = egs_internal_source

my_deps = $(common_source_deps)
extra_dep = $(addprefix $(DSOLIBS), $(my_deps))

link2_libs = egspp egs_glib

include $(SPEC_DIR)egspp_libs.spec

$(make_depend)

$(DSO2)egs_linternal.$(obje): ..$(DSEP)..$(DSEP)geometry$(DSEP)egs_glib$(DSEP)egs_glib.h