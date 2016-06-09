
###############################################################################
#
#  EGSnrc egs++ sample glib geometry
#  Copyright (C) 2015 National Research Council Canada
#
#  This file is part of EGSnrc.
#
#  EGSnrc is free software: you can redistribute it and/or modify it under
#  the terms of the GNU Affero General Public License as published by the
#  Free Software Foundation, either version 3 of the License, or (at your
#  option) any later version.
#
#  EGSnrc is distributed in the hope that it will be useful, but WITHOUT ANY
#  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
#  FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for
#  more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with EGSnrc. If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################
#
#  Author:          Randle Taylor, 2016
#
#  Contributors:
#
###############################################################################
#
#  An example geometry input file for the egs++ geometry package.
#
#  This input demonstrates loading an external geometry using the glib
#  library and inscribing it in an envelope geometry
#
###############################################################################


:start geometry definition:

    :start geometry:
        library = egs_glib
        type = egsphant
        name = my_egsphant_geom
        egsphant file = glib_egsphant.egsphant
        density file = glib_egsphant.density
    :stop geometry:

    :start geometry:

        library      = egs_glib
        name         = seed
        include file = I6702.inp

    :stop geometry:

    :start geometry:

        library              = egs_genvelope
        base geometry        = my_egsphant_geom
        inscribed geometries = seed
        name = phantom_w_seed

    :stop geometry:

    simulation geometry = phantom_w_seed

:stop geometry definition:

