/*
###############################################################################
#
#  EGSnrc egs++ glib geometry headers
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
# egs_glib was developed for Carleton Laboratory for Radiotherapy
# Physics (Rowan Thomson, Dave Rogers).
#
###############################################################################
*/

/*! \file egs_glib.h
 *
 *  \brief A small shim for loading geometries from files.
 *  \author Randle Taylor (randle.taylor@gmail.com)
 */

// need EXPLICIT_XYZ to get nd_geometry to include XYZ geometry
#define EXPLICIT_XYZ
#include <string>
#include <algorithm>
#include <fstream>
#include "egs_base_geometry.h"
#include "../egs_nd_geometry/egs_nd_geometry.h"

/*!

\ingroup Geometry


egs_glib is a shim to be used in conjunction with the EGS Input `include file`
directive for creating geometries defined in external files.

This is useful when you have a geometry defined in an external file:
\verbatim

 # your egsinp/geom file
:start geometry definition:

    :start geometry:
        library = egs_glib
        name = my_external_geom
        include file = /path/to/some/external/geom
    :stop geometry:

    :start geometry:
        library = egs_ndgeometry
        type = EGS_XYZGeometry
        name = my_base_geom
        x-planes =  -10 -5 0 5 10
        y-planes =  -10 -5 0 5 10
        z-planes =  -10 -5 0 5 10

        :start media input:
            media = WATER
        :stop media input:

    :stop geometry:

    :start geometry:
        library = egs_genvelope
        name = my_envelope_geom
        base geometry = my_base_geom
        inscribed geometries = my_external_geom
    :stop geometry:

    simulation geometry = my_envelope_geom

:stop geometry definition:

\endverbatim

and then in your external file you would have something like:

\verbatim

 #The external geom (/path/to/some/external/geom)
:start geometry definition:

    :start geometry:
        library   = egs_spheres
        midpoint  = 0 0 0
        radii     = 4
        name      = the_sphere
        :start media input:
            media = medium1
        :stop media input:
    :stop geometry:

    simulation geometry = the_sphere

:stop geometry definition:

\endverbatim

The above two definitions would result in a geometry equivalent to:

\verbatim

:start geometry definition:

    :start geometry:
        library   = egs_spheres
        midpoint  = 0 0 0
        radii     = 4
        name      = the_sphere
        :start media input:
            media = medium1
        :stop media input:
    :stop geometry:

    :start geometry:
        library = egs_ndgeometry
        type = XYZGeometry
        name = my_base_geom
        x-planes =  -10 5 0 5 10
        y-planes =  -10 5 0 5 10
        z-planes =  -10 5 0 5 10
    :stop geometry:

    :start geometry:
        library = egs_genvelope
        name = my_envelope_geom
        base geometry = my_base_geom
        inscribed geometries =  the_sphere
    :stop geometry:

    simulation geometry = my_envelope_geom

:stop geometry definition:

\endverbatim

Creating geometries based on egsphant files
-------------------------------------------
 egs_glib can also be used to construct an EGS_XYZGeometry based on
 an existing egsphant file.

\verbatim

    :start geometry:
        library = egs_glib
        type = egsphant
        name = my_egsphant_geom
        egsphant file = /path/to/some/egsphant/file
        density file = /path/to/density/file
    :stop geometry:

\endverbatim

by specifying `type=egsphant` in your input, egs_glib will parse the
egsphant file (either plane .egsphant text file or gzipped .egsphant.gz file)
you specify and construct a matching geometry.  The density file is
required so that the geometry may set the relative density of each
region correctly. Typically you would just set the `density file`
key to point to the pegs4dat file that you will
be using for your simulation although any file with a format:

\verbatim

MEDIUM=WATER
RHO=1.000
MEDIUM=AIR
RHO=1.2E-3

\endverbatim

will also work fine.



Examples
--------

Examples of using the glib library are available in the glib.geom and
seeds_in_xyz_aenv.geom files.

*/

#ifndef EGS_GLIB_
#define EGS_GLIB_

#ifdef WIN32

    #define EGS_GLIB_EXPORT __declspec(dllexport)
    #define EGS_GLIB_LOCAL

#else

    #ifdef HAVE_VISIBILITY
        #define EGS_GLIB_EXPORT __attribute__ ((visibility ("default")))
        #define EGS_GLIB_LOCAL  __attribute__ ((visibility ("hidden")))
    #else
        #define EGS_GLIB_EXPORT
        #define EGS_GLIB_LOCAL
    #endif

#endif

#endif


