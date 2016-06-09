
###############################################################################
#
#  EGSnrc egs++ sample autoenvelope geometry
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
#  This input file defines an egs++ geometry for a series of brachytherapy
#  seeds put into a xyz-geometry filled with water. It is very similar to the
#  seeds_in_xyz.geom file but now the seeds are inscribed in a xyz geometry
#  using an Auto Envelope geometry instead of a traditional envelope.
#  It also makes use of the glib library to load a brachytherapy seed from
#  an external file.
#
###############################################################################


:start geometry definition:

    ######################## begin of seed definition #######################

    :start geometry:

        library      = egs_glib
        name         = seed
        include file = I6702.inp

    :stop geometry:

    ######################## define XYZ-geometry for envelope
    :start geometry:
        library   = egs_ndgeometry
        type      = EGS_XYZGeometry
        name      = xyz
        x-planes  = -5.5 -4.5 -3.5 -2.5 -1.5 -0.5 0.5 1.5 2.5 3.5 4.5 5.5
        y-planes  = -5.5 -4.5 -3.5 -2.5 -1.5 -0.5 0.5 1.5 2.5 3.5 4.5 5.5
        z-planes  = -5.5 -4.5 -3.5 -2.5 -1.5 -0.5 0.5 1.5 2.5 3.5 4.5 5.5
        :start media input:
            media = H2O700ICRU
        :stop media input:
    :stop geometry:


    ######################## define autoenvelope with inscribed seeds
    :start geometry:

        name = seeds_in_box
        library = egs_autoenvelope
        type = EGS_AEnvelope
        base geometry = xyz

        :start inscribed geometry:

            inscribed geometry name = seed

            :start transformations:

                :start transformation:
                    translation = -3 -3 0
                    rotation vector = 1 -1 0
                :stop transformation:

                :start transformation:
                    translation = -2 -2 0
                :stop transformation:

                :start transformation:
                    translation = -1 -1 0
                    rotation vector = 1 -1 0
                :stop transformation:

                :start transformation:
                    translation =  1  1 0
                    rotation vector = 1 -1 0
                :stop transformation:

                :start transformation:
                    translation =  2  2 0
                :stop transformation:

                :start transformation:
                    translation =  3  3 0
                    rotation vector = 1 -1 0
                :stop transformation:

                :start transformation:
                    translation = 3 -3 0
                :stop transformation:

                :start transformation:
                    translation = 2 -2 0
                :stop transformation:

                :start transformation:
                    translation = 1 -1 0
                :stop transformation:

                :start transformation:
                    translation = -1  1 0
                :stop transformation:

                :start transformation:
                    translation = -2  2 0
                :stop transformation:

                :start transformation:
                    translation = -3  3 0
                :stop transformation:

            :stop transformations:

            :start volume correction:

                correction type = correct
                density of random points (cm^-3) = 1E8

                :start shape:

                    type = cylinder
                    radius = 0.04
                    height = 0.45

                :stop shape:

            :stop volume correction:

        :stop inscribed geometry:
    :stop geometry:

    simulation geometry = seeds_in_box

:stop geometry definition:

