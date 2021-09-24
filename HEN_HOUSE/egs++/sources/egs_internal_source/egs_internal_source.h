/*
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
*/

#ifndef EGS_INTERNAL_SOURCE_
#define EGS_INTERNAL_SOURCE_

#include "egs_vector.h"
#include "egs_base_source.h"
#include "egs_rndm.h"
#include "egs_shapes.h"
#include "egs_base_geometry.h"
#include "geometry/egs_glib/egs_glib.h"
#include "../../egs_alias_table.h"
#include "egs_math.h"


#ifdef WIN32

    #ifdef BUILD_INTERNAL_SOURCE_DLL
        #define EGS_INTERNAL_SOURCE_EXPORT __declspec(dllexport)
    #else
        #define EGS_INTERNAL_SOURCE_EXPORT __declspec(dllimport)
    #endif
    #define EGS_INTERNAL_SOURCE_LOCAL

#else

    #ifdef HAVE_VISIBILITY
        #define EGS_INTERNAL_SOURCE_EXPORT __attribute__ ((visibility ("default")))
        #define EGS_INTERNAL_SOURCE_LOCAL  __attribute__ ((visibility ("hidden")))
    #else
        #define EGS_INTERNAL_SOURCE_EXPORT
        #define EGS_INTERNAL_SOURCE_LOCAL
    #endif

#endif

class EGS_INTERNAL_SOURCE_EXPORT EGS_InternalSource :
    public EGS_BaseSimpleSource {

public:
	EGS_InternalSource(EGS_Input *, EGS_ObjectFactory *f=0);
	
    ~EGS_InternalSource() {
        if (geom)
            delete geom;
        if (reg) {
			vector<int>().swap(*reg);
            delete reg;
		}
        if (wgt) {
			vector<EGS_Float>().swap(*wgt);
            delete wgt;
		}
        if (table)
            delete table;
    };

    void getPositionDirection(EGS_RandomGenerator *rndm,
                              EGS_Vector &x, EGS_Vector &u, EGS_Float &wt) {
		// Set position
		int ind = (*reg)[table->sample(rndm)], xi, yi, zi;
		
		xi = ind%nx;
		x.x = rndm->getUniform()*(xb[xi+1]-xb[xi])+xb[xi];
		
		yi = int(ind/nx)%ny;
		x.y = rndm->getUniform()*(yb[yi+1]-yb[yi])+yb[yi];
		
		zi = ind/nxy;
		x.z = rndm->getUniform()*(zb[zi+1]-zb[zi])+zb[zi];
		
		// Set isotropic direction
		EGS_Float theta = rndm->getUniform()*2.0*M_PI, phi = acos(1.0-2.0*rndm->getUniform());
		u.z = cos(phi);
		u.y = sin(phi)*cos(theta);
		u.x = sin(phi)*sin(theta);
		
		// Set weight (always 1)
        wt = 1;
    };

    EGS_Float getFluence() const {
        return count;
    };

    bool storeFluenceState(ostream &) const {
        return true;
    };

    bool setFluenceState(istream &) {
        return true;
    };

protected:
    EGS_XYZGeometry      *geom;
    EGS_Float            *xb; // This will be a vector pointers
    EGS_Float            *yb; // This will be a vector pointers
    EGS_Float            *zb; // This will be a vector pointers
	EGS_I64               nx, ny, nz, nxy; // Number of voxels
    vector<EGS_Float>    *wgt; // Holds all weights
    vector<int>          *reg; // Holds all regions
    EGS_SimpleAliasTable *table; // Alias table to go from weight to region
};

#endif
