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


/*! \file egs_internal_source.cpp
 *  \brief An internal source
 *  \IK
 */

#include "egs_internal_source.h"
#include "egs_input.h"
#include "egs_math.h"

EGS_InternalSource::EGS_InternalSource(EGS_Input *input,
        EGS_ObjectFactory *f) : EGS_BaseSimpleSource(input,f), geom(0),
    reg(0), wgt(0), nx(0), ny(0), nz(0), nxy(0) {
    string geom_name;
    int err1 = input->getInput("geometry",geom_name);
    if (!err1) {
        geom = (EGS_XYZGeometry*)EGS_BaseGeometry::getGeometry(geom_name);
        if (!geom)
			egsFatal("EGS_InternalSource: no geometry named %s\n", geom_name.c_str());
        if (geom->getType()!="EGS_XYZGeometry")
			egsFatal("EGS_InternalSource: geometry %s not type XYZ\n", geom_name.c_str());
    }
	else
		egsFatal("EGS_InternalSource: no \"geometry\" input found\n");
	xb = geom->getXPositions(); nx = geom->getNx();
    yb = geom->getYPositions(); ny = geom->getNy();
    zb = geom->getZPositions(); nz = geom->getNz();
	nxy = nx*ny;
	
	reg = new vector <int>;
	err1 = input->getInput("regions",*reg);
	
	wgt = new vector <EGS_Float>;
	int err2 = input->getInput("weights",*wgt);
	
    string activity_name;
	int err3 = input->getInput("activity file",activity_name);
	
    if (err1 && err3)
		egsFatal("EGS_InternalSource: no \"regions\" input found\n");
	
    if (!err1 && !err2) {
		if (reg->size() != wgt->size()) {
			egsFatal("EGS_InternalSource: mismatch in number of regions and weights\n");
		} else {
			table = new EGS_SimpleAliasTable(wgt->size(),wgt->data());
		}
    } else if (err2 && err3) {
		egsFatal("EGS_InternalSource: no \"weights\" input found\n"); // This was often too subtle a warning
		//table = new EGS_SimpleAliasTable(reg->size(),&vector<EGS_Float>(reg->size(),1.0)[0]);
	} else if (!err3) {
		ifstream actFile(activity_name.c_str());

        if (actFile) {
            int readInt;
			double readDub;
			
			while (actFile >> readInt) { // Get region
				reg->push_back(readInt);
				actFile >> readDub; // Get weight
				wgt->push_back(readDub);
			}
        }
		else {
			egsFatal("EGS_InternalSource: failed to open specified activity file\n");			
		}
		
		if (reg->size() != wgt->size()) {
			egsFatal("EGS_InternalSource: mismatch in number of regions and weights\n");
		} else {
			table = new EGS_SimpleAliasTable(wgt->size(),wgt->data());
		}
	} else {
		egsFatal("EGS_InternalSource: missing weight and region definition\n");		
	}
	
	if (!err1 && !err2 && !err3) {
		egsWarning("EGS_InternalSource: Found weight, region, and activity file input\n                    ignoring activity file\n");
	}
	
	// Output ---------------------------------------------------------------------------------
    otype = "EGS_InternalSource";
	description = "Internal source for ";
	description += geom->ref();
	description += " XYZ geometry";
	if (q == -1) {
		description += ", electrons";
	}
	else if (q == 0) {
		description += ", photons";
	}
	else if (q == 1) {
		description += ", positrons";
	}
	else {
		description += ", unknown particle type";
	}
}

extern "C" {

    EGS_INTERNAL_SOURCE_EXPORT EGS_BaseSource *createSource(EGS_Input *input,
            EGS_ObjectFactory *f) {
        return
            createSourceTemplate<EGS_InternalSource>(input,f,"internal source");
    }

}
