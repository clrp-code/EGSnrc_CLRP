/*
###############################################################################
#
#  EGSnrc egs++ envelope geometry headers
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
#  Author:          Randle Taylor (randle.taylor@gmail.com) on behalf of
#                   the Carleton Laboratory for Radiotherapy Physics
#
#  Contributors:
#
###############################################################################
*/

/*! \file egs_glib.cpp
 * \brief The createGeometry function for the library shim
 *  \author Randle Taylor (randle.taylor@gmail.com)
 **/

#include <map>
#include <iostream>
#include <istream>
#include <sstream>
#include <fstream>
#include "egs_input.h"
#include "egs_functions.h"
#include "egs_transformations.h"
#include "egs_glib.h"
#include "../egs_nd_geometry/egs_nd_geometry.h"
#include "../egs_planes/egs_planes.h"
#include "gzstream.h"

/*! \brief Split a string on input delimeter */
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


/*! \brief Split a string on input delimeter */
std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

/*! \brief trim whitespace from start of string */
static inline std::string &ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}

/*! \brief trim whitespace from end of string */
static inline std::string &rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
}

/*! \brief trim whitespace from both ends of string */
static inline std::string &trim(std::string &s) {
    return ltrim(rtrim(s));
}


/*! \brief parse density file */
map<string, EGS_Float> getMedRhos(ifstream &in) {
    const string med_delim = "MEDIUM=";
    const string rho_delim = "RHO= ";
    string cur_med_name;
    EGS_Float cur_density;
    map<string, EGS_Float> med_rhos;

    while (in) {
        string line;
        getline(in, line);

        bool new_med = line.find(med_delim) != string::npos;

        if (new_med) {
            cur_med_name = trim(split(trim(split(line, '=')[1]), ',')[0]);
            getline(in, line);
            cur_density = atof(trim(split(trim(split(line, '=')[1]), ',')[0]).c_str());
            med_rhos[cur_med_name] = cur_density;
        }
    }

    return med_rhos;
}

EGS_BaseGeometry *readEGSPhant(istream &data, map<string, EGS_Float> med_rhos) {

    int nmed;
    data >> nmed;

    // read in media names and create map from phant medium code character to egs++ name
    map<char, string> phant2egs;
    map<char, int> phant2egs_idx;
    string phant_meds_str = "123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    for (int i=0; i < nmed; i++) {
        string med;
        data >> med;
        char med_key = phant_meds_str.at(i);
        phant2egs_idx[med_key] = EGS_BaseGeometry::addMedium(med);
        phant2egs[med_key] = med;
    }

    // estepe is ignored!
    for (int i=0; i < nmed; i++) {
        EGS_Float estepe;
        data >> estepe;
    }

    /* read in all bounds and create the geometry */
    int nx, ny, nz, nreg;
    vector<EGS_Float> xbounds, ybounds,zbounds;
    EGS_Float bound;
    data >> nx >> ny >> nz;
    nreg = nx*ny*nz;

    for (int i=0; i < nx+1; i++) {
        data >> bound;
        xbounds.push_back(bound);
    }
    for (int i=0; i < ny+1; i++) {
        data >> bound;
        ybounds.push_back(bound);
    }
    for (int i=0; i < nz+1; i++) {
        data >> bound;
        zbounds.push_back(bound);
    }

    /* now we've got all geometry information so construct our geom */
    EGS_PlanesX *xp = new EGS_PlanesX(xbounds,"",EGS_XProjector("x-planes"));
    EGS_PlanesY *yp = new EGS_PlanesY(ybounds,"",EGS_YProjector("y-planes"));
    EGS_PlanesZ *zp = new EGS_PlanesZ(zbounds,"",EGS_ZProjector("z-planes"));
    EGS_XYZGeometry *result = new EGS_XYZGeometry(xp,yp,zp);

    // read in region media and set them in the geometry
    char cur_med;
    for (int i=0; i < nreg ; i++) {
        data >> cur_med;
        result->setMedium(i, i, phant2egs[cur_med]);
    }

    // read in region rhos and set the relative rho value if required
    EGS_Float cur_rho;
    for (int i=0; i < nreg ; i++) {

        data >> cur_rho;

        string reg_med = EGS_BaseGeometry::getMediumName(result->medium(i));
        if (med_rhos.find(reg_med) != med_rhos.end()) {
            if (cur_rho != med_rhos[reg_med]) {
                result->setRelativeRho(i, i, cur_rho/med_rhos[reg_med]);
            }
        }
        else {
            egsFatal("While construction geometry from egsphant file, the medium %s was not found in the density file\n", reg_med.c_str());
        }
    }

    return result;

}

bool isGZip(istream &vfile) {
    return (vfile.get() == 0x1f && vfile.get() == 0x8b);
}

/*! \brief read egsphant file and construct an XYZ Geometry from it */
EGS_BaseGeometry *parse_egsphant(string fpath, map<string, EGS_Float> med_rhos) {

    ifstream data(fpath.c_str(), ios::binary);

    if (!data) {
        return 0;
    }
    EGS_BaseGeometry *result;
    if (isGZip(data)) {
        data.close();
        igzstream gzf(fpath.c_str());
        result = readEGSPhant(gzf, med_rhos);
        gzf.close();
    }
    else {
        data.close();
        ifstream tdata(fpath.c_str());
        result = readEGSPhant(tdata, med_rhos);
    }

    return result;
}


extern "C" {

    /*! createGeometry function for glib shim */
    EGS_GLIB_EXPORT EGS_BaseGeometry *createGeometry(EGS_Input *input) {

        if (!input) {
            egsWarning("createGeometry(egs_glib): null input?\n");
            return 0;
        }


        vector<string> gtype_choices;
        gtype_choices.push_back("external");
        gtype_choices.push_back("egsphant");

        int gtype = input->getInput("type", gtype_choices, 0);
        EGS_BaseGeometry *final = 0;
        if (gtype == 0) {

            EGS_Input *egs_geom_input = input->takeInputItem("geometry definition");
            if (!egs_geom_input) {

                egsWarning("createGeometry(egs_glib): missing `geometry definition` input\n");
                return 0;
            }

            final = EGS_BaseGeometry::createGeometry(egs_geom_input);

            delete egs_geom_input;
        }
        else {
            string egsphant_file;
            input->getInput("egsphant file", egsphant_file);

            string density_file;
            bool missing_rho = input->getInput("density file", density_file) != 0;
            map<string, EGS_Float> default_med_rho;
            if (missing_rho) {
                egsWarning("\n\nEGS_Glib::EgsPhant input Missing 'density file' key. Default media densities will be used. \n\n");
            }
            else {

                ifstream rho_data(density_file.c_str());

                if (!rho_data) {
                    egsFatal("\n\nEGS_Glib::EgsPhant unable to open density file '%s'.", density_file.c_str());
                    return 0;
                }

                default_med_rho = getMedRhos(rho_data);
            }

            final = parse_egsphant(egsphant_file, default_med_rho);
        }

        if (!final) {
            egsWarning("createGeometry(egs_glib): unable to create geometry\n");
            return 0;
        }

        final->setName(input);

        return final;

    }

}
