/*
################################################################################
#
#  egs_mird egs_mird.cpp
#  Copyright (C) 2021 Ting Lee, Rowan Thomson, Martin Martinov
#
#  This file is part of egs_mird
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
#include <cstdlib>
#include <fstream>
#include <sstream>
// We derive from EGS_AdvancedApplication => need the header file.
#include "egs_advanced_application.h" 
// We use scoring objects provided by egspp => need the header file.
#include "egs_scoring.h"
// Every C++ EGSnrc application needs this header file
#include "egs_interface2.h"
// We use egsInformation() => need the egs_functions.h header file.
#include "egs_functions.h"
// We use the EGS_Input class
#include "egs_input.h"
// To get the maximum source energy
#include "egs_base_source.h"
// The random number generator
#include "egs_rndm.h"
// Transformations
#include "egs_transformations.h"
// Interpolators
#include "egs_interpolator.h"
// Parallel submission control
#include "egs_run_control.h"

// Coefficient to convert from Mev to Joules
#define MEV_TO_J 1.6021773e-13
#define MEV_TO_GY 1.6021773e-10

// Get XYZ geoms for building 3ddose arrays
#define EXPLICIT_XYZ
#include "geometry/egs_nd_geometry/egs_nd_geometry.h"

class APP_EXPORT egs_mird : public EGS_AdvancedApplication {
    EGS_ScoringArray *score;          // scoring energies deposited for 3ddose array 
    EGS_BaseGeometry *doseg;          // scoring geometry
    int               nreg;           // number of regions in the scoring geometry
	int               offset;         // swap from global region to dose region
	int				  discard;        // flag a particle to be discarded
	EGS_Float		  scaling;        // dose scaling factor
	
	bool                       score_tlen; // score tracklength
	vector <EGS_Interpolator*> muenDat;    // interpolators holding all the muen data
	bool                       regWeight_RR;  // internal source russian roulette
    vector<EGS_Float>    	  *regWeight_wgt; // holds all weights
    vector<int>          	  *regWeight_reg; // holds all regions
	EGS_Float				   regWeight_min; // lowest weight for RR
	EGS_Float				   regWeight_max; // highest weight for RR
	EGS_Float				   regWeight_mul; // multiplier for RR
	
	bool              binDose;  // option for outputting binary 3ddose files
	string            fileName; 
	string            doseGeom; 
	string            type;
    
public:
    // Constructor
    egs_mird(int argc, char **argv) :
        EGS_AdvancedApplication(argc,argv), score(0), nreg(0), binDose(false), scaling(1.0),
		fileName(""), type("") {};

    // Destructor
    ~egs_mird() {
        if (score) {
            delete score;
        }
        if (doseg) {
            delete doseg;
        }
		for (int i = 0; i < muenDat.size(); i++)
			delete muenDat[i];
		muenDat.clear();
		if (regWeight_wgt) {
			delete regWeight_wgt;
		}
    };

	// Initiliazatin functions
    void describeUserCode() const;
	void initTrackLengthScoring(EGS_Input*);
	void initRegionWeightRoulette(EGS_Input*);
    int initScoring();
	
	// Scoring functions
    int ausgab(int); 
	
	// Parallel submission functions
    int outputData();
    int readData();
    void resetCounter();
    int addState(istream &data);
    void getCurrentResult(EGS_Float &sum, EGS_Float &sum2, EGS_Float &norm, EGS_Float &count);

    // Final output functions
    void outputResults();
	int output_3ddose();
	int output_b3ddose();

protected:
    int simulateSingleShower();
    int startNewShower();
};

void egs_mird::describeUserCode() const {
    egsInformation(
        "\n               ***************************************************"
        "\n               *                                                 *"
        "\n               *                     egs_mird                    *"
        "\n               *                                                 *"
        "\n               ***************************************************"
        "\n\n");
}

int egs_mird::initScoring() {
    EGS_Input *options = input->takeInputItem("scoring options");
    if (options) {
		// Get scoring type and get output file name
		if (!options->getInput("type",type)) {
			if (type != "3ddose" && type != "b3ddose")
				egsFatal("initScoring(): %s does not name a valid type\n", type);
		}
		else {
			egsWarning("initScoring(): did not find 'type' input, assuming 3ddose output preferred\n");
			type = "3ddose";
		}
		if (options->getInput("scaling factor",scaling)) {
			egsWarning("initScoring(): did not find 'scaling factor' input, assuming no scaling\n");
		}
		if (options->getInput("file name",fileName)) {
			egsWarning("initScoring(): did not find 'file name' input, setting output to input filename\n");
		}
		string gName;
		if (options->getInput("scoring geometry",gName)) {
			egsWarning("initScoring(): did not find 'scoring geometry' input\n");
		}
		doseg = EGS_BaseGeometry::getGeometry(gName);
		if (!doseg) {
			egsWarning("createGeometry(gtransformed): no geometry named %s"
					   " is defined\n",gName.c_str());
			return 0;
		}
		
		// Check primary geometry for being type XYZ
		if (type == "3ddose" || type == "b3ddose") {
			if (doseg->getType() != "EGS_XYZGeometry")
				egsFatal("initScoring(): geometry %s is type %s, not EGS_XYZGeometry\n", doseg->getName().c_str(), doseg->getType().c_str());
			int checkEndian = 1;
			if (type == "b3ddose" && *((char*)(&checkEndian)) != 1)
				egsFatal("initScoring(): b3ddose file output is only supported on little endian machines,\n"
						 "please use 3ddose type for this installation\n");
		}
		
		nreg = doseg->regions(); // set nreg to geometry's regions
		score = new EGS_ScoringArray(nreg); // set scoring array to score every region
	}
	else {
		egsFatal("initScoring(): did not find scoring options block\n");
	}
	
	delete options;
	
	// A fairly not robust method, querying what geometry returns and doseg region numbers for voxel 0 in doseg;
	EGS_Vector offPoint(doseg->getBound(0,0)+(doseg->getBound(0,1)-doseg->getBound(0,0))/100.0,
	doseg->getBound(1,0)+(doseg->getBound(1,1)-doseg->getBound(1,0))/100.0,
	doseg->getBound(2,0)+(doseg->getBound(2,1)-doseg->getBound(2,0))/100.0);
	offset = geometry->isWhere(offPoint);
	
	// Variance reduction technique setup
    options = input->takeInputItem("variance reduction");
	egsInformation("Variance Reduction techniques:\n"
				   "==============================\n");
	if (options) {

		// Check for track length scoring
		initTrackLengthScoring(options);
		
		// Check for internal source roulette
		initRegionWeightRoulette(options);
	}
	else {
		egsInformation("no variance reduction\n\n");
		score_tlen = 0;
		regWeight_RR = 0;
	}
		
    egsInformation("==============================\n");
	delete options;
		
    return 0;
}

int egs_mird::ausgab(int iarg) {
	int ir = top_p.ir-offset;
	if (ir < 0 || nreg <= ir) {
		return 0;
	}

    if (score_tlen && !top_p.q && !iarg) {
		EGS_Float muen_val = muenDat[the_useful->medium-1]->interpolateFast(the_epcont->gle);
		
		EGS_Float aux = the_epcont->tvstep*top_p.E*muen_val*top_p.wt;
		EGS_Float vol = geometry->getVolume(ir);

		if (vol > 0 && aux > 0) {
			aux /= vol;
			score->score(ir, aux);
		}
    }
	else if (!score_tlen && iarg <= ExtraEnergy) {
		EGS_Float aux = the_epcont->edep*top_p.wt;
		score->score(ir, aux);
	}

    return 0;
}

int egs_mird::outputData() {
    // Invoke base class first
    int err = EGS_AdvancedApplication::outputData();
    if (err) {
        return err;
    }
	
	// Store score
    if (!score->storeState(*data_out)) {
        return 101;
    }
		
    return 0;
}

int egs_mird::readData() {
    // Invoke base class first
    int err = EGS_AdvancedApplication::readData();
    if (err) {
        return err;
    }
  
    // Use EGS_AdvancedApplication's data_in to load score
    if (!score->setState(*data_in)) {
        return 201;
    }
    return 0;
}

void egs_mird::resetCounter() {
    // Invoke base class first
    EGS_AdvancedApplication::resetCounter();
	
	// Reset dose scoring arrays
    score->reset();
}

int egs_mird::addState(istream &data) {
    int err = EGS_AdvancedApplication::addState(data);
    if( err ) return err;
	
    // Read stored score in data into tmp then add it to score
    EGS_ScoringArray tmp(nreg);
    if (!tmp.setState(data)) {
        return 301;
    }
    (*score) += tmp;
	
	return 0;
}

void egs_mird::outputResults() {
    egsInformation("\n\n last case = %lld fluence = %g\n\n", current_case, source->getFluence());
	EGS_Float norm = MEV_TO_GY*current_case/source->getFluence()*scaling;
	int nx = doseg->getNRegDir(0);
	int ny = doseg->getNRegDir(1);
	int nz = doseg->getNRegDir(2);
	int nxy = nx*ny;
	EGS_Float norm2  = the_media->rho[doseg->medium(nreg/2)]*doseg->getRelativeRho(nreg/2);
			  norm2 *= (doseg->getBound(0,int(nx/2)+1)-doseg->getBound(0,int(nx/2)));
			  norm2 *= (doseg->getBound(1,int(ny/2)+1)-doseg->getBound(1,int(ny/2)));
			  norm2 *= (doseg->getBound(2,int(nz/2)+1)-doseg->getBound(2,int(nz/2)));
	double sum, sumE;
	
	if (!fileName.compare("")) {
		fileName = final_output_file;
		fileName += ".";
		fileName += type;
	}
	
	if (!final_job && output_file != final_output_file) {
		// do nothing if we aren't fully finished the simulation
	}
	else if (type == "3ddose") {
		output_3ddose();
		egsInformation("3ddose file output at %s\n", fileName.c_str());
		score->currentResult(nreg/2, sum, sumE);
		egsInformation("Dose at (center) region %d - %12.6e +/- %.5f\%\n", nreg/2, sum*norm/(score_tlen?1.0:norm2), 100.0*double(sum==0?1.0:sumE/sum));
	}
	else if (type == "b3ddose") {
		output_b3ddose();
		egsInformation("b3ddose file output at %s\n", fileName.c_str());
		score->currentResult(nreg/2, sum, sumE);
		egsInformation("Dose at (center) region %d - %12.6e +/- %.5f\%\n", nreg/2, sum*norm/(score_tlen?1.0:norm2), 100.0*double(sum==0?1.0:sumE/sum));
	}
	else {
		egsInformation("outputResults(): Nothing was output\n");
		score->currentResult(nreg/2, sum, sumE);
		egsInformation("Dose at (center) region %d - %12.6e +/- %.5f\%\n", nreg/2, sum*norm/(score_tlen?1.0:norm2), 100.0*double(sum==0?1.0:sumE/sum));
	}
}

void egs_mird::getCurrentResult(EGS_Float &sum, EGS_Float &sum2, EGS_Float &norm, EGS_Float &count) {
	norm = MEV_TO_GY*current_case/source->getFluence()*scaling;
	count = current_case;
	if (!score_tlen) {
		int nx = doseg->getNRegDir(0);
		int ny = doseg->getNRegDir(1);
		int nz = doseg->getNRegDir(2);
		int nxy = nx*ny;
		EGS_Float norm2  = the_media->rho[doseg->medium(nreg/2)]*doseg->getRelativeRho(nreg/2);
				  norm2 *= (doseg->getBound(0,int(nx/2)+1)-doseg->getBound(0,int(nx/2)));
				  norm2 *= (doseg->getBound(1,int(ny/2)+1)-doseg->getBound(1,int(ny/2)));
				  norm2 *= (doseg->getBound(2,int(nz/2)+1)-doseg->getBound(2,int(nz/2)));
		norm /= norm2;
	}
	score->currentScore(nreg/2, sum, sum2);
}

int egs_mird::simulateSingleShower() {
    int ireg;
	int ntry = 0;
	discard = 0;
	last_case = current_case;
	setEdep(0);
	do {
		ntry++;
		if (ntry > 100000) {
			egsWarning("EGS_Application::simulateSingleShower(): no particle"
					   " from the source has entered the geometry after 100000"
					   " attempts\n");
			return 1;
		}
		current_case = source->getNextParticle(rndm,p.q,p.latch,p.E,p.wt,p.x,p.u);
		ireg = geometry->isWhere(p.x);
		if (ireg < 0) {
			EGS_Float t = veryFar;
			ireg = geometry->howfar(ireg,p.x,p.u,t);
			if (ireg >= 0) {
				p.x += p.u*t;
			}
		}
	} while (ireg < 0);
	
	p.ir = ireg; // global geometry region is set

	// if getEdep() is non-zero, then egs_radionuclide_source has invoked userScoring from egs_application
	// so we score the same here for our 3ddose array
	int ir = top_p.ir-offset;
	if (ir > -1 && nreg > ir && ireg && getEdep())
		score->score(ir, getEdep()*(!p.wt?1:p.wt)); // egs_radionuclide_source sets weight to 0, so use 1
													// for scoring here
	int err = startNewShower();
	if (err) {
		return err;
	}
	
	ireg = doseg->isWhere(p.x); // ireg now holds scoring geometry local region
	
	if (score_tlen && p.q) { // Tracklength scoring and a non-photon
		EGS_Float aux = p.E;
		if (ireg >= 0 && ireg < nreg && aux > 0) { // Score it
			aux *= p.wt; // Scale by particle weight
			EGS_Float norm2 = the_media->rho[doseg->medium(ireg)]*doseg->getRelativeRho(ireg);
			norm2 *= geometry->getVolume(ireg);
			aux /= norm2; // Normalize to dose to match track length
			score->score(ireg,aux);
		}
		discard++; // flag a discard
	}
	
	if (regWeight_RR && !discard) { // If RR and particle not discarded
		if (ireg >= 0 && ireg < nreg && regWeight_wgt->at(ireg) > 0) {
			// Sample against the region weight
			if (rndm->getUniform() < regWeight_wgt->at(ireg)) { // If it survives, increase weight appropriately
				p.wt /= regWeight_wgt->at(ireg);
			}
			else { // Else, kill it
				discard++; // flag a discard
			}
		}
	}
	
	if (discard)
		p.wt = 0;
	
	err = shower();
	if (err) {
		return err;
	}
	
	return finishShower();
}

int egs_mird::startNewShower() {
	int res = EGS_Application::startNewShower();
	if( res ) return res;
	
	if (current_case != last_case) {
		score->setHistory(current_case);
		last_case = current_case;		
	}
	
	return 0;
};

int egs_mird::output_3ddose() {
	EGS_Float norm = MEV_TO_GY*current_case/source->getFluence()*scaling;
	
	int nx = doseg->getNRegDir(0);
	int ny = doseg->getNRegDir(1);
	int nz = doseg->getNRegDir(2);
	int nxy = nx*ny;
	int ir;
	EGS_Float bound, r, dr, norm2, rE, drE;
	
	ofstream out;
	out.open(fileName, ios::out);
	if (out.is_open()) {
		// Output voxel count
        out << nx << " " << ny << " " << nz << "\n";
		
        // Output voxel boundaries
		out.precision(6); // set boundary precision to 6
        for (int i = 0; i <= nx; i++) {
            bound = doseg->getBound(0,i);
            out << double(abs(bound)<epsilon?0.0:bound) << " ";
        }
        out << "\n";
        for (int j = 0; j <= ny; j++) {
            bound = doseg->getBound(1,j);
            out << double(abs(bound)<epsilon?0.0:bound) << " ";
        }
        out << "\n";
        for (int k = 0; k <= nz; k++) {
            bound = doseg->getBound(2,k);
            out << double(abs(bound)<epsilon?0.0:bound) << " ";
        }
        out << "\n";
		
		// Output doses
		out.precision(8); // set dose precision to 8
        for (int k = 0; k < nz; k++)
			for (int j = 0; j < ny; j++)
				for (int i = 0; i < nx; i++) {
					ir = i+j*nx+k*nxy;
					score->currentResult(ir,r,dr);
					if (!score_tlen) {
						norm2 = the_media->rho[doseg->medium(ir)]*doseg->getRelativeRho(ir);
						norm2 *= (doseg->getBound(0,i+1)-doseg->getBound(0,i));
						norm2 *= (doseg->getBound(1,j+1)-doseg->getBound(1,j));
						norm2 *= (doseg->getBound(2,k+1)-doseg->getBound(2,k));
					}
					else {
						norm2 = 1;
					}
					out << double(r*norm/norm2) << " ";
				}
        out << "\n";
		
		// Output uncertainty
		out.precision(6); // set relative error precision to 6
        for (int k = 0; k < nz; k++)
			for (int j = 0; j < ny; j++)
				for (int i = 0; i < nx; i++) {
					ir = i+j*nx+k*nxy;
					score->currentResult(ir,r,dr);
					out << double(r==0?1.0:dr/r) << " ";
				}
        out << "\n";
        out.close();
		return 1;
    }
    else {
        egsFatal("\noutput_3ddose(): Failed to open file %s for writing.\n", fileName);
		return 0;
    }
}

int egs_mird::output_b3ddose() {
	EGS_Float norm = MEV_TO_GY*current_case/source->getFluence()*scaling;
	int nx = doseg->getNRegDir(0);
	int ny = doseg->getNRegDir(1);
	int nz = doseg->getNRegDir(2);
	int nxy = nx*ny;
	int ir;
	char* buffer, temp;
	double bound;
	int nb = sizeof(bound);
	EGS_Float r, dr, norm2, rE, drE;
	
	ofstream out;
	out.open(fileName, ios::out | ios::binary);
	if (out.is_open()) {
		// Output voxel count
		buffer = new char (1);
		out.write(buffer,1); // Output XYZ filetype
		delete buffer;
		
		// Output voxel count in x, y, and z
		buffer = (char*)(&nx);
		out.write(buffer,sizeof(nx)); 
		buffer = (char*)(&ny);
		out.write(buffer,sizeof(nx)); 
		buffer = (char*)(&nz);
		out.write(buffer,sizeof(nx)); 
		
        // Output voxel boundaries
		buffer = (char*)(&bound);
        for (int i = 0; i <= nx; i++) {
            bound = geometry->getBound(0,i);
			out.write(buffer,nb); 
        }
        for (int j = 0; j <= ny; j++) {
            bound = geometry->getBound(1,j);
			out.write(buffer,nb); 
        }
        for (int k = 0; k <= nz; k++) {
            bound = geometry->getBound(2,k);
			out.write(buffer,nb); 
        }
		
		// Output doses
        for (int k = 0; k < nz; k++)
			for (int j = 0; j < ny; j++)
				for (int i = 0; i < nx; i++) {
					ir = i+j*nx+k*nxy;
					score->currentResult(ir,r,dr);
					if (!score_tlen) {
						norm2 = the_media->rho[doseg->medium(ir)]*doseg->getRelativeRho(ir);
						norm2 *= (doseg->getBound(0,i+1)-doseg->getBound(0,i));
						norm2 *= (doseg->getBound(1,j+1)-doseg->getBound(1,j));
						norm2 *= (doseg->getBound(2,k+1)-doseg->getBound(2,k));
					}
					else {
						norm2 = 1; 
					}
					bound = r*norm/norm2;
					out.write(buffer,nb); 
				}
		
		// Output uncertainty
        for (int k = 0; k < nz; k++)
			for (int j = 0; j < ny; j++)
				for (int i = 0; i < nx; i++) {
					ir = i+j*nx+k*nxy;
					score->currentResult(ir,r,dr);
					bound = r==0?1.0:dr/r;
					out.write(buffer,nb); 
				}
		
        out.close();
		return 1;
    }
    else {
        egsFatal("\noutput_b3ddose(): Failed to open file %s for writing.\n", fileName);
		return 0;
    }	
}

// This input process mirrors that of egs_internal_source and is meant to be used
// primarily with it
void egs_mird::initRegionWeightRoulette(EGS_Input *VR_options) {
	vector<string> choices;
    choices.push_back("no");
    choices.push_back("yes");

    regWeight_RR = VR_options->getInput("region weight roulette", choices, 0);
    if (!regWeight_RR) {
		egsInformation("Region weight roulette is off\n");
		return;
	}
	egsInformation("Region weight roulette is on\n");
	
    regWeight_reg = new vector <int>;
	int err1 = VR_options->getInput("regions",*regWeight_reg);
	
	regWeight_wgt = new vector <EGS_Float>;
	int err2 = VR_options->getInput("weights",*regWeight_wgt);
	
    string file_name;
	int err3 = VR_options->getInput("region weight file",file_name);
	
    if (err1 && err3)
		egsFatal("initRegionWeightRoulette: no \"regions\" input found\n");
	
    if (!err1 && !err2) {
		if (regWeight_reg->size() != regWeight_wgt->size()) {
			egsFatal("initRegionWeightRoulette: mismatch in number of regions and weights\n");
		}
		egsInformation("Region weights loaded from inputs\n");
    } else if (err2 && err3) {
		egsFatal("initRegionWeightRoulette: no \"weights\" input found\n");
	} else if (!err3) {
		ifstream actFile(file_name);

        if (actFile) {
            int readInt;
			double readDub;
			
			while (actFile >> readInt) { // Get region
				regWeight_reg->push_back(readInt);
				actFile >> readDub; // Get weight
				regWeight_wgt->push_back(readDub);
			}
        }
		else {
			egsFatal("initRegionWeightRoulette: failed to open specified weight file\n");			
		}
		
		if (regWeight_reg->size() != regWeight_wgt->size()) {
			egsFatal("initRegionWeightRoulette: mismatch in number of regions and weights\n");
		}
		egsInformation("Region weights loaded from file\n");
	} else {
		egsFatal("initRegionWeightRoulette: missing weight and region definition\n");		
	}
	
	if (!err1 && !err2 && !err3) {
		egsWarning("initRegionWeightRoulette: found weight, region, and weight file input\n"
				   "                          ignoring weight file\n");
	}
	
	regWeight_min = regWeight_wgt->at(0);
	regWeight_max = regWeight_wgt->at(0);
	for (int i = 1; i < regWeight_wgt->size(); i++) {
		if (regWeight_wgt->at(i) < regWeight_min)
			regWeight_min = regWeight_wgt->at(i);
		if (regWeight_wgt->at(i) > regWeight_max)
			regWeight_max = regWeight_wgt->at(i);
	}
	
	regWeight_mul = 2.0;
    if (VR_options->getInput("relative roulette weight", regWeight_mul) != 0) {
        egsWarning("initRegionWeightRoulette(): relative roulette weight input not found\n"
		           "                             assuming weight of 2.\n");
    }
	
	// swap from selected regions to all regions for faster lookup
    vector<EGS_Float> *regWeight_wgtFinal = new vector<EGS_Float> (nreg, 0.0);
	for (int i = 0; i < regWeight_reg->size(); i++) {
		(*regWeight_wgtFinal)[regWeight_reg->at(i)] = 1.0-(regWeight_wgt->at(i)-regWeight_min)/(regWeight_max-regWeight_min)*1.0/(regWeight_mul);
	}
	delete regWeight_wgt;
	delete regWeight_reg;
	regWeight_wgt = regWeight_wgtFinal;
	
	egsInformation("Region weights set, RR ranging as high as %4.2f%% for areas with max activity\n",100.0-1.0/regWeight_mul*100.0);
}

/***************************************************************************************/
// The code below has been taken in large part from egs_brachy.cpp and muen.h
void egs_mird::initTrackLengthScoring(EGS_Input *VR_options) {
    vector<string> choices;
    choices.push_back("no");
    choices.push_back("yes");

    score_tlen = VR_options->getInput("score tracklength dose", choices, 0);
    if (!score_tlen) {
		egsInformation("Tracklength scoring is off\n\n");
		return;
	}
	egsInformation("Tracklength scoring is on\n");
	
    string muen_file;
    if (VR_options->getInput("muen file", muen_file) != 0) {
        egsFatal("\ninitTrackLengthScoring(): muen file input not found\n");
    }
	
	ifstream muen_data(muen_file.c_str());
	vector <string> mediumNames;
	vector <EGS_Interpolator*> allMuenDat;
	if (muen_data.is_open()) {
		string line;
		int num;
		double minEn, maxEn;
		string tempStr;
		double tempDouble;
		EGS_Float* muenInput;
		while (muen_data.is_open()) {
			std::getline(muen_data,line);
			if (line.size() < 33) // Length of the string before medium name
				break;
			
			mediumNames.push_back(line.substr(32,-1));
			
			std::getline(muen_data,line);
			std::getline(muen_data,line);
			num = std::stoi(line.substr(30,-1), NULL);
			
			if (num < 2)
				egsFatal("\ninitTrackLengthScoring(): error reading data for medium %s\n", mediumNames.back().c_str());
			
			muenInput = new EGS_Float[num];
			getline(muen_data,line);
			getline(muen_data,line);
			istringstream pre_iss(line);
			pre_iss >> minEn >> tempStr >> tempDouble;
			muenInput[0] = tempDouble;
			
			//std::cout << "|" << mediumNames.back().c_str() << "|\n"; std::cout.flush();
			//std::cout << "\t" << num << " datapoints.\n"; std::cout.flush();
			num--;
			for (int i = 1; i < num; i++) {
				getline(muen_data,line);
				istringstream iss(line);
				iss >> tempDouble >> tempStr;
				iss	>> tempDouble;
				muenInput[i] = tempDouble;
			}
			
			getline(muen_data,line);
			istringstream post_iss(line);
			post_iss >> maxEn >> tempStr >> tempDouble;
			muenInput[num] = tempDouble;
			
			allMuenDat.push_back(new EGS_Interpolator(num+1, log(minEn), log(maxEn), muenInput));
		}
		egsInformation("muendat file read\n");
	}
	else {
        egsFatal("\ninitTrackLengthScoring(): failed to open muendat file\n");
	}
	
	string mediumName;
	for (int i = 0; i < getnMedia(); i++) {
		mediumName = getMediumName(i);
		for (int j = 0; j < mediumNames.size(); j++)
			if (!mediumNames[j].compare(mediumName)) {
				muenDat.push_back(allMuenDat[j]);
				allMuenDat[j] = NULL;
			}
		if (muenDat.size() <= i) {
			egsFatal("\ninitTrackLengthScoring(): failed to find %s in muendat file\n", mediumName.c_str());
		}
	}
	egsInformation("muen data found for all simulation media\n\n");
	
	for (int i = 0; i < allMuenDat.size(); i++)
		delete allMuenDat[i];
	allMuenDat.clear();
}

#ifdef BUILD_APP_LIB
APP_LIB(egs_mird);
#else
APP_MAIN(egs_mird);
#endif
