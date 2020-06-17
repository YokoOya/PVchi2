#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "fitsio.h"

using namespace std;

#define EPS 1e-308
#define INF 1e+308

#define naxis 2
#define nPairsLim 6
#define delim " "
#define noval 0.0
#define calcRatioLim 100
#define lenfilename 512
#define lenline 512
#define DataDir "fits/"
#define OutDir "out/"

const bool f_skipNegative = true;
const bool logger = false;
const double db_erratio = 5e-15;


FILE *fp_out;

void printerror(int status) {
	if (status) {
		fits_report_error(stderr, status);
		exit(status);
	}
	return;
}

//////////
class PVdata {
public:
	long npix[naxis];
	char filename[lenfilename];
	fitsfile *fptr;
	double **valData;
	double valMax, valMin;
	double crpix[2], crval[2], cdelt[2];
	
	
	PVdata();
	//PVdata(const PVdata &);
	~PVdata();
	
	void dataInit();
	void pixInit(long *npix_in, double *crpix_in, double *crval_in, double *cdelt_in);
	
	void readFits(char* filename_in);
	
	void setMaxMin();
	void getMaxMin(double *max, double *min);
	
	bool skipCell(int ip, int iv, bool f_threshold, double threshold4skip);
	void setVal(int ip, int iv, double val);
	double getVal(int ip, int iv);
	
	double getP(int ip);
	double getIP(double p);
	double getV(int iv);
	double getIV(double v);
};


PVdata::PVdata() {
}

PVdata::~PVdata() {
	free(valData[0]);
	free(valData);
	valData[0] = NULL;
	valData = NULL;
}

void PVdata::dataInit() {
	valData = (double**) malloc(sizeof(double*) * npix[0]);
	valData[0] = (double*) malloc(sizeof(double) * npix[0] * npix[1]);
	
	for (int ip = 0; ip < npix[0]; ++ip) {
		valData[ip] = valData[0] + npix[1] * ip;
		for (int iv = 0; iv < npix[1]; ++ iv) {
			valData[ip][iv] = 0.;
		}
	}
	return ;
}

void PVdata::pixInit(long *npix_in, double *crpix_in, double *crval_in, double *cdelt_in) {
	for (int i = 0; i < 2; ++i) {
		npix[i] = npix_in[i];
		crpix[i] = crpix_in[i];
		crval[i] = crval_in[i];
		cdelt[i] = cdelt_in[i];
	}
	
	dataInit();
	return ;
}

void PVdata::readFits(char* filename_in) {
	sscanf(filename_in, "%s", filename);
	
	int status = 0, nfound, dumint, firstelem = 1;
	double nullval = 0;
	long npixels = 1;
	
	double *tempData;
	
		
	if (fits_open_file(&fptr, filename, READONLY, &status)) printerror(status);
	
	if (fits_read_keys_lng(fptr, "NAXIS", 1, naxis, npix, &nfound, &status)) printerror(status);
		   
	tempData = (double*) malloc(sizeof(double) * npix[0] * npix[1]);
	if (fits_read_keys_dbl(fptr, "CRPIX", 1, naxis, crpix, &dumint, &status)) printerror(status);
	if (fits_read_keys_dbl(fptr, "CRVAL", 1, naxis, crval, &dumint, &status)) printerror(status);
	if (fits_read_keys_dbl(fptr, "CDELT", 1, naxis, cdelt, &dumint, &status)) printerror(status);
	
	
	for (int i = 0; i < naxis; ++i) {
		npixels *= npix[i];
		crpix[i] -= 1; // 1-origin -> 0-origin
	}
	
	crval[0] *= 60.*60.;	//radian -> arcsec
	cdelt[0] *= 60.*60.;
	crval[1] *= 1.0e-3;		// m/s -> km/s
	cdelt[1] *= 1.0e-3;

	dataInit();
	
	if (fits_read_img(fptr, TDOUBLE, firstelem, npixels, &nullval, tempData, &dumint, &status)) printerror(status);
	for (int iv = 0; iv < npix[1]; ++iv) for (int ip = 0; ip < npix[0]; ++ip) valData[ip][iv] = tempData[iv * npix[0] + ip];
	free(tempData);
	
	if (fits_close_file(fptr, &status)) printerror(status);
	return ;
}


void PVdata::setMaxMin() {
	valMax = -INF;
	valMin = INF;
	double val;
	
	for (int ip = 0; ip < npix[0]; ++ip) for (int iv = 0; iv < npix[1]; ++iv) {
		val = valData[ip][iv];
		if (valMax < val) valMax = val;
		if (valMin > val) valMin = val;
	}
	return ;
}

void PVdata::getMaxMin(double *max, double *min) {
	setMaxMin();
	if (max) *max = valMax;
	if (min) *min = valMin;
	return ;
}


bool PVdata::skipCell(int ip, int iv, bool f_threshold = false, double threshold4skip = 0.) {
	if (ip < 0 || ip >= npix[0] || iv < 0 || iv >= npix[1]) return true;
	
	if (f_threshold) {
		if (abs(getVal(ip, iv)) < threshold4skip) return true;
		if (f_skipNegative && getVal(ip, iv) < threshold4skip) return true;
	}
	return false;
}


void PVdata::setVal(int ip, int iv, double val) {
	if (skipCell(ip, iv)) return ;
	valData[ip][iv] = val;
	return ;
}

double PVdata::getVal(int ip, int iv) {
	if (skipCell(ip, iv)) return noval;
	return valData[ip][iv];
}

double PVdata::getP(int ip) {
	return crval[0] + (ip - crpix[0]) * cdelt[0];
}

double PVdata::getIP(double p) {
	return crpix[0] + (p - crval[0]) / cdelt[0];
}

double PVdata::getV(int iv) {
	return crval[1] + (iv - crpix[1]) * cdelt[1];
}

double PVdata::getIV(double v) {
	return crpix[1] + (v - crval[1]) / cdelt[1];
}


//////////
class PVpair {
public:
	PVdata *PV1, *PV2, *PV2regrid;
	char filename1[lenfilename], filename2[lenfilename];
	
	double VsysOffset, weight4PV2regrid;
	
	double limPmin, limPmax, limVmin, limVmax, rPmin, rPmax, rVmin, rVmax;
	int rIPmin, rIPmax, rIVmin, rIVmax, nregion;
	
	
	PVpair();
	~PVpair();
	
	void init(PVdata *PV1_in, PVdata *PV2_in, PVdata *PV2_regrid, char* filename1_in, char* filename2_in, double VsysOffset_in);
	bool skipCell(int ip, int iv);
	bool skipCell_andThreshold(int ip, int iv, double threshold4skip_PV1_in, double threshold4skip_PV2regrid_in);
	bool skipCell_orThreshold(int ip, int iv, double threshold4skip_PV1_in, double threshold4skip_PV2regrid_in);
	bool skipCell_region4comparing(int ip, int iv, bool f_threshold, double threshold4skip_PV1_in, double threshold4skip_PV2regrid_in);
	void regrid();
	void setRegion(double limPmin_in, double limPmax_in, double limVmin_in, double limVmax_in);
	void setWeight(double weight_in);
};

PVpair::PVpair() {
}

PVpair::~PVpair() {
}

 
void PVpair::init(PVdata *PV1_in, PVdata *PV2_in, PVdata *PV2regrid_in, char* filename1_in, char* filename2_in, double VsysOffset_in = 0.) {
	PV1 = PV1_in;
	PV2 = PV2_in;
	PV2regrid = PV2regrid_in;
	
	sscanf(filename1_in, "%s", filename1);
	sscanf(filename2_in, "%s", filename2);
	VsysOffset = VsysOffset_in;
	
	PV1->readFits(filename1);
	PV2->readFits(filename2);
	regrid();
}

bool PVpair::skipCell(int ip, int iv) {
	return PV1->skipCell(ip, iv) || PV2regrid->skipCell(ip, iv);
}

bool PVpair::skipCell_andThreshold(int ip, int iv, double threshold4skip_PV1_in, double threshold4skip_PV2regrid_in) {
	skipCell(ip, iv);
	return PV1->skipCell(ip, iv, true, threshold4skip_PV1_in) && PV2regrid->skipCell(ip, iv, true, threshold4skip_PV2regrid_in);
}

bool PVpair::skipCell_orThreshold(int ip, int iv, double threshold4skip_PV1_in, double threshold4skip_PV2regrid_in) {
	skipCell(ip, iv);
	return PV1->skipCell(ip, iv, true, threshold4skip_PV1_in) || PV2regrid->skipCell(ip, iv, true, threshold4skip_PV2regrid_in);
}

bool PVpair::skipCell_region4comparing(int ip, int iv, bool f_threshold = false, double threshold4skip_PV1_in = 0., double threshold4skip_PV2regrid_in = 0.) {
	if (ip < rIPmin || ip >= rIPmax || iv < rIVmin || iv >= rIVmax) return true;
	
	if (f_threshold) {
		return skipCell_orThreshold(ip, iv, threshold4skip_PV1_in, threshold4skip_PV2regrid_in);
	} else {
		return skipCell(ip, iv);
	}
}


void PVpair::regrid() {
	PV2regrid->pixInit(PV1->npix, PV1->crpix, PV1->crval, PV1->cdelt);
	double p, ip_dbl, v, iv_dbl, dip, div;
	double val, valLL, valLU, valUL, valUU;
	int ip_low, ip_upp, iv_low, iv_upp;
	
	for (int ip = 0; ip < PV2regrid->npix[0]; ++ip) {
		p = PV2regrid->getP(ip);
		ip_dbl = PV2->getIP(p);
		ip_low = floor(ip_dbl * (1. + db_erratio));
		ip_upp = ceil(ip_dbl * (1. - db_erratio));
		
		if (ip_low < 0 || ip_upp >= PV2->npix[0]) continue;
		
		dip = ip_dbl - ip_low;
		
		
		for (int iv = 0; iv < PV2regrid->npix[1]; ++iv) {
			v = PV2regrid->getV(iv);
			iv_dbl = PV2->getIV(v - VsysOffset);
			iv_low = floor(iv_dbl * (1. + db_erratio));
			iv_upp = ceil(iv_dbl * (1. - db_erratio));
			
			if (iv_low < 0 || iv_upp >= PV2->npix[1]) continue;
			
			div = iv_dbl - iv_low;
			
			valLL = PV2->getVal(ip_low, iv_low);
			valLU = PV2->getVal(ip_low, iv_upp);
			valUL = PV2->getVal(ip_upp, iv_low);
			valUU = PV2->getVal(ip_upp, iv_upp);
			
			val = (valLL * (1. - div) + valLU * div) * (1. - dip) + (valUL * (1. - div) + valUU * div) * dip;
			
			PV2regrid->setVal(ip, iv, val);
			
		}
	}
	

	PV2regrid->setMaxMin();
		
	return ;
}

void PVpair::setRegion(double limPmin_in, double limPmax_in, double limVmin_in, double limVmax_in) {
	limPmin = limPmin_in;
	limPmax = limPmax_in;
	limVmin = limVmin_in;
	limVmax = limVmax_in;
	
	
	rIPmin = floor(PV1->getIP(limPmin));
	rIPmax = ceil(PV1->getIP(limPmax));
	rIVmin = floor(PV1->getIV(limVmin));
	rIVmax = ceil(PV1->getIV(limVmax));
	
	if (rIPmin < 0) rIPmin = 0;
	if (rIPmax >= PV1->npix[0]) rIPmax = PV1->npix[0] - 1;
	if (rIVmin < 0) rIVmin = 0;
	if (rIVmax >= PV1->npix[1]) rIVmax = PV1->npix[1] - 1;
	
	nregion = (rIPmax - rIPmin + 1) * (rIVmax - rIVmin + 1);
	
	rPmin = PV1->getP(rIPmin);
	rPmax = PV1->getP(rIPmax);
	rVmin = PV1->getV(rIVmin);
	rVmax = PV1->getV(rIVmax);

	fprintf(fp_out, "!\tRegion in PV 1: ip=[0, %ld]\t<-> p = [%.4f, %.4f],\tiv=[0, %ld]\t<-> v = [%.4f, %.4f]\n", PV1->npix[0] - 1, PV1->getP(0), PV1->getP(PV1->npix[0] - 1), PV1->npix[1] - 1, PV1->getV(0), PV1->getV(PV1->npix[1] - 1));
	fprintf(fp_out, "!\tRegion for Comparison: ip=[%d, %d]\t<-> p = [%.4f, %.4f],\tiv=[%d, %d]\t<-> v=[%.4f, %.4f], #Cell = %d\n", rIPmin, rIPmax, rPmin, rPmax, rIVmin, rIVmax, rVmin, rVmax, nregion);
	
	return ;
}

void PVpair::setWeight(double weight_in) {
	weight4PV2regrid = weight_in;
	
	for (int ip = 0; ip < PV2regrid->npix[0]; ++ip) for (int iv = 0; iv < PV2regrid->npix[1]; ++iv) {
		PV2regrid->setVal(ip, iv, PV2regrid->getVal(ip, iv) * weight4PV2regrid);
	}
	
	PV2regrid->setMaxMin();
	return ;
}

//////////
class PVpairSet {
public:
	PVpair *PVpairList[nPairsLim];
	int nPairs, nPairs4weight, pairs4weight[nPairsLim];
	double weight4PV2regrid, threshold_PV1, VsysOffset;
	char paramline4weight[lenline];
	
	PVpairSet(int nPairs_in, char *paramline4weight_in, double threshold_PV1_in, double VsysOffset_in);
	~PVpairSet();
	
	void initPair(int ipair, PVpair *PVpair_in);
	void calcRatio();
	void setRatio();
	void compare(double rms_PV1);
};

PVpairSet::PVpairSet(int nPairs_in, char *paramline4weight_in, double threshold_PV1_in, double VsysOffset_in = 0.) {
	nPairs = nPairs_in;
	strcpy(paramline4weight, paramline4weight_in);
	threshold_PV1 = threshold_PV1_in;
	VsysOffset = VsysOffset_in;
}

PVpairSet::~PVpairSet() {
}

void PVpairSet::initPair(int ipair, PVpair *PVpair_in) {
	PVpairList[ipair] = PVpair_in;
	return ;
}



void PVpairSet::calcRatio() {
	double val_PV1, val_PV2regrid;
	int ndatapoints4weight = 0;
	
	double sxx = 0., sxy = 0., sx = 0., sy = 0.;
	
	for (int i = 0; i < nPairs4weight; ++i) {
		int ipair = pairs4weight[i];
		for (int ip = PVpairList[ipair]->rIPmin; ip <= PVpairList[ipair]->rIPmax; ++ip) for (int iv = PVpairList[ipair]->rIVmin; iv <= PVpairList[ipair]->rIVmax; ++iv) {
			val_PV1 = PVpairList[ipair]->PV1->getVal(ip, iv);
			val_PV2regrid = PVpairList[ipair]->PV2regrid->getVal(ip, iv);
			
			
			if (PVpairList[ipair]->skipCell_orThreshold(ip, iv, threshold_PV1, threshold_PV1 / weight4PV2regrid)) continue;
			
			++ndatapoints4weight;
			sxx += val_PV1 * val_PV1;
			sxy += val_PV1 * val_PV2regrid;
			sx += val_PV1;
			sy += val_PV2regrid;
		}
	}
	
	if (ndatapoints4weight == 0) {
		printf("\n\tERROR :: There is no data point for calculating the normalization constant.  Check the threshold level. \n\n");
		exit(1);
	}
	
	weight4PV2regrid =  sxx / sxy;
	
	
	if (abs(weight4PV2regrid) < EPS) {
		printf("\n\tWARN :: The normalization constant is 0. \n");
		printf("\t\t\tweight = %le\n\n", weight4PV2regrid);
		exit(1);
	}
	
	if (weight4PV2regrid < 0.) {
		printf("\n\tWARN :: The normalization constant is negative. \n");
		printf("\t\t\tweight = %le\n", weight4PV2regrid);
	}
	
	return ;
}


void PVpairSet::setRatio() {
	char *paramline4weight_delim;
	double weight4PV2regrid_in;
	
	paramline4weight_delim = strtok(paramline4weight, delim);
	sscanf(paramline4weight_delim, "%d", &nPairs4weight);

	if (nPairs4weight == 0) {
		paramline4weight_delim = strtok(NULL, delim);
		sscanf(paramline4weight_delim, "%lf", &weight4PV2regrid);
		return ;
	}
	
	for (int i = 0; i < nPairs4weight; ++i) {
		paramline4weight_delim = strtok(NULL, delim);
		sscanf(paramline4weight_delim, "%d", &pairs4weight[i]);
	}
	
	
	double maxval_PV1, maxval_PV2regrid, dumdbl;
	PVpairList[pairs4weight[0]]->PV1->getMaxMin(&maxval_PV1, &dumdbl);
	PVpairList[pairs4weight[0]]->PV2regrid->getMaxMin(&maxval_PV2regrid, &dumdbl);
	weight4PV2regrid = maxval_PV1 / maxval_PV2regrid;
	
	weight4PV2regrid_in = weight4PV2regrid;
	for (int i = 0; i < calcRatioLim; ++i) {
		calcRatio();
		if (abs(weight4PV2regrid - weight4PV2regrid_in) < EPS) break;
		weight4PV2regrid_in = weight4PV2regrid;
	}
	
	return ;
}


void PVpairSet::compare(double rms_PV1) {
	int ncount_region[nPairs], ncount_region_allPairs = 0, ncount_region_chi2[nPairs], ncount_region_chi2_allPairs = 0, ncount_threshold[nPairs], ncount_threshold_allPairs = 0;
	int npix_max[2] = {};
	
	double sumres_region[nPairs], sumres_region_allPairs = 0., sumres_threshold[nPairs], sumres_threshold_allPairs = 0.;
	double sumres_region_chi2[nPairs], sumres_region_chi2_allPairs = 0., sumres_threshold_chi2[nPairs], sumres_threshold_chi2_allPairs = 0.;
	double p, v, val1, val2regrid;
	
	
	setRatio();
	
	for (int i = 0; i < nPairs; ++i) {
		PVpairList[i]->setWeight(weight4PV2regrid);
		for (int j = 0; j < 2; ++j) if (npix_max[j] < PVpairList[i]->PV1->npix[j]) npix_max[j] = PVpairList[i]->PV1->npix[j];

		ncount_region[i] = 0;
		ncount_region_chi2[i] = 0;
		ncount_threshold[i] = 0;
		sumres_region[i] = 0.;
		sumres_threshold[i] = 0.;
		sumres_region_chi2[i] = 0.;
		sumres_threshold_chi2[i] = 0.;
	}
	
	printf("Normalization constant for regridded PV data = %.4e\n", weight4PV2regrid);
	fprintf(fp_out, "!Normalization constant for regridded PV data = %.4e\n", weight4PV2regrid);
	fprintf(fp_out, "!\n");
	
    fprintf(fp_out, "!  \t  ");
    for (int i = 0; i < nPairs; ++i) fprintf(fp_out, "\t|Pair%d\t\t\t\t\t\t\t\t\t\t\t\t", i);
    fprintf(fp_out, "\n");
    fprintf(fp_out, "!ip\tiv");
    for (int i = 0; i < nPairs; ++i) fprintf(fp_out, "\t|used\tp\t\t\tv\t\t\tPV1\t\tPV2regrid\t");
    fprintf(fp_out, "\n");
    
	for (int ip = 0; ip < npix_max[0]; ++ip) for (int iv = 0; iv < npix_max[1]; ++iv) {
		fprintf(fp_out, "%d\t%d", ip, iv);
		for (int ipair = 0; ipair < nPairs; ++ipair) {
			p = PVpairList[ipair]->PV1->getP(ip);
			v = PVpairList[ipair]->PV1->getV(iv);
			
			if (PVpairList[ipair]->skipCell(ip, iv)) {
				fprintf(fp_out, "\t0\t%8.4f\t%8.4f\t%.4e\t%.4e\t", p, v, noval, noval);
				continue;
			}
			
			val1 = PVpairList[ipair]->PV1->getVal(ip, iv);
			val2regrid = PVpairList[ipair]->PV2regrid->getVal(ip, iv);
			
			
			if (PVpairList[ipair]->skipCell_region4comparing(ip, iv)) {
				fprintf(fp_out, "\t0\t%8.4f\t%8.4f\t%.4e\t%.4e\t", p, v, val1, val2regrid);
				continue;
			}

			if (PVpairList[ipair]->skipCell_region4comparing(ip, iv, true, threshold_PV1, threshold_PV1)) {
				fprintf(fp_out, "\t1\t%8.4f\t%8.4f\t%.4e\t%.4e\t", p, v, val1, val2regrid);
				sumres_region[ipair] += (val1 - val2regrid) * (val1 - val2regrid);
				++ncount_region[ipair];
				
				if (abs(val2regrid) > EPS) {
					++ncount_region_chi2[ipair];
					sumres_region_chi2[ipair] += (val1 - val2regrid) * (val1 - val2regrid) / abs(val2regrid);
				}
				continue;
			}
			
			++ncount_threshold[ipair];
			fprintf(fp_out, "\t2\t%8.4f\t%8.4f\t%.4e\t%.4e\t", p, v, val1, val2regrid);
			
			sumres_threshold[ipair] += (val1 - val2regrid) * (val1 - val2regrid);
			sumres_threshold_chi2[ipair] += (val1 - val2regrid) * (val1 - val2regrid) / abs(val2regrid);
		}
		fprintf(fp_out, "\n");
	}
	
	for (int i = 0; i < nPairs; ++i) {
		ncount_region[i] += ncount_threshold[i];
		ncount_region_chi2[i] += ncount_threshold[i];
		sumres_region[i] += sumres_threshold[i];
		sumres_region_chi2[i] += sumres_threshold_chi2[i];
		
		ncount_region_allPairs += ncount_region[i];
		ncount_region_chi2_allPairs += ncount_region_chi2[i];
		ncount_threshold_allPairs += ncount_threshold[i];
		
		sumres_region_allPairs += sumres_region[i];
		sumres_region_chi2_allPairs += sumres_region_chi2[i];
		sumres_threshold_allPairs += sumres_threshold[i];
		sumres_threshold_chi2_allPairs += sumres_threshold_chi2[i];
	}
	
	
	
    if (logger) printf("\n                                                          \t");
    fprintf(fp_out, "!\n!                                                          \t");
    for (int i = 0; i < nPairs; ++i) {
        if (logger) printf("|Pair%2d\t\t", i);
        fprintf(fp_out, "|Pair%2d\t\t\t", i);
    }
    if (logger) printf("\n#Cell in Region                                           \t");
    fprintf(fp_out, "\n!#Cell in Region                                            \t");
    for (int i = 0; i < nPairs; ++i) {
        if (logger) printf("|%15d", ncount_region[i]);
        fprintf(fp_out, "|%15d", ncount_region[i]);
    }
    if (logger) printf("\n#Cell in Region > threshold                               \t");
    fprintf(fp_out, "\n!#Cell in Region > threshold                                \t");
    for (int i = 0; i < nPairs; ++i) {
        if (logger) printf("|%15d", ncount_threshold[i]);
        fprintf(fp_out, "|%15d", ncount_threshold[i]);
    }
    if (logger) printf("\nSum of the squared residual (in region)                   \t");
    fprintf(fp_out, "\n!Sum of the squared residual (in region)                    \t");
    for (int i = 0; i < nPairs; ++i) {
        if (logger) printf("|%15.4e", sumres_region[i]);
        fprintf(fp_out, "|%15.4e", sumres_region[i]);
    }
    if (logger) printf("\nSum of the squared residual (> threshold)                 \t");
    fprintf(fp_out, "\n!Sum of the squared residual (> threshold)                  \t");
    for (int i = 0; i < nPairs; ++i) {
        if (logger) printf("|%15.4e", sumres_threshold[i]);
        fprintf(fp_out, "|%15.4e", sumres_threshold[i]);
    }
	if (logger) printf("\n#Cell in Region (nonZero in PV2)                           \t");
	fprintf(fp_out, "\n!#Cell in Region (nonZero in PV2)                           \t");
	for (int i = 0; i < nPairs; ++i) {
		if (logger) printf("|%15d", ncount_region_chi2[i]);
		fprintf(fp_out, "|%15d", ncount_region_chi2[i]);
	}
	if (logger) printf("\n#Cell > threshold                                          \t");
	fprintf(fp_out, "\n!#Cell > threshold                                          \t");
	for (int i = 0; i < nPairs; ++i) {
		if (logger) printf("|%15d", ncount_threshold[i]);
		fprintf(fp_out, "|%15d", ncount_threshold[i]);
	}
	
    if (logger) printf("\nSum of (PV1 - PV2)^2 / rms^2 (in region)                   \t");
    fprintf(fp_out, "\n!Sum of (PV1 - PV2)^2 / rms^2 (in region)                   \t");
    for (int i = 0; i < nPairs; ++i) {
        if (logger) printf("|%15.4e", sumres_region[i] / threshold_PV1);
        fprintf(fp_out, "|%15.4e", sumres_region[i] / threshold_PV1);
    }
    if (logger) printf("\nSum of (PV1 - PV2)^2 / rms^2 (> threshold)                 \t");
    fprintf(fp_out, "\n!Sum of (PV1 - PV2)^2 / rms^2 (> threshold)                  \t");
    for (int i = 0; i < nPairs; ++i) {
        if (logger) printf("|%15.4e", sumres_threshold[i] / threshold_PV1);
        fprintf(fp_out, "|%15.4e", sumres_threshold[i] / threshold_PV1);
    }
	/*
	if (logger) printf("\nSum of (PV1 - PV2)^2 / PV2 (in region && nonZero in PV2)   \t");
	fprintf(fp_out, "\n!Sum of (PV1 - PV2)^2 / PV2 (in region && nonZero in PV2)   \t");
	for (int i = 0; i < nPairs; ++i) {
		if (logger) printf("|%15.4e", sumres_region_chi2[i]);
		fprintf(fp_out, "|%15.4e", sumres_region_chi2[i]);
	}
	if (logger) printf("\nSum of (PV1 - PV2)^2 / PV2 (> threshold)                   \t");
	fprintf(fp_out, "\n!Sum of (PV1 - PV2)^2 / PV2 (> threshold)                   \t");
	for (int i = 0; i < nPairs; ++i) {
		if (logger) printf("|%15.4e", sumres_threshold_chi2[i]);
		fprintf(fp_out, "|%15.4e", sumres_threshold_chi2[i]);
	}
	//*/
	
	if (true) printf("\nSum of (PV1 - PV2)^2 / rms^2 (in region)   for all Pairs \t|%15.4e| #Cell = %d \t-> Chi^2 / DOF = %lf", sumres_region_allPairs / rms_PV1 / rms_PV1, ncount_region_allPairs, sumres_region_allPairs / rms_PV1 / rms_PV1 / ncount_region_allPairs);
	fprintf(fp_out, "\n!Sum of (PV1 - PV2)^2 / rms^2 (in region)   for all Pairs \t|%15.4e| #Cell = %d \t-> Chi^2 / DOF = %lf", sumres_region_allPairs / rms_PV1 / rms_PV1, ncount_region_allPairs, sumres_region_allPairs / rms_PV1 / rms_PV1 / ncount_region_allPairs);
	if (logger) printf("\nSum of (PV1 - PV2)^2 / rms^2 (> threshold) for all Pairs \t|%15.4e| #Cell = %d \t-> Chi^2 / DOF = %lf", sumres_threshold_allPairs / rms_PV1 / rms_PV1, ncount_threshold_allPairs, sumres_threshold_allPairs / rms_PV1 / rms_PV1 / ncount_region_allPairs);
	fprintf(fp_out, "\n!Sum of (PV1 - PV2)^2 / rms^2 (> threshold) for all Pairs \t|%15.4e| #Cell = %d \t-> Chi^2 / DOF = %lf", sumres_threshold_allPairs / rms_PV1 / rms_PV1, ncount_threshold_allPairs, sumres_threshold_allPairs / rms_PV1 / rms_PV1 / ncount_region_allPairs);

	/*
	if (logger) printf("\nSum of (PV1 - PV2)^2 / PV2^2 (in region)   for all Pairs \t|%15.4e| #Cell = %d \t-> Chi^2 / DOF = %lf", sumres_region_chi2_allPairs, ncount_region_chi2_allPairs, sumres_region_chi2_allPairs / ncount_region_chi2_allPairs);
	fprintf(fp_out, "\n!Sum of (PV1 - PV2)^2 / PV2^2 (in region)   for all Pairs \t|%15.4e| #Cell = %d \t-> Chi^2 / DOF = %lf", sumres_region_chi2_allPairs, ncount_region_chi2_allPairs, sumres_region_chi2_allPairs / ncount_region_chi2_allPairs);
	if (logger) printf("\nSum of (PV1 - PV2)^2 / PV2^2 (> threshold) for all Pairs \t|%15.4e| #Cell = %d \t-> Chi^2 / DOF = %lf", sumres_threshold_chi2_allPairs, ncount_threshold_allPairs, sumres_threshold_chi2_allPairs / ncount_threshold_allPairs);
	fprintf(fp_out, "\n!Sum of (PV1 - PV2)^2 / PV2^2 (> threshold) for all Pairs \t|%15.4e| #Cell = %d \t-> Chi^2 / DOF = %lf", sumres_threshold_chi2_allPairs, ncount_threshold_allPairs, sumres_threshold_chi2_allPairs / ncount_threshold_allPairs);
	//*/

    printf("\n");
    fprintf(fp_out, "\n!\n");
    

	
	return ;
}



//////////


int main() {
    
    int nPairs;
	char paramline[lenline];
	
    while (~scanf("%d", &nPairs)) {
		fgets(paramline, lenline, stdin);
		
        if (nPairs < 1) break;
        puts("---------------------------------------------");
        char outfilename[lenfilename], PVfilename1[lenfilename], PVfilename2[lenfilename], tempfilename[lenfilename];
        double Vsys, rms_PV1, threshold, limPmin, limPmax, limVmin, limVmax;
		
		fgets(paramline, lenline, stdin);
		sscanf(paramline, "%lf ", &Vsys);
		
        fgets(paramline, lenline, stdin);
        sscanf(paramline, "%lf %lf %lf %lf %lf %lf ", &rms_PV1, &threshold, &limPmin, &limPmax, &limVmin, &limVmax);
		
        fgets(paramline, lenline, stdin);
		
		
		fgets(tempfilename, lenline, stdin);
		sscanf(tempfilename, "%s", tempfilename);
        sscanf(OutDir, "%s", outfilename);
        strcat(outfilename, tempfilename);
        
        fp_out = fopen(outfilename, "w");
        
        fprintf(fp_out, "!Number of Pairs   : %d\n", nPairs);
        fprintf(fp_out, "!Systemic Velocity: %.4f km/s\n", Vsys);
        fprintf(fp_out, "!rms in PV1: %.4f km/s\n", rms_PV1);
        fprintf(fp_out, "!Region to compare: threshold = %lf Jy/beam, xrange=[%lf, %lf], yrange=[%lf, %lf]\n", threshold, limPmin, limPmax, limVmin, limVmax);
        fprintf(fp_out, "!Other Parameters: %s!---------------\n", paramline);
		

		PVdata PVdataList[nPairs][3];
		PVpair PVpairList[nPairs];
		PVpairSet PVpairSet_in(nPairs, paramline, threshold, Vsys);
		
        for (int i = 0; i < nPairs; ++i) {
			fgets(tempfilename, lenline, stdin);
			sscanf(tempfilename, "%s", tempfilename);
            sscanf(DataDir, "%s", PVfilename1);
            strcat(PVfilename1, tempfilename);
			
			fgets(tempfilename, lenline, stdin);
			sscanf(tempfilename, "%s", tempfilename);
            sscanf(DataDir, "%s", PVfilename2);
            strcat(PVfilename2, tempfilename);
			
			PVpairList[i].init(&PVdataList[i][0], &PVdataList[i][1], &PVdataList[i][2], PVfilename1, PVfilename2);
			PVpairSet_in.initPair(i, &PVpairList[i]);
		}
		
        for (int i = 0; i < nPairs; ++i) {
            fprintf(fp_out, "!Pair %d----------\n", i);
			
            fprintf(fp_out, "!\tPV 1: %s\n", PVfilename1);
            fprintf(fp_out, "!\t          p: crpix=%.4e, crval=%.4e arcsec, cdelt=%.4e arcsec\n", PVpairSet_in.PVpairList[i]->PV1->crpix[0], PVpairSet_in.PVpairList[i]->PV1->crval[0], PVpairSet_in.PVpairList[i]->PV1->cdelt[0]);
            fprintf(fp_out, "!\t          v: crpix=%.4e, crval=%.4e km/s, cdelt=%.4e km/s\n", PVpairSet_in.PVpairList[i]->PV1->crpix[1], PVpairSet_in.PVpairList[i]->PV1->crval[1], PVpairSet_in.PVpairList[i]->PV1->cdelt[1]);

            fprintf(fp_out, "!\tPV 2: %s\n", PVfilename2);
            fprintf(fp_out, "!\t          p: crpix=%.4e, crval=%.4e arcsec, cdelt=%.4e arcsec\n", PVpairSet_in.PVpairList[i]->PV2->crpix[0], PVpairSet_in.PVpairList[i]->PV2->crval[0], PVpairSet_in.PVpairList[i]->PV2->cdelt[0]);
            fprintf(fp_out, "!\t          v: crpix=%.4e, crval=%.4e km/s, cdelt=%.4e km/s\n\n", PVpairSet_in.PVpairList[i]->PV2->crpix[1], PVpairSet_in.PVpairList[i]->PV2->crval[1], PVpairSet_in.PVpairList[i]->PV2->cdelt[1]);

			PVpairList[i].setRegion(limPmin, limPmax, limVmin, limVmax);

            if (logger) printf("\nPair (%d/%d)\n", i + 1, nPairs);
            if (logger) printf("\tPV 1: %s\n", PVfilename1);
            if (logger) printf("\tPV 2: %s\n", PVfilename2);
        }
        fprintf(fp_out, "!--------------\n");
        
		
		PVpairSet_in.compare(rms_PV1);
        
        fclose(fp_out);
        printf("\nSaved: %s\n\n", outfilename);
		
    }
	
    
    return 0;
}


