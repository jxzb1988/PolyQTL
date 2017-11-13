#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <math.h>
#include <unistd.h> 
#include <armadillo>
#ifndef _UTIL_H
#define _UTIL_H
#include "Util.h"
#endif

#ifndef _POSTCAL_H
#define _POSTCAL_H
#include "PostCal.h"
#endif

#ifndef _TOPKSNP_H
#define _TOPKSNP_H
#include "TopKSNP.h"
#endif

#ifndef _MEQTLPOLYGMODEL_H
#define _MEQTLPOLYGMODEL_H
#include "MeQTLPolyGModel.h"
#endif

#include "conditional_function.h"

#include <sys/stat.h>

using namespace std;


int main( int argc, char *argv[]  ){
        
        
        if (argc <= 1) {
                cout<<"Error, parameter(s) needed"<<endl;
                return EXIT_SUCCESS;
        }
        if (argc==2 && argv[1][0] == '-' && argv[1][1] == 'c') {
                
                return EXIT_SUCCESS;
        }
        int totalCausalSNP = 1;
	double gamma = 0.01;
	float rho = 0.95;
	bool histFlag = false;
	int oc = 0;	
	string yFile  = "";
	string outputFileName = "";
	string geneMapFile = "";	
        string weight = "";
        string covariate = "";
        string grm_test ="";
        string gene="";
        vector<string> grm_file;
        string X_file="";
        string geno_file="";        
        int number;
        int nthread=1;
        int mode=-1;
        string input="";
	while ((oc = getopt(argc, argv, "vhl:t:o:x:Z:p:n:w:g:r:c:G:w:f:m:P:T:")) != -1) {
		switch (oc) {
			case 'v':
				cout << "version 1.0:" << endl;
			case 'h':
				cout << "Options: " << endl;
  				cout << "-h, --help            		show this help message and exit " << endl;
  				cout << "-o OUTFILE, --out=OUTFILE 	specify the output file" << endl;
  				cout << "-l LDFILE, --ld_file=LDFILE  	the ld input file" << endl;
  				cout << "-p yFILE, --y_file=yFILE	phenotype" << endl;
  				cout << "-r RHO, --rho-prob=RHO		set $pho$ probability (default 0.95)" << endl;
				cout << "-g GAMMA, --gamma		set $gamma$ the prior of a SNP being causal (default 0.01)" << endl;
				cout << "-c causal			set the maximum number of causal SNPs" << endl;
                                cout << "-C covariate, --covariate      set the covariate matrix "<<endl;
                                cout << "-t gene whose abundance, --target      gene name "<<endl;
                                cout << "-G genetic relatedness matrix, --GRM  set the genetic relatedness matrix "<<endl;
				cout << "-f 1				to out the probaility of different number of causal SNP" << endl;
                                cout << "-w Weight file, --weight       set the biological annotation to use" << endl;
                                cout << "-n Number of samples, --number       set the biological annotation to use" << endl;
                                cout << "-x genotype information,       genotype file for explored variants" << endl;
                                cout << "-t threads to use, --nthread       set the threads to use" << endl;
				exit(0);
                        case 'n':
                                number = atoi(optarg);
                                break;
//                        case 'm':
  //                              conditional = atoi(optarg);
    //                            break;
			case 'o':
				outputFileName = string(optarg);
				break;
                        case 'T':
                                gene           =string(optarg);
                                break;
                        case 'P':
                                input          = string(optarg);
                                break;
			case 'p':
				yFile = string(optarg);
				break;
			case 'r':
				rho = atof(optarg);
				break;
                        case 'x':
                                X_file = string(optarg);
                                break;
			case 'c':
				totalCausalSNP = atoi(optarg);
				break;
			case 'g':
				gamma = atof(optarg);
				break;
                        case 'w':
                                weight=string(optarg);
                                break;
                        case 't':
                                nthread=atoi(optarg);
			case 'f':
                                histFlag = true;
                                break;
                        case 'G':
                                grm_test     =string(optarg);
                                grm_file.push_back(grm_test);
                                break;
                        case 'C':
                                covariate=string(optarg);
                                break;
                        case 'Z':
                                geno_file=string(optarg);
                                break;
			case ':':
			case '?':
			default:
				cout << "Strange" << endl;
				break;
		}
	}
        cout<<"Getting parameter information is over"<<endl;
        cout<<"input is: "<<input<<endl;
        cout<<"gene is: "<<gene<<endl;
        
        ifstream check_dir(gene);
        if (!check_dir) {
                mkdir(gene.c_str(), S_IRWXU|S_IRGRP|S_IROTH);
        }

        ifstream check_dir2("output/");
        if (!check_dir2) {
                mkdir("output", S_IRWXU|S_IRGRP|S_IROTH);
        }


//        return 0;
        if(input!="")
         {
           cout<<"input is: "<<input<<endl;
           cout<<"gene is: "<<gene<<endl;
           cout<<"Perform conditional analysis"<<endl;
//           conditional_analysis(input,gene, totalCausalSNP,rho, histFlag, gamma,weight, nthread,covariate, grm_file);
           conditional_analysis(input,gene, totalCausalSNP,rho, histFlag, gamma,weight, nthread,covariate, grm_file,outputFileName, geno_file);                     
           //  conditional_analysis(input,gene, totalCausalSNP,rho, histFlag);
           return 0;
         } else
         {
           MeQTLPolyGModel  MeQTLPoly(yFile, outputFileName, totalCausalSNP,rho, histFlag, gamma,weight, nthread,covariate, grm_file,X_file);
	   MeQTLPoly.run();
	   MeQTLPoly.finishUp();		
	   return 0;
         }
}
