#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sys/stat.h>
#include <cmath>

#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_eigen.h"
#include "gsl/gsl_cdf.h"

#include "gemma_lapack.h"  //for functions EigenDecomp

#ifdef FORCE_FLOAT
#include "io_float.h"   //for function ReadFile_kin
#include "gemma_float.h"
#include "lmm_float.h"  //for LMM class, and functions CalcLambda, CalcPve, CalcVgVe
#include "mathfunc_float.h"	//for a few functions
#else
#include "gemma_io.h"
#include "gemma.h"
#include "gemma_lmm.h"
#include "gemma_mathfunc.h"
#endif
using namespace std;
GEMMA::GEMMA(void):	
version("0.94"), date("01/12/2014"), year("2011")
{}
void GEMMA::PrintHeader (void)
{
}


void GEMMA::PrintLicense (void)
{
}



void GEMMA::PrintHelp(size_t option)
{
	if (option==0) {
		cout<<endl; 
		cout<<" GEMMA version "<<version<<", released on "<<date<<endl;
		cout<<" implemented by Xiang Zhou"<<endl; 
		cout<<endl;
		cout<<" type ./gemma -h [num] for detailed helps"<<endl;
		cout<<" options: " << endl;
		cout<<" 1: quick guide"<<endl;
		cout<<" 2: file I/O related"<<endl;
		cout<<" 3: SNP QC"<<endl;
		cout<<" 4: calculate relatedness matrix"<<endl;
		cout<<" 5: perform eigen decomposition"<<endl;
		cout<<" 6: fit a linear model"<<endl;
		cout<<" 7: fit a linear mixed model"<<endl;
		cout<<" 8: fit a multivariate linear mixed model"<<endl;
		cout<<" 9: fit a Bayesian sparse linear mixed model"<<endl;
		cout<<" 10: obtain predicted values"<<endl;
		cout<<" 11: note"<<endl;
		cout<<endl;
	}	
	
	if (option==1) {
		cout<<" QUICK GUIDE" << endl;
		cout<<" to generate a relatedness matrix: "<<endl;
		cout<<"         ./gemma -bfile [prefix] -gk [num] -o [prefix]"<<endl;
		cout<<"         ./gemma -g [filename] -p [filename] -gk [num] -o [prefix]"<<endl;
		cout<<" to perform eigen decomposition of the relatedness matrix: "<<endl;
		cout<<"         ./gemma -bfile [prefix] -k [filename] -eigen -o [prefix]"<<endl;
		cout<<"         ./gemma -g [filename] -p [filename] -k [filename] -eigen -o [prefix]"<<endl;
		cout<<" to fit a linear mixed model: "<<endl;
		cout<<"         ./gemma -bfile [prefix] -k [filename] -lmm [num] -o [prefix]"<<endl;
		cout<<"         ./gemma -g [filename] -p [filename] -a [filename] -k [filename] -lmm [num] -o [prefix]"<<endl;	
		cout<<" to fit a multivariate linear mixed model: "<<endl;
		cout<<"         ./gemma -bfile [prefix] -k [filename] -lmm [num] -n [num1] [num2] -o [prefix]"<<endl;
		cout<<"         ./gemma -g [filename] -p [filename] -a [filename] -k [filename] -lmm [num] -n [num1] [num2] -o [prefix]"<<endl;	
		cout<<" to fit a Bayesian sparse linear mixed model: "<<endl;
		cout<<"         ./gemma -bfile [prefix] -bslmm [num] -o [prefix]"<<endl;
		cout<<"         ./gemma -g [filename] -p [filename] -a [filename] -bslmm [num] -o [prefix]"<<endl;
		cout<<" to obtain predicted values: "<<endl;
		cout<<"         ./gemma -bfile [prefix] -epm [filename] -emu [filename] -ebv [filename] -k [filename] -predict [num] -o [prefix]"<<endl;
		cout<<"         ./gemma -g [filename] -p [filename] -epm [filename] -emu [filename] -ebv [filename] -k [filename] -predict [num] -o [prefix]"<<endl;
		cout<<endl;
	}
	
	
	
	if (option==7) {
		cout<<" LINEAR MIXED MODEL OPTIONS" << endl;		
		cout<<" -lmm      [num]         "<<" specify analysis options (default 1)."<<endl;
		cout<<"          options: 1: Wald test"<<endl;		
		cout<<"                   2: Likelihood ratio test"<<endl;
		cout<<"                   3: Score test"<<endl;
		cout<<"                   4: 1-3"<<endl;
		cout<<"                   5: Parameter estimation in the null model only"<<endl;
		cout<<" -lmin     [num]          "<<" specify minimal value for lambda (default 1e-5)" << endl; 
		cout<<" -lmax     [num]          "<<" specify maximum value for lambda (default 1e+5)" << endl; 
		cout<<" -region   [num]          "<<" specify the number of regions used to evaluate lambda (default 10)" << endl; 
		cout<<endl;
	}
	
	
	if (option==11) {
		cout<<" NOTE"<<endl;
		cout<<" 1. Only individuals with non-missing phenotoypes and covariates will be analyzed."<<endl;
		cout<<" 2. Missing genotoypes will be repalced with the mean genotype of that SNP."<<endl;
		cout<<" 3. For lmm analysis, memory should be large enough to hold the relatedness matrix and to perform eigen decomposition."<<endl;
		cout<<" 4. For multivariate lmm analysis, use a large -pnr for each snp will increase computation time dramatically."<<endl;
		cout<<" 5. For bslmm analysis, in addition to 3, memory should be large enough to hold the whole genotype matrix."<<endl;
		cout<<endl;
	}
	
	return;
}



void GEMMA::Assign(int argc, char ** argv, PARAM &cPar)
{
	string str;
	
	for(int i = 1; i < argc; i++) {		
		if (strcmp(argv[i], "-bfile")==0 || strcmp(argv[i], "--bfile")==0 || strcmp(argv[i], "-b")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.file_bfile=str;
		}
		else if (strcmp(argv[i], "-silence")==0) {
			cPar.mode_silence=true;
		}
		else if (strcmp(argv[i], "-g")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.file_geno=str;
		}
		else if (strcmp(argv[i], "-p")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.file_pheno=str;
		}
		else if (strcmp(argv[i], "-a")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.file_anno=str;
		}
		else if (strcmp(argv[i], "-k")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.file_kin=str;
		}
		else if (strcmp(argv[i], "-c")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.file_cvt=str;
		}
		else if (strcmp(argv[i], "-r")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.file_read=str;
		}
		else if (strcmp(argv[i], "-snps")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.file_snps=str;
		}
		else if (strcmp(argv[i], "-km")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.k_mode=atoi(str.c_str());
		}		
		else if (strcmp(argv[i], "-n")==0) {
			(cPar.p_column).clear();
			while (argv[i+1] != NULL && argv[i+1][0] != '-') {
				++i;
				str.clear();
				str.assign(argv[i]);
				(cPar.p_column).push_back(atoi(str.c_str()));
			}
		}
		else if (strcmp(argv[i], "-pace")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.d_pace=atoi(str.c_str());
		}
		else if (strcmp(argv[i], "-o")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.file_out=str;
		}		
		else if (strcmp(argv[i], "-miss")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.miss_level=atof(str.c_str());
		}
		else if (strcmp(argv[i], "-maf")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			if (cPar.maf_level!=-1) {cPar.maf_level=atof(str.c_str());}
		}
		else if (strcmp(argv[i], "-hwe")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.hwe_level=atof(str.c_str());
		}
		else if (strcmp(argv[i], "-r2")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.r2_level=atof(str.c_str());
		}
		else if (strcmp(argv[i], "-notsnp")==0) {
			cPar.maf_level=-1;
		}
		else if (strcmp(argv[i], "-gk")==0) {
			if (cPar.a_mode!=0) {cPar.error=true; cout<<"error! only one of -gk -eigen -lm -lmm -bslmm -predict options is allowed."<<endl; break;}
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {cPar.a_mode=21; continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.a_mode=20+atoi(str.c_str());
		}	
		else if (strcmp(argv[i], "-fa")==0 || strcmp(argv[i], "-lmm")==0) {
			if (cPar.a_mode!=0) {cPar.error=true; cout<<"error! only one of -gk -eigen -lm -lmm -bslmm -predict options is allowed."<<endl; break;}
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {cPar.a_mode=1; continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.a_mode=atoi(str.c_str());
		}
		else if (strcmp(argv[i], "-lmin")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.l_min=atof(str.c_str());
		}
		else if (strcmp(argv[i], "-lmax")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.l_max=atof(str.c_str());
		}
		else if (strcmp(argv[i], "-region")==0) {
			if(argv[i+1] == NULL || argv[i+1][0] == '-') {continue;}
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.n_region=atoi(str.c_str());
		}
		else {cout<<"error! unrecognized option: "<<argv[i]<<endl; cPar.error=true; continue;}
	}
	
	//change prediction mode to 43, if the epm file is not provided
	if (cPar.a_mode==41 && cPar.file_epm.empty()) {cPar.a_mode=43;}
	
	return;
}



void GEMMA::BatchRun (PARAM &cPar) 
{
	clock_t time_begin, time_start;
	time_begin=clock();	
	
	//Read Files
	cout<<"Reading Files ... "<<endl;
	cPar.ReadFiles();
	if (cPar.error==true) {cout<<"error! fail to read files. "<<endl; return;}

	cPar.CheckData();
	if (cPar.error==true) {cout<<"error! fail to check data. "<<endl; return;}

	if (cPar.a_mode==1 || cPar.a_mode==2 || cPar.a_mode==3 || cPar.a_mode==4 || cPar.a_mode==5 || cPar.a_mode==31) {  //Fit LMM or mvLMM or eigen
		gsl_matrix *Y=gsl_matrix_alloc (cPar.ni_test, cPar.n_ph);
		gsl_matrix *W=gsl_matrix_alloc (Y->size1, cPar.n_cvt);
		gsl_matrix *B=gsl_matrix_alloc (Y->size2, W->size2);	//B is a d by c matrix
		gsl_matrix *se_B=gsl_matrix_alloc (Y->size2, W->size2);
		gsl_matrix *G=gsl_matrix_alloc (Y->size1, Y->size1);
		gsl_matrix *U=gsl_matrix_alloc (Y->size1, Y->size1); 
		gsl_matrix *UtW=gsl_matrix_alloc (Y->size1, W->size2);
		gsl_matrix *UtY=gsl_matrix_alloc (Y->size1, Y->size2);
		gsl_vector *eval=gsl_vector_alloc (Y->size1);
		cPar.CopyCvtPhen (W, Y, 0);
		if (!(cPar.file_kin).empty()) {
			ReadFile_kin (cPar.file_kin, cPar.indicator_idv, cPar.mapID2num, cPar.k_mode, cPar.error, G);
			if (cPar.error==true) {cout<<"error! fail to read kinship/relatedness file. "<<endl; return;}
			CenterMatrix (G);
			cout<<"Start Eigen-Decomposition..."<<endl;
			time_start=clock();	
	
			if (cPar.a_mode==31) {
				cPar.trace_G=EigenDecomp (G, U, eval, 1);
			} else {
				cPar.trace_G=EigenDecomp (G, U, eval, 0);
			}

			cPar.trace_G=0.0;
			for (size_t i=0; i<eval->size; i++) {
				if (gsl_vector_get (eval, i)<1e-10) {gsl_vector_set (eval, i, 0);}
				cPar.trace_G+=gsl_vector_get (eval, i);
			}
			cPar.trace_G/=(double)eval->size;

			cPar.time_eigen=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);	
		} else {
			ReadFile_eigenU (cPar.file_ku, cPar.error, U);
			if (cPar.error==true) {cout<<"error! fail to read the U file. "<<endl; return;}
			
			ReadFile_eigenD (cPar.file_kd, cPar.error, eval);			
			if (cPar.error==true) {cout<<"error! fail to read the D file. "<<endl; return;}
			
			cPar.trace_G=0.0;
			for (size_t i=0; i<eval->size; i++) {
				if (gsl_vector_get(eval, i)<1e-10) {gsl_vector_set(eval, i, 0);}
			  	cPar.trace_G+=gsl_vector_get(eval, i);
			}
			cPar.trace_G/=(double)eval->size;
		}
		
		if (cPar.a_mode==31) {
			cPar.WriteMatrix(U, "eigenU");
			cPar.WriteVector(eval, "eigenD");
		} else {
			CalcUtX (U, W, UtW);
			CalcUtX (U, Y, UtY);			
			if (cPar.n_ph==1) {
				gsl_vector_view beta=gsl_matrix_row (B, 0);
				gsl_vector_view se_beta=gsl_matrix_row (se_B, 0);
				gsl_vector_view UtY_col=gsl_matrix_column (UtY, 0);

				CalcLambda ('L', eval, UtW, &UtY_col.vector, cPar.l_min, cPar.l_max, cPar.n_region, cPar.l_mle_null, cPar.logl_mle_H0);
				CalcLmmVgVeBeta (eval, UtW, &UtY_col.vector, cPar.l_mle_null, cPar.vg_mle_null, cPar.ve_mle_null, &beta.vector, &se_beta.vector);

				cPar.beta_mle_null.clear();
				cPar.se_beta_mle_null.clear();
				for (size_t i=0; i<B->size2; i++) {
					cPar.beta_mle_null.push_back(gsl_matrix_get(B, 0, i) );
					cPar.se_beta_mle_null.push_back(gsl_matrix_get(se_B, 0, i) );
				}

				CalcLambda ('R', eval, UtW, &UtY_col.vector, cPar.l_min, cPar.l_max, cPar.n_region, cPar.l_remle_null, cPar.logl_remle_H0);
				CalcLmmVgVeBeta (eval, UtW, &UtY_col.vector, cPar.l_remle_null, cPar.vg_remle_null, cPar.ve_remle_null, &beta.vector, &se_beta.vector);
				cPar.beta_remle_null.clear();
				cPar.se_beta_remle_null.clear();
				for (size_t i=0; i<B->size2; i++) {
					cPar.beta_remle_null.push_back(gsl_matrix_get(B, 0, i) );
					cPar.se_beta_remle_null.push_back(gsl_matrix_get(se_B, 0, i) );
				}
				
				CalcPve (eval, UtW, &UtY_col.vector, cPar.l_remle_null, cPar.trace_G, cPar.pve_null, cPar.pve_se_null);
				cPar.PrintSummary();
				if (cPar.a_mode==5) {
					gsl_vector *Utu_hat=gsl_vector_alloc (Y->size1);
					gsl_vector *Ute_hat=gsl_vector_alloc (Y->size1);
					gsl_vector *u_hat=gsl_vector_alloc (Y->size1);
					gsl_vector *e_hat=gsl_vector_alloc (Y->size1);
					gsl_vector *y_hat=gsl_vector_alloc (Y->size1);
					gsl_vector_memcpy (y_hat, &UtY_col.vector);
					gsl_blas_dgemv (CblasNoTrans, -1.0, UtW, &beta.vector, 1.0, y_hat);
					double d, u, e;
					for (size_t i=0; i<eval->size; i++) {
						d=gsl_vector_get (eval, i);
						u=cPar.l_remle_null*d/(cPar.l_remle_null*d+1.0)*gsl_vector_get(y_hat, i);
						e=1.0/(cPar.l_remle_null*d+1.0)*gsl_vector_get(y_hat, i);
						gsl_vector_set (Utu_hat, i, u);
						gsl_vector_set (Ute_hat, i, e);
					}
					gsl_blas_dgemv (CblasNoTrans, 1.0, U, Utu_hat, 0.0, u_hat);
					gsl_blas_dgemv (CblasNoTrans, 1.0, U, Ute_hat, 0.0, e_hat);
					cPar.WriteVector(u_hat, "residU");
					cPar.WriteVector(e_hat, "residE");
					gsl_vector_free(u_hat);
					gsl_vector_free(e_hat);
					gsl_vector_free(y_hat);
				}							
			} 
			if (cPar.a_mode==1 || cPar.a_mode==2 || cPar.a_mode==3 || cPar.a_mode==4) {
				if (cPar.n_ph==1) {			
					LMM cLmm;
					cLmm.CopyFromParam(cPar);
					
					gsl_vector_view Y_col=gsl_matrix_column (Y, 0);
					gsl_vector_view UtY_col=gsl_matrix_column (UtY, 0);
					
					if (!cPar.file_gene.empty()) {		
						cLmm.AnalyzeGene (U, eval, UtW, &UtY_col.vector, W, &Y_col.vector); //y is the predictor, not the phenotype
					} else if (!cPar.file_bfile.empty()) {
						cLmm.AnalyzePlink (U, eval, UtW, &UtY_col.vector, W, &Y_col.vector);
					} else {
						cLmm.AnalyzeBimbam (U, eval, UtW, &UtY_col.vector, W, &Y_col.vector);
					}	
					
					cLmm.WriteFiles();
					cLmm.CopyToParam(cPar);
				} else {			 
//					MVLMM cMvlmm;
//					cMvlmm.CopyFromParam(cPar);			
//					
//					if (!cPar.file_bfile.empty()) {
//						cMvlmm.AnalyzePlink (U, eval, UtW, UtY);
//					} else {
//						cMvlmm.AnalyzeBimbam (U, eval, UtW, UtY);
//					}
//					
//					cMvlmm.WriteFiles();
//					cMvlmm.CopyToParam(cPar);
				}
			} 
		}
		gsl_matrix_free (Y);
		gsl_matrix_free (W);
		gsl_matrix_free(B);
		gsl_matrix_free(se_B);
		gsl_matrix_free (G);	
		gsl_matrix_free (U);
		gsl_matrix_free (UtW);
		gsl_matrix_free (UtY);
		gsl_vector_free (eval);
	} 
	cPar.time_total=(clock()-time_begin)/(double(CLOCKS_PER_SEC)*60.0);
	return;
}
void GEMMA::WriteLog (int argc, char ** argv, PARAM &cPar) 
{
	string file_str;
	file_str="./output/"+cPar.file_out;
	file_str+=".log.txt";
	
	ofstream outfile (file_str.c_str(), ofstream::out);
	if (!outfile) {cout<<"error writing log file: "<<file_str.c_str()<<endl; return;}
	
	outfile<<"##"<<endl;
	outfile<<"## GEMMA Version = "<<version<<endl;
	
	outfile<<"##"<<endl;
	outfile<<"## Command Line Input = ";
	for(int i = 1; i < argc; i++) {	
		outfile<<argv[i]<<" ";
	}
	outfile<<endl;
	
	outfile<<"##"<<endl;
	outfile<<"## Summary Statistics:"<<endl;
	outfile<<"## number of total individuals = "<<cPar.ni_total<<endl;	
	if (cPar.a_mode==43) {
		outfile<<"## number of analyzed individuals = "<<cPar.ni_cvt<<endl;
		outfile<<"## number of individuals with full phenotypes = "<<cPar.ni_test<<endl;
	} else {
		outfile<<"## number of analyzed individuals = "<<cPar.ni_test<<endl;
	}
	outfile<<"## number of covariates = "<<cPar.n_cvt<<endl;
	outfile<<"## number of phenotypes = "<<cPar.n_ph<<endl;
	if (cPar.a_mode==43) {
		outfile<<"## number of observed data = "<<cPar.np_obs<<endl;
		outfile<<"## number of missing data = "<<cPar.np_miss<<endl;
	}
		
	if (!(cPar.file_gene).empty()) {
		outfile<<"## number of total genes = "<<cPar.ng_total<<endl;
		outfile<<"## number of analyzed genes = "<<cPar.ng_test<<endl;		
	} else if (cPar.file_epm.empty()) {	
		outfile<<"## number of total SNPs = "<<cPar.ns_total<<endl;	
		outfile<<"## number of analyzed SNPs = "<<cPar.ns_test<<endl;
	} else {
		outfile<<"## number of analyzed SNPs = "<<cPar.ns_test<<endl;
	}
	
	if (cPar.a_mode==13) {
		outfile<<"## number of cases = "<<cPar.ni_case<<endl;
		outfile<<"## number of controls = "<<cPar.ni_control<<endl;
	}
	
	if (cPar.a_mode==1 || cPar.a_mode==2 || cPar.a_mode==3 || cPar.a_mode==4 || cPar.a_mode==5 || cPar.a_mode==11 || cPar.a_mode==12 || cPar.a_mode==13) {
		outfile<<"## REMLE log-likelihood in the null model = "<<cPar.logl_remle_H0<<endl;
		outfile<<"## MLE log-likelihood in the null model = "<<cPar.logl_mle_H0<<endl;
		if (cPar.n_ph==1) {
			//outfile<<"## lambda REMLE estimate in the null (linear mixed) model = "<<cPar.l_remle_null<<endl;
			//outfile<<"## lambda MLE estimate in the null (linear mixed) model = "<<cPar.l_mle_null<<endl;	
			outfile<<"## pve estimate in the null model = "<<cPar.pve_null<<endl;
			outfile<<"## se(pve) in the null model = "<<cPar.pve_se_null<<endl;	
			outfile<<"## vg estimate in the null model = "<<cPar.vg_remle_null<<endl;
			outfile<<"## ve estimate in the null model = "<<cPar.ve_remle_null<<endl;	
			outfile<<"## beta estimate in the null model = ";
			for (size_t i=0; i<cPar.beta_remle_null.size(); i++) {
				outfile<<"  "<<cPar.beta_remle_null[i];
			}
			outfile<<endl;
			outfile<<"## se(beta) = ";
			for (size_t i=0; i<cPar.se_beta_remle_null.size(); i++) {
				outfile<<"  "<<cPar.se_beta_remle_null[i];
			}
			outfile<<endl;
			
		} else {
			size_t c;
			outfile<<"## REMLE estimate for Vg in the null model: "<<endl;			
			for (size_t i=0; i<cPar.n_ph; i++) {
				for (size_t j=0; j<=i; j++) {
					c=(2*cPar.n_ph-min(i,j)+1)*min(i,j)/2+max(i,j)-min(i,j);
					outfile<<cPar.Vg_remle_null[c]<<"\t";
				}
				outfile<<endl;
			}
			outfile<<"## se(Vg): "<<endl;	
			for (size_t i=0; i<cPar.n_ph; i++) {
				for (size_t j=0; j<=i; j++) {
					c=(2*cPar.n_ph-min(i,j)+1)*min(i,j)/2+max(i,j)-min(i,j);
					outfile<<sqrt(cPar.VVg_remle_null[c])<<"\t";
				}
				outfile<<endl;
			}
			outfile<<"## REMLE estimate for Ve in the null model: "<<endl;	
			for (size_t i=0; i<cPar.n_ph; i++) {
				for (size_t j=0; j<=i; j++) {
					c=(2*cPar.n_ph-min(i,j)+1)*min(i,j)/2+max(i,j)-min(i,j);
					outfile<<cPar.Ve_remle_null[c]<<"\t";
				}
				outfile<<endl;
			}
			outfile<<"## se(Ve): "<<endl;	
			for (size_t i=0; i<cPar.n_ph; i++) {
				for (size_t j=0; j<=i; j++) {
					c=(2*cPar.n_ph-min(i,j)+1)*min(i,j)/2+max(i,j)-min(i,j);
					outfile<<sqrt(cPar.VVe_remle_null[c])<<"\t";
				}
				outfile<<endl;
			}
			
			outfile<<"## MLE estimate for Vg in the null model: "<<endl;
			for (size_t i=0; i<cPar.n_ph; i++) {
				for (size_t j=0; j<cPar.n_ph; j++) {
					c=(2*cPar.n_ph-min(i,j)+1)*min(i,j)/2+max(i,j)-min(i,j);
					outfile<<cPar.Vg_mle_null[c]<<"\t";
				}
				outfile<<endl;
			}
			outfile<<"## se(Vg): "<<endl;	
			for (size_t i=0; i<cPar.n_ph; i++) {
				for (size_t j=0; j<=i; j++) {
					c=(2*cPar.n_ph-min(i,j)+1)*min(i,j)/2+max(i,j)-min(i,j);
					outfile<<sqrt(cPar.VVg_mle_null[c])<<"\t";
				}
				outfile<<endl;
			}
			outfile<<"## MLE estimate for Ve in the null model: "<<endl;	
			for (size_t i=0; i<cPar.n_ph; i++) {
				for (size_t j=0; j<cPar.n_ph; j++) {
					c=(2*cPar.n_ph-min(i,j)+1)*min(i,j)/2+max(i,j)-min(i,j);
					outfile<<cPar.Ve_mle_null[c]<<"\t";
				}
				outfile<<endl;
			}
			outfile<<"## se(Ve): "<<endl;	
			for (size_t i=0; i<cPar.n_ph; i++) {
				for (size_t j=0; j<=i; j++) {
					c=(2*cPar.n_ph-min(i,j)+1)*min(i,j)/2+max(i,j)-min(i,j);
					outfile<<sqrt(cPar.VVe_mle_null[c])<<"\t";
				}
				outfile<<endl;
			}
			outfile<<"## estimate for B (d by c) in the null model (columns correspond to the covariates provided in the file): "<<endl;
			for (size_t i=0; i<cPar.n_ph; i++) {
				for (size_t j=0; j<cPar.n_cvt; j++) {
					c=i*cPar.n_cvt+j;
					outfile<<cPar.beta_remle_null[c]<<"\t";
				}
				outfile<<endl;
			}
			outfile<<"## se(B): "<<endl;
			for (size_t i=0; i<cPar.n_ph; i++) {
				for (size_t j=0; j<cPar.n_cvt; j++) {
					c=i*cPar.n_cvt+j;
					outfile<<cPar.se_beta_remle_null[c]<<"\t";
				}
				outfile<<endl;
			}
		}
	}
	
	
	if (cPar.a_mode==11 || cPar.a_mode==12 || cPar.a_mode==13) {
		outfile<<"## estimated mean = "<<cPar.pheno_mean<<endl;
	}
	
	if (cPar.a_mode==11 || cPar.a_mode==13) {	
		outfile<<"##"<<endl;
		outfile<<"## MCMC related:"<<endl;	
		outfile<<"## initial value of h = "<<cPar.cHyp_initial.h<<endl;
		outfile<<"## initial value of rho = "<<cPar.cHyp_initial.rho<<endl;
		outfile<<"## initial value of pi = "<<exp(cPar.cHyp_initial.logp)<<endl;
		outfile<<"## initial value of |gamma| = "<<cPar.cHyp_initial.n_gamma<<endl;
		outfile<<"## random seed = "<<cPar.randseed<<endl;
		outfile<<"## acceptance ratio = "<<(double)cPar.n_accept/(double)((cPar.w_step+cPar.s_step)*cPar.n_mh)<<endl;
	}
	
	outfile<<"##"<<endl;
	outfile<<"## Computation Time:"<<endl;
	outfile<<"## total computation time = "<<cPar.time_total<<" min "<<endl;
	outfile<<"## computation time break down: "<<endl;
	if (cPar.a_mode==1 || cPar.a_mode==2 || cPar.a_mode==3 || cPar.a_mode==4 || cPar.a_mode==5 || cPar.a_mode==11 || cPar.a_mode==12 || cPar.a_mode==13) {
		outfile<<"##      time on eigen-decomposition = "<<cPar.time_eigen<<" min "<<endl;
		outfile<<"##      time on calculating UtX = "<<cPar.time_UtX<<" min "<<endl;		
	}
	if (cPar.a_mode==1 || cPar.a_mode==2 || cPar.a_mode==3 || cPar.a_mode==4) {
		outfile<<"##      time on optimization = "<<cPar.time_opt<<" min "<<endl;
	}
	if (cPar.a_mode==11 || cPar.a_mode==13) {
		outfile<<"##      time on mcmc = "<<cPar.time_opt<<" min "<<endl;
		outfile<<"##      time on Omega = "<<cPar.time_Omega<<" min "<<endl;
	}
	if (cPar.a_mode==41 || cPar.a_mode==42) {
		outfile<<"##      time on eigen-decomposition = "<<cPar.time_eigen<<" min "<<endl;
	}
	if (cPar.a_mode==43) {
		outfile<<"##      time on eigen-decomposition = "<<cPar.time_eigen<<" min "<<endl;
		outfile<<"##      time on predicting phenotypes = "<<cPar.time_opt<<" min "<<endl;
	}
	outfile<<"##"<<endl;
	
	outfile.close();
	outfile.clear();
	return;
}


