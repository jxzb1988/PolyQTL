#ifndef _CONDITIONAL_ANALYSIS_H
#define _CONDITIONAL_ANALYSIS_H

#include<iostream>
#include<fstream>
#include <armadillo>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sstream>
#include <string>
#include <vector>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <regex>
#include <iomanip>


using namespace std;
using namespace arma;

void read_table(string data, vector< vector<double> > &test,  vector<string> &ind, vector<string> &variant);
int match(vector<string> &data, string target);
double calculate_VIF(mat &input,vector<string> &variant,string target,vector<string> &detected);
void readgemma(string input, vector < vector<double> > &out_gemma, vector<string> & variant);
int  updatefam(vector< vector<string> > &input , mat &value,string output);
void copybim(string input, string output);
void copybim(string input, string output, vector<int> &index);
void copyfam(string input, string output);
void  outputbed(mat & _snp_2, mat & _snp_1, string &output);
void readbim(string file,map<string, int> &order);
int readfam(string file);
int readfam(string file,vector<double> &pheno, vector< vector<string> > &Fam_infor);
int readbim(string file);
void readbim(string file,vector<string>&variant);
void extractvariant(string input, string output, double r2_cutoff, string target);
void copybed(string input, string output);
int conditional_analysis(string input,string gene, int totalCausalSNP, float rho, bool histFlag, double gamma, string weight, int nthread, string covariate, vector<string> &grm_file, string outputFile, string geno_file);
int conditional_analysis(string input,string gene, int totalCausalSNP, float rho, bool histFlag);

#endif
