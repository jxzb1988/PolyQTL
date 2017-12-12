
//The part of the function in the code from adapted from CAVIAR developped by Hormozdiari, F et al.



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream> 
#include <sstream>
#include <string>
#include <vector>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

long int fact(int n) {
        if(n==0)
                return 1;
        return n* fact(n-1);
}

void copyConfigure(double *dest, double *src, int size) {
	for(int i = 0; i < size; i++) 
		dest[i] = src[i];
}

double min(double a, double b) {
	if(a>b)
		return b;
	else
		return a;
}

long int nCr(int n, int r) {
        long int result = 1;
        for(int i = n; i > n-r; i--)
                result *= i;
        return result/fact(r);
}

void printVector(char * data, int size) {
        for(int i = 0; i < size; i++)
                printf("%c, ", data[i]);
}

void printVector(int * data, int size) {
        for(int i = 0; i < size; i++)
                printf("%d, ", (int)data[i]);
}

void printVector(double * data, int size) {
        for(int i = 0; i < size; i++)
                printf("%lf, ", data[i]);
}

void diffVector(double * data1, double * data2, int size, double * result) {
	for(int i = 0; i < size; i++ ) 
		result[i] = data1[i] - data2[i];
}

void sumVector(double * data1, double * data2, int size, double * result) {
        for(int i = 0; i < size; i++ ) 
                result[i] = data1[i] + data2[i];
}

double multVector(double * data1, double * data2, int size) {
	double res = 0;
	for(int i = 0; i < size; i++ ) 
                res += data1[i] * data2[i];
	return res;
}

void dotVector(double * data1, double * data2, int size, double * result) {
	for(int i = 0; i < size; i++ ) 
                result[i] = data1[i] * data2[i];
}

void multVectorMatrix(double *vector, double * matrix, int size, double * result) {
	double total_row = 0;
	for(int i = 0; i < size; i++) {
		total_row = 0;
		for(int j = 0; j < size; j++) {
			total_row += vector[j] * matrix[i + j * size];
		}
		result[i]= total_row;
	}
}

void importData(string fileName, double * vector) {
	int index = 0;
	double data = 0;
        ifstream fin(fileName.c_str(), std::ifstream::in);
        fin >> data;
        while (fin.good()) {
                vector[index] = data;
                index++;
                fin >> data;
        }
        fin.close();
}

void importData(string fileName, int * vector) {
        int index = 0;
        double data = 0;
        ifstream fin(fileName.c_str(), std::ifstream::in);
        while( fin.good()  ){
                fin >> data;
                vector[index] = (int)data;
                index++;
        }
        fin.close();
}

/*
	The column index starts by 1 in this implemenation
*/
void importDataSecondColumn(string fileName, double * vector) {
	int index = 0;
	string line = "";
	string dataS = "";
	double data = 0.0;	
	ifstream fin(fileName.c_str(), std::ifstream::in);
	while( getline(fin, line) ){
		istringstream iss(line);
		iss >> dataS;
		iss >> data;
	        vector[index] = (double)data;
                index++;
        }
//	cout << "reach=" << index << endl;
        fin.close();
}

/*
	The column index starts by 1 in this implemenation
*/
void importDataNthColumn(string fileName, double * vector, int colNum) {
        int index = 0;
        string line = "";
        string dataS = "";
        double data = 0.0;
	ifstream fin(fileName.c_str(), std::ifstream::in);
        while( getline(fin, line) ){
                istringstream iss(line);
                iss >> dataS;
                for(int i = 0; i < colNum-1;i++)
			iss >> data;
                vector[index] = (double)data;
                index++;
        }
 //       cout << "reach=" << index << endl;
        fin.close();
}

void importDataFirstColumn(string fileName, double * list) {
   //     cout<<"Come into to extract variants names"<<endl;
 	int index = 0;
        double data = 0.0;
        string line = "";
	ifstream fin(fileName.c_str(), std::ifstream::in);
   //     cout<<"I am safe in extracting variants"<<endl;
        while( getline(fin, line) ){
   //             cout<<"Come into"<<endl;
   //             cout<<"line is: "<<line<<endl;
//		istringstream iss(line);
                istringstream iss;
   //             string test_b="Biao Zeng is the best";
          //      istringstream iss_test;
   //             cout<<"Tested string is: "<<test_b<<endl;
        //        iss_test.str(test_b);
   //             cout<<"It was changed"<<endl;
                iss.str(line);
   //             cout<<"Still safe 1"<<endl;
//                cout<<"Come into"<<endl;
                iss >> data;
   //             cout<<"Still safe 2"<<endl;
                list[index] = (double) data ;
   //             cout<<"Still safe 3"<<endl;
		index++;
          }
        double mean=0;
        for(int i=0;i<index;i++)
         {
           mean+=list[i];
         }
        mean/=index;
        for(int i=0;i<index;i++)
         {
           list[i]-=mean;
         }
     //   double mean=0;
     //   for(int i=0;i<index;i++)
     //    {
     //      mean+=list[index]; 
     //    }
     //   mean=mean/index;
     //   for(int i=0;i<index;i++)
     //    {
     //      list[index]-=mean;
     //    }
	cout << "FINISH" << endl;
        fin.close();
}

void fileSize(string fileName, int & size) {
    //    cout<<"Come into to extract variant number"<<endl;
	size = 0;
	double data = 0;	
	ifstream fin(fileName.c_str(), std::ifstream::in);
	if( fin.good()  ){

		fin >> data;
      //          cout<<" "<<data;
		size++;
	}
     //   cout<<endl;
     //   cout<<"size is: "<<size<<endl;
	fin.close();
}

string convertInt(int number) {
   stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}

void resetVector(char *data, int size){
	for(int i = 0; i < size; i++)
		data[i] = '0';
}

void resetVector(int * data, int size) {
	for(int i = 0; i < size; i++)
		data[i] = 0;
}

void resetVector(double * data, int size) {
        for(int i = 0; i < size; i++)
                data[i] = 0;
}

void exportVector2File(string fileName, char * data, int size) {
	ofstream outfile(fileName.c_str(), ios::out );
	for (int i = 0; i < size; i++)
		outfile << data[i] << " ";
	outfile << endl;
	outfile.close();
}

void exportVector2File(string fileName, double * data, int size) {
        ofstream outfile(fileName.c_str(), ios::out );
        for (int i = 0; i < size; i++)
                outfile << data[i] << " ";
        outfile << endl;
        outfile.close();
}

void exportVector2File(string fileName, int * data, int size) {
        ofstream outfile(fileName.c_str(), ios::out );
        for (int i = 0; i < size; i++)
                outfile << data[i] << " ";
        outfile << endl;
        outfile.close();
}

void export2File(string fileName, int data) {
        ofstream outfile(fileName.c_str(), ios::out );
        outfile << data << endl;
        outfile.close();
}

int snp2Gene(int * G, int snpId, int snpCount, int geneCount) {
	for(int i = 0; i < geneCount; i++) {
		if(G[snpId*geneCount + i] == 1)
			return i;
	}
	return -1;
}

void setIdentitymatrix(int * G, int snpCount, int geneCount) {
	for(int i = 0; i < snpCount; i++) {
		for(int j = 0; j < geneCount; j++) {
			G[i*geneCount + j] = 0;
		}
		G[i*geneCount + (i/(snpCount/geneCount))] = 1;
	}
	  for(int i = 0; i < snpCount; i++) {
                for(int j = 0; j < geneCount; j++) {
                        printf("%d ", G[i*geneCount+j]);
                }
		printf("\n");
        }
}

void makeSigmaPositiveSemiDefinite(double * sigma, int size) {
	int gsl_tmp = 0;
	double matDet  = 0;
        double addDiag = 0;
        bool positive = false;
 //       cout<<"Come into the code of makeSigmaPositiveSemiDefinite"<<endl;
 //       cout<<"Value of size is: "<<size<<endl;	
	//gsl_set_error_handler_off();	
	gsl_matrix * tmpResultMatrix = gsl_matrix_calloc (size, size);	
 //       cout<<"Come into the code of makeSigmaPositiveSemiDefinite"<<endl;
	gsl_permutation *p = gsl_permutation_alloc(size);
 //       cout<<"Come into the code of makeSigmaPositiveSemiDefinite"<<endl;
        do{
		for(int i = 0; i < size; i++) {
                	for (int j = 0; j < size; j++) {
                      		if(i==j)
					gsl_matrix_set(tmpResultMatrix,i,j,sigma[i*size+j]+addDiag);
				else
					gsl_matrix_set(tmpResultMatrix,i,j,sigma[i*size+j]);
			}
        	}
	
		gsl_linalg_LU_decomp(tmpResultMatrix, p, &gsl_tmp);
       		matDet = gsl_linalg_LU_det(tmpResultMatrix,gsl_tmp);	
//		cout << matDet << "\t" << addDiag << endl;
		if(matDet > 0 ) 
			positive = true;
		else {
//			cout << "add" << endl;
			addDiag+=0.1;		
		}
	} while(!positive);
 //       cout<<"Come into the code of makeSigmaPositiveSemiDefinite"<<endl;
	for(int i = 0; i < size*size; i++){
                if(i%(size+1) == 0)
                        sigma[i] = sigma[i] + addDiag;
        }
 //       cout<<"Come into the code of makeSigmaPositiveSemiDefinite"<<endl;
}
