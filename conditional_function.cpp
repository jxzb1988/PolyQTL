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
#include "gemma_param.h"
#include "gemma.h"
#include <iomanip>
#ifndef _MEQTLPOLYGMODEL_H
#define _MEQTLPOLYGMODEL_H
#include "MeQTLPolyGModel.h"
#endif

#include "conditional_function.h"


using namespace std;
using namespace arma;


//Read phe+genotype table
void read_table(string data, vector< vector<double> > &test,  vector<string> &ind, vector<string> &variant)
  {
//    cout<<"Come into input file to read information"<<endl;
    ifstream FHIN1(data.c_str(),ifstream::in);
    string a;
    int b=0;
    while(FHIN1 &&getline(FHIN1,a))
     {
       b++;
       if(b==1)
        {

          stringstream ss(a);
          string word;
          ss>>word;
          while(ss && ss>>word)
           {
             variant.push_back(word);
           }
         } else
         {

           double c=0;
           stringstream ss(a);
           string waste;
           ss>>waste;
           ind.push_back(waste);
           vector<double> T;
           while(ss && ss>>c)
            {

              T.push_back(c);
            }
           test.push_back(T);
         }
      }
//    cout<<"Read phenotype and genotype finished"<<endl;

    FHIN1.clear();
    FHIN1.close();
  }



//Find the index of a specific string in a string vector
int match(vector<string> &data, string target)
 {
   cout<<"Come into function match"<<endl;
//   cout<<"Target data is: "<<endl;
   
   int j=-1;
   vector<string>::iterator Col;

   for(Col=data.begin();Col!=data.end();Col++)
    {
      j++;
      if((*Col)==target)
       {
         cout<<"Good, the string is detected"<<endl;
         cout<<"Index for the target string is: "<<j<<endl;
     
         return j;
       }
    }
   return j;
 }

//Calculate the variance inflation factor to determine whether include the explored variant or not
double calculate_VIF(mat &input,vector<string> &variant,string target,vector<string> &detected)
 {
   cout<<"Come to function calculate_VIF"<<endl;
   int res=match(variant,target);
   mat X_chosen=mat(input.n_rows,detected.size(),fill::zeros);
   for(int i=0;i<detected.size();i++)
    {
      int ind=match(variant, detected[i]);
   
      X_chosen.col(i)=input.col(ind);
    }
   cout<<"Genotype get, and start to calculate VIF"<<endl;

   mat stat_test=input.col(res);
   mat XtX=trans(X_chosen)*X_chosen;
   mat inv_XtX=pinv(XtX);

   mat res_test=inv_XtX*trans(X_chosen)*stat_test;
   mat resi=stat_test-X_chosen*res_test;
   double  phe_variance=stddev(resi.col(0));
   colvec residuals = stat_test - X_chosen * res_test;
   double s2=0;
   double s=0;
   for(int i_x=0;i_x<input.n_rows;i_x++)
    {
      s2+=residuals(i_x)*residuals(i_x);
      s+= stat_test(i_x,0)*stat_test(i_x,0);
    }
   double r2=1-s2/s;
   double vif=1/r2;
   return vif;
 }




//read GEMMA result output
void readgemma(string input, vector < vector<double> > &out_gemma, vector<string> & variant)
 {
   ifstream FHIN(input.c_str(),ios::in);
   if(!FHIN)
    {
//      cout<<"There is some error for the input, please have a check"<<endl;
      return ;
    }
   string A;
   getline(FHIN,A);
   while(FHIN && getline(FHIN, A))
    {
      stringstream ss(A);
      vector<string> X;
      string B;
      vector<double> Y;
      ss>>B;
      ss>>B;
      variant.push_back(B);
      for(int i=3;i<8;i++)
       {
         ss>>B;
       }
      double test;
      ss>>test;
      Y.push_back(test);
      ss>>test;
      ss>>test;
      ss>>test;
      Y.push_back(test);
      out_gemma.push_back(Y);
    }
   FHIN.clear();
   FHIN.close();

 }


//Attached new residual phenotype, removing effect of detected SNPs, to plink fam file.
int  updatefam(vector< vector<string> > &input , mat &value,string output)
 {
//   cout<<"Come into updatefam fucntion"<<endl;
   string A;
   fstream  FHOU(output.c_str(),ios::out);
   if(!FHOU)
    {
      cout<<"Error for input or output"<<endl;
      return (1);
    }
   if(value.n_rows!=input.size())
    {
//      cout<<"The number of individuals in phenotype doesn't match that in fam file "<<endl;
      return (1);
    }
   for(int i=0;i<input.size();i++)
    {
      FHOU<<input[i][0];
      for(int j=1;j<input[i].size();j++)
       {
         FHOU<<"        "<<input[i][j];
       }
      FHOU<<"   "<<value(i,0)<<endl;
    }
   FHOU.clear();
   FHOU.clear();
   return (0);
 }



//Make a copy of plink bim file
void copybim(string input, string output)
 {
   ifstream FHIN(input.c_str(),ios::in);
   string A;
   fstream  FHOU(output.c_str(),ios::out);
   if(!FHIN || !FHOU)
    {
      cout<<"Error for input or output"<<endl;
      return;
    }
   while(FHIN && getline(FHIN,A))
    {
      FHOU<<A;
      FHOU<<endl;
    }
   FHIN.clear();
   FHOU.clear();
   FHIN.close();
   FHOU.clear();
 }

//Extract a specific set of variants, and output a new bim file
void copybim(string input, string output, vector<int> &index)
 {
   ifstream FHIN(input.c_str(),ios::in);
   string A;
   fstream  FHOU(output.c_str(),ios::out);
   if(!FHIN || !FHOU)
    {
      cout<<"Error for input or output"<<endl;
      return;
    }
   int val=-1;
   while(FHIN && getline(FHIN,A))
    {
      val++;
      if(find(index.begin(), index.end(), val) != index.end())
       {       
         FHOU<<A;
         FHOU<<endl;
       }
    }
   FHIN.clear();
   FHOU.clear();
   FHIN.close();
   FHOU.clear();
 }


//Make a copy of fam file
void copyfam(string input, string output)
 {
   ifstream FHIN(input.c_str(),ios::in);
   string A;
   fstream FHOU(output.c_str(),ios::out);
   if(!FHIN || !FHOU)
    {
      cout<<"Error for input or output"<<endl;
    }
   while(FHIN && getline(FHIN,A))
    {
      stringstream ss(A);
      vector<string> B;
      string Test;
      while(ss && ss>>Test)
       {
         FHOU<<Test<<"	"; 
       }
      FHOU<<endl;
    }
   FHIN.clear();
   FHIN.close();
   FHOU.clear();
   FHOU.close();
 }

//Make a copy of bed file
void  outputbed(mat & _snp_2, mat & _snp_1, string &output)
 {
   int i=0, pos=0, j=0;
   string OutBedFile=output;
   fstream OutBed(OutBedFile.c_str(), ios::out|ios::binary);
   if(!OutBed) throw("Error: can not open the file ["+OutBedFile+"] to write.");
//   cout<<"Writing genotypes to PLINK BED file ["+OutBedFile+"] ..."<<endl;
   bitset<8> b;
   char ch[1];
   b.reset();
   b.set(2);  b.set(3);  b.set(5);  b.set(6);
   ch[0] = (char)b.to_ulong();
   OutBed.write(ch,1);
   b.reset();
   b.set(0);  b.set(1);  b.set(3);  b.set(4);
   ch[0] = (char)b.to_ulong();
   OutBed.write(ch,1);
   b.reset();
   b.set(0);
   ch[0] = (char)b.to_ulong();
   OutBed.write(ch,1);
   int n_sid=_snp_2.n_rows;
   int n_ind=_snp_2.n_cols;
//   cout<<"Number of individuals is: "<<n_ind<<endl;
//   cout<<"Number of variants is: "<<n_sid;
//   cout<<"Right now, it is OK"<<endl;
   for(i=0; i<n_sid; i++){
      pos=0;
      b.reset();
      for(j=0; j<n_ind; j++){
        b[pos++]=(!_snp_2(i,j));
        b[pos++]=(!_snp_1(i,j));
        if(pos>7 || j==n_ind-1){
          ch[0]=(char)b.to_ulong();
          OutBed.write(ch,1);
          pos=0;
          b.reset();
        }
      }
    }
   OutBed.close();
 }

//Extract a specific set of variants from a bim file 
void readbim(string file,map<string, int> &order)
 {

   int ibuf=-1;
   string A;
   string cbuf="0";
   double dbuf=0.0;
   string str_buf;
   ifstream Bim(file.c_str());
   if(!Bim) throw("Error: can not open the file ["+file+"] to read.");
//   cout<<"Reading PLINK BIM file from ["+file+"]."<<endl;
   while(Bim && getline(Bim, A))
    {
      ibuf++;
      stringstream ss(A);
      string waste;
      ss>>waste;
      ss>>waste;
//      cout<<"waste is: "<<waste<<", and ibuf is: "<<ibuf<<endl;
      order.insert(pair<string,int>(waste,ibuf));
    }
//   cout<<"Read bim file is over"<<endl;
   Bim.close();
 }

//Read a bim file
void readbim(string file,vector<string>&variant)
 {

   int ibuf=-1;
   string A;
   string cbuf="0";
   double dbuf=0.0;
   string str_buf;
   ifstream Bim(file.c_str());
   if(!Bim) throw("Error: can not open the file ["+file+"] to read.");

   while(Bim && getline(Bim, A))
    {
      ibuf++;
      stringstream ss(A);
      string waste;
      ss>>waste;
      ss>>waste;
      variant.push_back(waste);
    }
   Bim.close();
 }


//Read a fam file, and calcuate the individual number 
int readfam(string file)
 {
   cout<<"Come to estimate individual number"<<endl;
   ifstream Fam(file.c_str());
   string A;
   if(!Fam) throw("Error: can not open the file ["+file+"] to read.");
   cout<<"Reading PLINK FAM file from ["+file+"]."<<endl;
   int i=0;
   while(Fam && getline(Fam, A))
    {
      i++;
    }
   Fam.clear();
   Fam.close();
   cout<<"Finish estimation"<<endl;
   return i;
 }

//Read a fam file and extract phenotype 
int readfam(string file,vector<double> &pheno, vector< vector<string> > &Fam_infor)
 {
   ifstream Fam(file.c_str());
   string A;
   if(!Fam) throw("Error: can not open the file ["+file+"] to read.");
   cout<<"Reading PLINK FAM file from ["+file+"]."<<endl;
   while(Fam && getline(Fam, A))
    {
      stringstream ss(A);
      vector<string> T;
      string TT;
      for(int j=0;j<5;j++)
       {
         ss>>TT;
         T.push_back(TT);
       }
      Fam_infor.push_back(T);
      double TTT;
      ss>>TTT;
      pheno.push_back(TTT);
    }
   Fam.clear();
   Fam.close();
 }

//Read a bim file, and calculate the number of variants
int readbim(string file)
 {
   int ibuf=0;
   string A;
   string cbuf="0";
   double dbuf=0.0;
   string str_buf;
   ifstream Bim(file.c_str());
   if(!Bim) throw("Error: can not open the file ["+file+"] to read.");
   cout<<"Reading PLINK BIM file from ["+file+"]."<<endl;
   while(Bim && getline(Bim, A))
    {
      ibuf++;
    }
   Bim.close();
   return ibuf;
 }


//Given a r2 cutoff, extract the variants locating in high LD with a specific varaint
void extractvariant(string input, string output, double r2_cutoff, string target)
 {
    int i=0, j=0, k=0;
    string strfam="fam";
    string strbim="bim";
    string famfile=input;
    string bimfile=input;
    famfile.replace(famfile.end()-3,famfile.end(),strfam);
    bimfile.replace(bimfile.end()-3,bimfile.end(),strbim);
    string famout=output;
    string bimout=output;
    famout.replace(famout.end()-3,famout.end(),strfam);
    bimout.replace(bimout.end()-3,bimout.end(),strbim);
    int nsnp=readbim(bimfile);
    int nind=readfam(famfile);
    mat _snp_2 =mat(nsnp,nind, fill::zeros);
    mat _snp_1 =mat(nsnp,nind, fill::zeros);
    char ch[1];
    bitset<8> b;
    fstream BIT(input.c_str(), ios::in|ios::binary);
    if(!BIT) throw("Error: can not open the file ["+input+"] to read.");
    for(i=0; i<3; i++) BIT.read(ch,1); // skip the first three bytes
    int snp_indx=0, indi_indx=0;
    for(j=0, snp_indx=0; j<nsnp; j++)
     {
       for(i=0, indi_indx=0; i<nind;)
        {
          BIT.read(ch,1);
          if(!BIT) throw("Error: problem with the BED file ... has the FAM/BIM file been changed?");
          b=ch[0];
          k=0;
          while(k < 7 && i < nind)
           {
             _snp_2(snp_indx,indi_indx)=(!b[k++]);
             _snp_1(snp_indx,indi_indx)=(!b[k++]);
             indi_indx++;
             i++;
           }
        }
       snp_indx++;
     }
    cout<<"Reading PLINK BED file is over"<<endl;
    mat dosage=mat(nsnp,nind, fill::zeros);
    for(int i=0;i<_snp_2.n_rows;i++)
     {
       for(int j=0;j<_snp_2.n_cols;j++)
        {
          if(_snp_2(i,j)==0 && _snp_1(i,j)==1)
           {
             dosage(i,j)=-9;
           } else
           {
             dosage(i,j)=_snp_2(i,j)+_snp_1(i,j);
           }
        }
     }
    for(int i=0;i<dosage.n_rows;i++)
     {
       double Mean=0;
       int num=0;
       for(int j=0;j<dosage.n_cols;j++)
        {
          if(dosage(i,j)!=-9)
           {
             Mean+=dosage(i,j);
             num++;
           }
        }
       Mean/=num;
       for(int j=0;j<dosage.n_cols;j++)
        {
          if(dosage(i,j)==-9)
           {
             dosage(i,j)=Mean;
           }
        }
     }

    map<string,int> order;
    readbim(bimfile,order);
    map<string,int>::iterator Ite; 
    Ite=order.find(target);
    if(Ite==order.end())
     {
       return ;
     }
    int T_index=Ite->second;
//    cout<<"T_index is: "<<T_index<<endl;
    vector<int> index;
    vector< vector<int> > _Snp_2;
    vector< vector<int> > _Snp_1;
    vector <int> C;
    for(int x=0;x<_snp_2.n_rows;x++)
     {
       mat r2=cor(dosage.row(x),dosage.row(T_index));
 //      cout<<"r2 is: "<<r2<<endl;
       if(r2(0,0)*r2(0,0)>=r2_cutoff)
        {
          C.push_back(x);
          vector<int> Test1;
          vector<int> Test2;
          for(int y=0;y<_snp_2.n_cols;y++)
           {
             Test1.push_back(_snp_1(x,y));
             Test2.push_back(_snp_2(x,y));
           }
          _Snp_2.push_back(Test2);
          _Snp_1.push_back(Test1);
        }
     }

    mat _snp_2_target=mat(_Snp_2.size(),_Snp_2[0].size(),fill::zeros);
    mat _snp_1_target=mat(_Snp_1.size(),_Snp_1[0].size(),fill::zeros);
    for(int i=0;i<_Snp_1.size();i++)
     {
       for(int j=0;j<_Snp_1[0].size();j++)
        {
          _snp_1_target(i,j)=_Snp_1[i][j];
          _snp_2_target(i,j)=_Snp_2[i][j];
        }
     }
    cout<<"Variants locating in high LD with target variant is over"<<endl;
    outputbed(_snp_2_target,_snp_1_target,output);
    copyfam(famfile,famout);
    copybim(bimfile,bimout,C);
    BIT.clear();
    BIT.close();
 }

//Make a copy of plink bed file
void copybed(string input, string output)
 {
    cout<<"Come into copybed function"<<endl;
    int i=0, j=0, k=0;
    string strfam="fam";
    string strbim="bim";
    string famfile=input;
    string bimfile=input;

    famfile.replace(famfile.end()-3,famfile.end(),strfam);
    bimfile.replace(bimfile.end()-3,bimfile.end(),strbim);
    int nsnp=readbim(bimfile);
    int nind=readfam(famfile);
    cout<<"Number of variants is: "<<nsnp<<endl;
    cout<<"Number of individuals is: "<<nind<<endl;
    mat _snp_2 =mat(nsnp,nind, fill::zeros);
    mat _snp_1 =mat(nsnp,nind, fill::zeros);
    char ch[1];
    bitset<8> b;
    fstream BIT(input.c_str(), ios::in|ios::binary);
    if(!BIT) throw("Error: can not open the file ["+input+"] to read.");
//    cout<<"Reading PLINK BED file from ["+input+"] in SNP-major format ..."<<endl;
    for(i=0; i<3; i++) BIT.read(ch,1); // skip the first three bytes
    int snp_indx=0, indi_indx=0;
    for(j=0, snp_indx=0; j<nsnp; j++)
     {
       for(i=0, indi_indx=0; i<nind;)
        {
          BIT.read(ch,1);
          if(!BIT) throw("Error: problem with the BED file ... has the FAM/BIM file been changed?");
          b=ch[0];
          k=0;
          while(k < 7 && i < nind)
           {
             _snp_2(snp_indx,indi_indx)=(!b[k++]);
             _snp_1(snp_indx,indi_indx)=(!b[k++]);
             indi_indx++;
             i++;
           }
        }
       snp_indx++;
     }
    outputbed(_snp_2,_snp_1,output);
    BIT.clear();
    BIT.close();
 }





//If the phenotype and genotype are in plink format, convert it to dosage manner
void prepare_geno_phe(string input, string X_file, string YFile)
 {
   cout<<"Come into function to prepare genotype and phenotype"<<endl;
   int i,j,k;
   string bedfile=input+".bed";
   string famfile=input+".fam";
   string bimfile=input+".bim";
   cout<<"bedfile is: "<<bedfile<<endl;
   cout<<"famfile is: "<<famfile<<endl;
   cout<<"bimfile is: "<<bimfile<<endl;
   int nsnp=readbim(bimfile);
   int nind=readfam(famfile);
   vector<string> variant;
   readbim(bimfile,variant);
   mat _snp_2 =mat(nsnp,nind, fill::zeros);
   mat _snp_1 =mat(nsnp,nind, fill::zeros);
   cout<<"nsnp is: "<<nsnp<<endl;
   cout<<"nind is: "<<nind<<endl;
   char ch[1];
   bitset<8> b;
   fstream BIT(bedfile.c_str(), ios::in|ios::binary);
   if(!BIT) throw("Error: can not open the file ["+input+"] to read.");
   for(i=0; i<3; i++) BIT.read(ch,1); // skip the first three bytes
   int snp_indx=0, indi_indx=0;
   for(j=0, snp_indx=0; j<nsnp; j++)
    {
      for(i=0, indi_indx=0; i<nind;)
       {
         BIT.read(ch,1);
         if(!BIT) throw("Error: problem with the BED file ... has the FAM/BIM file been changed?");
         b=ch[0];
         k=0;
         while(k < 7 && i < nind)
          {
            _snp_2(snp_indx,indi_indx)=(!b[k++]);
            _snp_1(snp_indx,indi_indx)=(!b[k++]);
            indi_indx++;
            i++;
          }
       }
      snp_indx++;
    }
   cout<<"Reading PLINK BED file is over"<<endl;
   mat dosage=mat(nsnp,nind, fill::zeros);
   for(int i=0;i<_snp_2.n_rows;i++)
    {
      for(int j=0;j<_snp_2.n_cols;j++)
       {
         if(_snp_2(i,j)==0 && _snp_1(i,j)==1)
          {
            dosage(i,j)=-9;
          } else
          {
            dosage(i,j)=_snp_2(i,j)+_snp_1(i,j);
          }
       }
    }
   for(int i=0;i<dosage.n_rows;i++)
    {
      double Mean=0;
      int num=0;
      for(int j=0;j<dosage.n_cols;j++)
       {
         if(dosage(i,j)!=-9)
          {
            Mean+=dosage(i,j);
            num++;
          }
       }
      Mean/=num;
      for(int j=0;j<dosage.n_cols;j++)
       {
         if(dosage(i,j)==-9)
          {
            dosage(i,j)=Mean;
          }
       }
    }
   ofstream outfile1(X_file.c_str(), ios::out);
   outfile1<<variant[0];
   for(int x=1;x<variant.size();x++)
    {
      outfile1<<" "<<variant[x];
    }
   outfile1<<endl;
   for(int i=0;i<dosage.n_cols;i++)
    {
      outfile1<<dosage(0,i);
      for(int j=1;j<dosage.n_rows;j++)
       {
         outfile1<<" "<<dosage(j,i);
       }
      outfile1<<endl;
    }
   outfile1.close();

   vector<double>  pheno;
   vector < vector<string> > fam_infor;
   readfam(famfile,pheno,fam_infor);
   ofstream outfile2(YFile.c_str(), ios::out);
   for (int i = 0; i < pheno.size(); i++)
     outfile2 << pheno[i] << endl;
   outfile2.close();

 }

//Just for test
int conditional_analysis(string input,string gene, int totalCausalSNP,float rho, bool histFlag)
 {
   cout<<"test"<<endl;
 }


//Perform conditional analysis to detect peak signal
int conditional_analysis(string file, string gene, int totalCausalSNP,float rho, bool histFlag, double gamma,string weight, int nthread,string covariate, vector <string> &grm_file,string outputFile,string geno_file)
 {

   vector< vector<double> >  Data;
   vector<string> ind;
   vector<string> variant;
   cout<<"Read input file"<<endl;
   cout<<"file is: "<<file<<endl;
   cout<<"explored gene is: "<<gene<<endl;
//   return (0);
   read_table(file,Data,ind,variant);
   cout<<"Reading input file is finished"<<endl;
   mat data=mat(Data.size(),Data[0].size(),fill::zeros);
   string out_test="";
   string fam_file_test="";
   for(int i=0;i<Data.size();i++)
    {
      for(int j=0;j<Data[i].size();j++)
       {
         data(i,j)=Data[i][j];
       }
    }
   cout<<"Convert input vector into mat finished"<<endl;
   int len=data.n_cols;
   int dim=data.n_rows;
   string out_text="";
   vector<string> col;
   for(int i=0;i<variant.size();i++)
    {
      col.push_back(variant[i]);
    }
   int j=0;
   for(int fea=0;fea<col.size();fea++)
    {
      if (regex_match (col[fea],regex("(IL)(.*)") ) || regex_match (col[fea],regex("(phe)(.*)") ) )
       {
         j=j+1;
       }
    }
   double cutoff=0.00001;
   cout<<"It is safe"<<endl;
   cout<<"Output individuals"<<endl;

   cout<<"j is: "<<j<<endl;
   for(int t=0;t<j;t++)
    {
      cout<<"Step1"<<endl;
      cout<<"I am a"<<endl;
      int x=0;
      int bi=0;
      vector<string> left;
      int num_rev=0;
      vector<string> removal;
      int target=0;
      double p_target;
      mat data1=data;
      vector<string> variant1;
      for(int i=0;i<variant.size();i++)
       {
         variant1.push_back(variant[i]);
       }
      string bed_file=string("probe")+to_string(t);

      cout<<"Copy genotype to each probe"<<endl;
      cout<<"Copy bed file"<<endl;
      copybed(geno_file+".bed",gene+"/"+bed_file+".bed");
      cout<<"Copy fam file"<<endl;
      copyfam(geno_file+".fam",gene+"/"+bed_file+".fam");
      cout<<"Copy bim file"<<endl;
      copybim(geno_file+".bim",gene+"/"+bed_file+".bim");

      string fam=gene+"/"+bed_file+".fam";
      vector<double>  pheno;
      vector < vector<string> > fam_infor;
      readfam(fam,pheno,fam_infor);
      mat fam_file=mat(pheno.size(),1,fill::zeros);
      for(int i=0;i<pheno.size();i++)
       {
         fam_file(i,0)=pheno[i];
       }
//      cout<<"fam_file is: "<<fam_file<<endl;
      int BIAO=fam_file.n_rows;
//      cout<<"BIAO is: "<<BIAO<<endl;
      BIAO=data.n_rows;
//      cout<<"BIAO is: "<<BIAO<<endl;
      fam_file.col(0)=data.col(t);
      string fam_file_out=fam+"_adjust";
      cout<<"updatefam"<<endl;
//      updatefam(fam_infor,fam_file,fam_file_out);
      updatefam(fam_infor,fam_file,fam);
      cout<<"I am Ok here"<<endl;
      for(int k=j;k<len;k++)
       {         
//         cout<<"Step2"<<endl;
//         cout<<"I am in"<<endl;
         vector<string> col1;
         for(int i=0;i<variant.size();i++)
          {
            col1.push_back(variant1[i]);
          }
         int len1=col1.size();
         data1.col(t)=data.col(t);
         target=0;
         p_target=1;
         int alle_target=0;
         string file_gemma_out=string("output_of_GEMMA_")+gene;
         PARAM cPar;
         GEMMA cGemma;
         cout<<"Inititation for GEMMA parameter and GEMMA core Over"<<endl;
//         cGemma.Assign(argc, argv, cPar);
         string str=gene+"/"+bed_file;
         cPar.file_bfile=str;

         str.assign(grm_file[0]);
         cPar.file_kin=str;
         str.assign(file_gemma_out);
         cPar.file_out=str;
         cPar.a_mode=1;
         cout<<"Perform parameter checking"<<endl;
         cPar.CheckParam();
         cout<<"Parameter check is over"<<endl;
         cGemma.BatchRun(cPar);
         cout<<"GEMMA running is over"<<endl;



         string file_gemma=string("./output")+"/output_of_GEMMA_"+gene+".assoc.txt";
         cout<<"Output of gemma is: "<<file_gemma<<endl;
         vector< vector<double> > Out_gemma;
         vector<string> Out_gemma_variant;
         readgemma(file_gemma,Out_gemma,Out_gemma_variant);
         mat out_gemma=mat(Out_gemma.size(),Out_gemma[0].size(),fill::zeros);
         for(int i=0;i<Out_gemma.size();i++)
          {
            for(int j=0;j<Out_gemma[i].size();j++)
             {
               out_gemma(i,j)=Out_gemma[i][j];
             }
          }

         int target_gemma  =index_min(out_gemma.col(1));
         cout<<"position for detected SNP is: "<<target_gemma<<endl;
         cout<<"detected SNP in GEMMA is: "<<Out_gemma_variant[target_gemma]<<endl;
         target=match(col1,Out_gemma_variant[target_gemma]);
         p_target=out_gemma(target_gemma,1);
         double beta_target =out_gemma(target_gemma,0);
         cout<<"beta value is: "<<beta_target<<endl;
         for(int X=0;X<removal.size();X++)
          {
            cout<<"removal is: "<<removal[X]<<endl;
          }
         cout<<"p_target is: "<<p_target<<", target is: "<<target<<" and SNP is: "<<Out_gemma_variant[target_gemma]<<endl;
         cout<<"If there are some peak signals detected, remove"<<endl;
         while(removal.size()>1&& match(removal,Out_gemma_variant[target_gemma])!=-1)
          {
             cout<<"Please skip this SNP "<<Out_gemma_variant[target_gemma]<<endl;
             out_gemma.shed_row(target_gemma);
             Out_gemma_variant.erase(Out_gemma_variant.begin()+target_gemma);
             cout<<"Try to find the index of minimum p value"<<endl;
             target_gemma  =index_min(out_gemma.col(1));
             target=match(col1,Out_gemma_variant[target_gemma]);
             p_target=out_gemma(target_gemma,1);
             if(p_target>cutoff)
              {
                break;
              }
             beta_target =out_gemma(target_gemma,0);
          }
         cout<<"p_target is:"<<p_target<<", target is: "<<target<<" and SNP is: "<<Out_gemma_variant[target_gemma]<<endl;
         if(p_target>cutoff)
          {
            cout<<"No significant variant found"<<endl;
            break ;
          }
         if(p_target<=cutoff)
          {
             cout<<"bi is: "<<bi<<endl;
             if(bi>0)
              {
                cout<<"Some peak signals are already detected"<<endl;
                int y=0;
                vector <string> left_test;
                for(int test=0;test<left.size();test++)
                 {
                   left_test.push_back(left[test]);
                 }
                left_test.push_back(col1[target]);
//                cout<<"left_test is: "<<left_test<<endl;
                vector<string>::iterator z;
                for(z=left_test.begin();z!=left_test.end();z++)
                 {
                   vector<string> col_test;
                   left_test.erase(z);
                   for(int i=0;i<left_test.size();i++)
                     {
                       col_test.push_back(left_test[i]);
                     }
                   cout<<"Come to calculate VIF"<<endl;
                   double vif=calculate_VIF(data,variant,(*z),col_test);
                   cout<<"VIF value is: "<<vif<<endl;
                   if(vif>=10)
                    {
                      y=1;
                    }
                 }
              
          
                if(y==1)
                 {
                   cout<<"Sorry, skip this SNP, as it locates in high LD with detected variants"<<endl;
                   num_rev++;
                   removal.push_back(col1[target]);
                   data1.shed_row(target);
                   col1.erase(col1.begin()+target);
                   int z=0;
                   while(y==1)
                    {
                      y=0;
                      int test=index_min(out_gemma.col(1));
                      out_gemma.shed_row(test);
                      int target=match(col1,Out_gemma_variant[target_gemma]);
                      while(target<0)
                       {
                         out_gemma.shed_row(target_gemma);
                         target_gemma=index_min(out_gemma.col(1));
                         target=match(col1,Out_gemma_variant[target_gemma]);
                         cout<<"Please skip this SNP, as it is in the list of removal"<<endl;
                       }
                      p_target=out_gemma(target_gemma,1);
                      beta_target=out_gemma(target_gemma,0);
                      if(p_target>cutoff)
                       {
                         z=1;
                         break;
                       }
                      if(p_target<=cutoff)
                       {
                            
                            y=0;
                            vector<string>::iterator z;
                            for(z=left_test.begin();z!=left_test.end();z++)
                             {
                               vector<string> col_test;
                               left_test.erase(z);
                               for(int i=0;i<left_test.size();i++)
                                {
                                  col_test.push_back(left_test[i]);
                                }
                               double vif=calculate_VIF(data,variant,(*z),col_test);
                               cout<<"VIF value is: "<<vif<<endl;
                               if(vif>=10)
                                {
                                  y=1;
                                }
                             }

                            if(y==1)
                             {
                               cout<<"Sorry, skip this SNP, as it locates in high LD with detected variants"<<endl;
                               num_rev=num_rev+1;
                               removal.push_back(col1[target]);
                               data1.shed_row(target);
                               col1.erase(col1.begin()+target);
                             }
                       }
                    }
                  if(z==1)
                   {
                     continue;
                   }
                  }
                int x=1;              
                if(x==1)
                 {
                  cout<<"Conglatulations, another e-signal detected"<<endl;
                  num_rev=num_rev+1;
                  removal.push_back(col1[target]);
                 
                  ostringstream streamObj;
                  streamObj << p_target;
                  string str_p_target = streamObj.str();

                  string res=col1[t]+" "+col1[target]+" "+str_p_target+" "+to_string(beta_target);
                  extractvariant(gene+"/"+bed_file+".bed", gene+"/"+bed_file+"_high_LD.bed", 0.3, col1[target]);
                  string yFile=gene+"/"+bed_file+"_high_LD.phe";
                  string outputFileName=gene+"/"+col1[t]+"_"+to_string(bi)+"_high_LD.PolyQTL.output";
                  string X_file=gene+"/"+bed_file+"_high_LD.geno";
                  prepare_geno_phe(gene+"/"+bed_file+"_high_LD",X_file,yFile);
                  MeQTLPolyGModel  MeQTLPoly(yFile, outputFileName, totalCausalSNP,rho, histFlag, gamma,weight, nthread,covariate, grm_file,X_file);
                  MeQTLPoly.run();
                  MeQTLPoly.finishUp();


                  out_test=out_test+"\n"+res;
                  cout<<"Get residuals and save it to fam file"<<endl;
                  cout<<"beta value is: "<<beta_target<<endl;
                  fam_file.col(0)=fam_file.col(0)- beta_target*data1.col(target);
//                  fam_file.col(0)=fam_file.col(0)- beta_target*data.col(target+j-1);
                  updatefam(fam_infor,fam_file,fam);
                  string fam_file_out=fam+"_adjust";
                  fam_file_test=fam+"_"+to_string(num_rev);
                  data1.shed_col(target);
                  bi++;
                  left.push_back(col1[target]);
                  col1.erase(col1.begin()+target);

                  len--;
                  if(len<=0)
                   {
                     break;
                   }                  
                }
              cout<<"Next running"<<endl;
              } else
              {
                  num_rev=num_rev+1;
                  removal.push_back(col1[target]);
                  ostringstream streamObj;
                  streamObj << p_target;
                  string str_p_target = streamObj.str();
                  
                  string res=col1[t]+" "+col1[target]+" "+str_p_target+" "+to_string(beta_target);
                  extractvariant(gene+"/"+bed_file+".bed", gene+"/"+bed_file+"_high_LD.bed", 0.3, col1[target]);
                  string yFile=gene+"/"+bed_file+"_high_LD.phe";
                  string outputFileName=gene+"/"+col1[t]+"_"+to_string(bi)+"_high_LD.PolyQTL.output";
                  string X_file=gene+"/"+bed_file+"_high_LD.geno";
                  prepare_geno_phe(gene+"/"+bed_file+"_high_LD",X_file,yFile);
                  MeQTLPolyGModel  MeQTLPoly(yFile, outputFileName, totalCausalSNP,rho, histFlag, gamma,weight, nthread,covariate, grm_file,X_file);
                  MeQTLPoly.run();
                  MeQTLPoly.finishUp(); 
                  out_test=res;
                  cout<<"Get residuals and save it to fam file"<<endl;
                  cout<<"beta value is: "<<beta_target<<endl;

                  fam_file.col(0)=fam_file.col(0)- beta_target*data1.col(target);
                  updatefam(fam_infor,fam_file,fam);
                  string fam_file_out=fam+"_adjust";
                  fam_file_test=fam+"_"+to_string(num_rev);

                  data1.shed_col(target);
                  bi++;
                  left.push_back(col1[target]);
                  col1.erase(col1.begin()+target);
//                  cout<<"left is: "<<left<<endl;
                  len--;
                  if(len<=0)
                   {
                     break;
                   }
//                  return (0);
              }                     
           }
        }           
      cout<<"conditional analysis result is: "<<out_test<<endl;
      fstream  FHOU(outputFile.c_str(),ios::out);
      if(!FHOU)
       {
         cout<<"Error for input or output"<<endl;
         return (1);
       }
      FHOU<<out_test<<endl;
      FHOU.clear();
      FHOU.close();
      
//      cout<<"left is: "<<left<<endl
      
   }
 }
