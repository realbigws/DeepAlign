#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <string.h>
using namespace std;


//======================= I/O related ==========================//
//-------- utility ------//
void getBaseName(string &in,string &out,char slash,char dot)
{
	int i,j;
	int len=(int)in.length();
	for(i=len-1;i>=0;i--)
	{
		if(in[i]==slash)break;
	}
	i++;
	for(j=len-1;j>=0;j--)
	{
		if(in[j]==dot)break;
	}
	if(j==-1)j=len;
	out=in.substr(i,j-i);
}
void getRootName(string &in,string &out,char slash)
{
	int i;
	int len=(int)in.length();
	for(i=len-1;i>=0;i--)
	{
		if(in[i]==slash)break;
	}
	if(i<=0)out=".";
	else out=in.substr(0,i);
}


//------- read EVD parameter: miu and beta ---------//
//-> example
/*
miu = %lf , beta = %lf
*/
void Read_Miu_Beta(string &infile, double &miu,double &beta)
{
	ifstream fin;
	string buf;
	//load
	fin.open(infile.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"file %s not found!!\n",infile.c_str());
		exit(-1);
	}
	//read
	if(!(fin>>buf))
	{
		fprintf(stderr,"file %s bad!!\n",infile.c_str());
		exit(-1);
	}
	if(!(fin>>buf))
	{
		fprintf(stderr,"file %s bad!!\n",infile.c_str());
		exit(-1);
	}
	if(!(fin>>buf))
	{
		fprintf(stderr,"file %s bad!!\n",infile.c_str());
		exit(-1);
	}
	miu=atof(buf.c_str());
	if(!(fin>>buf))
	{
		fprintf(stderr,"file %s bad!!\n",infile.c_str());
		exit(-1);
	}
	if(!(fin>>buf))
	{
		fprintf(stderr,"file %s bad!!\n",infile.c_str());
		exit(-1);
	}
	if(!(fin>>buf))
	{
		fprintf(stderr,"file %s bad!!\n",infile.c_str());
		exit(-1);
	}
	if(!(fin>>buf))
	{
		fprintf(stderr,"file %s bad!!\n",infile.c_str());
		exit(-1);
	}
	beta=atof(buf.c_str());
}

//----------- calculate Pvalue related -----------//
double Pvalue(double x, double lamda,double miu)
{
	double h = lamda * (x - miu);
	return (h > 10) ? exp(-h) : double(1.0) - exp(-exp(-h));
}
int Get_Score_List(string &file,vector <double> &score_list)
{
	ifstream fin;
	string buf;
	fin.open(file.c_str(), ios::in);
	if(fin.fail()!=0)return 0;
	int count=0;
	score_list.clear();
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		double score=atof(buf.c_str());
		score_list.push_back(score);
		count++;
	}
	return count;
}
int Get_Cutoff_Given_EVD(vector <double> &score, double miu, double beta, double thres)
{
	int i;
	int size=(int)score.size();
	for(i=0;i<size;i++)
	{
		double ws_pvalue_rank=Pvalue(score[i],1/beta,miu);
		if(ws_pvalue_rank>thres)return i;
	}
}


//---------- main ---------//
int main(int argc,char **argv)
{
	//---- determine TopK by pValue ----//
	{
		if(argc<4)
		{
			printf("TopK_by_EVD <sorted_score> <evd_param> <pval_cut> \n");
			exit(-1);
		}
		string sorted_score=argv[1];
		string evd_param=argv[2];
		double pval_cut=atof(argv[3]);
		//process
		int TopK=-1;
		if(pval_cut>0)
		{
			//-> read sorted score
			vector <double> pvalue_score;
			Get_Score_List(sorted_score,pvalue_score);
			//-> read EVD parameter
			double miu=0;
			double beta=1;
			Read_Miu_Beta(evd_param,miu,beta);
			//-> determine TopK by pVal cutoff
			TopK=Get_Cutoff_Given_EVD(pvalue_score,miu,beta,pval_cut);
		}
		printf("%d\n",TopK);
		//exit
		exit(0);
	}
}

