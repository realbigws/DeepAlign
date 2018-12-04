#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstring>
#include <time.h>
#include <math.h>
#include <algorithm>
using namespace std;



//--------- read score ------//
//only read the first part
void Read_Score(string &infile,vector <double> &out_score)
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(infile.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"list %s not found!!\n",infile.c_str());
		exit(-1);
	}
	out_score.clear();
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		istringstream www(buf);
		www>>temp;
		double value=atof(temp.c_str());
		if(value<0)value=0;
		out_score.push_back(value);
	}
}

//--------- Pvalue related ------//
inline double logPvalue(double x, double lamda,double miu) 
{
	double h = lamda * (x - miu);
	//return log((double(1.0) - exp(-exp(-h))));
	//cout<<(double(1.0)-exp(-exp(-h)))<<" ";
	return (h > 10) ? -h : (h < -2.5) ? -exp(-exp(-h)) : log((double(1.0)-exp(-exp(-h))));
}
inline double Pvalue(double x, double lamda,double miu) 
{
	double h = lamda * (x - miu);
	return (h > 10) ? exp(-h) : double(1.0) - exp(-exp(-h));
}

//------- fit gumbel extreme value distribution --------//
int ws_function(vector<double> & score_list,double lamda,double mean,double &lamda_new)
{
	double zero_momen=0;
	double fist_momen=0;
	double send_momen=0;
	int size=(int)score_list.size();
	for(int j=0;j<size;j++)
	{
		double x=score_list[j];
		double exp_value = exp(-1.0*lamda*x);
		zero_momen += exp_value;
		fist_momen += x*exp_value;
		send_momen += x*x*exp_value;
	}
	double tmp=fist_momen/zero_momen;
	double f=1.0/lamda-mean+tmp;
	double g=tmp*tmp-send_momen/zero_momen-1/lamda/lamda;
	lamda_new=lamda-f/g;
	//--- return check ---//
	if(fabs(f)<0.001)return 1;  //converged
	return 0;
}

pair<double,double> fit_gumbel_extreme(vector<double> & score_list)
{
	int size=(int)score_list.size();
	if(size<=0)
	{
		pair<double,double> result(1,0);
		return result;
	}
	//--- calculate mean & vari ---//
	double mean = 0;
	for(int i=0;i<size;i++)
		mean += score_list[i];
	mean = mean/size;
	double vari = 0;
	for(int i=0;i<size;i++)
		vari += (score_list[i]-mean)*(score_list[i]-mean);
	vari= sqrt(1.0*vari/score_list.size());

	//--- fitting parameter ---//
	double lamda=1.0/vari;
	double lamda_new=lamda;
	for(int i=0;i<1000;i++)
	{
		int retv=ws_function(score_list,lamda,mean,lamda_new);
		if(retv==1)break;
		lamda=lamda_new;
	}
	double zero_momen=0;
	for(int j=0;j<size;j++)
	{
		double x=score_list[j];
		double exp_value = exp(-1.0*lamda*x);
		zero_momen += exp_value;
	}
	double miu = 1.0/lamda_new * log(size/zero_momen);
	double beta = 1.0/lamda_new;
	
	//--- final check ----//
	if(fabs(miu-mean)>1.5*vari)
	{
		fprintf(stderr,"Warning: evd fitting not converged !!, mean=%lf , vari=%lf ; miu=%lf , beta=%lf \n",
			mean,vari,miu,beta);
		miu = mean;
		beta = vari;
	}
	pair<double,double> result(miu,beta);
	return result;
}


//--------- main ----------//
int main(int argc,char **argv)
{
	//---- Fitting EVD parameter ----//
	{
		if(argc<2)
		{
			fprintf(stderr,"Version: 1.32 \n");
			fprintf(stderr,"Fitting_EVD <score_list> \n");
			exit(-1);
		}
		string infile=argv[1];
		vector <double> score_list;
		Read_Score(infile,score_list);
		pair<double,double> param = fit_gumbel_extreme (score_list);
		printf("miu = %lf , beta = %lf \n",param.first,param.second);
		exit(0);
	}
}

