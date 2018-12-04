#include <string>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <map>
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

//------- calc mean and vari -----//
pair<double,double> Calc_Mean_Vari(vector<double> & score_list)
{
	int count=(int)score_list.size();
	if(count<=0)
	{
		pair<double,double> result(1,0);
		return result;
	}
	//final stat
	int i;
	double mean,vari;
	mean=0;
	for(i=0;i<count;i++)mean+=score_list[i];
	mean/=count;
	vari=0;
	for(i=0;i<count;i++)vari+=(score_list[i]-mean)*(score_list[i]-mean);
	vari=1.0*sqrt(vari/count);
	//return
	pair<double,double> result(mean,vari);
	return result;
}


//------ main ------//
int main(int argc,char **argv)
{
	//----- Stat List ---//
	{
		if(argc<2)
		{
			printf("Version: 1.02 \n");
			printf("Stat_List <score_list> \n");
			exit(-1);
		}
		string infile=argv[1];
		vector <double> score_list;
		Read_Score(infile,score_list);
		pair<double,double> param = Calc_Mean_Vari(score_list);
		printf("mean = %lf , vari = %lf \n",param.first,param.second);
		exit(0);

	}
}
