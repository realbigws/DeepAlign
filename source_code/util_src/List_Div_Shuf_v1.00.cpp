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
#include <algorithm>
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



//------- shuffle and cut -------//
void Shuffle_And_Cut_List(string &list,string &outnam,int num, string &out_root)
{
	ifstream fin;
	string buf,temp;
	fin.open(list.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"list file not found [%s] !!!\n",list.c_str());
		exit(-1);
	}
	vector <string> shuffle_list;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		shuffle_list.push_back(buf);
	}
	//shuffle list
	random_shuffle(shuffle_list.begin(),shuffle_list.end());
	//cut list
	int i;
	int totnum=(int)shuffle_list.size();
	int bignum=totnum/num;
	int minnum=totnum%num;
	int count=0;
	int iter=0;
	FILE *fp;
	string name;
	char command[30000];
	sprintf(command,"%s/%s.%d",out_root.c_str(),outnam.c_str(),iter);
	name=command;
	fp=fopen(name.c_str(),"wb");
	for(i=0;i<totnum;i++)
	{
		if(count==bignum && iter<num-1)
		{
			fclose(fp);
			iter++;
			sprintf(command,"%s/%s.%d",out_root.c_str(),outnam.c_str(),iter);
			name=command;
			fp=fopen(name.c_str(),"wb");
			count=0;
		}
		fprintf(fp,"%s\n",shuffle_list[i].c_str());
		count++;
	}
	fclose(fp);
	//judge
	if(iter!=num-1)
	{
		for(i=iter+1;i<=num;i++)
		{
			sprintf(command,"%s/%s.%d",out_root.c_str(),outnam.c_str(),i);
			name=command;
			fp=fopen(name.c_str(),"wb");
			fclose(fp);
		}
	}
}


//---------- main ---------//
int main(int argc,char **argv)
{
	//---- shuffle and divide list ----//
	{
		if(argc<4)
		{
			printf("List_Div_Shuf <list> <div_num> <out_root> \n");
			exit(-1);
		}
		string list=argv[1];
		int div_num=atoi(argv[2]);
		string out_root=argv[3];
		//process
		string name;
		getBaseName(list,name,'/','.');
		Shuffle_And_Cut_List(list,name,div_num, out_root);
		//exit
		exit(0);
	}
}

