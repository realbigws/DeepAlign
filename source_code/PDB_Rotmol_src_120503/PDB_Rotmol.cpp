#include <string>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include "Kabsch.h"
using namespace std;

//----- read rotmat -----//
void Load_Rotmat(string &fn,double *rotmat,int skip=0)
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(fn.c_str(),ios::in);
	if(fin.fail()!=0)
	{
		printf("no such file![%s]\n",fn.c_str());
		exit(-1);
	}
	//skip
	int i;
	for(i=0;i<skip;i++)
	{
		if(!getline(fin,buf,'\n'))
		{
			printf("file error![%s]\n",fn.c_str());
			exit(-1);
		}
	}
	//load
	string l1,l2,l3;
	if(!getline(fin,l1,'\n'))
	{
		printf("file error![%s]\n",fn.c_str());
		exit(-1);
	}
	if(!getline(fin,l2,'\n'))
	{
		printf("file error![%s]\n",fn.c_str());
		exit(-1);
	}
	if(!getline(fin,l3,'\n'))
	{
		printf("file error![%s]\n",fn.c_str());
		exit(-1);
	}
	//process
	temp=l1.substr(0,9);
	rotmat[0]=atof(temp.c_str());
	temp=l1.substr(10,9);
	rotmat[1]=atof(temp.c_str());
	temp=l1.substr(20,9);
	rotmat[2]=atof(temp.c_str());
	temp=l1.substr(30,12);
	rotmat[9]=atof(temp.c_str());
	temp=l2.substr(0,9);
	rotmat[3]=atof(temp.c_str());
	temp=l2.substr(10,9);
	rotmat[4]=atof(temp.c_str());
	temp=l2.substr(20,9);
	rotmat[5]=atof(temp.c_str());
	temp=l2.substr(30,12);
	rotmat[10]=atof(temp.c_str());
	temp=l3.substr(0,9);
	rotmat[6]=atof(temp.c_str());
	temp=l3.substr(10,9);
	rotmat[7]=atof(temp.c_str());
	temp=l3.substr(20,9);
	rotmat[8]=atof(temp.c_str());
	temp=l3.substr(30,12);
	rotmat[11]=atof(temp.c_str());
}

//----- Process_PDB ----//
void Process_PDB(string &fn,double *rotmat,FILE *fp)
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(fn.c_str(),ios::in);
	if(fin.fail()!=0)
	{
		printf("no such file![%s]\n",fn.c_str());
		exit(-1);
	}
	//process
	int len;
	Kabsch kabsch;
	XYZ ori,nxt;
	string head,tail;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		len=(int)buf.length();
		if(len>=54)
		{
			temp=buf.substr(0,4);
			if(temp=="ATOM"||temp=="HETA")
			{
				head="";
				head=buf.substr(0,30);
				tail="";
				tail=buf.substr(54,len-54);
				temp=buf.substr(30,8);
				ori.X=atof(temp.c_str());
				temp=buf.substr(38,8);
				ori.Y=atof(temp.c_str());
				temp=buf.substr(46,8);
				ori.Z=atof(temp.c_str());
				kabsch.rot_point(ori,nxt,rotmat);
				fprintf(fp,"%s%8.3f%8.3f%8.3f%s\n",head.c_str(),nxt.X,nxt.Y,nxt.Z,tail.c_str());
			}
			else fprintf(fp,"%s\n",buf.c_str());
		}
		else fprintf(fp,"%s\n",buf.c_str());
	}
}


//------- main -------//
int main(int argc,char **argv)
{
	//---- real_process ----//
	{
		if(argc<5)
		{
			fprintf(stderr,"PDB_Rotmol <pdb_file> <rotmat> <out> <skip>\n");
			exit(-1);
		}
		//arguments
		string pdb_file=argv[1];
		string rot_file=argv[2];
		string out_file=argv[3];
		int skip=atoi(argv[4]);
		//process
		double rotmat[12];
		Load_Rotmat(rot_file,rotmat,skip);
		FILE *fp=fopen(out_file.c_str(),"wb");
		Process_PDB(pdb_file,rotmat,fp);
		fclose(fp);
	}
}

