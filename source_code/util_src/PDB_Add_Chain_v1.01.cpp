#include <string>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>
using namespace std;


//------- int2string -----//
string convertInt(int number)
{
	stringstream ss;  //create a stringstream
	ss << number;     //add number to the stream
	return ss.str();  //return a string with the contents of the stream
}


//========== process PDB ============//
void PDB_Add_Chain(string &pdbfile,FILE *fp,char chain) //->from .pdb file
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(pdbfile.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"pdbfile %s not found!!\n",pdbfile.c_str());
		exit(-1);
	}
	int pos;
	int len;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		len=(int)buf.length();
		if(len>=4)
		{
			temp=buf.substr(0,4);
			if(temp=="ATOM"||temp=="HETA"||temp=="MISS")
			{
				buf[21]=chain;
				fprintf(fp,"%s\n",buf.c_str());
			}
			else fprintf(fp,"%s\n",buf.c_str());
		}
		else fprintf(fp,"%s\n",buf.c_str());
	}
}

//---------- main ----------//
int main(int argc,char **argv)
{
	//---- PDB_Residue_Trans ----//
	{
		if(argc<4)
		{
			fprintf(stderr,"Version 1.01 \n");
			fprintf(stderr,"PDB_Add_Chain <pdbfile> <chain> <outfile> \n");
			fprintf(stderr,"[note]: set chain to -1 to remove \n");
			exit(-1);
		}
		string infile=argv[1];
		string chain=argv[2];
		string outfile=argv[3];
		//--- check chain ---//
		char chain_char;
		if(chain=="-1")chain_char=' ';
		else chain_char=chain[0];
		//--- process chain ---//
		FILE *fp=fopen(outfile.c_str(),"wb");
		PDB_Add_Chain(infile,fp,chain_char);
		fclose(fp);
		exit(0);
	}
}
