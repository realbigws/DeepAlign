#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
using namespace std;


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

//------- int to string --------//
template <typename T>
string NumberToString( T Number )
{
	ostringstream ss;
	ss << Number;
	return ss.str();
}


//-------- PDB_Residue_Copy ---------//
//-> input PDB_1 (source)
/*
ATOM      1  N   CYS A 162      -4.982  26.111  19.390  1.00 25.54      A    N
ATOM      2  CA  CYS A 162      -5.505  27.238  20.155  1.00 29.50      A    C
ATOM      3  C   CYS A 162      -4.430  28.295  20.422  1.00 29.50      A    C
ATOM      4  O   CYS A 162      -4.750  29.450  20.711  1.00 31.66      A    O
ATOM      5  CB  CYS A 162      -6.103  26.756  21.475  1.00 25.09      A    C
ATOM      6  SG  CYS A 162      -4.911  25.951  22.567  1.00 20.61      A    S
ATOM      7  N   ASN A 163      -3.160  27.887  20.343  1.00 22.52      A    N
ATOM      8  CA  ASN A 163      -2.029  28.816  20.425  1.00 20.46      A    C
...
*/

//-> input PDB_2 (target)
/*
ATOM      1  N   CYS A   1      -4.982  26.111  19.390  1.00 25.54      A    N
ATOM      2  CA  CYS A   1      -5.505  27.238  20.155  1.00 29.50      A    C
ATOM      3  CB  CYS A   1      -6.103  26.756  21.475  1.00 25.09      A    C
ATOM      4  SG  CYS A   1      -4.911  25.951  22.567  1.00 20.61      A    S
ATOM      5  C   CYS A   1      -4.430  28.295  20.422  1.00 29.50      A    C
ATOM      6  O   CYS A   1      -4.750  29.450  20.711  1.00 31.66      A    O
ATOM      7  N   ASN A   2      -3.160  27.887  20.343  1.00 22.52      A    N
ATOM      8  CA  ASN A   2      -2.029  28.816  20.425  1.00 20.46      A    C
...
*/

//our purpose is to COPY the residue number from PDB_1 to PDB_2

//----- read source PDB --------//
int PDB_Residue_Read(string &pdb,vector <string> &resi)
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(pdb.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"list %s not found!!\n",pdb.c_str());
		exit(-1);
	}
	int len;
	int first=1;
	string prev="";
	string curr;
	resi.clear();
	int count=0;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		len=(int)buf.length();
		if(len<3)continue;
		//check TER
		temp=buf.substr(0,3);
		if(temp=="TER"||temp=="END")break;
		//check ATOM
		if(len<4)continue;
		temp=buf.substr(0,4);
		if(temp!="ATOM"&&temp!="HETA")continue;
		//resi
		curr=buf.substr(21,6);
		if(first==1)
		{
			first=0;
			prev=curr;
		}
		if(curr!=prev)
		{
			resi.push_back(prev);
			count++;
			prev=curr;
		}
	}
	if(prev!="")
	{
		resi.push_back(curr);
		count++;
	}
	return count;
}


//----- process target PDB --------//
void PDB_Residue_Copy(string &pdb,FILE *fp,vector <string> &resi)
{
	ifstream fin;
	string buf,temp,name;
	//read
	fin.open(pdb.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"list %s not found!!\n",pdb.c_str());
		exit(-1);
	}
	int len;
	int count=-1;
	string orirec="";
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		len=(int)buf.length();
		if(len>=3)
		{
			temp=buf.substr(0,3);
			if(temp=="TER")continue;
		}
		if(len>=4)
		{
			temp=buf.substr(0,4);
			if(temp=="ATOM"||temp=="HETA"||temp=="MISS")
			{
				temp=buf.substr(21,6);
				if(temp!=orirec)
				{
					orirec=temp;
					count++;
				}
				string pos_str=resi[count];
				string part1=buf.substr(0,21);
				string part2=buf.substr(27,len-27);
				fprintf(fp,"%s%6s%s\n",part1.c_str(),pos_str.c_str(),part2.c_str());
			}
			else fprintf(fp,"%s\n",buf.c_str());
		}
		else fprintf(fp,"%s\n",buf.c_str());
	}
}


//-------- main ---------//
int main(int argc,char **argv)
{
	//---- PDB_Residue_Trans ----//
	{
		if(argc<4)
		{
			fprintf(stderr,"Version 1.00 \n");
			fprintf(stderr,"PDB_Resi_Copy <source_pdb> <target_pdb> <output_pdb> \n");
			fprintf(stderr,"[note]: copy the residue index from source_pdb to target_pdb. \n");
			fprintf(stderr,"        the length of the two input PDBs must be the same. \n");
			exit(-1);
		}
		string source_pdb=argv[1];
		string target_pdb=argv[2];
		string output_pdb=argv[3];
		//read residue
		vector <string> source_residue;
		int source_len=PDB_Residue_Read(source_pdb,source_residue);
		vector <string> target_residue;
		int target_len=PDB_Residue_Read(target_pdb,target_residue);
		//check length
		if(source_len!=target_len)
		{
			fprintf(stderr,"source_len %d not equal to target_len %d \n",
				source_len,target_len);
			exit(-1);
		}
		//copy residue
		FILE *fp=fopen(output_pdb.c_str(),"wb");
		PDB_Residue_Copy(target_pdb,fp,source_residue);
		fclose(fp);
		exit(0);
	}
}

