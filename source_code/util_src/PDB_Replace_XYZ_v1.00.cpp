#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <map>
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


//------------- load Rosetta model -----------------//
//-> example
/*
ATOM      1  N   MET A   1     -13.432  56.028  15.302  1.00  0.00           N
ATOM      2  CA  MET A   1     -13.637  57.376  15.819  1.00  0.00           C
ATOM      3  C   MET A   1     -12.712  58.373  15.129  1.00  0.00           C
ATOM      4  O   MET A   1     -12.077  58.051  14.124  1.00  0.00           O
ATOM      5  CB  MET A   1     -13.416  57.392  17.333  1.00  0.00           C
ATOM      6  CG  MET A   1     -14.341  56.471  18.114  1.00  0.00           C
ATOM      7  SD  MET A   1     -16.075  56.944  17.965  1.00  0.00           S
ATOM      8  CE  MET A   1     -16.099  58.448  18.937  1.00  0.00           C
.....
ATOM   5338 1HH2 ARG A 341      30.826  78.687 -16.915  1.00  0.00           H
ATOM   5339 2HH2 ARG A 341      32.426  79.364 -16.709  1.00  0.00           H
TER
# All scores below are weighted scores, not raw scores.
#BEGIN_POSE_ENERGIES_TABLE S_0001.pdb
label fa_atr fa_rep fa_sol fa_intra_rep pro_close fa_pair hbond_sr_bb hbond_lr_bb hbond_bb_sc hbond_sc dslf_ss_dst dslf_cs_ang dslf_ss_dih dslf_ca_dih atom_pair_constraint rama omega fa_dun p_aa_pp ref cart_bonded total
weights 0.8 0.44 0.65 0.004 1 0.49 0.585 1.17 1.17 1.1 0.5 2 5 5 0.5 0.2 0.5 0.56 0.32 1 0.5 NA

*/
int Load_Rosetta_Model(string &pdbfile,
	vector <vector <string> > &head_rec,
	vector <vector <string> > &xyz_rec) 
{
	map<string, int>::iterator iter;
	map<string, int> mapping;
	ifstream fin;
	string buf,temp,name;
	//read
	fin.open(pdbfile.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"pdbfile %s not found!!\n",pdbfile.c_str());
		exit(-1);
	}
	int count=0;
	head_rec.clear();
	xyz_rec.clear();
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		int len=(int)buf.length();
		if(len<3)continue;
		//check TER
		temp=buf.substr(0,3);
		if(temp=="TER"||temp=="END")break;
		//check ATOM
		if(len<4)continue;
		temp=buf.substr(0,4);
		if(temp!="ATOM" && temp!="HETA")continue;
		//record name
		name=buf.substr(21,6);
		iter = mapping.find(name);
		if(iter != mapping.end())    //-> add additional residue
		{
			int key=mapping[name]-1;
			//-> atom content
			string atomtype=buf.substr(0,30);
			head_rec[key].push_back(atomtype);
			string coordinate=buf.substr(30,24);
			xyz_rec[key].push_back(coordinate);
		}
		else                         //-> add a new residue
		{
			count++;
			mapping.insert(map < string, int >::value_type(name, count));
			//add new residue
			//-> atom content
			vector <string> atomrec;
			string atomtype=buf.substr(0,30);
			atomrec.push_back(atomtype);
			head_rec.push_back(atomrec);
			//-> coordinate
			vector <string> tmprec;
			string coordinate=buf.substr(30,24);
			tmprec.push_back(coordinate);
			xyz_rec.push_back(tmprec);
		}
	}
	//return
	return count;
}


//---------------- replcae coordinate from 'start' position -----------------//
//-> start_ is 0-base
//-> example
/*
ATOM   2768  N   SER   338     -20.777   0.877  18.287  1.00 18.16           N
ATOM   2769  CA  SER   338     -20.011   0.145  17.319  1.00 18.16           C
ATOM   2770  CB  SER   338     -20.014   0.808  15.932  1.00 18.16           C
ATOM   2771  OG  SER   338     -19.439   2.103  16.003  1.00 18.16           O
ATOM   2772  C   SER   338     -18.594   0.029  17.795  1.00 18.16           C
ATOM   2773  O   SER   338     -17.951  -1.002  17.605  1.00 18.16           O
*/
void Replace_Coordinate(string &pdbfile, FILE *fp, int start_, 
	vector <vector <string> > &head_rec,
	vector <vector <string> > &xyz_rec)
{
	map<int, int>::iterator iter;
	map<int, int> mapping;
	ifstream fin;
	string buf,temp,name;
	//read
	fin.open(pdbfile.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"pdbfile %s not found!!\n",pdbfile.c_str());
		exit(-1);
	}
	int start=start_+1;
	int end=start_+xyz_rec.size();
	int count=0;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		int len=(int)buf.length();
		if(len>=54)
		{
			temp=buf.substr(0,4);
			if(temp=="ATOM"||temp=="HETA")
			{
				//----- find start -----//
				int pos=atoi(buf.substr(21,5).c_str());
				if(pos<start || pos>end)fprintf(fp,"%s\n",buf.c_str());
				else
				{
					iter = mapping.find(pos);
					if(iter == mapping.end())    //-> add a new residue
					{
						count++;
						mapping.insert(map < int, int >::value_type(pos, count));
						//get head and tail
						string head=buf.substr(0,30);
						string tail=buf.substr(54,len-54);
						//replace coordinate
						for(int i=0;i<(int)xyz_rec[count-1].size();i++)
							fprintf(fp,"%s%s%s\n",head_rec[count-1][i].c_str(),xyz_rec[count-1][i].c_str(),tail.c_str());
					}
				}
			}
			else fprintf(fp,"%s\n",buf.c_str());
		}
		else fprintf(fp,"%s\n",buf.c_str());
	}
	//---- final check -----//
	if(count!=xyz_rec.size())
	{
		fprintf(stderr,"count %d not equal to xyz_rec.size() %d \n",
			count,xyz_rec.size());
		exit(-1);
	}
}

//---------- main ----------//
int main(int argc,char **argv)
{
	//---- PDB_Replace_XYZ ----//
	{
		if(argc<5)
		{
			fprintf(stderr,"Version 1.00 \n");
			fprintf(stderr,"PDB_Replace_XYZ <pdbfile> <start_pos> <replcae_pdb> <outfile> \n");
			fprintf(stderr,"[note]: start_pos shall be 0-base \n");
			exit(-1);
		}
		string pdbfile=argv[1];
		int start_pos=atoi(argv[2]);
		string replace_file=argv[3];
		string outfile=argv[4];
		//--- load replace_file ---//
		vector <vector <string> > head_rec;
		vector <vector <string> > xyz_rec;
		Load_Rosetta_Model(replace_file,head_rec,xyz_rec);
		//--- replace coordinate ---//
		FILE *fp=fopen(outfile.c_str(),"wb");
		Replace_Coordinate(pdbfile,fp,start_pos,head_rec,xyz_rec);
		fclose(fp);
		exit(0);
	}
}

