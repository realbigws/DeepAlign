#include <time.h>
#include <sstream>
#include <map>
#include "getopt.h"
#include "Mol_Load.h"
#include "Mol_Out.h"
#include "CLEFAPS_Main.h"
#include "CLEPAPS_Out.h"
#include "reference_nam.h"
#include "reference_mca.h"
#include "reference_mcb.h"
#include "reference_ami.h"
#include "reference_cle.h"
#include "Fast_Sort.h"
#include "Confo_Beta.h"
#define REFERENCE_NUM 1814
using namespace std;

Mol_Load mol_input;   //-> for database search
Mol_Load mol_input1;  //-> for load first mol
Mol_Load mol_input2;  //-> for load second mol
Mol_Out mol_out;
Confo_Lett confo_lett;

//-> [1] for output PDB in linear co-ordinate, otherwise not output
int Out_PDB;          //-> output PDB



//--------- FASTA I/O ------------//
//FASTA
int ReadToFile_FASTA(string &fn,vector<pair<int, int> > &alignment,
	string &nam1_content,string &nam2_content,
	string &nam1_full,string &nam2_full,
	string &nam1,string &nam2)
{
	int i;
	int cur1=0;
	int cur2=0;
	int len;
	int len1,len2;
	alignment.clear();
	//init
	string seq="";  //sequence
	string tmp="";  //template
	//load
	ifstream fin;
	string buf,temp;
	fin.open(fn.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"alignment file not found [%s] !!!\n",fn.c_str());
		return -1;
	}
	//read tmp
	for(;;)
	{
		if(!getline(fin,buf,'\n'))goto badend;
		len=(int)buf.length();
		if(len>1)
		{
			if(buf[0]=='>')
			{
				istringstream www(buf);
				www>>temp;
				len=(int)temp.length();
				nam1=temp.substr(1,len-1);
				break;
			}
		}
	}
	for(;;)
	{
		if(!getline(fin,buf,'\n'))goto badend;
		len=(int)buf.length();
		if(len==0)continue;
		if(len>1)
		{
			if(buf[0]=='>')
			{
				istringstream www(buf);
				www>>temp;
				len=(int)temp.length();
				nam2=temp.substr(1,len-1);
				break;
			}
		}
		tmp+=buf;
	}
	//read seq
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		len=(int)buf.length();
		if(len==0)continue;
		seq+=buf;
	}
	//process
	len1=(int)seq.length();
	len2=(int)tmp.length();
	if(len1!=len2)
	{
		fprintf(stderr,"alignment len not equal [%s] !!!\n",fn.c_str());
		return -1;
	}
	len=len1;
	nam1_content.clear();
	nam2_content.clear();
	for(i=0;i<len;i++)
	{
		if(tmp[i]!='-' && seq[i]!='-') //match
		{
			nam1_content.push_back(tmp[i]);
			nam2_content.push_back(seq[i]);
			cur1++;
			cur2++;
			alignment.push_back(pair<int,int>(cur1,cur2));
		}
		else
		{
			if(tmp[i]!='-') //Ix
			{
				nam1_content.push_back(tmp[i]);
				cur1++;
				alignment.push_back(pair<int,int>(cur1,-cur2));
			}
			if(seq[i]!='-') //Iy
			{
				nam2_content.push_back(seq[i]);
				cur2++;
				alignment.push_back(pair<int,int>(-cur1,cur2));
			}
		}
	}
	//return
	nam1_full=tmp;
	nam2_full=seq;
	return 1; //success
	
	badend:
	fprintf(stderr,"alignment file format bad [%s] !!!\n",fn.c_str());
	return -1;
}


//============= Exttract Alignment ========//
//-> get mapping alignment detail
void process_oriami_record_simp(const char *seq_,const char *ami_,
	vector<pair<int,int> > &WWW_alignment)
{
	int i,j;
	int n1,n2;
	//--[1]dynamic_programming	
	n1=(int)strlen(seq_);
	n2=(int)strlen(ami_);
	double *WWW_score=new double[n1*n2];
	for(i=0;i<n1;i++)
	{
		for(j=0;j<n2;j++)
		{
			if(seq_[i]==ami_[j])
			{
				if(seq_[i]=='X'||seq_[i]=='Z'||seq_[i]=='.')WWW_score[i*n2+j]=0;
				else if(ami_[j]=='X'||ami_[j]=='Z'||ami_[j]=='.')WWW_score[i*n2+j]=0;
				else WWW_score[i*n2+j]=10;
			}
			else
			{
				if(seq_[i]=='X'||seq_[i]=='Z'||seq_[i]=='.')WWW_score[i*n2+j]=0;
				else if(ami_[j]=='X'||ami_[j]=='Z'||ami_[j]=='.')WWW_score[i*n2+j]=0;
				else WWW_score[i*n2+j]=-15;
			}
		}
	}
	double sco;
	int matchs;
	matchs=Advance_Align_Dyna_Prog_Double(n1,n2,WWW_score,-11,-1,-11,-1,0,0,0,0,
		WWW_alignment,sco);
	delete [] WWW_score;
}
//-> extract alignment
//[note]: we regard namX_pdb as alignment_out
//        the namX_seq as alignment_in
//[example]:
/*
//         AAAAABBBBBCCCCDDDDD  -> 1 (nam1_pdb) \ -> ali1 \
//              BBBBBCCCC       -> 2 (nam1_seq)  \         \
//               ||||||                          -> ali2   -> ali4
//               BBBBCC         -> 3 (nam2_seq)  /         /
//          AAAABBBBBCCCCDD     -> 4 (nam2_pdb) / -> ali3 /
*/
//[note]: ali2 is from alignment_in, and ali4 is from alignment_out
void Extract_Alignment(string &nam1_seq,string &nam1_pdb,string &nam2_seq,string &nam2_pdb,
	vector<pair<int,int> > &alignment_in,vector<pair<int,int> > &alignment_out)
{
	//init
	int i;
	int l1,l2,l3,l4;
	l1=(int)nam1_pdb.length();
	l2=(int)nam1_seq.length();
	l3=(int)nam2_seq.length();
	l4=(int)nam2_pdb.length();
	int *ali1=new int[l1];
	int *ali2=new int[l2];
	int *ali3=new int[l3];
	int *ali4=new int[l1];
	vector<pair<int,int> > WWW_alignment;
	//first alignment
	process_oriami_record_simp(nam1_pdb.c_str(),nam1_seq.c_str(),WWW_alignment);
	for(i=0;i<l1;i++)ali1[i]=-1;
	for(i=0;i<(int)WWW_alignment.size();i++)
	{
		int pos1=WWW_alignment[i].first;
		int pos2=WWW_alignment[i].second;
		if(pos1>0 && pos2>0)ali1[pos1-1]=pos2-1;
	}
	//second alignment
	for(i=0;i<l2;i++)ali2[i]=-1;
	for(i=0;i<(int)alignment_in.size();i++)
	{
		int pos1=alignment_in[i].first;
		int pos2=alignment_in[i].second;
		if(pos1>0 && pos2>0)ali2[pos1-1]=pos2-1;
	}
	//third alignment
	process_oriami_record_simp(nam2_seq.c_str(),nam2_pdb.c_str(),WWW_alignment);
	for(i=0;i<l3;i++)ali3[i]=-1;
	for(i=0;i<(int)WWW_alignment.size();i++)
	{
		int pos1=WWW_alignment[i].first;
		int pos2=WWW_alignment[i].second;
		if(pos1>0 && pos2>0)ali3[pos1-1]=pos2-1;
	}
	//fourth alignment
	for(i=0;i<l1;i++)ali4[i]=-1;
	for(i=0;i<l1;i++)
	{
		int pos1=ali1[i];
		if(pos1==-1)continue;
		int pos2=ali2[pos1];
		if(pos2==-1)continue;
		int pos3=ali3[pos2];
		if(pos3==-1)continue;
		ali4[i]=pos3;
	}
	//final
	{
		alignment_out.clear();
		//start
		int j;
		int ii,jj;
		int wlen;
		int pre_ii=0;
		int pre_jj=0;
		for(i=1;i<=l1;i++)
		{
			ii=i;
			jj=ali4[i-1];  //ali1 starts from 0, correspondence also from 0
			if(jj==-1)
			{
				continue;
			}
			else
			{
				jj++;
				//previous_path
				wlen=ii-pre_ii;
				for(j=1;j<wlen;j++)
				{
					pre_ii++;
					alignment_out.push_back (pair<int,int>(pre_ii, -pre_jj)); //Ix
				}
				wlen=jj-pre_jj;
				for(j=1;j<wlen;j++)
				{
					pre_jj++;
					alignment_out.push_back (pair<int,int>(-pre_ii, pre_jj)); //Iy
				}
				//current_path
				alignment_out.push_back (pair<int,int>(ii, jj)); //Match
				//update
				pre_ii=ii;
				pre_jj=jj;
			}
		}
		//termi
		pre_ii++;
		for(i=pre_ii;i<=l1;i++)alignment_out.push_back (pair<int,int>(i, -pre_jj)); //Ix
		pre_jj++;
		for(i=pre_jj;i<=l4;i++)alignment_out.push_back (pair<int,int>(-l1, i));  //Iy
	}
	//delete
	delete [] ali1;
	delete [] ali2;
	delete [] ali3;
	delete [] ali4;
}



//=========== blosum calculate ============//
double Local_Para=10;
//=================//
//--Ori_BLOSUM----//
//===============//
int Ori_BLOSUM_62[21][21]={
{  4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -5 },  //A
{ -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -5 },  //R
{ -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3, -5 },  //N
{ -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3, -5 },  //D
{  0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -5 },  //C
{ -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2, -5 },  //Q
{ -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2, -5 },  //E
{  0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -5 },  //G
{ -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3, -5 },  //H
{ -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -5 },  //I
{ -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -5 },  //L
{ -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2, -5 },  //K
{ -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -5 },  //M
{ -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -5 },  //F
{ -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -5 },  //P
{  1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2, -5 },  //S
{  0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -5 },  //T
{ -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -5 },  //W
{ -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -5 },  //Y
{  0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -5 },  //V
{ -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5 }}; //Z
// A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   Z

int Ori_BLOSUM_45[21][21]={
{  5, -2, -1, -2, -1, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -2, -2,  0, -5 },  //A
{ -2,  7,  0, -1, -3,  1,  0, -2,  0, -3, -2,  3, -1, -2, -2, -1, -1, -2, -1, -2, -5 },  //R
{ -1,  0,  6,  2, -2,  0,  0,  0,  1, -2, -3,  0, -2, -2, -2,  1,  0, -4, -2, -3, -5 },  //N
{ -2, -1,  2,  7, -3,  0,  2, -1,  0, -4, -3,  0, -3, -4, -1,  0, -1, -4, -2, -3, -5 },  //D
{ -1, -3, -2, -3, 12, -3, -3, -3, -3, -3, -2, -3, -2, -2, -4, -1, -1, -5, -3, -1, -5 },  //C
{ -1,  1,  0,  0, -3,  6,  2, -2,  1, -2, -2,  1,  0, -4, -1,  0, -1, -2, -1, -3, -5 },  //Q
{ -1,  0,  0,  2, -3,  2,  6, -2,  0, -3, -2,  1, -2, -3,  0,  0, -1, -3, -2, -3, -5 },  //E
{  0, -2,  0, -1, -3, -2, -2,  7, -2, -4, -3, -2, -2, -3, -2,  0, -2, -2, -3, -3, -5 },  //G
{ -2,  0,  1,  0, -3,  1,  0, -2, 10, -3, -2, -1,  0, -2, -2, -1, -2, -3,  2, -3, -5 },  //H
{ -1, -3, -2, -4, -3, -2, -3, -4, -3,  5,  2, -3,  2,  0, -2, -2, -1, -2,  0,  3, -5 },  //I
{ -1, -2, -3, -3, -2, -2, -2, -3, -2,  2,  5, -3,  2,  1, -3, -3, -1, -2,  0,  1, -5 },  //L
{ -1,  3,  0,  0, -3,  1,  1, -2, -1, -3, -3,  5, -1, -3, -1, -1, -1, -2, -1, -2, -5 },  //K
{ -1, -1, -2, -3, -2,  0, -2, -2,  0,  2,  2, -1,  6,  0, -2, -2, -1, -2,  0,  1, -5 },  //M
{ -2, -2, -2, -4, -2, -4, -3, -3, -2,  0,  1, -3,  0,  8, -3, -2, -1,  1,  3,  0, -5 },  //F
{ -1, -2, -2, -1, -4, -1,  0, -2, -2, -2, -3, -1, -2, -3,  9, -1, -1, -3, -3, -3, -5 },  //P
{  1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -3, -1, -2, -2, -1,  4,  2, -4, -2, -1, -5 },  //S
{  0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -1, -1,  2,  5, -3, -1,  0, -5 },  //T
{ -2, -2, -4, -4, -5, -2, -3, -2, -3, -2, -2, -2, -2,  1, -3, -4, -3, 15,  3, -3, -5 },  //W
{ -2, -1, -2, -2, -3, -1, -2, -3,  2,  0,  0, -1,  0,  3, -3, -2, -1,  3,  8, -1, -5 },  //Y
{  0, -2, -3, -3, -1, -3, -3, -3, -3,  3,  1, -2,  1,  0, -3, -1,  0, -3, -1,  5, -5 },  //V
{ -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5 }}; //Z
// A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   Z  

int Ori_BLOSUM_80[21][21]={
{  5, -2, -2, -2, -1, -1, -1,  0, -2, -2, -2, -1, -1, -3, -1,  1,  0, -3, -2,  0, -5 },  //A
{ -2,  6, -1, -2, -4,  1, -1, -3,  0, -3, -3,  2, -2, -4, -2, -1, -1, -4, -3, -3, -5 },  //R
{ -2, -1,  6,  1, -3,  0, -1, -1,  0, -4, -4,  0, -3, -4, -3,  0,  0, -4, -3, -4, -5 },  //N
{ -2, -2,  1,  6, -4, -1,  1, -2, -2, -4, -5, -1, -4, -4, -2, -1, -1, -6, -4, -4, -5 },  //D
{ -1, -4, -3, -4,  9, -4, -5, -4, -4, -2, -2, -4, -2, -3, -4, -2, -1, -3, -3, -1, -5 },  //C
{ -1,  1,  0, -1, -4,  6,  2, -2,  1, -3, -3,  1,  0, -4, -2,  0, -1, -3, -2, -3, -5 },  //Q
{ -1, -1, -1,  1, -5,  2,  6, -3,  0, -4, -4,  1, -2, -4, -2,  0, -1, -4, -3, -3, -5 },  //E
{  0, -3, -1, -2, -4, -2, -3,  6, -3, -5, -4, -2, -4, -4, -3, -1, -2, -4, -4, -4, -5 },  //G
{ -2,  0,  0, -2, -4,  1,  0, -3,  8, -4, -3, -1, -2, -2, -3, -1, -2, -3,  2, -4, -5 },  //H
{ -2, -3, -4, -4, -2, -3, -4, -5, -4,  5,  1, -3,  1, -1, -4, -3, -1, -3, -2,  3, -5 },  //I
{ -2, -3, -4, -5, -2, -3, -4, -4, -3,  1,  4, -3,  2,  0, -3, -3, -2, -2, -2,  1, -5 },  //L
{ -1,  2,  0, -1, -4,  1,  1, -2, -1, -3, -3,  5, -2, -4, -1, -1, -1, -4, -3, -3, -5 },  //K
{ -1, -2, -3, -4, -2,  0, -2, -4, -2,  1,  2, -2,  6,  0, -3, -2, -1, -2, -2,  1, -5 },  //M
{ -3, -4, -4, -4, -3, -4, -4, -4, -2, -1,  0, -4,  0,  6, -4, -3, -2,  0,  3, -1, -5 },  //F
{ -1, -2, -3, -2, -4, -2, -2, -3, -3, -4, -3, -1, -3, -4,  8, -1, -2, -5, -4, -3, -5 },  //P
{  1, -1,  0, -1, -2,  0,  0, -1, -1, -3, -3, -1, -2, -3, -1,  5,  1, -4, -2, -2, -5 },  //S
{  0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -2, -1, -1, -2, -2,  1,  5, -4, -2,  0, -5 },  //T
{ -3, -4, -4, -6, -3, -3, -4, -4, -3, -3, -2, -4, -2,  0, -5, -4, -4, 11,  2, -3, -5 },  //W
{ -2, -3, -3, -4, -3, -2, -3, -4,  2, -2, -2, -3, -2,  3, -4, -2, -2,  2,  7, -2, -5 },  //Y
{  0, -3, -4, -4, -1, -3, -3, -4, -4,  3,  1, -3,  1, -1, -3, -2,  0, -3, -2,  4, -5 },  //V
{ -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5 }}; //Z
// A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   Z  


//BLOSUM_Mapping//--------------ARNDCQEGHILKMFPSTWYVZ
int Blo_AA_Map_WS[21]=
{ 0,19, 4, 3, 6, 13,7, 8, 9, 17,11,10,12,2, 18,14,5, 1, 15,16,20};
//A  V  C  D  E  F  G  H  I  W  K  L  M  N  Y  P  Q  R   S  T  Z
//0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17  18 19 20
//Ori_Mapping//-----------------AVCDEFGHIWKLMNYPQRSTZ
int Ori_AA_Map_WS[26]=
{ 0,20,2,3,4,5,6,7,8,20,10,11,12,13,20,15,16,17,18,19,20, 1, 9,20,14,20};
// A B C D E F G H I  J  K  L  M  N  O  P  Q  R  S  T  U  V  W  X  Y  Z
// 0 1 2 3 4 5 6 7 8  9 10 11 12 14 14 15 16 17 18 19 20 21 22 23 24 25

//------ calculate -------//
int BLOSUM62_Calc(char a,char b)
{
	int ii,jj;
	if(a<'A' || a>'Z')a='Z';
	ii=Blo_AA_Map_WS[Ori_AA_Map_WS[a-'A']];
	if(b<'A' || b>'Z')b='Z';
	jj=Blo_AA_Map_WS[Ori_AA_Map_WS[b-'A']];
	return Ori_BLOSUM_62[ii][jj];
}
int BLOSUM45_Calc(char a,char b)
{
	int ii,jj;
	if(a<'A' || a>'Z')a='Z';
	ii=Blo_AA_Map_WS[Ori_AA_Map_WS[a-'A']];
	if(b<'A' || b>'Z')b='Z';
	jj=Blo_AA_Map_WS[Ori_AA_Map_WS[b-'A']];
	return Ori_BLOSUM_45[ii][jj];
}
int BLOSUM80_Calc(char a,char b)
{
	int ii,jj;
	if(a<'A' || a>'Z')a='Z';
	ii=Blo_AA_Map_WS[Ori_AA_Map_WS[a-'A']];
	if(b<'A' || b>'Z')b='Z';
	jj=Blo_AA_Map_WS[Ori_AA_Map_WS[b-'A']];
	return Ori_BLOSUM_80[ii][jj];
}

//================ clesum calculate =============//
int Ori_CLESUM_WS[18][18]={
{ 73, 20, 13, -17, -25, -20, -6, -45, -31,-23,-19,-11, -2, 10, 25, 35, 16,0}, //A 
{ 20, 51,  7,  13,  15,   7, 13, -96, -74,-57,-50,-12,-13,-11,-12, 42, 12,0}, //B 
{ 13,  7, 53,  21,   3,  20, -4, -77, -56,-43,-33,  0,-12, -5,  3,  4, 29,0}, //C 
{-17, 13, 21,  52,  22,  22,-31,-124,-105,-88,-81,-22,-49,-44,-42,-10, 14,0}, //D 
{-25, 15,  3,  22,  36,  26,-22,-127,-108,-93,-84,-21,-47,-43,-48, -5, -6,0}, //E 
{-20,  7, 20,  22,  26,  50, -5,-107, -88,-73,-69,-16,-33,-32,-30,  0,  3,0}, //F 
{ -6, 13, -4, -31, -22,  -5, 69, -51, -34,-21,-13, 29, 21, -8, -1,  5,  8,0}, //G 
{-45,-96,-77,-124,-127,-107,-51,  23,  18, 13,  5,-62, -4,-34,-55,-60,-87,0}, //H 
{-31,-74,-56,-105,-108, -88,-34,  18,  23, 16, 21,-41,  1,-11,-34,-49,-62,0}, //I 
{-23,-57,-43, -88, -93, -73,-21,  13,  16, 37, 13,-32, 16, -2,-24,-34,-44,0}, //J 
{-19,-50,-33, -81, -84, -69,-13,   5,  21, 13, 49, -1, 12, 28,  5,-36,-24,0}, //K 
{-11,-12,  0, -22, -21, -16, 29, -62, -41,-32, -1, 74,  5,  8, -4,-12, 26,0}, //L 
{ -2,-13,-12, -49, -47, -33, 21,  -4,   1, 16, 12,  5, 61,  7,  5,  8, -7,0}, //M 
{ 10,-11, -5, -44, -43, -32, -8, -34, -11, -2, 28,  8,  7, 90, 15, -3, 32,0}, //N 
{ 25,-12,  3, -42, -48, -30, -1, -55, -34,-24,  5, -4,  5, 15,104,  4,-13,0}, //O 
{ 35, 42,  4, -10,  -5,   0,  5, -60, -49,-34,-36,-12,  8, -3,  4, 66,  7,0}, //P 
{ 16, 12, 29,  14,  -6,   3,  8, -87, -62,-44,-24, 26, -7, 32,-13,  7, 90,0}, //Q 
{  0,  0,  0,   0,   0,   0,  0,   0,   0,  0,  0,  0,  0,  0,  0,  0,  0,0}};//R  
// A   B   C    D    E    F   G    H    I   J   K   L   M   N   O   P   Q R

//------ calculate -------//
int CLESUM_Calc(char a,char b)
{
	if(a<'A' || a>'R')a='R';
	if(b<'A' || b>'Z')b='R';
	return Ori_CLESUM_WS[a-'A'][b-'A'];
}



//-----------------------------------------------------------------------------------------------------------//
//-------------- Single_Output -------------//
//-> protein structure -------//
char *TM_IND1,*TM_IND2;
char *TM_AMI1,*TM_AMI2;
char *TM_CLE1,*TM_CLE2;
XYZ *TM_MOL1,*TM_MOL2,*TM_MOL_TMP;
PDB_Residue *TM_PDB1,*TM_PDB2,*TM_PDB_TMP;
int TM_MOLN1,TM_MOLN2;
//-> protein alignment ------//
int *TM_ALIGNMENT;
int *TM_ALI1;
int *TM_AFP;
double *TM_ROTMAT;
//-> CB_additional ----//__110830__//
char *tmp_ami;
XYZ *tmp_mcb;
int *ws_reco;
XYZ *TM_MCB1,*TM_MCB2;
//-> weighted matrix ---//
int *wwa1;
int *wwa2;
int *wwc1;
int *wwc2;
double *wswei;

//-------- Create -------//
void TM_Align_Init_WS(int maxlen1,int maxlen2)
{
	//-- protein structure --//
	//-> local
	TM_IND1=new char[6*maxlen1];
	TM_IND2=new char[6*maxlen2];
	TM_AMI1=new char[maxlen1+1];
	TM_AMI2=new char[maxlen2+1];
	TM_CLE1=new char[maxlen1+1];
	TM_CLE2=new char[maxlen2+1];
	//-> global
	TM_MOL1=new XYZ[maxlen1];
	TM_MOL2=new XYZ[maxlen2];
	TM_MOL_TMP=new XYZ[maxlen1+maxlen2];
	TM_PDB1=new PDB_Residue[maxlen1];
	TM_PDB2=new PDB_Residue[maxlen2];
	TM_PDB_TMP=new PDB_Residue[maxlen1+maxlen2];
	//-- protein alignment --//
	TM_ALIGNMENT=new int[maxlen1+maxlen2];
	TM_ALI1=new int[maxlen1+maxlen2];
	TM_AFP=new int[(maxlen1+maxlen2)*4];
	TM_ROTMAT=new double[12];
	//---- CB_additional ----//__110830__//
	tmp_ami=new char[maxlen1+maxlen2+1];
	tmp_mcb=new XYZ[maxlen1+maxlen2];
	ws_reco=new int[maxlen1+maxlen2];
	TM_MCB1=new XYZ[maxlen1];
	TM_MCB2=new XYZ[maxlen2];
	//---- others ----//
	wwa1=new int[maxlen1];
	wwa2=new int[maxlen2];
	wwc1=new int[maxlen1];
	wwc2=new int[maxlen2];
	wswei=new double[maxlen1*maxlen2];
}
void TM_Align_dele_WS(void)
{
	//---- protein structure ---//
	//-> local 
	delete [] TM_IND1;
	delete [] TM_IND2;
	delete [] TM_AMI1;
	delete [] TM_AMI2;
	delete [] TM_CLE1;
	delete [] TM_CLE2;
	//-> global
	delete [] TM_MOL1;
	delete [] TM_MOL2;
	delete [] TM_MOL_TMP;
	delete [] TM_PDB1;
	delete [] TM_PDB2;
	delete [] TM_PDB_TMP;
	//---- protein alignment ---//
	delete [] TM_ALIGNMENT;
	delete [] TM_ALI1;
	delete [] TM_AFP;
	delete [] TM_ROTMAT;
	//---- CB_additional ----//__110830__//
	delete [] tmp_ami;
	delete [] tmp_mcb;
	delete [] ws_reco;
	delete [] TM_MCB1;
	delete [] TM_MCB2;
	//---- others ----//
	delete [] wwa1;
	delete [] wwa2;
	delete [] wwc1;
	delete [] wwc2;
	delete [] wswei;
}

//=========== Fasta_Output related ===========//
void FASTA_Output(FILE *fp,string &nam1,string &nam2,char *ami1,char *ami2,int *ali1,int *ali2,int moln1,int moln2)
{
	//ali2->ali1
	int i;
	int ii,jj;
	for(i=0;i<moln1;i++)ali1[i]=-1;
	for(i=0;i<moln2;i++)
	{
		if(ali2[i]>=0)
		{
			ii=ali2[i];
			jj=i;
			ali1[ii]=jj;
		}
	}

	//ali1->alignment
	//init
	vector<pair<int,int> > alignment;
	alignment.clear();
	//start
	int j;
	int wlen;
	int pre_ii=0;
	int pre_jj=0;
	for(i=1;i<=moln1;i++)
	{
		ii=i;
		jj=ali1[i-1];  //ali1 starts from 0, correspondence also from 0
		if(jj==-1)
		{
			continue;
		}
		else
		{
			jj++;
			//previous_path
			wlen=ii-pre_ii;
			for(j=1;j<wlen;j++)
			{
				pre_ii++;
				alignment.push_back (pair<int,int>(pre_ii, -pre_jj)); //Ix
			}
			wlen=jj-pre_jj;
			for(j=1;j<wlen;j++)
			{
				pre_jj++;
				alignment.push_back (pair<int,int>(-pre_ii, pre_jj)); //Iy
			}
			//current_path
			alignment.push_back (pair<int,int>(ii, jj)); //Match
			//update
			pre_ii=ii;
			pre_jj=jj;
		}
	}
	//termi
	pre_ii++;
	for(i=pre_ii;i<=moln1;i++)alignment.push_back (pair<int,int>(i, -pre_jj)); //Ix
	pre_jj++;
	for(i=pre_jj;i<=moln2;i++)alignment.push_back (pair<int,int>(-moln1, i));  //Iy

	//alignment->output
	//output
	char c;
	int size=(int)alignment.size();
	fprintf(fp,">%s\n",nam1.c_str());
	for(i=0;i<size;i++)
	{
		ii=alignment[i].first;
		if(ii<=0)fprintf(fp,"-");
		else
		{
			c=ami1[ii-1];
			fprintf(fp,"%c",c);
//			if(c=='X'||c=='Z')fprintf(fp,".");
//			else fprintf(fp,"%c",c);
		}
	}
	fprintf(fp,"\n");
	fprintf(fp,">%s\n",nam2.c_str());
	for(i=0;i<size;i++)
	{
		jj=alignment[i].second;
		if(jj<=0)fprintf(fp,"-");
		else
		{
			c=ami2[jj-1];
			fprintf(fp,"%c",c);
//			if(c=='X'||c=='Z')fprintf(fp,".");
//			else fprintf(fp,"%c",c);
		}
	}
	fprintf(fp,"\n");
}
void FASTA_Output_More(string &ws_output_tot,string &nam1_,string &nam2_,
		char *ami1_,char *ami2_, char *cle1_,char *cle2_,XYZ *mol1,XYZ *mol2,
		int *ali1,int *ali2,int moln1,int moln2,double d0,int &seqid,
		int &blosum,int &clesum,vector <double> &match_wei,vector <double> &rmsd_out,int CLEorNOT, int threscut=80)
{
	ws_output_tot="";
	//ali2->ali1
	int i;
	int ii,jj;
	for(i=0;i<moln1;i++)ali1[i]=-1;
	for(i=0;i<moln2;i++)
	{
		if(ali2[i]>=0)
		{
			ii=ali2[i];
			jj=i;
			ali1[ii]=jj;
		}
	}

	//ali1->alignment
	//init
	vector<pair<int,int> > alignment;
	alignment.clear();
	//start
	int j;
	int wlen;
	int pre_ii=0;
	int pre_jj=0;
	for(i=1;i<=moln1;i++)
	{
		ii=i;
		jj=ali1[i-1];  //ali1 starts from 0, correspondence also from 0
		if(jj==-1)
		{
			continue;
		}
		else
		{
			jj++;
			//previous_path
			wlen=ii-pre_ii;
			for(j=1;j<wlen;j++)
			{
				pre_ii++;
				alignment.push_back (pair<int,int>(pre_ii, -pre_jj)); //Ix
			}
			wlen=jj-pre_jj;
			for(j=1;j<wlen;j++)
			{
				pre_jj++;
				alignment.push_back (pair<int,int>(-pre_ii, pre_jj)); //Iy
			}
			//current_path
			alignment.push_back (pair<int,int>(ii, jj)); //Match
			//update
			pre_ii=ii;
			pre_jj=jj;
		}
	}
	//termi
	pre_ii++;
	for(i=pre_ii;i<=moln1;i++)alignment.push_back (pair<int,int>(i, -pre_jj)); //Ix
	pre_jj++;
	for(i=pre_jj;i<=moln2;i++)alignment.push_back (pair<int,int>(-moln1, i));  //Iy


	//-------- get start and end ---------//
	int start1,start2;
	int end1,end2;
	int wsstart=-1;
	int wsend=-1;
	int x=0;
	int y=0;
	rmsd_out.clear();
	int size=(int)alignment.size();
	{
		//omit head gap
		for(i=0;i<size;i++)
		{
			x = alignment[i].first;
			y = alignment[i].second;
			wsstart=i;
			if(x>0 && y>0)break;
			rmsd_out.push_back(-1);
		}
		if(x==0 || y==0)return;
		x--;
		y--;
		start1=x;
		start2=y;
		end1=x;
		end2=y;
		//omit tail gap
		for(i=size-1;i>=0;i--)
		{
			x = alignment[i].first;
			y = alignment[i].second;
			wsend=i;
			if(x>0 && y>0)break;
		}
	}

	//---------- calc -----------//
	{
		//get name
		string nam1;
		if(nam1_.length()<14)nam1=nam1_;
		else nam1=nam1_.substr(0,14);
		string nam2;
		if(nam2_.length()<14)nam2=nam2_;
		else nam2=nam2_.substr(0,14);
		//others
		double wscur1,wscur2,wscur;
		int first=1;
		int lali=0;
		char wstmp;
		string cle1,ami1,cle2,ami2;
		string score_out;
		//string clear
		cle1.clear();
		ami1.clear();
		cle2.clear();
		ami2.clear();
		score_out.clear();
		//output related
		char ws_output[300000];
		string ws_output_str;
		ws_output_str.clear();
		//[real process]
		int num=0;
		blosum=0;
		clesum=0;
		seqid=0;
		match_wei.clear();
		for(i=wsstart; i <= wsend; i++)
		{
			x = alignment[i].first;
			y = alignment[i].second;
			x--,y--;
			//output
			if(num%threscut==0 && first==0)
			{
				//fragment output
				ws_output_str.clear();
				sprintf(ws_output,"\n");
				ws_output_str+=ws_output;
				if(CLEorNOT==1)
				{
					sprintf(ws_output,"T cle                 %s           \n",cle1.c_str());
					ws_output_str+=ws_output;
				}
				sprintf(ws_output,"T %-14s %4d %s %4d (%d)\n",
					nam1.c_str(),start1+1,ami1.c_str(),end1+1,moln1);
				ws_output_str+=ws_output;
				sprintf(ws_output,"RMSD                  %s           \n",score_out.c_str());
				ws_output_str+=ws_output;
				sprintf(ws_output,"S %-14s %4d %s %4d (%d)\n",
					nam2.c_str(),start2+1,ami2.c_str(),end2+1,moln2);
				ws_output_str+=ws_output;
				if(CLEorNOT==1)
				{
					sprintf(ws_output,"S cle                 %s           \n",cle2.c_str());
					ws_output_str+=ws_output;
				}
				sprintf(ws_output,"\n");
				ws_output_str+=ws_output;
				//string clear
				cle1.clear();
				ami1.clear();
				cle2.clear();
				ami2.clear();
				score_out.clear();
				//position update
				start1=end1+1;
				start2=end2+1;
				//final output
				ws_output_tot+=ws_output_str;
			}
			num++;
			//lali
			if(x>=0 && y>=0)lali++;
			if(x>=0)end1=x;
			if(y>=0)end2=y;
			first=0;
			//record
			//[temp]
			if(x>=0)
			{
				cle1=cle1+cle1_[x];
				ami1=ami1+ami1_[x];
			}
			else
			{
				cle1=cle1+'-';
				ami1=ami1+'-';
			}
			//[targ]
			if(y>=0)
			{
				cle2=cle2+cle2_[y];
				ami2=ami2+ami2_[y];
			}
			else
			{
				cle2=cle2+'-';
				ami2=ami2+'-';
			}
			//[score]
			if(x>=0 && y>=0)
			{
				//calculate local RMSD
				double dist2=mol1[x].distance_square(mol2[y]);
//				double tmscore=1.0/(1.0+dist2/d0/d0);
				double tmscore=sqrt(1.0*dist2);
				int wssco=(int)(tmscore);
				if(wssco<0)wssco=0;
				if(wssco>9)wssco=9;
				wstmp=wssco+'0';
				if(wstmp=='9')wstmp='.';
				score_out=score_out+wstmp;
				rmsd_out.push_back(tmscore);
				//calculate BLOSUM and CLESUM
				wscur1=BLOSUM62_Calc(ami1_[x],ami2_[y]);
				wscur2=CLESUM_Calc(cle1_[x],cle2_[y]);
				blosum+=(int)wscur1;
				clesum+=(int)wscur2;
				//calculate sequence identity
				if(ami1_[x]==ami2_[y])
				{
					if(ami1_[x]!='.' && ami1_[x]!='X' && ami1_[x]!='Z')seqid++;
				}
				//add to match_wei
				if(wscur1<0)wscur1=0;     //max(BLOSUM,0)
				wscur=wscur1+0.1*wscur2;   //local parameter
				match_wei.push_back(wscur);
			}
			else
			{
				score_out=score_out+" ";
				rmsd_out.push_back(-1);
			}
		}

		//rmsd_termi
		for(i=wsend+1;i<size;i++)rmsd_out.push_back(-1);

		//termi_iprocess
	//	if(i%threscut==0 && first==0)
		{
			//fragment output
			ws_output_str.clear();
			sprintf(ws_output,"\n");
			ws_output_str+=ws_output;
			if(CLEorNOT==1)
			{
				sprintf(ws_output,"T cle                 %s           \n",cle1.c_str());
				ws_output_str+=ws_output;
			}
			sprintf(ws_output,"T %-14s %4d %s %4d (%d)\n",
				nam1.c_str(),start1+1,ami1.c_str(),end1+1,moln1);
			ws_output_str+=ws_output;
			sprintf(ws_output,"RMSD                  %s           \n",score_out.c_str());
			ws_output_str+=ws_output;
			sprintf(ws_output,"S %-14s %4d %s %4d (%d)\n",
				nam2.c_str(),start2+1,ami2.c_str(),end2+1,moln2);
			ws_output_str+=ws_output;
			if(CLEorNOT==1)
			{
				sprintf(ws_output,"S cle                 %s           \n",cle2.c_str());
				ws_output_str+=ws_output;
			}
			sprintf(ws_output,"\n");
			ws_output_str+=ws_output;
			//string clear
			cle1.clear();
			ami1.clear();
			cle2.clear();
			ami2.clear();
			score_out.clear();
			//position update
			start1=end1+1;
			start2=end2+1;
			//final output
			ws_output_tot+=ws_output_str;
		}
	}
	clesum/=10;
}
void Molecular_Output(FILE *fp,XYZ *mol1,XYZ *mol2,char *ami1,char *ami2,int moln1,int moln2)
{
	string TER="TER                                                                             ";
	mol_out.Output_PDB(fp,moln1,mol1,ami1,0,'A');
	fprintf(fp,"%s\n",TER.c_str());
	mol_out.Output_PDB(fp,moln2,mol2,ami2,0,'B');
	fprintf(fp,"%s\n",TER.c_str());
}

//=============== single_pair_process ================//
/* an input file example
#start
!alignment file between sequence F0299 and template 1akea
!this alignment is generated by RAPTOR
!zScore=59.429024
>P1;1ake
structure:1ake:1:A:214:A::::
--MRIILLGAPGAGKGTQAQFIMEKYGIPQISTGDMLRAAVKSGSELGKQAKDIMDAGKL
V-----TDELVIALVKERIAQEDCRNGFLLDGFPRTIPQADAMKEAGINVD------YVL
EFD--VPDELIVDRIV-----GRRVHAPSGRVYHVKFNPPKVEGKDDVTGEELTTRKDDQ
EETVRKRLVEYHQMTAPLIGYYSKEAEAGNTKYAKVDGTKPVAEVRADLEKILG-*

>P1;F0299
sequence:F0299::::::::
MTRYALLVRGINVG------------GKNKVVMAE-LRQELTNLG--LEKVESYINSGNI
FFTSIDSKAQLVEKLETFFAVHYPFIQSFSLL--SLEDFEAELENLPAWWSRDLARKDFL
FYTEGLDVDQVIATVESLELKDEVLYFGKLGIFWGKFSEESYSKT---------------
-------------------AYHKYLLKVPFYRHITI---RNAKTF-DKIGQMLKK*
#end
*/
// PIR_Type -> [0] both structure, [1] first structure, [2] second structure
void PIR_Output(FILE *fp,char *ami1,char *ami2,int *ali1,int *ali2,int moln1,int moln2,
				   string &nam1,string &nam2,char c1,char c2,int PIR_Type)
{
	//ali2->ali1
	int i;
	int ii,jj;
	for(i=0;i<moln1;i++)ali1[i]=-1;
	for(i=0;i<moln2;i++)
	{
		if(ali2[i]>=0)
		{
			ii=ali2[i];
			jj=i;
			ali1[ii]=jj;
		}
	}

	//ali1->alignment
	//init
	vector<pair<int,int> > alignment;
	alignment.clear();
	//start
	int j;
	int wlen;
	int pre_ii=0;
	int pre_jj=0;
	for(i=1;i<=moln1;i++)
	{
		ii=i;
		jj=ali1[i-1];  //ali1 starts from 0, correspondence also from 0
		if(jj==-1)
		{
			continue;
		}
		else
		{
			jj++;
			//previous_path
			wlen=ii-pre_ii;
			for(j=1;j<wlen;j++)
			{
				pre_ii++;
				alignment.push_back (pair<int,int>(pre_ii, -pre_jj)); //Ix
			}
			wlen=jj-pre_jj;
			for(j=1;j<wlen;j++)
			{
				pre_jj++;
				alignment.push_back (pair<int,int>(-pre_ii, pre_jj)); //Iy
			}
			//current_path
			alignment.push_back (pair<int,int>(ii, jj)); //Match
			//update
			pre_ii=ii;
			pre_jj=jj;
		}
	}
	//termi
	pre_ii++;
	for(i=pre_ii;i<=moln1;i++)alignment.push_back (pair<int,int>(i, -pre_jj)); //Ix
	pre_jj++;
	for(i=pre_jj;i<=moln2;i++)alignment.push_back (pair<int,int>(-moln1, i));  //Iy

	//alignment->output
	//mol1 output
	fprintf(fp,">P1;%s\n",nam1.c_str());
	if(PIR_Type==0 || PIR_Type==1)fprintf(fp,"structure:%s::%c::%c::::",nam1.c_str(),c1,c1);
	else fprintf(fp,"sequence:%s::::::::",nam1.c_str());
	//output
	char c;
	int size=(int)alignment.size();
	for(i=0;i<size;i++)
	{
		if(i%60==0)fprintf(fp,"\n");
		ii=alignment[i].first;
		if(ii<=0)fprintf(fp,"-");
		else
		{
			c=ami1[ii-1];
			if(c=='X'||c=='Z')fprintf(fp,".");
			else fprintf(fp,"%c",c);
		}
	}
	fprintf(fp,"*\n\n\n");
	//mol2 output
	fprintf(fp,">P1;%s\n",nam2.c_str());
	if(PIR_Type==0 || PIR_Type==2)fprintf(fp,"structure:%s::%c::%c::::",nam2.c_str(),c2,c2);
	else fprintf(fp,"sequence:%s::::::::",nam2.c_str());
	//output
	for(i=0;i<size;i++)
	{
		if(i%60==0)fprintf(fp,"\n");
		jj=alignment[i].second;
		if(jj<=0)fprintf(fp,"-");
		else
		{
			c=ami2[jj-1];
			if(c=='X'||c=='Z')fprintf(fp,".");
			else fprintf(fp,"%c",c);
		}
	}
	fprintf(fp,"*\n\n");
}
void PIR_Output_II(FILE *fp,char *seq1,char *seq2,
				   string &nam1,string &nam2,char c1,char c2,int PIR_Type)
{
	//alignment->output
	//mol1 output
	fprintf(fp,">P1;%s\n",nam1.c_str());
	if(PIR_Type==0 || PIR_Type==1)fprintf(fp,"structure:%s::%c::%c::::",nam1.c_str(),c1,c1);
	else fprintf(fp,"sequence:%s::::::::",nam1.c_str());
	//output
	int i;
	char c;
	int size=(int)strlen(seq1);
	for(i=0;i<size;i++)
	{
		if(i%60==0)fprintf(fp,"\n");
		c=seq1[i];
		if(c=='X'||c=='Z')fprintf(fp,".");
		else fprintf(fp,"%c",c);
	}
	fprintf(fp,"*\n\n\n");
	//mol2 output
	fprintf(fp,">P1;%s\n",nam2.c_str());
	if(PIR_Type==0 || PIR_Type==2)fprintf(fp,"structure:%s::%c::%c::::",nam2.c_str(),c2,c2);
	else fprintf(fp,"sequence:%s::::::::",nam2.c_str());
	//output
	for(i=0;i<size;i++)
	{
		if(i%60==0)fprintf(fp,"\n");
		c=seq2[i];
		if(c=='X'||c=='Z')fprintf(fp,".");
		else fprintf(fp,"%c",c);
	}
	fprintf(fp,"*\n\n");
}
/*
//============== DALI_Score Calc =============//
double DALI_Score_Calc(XYZ *mol1,XYZ *mol2,int moln)
{
	int i,j;
	double dist1,dist2,dist_ave;
	double dali_cur_score;
	double dali_tot_score=0.0;
	double Phi_E=0.2;
	double Alpha=20.0;
	for(i=0;i<moln;i++)
	{
		for(j=0;j<moln;j++)
		{
			if(i==j)dali_cur_score=Phi_E;
			else
			{
				dist1=mol1[i].distance(mol1[j]);
				dist2=mol2[i].distance(mol2[j]);
				dist_ave=0.5*(dist1+dist2);
				dali_cur_score=(Phi_E-1.0*fabs((dist1-dist2)/dist_ave))*exp(-1.0*dist_ave*dist_ave/Alpha/Alpha);
			}
			dali_tot_score+=dali_cur_score;
		}
	}
	return dali_tot_score;
}
double DALI_Score_Total(XYZ *mol1,XYZ *mol2,int moln1,int moln2,int *ali2)
{
	XYZ *ww1=new XYZ[moln1];
	XYZ *ww2=new XYZ[moln2];
	int i;
	int ii,jj;
	int count=0;
	for(i=0;i<moln2;i++)
	{
		ii=ali2[i];
		jj=i;
		if(ii<0)continue;
		ww1[count]=mol1[ii];
		ww2[count]=mol2[jj];
		count++;
	}
	double score=DALI_Score_Calc(ww1,ww2,count);
	delete [] ww1;
	delete [] ww2;
	return score;
}
*/

//============= Weight_Matrix_Built ===============//__110815__//
void Linear_Transform(double *in_mat,int moln1,int moln2,double *out_mat)
{
	int i;
	int totnum=moln1*moln2;
	double wscur;
	double wsmax,wsmin;
	wsmax=-99999;
	wsmin=99999;
	for(i=0;i<totnum;i++)
	{
		wscur=in_mat[i];
		if(wscur>wsmax)wsmax=wscur;
		if(wscur<wsmin)wsmin=wscur;
	}
	for(i=0;i<totnum;i++)
	{
		wscur=in_mat[i];
		out_mat[i]=1.0*(wscur-wsmin)/(wsmax-wsmin);
	}
}
void Weight_Matrix_Built_WS(char *ami1,char *cle1,char *ami2,char *cle2,double *wei_out,
							   double ws1,double ws2,double Local_Para)
{
	Bioinfo_Code bioinfo_code(2);
	int len1,len2,w1,w2;
	//get GLO
	len1=(int)strlen(ami1);
	len2=(int)strlen(cle1);
	if(len1!=len2)
	{
		fprintf(stderr,"ami1 != cle1 [%d!=%d]\n",len1,len2);
		exit(-1);
	}
	w1=len1;
	len1=(int)strlen(ami2);
	len2=(int)strlen(cle2);
	if(len1!=len2)
	{
		fprintf(stderr,"ami2 != cle2 [%d!=%d]\n",len1,len2);
		exit(-1);
	}
	w2=len2;

	//AMI only
	bioinfo_code.AMI_transform(ami1,wwa1);
	bioinfo_code.AMI_transform(ami2,wwa2);
	//CLE only
	bioinfo_code.CLE_transform(cle1,wwc1);
	bioinfo_code.CLE_transform(cle2,wwc2);

	//gen WEI
	int i,j,k;
	int range=3;
	int ii,jj;
	int len;
	double wscur;
	double wscur1,wscur2;

	//calc
	for(i=0;i<w1;i++)
	{
		for(j=0;j<w2;j++)
		{
			ii=i;
			jj=j;
			len=1;
			for(k=1;k<=range;k++)
			{
				if(ii-k<0 || jj-k<0)break;
				ii--;
				jj--;
				len++;
			}
			for(k=1;k<=range;k++)
			{
				if(i+k>=w1 || j+k>=w2)break;
				len++;
			}
			wscur1=bioinfo_code.Universal_Calc(ii,jj,len,wwa1,wwa2,w1,w2,3);  //-> BLOSUM
			wscur2=bioinfo_code.Universal_Calc(ii,jj,len,wwc1,wwc2,w1,w2,1);  //-> CLESUM
			//calc_score
			if(wscur1<0)wscur1=0;                   //max(BLOSUM,0)
			wscur=ws1*Local_Para*wscur1+ws2*wscur2; //local parameter
			wei_out[i*w2+j]=wscur;
		}
	}
}

//---------------- superimpose_fullatom --------------------//
void Superimpose_FullAtom(Kabsch &kabsch,PDB_Residue *in,int moln,PDB_Residue *out,double *rotmat)
{
	int i,k;
	int num;
	PDB_Residue pdb;
	XYZ xyz;
	int numb;
	double R,tmpr;
	for(i=0;i<moln;i++)
	{
		pdb=in[i];
		//rotate backbone
		num=pdb.get_backbone_totnum();
		for(k=0;k<num;k++)
		{
			if(pdb.get_backbone_part_index(k)==0)continue;
			pdb.get_backbone_atom(k,xyz,numb,R,tmpr);
			kabsch.rot_point(xyz,xyz,rotmat);
			pdb.set_backbone_atom(k,xyz,numb,R,tmpr);
		}
		//rotate sidechain
		num=pdb.get_sidechain_totnum();
		for(k=0;k<num;k++)
		{
			if(pdb.get_sidechain_part_index(k)==0)continue;
			pdb.get_sidechain_atom(k,xyz,numb,R,tmpr);
			kabsch.rot_point(xyz,xyz,rotmat);
			pdb.set_sidechain_atom(k,xyz,numb,R,tmpr);
		}
		out[i]=pdb;
	}
}

//------------- generate CB --------------//
void Generate_CB(Confo_Beta &confo_beta,int moln,PDB_Residue *pdb,char *ami,char *cle,XYZ *mol,XYZ *mcb)
{
	int k;
	int ws_correct; //default: OK
	PDB_Residue PDB;
	double dist;
	//mol1 process
	ws_correct=1;
	for(k=0;k<moln;k++)
	{
		PDB=pdb[k];
		if(PDB.get_sidechain_part_index(0)==0)
		{
			ws_correct=0;
			ws_reco[k]=0;
		}
		else
		{
			PDB.get_sidechain_atom(0,mcb[k]);
			dist=mcb[k].distance_square(mol[k]);
			if(dist<1.0||dist>4.0)
			{
				ws_correct=0;
				ws_reco[k]=0;
			}
			else ws_reco[k]=1;
		}
	}
	if(ws_correct==0)
	{
		strcpy(tmp_ami,ami);
		for(k=0;k<moln;k++)if(tmp_ami[k]=='G')tmp_ami[k]='A';
		confo_beta.Recon_Beta_21(mol,tmp_mcb,moln,tmp_ami,cle);
		for(k=0;k<moln;k++)if(ws_reco[k]==0)mcb[k]=tmp_mcb[k];
	}
}

//--------- BLOSUM_CLESUM_Output --------//
void BLOSUN_CLESUM_Output(char *ami1,char *ami2,char *cle1,char *cle2,
	int *ali1,int *ali2,int moln1,int moln2,int &seqid,int &blosum,int &clesum,vector <double> &match_wei)
{
	//ali2->ali1
	int i;
	int ii,jj;
	for(i=0;i<moln1;i++)ali1[i]=-1;
	for(i=0;i<moln2;i++)
	{
		if(ali2[i]>=0)
		{
			ii=ali2[i];
			jj=i;
			ali1[ii]=jj;
		}
	}

	//ali1->alignment
	//init
	vector<pair<int,int> > alignment;
	alignment.clear();
	//start
	int j;
	int wlen;
	int pre_ii=0;
	int pre_jj=0;
	for(i=1;i<=moln1;i++)
	{
		ii=i;
		jj=ali1[i-1];  //ali1 starts from 0, correspondence also from 0
		if(jj==-1)
		{
			continue;
		}
		else
		{
			jj++;
			//previous_path
			wlen=ii-pre_ii;
			for(j=1;j<wlen;j++)
			{
				pre_ii++;
				alignment.push_back (pair<int,int>(pre_ii, -pre_jj)); //Ix
			}
			wlen=jj-pre_jj;
			for(j=1;j<wlen;j++)
			{
				pre_jj++;
				alignment.push_back (pair<int,int>(-pre_ii, pre_jj)); //Iy
			}
			//current_path
			alignment.push_back (pair<int,int>(ii, jj)); //Match
			//update
			pre_ii=ii;
			pre_jj=jj;
		}
	}
	//termi
	pre_ii++;
	for(i=pre_ii;i<=moln1;i++)alignment.push_back (pair<int,int>(i, -pre_jj)); //Ix
	pre_jj++;
	for(i=pre_jj;i<=moln2;i++)alignment.push_back (pair<int,int>(-moln1, i));  //Iy

	//middle line
	int size=alignment.size();
	blosum=0;
	clesum=0;
	seqid=0;
	double wscur1,wscur2,wscur;
	match_wei.clear();
	for(i=0;i<size;i++)
	{
		ii=alignment[i].first;
		jj=alignment[i].second;
		if(ii>0 && jj>0)
		{
			//calculate BLOSUM and CLESUM
			wscur1=BLOSUM62_Calc(ami1[ii-1],ami2[jj-1]);
			wscur2=CLESUM_Calc(cle1[ii-1],cle2[jj-1]);
			blosum+=(int)wscur1;
			clesum+=(int)wscur2;
			//calculate seqid
			if(ami1[ii-1]==ami2[jj-1])
			{
				if(ami1[ii-1]!='.' && ami1[ii-1]!='X' && ami1[ii-1]!='Z')seqid++;
			}
			//add to match_wei
			if(wscur1<0)wscur1=0;     //max(BLOSUM,0)
			wscur=wscur1+0.1*wscur2;   //local parameter
			match_wei.push_back(wscur);
		}
	}
	clesum/=10;
}

//-------------- mask region strict check -----------------//
//-> we enforce that all mask regions should be aligned
//[note]: regions should be range from 1 to 9 (0 indicate NULL region)
//        cover_rate is the threshold that the mask region is covered by alignment
int Mask_Region_StrictCheck(CLEFAPS_Main &clepaps,double cover_rate)
{
	//-- init check --//
	if(clepaps.mas1==0 || clepaps.mas2==0)return 1;

	//-- create mask region --//
	vector <int> check_mask1(10,0);
	vector <int> check_mask2(10,0);
	for(int wsi=0;wsi<TM_MOLN1;wsi++)
	{
		int mask=clepaps.mas1[wsi];
		if(mask<0 || mask>9)
		{
			fprintf(stderr,"BAD MASK in mask file 1 \n");
			exit(-1);
		}
		check_mask1[mask]++;
	}
	for(int wsi=0;wsi<TM_MOLN2;wsi++)
	{
		int mask=clepaps.mas2[wsi];
		if(mask<0 || mask>9)
		{
			fprintf(stderr,"BAD MASK in mask file 2 \n");
			exit(-1);
		}
		check_mask2[mask]++;
	}

	//-- check each alignment solution --//
	vector <Align_Record> FM_align_tot;
	FM_align_tot.clear();
	int count=0;
	for(int i=0;i<(int)clepaps.FM_align_tot.size();i++)
	{
		//check mask region
		vector <int> exist_mask1(10,0);
		vector <int> exist_mask2(10,0);
		//get alignment
		for(int wsi=0;wsi<TM_MOLN2;wsi++)TM_ALIGNMENT[wsi]=clepaps.FM_align_tot[i].alignment[wsi];
		//check alignment
		for(int wsi=0;wsi<TM_MOLN2;wsi++)
		{
			if(TM_ALIGNMENT[wsi]>=0)
			{
				int ii=TM_ALIGNMENT[wsi];
				int jj=wsi;
				exist_mask1[clepaps.mas1[ii]]++;
				exist_mask2[clepaps.mas2[jj]]++;
			}
		}
		//check mask
		int correct=1;
		for(int m=1;m<10;m++)
		{
			int exist1=exist_mask1[m];
			int exist2=exist_mask2[m];
			int orig1=check_mask1[m];
			int orig2=check_mask2[m];
			if(orig1==0 || orig2==0)continue;
			int smaller_orig=orig1<orig2?orig1:orig2;
			if(1.0*exist1/smaller_orig<cover_rate || 1.0*exist2/smaller_orig<cover_rate)
			{
				correct=0;
				break;
			}
		}
		if(correct==1)
		{
			FM_align_tot.push_back(clepaps.FM_align_tot[i]);
			count++;
		}
	}
	clepaps.FM_align_tot=FM_align_tot;
	return count;
}


//====================== DeepAlign_main function ===================//__110630__//
int DeepAlign_main(CLEFAPS_Main &clepaps,string &wsnam1,string &wsnam2,string &out,
		string &ali_file,int Out_More,int Out_Screen,int Normalize,int Kill_Gaps,int Score_Func,int QualityLevel,int QualityDegree,
		string &output_root,string &out_str,double &ret_val,int SIMPLY_LOAD,int Kill_Frag,double Mask_StrickCheck,
		double Distance_Cutoff,double Multi_Cut)
{
	//-------- ws_weight_process ------------//__110810__//
	if(Score_Func>=4) //Wei_Score
	{
		double rate1=1.0;
		double rate2=0.1;
		Weight_Matrix_Built_WS(TM_AMI1,TM_CLE1,TM_AMI2,TM_CLE2,wswei,rate1,rate2,Local_Para);
		Linear_Transform(wswei,TM_MOLN1,TM_MOLN2,wswei);
		clepaps.TMali_Weight=wswei;
		if(Score_Func==4)clepaps.TM_Wei_Score=1;
		if(Score_Func==5)clepaps.TM_Vect_Score=0;
		if(Score_Func==6)clepaps.TM_Vect_Score=1;
	}
	else
	{
		if(Score_Func==1)clepaps.TM_Vect_Score=0;
		if(Score_Func==2)clepaps.TM_Vect_Score=1;
	}

	//---- ws_additional_set -------//
	if(Kill_Gaps==0)
	{
		clepaps.FM_Kill_Gap=0;  //don't kill gaps~~~
	}
	else
	{
		clepaps.FM_Kill_Gap=1;  //kill gaps~~~
	}
	if(Kill_Frag==0)
	{
		clepaps.Kill_Frag=0;    //don't kill frag~~~
	}
	else
	{
		clepaps.Kill_Frag=1;    //kill frag~~~
	}

	//------ process ------//
	clepaps.CLEFAPS_Input_Func(TM_CLE1,TM_CLE2,TM_AMI1,TM_AMI2,TM_MOL1,TM_MOL2,TM_MOLN1,TM_MOLN2);
	clepaps.TM_Align_Init(TM_MOLN1,TM_MOLN2);
	//normalize //----------------> need to be rewritten !! //__110630__//
	int norm_len=TM_MOLN2; //normal, using the second to normalize
	double norm_d0=clepaps.d0;
	{
		if(Normalize==-2)norm_len=TM_MOLN2;
		else if(Normalize==-1)norm_len=TM_MOLN1;
		else if(Normalize==0)norm_len=min(TM_MOLN1,TM_MOLN2);
		else if(Normalize>0)norm_len=Normalize;
		else
		{
			fprintf(stderr,"WARNING: incorrect normalize length [%d]!! using min as default !!\n",Normalize);
			norm_len=min(TM_MOLN1,TM_MOLN2);
		}
		norm_d0=clepaps.Calc_TM_d0_Simp(norm_len);
	}

	//---- distance cutoff setup ------//__2016.05.20__//
	if(Distance_Cutoff==0)
	{
		clepaps.TM_DIST_CUT=0;
	}
	else if(Distance_Cutoff<0)
	{
		clepaps.TM_DIST_CUT=1;
		Distance_Cutoff=clepaps.d8;
	}
	else
	{
		clepaps.TM_DIST_CUT=1;
		if(Distance_Cutoff<clepaps.d0+1.0)
		{
			fprintf(stderr,"WARNING: user assigned distance cutoff %lf is smaller than d0+1.0 %lf !!\n",
				Distance_Cutoff,clepaps.d0+1.0);
			clepaps.d8=clepaps.d0+1.0;
		}
		else
		{
			clepaps.d8=Distance_Cutoff;
		}
		clepaps.TM_DistCut=Distance_Cutoff;
	}

	//-------- speed up -------//
	double ws_ret;
	if(ali_file=="") // no initial alignment file input
	{
		if(QualityLevel!=1) //extreme_fast
		{
			clepaps.TM_ADDITION=0;     //no additional TMalign procedure !! //__110720__//
			clepaps.REFINE_ITER=3;     //only three iteration
			clepaps.Refine_Cut=0.5; //50% cutoff to the maximal
			clepaps.ZM_Upper_Max=50*QualityDegree;   //only consider initial 50 SFP_H as TopJ
			clepaps.ZM_Lower_Max=10*QualityDegree;   //only consider initial 10 SFP_H as TopJ
			clepaps.ZM_TopK=10+QualityDegree-1;      //only consider TopK as 10
			if(Out_More>=2)clepaps.ZM_TopL=5;  //only consider top5
			else clepaps.ZM_TopL=2;            //only consider top2
			ws_ret=clepaps.FM_Align_Normal(TM_MOL1,TM_MOL2,TM_MOLN1,TM_MOLN2,TM_ALIGNMENT,norm_len,norm_d0); //fast version
		}
		else //normal version
		{
			ws_ret=clepaps.FM_Align_Normal(TM_MOL1,TM_MOL2,TM_MOLN1,TM_MOLN2,TM_ALIGNMENT,norm_len,norm_d0); //normal version
		}
//		else                     //slow version
//		{
//			ws_ret=clepaps.FM_Align_Total(TM_MOL1,TM_MOL2,TM_MOLN1,TM_MOLN2,TM_ALIGNMENT,norm_len,norm_d0);  //slow version
//		}
	}
	else             // initial alignment file input
	{
		//load alignment file
		vector<pair<int, int> > alignment;
		string nam1_content,nam2_content;
		string nam1_full,nam2_full;
		string wnam1,wnam2;
		int retv=ReadToFile_FASTA(ali_file,alignment,nam1_content,nam2_content,nam1_full,nam1_full,wnam1,wnam2);
		if(retv==-1)exit(-1);
		//extract alignment
		string nam1_pdb=TM_AMI1;
		string nam2_pdb=TM_AMI2;
		vector<pair<int, int> > alignment_out;
		Extract_Alignment(nam1_content,nam1_pdb,nam2_content,nam2_pdb,alignment,alignment_out);
		//process alignment
		for(int wwi=0;wwi<TM_MOLN2;wwi++)TM_ALIGNMENT[wwi]=-1;
		int ww_lali=0;
		for(int wwi=0;wwi<(int)alignment_out.size();wwi++)
		{
			int pos1=alignment_out[wwi].first;
			int pos2=alignment_out[wwi].second;
			if(pos1>0 && pos2>0)
			{
				TM_ALIGNMENT[pos2-1]=pos1-1;
				ww_lali++;
			}
		}
		//refine alignment
		clepaps.CLEPAPS_Input(TM_MOL1,TM_MOL2,TM_MOLN1,TM_MOLN2);
		ws_ret=clepaps.FM_Align_WithAli(TM_MOL1,TM_MOL2,TM_MOLN1,TM_MOLN2,TM_ALIGNMENT,norm_len,norm_d0);  //align version
	}
	
	//--- we do realign process ----//
	if(ws_ret==0 && QualityLevel>0)
	{
		ws_ret=clepaps.FM_Align_Total(TM_MOL1,TM_MOL2,TM_MOLN1,TM_MOLN2,TM_ALIGNMENT,norm_len,norm_d0);  //slow version
	}
	if(ws_ret==0)
	{
		char ws_command[300000];
		double zerod=0;
		int zeroi=0;
		sprintf(ws_command,"%s %s %4d %4d -> %6d %6d %9.2f -> %4d %7.3f %7.3f -> %6.3f %6.3f %6.3f -> %5d %5d %6.3f %.3f\n",
			wsnam1.c_str(),wsnam2.c_str(),TM_MOLN1,TM_MOLN2,zeroi,zeroi,zerod,zeroi,zerod,zerod,zerod,zerod,zerod,zeroi,norm_len,Distance_Cutoff,zerod);
		out_str=ws_command;
		ret_val=zerod;
		return 0;
	}

	//--- strict check mask region alignment ----//__180820__//
	if(Mask_StrickCheck>0)
	{
		int mask_ret=Mask_Region_StrictCheck(clepaps,Mask_StrickCheck);
		if(mask_ret==0)
		{
			char ws_command[300000];
			double zerod=0;
			int zeroi=0;
			sprintf(ws_command,"%s %s %4d %4d -> %6d %6d %9.2f -> %4d %7.3f %7.3f -> %6.3f %6.3f %6.3f -> %5d %5d %6.3f %.3f\n",
			wsnam1.c_str(),wsnam2.c_str(),TM_MOLN1,TM_MOLN2,zeroi,zeroi,zerod,zeroi,zerod,zerod,zerod,zerod,zerod,zeroi,norm_len,Distance_Cutoff,zerod);
			out_str=ws_command;
			ret_val=zerod;
			return 0;
		}
	}

	//------ output ------//
	FILE *fp;
	string name,file;
	name=out;
	double Ret_Sco[8];
	int fin_lali=0;
	int fin_blosum=0;
	int fin_clesum=0;
	int fin_seqid=0;
	double fin_tms,fin_rms,fin_wms,fin_gdt,fin_gdtha,fin_maxsub,fin_ugdt;
	fin_tms=fin_rms=fin_wms=fin_gdt=fin_gdtha=fin_maxsub=fin_ugdt=0;
	string final_record="";
	if(Out_More>=0)
	{
		//init
		int wsk,wsi;
		int wssize=(int)clepaps.FM_align_tot.size();
		int relout=wssize<10?wssize:10;
		char www_nam[300000];
		Kabsch kabsch;
		string TER="TER                                                                             ";
		string END="END                                                                             ";
		int lali;
		int seqid;
		int blosum,clesum;
		double tms,wms,rms;
		double gdt1,gdt2,maxsub,ugdt;
		vector <double> match_wei;
		vector <double> rmsd_out;
		//--------- single solution -----------//
		wsk=0;
		for(wsi=0;wsi<TM_MOLN2;wsi++)TM_ALIGNMENT[wsi]=clepaps.FM_align_tot[0].alignment[wsi];
		for(wsi=0;wsi<12;wsi++)TM_ROTMAT[wsi]=clepaps.FM_align_tot[0].rotmat[wsi];
		if(Out_More>=0)
		{
			double d0=clepaps.d0;
			string ws_output_tot;
			kabsch.rot_mol(TM_MOL1,TM_MOL_TMP,TM_MOLN1,TM_ROTMAT);
			FASTA_Output_More(ws_output_tot,wsnam1,wsnam2,TM_AMI1,TM_AMI2,TM_CLE1,TM_CLE2,TM_MOL_TMP,TM_MOL2,TM_ALI1,
				TM_ALIGNMENT,TM_MOLN1,TM_MOLN2,d0,seqid,blosum,clesum,match_wei,rmsd_out,0);
			final_record=ws_output_tot;
		}
		//---- for file output ----//
		if(Out_More>=1)
		{
			//[1]fasta_out
			if(Out_More>=2)sprintf(www_nam,"%s/%s.%c.fasta",output_root.c_str(),name.c_str(),wsk+'0');
			else sprintf(www_nam,"%s/%s.fasta",output_root.c_str(),name.c_str());
			fp=fopen(www_nam,"wb");
			if(fp==0)
			{
				fprintf(stderr,"ERROR: file %s can't be opened. \n",www_nam);
			}
			else
			{
				FASTA_Output(fp,wsnam1,wsnam2,TM_AMI1,TM_AMI2,TM_ALI1,TM_ALIGNMENT,TM_MOLN1,TM_MOLN2);
				fclose(fp);
			}
			//[1-1]cle_fasta_out
//			if(Out_More>=2)sprintf(www_nam,"%s/%s.%c.fasta_cle",output_root.c_str(),name.c_str(),wsk+'0');
//			else sprintf(www_nam,"%s/%s.fasta_cle",output_root.c_str(),name.c_str());
//			fp=fopen(www_nam,"wb");
//			FASTA_Output(fp,wsnam1,wsnam2,TM_CLE1,TM_CLE2,TM_ALI1,TM_ALIGNMENT,TM_MOLN1,TM_MOLN2);
//			fclose(fp);
			//[2]pdb_out
			kabsch.rot_mol(TM_MOL1,TM_MOL_TMP,TM_MOLN1,TM_ROTMAT);
			if(Out_More>=2)sprintf(www_nam,"%s/%s.%c.pdb",output_root.c_str(),name.c_str(),wsk+'0');
			else sprintf(www_nam,"%s/%s.pdb",output_root.c_str(),name.c_str());
			fp=fopen(www_nam,"wb");
			if(fp==0)
			{
				fprintf(stderr,"ERROR: file %s can't be opened. \n",www_nam);
			}
			else
			{
				if(SIMPLY_LOAD==0)  //-> output Full_Atom PDB
				{
					Superimpose_FullAtom(kabsch,TM_PDB1,TM_MOLN1,TM_PDB_TMP,TM_ROTMAT);
					mol_out.Output_PDB_III(fp,TM_MOLN1,TM_PDB_TMP,'A',1);
					fprintf(fp,"%s\n",TER.c_str());
					mol_out.Output_PDB_III(fp,TM_MOLN2,TM_PDB2,'B',1);
					fprintf(fp,"%s\n",TER.c_str());
					fprintf(fp,"%s\n",END.c_str());
					//--> output for linear chain 
					if(Out_PDB==1)
					{
						// output file
						if(Out_More>=2)sprintf(www_nam,"%s/%s.%c.pdb_linear",output_root.c_str(),name.c_str(),wsk+'0');
						else sprintf(www_nam,"%s/%s.pdb_linear",output_root.c_str(),name.c_str());
						FILE *fq=fopen(www_nam,"wb");
						if(fq==0)
						{
							fprintf(stderr,"ERROR: file %s can't be opened. \n",www_nam);
						}
						// output
						mol_out.Output_PDB_III(fq,TM_MOLN1,TM_PDB_TMP,'A',-1);
						fprintf(fq,"%s\n",TER.c_str());
						mol_out.Output_PDB_III(fq,TM_MOLN2,TM_PDB2,'B',-1);
						fprintf(fq,"%s\n",TER.c_str());
						fprintf(fq,"%s\n",END.c_str());
						// close
						fclose(fq);
					}
				}
				else                //-> output CA_Only PDB
				{
					mol_out.Output_PDB(fp,TM_MOLN1,TM_MOL_TMP,TM_AMI1,0,'A');
					fprintf(fp,"%s\n",TER.c_str());
					mol_out.Output_PDB(fp,TM_MOLN2,TM_MOL2,TM_AMI2,0,'B');
					fprintf(fp,"%s\n",TER.c_str());
					fprintf(fp,"%s\n",END.c_str());
				}
				fclose(fp);
			}
			//[3]local_out
			if(Out_More>=2)sprintf(www_nam,"%s/%s.%c.local",output_root.c_str(),name.c_str(),wsk+'0');
			else sprintf(www_nam,"%s/%s.local",output_root.c_str(),name.c_str());
			fp=fopen(www_nam,"wb");
			if(fp==0)
			{
				fprintf(stderr,"ERROR: file %s can't be opened. \n",www_nam);
			}
			else
			{
				double d0=clepaps.d0;
				string ws_output_tot;
				FASTA_Output_More(ws_output_tot,wsnam1,wsnam2,TM_AMI1,TM_AMI2,TM_CLE1,TM_CLE2,TM_MOL_TMP,TM_MOL2,TM_ALI1,
					TM_ALIGNMENT,TM_MOLN1,TM_MOLN2,d0,seqid,blosum,clesum,match_wei,rmsd_out,1);
				fprintf(fp,"%s",ws_output_tot.c_str());
				fclose(fp);
			}
			//[4]rmsd_out
			//if(Out_More>=2)sprintf(www_nam,"%s/%s.%c.rmsd",output_root.c_str(),name.c_str(),wsk+'0');
			//else sprintf(www_nam,"%s/%s.rmsd",output_root.c_str(),name.c_str());
			//fp=fopen(www_nam,"wb");
			//fprintf(fp,">Alignment_RMSD\n");
			//for(int kk=0;kk<(int)rmsd_out.size();kk++)
			//{
			//	char wstmp;
			//	if(rmsd_out[kk]<0)wstmp='-';
			//	else if(rmsd_out[kk]>9)wstmp='9';
			//	else
			//	{
			//		int wssco=(int)rmsd_out[kk];
			//		wstmp=wssco+'0';
			//	}
			//	fprintf(fp,"%c",wstmp);
			//}
			//fprintf(fp,"\n");
			//fclose(fp);
		}
		else  //-> no file output 
		{
			BLOSUN_CLESUM_Output(TM_AMI1,TM_AMI2,TM_CLE1,TM_CLE2,TM_ALI1,TM_ALIGNMENT,TM_MOLN1,TM_MOLN2,seqid,blosum,clesum,match_wei);
		}
		//[4]score_out
		tms=clepaps.TM_Align_TM_Score(TM_MOL1,TM_MOL2,TM_MOLN1,TM_MOLN2,TM_ALIGNMENT,norm_len,norm_d0,rms,lali,Ret_Sco);  //tm-score
		wms=clepaps.TM_Align_Get_Score_Simp_MatchWei(TM_MOL1,TM_MOL2,TM_ROTMAT,TM_MOLN1,TM_MOLN2,TM_ALIGNMENT,match_wei); //DeepAlign-score
		gdt1=1.0*(Ret_Sco[1]+Ret_Sco[2]+Ret_Sco[3]+Ret_Sco[4])/(4.0*norm_len);  //ori_GDT
		gdt2=1.0*(Ret_Sco[0]+Ret_Sco[1]+Ret_Sco[2]+Ret_Sco[3])/(4.0*norm_len);  //ha_GDT
		ugdt=0.25*(Ret_Sco[1]+Ret_Sco[2]+Ret_Sco[3]+Ret_Sco[4]);                //uGDT
		maxsub=1.0*Ret_Sco[5]/norm_len;                                         //maxsub
		if(Out_More>=1)
		{
			if(Out_More>=2)sprintf(www_nam,"%s/%s.%c.score",output_root.c_str(),name.c_str(),wsk+'0');
			else sprintf(www_nam,"%s/%s.score",output_root.c_str(),name.c_str());
			fp=fopen(www_nam,"wb");
			if(fp==0)
			{
				fprintf(stderr,"ERROR: file %s can't be opened. \n",www_nam);
			}
			else
			{
				fprintf(fp,"#name1 name2 len1 len2 -> BLOSUM CLESUM DeepScore -> LALI RMSDval TMscore -> MAXSUB GDT_TS GDT_HA -> SeqID nLen dCutoff uGDT\n");
				fprintf(fp," %s %s %4d %4d -> %6d %6d %9.2f -> %4d %7.3f %7.3f -> %6.3f %6.3f %6.3f -> %5d %5d %6.3f %.3f\n",
					wsnam1.c_str(),wsnam2.c_str(),TM_MOLN1,TM_MOLN2,blosum,clesum,wms,lali,rms,tms,maxsub,gdt1,gdt2,seqid,norm_len,Distance_Cutoff,ugdt);
				fprintf(fp,"#---------------- transformation to superpose 1st structure onto the 2nd ------------------------\n");
				fprintf(fp,"%9.6f %9.6f %9.6f %12.6f\n",TM_ROTMAT[0],TM_ROTMAT[1],TM_ROTMAT[2],TM_ROTMAT[9]);
				fprintf(fp,"%9.6f %9.6f %9.6f %12.6f\n",TM_ROTMAT[3],TM_ROTMAT[4],TM_ROTMAT[5],TM_ROTMAT[10]);
				fprintf(fp,"%9.6f %9.6f %9.6f %12.6f\n",TM_ROTMAT[6],TM_ROTMAT[7],TM_ROTMAT[8],TM_ROTMAT[11]);
				fclose(fp);
			}
		}
		//record
		fin_tms=tms;
		fin_rms=rms;
		fin_wms=wms;
		fin_gdt=gdt1;
		fin_gdtha=gdt2;
		fin_ugdt=ugdt;
		fin_maxsub=maxsub;
		fin_lali=lali;
		fin_blosum=blosum;
		fin_clesum=clesum;
		fin_seqid=seqid;
		//--------- multiple solution -----------//
		if(Out_More>=2)
		{
			double *temp_score=new double[relout];
			int *temp_index=new int[relout];
			//----- get all scores ----//
			temp_score[0]=-1;
			for(int wsk_=1;wsk_<relout;wsk_++)
			{
				temp_score[wsk_]=-1;
				int pos=wsk_;
				for(wsi=0;wsi<TM_MOLN2;wsi++)TM_ALIGNMENT[wsi]=clepaps.FM_align_tot[pos].alignment[wsi];
				for(wsi=0;wsi<12;wsi++)TM_ROTMAT[wsi]=clepaps.FM_align_tot[pos].rotmat[wsi];
				//score quality check
				tms=clepaps.TM_Align_TM_Score(TM_MOL1,TM_MOL2,TM_MOLN1,TM_MOLN2,TM_ALIGNMENT,norm_len,norm_d0,rms,lali,Ret_Sco);  //tm-score
				gdt1=1.0*(Ret_Sco[1]+Ret_Sco[2]+Ret_Sco[3]+Ret_Sco[4])/(4.0*norm_len);  //ori_GDT
				gdt2=1.0*(Ret_Sco[0]+Ret_Sco[1]+Ret_Sco[2]+Ret_Sco[3])/(4.0*norm_len);  //ha_GDT
				ugdt=0.25*(Ret_Sco[1]+Ret_Sco[2]+Ret_Sco[3]+Ret_Sco[4]);                //uGDT
				maxsub=1.0*Ret_Sco[5]/norm_len;
//				if(tms < 0.8*fin_tms || gdt1 < 0.8*fin_gdt)  // we use 0.8 cutoff here !!
//				{
//					temp_score[wsk_]=-1;
//					continue;
//				}
//				if(tms < 0.35)  //-> we show multi solutions for TMscore > 0.35 !!
				if(tms < Multi_Cut )
				{
					temp_score[wsk_]=-1;
					continue;
				}
				//record
				temp_score[wsk_]=gdt1;
			}
			Fast_Sort <double> fast_sort;
			fast_sort.fast_sort_1(temp_score,temp_index,relout);

			//----- all the others ----//
			wsk=1;
			for(int wsk_=0;wsk_<relout;wsk_++)
			{
				int pos=temp_index[wsk_];
				if(temp_score[pos]<0)break;
				for(wsi=0;wsi<TM_MOLN2;wsi++)TM_ALIGNMENT[wsi]=clepaps.FM_align_tot[pos].alignment[wsi];
				for(wsi=0;wsi<12;wsi++)TM_ROTMAT[wsi]=clepaps.FM_align_tot[pos].rotmat[wsi];
				//[1]fasta_out
				sprintf(www_nam,"%s/%s.%c.fasta",output_root.c_str(),name.c_str(),wsk+'0');
				fp=fopen(www_nam,"wb");
				if(fp==0)
				{
					fprintf(stderr,"ERROR: file %s can't be opened. \n",www_nam);
				}
				else
				{
					FASTA_Output(fp,wsnam1,wsnam2,TM_AMI1,TM_AMI2,TM_ALI1,TM_ALIGNMENT,TM_MOLN1,TM_MOLN2);
					fclose(fp);
				}
				//[1-1]cle_fasta_out
//				sprintf(www_nam,"%s/%s.%c.fasta_cle",output_root.c_str(),name.c_str(),wsk+'0');
//				fp=fopen(www_nam,"wb");
//				FASTA_Output(fp,wsnam1,wsnam2,TM_CLE1,TM_CLE2,TM_ALI1,TM_ALIGNMENT,TM_MOLN1,TM_MOLN2);
//				fclose(fp);
				//[2]pdb_out
				kabsch.rot_mol(TM_MOL1,TM_MOL_TMP,TM_MOLN1,TM_ROTMAT);
				sprintf(www_nam,"%s/%s.%c.pdb",output_root.c_str(),name.c_str(),wsk+'0');
				fp=fopen(www_nam,"wb");
				if(fp==0)
				{
					fprintf(stderr,"ERROR: file %s can't be opened. \n",www_nam);
				}
				else
				{
					if(SIMPLY_LOAD==0)  //-> output Full_Atom PDB
					{
						Superimpose_FullAtom(kabsch,TM_PDB1,TM_MOLN1,TM_PDB_TMP,TM_ROTMAT);
						mol_out.Output_PDB_III(fp,TM_MOLN1,TM_PDB_TMP,'A',1);
						fprintf(fp,"%s\n",TER.c_str());
						mol_out.Output_PDB_III(fp,TM_MOLN2,TM_PDB2,'B',1);
						fprintf(fp,"%s\n",TER.c_str());
						fprintf(fp,"%s\n",END.c_str());
						//--> output for linear chain
						if(Out_PDB==1)
						{
							// output file
							if(Out_More>=2)sprintf(www_nam,"%s/%s.%c.pdb_linear",output_root.c_str(),name.c_str(),wsk+'0');
							else sprintf(www_nam,"%s/%s.pdb_linear",output_root.c_str(),name.c_str());
							FILE *fq=fopen(www_nam,"wb");
							if(fq==0)
							{
								fprintf(stderr,"ERROR: file %s can't be opened. \n",www_nam);
							}
							// output
							mol_out.Output_PDB_III(fq,TM_MOLN1,TM_PDB_TMP,'A',-1);
							fprintf(fq,"%s\n",TER.c_str());
							mol_out.Output_PDB_III(fq,TM_MOLN2,TM_PDB2,'B',-1);
							fprintf(fq,"%s\n",TER.c_str());
							fprintf(fq,"%s\n",END.c_str());
							// close
							fclose(fq);
						}
					}
					else                //-> output CA_Only PDB
					{
						mol_out.Output_PDB(fp,TM_MOLN1,TM_MOL_TMP,TM_AMI1,0,'A');
						fprintf(fp,"%s\n",TER.c_str());
						mol_out.Output_PDB(fp,TM_MOLN2,TM_MOL2,TM_AMI2,0,'B');
						fprintf(fp,"%s\n",TER.c_str());
						fprintf(fp,"%s\n",END.c_str());
					}
					fclose(fp);
				}
				//[3]local_out
				sprintf(www_nam,"%s/%s.%c.local",output_root.c_str(),name.c_str(),wsk+'0');
				fp=fopen(www_nam,"wb");
				if(fp==0)
				{
					fprintf(stderr,"ERROR: file %s can't be opened. \n",www_nam);
				}
				else
				{
					double d0=clepaps.d0;
					string ws_output_tot;
					FASTA_Output_More(ws_output_tot,wsnam1,wsnam2,TM_AMI1,TM_AMI2,TM_CLE1,TM_CLE2,TM_MOL_TMP,TM_MOL2,TM_ALI1,
						TM_ALIGNMENT,TM_MOLN1,TM_MOLN2,d0,seqid,blosum,clesum,match_wei,rmsd_out,1);
					fprintf(fp,"%s",ws_output_tot.c_str());
					fclose(fp);
				}
				//[4]rmsd_out
				//sprintf(www_nam,"%s/%s.%c.rmsd",output_root.c_str(),name.c_str(),wsk+'0');
				//fp=fopen(www_nam,"wb");
				//fprintf(fp,">Alignment_RMSD\n");
				//for(int kk=0;kk<(int)rmsd_out.size();kk++)
				//{
				//	char wstmp;
				//	if(rmsd_out[kk]<0)wstmp='-';
				//	else if(rmsd_out[kk]>9)wstmp='9';
				//	else
				//	{
				//		int wssco=(int)rmsd_out[kk];
				//		wstmp=wssco+'0';
				//	}
				//	fprintf(fp,"%c",wstmp);
				//}
				//fprintf(fp,"\n");
				//fclose(fp);
				//[5]score_out
				kabsch.rot_mol(TM_MOL1,TM_MOL_TMP,TM_MOLN1,TM_ROTMAT);
				tms=clepaps.TM_Align_TM_Score(TM_MOL1,TM_MOL2,TM_MOLN1,TM_MOLN2,TM_ALIGNMENT,norm_len,norm_d0,rms,lali,Ret_Sco);  //tm-score
				wms=clepaps.TM_Align_Get_Score_Simp_MatchWei(TM_MOL1,TM_MOL2,TM_ROTMAT,TM_MOLN1,TM_MOLN2,TM_ALIGNMENT,match_wei); //DeepAlign-score
				gdt1=1.0*(Ret_Sco[1]+Ret_Sco[2]+Ret_Sco[3]+Ret_Sco[4])/(4.0*norm_len);  //ori_GDT
				gdt2=1.0*(Ret_Sco[0]+Ret_Sco[1]+Ret_Sco[2]+Ret_Sco[3])/(4.0*norm_len);  //ha_GDT
				ugdt=0.25*(Ret_Sco[1]+Ret_Sco[2]+Ret_Sco[3]+Ret_Sco[4]);                //uGDT
				maxsub=1.0*Ret_Sco[5]/norm_len;
				sprintf(www_nam,"%s/%s.%c.score",output_root.c_str(),name.c_str(),wsk+'0');
				fp=fopen(www_nam,"wb");
				if(fp==0)
				{
					fprintf(stderr,"ERROR: file %s can't be opened. \n",www_nam);
				}
				else
				{
					fprintf(fp,"#name1 name2 len1 len2 -> BLOSUM CLESUM DeepScore -> LALI RMSDval TMscore -> MAXSUB GDT_TS GDT_HA -> SeqID nLen dCutoff uGDT\n");
					fprintf(fp," %s %s %4d %4d -> %6d %6d %9.2f -> %4d %7.3f %7.3f -> %6.3f %6.3f %6.3f -> %5d %5d %6.3f %.3f\n",
						wsnam1.c_str(),wsnam2.c_str(),TM_MOLN1,TM_MOLN2,blosum,clesum,wms,lali,rms,tms,maxsub,gdt1,gdt2,seqid,norm_len,Distance_Cutoff,ugdt);
					fprintf(fp,"#---------------- transformation to superpose 1st structure onto the 2nd ------------------------\n");
					fprintf(fp,"%9.6f %9.6f %9.6f %12.6f\n",TM_ROTMAT[0],TM_ROTMAT[1],TM_ROTMAT[2],TM_ROTMAT[9]);
					fprintf(fp,"%9.6f %9.6f %9.6f %12.6f\n",TM_ROTMAT[3],TM_ROTMAT[4],TM_ROTMAT[5],TM_ROTMAT[10]);
					fprintf(fp,"%9.6f %9.6f %9.6f %12.6f\n",TM_ROTMAT[6],TM_ROTMAT[7],TM_ROTMAT[8],TM_ROTMAT[11]);
					fclose(fp);
				}
				//increase
				wsk++;
			}
			//delete
			delete [] temp_score;
			delete [] temp_index;
		}
	}

	//--- print out ----//
	char ws_command[300000];
	if(Out_Screen==0)  //simplest screen-out
	{
		sprintf(ws_command,"%s %s %4d %4d -> %6d %6d %9.2f -> %4d %7.3f %7.3f -> %6.3f %6.3f %6.3f -> %5d %5d %6.3f %.3f\n",
			wsnam1.c_str(),wsnam2.c_str(),TM_MOLN1,TM_MOLN2,fin_blosum,fin_clesum,fin_wms,fin_lali,fin_rms,
			fin_tms,fin_maxsub,fin_gdt,fin_gdtha,fin_seqid,norm_len,Distance_Cutoff,fin_ugdt);
		out_str=ws_command;
	}
	else //-> screen out new
	{
		string fin_tmp_str="";
		//--- header out ---//
		sprintf(ws_command,"#------------------------------------------------------#\n");
		fin_tmp_str+=ws_command;
		sprintf(ws_command,"#                    DeepAlign                         #\n");
		fin_tmp_str+=ws_command;
		sprintf(ws_command,"# Protein structure alignment beyond spatial proximity #\n");
		fin_tmp_str+=ws_command;
		sprintf(ws_command,"# Reference:     S. Wang, J. Ma, J. Peng and J. Xu.    #\n");
		fin_tmp_str+=ws_command;
		sprintf(ws_command,"#                Scientific Reports, 3, 1448, (2013)   #\n");
		fin_tmp_str+=ws_command;
		sprintf(ws_command,"#------------------------------------------------------#\n");
		fin_tmp_str+=ws_command;
		sprintf(ws_command,"\n");
		fin_tmp_str+=ws_command;
		sprintf(ws_command,"1st input protein: %s    length= %4d \n",wsnam1.c_str(),TM_MOLN1);
		fin_tmp_str+=ws_command;
		sprintf(ws_command,"2nd input protein: %s    length= %4d \n",wsnam2.c_str(),TM_MOLN2);
		fin_tmp_str+=ws_command;
		sprintf(ws_command,"\n");
		fin_tmp_str+=ws_command;
		sprintf(ws_command,"#----- Normalization length for TMscore, MAXSUB, GDT_TS, and GDT_HA is %4d\n",norm_len);
		fin_tmp_str+=ws_command;
		sprintf(ws_command,"#----- Distance cutoff value is %lf\n",Distance_Cutoff);
		fin_tmp_str+=ws_command;
		sprintf(ws_command,"# BLOSUM CLESUM DeepScore SeqID LALI RMSD(A) TMscore MAXSUB GDT_TS GDT_HA uGDT\n");
		fin_tmp_str+=ws_command;
		sprintf(ws_command,"  %6d %6d %9.2f %5d %4d %7.3f %7.3f %6.3f %6.3f %6.3f %.3f\n",
			fin_blosum,fin_clesum,fin_wms,fin_seqid,fin_lali,fin_rms,fin_tms,fin_maxsub,fin_gdt,fin_gdtha,fin_ugdt);
		fin_tmp_str+=ws_command;
		sprintf(ws_command,"#----- Please see http://raptorx.uchicago.edu/DeepAlign/documentation/ for explanation of these scores \n\n");
		fin_tmp_str+=ws_command;
		sprintf(ws_command,"#----- Transformation to superpose the 1st protein onto the 2nd ---#\n");
		fin_tmp_str+=ws_command;
		sprintf(ws_command,"# i          t(i)         u(i,1)         u(i,2)         u(i,3) \n");
		fin_tmp_str+=ws_command;
		sprintf(ws_command,"  1      %12.6f   %12.9f   %12.9f   %12.9f\n",TM_ROTMAT[9], TM_ROTMAT[0],TM_ROTMAT[1],TM_ROTMAT[2]);
		fin_tmp_str+=ws_command;
		sprintf(ws_command,"  2      %12.6f   %12.9f   %12.9f   %12.9f\n",TM_ROTMAT[10],TM_ROTMAT[3],TM_ROTMAT[4],TM_ROTMAT[5]);
		fin_tmp_str+=ws_command;
		sprintf(ws_command,"  3      %12.6f   %12.9f   %12.9f   %12.9f\n",TM_ROTMAT[11],TM_ROTMAT[6],TM_ROTMAT[7],TM_ROTMAT[8]);
		fin_tmp_str+=ws_command;
		sprintf(ws_command,"\n");
		fin_tmp_str+=ws_command;
		//--- local out ---//
		sprintf(ws_command,"%s",final_record.c_str());
		fin_tmp_str+=ws_command;
		//-- final return --//
		out_str=fin_tmp_str;
	}
	ret_val=fin_wms;
	return 1;
}

//-- screen-out new ---// according to TMalign
/*
 **************************************************************************
 *                               TM-align                                 *
 * A protein structural alignment algorithm based on TM-score             *
 * Reference: Y. Zhang and J. Skolnick, Nucl. Acids Res. 2005 33, 2302-9  *
 * Comments on the program, please email to: yzhang@ku.edu                *
 **************************************************************************

Chain 1:../../DATA  Size= 120
Chain 2:../../DATA  Size= 120 (TM-score is normalized by  120)

Aligned length= 120, RMSD=  0.00, TM-score=1.00000, ID=1.000

 -------- rotation matrix to rotate Chain-1 to Chain-2 ------
 i          t(i)         u(i,1)         u(i,2)         u(i,3)
 1      0.0000000000   1.0000000000   0.0000000000   0.0000000000
 2      0.0000000000   0.0000000000   1.0000000000   0.0000000000
 3      0.0000000000   0.0000000000   0.0000000000   1.0000000000

(":" denotes the residue pairs of distance < 5.0 Angstrom)
ENIEVHMLNKGAEGAMVFEPAYIKANPGDTVTFIPVDKGHNVESIKDMIPEGAEKFKSKINENYVLTVTQPGAYLVKCTPHYAMGMIALIAVGDSPANLDQIVSAKKPKIVQERLEKVIA
111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111
ENIEVHMLNKGAEGAMVFEPAYIKANPGDTVTFIPVDKGHNVESIKDMIPEGAEKFKSKINENYVLTVTQPGAYLVKCTPHYAMGMIALIAVGDSPANLDQIVSAKKPKIVQERLEKVIA
*/


//====================== DeepAlign_search ===================//__121130__//
int DeepAlign_search(CLEFAPS_Main &clepaps,string &wsnam1,string &wsnam2,string &output,double &ret_val,double tms_thres)
{
	ret_val=-1.0;
	//process
	clepaps.CLEFAPS_Input_Func(TM_CLE1,TM_CLE2,TM_AMI1,TM_AMI2,TM_MOL1,TM_MOL2,TM_MOLN1,TM_MOLN2);
	clepaps.TM_Align_Init(TM_MOLN1,TM_MOLN2);
	//normalize //----------------> need to be rewritten !! //__110630__//
	int norm_len=TM_MOLN2; //normal, using the second to normalize
	double norm_d0=clepaps.d0;
	int Normalize=0;
	{
		if(Normalize==-2)norm_len=TM_MOLN2;
		else if(Normalize==-1)norm_len=TM_MOLN1;
		else if(Normalize==0)norm_len=min(TM_MOLN1,TM_MOLN2);
		else if(Normalize>0)norm_len=Normalize;
		else
		{
			fprintf(stderr,"WARNING: incorrect normalize length [%d]!! using min as default !!\n",Normalize);
			norm_len=min(TM_MOLN1,TM_MOLN2);
		}
		norm_d0=clepaps.Calc_TM_d0_Simp(norm_len);
	}

	//-------- record original parameters ------//__181204__//
	int FAST_Chk_=clepaps.FAST_Chk;
	int REF_Strategy_=clepaps.REF_Strategy;
	double CUR_MaxJ_thres_=clepaps.CUR_MaxJ_thres;
	//------------------------------------------//__181204__//

	//-------- speed up -------//
	double ws_ret;
	int ws_lali;
	clepaps.FAST_Chk=1;
	clepaps.REF_Strategy=0;
	tms_thres=tms_thres>0.17?tms_thres:0.17;
	tms_thres=tms_thres<0.25?tms_thres:0.25;
	clepaps.CUR_MaxJ_thres=tms_thres;    //tms_thres
//	if(QualityLevel==0) //extreme_fast
	{
		clepaps.TM_ADDITION=0;     //no additional TMalign procedure !! //__110720__//
		clepaps.REFINE_ITER=3;     //only three iteration
		clepaps.Refine_Cut=0.5; //50% cutoff to the maximal
		clepaps.ZM_Upper_Max=50;   //only consider initial 50 SFP_H as TopJ upper bound
		clepaps.ZM_Lower_Max=10;   //only consider initial 10 SFP_H as TopJ lower bound
		clepaps.ZM_TopK=10;        //only consider TopK as 10
		clepaps.ZM_TopL=2;         //only consider top2
		ws_ret=clepaps.FM_Align_Lite(TM_MOL1,TM_MOL2,TM_MOLN1,TM_MOLN2,TM_ALIGNMENT,norm_len,norm_d0); //fast version
	}

	//-------- return original parameters ------//__181204__//
	clepaps.FAST_Chk=FAST_Chk_;
	clepaps.REF_Strategy=REF_Strategy_;
	clepaps.CUR_MaxJ_thres=CUR_MaxJ_thres_;
	//------------------------------------------//__181204__//

	//------- result analysis -------//
	if(ws_ret<=0)
	{
		return 0;
	}
	ws_lali=clepaps.FM_align_tot[0].lali;
	char ws_command[300000];
	sprintf(ws_command,"%s %s %4d %4d -> %5.3f %4d -> %lf \n",wsnam1.c_str(),wsnam2.c_str(),TM_MOLN1,TM_MOLN2,ws_ret,ws_lali,ws_ret*ws_lali);
	output=ws_command;
	ret_val=ws_ret;

	//----- return -----//
	return 1;
}

//========================= list process =========================//
//----- load pair_list -----//__130630__//
//-> in format <temp targ> [range1;range2;]
int Get_PairList(string list,vector <string> &temp_rec,vector <string> &targ_rec,
	vector <string> &temp_rec_addi,vector <string> &targ_rec_addi)
{
	temp_rec.clear();
	targ_rec.clear();
	temp_rec_addi.clear();
	targ_rec_addi.clear();
	ifstream fin;
	string buf,temp1,temp2,rang1,rang2;
	fin.open(list.c_str(), ios::in);
	if(fin.fail()!=0)return 0;
	int count=0;
	string def_range="_";
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		istringstream www(buf);
		//-> read pair
		if(! (www>>temp1))continue;
		if(! (www>>temp2))continue;
		temp_rec.push_back(temp1);
		targ_rec.push_back(temp2);
		temp_rec_addi.push_back(def_range);
		targ_rec_addi.push_back(def_range);
		count++;
		//-> read range
		rang1="";
		rang2="";
		if(!getline(www,rang1,';'))continue;
		if(!getline(www,rang2,';'))continue;
		if(rang1!="")temp_rec_addi[count-1]=rang1;
		if(rang2!="")targ_rec_addi[count-1]=rang2;
	}
	return count;
}

//----- load single_list -----//__130630__//
//-> in format <temp> [range;]
int Get_TempList(string list,vector <string> &temp_rec,vector <string> &temp_rec_addi)
{
	temp_rec.clear();
	temp_rec_addi.clear();
	ifstream fin;
	string buf,temp1,rang1;
	fin.open(list.c_str(), ios::in);
	if(fin.fail()!=0)return 0;
	int count=0;
	string def_range="_";
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		istringstream www(buf);
		//-> read pair
		if(! (www>>temp1))continue;
		temp_rec.push_back(temp1);
		temp_rec_addi.push_back(def_range);
		count++;
		//-> read range
		int start=(int)temp1.length();
		int len=(int)buf.length();
		int success=0;
		for(int k=start;k<len;k++)
		{
			if(buf[k]==';')
			{
				success=1;
				break;
			}
		}
		if(success==1)
		{
			rang1="";
			if(!getline(www,rang1,';'))continue;
			if(rang1!="")temp_rec_addi[count-1]=rang1;
		}
	}
	return count;
}

//---- transfer float[][] to XYZ format ----//
void Transfer_Float_To_XYZ(float *input,XYZ *output,int moln)
{
	int i;
	for(i=0;i<moln;i++)
	{
		output[i].X=input[3*i+0];
		output[i].Y=input[3*i+1];
		output[i].Z=input[3*i+2];
	}
}

//---- calculate reference score -----//
int Calculate_Reference_Score(CLEFAPS_Main &clepaps,double tms_cut,
	string &wsnam1,string &wsnam2,int ref_id,string &ret_str,double &ret_val)
{
	//assign mol1
	TM_MOLN1=(int)reference_ami[ref_id].length();
	strcpy(TM_AMI1,reference_ami[ref_id].c_str());
	strcpy(TM_CLE1,reference_cle[ref_id].c_str());
	Transfer_Float_To_XYZ(reference_mca[ref_id],TM_MOL1,TM_MOLN1);
	Transfer_Float_To_XYZ(reference_mcb[ref_id],TM_MCB1,TM_MOLN1);
	//calculate score
	int retv=DeepAlign_search(clepaps,wsnam1,wsnam2,ret_str,ret_val,tms_cut);
	return retv;
}


//----- simply load PDB file -----//
//========== process PDB ============//
char WWW_Three2One_III(const char *input)
{
	int i;
	int len;
	int result;
	//encoding
	len=(int)strlen(input);
	if(len!=3)return 'X';
	result=0;
	for(i=0;i<len;i++)result+=(input[i]-'A')*(int)pow(1.0*26,1.0*i);
	//switch
	switch(result)
	{
		case 286:return 'A';
		case 4498:return 'R';
		case 9256:return 'N';
		case 10608:return 'D';
		case 12794:return 'C';
		case 9080:return 'Q';
		case 13812:return 'E';
		case 16516:return 'G';
		case 12383:return 'H';
		case 2998:return 'I';
		case 13635:return 'L';
		case 12803:return 'K';
		case 12960:return 'M';
		case 2901:return 'F';
		case 9921:return 'P';
		case 11614:return 'S';
		case 11693:return 'T';
		case 10601:return 'W';
		case 12135:return 'Y';
		case 7457:return 'V';
		default:return 'X';
	}
}
int Simply_Load_PDB(string &pdbfile,XYZ *mca,XYZ *mcb,char *ami) //->from .pdb file
{
	//--- list for mapping ---//
	map<string, int > ws_mapping;
	map<string, int>::iterator iter;
	ws_mapping.clear();
	ifstream fin;
	string buf,temp,name;
	//read
	fin.open(pdbfile.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		return 0;
	}
	int len;
	int count=0;
	int cb_count=0;
	string wstmp;
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
		if(temp!="ATOM" && temp!="HETA")continue;
		//check CA
		temp=buf.substr(13,2);
		wstmp=temp;
		if(temp!="CA" && temp!="CB")continue;
		//record name
		name=buf.substr(21,6);
		iter = ws_mapping.find(name);
		if(iter != ws_mapping.end())
		{
			if( wstmp=="CB"  )goto next;
			continue;
		}
		else
		{
//			if(first==1)
			{
				count++;
				ws_mapping.insert(map < string, int >::value_type(name, count));
				//record amino acid
				temp=buf.substr(17,3);
				char c=WWW_Three2One_III(temp.c_str());
				ami[count-1]=c;
				mcb[count-1].X=-99999;
				mcb[count-1].Y=-99999;
				mcb[count-1].Z=-99999;
			}
			//record coordinate
			if(wstmp=="CA")
			{
				temp=buf.substr(30,8);
				mca[count-1].X=atof(temp.c_str());
				temp=buf.substr(38,8);
				mca[count-1].Y=atof(temp.c_str());
				temp=buf.substr(46,8);
				mca[count-1].Z=atof(temp.c_str());
			}
next:
			if(wstmp=="CB")
			{
				temp=buf.substr(30,8);
				mcb[count-1].X=atof(temp.c_str());
				temp=buf.substr(38,8);
				mcb[count-1].Y=atof(temp.c_str());
				temp=buf.substr(46,8);
				mcb[count-1].Z=atof(temp.c_str());
				cb_count++;
			}
		}
	}
	ami[count]='\0';
	if(count==cb_count)return count;
	else return -1*count;
}
void Generate_CB_Simp(Confo_Beta &confo_beta,int moln,char *ami,char *cle,XYZ *mol,XYZ *mcb)
{
	int k;
	double dist;
	int ws_correct=1; //default: OK
	for(k=0;k<moln;k++)
	{
		if( IsZero( mcb[k].X + 99999 ) )
		{
			ws_correct=0;
			ws_reco[k]=0;
		}
		else
		{
			dist=mcb[k].distance_square(mol[k]);
			if(dist<1.0||dist>4.0)
			{
				ws_correct=0;
				ws_reco[k]=0;
			}
			else ws_reco[k]=1;
		}
	}
	if(ws_correct==0)
	{
		strcpy(tmp_ami,ami);
		for(k=0;k<moln;k++)if(tmp_ami[k]=='G')tmp_ami[k]='A';
		confo_beta.Recon_Beta_21(mol,tmp_mcb,moln,tmp_ami,cle);
		for(k=0;k<moln;k++)if(ws_reco[k]==0)mcb[k]=tmp_mcb[k];
	}
}

//============== get PDB_File length ============//__130830__//
//[note]: only pdb_BC100 file format is accepted !!
int Get_PDB_File_Len(string &pdbfile) //-> only suitable for pdb_BC100 pdb_file
{
	//--- list for mapping ---//
	map<string, int > ws_mapping;
	map<string, int>::iterator iter;
	//read
	ifstream fin;
	string buf,temp;
	fin.open(pdbfile.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"pdbfile %s not found!!\n",pdbfile.c_str());
		exit(-1);
	}
	//process
	int len;
	int count=0;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		len=(int)buf.length();
		if(len<3)continue;
		//check TER
//		temp=buf.substr(0,3);
//		if(temp=="TER"||temp=="END")break;
		//check ATOM
		if(len<4)continue;
		temp=buf.substr(0,4);
		if(temp!="ATOM" && temp!="HETA")continue;
		//check CA
		temp=buf.substr(13,2);
		if(temp!="CA")continue;
		//record name
		temp=buf.substr(21,6);
		iter = ws_mapping.find(temp);
		if(iter != ws_mapping.end())continue;
		count++;
		ws_mapping.insert(map < string, int >::value_type(temp, count));
	}
	//return
	return count;
}
int Return_Max_Length(string &list,string &root)
{
	//read
	ifstream fin;
	string buf,temp;
	fin.open(list.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"list %s not found!!\n",list.c_str());
		exit(-1);
	}
	//process
	vector <string> nam_rec;
	int count=0;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		nam_rec.push_back(buf);
		count++;
	}
	fin.clear();
	fin.close();
	//record len
	vector <int> len_rec;
	len_rec.resize(count);
	//#pragma omp parallel for schedule(dynamic)
	for(int i=0;i<count;i++)
	{
		string file=root+"/"+nam_rec[i]+".pdb";
		len_rec[i]=Get_PDB_File_Len(file);
	}
	//get max_len
	int maxlen=0;
	for(int i=0;i<count;i++)
	{
		if(len_rec[i]>maxlen)maxlen=len_rec[i];
	}
	//return
	return maxlen;
}

//-------------- read region mask file -----------------//
//-> example:   1111111111122222222222222
int Read_Region_Mask(string &in, vector <int> &out)
{
	//read
	ifstream fin;
	string buf,temp;
	fin.open(in.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"in file %s not found!!\n",in.c_str());
		exit(-1);
	}
	//process
	out.clear();
	if(!getline(fin,buf,'\n'))
	{
		fprintf(stderr,"in file %s format bad!!\n",in.c_str());
		exit(-1);
	}
	//calc
	out.resize(buf.length());
	for(int i=0;i<(int)buf.length();i++)out[i]=buf[i]-'0';
	//return
	return buf.length();
}





//------- usage ------//
void Usage(void)
{
/*
	fprintf(stderr,"DeepAlign V1.10 [Dec-30-2013] \n");
	fprintf(stderr,"-------------- REFERENCE ------------------------------\n");
	fprintf(stderr,"Sheng Wang, Jianzhu Ma, Jian Peng and Jinbo Xu.\n");
	fprintf(stderr,"   PROTEIN STRUCTURE ALIGNMENT BEYOND SPATIAL PROXIMITY \n");
	fprintf(stderr,"                Scientific Reports, 3, 1448, (2013) \n");
	fprintf(stderr,"============== USAGE:1 pairwise alignment =============\n");
	fprintf(stderr,"./DeepAlign <-t protein_A> <-q protein_B> [-o out_name] \n");
	fprintf(stderr,"    [-a alignment] [-x range_A] [-y range_B] \n");
	fprintf(stderr,"Example: ./DeepAlign -t 1col.pdb -q 1cpc.pdb -o testout \n");
	fprintf(stderr,"             -x A:5-121 -y L:1-132 \n");
	fprintf(stderr,"[range_A or range_B] -> the optional residue range. \n"); 
	fprintf(stderr,"if not specified, all residues in the 1st chain are used.\n");
	fprintf(stderr,"[alignment] -> the alignment for initial superposition. \n");
	fprintf(stderr,"[out_name]  -> the output name for the alignment files. \n");
	fprintf(stderr,"    Default is the basenames of the two input files. \n");
	fprintf(stderr,"============== USAGE:2 database search ================\n");
	fprintf(stderr,"./DeepAlign <-f search_list> <-r root> <-q protein_B> \n");
//	fprintf(stderr,"    [-m f_on] [-b f_out] [-c f_thres] [-i simp_in] \n");
	fprintf(stderr,"Example: ./DeepAlign -f bc100_list -r pdb_BC100/ -q 1col.pdb\n");
//	fprintf(stderr,"             -m 1 -b filter_out -c 0.5 \n");
//	fprintf(stderr,"[f_on]-> CLE filter. 0: no filter, [1]: apply filter, \n");
//	fprintf(stderr,"    2: apply reference-based filter, 3: only filter \n");
//	fprintf(stderr,"[f_out]-> CLE filter temporary file out, default==NULL \n");
//	fprintf(stderr,"[f_thres]-> TMscore for CLE filter, default==0.5 \n");
//	fprintf(stderr,"[simp_in]    -> 0: complex pdb format, [1]: simple format\n");
	fprintf(stderr,"============== OPTIONS ================================\n");
//	fprintf(stderr,"   [-n normal_len] [-j kill_frag] [-k kill_gap] [-l loc_para]\n");
	fprintf(stderr,"   [-p out_option] [-s score_func] [-u quality] \n");
//	fprintf(stderr,"   [-d out_root] [-w rank_file] [-v verbose] [-e maxsize]\n");
	fprintf(stderr,"-------------------------------------------------------\n");
	fprintf(stderr,"[out_option] -> 0: don't output detail information,\n");
	fprintf(stderr,"    [1]: single solution, 2: multiple solutions \n");
//	fprintf(stderr,"[normal_len] -> score normalized by the following length,\n");
//	fprintf(stderr,"   -2: second, -1: first, [0]: min, or other given length\n");
//	fprintf(stderr,"[kill_gap]   -> 0: don't kill gaps, [1]: kill gaps \n");
	fprintf(stderr,"[score_func] -> 1:dist-score, 2:vect-score, 4:local-score\n");
	fprintf(stderr,"   these scores could be combined, e.g., 3 == dist * vect\n");
	fprintf(stderr,"   and the default score_func is [7], i.e., using all \n");
	fprintf(stderr,"[quality]    -> [0]: fastest, low quality \n");
	fprintf(stderr,"    1: fast, normal quality, 2: slow, high quality \n");
//	fprintf(stderr,"[loc_para]   -> local-score function's parameter p,\n");
//	fprintf(stderr,"   defined as p*BLOSUM + CLESUM, default is [10] \n");
//	fprintf(stderr,"[kill_frag]  -> 0: don't kill small frags, [1]: kill frags\n");
	fprintf(stderr,"-------------------------------------------------------\n");
*/

	fprintf(stderr,"DeepAlign v1.4 [Aug-20-2018] \n");
	fprintf(stderr,"Sheng Wang, Jianzhu Ma, Jian Peng and Jinbo Xu.\n");
	fprintf(stderr,"   PROTEIN STRUCTURE ALIGNMENT BEYOND SPATIAL PROXIMITY\n");
	fprintf(stderr,"                Scientific Reports, 3, 1448, (2013) \n\n");
	fprintf(stderr,"Usage: \n");
	fprintf(stderr,"./DeepAlign protein_1 protein_2 [-x range_1] [-y range_2] [-a alignment] [-o out_name]\n");
	fprintf(stderr,"                                [-p out_option] [-u quality] [-P screenout] [-n normalize_len] \n");
	fprintf(stderr,"                                [-s score_func] [-C distance_cut] [-M multi_cut] \n");
	fprintf(stderr,"                                [-A mask_1] [-B mask_2] \n\n");
	fprintf(stderr,"Required input: \n");
	fprintf(stderr," protein_1:             The 1st input protein file in PDB format. \n");
	fprintf(stderr," protein_2:             The 2nd input protein file in PDB format. \n\n");
	fprintf(stderr,"Options: \n\n");
	fprintf(stderr,"-x range_1:             The residue range for the 1st input protein, (e.g., A:1-150) \n");
	fprintf(stderr,"-y range_2:             The residue range for the 2nd input protein, \n");
	fprintf(stderr,"                        If not specified, all residues in the first chain will be used. See README for more details.\n\n");
	fprintf(stderr,"-a alignment:           Specify an initial alignment in FASTA format from which to optimize it. \n\n");
	fprintf(stderr,"-o out_name:            Specify an output file name for the alignment. If not specified, \n");
	fprintf(stderr,"                        the output file name is derived from the basenames of the two inputs. \n\n");
	fprintf(stderr,"-p out_option:         [0], do not output the alignment files, just screenout. (Set as default) \n");
	fprintf(stderr,"                        1,  output alignment files for single solution. \n");
	fprintf(stderr,"                        2,  output alignment files for multiple solutions.\n\n");
	fprintf(stderr,"-u quality:             0,  fast, low quality. \n");
	fprintf(stderr,"                       [1], slow, high quality. (Set as default) \n");
	fprintf(stderr,"-U quality_degree:      Quality degree for -u 0 only, ranges from 1 to 10. (set 1 as default) \n\n");
	fprintf(stderr,"-P screenout:           0,  simple screenout alignment scores. \n");
	fprintf(stderr,"                       [1], detailed screenout alignment scores. (Set as default) \n\n");
	fprintf(stderr,"-n normalize_len:       Specify a normalization length for the calculation of TMscore,MAXSUB,GDT_TS/HA. In particular,\n");
	fprintf(stderr,"                       [0], the minimal length of the 1st and 2nd input protein. (Set as default)\n");
	fprintf(stderr,"                       -1,  the length of the first input protein. \n");
	fprintf(stderr,"                       -2,  the length of the second input protein. \n\n");
	fprintf(stderr,"-s score_func:          1:distance-score, 2:vector-score, 4:evolution-score; Note that these scores could be combined,\n");
	fprintf(stderr,"                       [7], using all score combinations. (Set as default) \n\n");
	fprintf(stderr,"-C distance_cut:        Specify a distance cutoff to remove residue pairs whose distance exceeds the threshold. \n");
	fprintf(stderr,"                        0,  keep all residue pairs. \n");
	fprintf(stderr,"                       [-1],automatically assign a distance cutoff value according to d0 in TMscore. (Set as default) \n\n");
	fprintf(stderr,"-M multi_cut:           if multiple solution is specified, set a TMscore cutoff for minimal quality of the solution/ \n");
	fprintf(stderr,"                        (default is 0.35) \n\n");
	fprintf(stderr,"-A mask_1:              The residue mask region for the 1st input protein, (e.g., 11111112222222) \n");
	fprintf(stderr,"-B mask_2:              The residue mask region for the 2nd input protein. \n");
	fprintf(stderr,"-K strict_check:        Enforce that all mask regions shall be aligned above a ratio. (set 0 as NOT to check by default) \n\n");
	fprintf(stderr,"Simple screenout description (please refer to README file for more details):\n");
	fprintf(stderr,"   name1 name2 len1 len2 -> BLOSUM CLESUM DeepScore -> LALI RMSDval TMscore -> MAXSUB GDT_TS GDT_HA -> SeqID nLen dCut uGDT \n");

}

//--- check argument --//__140530__//
int Check_Argument(int argc,char **argv)
{
	int retv=0;
	int hast=0;
	int hasq=0;
	//---- process argument ----//
	for(int i=0;i<argc;i++)
	{
		string tmp=argv[i];
		if(tmp=="-t")hast=1;
		if(tmp=="-q")hasq=1;
	}
	//---- final check -----//
	retv=hast+hasq;
	return retv;
}

//============== main ===============//
int main(int argc,char **argv)
{
	Out_PDB=0;     //-> default, not output linear PDB
	//---- argument check ----//
	int NEWorOLD=0;
	{
		int retv;
		//check input style
		retv=Check_Argument(argc,argv);
		if(retv==0) //-> new style: i.e., ./DeepAlign pdb1 pdb2
		{
			NEWorOLD=1;
		}
		else        //-> old style: i.e., ./DeepAlign -t pdb1 -q pdb2
		{
			NEWorOLD=0;
		}
	}
	
	//================= DeepAlign ===================//__110720__//
	{
		if(argc<3)
		{
			Usage();
			exit(-1);
		}
		//maxsize
		int DEEPALIGN_MAXSIZE=PROT_MAX_NUM;
		int SIMPLY_LOAD=1;  // [1]: simply load; 0: complex load, ( only for Database_Search )
		//job style
		int job_style_search=0;
		int job_style_proc=0;
		int job_style_pair=0;
		//cle filter
		//-> 0: no filter;  1: filter,  2: add reference,  3: only filter
		int CLE_Filter=1;            //default: use CLE_Filter to haste the search process
		string CLE_Filter_out="";    //default: don't output CLE_Filter temporary result
		double CLE_Filter_thres=0.5; //default: use TMscore=0.5 for filter
		//basic input
		string template_list="";
		string template_root="./";
		string name1="";
		string name2="";
		string range1="_";
		string range2="_";
		string ali_file="";
		string out="";
		//verbose
		int verbose=0;
		int rankbose=0;
		string rank_file="";
		string out_root="./";
		//options
		int Out_More=0;   // screen out
		int Out_Screen=1; // detailed screen out
		int Normalize=0;
		int Kill_Gaps=1;
		int Score_Func=7;
		int quality_level=1;
		int quality_degree=1;
		double loc_para=10;
		int Kill_Frag=1;   // we kill frags by default
		double Distance_Cutoff=-1; // we cut the previous aligned positions exceeding this distance cutoff
		double Multi_Cut = 0.35;   // we kill multi-solution below this value
		//reference z_score
		double refer_z_mean=-1;
		double refer_z_vari=-1;
		//mask region file
		string mask_file1="";
		string mask_file2="";
		double Strict_Check=-1;

		//---- process argument ----//
		extern char* optarg;
		int c=0;
		if(NEWorOLD==0) //old style
		{
			while((c=getopt(argc,argv,"g:h:e:i:f:z:r:m:b:c:t:q:a:x:y:o:O:p:P:n:k:s:u:U:l:j:d:w:C:M:A:B:K:v"))!=EOF)
			{
				switch(c) 
				{
					//----- reference z_score ---//
					case 'g':
						refer_z_mean = atof(optarg);
						break;
					case 'h':
						refer_z_vari = atof(optarg);
						break;
	
					//----- max size -----//
					case 'e':
						DEEPALIGN_MAXSIZE = atoi(optarg);
						break;
					case 'i':
						SIMPLY_LOAD = atoi(optarg);
						break;
	
					//----- list job -----//
					case 'f':
						template_list = optarg;
						job_style_search=1; //-> do search job against a template list
						break;
					case 'z':
						template_list = optarg;
						job_style_proc=1;   //-> do <target query> process list
						break;
					case 'r':
						template_root = optarg;
						break;
	
					//----- filter job ----//
					case 'm':
						CLE_Filter = atoi(optarg);  //0, 1, 2 or 3
						break;
					case 'b':
						CLE_Filter_out = optarg;    //default ""
						break;
					case 'c':
						CLE_Filter_thres = atof(optarg); //default 0.5
						break;
						
					//---- pairsie job ---//
					case 't':
						name1 = optarg;
						job_style_pair=1;   //-> do pairwise alignment
						break;
					case 'q':
						name2 = optarg;
						break;
					case 'a':
						ali_file = optarg;
						break;
					case 'x':
						range1 = optarg;
						break;
					case 'y':
						range2 = optarg;
						break;
					case 'o':
						out = optarg;
						break;
					case 'O':
						//-> [1] for output PDB in linear co-ordinate, otherwise not output
						Out_PDB = atoi(optarg);  
						break;
						
					//---- options ----//
					case 'p':
						Out_More = atoi(optarg);
						break;
					case 'P':
						Out_Screen = atoi(optarg);
						break;
					case 'n':
						Normalize = atoi(optarg);
						break;
					case 'k':
						Kill_Gaps = atoi(optarg);
						break;
					case 's':
						Score_Func = atoi(optarg);
						break;
					case 'u':
						quality_level = atoi(optarg);
						break;
					case 'U':
						quality_degree = atoi(optarg);
						break;
					case 'l':
						loc_para = atof(optarg);
						break;
					case 'j':
						Kill_Frag = atoi(optarg);
						break;
					case 'C':
						Distance_Cutoff = atof(optarg);
						break;
					case 'M':
						Multi_Cut = atof(optarg);
						break;
	
					//---- verbose related ---//
					case 'd':
						out_root = optarg;
						break;
					case 'w':
						rank_file = optarg;
						rankbose=1;
						break;
					case 'v':
						verbose=1;
						break;

					//---- region mask file ---//
					case 'A':
						mask_file1 = optarg;
						break;
					case 'B':
						mask_file2 = optarg;
						break;
					case 'K':
						Strict_Check = atof(optarg);
						break;

					//----- default ----//
					default:
						Usage();
						exit(-1);
				}
			}
		}
		else //new style
		{
			job_style_pair=1;
			name1 = argv[1];
			name2 = argv[2];
			while((c=getopt(argc,argv,"e:a:x:y:o:O:p:P:n:k:s:u:U:l:j:C:M:A:B:K:"))!=EOF)
			{
				switch(c) 
				{
					//----- max size -----//
					case 'e':
						DEEPALIGN_MAXSIZE = atoi(optarg);
						break;
						
					//---- pairsie job ---//
					case 'a':
						ali_file = optarg;
						break;
					case 'x':
						range1 = optarg;
						break;
					case 'y':
						range2 = optarg;
						break;
					case 'o':
						out = optarg;
						break;
					case 'O':
						//-> [1] for output PDB in linear co-ordinate, otherwise not output
						Out_PDB = atoi(optarg);
						break;
						
					//---- options ----//
					case 'p':
						Out_More = atoi(optarg);
						break;
					case 'P':
						Out_Screen = atoi(optarg);
						break;
					case 'n':
						Normalize = atoi(optarg);
						break;
					case 'k':
						Kill_Gaps = atoi(optarg);
						break;
					case 's':
						Score_Func = atoi(optarg);
						break;
					case 'u':
						quality_level = atoi(optarg);
						break;
					case 'U':
						quality_degree = atoi(optarg);
						break;
					case 'l':
						loc_para = atof(optarg);
						break;
					case 'j':
						Kill_Frag = atoi(optarg);
						break;
					case 'C':
						Distance_Cutoff = atof(optarg);
						break;
					case 'M':
						Multi_Cut = atof(optarg);
						break;

					//---- region mask file ---//
					case 'A':
						mask_file1 = optarg;
						break;
					case 'B':
						mask_file2 = optarg;
						break;
					case 'K':
						Strict_Check = atof(optarg);
						break;

					//----- default ----//
					default:
						Usage();
						exit(-1);
				}
			}
		}

		//----- check argument ----//
		//-> check refer_z_score
		int refer_z_score=0; //don't apply input refer_z_score
		if(refer_z_mean>0 || refer_z_vari>0)refer_z_score=1;
		//-> check verbose
		int verb_total=rankbose+verbose;
		//-> check job_style
		int job_total=job_style_search + job_style_proc + job_style_pair;
		if(job_total!=1)
		{
			fprintf(stderr,"ERROR: must input a valid job_style, either search_list, proc_list or pairwise \n"); 
			exit(-1);
		}
		//-> check query_mode
		if(job_style_search==1 || job_style_pair==1) //-> search_list or pairwise, then check query_file
		{
			ifstream fin;
			fin.open(name2.c_str(), ios::in);
			if(fin.fail()!=0)
			{
				fprintf(stderr,"ERROR: query_file %s must not be empty \n",name2.c_str());
				exit(-1);
			}
		}
		//-> check list_mode
		if(job_style_search==1 || job_style_proc==1) //-> search_list or proc_list, then check template_list
		{
			ifstream fin;
			string buf;
			fin.open(template_list.c_str(), ios::in);
			if(fin.fail()!=0)
			{
				fprintf(stderr,"ERROR: template_list %s must not be empty \n",template_list.c_str());
				exit(-1);
			}
			int count=0;
			for(;;)
			{
				if(!getline(fin,buf,'\n'))break;
				count++;
			}
			if(count==0)
			{
				exit(-1);
			}
		}
		//-> check DEEPALIGN_MAXSIZE
		if(job_style_search==1)
		{
			//check search list to get maxsize
			if(DEEPALIGN_MAXSIZE<0)
			{
				DEEPALIGN_MAXSIZE=Return_Max_Length(template_list,template_root);
			}
		}
		//-> check CLE filter
		if(CLE_Filter<0 || CLE_Filter>3)
		{
			fprintf(stderr,"ERROR: CLE_Filter must be 0 to 3 \n");
			exit(-1);
		}
		if(CLE_Filter_thres<0 || CLE_Filter_thres>1)
		{
			fprintf(stderr,"ERROR: CLE_Filter_thres must be a real number between 0 to 1 \n");
			exit(-1);
		}
		//-> check other arguments
		if(Out_More<0 || Out_More>2)
		{
			fprintf(stderr,"ERROR: Out_More must be integer between 0 to 2 \n");
			exit(-1);
		}
		if(Out_Screen<0 || Out_Screen>1)
		{
			fprintf(stderr,"ERROR: Out_Screen must be integer between 0 to 1 \n");
			exit(-1);
		}
		if(Normalize<-2)
		{
			fprintf(stderr,"ERROR: normal_len must be interger >= -2 \n");
			exit(-1);
		}
		if(Kill_Gaps<0 || Kill_Gaps>1)
		{
			fprintf(stderr,"ERROR: kill_gap must be integer between 0 to 1 \n");
			exit(-1);
		}
		if(Score_Func<1 || Score_Func>7)
		{
			fprintf(stderr,"ERROR: score_func must be integer between 1 to 7 \n");
			exit(-1);
		}
		if(quality_level<0 || quality_level>1)
		{
			fprintf(stderr,"ERROR: quality must be integer between 0 to 1 \n");
			exit(-1);
		}
		if(quality_degree<1 || quality_degree>10)
		{
			fprintf(stderr,"ERROR: quality_degree must be integer between 1 to 10 \n");
			exit(-1);
		}
		if(loc_para<0)
		{
			fprintf(stderr,"ERROR: loc_para must > 0 \n");
			exit(-1);
		}
		if(Kill_Frag<0 || Kill_Frag>1)
		{
			fprintf(stderr,"ERROR: kill_frag must be integer between 0 to 1 \n");
			exit(-1);
		}
		if(SIMPLY_LOAD<0 || SIMPLY_LOAD>1)
		{
			fprintf(stderr,"ERROR: SIMPLY_LOAD must be integer between 0 to 1 \n");
			exit(-1);
		}
		if(Multi_Cut<0 || Multi_Cut>1)
		{
			fprintf(stderr,"ERROR: Multi_Cut must be float between 0 to 1 \n");
			exit(-1);
		}
		//-> check mask file
		if( (mask_file1=="" && mask_file2!="") || (mask_file2=="" && mask_file1!="") )
		{
			fprintf(stderr,"mask_file1 or mask_file2 shall not exist single \n");
			exit(-1);
		}
		if(Strict_Check>1)
		{
			fprintf(stderr,"ERROR: Strict_Check must be float between 0 to 1 \n");
			exit(-1);
		}


		//=========== DeepAlign process =========//
		int retv;
		string wsnam1,wsnam2;
		//rank file
		FILE *fp=0;
		if(rank_file.length()!=0)
		{
			fp=fopen(rank_file.c_str(),"wb");
			if(fp==0)
			{
				fprintf(stderr,"ERROR: file %s can't be opened. \n",rank_file.c_str());
				exit(-1);
			}
		}

		//----- pairwise DeepAlign ----//
		if(job_style_pair==1)
		{
			//-> load mol1
			retv=mol_input1.XYZ_Input(name1,range1,0,TM_MOLN1,0,0,0,0,0);
			if(retv!=1)
			{
				fprintf(stderr,"FILE1 NOT FOUND!!!![%s]\n",name1.c_str());
				exit(-1);
			}
			if(TM_MOLN1<=0)
			{
				fprintf(stderr,"FILE1 IS EMPTY!!!![%s]\n",name1.c_str());
				exit(-1);
			}
			//-> load mol2
			retv=mol_input2.XYZ_Input(name2,range2,0,TM_MOLN2,0,0,0,0,0);
			if(retv!=1)
			{
				fprintf(stderr,"FILE2 NOT FOUND!!!![%s]\n",name2.c_str());
				exit(-1);
			}
			if(TM_MOLN2<=0)
			{
				fprintf(stderr,"FILE2 IS EMPTY!!!![%s]\n",name2.c_str());
				exit(-1);
			}
			//-------------- load region mask ------------------//
			int *reg_mask1=0;
			int *reg_mask2=0;
			if(mask_file1!="" && mask_file2!="")
			{
				vector <int> out1;
				int mask1_len=Read_Region_Mask(mask_file1,out1);
				vector <int> out2;
				int mask2_len=Read_Region_Mask(mask_file2,out2);
				//judge
				if(mask1_len!=TM_MOLN1 || mask2_len!=TM_MOLN2)
				{
					fprintf(stderr,"mask1_len %d not equal to TM_MOLN1 %d || mask2_len %d not equal to TM_MOLN2 %d for pair %s %s \n",
						mask1_len,TM_MOLN1,mask2_len,TM_MOLN2,mask_file1.c_str(),mask_file2.c_str());
					exit(-1);
				}
				//new
				reg_mask1=new int[mask1_len];
				for(int z=0;z<mask1_len;z++)reg_mask1[z]=out1[z];
				reg_mask2=new int[mask2_len];
				for(int z=0;z<mask2_len;z++)reg_mask2[z]=out2[z];
			}

			//-> initialization
			TM_Align_Init_WS(TM_MOLN1,TM_MOLN2);
			int wsmax=TM_MOLN1>TM_MOLN2?TM_MOLN1:TM_MOLN2;
			Confo_Beta confo_beta(wsmax);
			CLEFAPS_Main clepaps(wsmax);
			if(reg_mask1!=0 && reg_mask2!=0)
			{
				clepaps.mas1=reg_mask1;
				clepaps.mas2=reg_mask2;
			}
			//-> load structure
			int PRE_LOAD_1=mol_input1.PRE_LOAD;
			int WARNING_out_1=mol_input1.WARNING_out;
			mol_input1.PRE_LOAD=1;
			mol_input1.WARNING_out=0;
			mol_input1.XYZ_Input(name1,range1,0,TM_MOLN1,TM_MOL1,TM_AMI1,TM_CLE1,0,TM_PDB1);
			mol_input1.PRE_LOAD=PRE_LOAD_1;
			mol_input1.WARNING_out=WARNING_out_1;
			int PRE_LOAD_2=mol_input2.PRE_LOAD;
			int WARNING_out_2=mol_input2.WARNING_out;
			mol_input2.PRE_LOAD=1;
			mol_input2.WARNING_out=0;
			mol_input2.XYZ_Input(name2,range2,0,TM_MOLN2,TM_MOL2,TM_AMI2,TM_CLE2,0,TM_PDB2);
			mol_input2.PRE_LOAD=PRE_LOAD_2;
			mol_input2.WARNING_out=WARNING_out_2;
			//-> get_base_name
			getBaseName(name1,wsnam1,'/','.');
			getBaseName(name2,wsnam2,'/','.');
			//-> process CB
			Generate_CB(confo_beta,TM_MOLN1,TM_PDB1,TM_AMI1,TM_CLE1,TM_MOL1,TM_MCB1);
			Generate_CB(confo_beta,TM_MOLN2,TM_PDB2,TM_AMI2,TM_CLE2,TM_MOL2,TM_MCB2);
			clepaps.TM_cb1=TM_MCB1;
			clepaps.TM_cb2=TM_MCB2;
			//-> real process
			string out_str;
			double ret_val;
			//-> check out
			if(out!="" && Out_More<=0)Out_More=1;
			if(out=="" && Out_More>0)out=wsnam1+"-"+wsnam2;
			DeepAlign_main(clepaps,wsnam1,wsnam2,out,ali_file,Out_More,Out_Screen,
				Normalize,Kill_Gaps,Score_Func,quality_level,quality_degree,
				out_root,out_str,ret_val,0,Kill_Frag,Strict_Check,Distance_Cutoff,Multi_Cut);
			//-> output
			if(rankbose==1)
			{
				fprintf(fp,"%s",out_str.c_str());
			}
			else
			{
				printf("%s",out_str.c_str());
			}
			//-> delete
			if(reg_mask1!=0 && reg_mask2!=0)
			{
				delete [] reg_mask1;
				delete [] reg_mask2;
			}
		}
		//----- proc_list DeepAlign ----//
		if(job_style_proc==1)
		{
			//init
			vector <string> temp_list;
			vector <string> targ_list;
			vector <string> temp_list_addi;
			vector <string> targ_list_addi;
			int size=Get_PairList(template_list,temp_list,targ_list,temp_list_addi,targ_list_addi);
			TM_Align_Init_WS(DEEPALIGN_MAXSIZE,DEEPALIGN_MAXSIZE);
			CLEFAPS_Main clepaps(DEEPALIGN_MAXSIZE);
			Confo_Beta confo_beta(DEEPALIGN_MAXSIZE);
			//proc templist
			int outsize=size/100;
			if(outsize==0)outsize=1;
			for(int i=0;i<size;i++)
			{
				//-> load mol1
				name1=template_root+"/"+temp_list[i]+".pdb";
				range1=temp_list_addi[i];
				retv=mol_input.XYZ_Input(name1,range1,0,TM_MOLN1,TM_MOL1,TM_AMI1,TM_CLE1,0,TM_PDB1);
				if(retv!=1)
				{
					fprintf(stderr,"FILE1 NOT FOUND!!!![%s]\n",name1.c_str());
					continue;
				}
				if(TM_MOLN1<=0 || TM_MOLN1>DEEPALIGN_MAXSIZE)
				{
					fprintf(stderr,"FILE1 IS EMPTY or OVERANGE!!!![%s], <%d,%d>; if overange, try -e size to increase \n",
						name1.c_str(),TM_MOLN1,DEEPALIGN_MAXSIZE);
					continue;
				}
				//-> load mol2
				name2=template_root+"/"+targ_list[i]+".pdb";
				range2=targ_list_addi[i];
				retv=mol_input.XYZ_Input(name2,range2,0,TM_MOLN2,TM_MOL2,TM_AMI2,TM_CLE2,0,TM_PDB2);
				if(retv!=1)
				{
					fprintf(stderr,"FILE2 NOT FOUND!!!![%s]\n",name2.c_str());
					continue;
				}
				if(TM_MOLN2<=0 || TM_MOLN2>DEEPALIGN_MAXSIZE)
				{
					fprintf(stderr,"FILE2 IS EMPTY or OVERANGE!!!![%s], <%d,%d>; if overange, try -e size to increase \n",
						name2.c_str(),TM_MOLN2,DEEPALIGN_MAXSIZE);
					continue;
				}
				//-> get_base_name
				getBaseName(name1,wsnam1,'/','.');
				getBaseName(name2,wsnam2,'/','.');
				//-> process CB
				Generate_CB(confo_beta,TM_MOLN1,TM_PDB1,TM_AMI1,TM_CLE1,TM_MOL1,TM_MCB1);
				Generate_CB(confo_beta,TM_MOLN2,TM_PDB2,TM_AMI2,TM_CLE2,TM_MOL2,TM_MCB2);
				clepaps.TM_cb1=TM_MCB1;
				clepaps.TM_cb2=TM_MCB2;
				//-> real process
				string out_str;
				string ali_null="";
				double ret_val;
				string wsout=wsnam1+"-"+wsnam2;
				DeepAlign_main(clepaps,wsnam1,wsnam2,wsout,ali_null,Out_More,Out_Screen,
					Normalize,Kill_Gaps,Score_Func,quality_level,quality_degree,
					out_root,out_str,ret_val,0,Kill_Frag,-1,Distance_Cutoff,Multi_Cut);
				//-> output
				if(rankbose==1)
				{
					fprintf(fp,"%s",out_str.c_str());
				}
				else
				{
					printf("%s",out_str.c_str());
				}
				//verbose output
				if(verb_total==2)
				{
					if(i%outsize==0)fprintf(stderr,".");
				}
			}
			if(verb_total==2)fprintf(stderr,"\n");
		}
		//----- search_list DeepAlign ----//
		if(job_style_search==1)
		{
			//-> load mol2 (query)
			retv=mol_input.XYZ_Input(name2,range2,0,TM_MOLN2,0,0,0,0,0);
			if(retv!=1)
			{
				fprintf(stderr,"QUERY NOT FOUND!!!![%s]\n",name2.c_str());
				exit(-1);
			}
			if(TM_MOLN2<=0)
			{
				fprintf(stderr,"QUERY IS EMPTY!!!![%s]\n",name2.c_str());
				exit(-1);
			}
			if(TM_MOLN2>DEEPALIGN_MAXSIZE)DEEPALIGN_MAXSIZE=TM_MOLN2;
			//-> initialize
			TM_Align_Init_WS(DEEPALIGN_MAXSIZE,TM_MOLN2);
			CLEFAPS_Main clepaps(DEEPALIGN_MAXSIZE);
			Confo_Beta confo_beta(DEEPALIGN_MAXSIZE);
			//-> load structure
			int PRE_LOAD_=mol_input.PRE_LOAD;
			int WARNING_out_=mol_input.WARNING_out;
			mol_input.PRE_LOAD=1;
			mol_input.WARNING_out=0;
			mol_input.XYZ_Input(name2,range2,0,TM_MOLN2,TM_MOL2,TM_AMI2,TM_CLE2,0,TM_PDB2);
			mol_input.PRE_LOAD=PRE_LOAD_;
			mol_input.WARNING_out=WARNING_out_;
			//-> get base_name
			getBaseName(name2,wsnam2,'/','.');
			//-> process CB
			Generate_CB(confo_beta,TM_MOLN2,TM_PDB2,TM_AMI2,TM_CLE2,TM_MOL2,TM_MCB2);
			clepaps.TM_cb2=TM_MCB2;
			//-> load list
			vector <string> temp_list;
			vector <string> temp_list_addi;
			int size=0;
			
			//---- run CLE_Filter or not? ----//__130202__//
			if(CLE_Filter>=1)  //-> apply CLE_Filter 
			{
				//-> scan reference first
				int CLE_Filter_add;
				double z_mean,z_vari;
				int z_count;
				z_mean=0;
				z_vari=1;
				z_count=0;
				if(CLE_Filter>=2 && refer_z_score==0) //1814 reference filter by z_score
				{
					int cur_count=0;
					double cur_mean=0;
					double cur_vari=0;
					string ret_str;
					double ret_val;
					int retv;
					double tmp_sss[REFERENCE_NUM];
					//-> pre2-process list
					int outsize=REFERENCE_NUM/100;
					if(outsize<=2)outsize=2;
					if(verb_total==2)fprintf(stderr,"process@: ");
					for(int i=0;i<REFERENCE_NUM;i++)
					{
						retv=Calculate_Reference_Score(clepaps,CLE_Filter_thres/2,
							reference_nam[i],wsnam2,i,ret_str,ret_val);
						if(retv==1)
						{
							tmp_sss[cur_count]=ret_val;
							cur_count++;
						}
						//verbose output
						if(verb_total==2)
						{
							if(i%outsize==0)fprintf(stderr,".");
						}
					}
					//-> sort
					if(cur_count>3) //omit top3
					{
						//sort
						Fast_Sort <double> fast_sort;
						double *sort_sco=new double[cur_count];
						int *sort_idx=new int[cur_count];
						for(int i=0;i<cur_count;i++)sort_sco[i]=tmp_sss[i];
						fast_sort.fast_sort_1(sort_sco,sort_idx,cur_count);
						//assign
						for(int i=0;i<cur_count;i++)tmp_sss[i]=sort_sco[sort_idx[i]];
						//delete
						delete [] sort_sco;
						delete [] sort_idx;
					}
					//-> estimate z_score
					if(cur_count>3)
					{
						for(int i=3;i<cur_count;i++)
						{
							cur_mean+=tmp_sss[i];
						}
						cur_mean/=(cur_count-3);
						for(int i=3;i<cur_count;i++)
						{
							cur_vari+=(tmp_sss[i]-cur_mean)*(tmp_sss[i]-cur_mean);
						}
						cur_vari=sqrt(cur_vari/(cur_count-3));
						//assign
						z_mean=cur_mean;
						z_vari=cur_vari;
						z_count=cur_count-3;
					}
					if(verb_total==2)fprintf(stderr,"\n");
				}
				if(CLE_Filter>=2 && refer_z_score==1)
				{
					z_mean=refer_z_mean;
					z_vari=refer_z_vari;
					z_count=1;
				}
				double CLE_Filter_thres_real=CLE_Filter_thres>z_mean?CLE_Filter_thres:z_mean;
				//-> scan reference first over

				//output tempfile
				FILE *filter_out=0;
				if(CLE_Filter_out!="")filter_out=fopen(CLE_Filter_out.c_str(),"wb");
				//load list first
				vector <string> temp_list_ori;
				vector <string> temp_list_addi_ori;
				int totnum=Get_TempList(template_list,temp_list_ori,temp_list_addi_ori);
				//-> pre-process list
				int outsize=totnum/100;
				if(outsize<=2)outsize=2;
				if(verb_total==2)fprintf(stderr,"process0: ");
				for(int i=0;i<totnum;i++)
				{
					name1=template_root+"/"+temp_list_ori[i]+".pdb";
					range1=temp_list_addi_ori[i];
					if(SIMPLY_LOAD==0)
					{
						retv=mol_input.XYZ_Input(name1,range1,0,TM_MOLN1,TM_MOL1,TM_AMI1,TM_CLE1,0,TM_PDB1);
					}
					else
					{
						retv=Simply_Load_PDB(name1,TM_MOL1,TM_MCB1,TM_AMI1);
						TM_MOLN1=abs(retv);
						confo_lett.btb_ori(0,0,0,TM_MOLN1,TM_MOL1,TM_CLE1);
						TM_CLE1[TM_MOLN1]='\0';
					}
					if(retv==0)
					{
						fprintf(stderr,"FILE1 NOT FOUND!!!![%s]\n",name1.c_str());
						continue;
					}
					if(TM_MOLN1<=0 || TM_MOLN1>DEEPALIGN_MAXSIZE)
					{
						fprintf(stderr,"FILE1 IS EMPTY or OVERANGE!!!![%s], <%d,%d>; if overange, try -e size to increase \n",
							name1.c_str(),TM_MOLN1,DEEPALIGN_MAXSIZE);
						continue;
					}
					//-> get base_name
					getBaseName(name1,wsnam1,'/','.');
					//-> apply CLE_Filter DeepAlign_search
					string out_str;
					double ret_val;
					retv=DeepAlign_search(clepaps,wsnam1,wsnam2,out_str,ret_val,CLE_Filter_thres_real/2);
					//-> out check
					if(retv==1) 
					{
						CLE_Filter_add=1;
						if(CLE_Filter>=2 && z_count>0)
						{
							if( (ret_val-z_mean)/z_vari < CLE_Filter_thres_real*2 )CLE_Filter_add=0;
						}
						if(ret_val > CLE_Filter_thres_real && CLE_Filter_add==1) //-> record this template
						{
							//output
							if(filter_out!=0)fprintf(filter_out,"%s",out_str.c_str());
							//record
							temp_list.push_back(temp_list_ori[i]);
							temp_list_addi.push_back(temp_list_addi_ori[i]);
							size++;
						}
					}
					//verbose output
					if(verb_total==2)
					{
						if(i%outsize==0)fprintf(stderr,".");
					}
				}
				if(verb_total==2)fprintf(stderr,"\n");
				//break check
				if(CLE_Filter>=3)  //only filter
				{
					exit(-1);
				}
			}
			else
			{
				size=Get_TempList(template_list,temp_list,temp_list_addi);
			}
			//proc templist
			int outsize=size/100;
			if(outsize<=2)outsize=2;
			if(verb_total==2)fprintf(stderr,"process1: ");
			vector <string> out_rec; //-> record out_string
			vector <double> out_sco; //-> record out_score
			int sort_size=0;         //-> record size
			for(int i=0;i<size;i++)
			{
				//-> load mol1
				name1=template_root+"/"+temp_list[i]+".pdb";
				range1=temp_list_addi[i];
				if(SIMPLY_LOAD==0)
				{
					retv=mol_input.XYZ_Input(name1,range1,0,TM_MOLN1,TM_MOL1,TM_AMI1,TM_CLE1,0,TM_PDB1);
				}
				else
				{
					retv=Simply_Load_PDB(name1,TM_MOL1,TM_MCB1,TM_AMI1);
					TM_MOLN1=abs(retv);
				}
				if(retv==0)
				{
					fprintf(stderr,"FILE1 NOT FOUND!!!![%s]\n",name1.c_str());
					continue;
				}
				if(TM_MOLN1<=0 || TM_MOLN1>DEEPALIGN_MAXSIZE)
				{
					fprintf(stderr,"FILE1 IS EMPTY or OVERANGE!!!![%s], <%d, %d>; if overange, try -e size to increase \n",
						name1.c_str(),TM_MOLN1,DEEPALIGN_MAXSIZE);
					continue;
				}
				//-> get_base_name
				getBaseName(name1,wsnam1,'/','.');
				//-> process CB
				if(SIMPLY_LOAD==0)
				{
					Generate_CB(confo_beta,TM_MOLN1,TM_PDB1,TM_AMI1,TM_CLE1,TM_MOL1,TM_MCB1);
				}
				else
				{
					confo_lett.btb_ori(0,0,0,TM_MOLN1,TM_MOL1,TM_CLE1);
					TM_CLE1[TM_MOLN1]='\0';
					Generate_CB_Simp(confo_beta,TM_MOLN1,TM_AMI1,TM_CLE1,TM_MOL1,TM_MCB1);

				}
				clepaps.TM_cb1=TM_MCB1;
				//-> real process
				string out_str;
				string ali_null="";
				double ret_val;
				string wsout=wsnam1+"-"+wsnam2;
				DeepAlign_main(clepaps,wsnam1,wsnam2,wsout,ali_null,Out_More,Out_Screen,
					Normalize,Kill_Gaps,Score_Func,quality_level,quality_degree,
					out_root,out_str,ret_val,SIMPLY_LOAD,Kill_Frag,-1,Distance_Cutoff,Multi_Cut);
				//-> output record
				out_rec.push_back(out_str);
				out_sco.push_back(ret_val);
				sort_size++;
				//verbose output
				if(verb_total==2)
				{
					if(i%outsize==0)fprintf(stderr,".");
				}
			}
			//verbose output
			if(verb_total==2)fprintf(stderr,"\n");
			//-> sort and output
			{
				//sort
				Fast_Sort <double> fast_sort;
				double *sort_sco=new double[sort_size];
				int *sort_idx=new int[sort_size];
				for(int i=0;i<sort_size;i++)sort_sco[i]=out_sco[i];
				fast_sort.fast_sort_1(sort_sco,sort_idx,sort_size);
				//output
				for(int i=0;i<sort_size;i++)
				{
					int index=sort_idx[i];
					if(rankbose==1)
					{
						fprintf(fp,"%s",out_rec[index].c_str());
					}
					else
					{
						printf("%s",out_rec[index].c_str());
					}
				}
				//delete
				delete [] sort_sco;
				delete [] sort_idx;
			}
		}

		//=========== final delete =========//
		TM_Align_dele_WS();
	}

	//return 0
	exit(0);
}
