#include <string>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <map>
#include "getopt.h"
#include "TM_align.h"
#include "Confo_Beta.h"
#include "Mol_Load.h"
#include "Mol_Out.h"
using namespace std;

Mol_Load mol_input;
Mol_Out mol_out;

//-> [1] for output PDB in linear co-ordinate, otherwise not output
int Out_PDB;          //-> output PDB


//--------- get min and max residue number from a given PDB file --------//
int Extract_MinMax_ResNum(string &pdbfile, int &minres,int &maxres)
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(pdbfile.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"pdbfile %s not found \n",pdbfile.c_str());
		exit(-1);
	}
	int resnum;
	minres=9999999;
	maxres=-9999999;
	int len;
	string wstmp;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		len=(int)buf.length();
		if(len<3)continue;
		//check ATOM
		if(len<4)continue;
		temp=buf.substr(0,4);
		if(temp!="ATOM" && temp!="HETA")continue;
		//check CA
		temp=buf.substr(13,2);
		if(temp!="CA")continue;
		//extract resnum
		temp=buf.substr(22,4);
		resnum=atoi(temp.c_str());
		if(resnum<minres)minres=resnum;
		if(resnum>maxres)maxres=resnum;
	}
	//check
	if(minres>maxres)
	{
		fprintf(stderr,"minres %d larger than maxres %d\n",minres,maxres);
		exit(-1);
	}
	//return
	return maxres-minres+1;
}


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

//---------- dynamic programming ----------//
int WWW_Advance_Align_Dyna_Prog_Double(int n1,int n2,const vector<double> &score,
								   double GAP_OPEN1,double GAP_EXT1,double GAP_OPEN2,double GAP_EXT2,
								   double GAP_HEAD1,double GAP_TAIL1,double GAP_HEAD2,double GAP_TAIL2,
								   vector<pair<int,int> > & alignment,double &ali_sco)
{
	int i,j;
	//input
	int m = n1 + 1;  // +1 to account for the extra row,col in
	int n = n2 + 1;  // the DP matrices corresponding to gaps
	int DP_maximal=n;
	int IN_maximal=n2;
	//const value
	const int _H_  = 0;
	const int _S_  = 1;
	const int _V_  = 2;

	//create D and M
	vector <int> D[3];      // the path (directions) matrix
	vector <double> M[3];   // the current scores (values) matrix
	//resize(m,n)
	for (i = 0; i < 3; ++i) 
	{
		D[i].resize(m*n);
		M[i].resize(m*n);
	}
	//init()
	double WS_MIN=-1000000;
	D[_S_][0*DP_maximal+ 0] = -1;
	D[_H_][0*DP_maximal+ 0] = -1;
	D[_V_][0*DP_maximal+ 0] = -1;
	M[_S_][0*DP_maximal+ 0] = 0;
	M[_H_][0*DP_maximal+ 0] = WS_MIN;
	M[_V_][0*DP_maximal+ 0] = WS_MIN;
	for (i = 1; i < m; i++) 
	{
		D[_S_][i*DP_maximal+ 0] = _V_;
		D[_H_][i*DP_maximal+ 0] = _V_;
		D[_V_][i*DP_maximal+ 0] = _V_;
		M[_S_][i*DP_maximal+ 0] = WS_MIN;
		M[_H_][i*DP_maximal+ 0] = WS_MIN;
		M[_V_][i*DP_maximal+ 0] = i*GAP_HEAD1; //-(Params::GAP_OPEN + (i-1)*Params::GAP_EXT);
	}
	for (j = 1; j < n; j++) 
	{
		D[_S_][0*DP_maximal+ j] = _H_;
		D[_H_][0*DP_maximal+ j] = _H_;
		D[_V_][0*DP_maximal+ j] = _H_;
		M[_S_][0*DP_maximal+ j] = WS_MIN;
		M[_H_][0*DP_maximal+ j] = j*GAP_HEAD2; //-(Params::GAP_OPEN + (j-1)*Params::GAP_EXT);
		M[_V_][0*DP_maximal+ j] = WS_MIN;
	}
	//fill(firstSeq, secondSeq, distFunc);
	double gap_open;
	double gap_ext;
	double v1,v2,v3;
	double dist;
	for (i = 1; i < m; i++) 
	{
		for (j = 1; j < n; j++) 
		{
			//condition upper
			if(j==n-1)
			{
				gap_open=GAP_TAIL1;
				gap_ext=GAP_TAIL1;
			}
			else
			{
				gap_open=GAP_OPEN1;
				gap_ext=GAP_EXT1;
			}
			v1 = M[_V_][(i-1)*DP_maximal+ j] + gap_ext;
			v2 = M[_S_][(i-1)*DP_maximal+ j] + gap_open;
			v3 = M[_H_][(i-1)*DP_maximal+ j] + gap_open;
			M[_V_][i*DP_maximal+ j] = std::max(v1, std::max(v2, v3));
			if (M[_V_][i*DP_maximal+ j] == v1) D[_V_][i*DP_maximal+ j] = _V_;
			else if(M[_V_][i*DP_maximal+ j] == v2) D[_V_][i*DP_maximal+ j] = _S_;
			else D[_V_][i*DP_maximal+ j] = _H_;
			//condition left
			if(i==m-1)
			{
				gap_open=GAP_TAIL2;
				gap_ext=GAP_TAIL2;
			}
			else
			{
				gap_open=GAP_OPEN2;
				gap_ext=GAP_EXT2;
			}
			v1 = M[_H_][i*DP_maximal+ j-1] + gap_ext;
			v2 = M[_S_][i*DP_maximal+ j-1] + gap_open;
			v3 = M[_V_][i*DP_maximal+ j-1] + gap_open;
			M[_H_][i*DP_maximal+ j] = std::max(v1, std::max(v2, v3));
			if (M[_H_][i*DP_maximal+ j] == v1) D[_H_][i*DP_maximal+ j] = _H_;
			else if(M[_H_][i*DP_maximal+ j] == v2) D[_H_][i*DP_maximal+ j] = _S_;
			else D[_H_][i*DP_maximal+ j] = _V_;
			//condition diag
			dist = score.at((i-1)*IN_maximal+ j-1);  //Params::K - distFunc(firstSeq[i-1], secondSeq[j-1]);
			v1 = M[_V_][(i-1)*DP_maximal+ j-1] + dist;
			v2 = M[_H_][(i-1)*DP_maximal+ j-1] + dist;
			v3 = M[_S_][(i-1)*DP_maximal+ j-1] + dist;
			M[_S_][i*DP_maximal+ j] = std::max(v1, std::max(v2, v3));
			if (M[_S_][i*DP_maximal+ j] == v3) D[_S_][i*DP_maximal+ j] = _S_;
			else if (M[_S_][i*DP_maximal+ j] == v1) D[_S_][i*DP_maximal+ j] = _V_;
			else D[_S_][i*DP_maximal+ j] = _H_;
		}
	}
	//build(ali, firstSeq, secondSeq, distFunc);
	i = m-1;
	j = n-1;
	v1=M[_V_][i*DP_maximal+ j];
	v2=M[_H_][i*DP_maximal+ j];
	v3=M[_S_][i*DP_maximal+ j];
	double maximal = std::max(v1, std::max(v2, v3));
	int k = -1;
	if(v3==maximal)k = _S_;
	else if(v2==maximal)k = _H_;
	else k = _V_;
	//trace_back
	alignment.clear();
	int count = 0;
	int matches = 0;
	int cur_case=k;
	int pre_case;
	for(;;)
	{
		if(i==0||j==0)break;
		pre_case=D[cur_case][i*DP_maximal+ j];
		switch (cur_case)
		{
			case _S_:
				alignment.push_back(pair<int,int>(i,j)); 
				i--;
				j--;
				++matches;
				break;
			case _V_:
				alignment.push_back(pair<int,int>(i,-j)); 
				i--;
				break;
			case _H_:
				alignment.push_back(pair<int,int>(-i,j)); 
				j--;
				break;
			default:
				cout << "ERROR!! -> advance_global: invalid direction D[" << k << "](" << i << ", " << j << ") = " 
				<< D[k][i*DP_maximal+ j] << endl;
				exit(-1);
		}
		cur_case=pre_case;
		count++;
	}
	while (j> 0) alignment.push_back(pair<int,int>(-i,j)),j--;
	while (i> 0) alignment.push_back(pair<int,int>(i,0)), i--;
	reverse(alignment.begin(), alignment.end());
	ali_sco=maximal;
	return matches;
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
	vector <double> WWW_score;
	WWW_score.resize(n1*n2);
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
	matchs=WWW_Advance_Align_Dyna_Prog_Double(n1,n2,WWW_score,-11,-1,-11,-1,0,0,0,0,
		WWW_alignment,sco);
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
	vector<pair<int,int> > &alignment_in,vector<pair<int,int> > &alignment_out,
	vector <int> &ali1_seq_pdb, vector <int> &ali2_seq_pdb)
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
	ali1_seq_pdb.resize(l2);
	process_oriami_record_simp(nam1_pdb.c_str(),nam1_seq.c_str(),WWW_alignment);
	for(i=0;i<l1;i++)ali1[i]=-1;
	for(i=0;i<l2;i++)ali1_seq_pdb[i]=-1;
	for(i=0;i<(int)WWW_alignment.size();i++)
	{
		int pos1=WWW_alignment[i].first;
		int pos2=WWW_alignment[i].second;
		if(pos1>0 && pos2>0)ali1[pos1-1]=pos2-1;
		if(pos1>0 && pos2>0)ali1_seq_pdb[pos2-1]=pos1-1;
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
	ali2_seq_pdb.resize(l3);
	process_oriami_record_simp(nam2_seq.c_str(),nam2_pdb.c_str(),WWW_alignment);
	for(i=0;i<l3;i++)ali3[i]=-1;
	for(i=0;i<l3;i++)ali2_seq_pdb[i]=-1;
	for(i=0;i<(int)WWW_alignment.size();i++)
	{
		int pos1=WWW_alignment[i].first;
		int pos2=WWW_alignment[i].second;
		if(pos1>0 && pos2>0)ali3[pos1-1]=pos2-1;
		if(pos1>0 && pos2>0)ali2_seq_pdb[pos1-1]=pos2-1;
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

//--------------- residue number alginment -------------//
int Residue_Number_Alignment(PDB_Residue *pdb1, PDB_Residue *pdb2,int l1,int l2,
	vector<pair<int,int> > &alignment_out)
{
	//create align
	int *ali1=new int[l1];
	int *ali2=new int[l2];
	for(int i=0;i<l1;i++)ali1[i]=-1;
	for(int i=0;i<l2;i++)ali2[i]=-1;
	//record position
	map<string, int > ws_mapping1;   //M1, mapping the PDB's name
	map<string, int > ws_mapping2;   //M2, mapping the PDB's name
	map<string, int >::iterator iter;
	//record pdb1
	for(int i=0;i<l1;i++)
	{
		//-> record
		string res1;
		if(pdb1[i].get_PDB_residue_number(res1)!=0)continue;
		iter = ws_mapping1.find(res1);
		if(iter != ws_mapping1.end())
		{
			fprintf(stderr,"duplicated mapping1 !! %s \n",res1.c_str());
			continue;
		}
		ws_mapping1.insert(map < string, int >::value_type(res1, i));
	}
	//extract pdb2
	int match=0;
	for(int i=0;i<l2;i++)
	{
		//-> record
		string res2;
		if(pdb2[i].get_PDB_residue_number(res2)!=0)continue;
		iter = ws_mapping2.find(res2);
		if(iter != ws_mapping2.end())
		{
			fprintf(stderr,"duplicated mapping2 !! %s \n",res2.c_str());
			continue;
		}
		ws_mapping2.insert(map < string, int >::value_type(res2, i));
		//-> mapping
		iter = ws_mapping1.find(res2);
		if(iter != ws_mapping1.end())
		{
			int key=ws_mapping1[res2];
			ali1[key]=i;
			ali2[i]=key;
			match++;
		}
	}
	//dump to alignment
	//final
	{
		alignment_out.clear();
		//start
		int i,j;
		int ii,jj;
		int wlen;
		int pre_ii=0;
		int pre_jj=0;
		for(i=1;i<=l1;i++)
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
		for(i=pre_jj;i<=l2;i++)alignment_out.push_back (pair<int,int>(-l1, i));  //Iy
	}
	//delete
	delete [] ali1;
	delete [] ali2;
	//return
	return match;
}




//=========== blosum calculate ============//
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
int WS_Simply_Load_PDB(string &pdbfile,XYZ *mca,XYZ *mcb,char *ami) //->from .pdb file
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
//		fprintf(stderr,"pdbfile %s not found!!\n",pdbfile.c_str());
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

//=============== Generate CB related ==============//
int *ws_reco;
XYZ *tmp_mcb;
char *tmp_ami;
void WS_Generate_CB_Simp(Confo_Beta &confo_beta,int moln,char *ami,char *cle,XYZ *mol,XYZ *mcb)
{
	int k;
	double dist;
	int ws_correct=1;
	for(k=0;k<moln;k++)
	{
		if( IsZero(mcb[k].X+99999) )
		{
			ws_correct=0;
			ws_reco[k]=0;
		}
		else
		{
			if(ami[k]!='G')ws_reco[k]=1;
			else ws_reco[k]=0;
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
		confo_beta.WS_Recon_Beta_21(mol,tmp_mcb,moln,tmp_ami,cle);
		for(k=0;k<moln;k++)if(ws_reco[k]==0)mcb[k]=tmp_mcb[k];
	}
}
void WS_Generate_CB(Confo_Beta &confo_beta,int moln,PDB_Residue *pdb,char *ami,char *cle,XYZ *mol,XYZ *mcb)
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
		confo_beta.WS_Recon_Beta_21(mol,tmp_mcb,moln,tmp_ami,cle);
		for(k=0;k<moln;k++)if(ws_reco[k]==0)mcb[k]=tmp_mcb[k];
	}
}

//============== Calc_Local_Score =============//
//[note]: get #Gap, #SeqID, LALI, blosum and clesum
int Calc_Local_Score(char *ami1,char *cle1,char *ami2,char *cle2,
	int moln1,int moln2,vector<pair<int, int> > &alignment,
	int &blos_out,int &cles_out,int &lali_out,int &seqid_out,
	int &gapo_out,int &gape_out,vector <double> &match_wei,
	vector <int> &blos_pos,vector <int> &cles_pos,vector <int> &seqid_pos)
{
	int i;
	int size=(int)alignment.size();
	//omit head gap
	int wsstart=0;
	for(i=0;i<size;i++)
	{
		int x = alignment[i].first;
		int y = alignment[i].second;
		wsstart=i;
		if(x>0 && y>0)break;
	}
	//omit tail gap
	int wsend=size-1;
	for(i=size-1;i>=0;i--)
	{
		int x = alignment[i].first;
		int y = alignment[i].second;
		wsend=i;
		if(x>0 && y>0)break;
	}
	//-> gap
	int gapo=0;
	int gape=0;
	int wlen1=0;
	int wlen2=0;
	//-> match
	int lali=0;
	int seqid=0;
	int blosum=0;
	int clesum=0;
	//count gaps and match
	int ii,jj;
	double wscur1,wscur2,wscur;
	match_wei.clear();
	blos_pos.clear();
	cles_pos.clear();
	seqid_pos.clear();
	for(i=wsstart;i<=wsend;i++)
	{
		ii=alignment[i].first;
		jj=alignment[i].second;
		//count the gap
		if(jj>0)
		{
			if(wlen1>0) gapo++,gape+=wlen1;
			wlen1 = 0;
		}
		else wlen1++;
		if(ii>0)
		{
			if(wlen2>0) gapo++,gape+=wlen2;
			wlen2 = 0;
		}
		else wlen2++;
		//count match
		if(ii>0 && jj>0)
		{
			lali++;
			if(ami1[ii-1]==ami2[jj-1])
			{
				seqid++;
				seqid_pos.push_back(1);
			}
			else seqid_pos.push_back(0);
			//calc blosum and clesum
			wscur1=BLOSUM62_Calc(ami1[ii-1],ami2[jj-1]);
			wscur2=CLESUM_Calc(cle1[ii-1],cle2[jj-1]);
			blosum+=(int)wscur1;
			clesum+=(int)wscur2;
			blos_pos.push_back((int)(wscur1));
			cles_pos.push_back((int)(0.1*wscur2));
			//add to match_wei
			if(wscur1<0)wscur1=0;      //max(BLOSUM,0)
			wscur=wscur1+0.1*wscur2;   //local parameter
			match_wei.push_back(wscur);
		}
	}
	//final process
	blos_out=blosum;
	cles_out=(int)(1.0*clesum/10);
	lali_out=lali;
	seqid_out=seqid;
	gapo_out=gapo;
	gape_out=gape;
	//return
	return lali;
}

//------ rotate full atom ------//
void Rotate_FullAtom(Kabsch &kabsch, PDB_Residue &in, PDB_Residue &out, double *rotmat)
{
	int i,k;
	int num;
	PDB_Residue pdb;
	XYZ xyz;
	int numb;
	double R,tmpr;
	//input
	pdb=in;
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
	//output
	out=pdb;
}


//------ calculate distance ---------//
//-> CA version
int Calc_Global_Score(XYZ *mol1,XYZ *mol2,int moln1,int moln2,
	double *rotmat,int *ali2,double d0,
	vector <double> &distance, vector <double> &tmsco)
{
	int i;
	int pos;
	double ori_d=d0*d0;
	double dist2;
	double dist;
	double tms;
	int lali=0;
	XYZ m1,m2,m3;
	distance.clear();
	tmsco.clear();
	for(i=0;i<moln2;i++)
	{
		pos=ali2[i];
		if(pos<0 || pos>=moln1)continue;
		//get m1,m2
		m1=mol1[pos];
		m2=mol2[i];
		m1.transform(m3,rotmat);
		dist2=m2.distance_square(m3);
		dist=sqrt(dist2);
		tms=1.0/(1.0+dist2/ori_d);
		//add
		distance.push_back(dist);
		tmsco.push_back(tms);
		lali++;
	}
	//return
	return lali;
}

//-> Full version
int Calc_Global_Score_FULL(Kabsch &kabsch,
	PDB_Residue *pdb1,PDB_Residue *pdb2,
	int moln1,int moln2,
	double *rotmat,int *ali2,double d0,
	vector <vector <string> > &recname,
	vector <vector <double> > &distance, 
	vector <vector <double> > &tmsco)
{
	int i;
	int pos;
	double ori_d=d0*d0;
	double dist2;
	double dist;
	double tms;
	int lali=0;
	PDB_Residue m1,m2,m3;
	XYZ x1,x2,x3;
	//init
	recname.clear();
	distance.clear();
	tmsco.clear();
	//proc
	for(i=0;i<moln2;i++)
	{
		pos=ali2[i];
		if(pos<0 || pos>=moln1)continue;
		//---- init ----//
		vector <string> nam;
		vector <double> dis;
		vector <double> ttm;
		//----- get m1,m2,m3 ----//
		m1=pdb1[pos];
		m2=pdb2[i];
		Rotate_FullAtom(kabsch, m1, m3, rotmat);
		//----- check ----//
		char amino1=m1.get_AA();
		char amino2=m2.get_AA();
		if(amino1 != amino2)
		{
			fprintf(stderr,"amino1 %c not equal to amino2 %c \n",amino1,amino2);
			//check add
			recname.push_back(nam);
			distance.push_back(dis);
			tmsco.push_back(ttm);
			lali++;
			continue;
		}
		char amino=amino1;
		int num;

		//-> check backbone
		num=m1.get_backbone_totnum();
		for(int k=0;k<num;k++)
		{
			//check index
			if(m2.get_backbone_part_index(k)==0)continue;
			if(m3.get_backbone_part_index(k)==0)continue;
			string atomname=backbone_atom_name_decode(k);
			//get atom
			m2.get_backbone_atom(k,x2);
			m3.get_backbone_atom(k,x3);
			//dist
			dist2=x2.distance_square(x3);
			dist=sqrt(dist2);
			tms=1.0/(1.0+dist2/ori_d);
			//add
			nam.push_back(atomname);
			dis.push_back(dist);
			ttm.push_back(tms);
		}

		//-> check sidechain
		num=m1.get_sidechain_totnum();
		for(int k=0;k<num;k++)
		{
			//check index
			if(m2.get_sidechain_part_index(k)==0)continue;
			if(m3.get_sidechain_part_index(k)==0)continue;
			string atomname=sidechain_atom_name_decode(k,amino);
			//get atom
			m2.get_sidechain_atom(k,x2);
			m3.get_sidechain_atom(k,x3);
			//dist
			dist2=x2.distance_square(x3);
			dist=sqrt(dist2);
			tms=1.0/(1.0+dist2/ori_d);
			//add
			nam.push_back(atomname);
			dis.push_back(dist);
			ttm.push_back(tms);
		}

		//--- final add ---//
		recname.push_back(nam);
		distance.push_back(dis);
		tmsco.push_back(ttm);
		lali++;
	}
	//return
	return lali;
}

//------ calculate local matrix -------//
void Calc_Local_Matrix(char *ami1,char *cle1,char *ami2,char *cle2,
	int moln1,int moln2,vector <double> &matrix_wei)
{
	int i,j;
	double wscur1,wscur2,wscur;
	matrix_wei.resize(moln1*moln2);
	for(i=0;i<moln1;i++)
	{
		int cur_index=i*moln2;
		for(j=0;j<moln2;j++)
		{
			wscur1=BLOSUM62_Calc(ami1[i],ami2[j]);
			wscur2=CLESUM_Calc(cle1[i],cle2[j]);
			//add to match_wei
			if(wscur1<0)wscur1=0;      //max(BLOSUM,0)
			wscur=wscur1+0.1*wscur2;   //local parameter
			matrix_wei[cur_index+j]=wscur;
		}
	}
}

//============== output process ===========//
//-> fasta_output
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
//-> fasta_output_simp
void FASTA_Output_Simp(FILE *fp,string &nam1,string &nam2,char *ami1,char *ami2,
	vector<pair<int,int> > &alignment)
{
	//alignment->output
	//output
	char c;
	int i;
	int ii,jj;
	int size=(int)alignment.size();
	fprintf(fp,">%s\n",nam1.c_str());
	for(i=0;i<size;i++)
	{
		ii=alignment[i].first;
		if(ii<=0)fprintf(fp,"-");
		else
		{
			c=ami1[ii-1];
			if(c=='X'||c=='Z')fprintf(fp,".");
			else fprintf(fp,"%c",c);
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
			if(c=='X'||c=='Z')fprintf(fp,".");
			else fprintf(fp,"%c",c);
		}
	}
	fprintf(fp,"\n");
}
//-> fasta_output_more
void FASTA_Output_More(string &ws_output_tot,string &nam1_,string &nam2_,
		char *ami1_,char *ami2_, char *cle1_,char *cle2_,XYZ *mol1,XYZ *mol2,double d0,
		vector<pair<int,int> > &alignment,int CLEorNOT,int threscut=80)
{
	ws_output_tot="";
	int moln1=(int)strlen(ami1_);
	int moln2=(int)strlen(ami2_);
	//-------- get start and end ---------//
	int start1,start2;
	int end1,end2;
	int wsstart=-1;
	int wsend=-1;
	int x=0;
	int y=0;
	int i;
	int size=(int)alignment.size();
	{
		//omit head gap
		for(i=0;i<size;i++)
		{
			x = alignment[i].first;
			y = alignment[i].second;
			wsstart=i;
			if(x>0 && y>0)break;
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
		char ws_output[60000];
		string ws_output_str;
		ws_output_str.clear();
		//[real process]
		int num=0;
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
					sprintf(ws_output,"cle                 %s           \n",cle1.c_str());
					ws_output_str+=ws_output;
				}
				sprintf(ws_output,"%-15s%4d %s %4d (%d)\n",
					nam1.c_str(),start1+1,ami1.c_str(),end1+1,moln1);
				ws_output_str+=ws_output;
				sprintf(ws_output,"RMSD                %s           \n",score_out.c_str());
				ws_output_str+=ws_output;
				sprintf(ws_output,"%-15s%4d %s %4d (%d)\n",
					nam2.c_str(),start2+1,ami2.c_str(),end2+1,moln2);
				ws_output_str+=ws_output;
				if(CLEorNOT==1)
				{
					sprintf(ws_output,"cle                 %s           \n",cle2.c_str());
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
			}
			else
			{
				score_out=score_out+" ";
			}
		}

		//termi_iprocess
	//	if(i%threscut==0 && first==0)
		{
			//fragment output
			ws_output_str.clear();
			sprintf(ws_output,"\n");
			ws_output_str+=ws_output;
			if(CLEorNOT==1)
			{
				sprintf(ws_output,"cle                 %s           \n",cle1.c_str());
				ws_output_str+=ws_output;
			}
			sprintf(ws_output,"%-15s%4d %s %4d (%d)\n",
				nam1.c_str(),start1+1,ami1.c_str(),end1+1,moln1);
			ws_output_str+=ws_output;
			sprintf(ws_output,"RMSD                %s           \n",score_out.c_str());
			ws_output_str+=ws_output;
			sprintf(ws_output,"%-15s%4d %s %4d (%d)\n",
				nam2.c_str(),start2+1,ami2.c_str(),end2+1,moln2);
			ws_output_str+=ws_output;
			if(CLEorNOT==1)
			{
				sprintf(ws_output,"cle                 %s           \n",cle2.c_str());
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
}
//---------------- superimpose_fullatom --------------------//
void WS_Superimpose_FullAtom(Kabsch &kabsch,PDB_Residue *in,int moln,PDB_Residue *out,double *rotmat)
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

//-------------- output detailed -------------//
//for each aligned position with BLOSUM CLESUM DeepScore RMSD TMscore.
void Output_Detailed(FILE *fp,PDB_Residue *pdb1,PDB_Residue *pdb2,
	string &nam1_content,string &nam2_content,
	vector<pair<int, int> > &alignment,
	vector <int> &ali1_seq_pdb,vector <int> &ali2_seq_pdb,
	vector <int> &blos_pos,vector <int> &cles_pos,
	vector <int> &seqid_pos,vector <double> &match_wei,
	vector <double> &distance,vector <double> &tmsco)
{
	int i;
	int ii,jj;
	int size=(int)alignment.size();
	int rel_lali=0;
	for(i=0;i<size;i++)
	{
		ii=alignment[i].first;
		jj=alignment[i].second;
		if(ii>0 && jj>0)
		{
			//check mapping align
			if(ali1_seq_pdb[ii-1]<0 || ali2_seq_pdb[jj-1]<0)
			{
				fprintf(fp,"%c%c X\n",nam1_content[ii-1],nam2_content[jj-1]);
			}
			else
			{
				//-> get residue number
				string res1="";
				if(pdb1[ii-1].get_PDB_residue_number(res1)==0)res1=res1.substr(1,4);
				string res2="";
				if(pdb2[jj-1].get_PDB_residue_number(res2)==0)res2=res2.substr(1,4);
				//-> output
				fprintf(fp,"%c%c %4s %4s | %1d %4d %4d %6.1f | %6.2f %6.3f\n",
					nam1_content[ii-1],nam2_content[jj-1],res1.c_str(),res2.c_str(),
					seqid_pos[rel_lali],blos_pos[rel_lali],cles_pos[rel_lali],match_wei[rel_lali],
					distance[rel_lali],tmsco[rel_lali]);
				rel_lali++;
			}
		}
		else
		{
			if(ii>0)
			{
				fprintf(fp,"%c-\n",nam1_content[ii-1]);
			}
			if(jj>0)
			{
				fprintf(fp,"-%c\n",nam2_content[jj-1]);
			}
		}
	}
}

//-> Full version
void Output_Detailed_FULL(FILE *fp,PDB_Residue *pdb1,PDB_Residue *pdb2,
	string &nam1_content,string &nam2_content,
	vector<pair<int, int> > &alignment,
	vector <int> &ali1_seq_pdb,vector <int> &ali2_seq_pdb,
	vector <int> &blos_pos,vector <int> &cles_pos,
	vector <int> &seqid_pos,vector <double> &match_wei,
	vector <vector <string> > &recname,
	vector <vector <double> > &distance,
	vector <vector <double> > &tmsco)
{
	int i;
	int ii,jj;
	int size=(int)alignment.size();
	int rel_lali=0;
	for(i=0;i<size;i++)
	{
		ii=alignment[i].first;
		jj=alignment[i].second;
		if(ii>0 && jj>0)
		{
			//check mapping align
			if(ali1_seq_pdb[ii-1]<0 || ali2_seq_pdb[jj-1]<0)
			{
				fprintf(fp,"%c%c  CA  X\n",nam1_content[ii-1],nam2_content[jj-1]);
			}
			else
			{
				//---- FULL version ----//
				int num=(int)recname[rel_lali].size();
				for(int k=0;k<num;k++)
				{
					//-> get residue number
					string res1="";
					if(pdb1[ii-1].get_PDB_residue_number(res1)==0)res1=res1.substr(1,4);
					string res2="";
					if(pdb2[jj-1].get_PDB_residue_number(res2)==0)res2=res2.substr(1,4);
					//-> output
					fprintf(fp,"%c%c %4s %4s %4s | %1d %4d %4d %6.1f | %6.2f %6.3f\n",
						nam1_content[ii-1],nam2_content[jj-1],recname[rel_lali][k].c_str(),res1.c_str(),res2.c_str(),
						seqid_pos[rel_lali],blos_pos[rel_lali],cles_pos[rel_lali],match_wei[rel_lali],
						distance[rel_lali][k],tmsco[rel_lali][k]);
				}
				rel_lali++;
			}
		}
		else
		{
			if(ii>0)
			{
				fprintf(fp,"%c-  CA\n",nam1_content[ii-1]);
			}
			if(jj>0)
			{
				fprintf(fp,"-%c  CA\n",nam2_content[jj-1]);
			}
		}
	}
}


//===================== script related =======================//__2014_05_30__//
//----------------- Alignment_To_Script ------------//
//[note]: the size of AFP should be (moln1+moln2)*4
int Alignment_To_AFP(vector<pair<int, int> > &alignment,
	int moln1,int moln2,int *AFP_Cor)
{
	int i;
	int num;
	int ii,jj;
	int count;
	int isFirst;
	int isLast;
	int head1,head2;
	int index;
	int cur_ii;
	int INT_MIN_NUM_ = -99999;

	//alignment to ali2
	int *ali1=new int[moln1];
	int *ali2=new int[moln2];
	for(i=0;i<moln1;i++)ali1[i]=-1;
	for(i=0;i<moln2;i++)ali2[i]=-1;
	//--- extract alignment ---//
	int totnum=(int)alignment.size();
	for(i=0;i<totnum;i++)
	{
		ii=alignment[i].first;
		jj=alignment[i].second;
		if(ii>0&&jj>0)
		{
			ali1[ii-1]=jj-1;
			ali2[jj-1]=ii-1;
		}
	}

	//--- final return ----//
	//init
	AFP_Cor[0]=0;
	num=0;
	ii=-1;
	jj=-1;
	cur_ii=-1;
	head1=-1;
	head2=-1;
	isLast=0;
	isFirst=1;
	count=INT_MIN_NUM_;
	for(i=0;i<moln2;i++)
	{
		cur_ii=ali2[i];
		if(cur_ii==-1) //purely blank
		{
			if(isFirst==0)
			{
				if(count>0)
				{
					AFP_Cor[0]++;
					index=AFP_Cor[0];
					AFP_Cor[4*index+0] = 0;
					AFP_Cor[4*index+1] = head1;
					AFP_Cor[4*index+2] = head2;
					AFP_Cor[4*index+3] = count;
					num+=count;
				}
				count=0;
				isFirst=1;
			}
			continue;
		}
		if(isFirst==1)
		{
ws_init:
			isFirst=0;
			ii=cur_ii;
			jj=i;
			count=1;
			head1=ii;
			head2=jj;
			continue;
		}
		if(i==jj+1&&cur_ii==ii+1)
		{
			ii=cur_ii;
			jj=i;
			count++;
			continue;
		}

ws_end:
		if(count>0)
		{
			AFP_Cor[0]++;
			index=AFP_Cor[0];
			AFP_Cor[4*index+0] = 0;
			AFP_Cor[4*index+1] = head1;
			AFP_Cor[4*index+2] = head2;
			AFP_Cor[4*index+3] = count;
			num+=count;
		}
		if(isLast==1)goto end;
		else goto ws_init;
	}
	if(count==INT_MIN_NUM_)goto end;
	isLast=1;
	goto ws_end;
end:
	delete [] ali1;
	delete [] ali2;
	return num;
}


//=============== AFP -> Script ============//__130830__//
//-> JMol
void Output_JMol_Script(FILE *fws,int *AFP_Cor)
{
	//process
	int j,k;
	int jj;
	int totnum;
	int winlen;
//	fprintf(fws,"load %s\n",name.c_str());
	fprintf(fws,"rotate x, 180\n");
	fprintf(fws,"backbone only\n");
	fprintf(fws,"backbone 10\n");
//	fprintf(fws,"set background white\n");
//	fprintf(fws,"select :A\n");
//	fprintf(fws,"color blue\n");
//	fprintf(fws,"select :B\n");
//	fprintf(fws,"color red\n");
	totnum=AFP_Cor[0];
	for(k=1;k<=totnum;k++)
	{
		for(j=0;j<2;j++)
		{
			//record
			jj=AFP_Cor[k*4+j+1];
			if(jj==-1)continue;
			winlen=AFP_Cor[k*4+3];
			fprintf(fws,"select   %d-%d:%c\n", jj+1,jj+winlen,j+'A');
			fprintf(fws,"backbone off\n");
			fprintf(fws,"cartoon\n");
		}
	}
	fprintf(fws,"select all\n");
	fprintf(fws,"color chain\n");
}
//-> RasMol
void Output_RasMol_Script(FILE *fws,int *AFP_Cor)
{
	//process
	int j,k;
	int jj;
	int totnum;
	int winlen;
//	fprintf(fws,"load %s\n",name.c_str());
	fprintf(fws,"wireframe off\n");
	fprintf(fws,"backbone 10\n");
	fprintf(fws,"set ambient 20\n");
//	fprintf(fws,"set background white\n");
//	fprintf(fws,"select :A\n");
//	fprintf(fws,"color blue\n");
//	fprintf(fws,"select :B\n");
//	fprintf(fws,"color red\n");
	totnum=AFP_Cor[0];
	for(k=1;k<=totnum;k++)
	{
		for(j=0;j<2;j++)
		{
			//record
			jj=AFP_Cor[k*4+j+1];
			if(jj==-1)continue;
			winlen=AFP_Cor[k*4+3];
			fprintf(fws,"select   %d-%d:%c\n", jj+1,jj+winlen,j+'A');
			fprintf(fws,"backbone off\n");
			fprintf(fws,"cartoon\n");
		}
	}
	fprintf(fws,"select all\n");
	fprintf(fws,"color chain\n");
}
//-> PyMol
void Output_PyMol_Script(FILE *fws,int *AFP_Cor)
{
	//process
	int k;
	int jj;
	int totnum;
	int winlen;
//	fprintf(fws,"load %s\n",name.c_str());
	fprintf(fws,"bg_color white\n");
	fprintf(fws,"rotate x, 180\n");
	fprintf(fws,"hide everything, all\n");
	fprintf(fws,"set ribbon_width, 1.8\n");
	totnum=AFP_Cor[0];
	string wsa="select aa, chain A and ";
	string wsb="select bb, chain B and ";
	string wsan="select aan, chain A and not ";
	string wsbn="select bbn, chain B and not ";
	string temp;
	char command[30000];
	if(totnum>1)
	{
		for(k=1;k<=totnum;k++)
		{
			//record A
			jj=AFP_Cor[k*4+0+1];
			if(jj==-1)continue;
			winlen=AFP_Cor[k*4+3];
			//-> for aa
			sprintf(command,"resi %d-%d ", jj+1,jj+winlen);
			temp=command;
			if(k==1)wsa=wsa+"( "+temp;
			else if(k==totnum)wsa=wsa+"or "+temp+" )";
			else wsa=wsa+"or "+temp;
			//-> for aan
			sprintf(command,"resi %d-%d ", jj+1+1,jj+winlen-1);
			temp=command;
			if(k==1)wsan=wsan+"( "+temp;
			else if(k==totnum)wsan=wsan+"or "+temp+" )";
			else wsan=wsan+"or "+temp;
			//record B
			jj=AFP_Cor[k*4+1+1];
			if(jj==-1)continue;
			winlen=AFP_Cor[k*4+3];
			//-> for bb
			sprintf(command,"resi %d-%d ", jj+1,jj+winlen);
			temp=command;
			if(k==1)wsb=wsb+"( "+temp;
			else if(k==totnum)wsb=wsb+"or "+temp+" )";
			else wsb=wsb+"or "+temp;
			//-> for bbn
			sprintf(command,"resi %d-%d ", jj+1+1,jj+winlen-1);
			temp=command;
			if(k==1)wsbn=wsbn+"( "+temp;
			else if(k==totnum)wsbn=wsbn+"or "+temp+" )";
			else wsbn=wsbn+"or "+temp;
		}
	}
	else
	{
		k=1;
		//record A
		jj=AFP_Cor[k*4+0+1];
		winlen=AFP_Cor[k*4+3];
		//-> for aa
		sprintf(command,"resi %d-%d ", jj+1,jj+winlen);
		temp=command;
		wsa=wsa+"( "+temp + " )";
		//-> for aan
		sprintf(command,"resi %d-%d ", jj+1+1,jj+winlen-1);
		temp=command;
		wsan=wsan+"( "+temp + " )";
		//record B
		jj=AFP_Cor[k*4+1+1];
		winlen=AFP_Cor[k*4+3];
		//-> for bb
		sprintf(command,"resi %d-%d ", jj+1,jj+winlen);
		temp=command;
		wsb=wsb+"( "+temp + " )";
		//-> for bbn
		sprintf(command,"resi %d-%d ", jj+1+1,jj+winlen-1);
		temp=command;
		wsbn=wsbn+"( "+temp + " )";
	}
	//-> aa
	fprintf(fws,"%s\n",wsa.c_str());
	fprintf(fws,"color blue, aa\n");
	fprintf(fws,"hide lines, aa\n");
	fprintf(fws,"show cartoon, aa\n");
	//-> aan
	fprintf(fws,"%s\n",wsan.c_str());
	fprintf(fws,"color blue, aan\n");
	fprintf(fws,"hide lines, aan\n");
	fprintf(fws,"show ribbon, aan\n");
	//-> bb
	fprintf(fws,"%s\n",wsb.c_str());
	fprintf(fws,"color red, bb\n");
	fprintf(fws,"hide lines, bb\n");
	fprintf(fws,"show cartoon, bb\n");
	//-> bbn
	fprintf(fws,"%s\n",wsbn.c_str());
	fprintf(fws,"color red, bbn\n");
	fprintf(fws,"hide lines, bbn\n");
	fprintf(fws,"show ribbon, bbn\n");
}



//============== main process =============//
void Main_Process(string &file1,string &range1,string &file2,string &range2,
	string &ali_file,int &ali_res,string &out_file,string &detail_file,int &detail_full,
	int Out_Script,int Normalize,double Distance_Cutoff,int Out_Screen,int SIMPLY_LOAD,int TM_type)
{
	//data structure
	Confo_Lett confo_lett;
	XYZ *TM_MOL1;
	XYZ *TM_MOL_TMP;
	XYZ *TM_MOL2;
	XYZ *TM_MCB1;
	XYZ *TM_MCB2;
	PDB_Residue *TM_PDB1;
	PDB_Residue *TM_PDB_TMP;
	PDB_Residue *TM_PDB2;
	char *TM_AMI1;
	char *TM_AMI2;
	char *TM_CLE1;
	char *TM_CLE2;
	int TM_MOLN1;
	int TM_MOLN2;

	//load molecular
	//-> check length
	TM_MOLN1=Get_PDB_File_Len(file1);
	if(TM_MOLN1==0)
	{
		fprintf(stderr,"FILE1 BLANK [%s]\n",file1.c_str());
		exit(-1);
	}
	TM_MOLN2=Get_PDB_File_Len(file2);
	if(TM_MOLN2==0)
	{
		fprintf(stderr,"FILE2 BLANK [%s]\n",file2.c_str());
		exit(-1);
	}
	//-> get name
	string nam1,nam2;
	getBaseName(file1,nam1,'/','.');
	getBaseName(file2,nam2,'/','.');
	//-> create data
	//--| mol1
	TM_MOL1=new XYZ[TM_MOLN1];
	TM_MOL_TMP=new XYZ[TM_MOLN1];
	TM_MCB1=new XYZ[TM_MOLN1];
	TM_PDB1=new PDB_Residue[TM_MOLN1];
	TM_PDB_TMP=new PDB_Residue[TM_MOLN1];
	TM_AMI1=new char[TM_MOLN1+1];
	TM_CLE1=new char[TM_MOLN1+1];
	//--| mol2
	TM_MOL2=new XYZ[TM_MOLN2];
	TM_MCB2=new XYZ[TM_MOLN2];
	TM_PDB2=new PDB_Residue[TM_MOLN2];
	TM_AMI2=new char[TM_MOLN2+1];
	TM_CLE2=new char[TM_MOLN2+1];
	//->create temp
	int DEEPALIGN_MAXSIZE=TM_MOLN1>TM_MOLN2?TM_MOLN1:TM_MOLN2;
	ws_reco=new int[DEEPALIGN_MAXSIZE];
	tmp_mcb=new XYZ[DEEPALIGN_MAXSIZE];
	tmp_ami=new char[DEEPALIGN_MAXSIZE+1];
	//-> load file
	if(SIMPLY_LOAD==0)
	{
		mol_input.XYZ_Input(file1,range1,0,TM_MOLN1,TM_MOL1,TM_AMI1,TM_CLE1,0,TM_PDB1);
		mol_input.XYZ_Input(file2,range2,0,TM_MOLN2,TM_MOL2,TM_AMI2,TM_CLE2,0,TM_PDB2);
	}
	else
	{
		WS_Simply_Load_PDB(file1,TM_MOL1,TM_MCB1,TM_AMI1);
		WS_Simply_Load_PDB(file2,TM_MOL2,TM_MCB2,TM_AMI2);
	}

	//construct c_beta
	Confo_Beta confo_beta(DEEPALIGN_MAXSIZE);
	if(SIMPLY_LOAD==0)
	{
		WS_Generate_CB(confo_beta,TM_MOLN1,TM_PDB1,TM_AMI1,TM_CLE1,TM_MOL1,TM_MCB1);
		WS_Generate_CB(confo_beta,TM_MOLN2,TM_PDB2,TM_AMI2,TM_CLE2,TM_MOL2,TM_MCB2);
	}
	else
	{
		confo_lett.btb_ori(0,0,0,TM_MOLN1,TM_MOL1,TM_CLE1);
		TM_CLE1[TM_MOLN1]='\0';
		WS_Generate_CB_Simp(confo_beta,TM_MOLN1,TM_AMI1,TM_CLE1,TM_MOL1,TM_MCB1);
		confo_lett.btb_ori(0,0,0,TM_MOLN2,TM_MOL2,TM_CLE2);
		TM_CLE1[TM_MOLN1]='\0';
		WS_Generate_CB_Simp(confo_beta,TM_MOLN2,TM_AMI2,TM_CLE2,TM_MOL2,TM_MCB2);
	}

	//alignment process
	//-> load alignment
	vector<pair<int, int> > alignment;
	vector<pair<int, int> > alignment_out;
	string nam1_con,nam2_con;
	vector <int> ali1_seq_pdb;
	vector <int> ali2_seq_pdb;
	if(ali_file!="")  //-> load alignment file
	{
		//load alignment file
		string nam1_content,nam2_content;
		string nam1_full,nam2_full;
		string wnam1,wnam2;
		int retv=ReadToFile_FASTA(ali_file,alignment,nam1_content,nam2_content,nam1_full,nam1_full,wnam1,wnam2);
		if(retv==-1)exit(-1);
		nam1_con=nam1_content;
		nam2_con=nam2_content;
		//extract alignment
		string nam1_pdb=TM_AMI1;
		string nam2_pdb=TM_AMI2;
		Extract_Alignment(nam1_content,nam1_pdb,nam2_content,nam2_pdb,alignment,alignment_out,
			ali1_seq_pdb,ali2_seq_pdb);
	}
	else
	{
		//init
		string nam1_pdb=TM_AMI1;
		string nam2_pdb=TM_AMI2;
		nam1_con=nam1_pdb;
		nam2_con=nam2_pdb;
		//align
		if(ali_res==1)  //-> use residue number
			Residue_Number_Alignment(TM_PDB1,TM_PDB2,TM_MOLN1,TM_MOLN2,alignment_out);
		else            //-> generate alignment
			process_oriami_record_simp(nam1_pdb.c_str(),nam2_pdb.c_str(),alignment_out);
		alignment=alignment_out;
		//assign
		ali1_seq_pdb.resize(nam1_pdb.length());
		for(int i=0;i<(int)nam1_pdb.length();i++)ali1_seq_pdb[i]=i;
		ali2_seq_pdb.resize(nam2_pdb.length());
		for(int i=0;i<(int)nam2_pdb.length();i++)ali2_seq_pdb[i]=i;
	}


	//-> ali1 & ali2
	int *ali1=new int[TM_MOLN1];
	int *ali2=new int[TM_MOLN2];
	for(int i=0;i<TM_MOLN1;i++)ali1[i]=-1;
	for(int i=0;i<TM_MOLN2;i++)ali2[i]=-1;
	for(int i=0;i<(int)alignment_out.size();i++)
	{
		int ii=alignment_out[i].first;
		int jj=alignment_out[i].second;
		if(ii>0 && jj>0)
		{
			ali1[ii-1]=jj-1;
			ali2[jj-1]=ii-1;
		}
	}

	//calculate global and local score
	//-> get weight matrix	
	TM_align tm_align(DEEPALIGN_MAXSIZE);
	tm_align.TM_Align_Init(TM_MOLN1,TM_MOLN2);
	if(Distance_Cutoff==0)
	{
		tm_align.TM_DIST_CUT=0;
	}
	else if(Distance_Cutoff<0)
	{
		tm_align.TM_DIST_CUT=1;
		Distance_Cutoff=tm_align.d8;
	}
	else
	{
		tm_align.TM_DIST_CUT=1;
		tm_align.d8=Distance_Cutoff;
	}
	tm_align.TM_cb1=TM_MCB1;
	tm_align.TM_cb2=TM_MCB2;
	//-> calculate local score
	int blos_out,cles_out;
	int lali_out,seqid_out;
	int gapo_out,gape_out;
	vector <int> blos_pos;
	vector <int> cles_pos;
	vector <int> seqid_pos;
	vector <double> match_wei;
	Calc_Local_Score(TM_AMI1,TM_CLE1,TM_AMI2,TM_CLE2,TM_MOLN1,TM_MOLN2,alignment_out,
		blos_out,cles_out,lali_out,seqid_out,gapo_out,gape_out,match_wei,blos_pos,cles_pos,seqid_pos);
	//normalize //----------------> need to be rewritten !! //__110630__//
	int norm_len=TM_MOLN2; //normal, using the second to normalize
	double norm_d0=tm_align.d0;
	{
		//----- Jinbo_Case -----//
		int minres,maxres;
		if(Normalize==-4)norm_len=Extract_MinMax_ResNum(file2,minres,maxres);
		else if(Normalize==-3)norm_len=Extract_MinMax_ResNum(file1,minres,maxres);
		else if(Normalize==-2)norm_len=TM_MOLN2;
		else if(Normalize==-1)norm_len=TM_MOLN1;
		else if(Normalize==0)norm_len=min(TM_MOLN1,TM_MOLN2);
		else if(Normalize>0)norm_len=Normalize;
		else
		{
			fprintf(stderr,"WARNING: incorrect normalize length [%d]!! using min as default !!\n",Normalize);
			norm_len=min(TM_MOLN1,TM_MOLN2);
		}
		norm_d0=tm_align.Calc_TM_d0_Simp(norm_len);
	}
	//-> calculate global score
	int lali;
	double rms;
	double tms;
	double Ret_Sco[8];
	double TM_ROTMAT[12];
	//different types of TMscore calculation
	if(TM_type==1) //typical type
		tms=tm_align.TM_Align_TM_Score(TM_MOL1,TM_MOL2,TM_MOLN1,TM_MOLN2,ali2,norm_len,norm_d0,rms,lali,Ret_Sco);  //tm-score (typical)
	else if(TM_type==2) //simplified type
		tms=tm_align.TM_Align_TM_Score_Simp(TM_MOL1,TM_MOL2,TM_MOLN1,TM_MOLN2,ali2,norm_len,norm_d0,rms,lali,Ret_Sco);  //tm-score (simplified)
	else
		tms=tm_align.TM_Align_TM_Score_Simplest(TM_MOL1,TM_MOL2,TM_MOLN1,TM_MOLN2,ali2,norm_len,norm_d0,rms,lali,Ret_Sco);  //tm-score (direct superimpose)
	for(int i=0;i<12;i++)TM_ROTMAT[i]=tm_align.finmat[i];
	double wms=tm_align.TM_Align_Get_Score_Simp_MatchWei(TM_MOL1,TM_MOL2,TM_ROTMAT,TM_MOLN1,TM_MOLN2,ali2,match_wei); //DeepAlign-score
	double gdt1=1.0*(Ret_Sco[1]+Ret_Sco[2]+Ret_Sco[3]+Ret_Sco[4])/(4.0*norm_len);  //ori_GDT
	double gdt2=1.0*(Ret_Sco[0]+Ret_Sco[1]+Ret_Sco[2]+Ret_Sco[3])/(4.0*norm_len);  //ha_GDT
	double ugdt=0.25*(Ret_Sco[1]+Ret_Sco[2]+Ret_Sco[3]+Ret_Sco[4]);                //uGDT
	double maxsub=1.0*Ret_Sco[5]/norm_len;                                         //maxsub

	//-> check lali
	if( Distance_Cutoff==0 &&  lali_out!=lali)
	{
		fprintf(stderr,"LALI NOT EQUAL [%d]!=[%d] \n",lali_out,lali);
		exit(-1);
	}

	//output superimposed structure
	FILE *fp;
	string ws_output_tot;
	Kabsch kabsch;
	if(out_file!="")
	{
		char www_nam[30000];
		string TER="TER                                                                             ";
		string END="END                                                                             ";
		string outnam=out_file;
		//[1]ami_fasta_out
		sprintf(www_nam,"%s.fasta",outnam.c_str());
		fp=fopen(www_nam,"wb");
		if(fp==0)
		{
			fprintf(stderr,"ERROR: file %s can't be opened. \n",www_nam);
		}
		else
		{
			FASTA_Output_Simp(fp,nam1,nam2,TM_AMI1,TM_AMI2,alignment_out);
			fclose(fp);
		}
		//[2]cle_fasta_out
		//sprintf(www_nam,"%s.fasta_cle",outnam.c_str());
		//fp=fopen(www_nam,"wb");
		//FASTA_Output_Simp(fp,nam1,nam2,TM_CLE1,TM_CLE2,alignment_out);
		//fclose(fp);
		//[3]pdb_out
		kabsch.rot_mol(TM_MOL1,TM_MOL_TMP,TM_MOLN1,TM_ROTMAT);
		sprintf(www_nam,"%s.pdb",outnam.c_str());
		fp=fopen(www_nam,"wb");
		if(fp==0)
		{
			fprintf(stderr,"ERROR: file %s can't be opened. \n",www_nam);
		}
		else
		{
			if(SIMPLY_LOAD==0)  //-> output Full_Atom PDB
			{
				WS_Superimpose_FullAtom(kabsch,TM_PDB1,TM_MOLN1,TM_PDB_TMP,TM_ROTMAT);
				mol_out.Output_PDB_III(fp,TM_MOLN1,TM_PDB_TMP,'A',1);
				fprintf(fp,"%s\n",TER.c_str());
				mol_out.Output_PDB_III(fp,TM_MOLN2,TM_PDB2,'B',1);
				fprintf(fp,"%s\n",TER.c_str());
				fprintf(fp,"%s\n",END.c_str());
				//--> output for linear chain
				if(Out_PDB==1)
				{
					// output file
					sprintf(www_nam,"%s.pdb_linear",outnam.c_str());
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
		//[4]local_out
		sprintf(www_nam,"%s.local",outnam.c_str());
		fp=fopen(www_nam,"wb");
		if(fp==0)
		{
			fprintf(stderr,"ERROR: file %s can't be opened. \n",www_nam);
		}
		else
		{
			FASTA_Output_More(ws_output_tot,nam1,nam2,TM_AMI1,TM_AMI2,TM_CLE1,TM_CLE2,TM_MOL_TMP,TM_MOL2,norm_d0,alignment_out,0);
			fprintf(fp,"%s",ws_output_tot.c_str());
			fclose(fp);
		}
		//[5]score_out
		sprintf(www_nam,"%s.score",outnam.c_str());
		fp=fopen(www_nam,"wb");
		if(fp==0)
		{
			fprintf(stderr,"ERROR: file %s can't be opened. \n",www_nam);
		}
		else
		{
			fprintf(fp,"#name1 name2 len1 len2 -> BLOSUM CLESUM DeepScore -> LALI RMSDval TMscore -> MAXSUB GDT_TS GDT_HA -> SeqID nLen dCutoff uGDT\n");
			fprintf(fp," %s %s %4d %4d -> %6d %6d %9.2f -> %4d %7.3f %7.3f -> %6.3f %6.3f %6.3f -> %5d %5d %6.3f %.3f\n",
				nam1.c_str(),nam2.c_str(),TM_MOLN1,TM_MOLN2,blos_out,cles_out,wms,lali,rms,tms,maxsub,gdt1,gdt2,seqid_out,norm_len,Distance_Cutoff,ugdt);
			fprintf(fp,"#---------------- transformation to superpose 1st structure onto the 2nd ------------------------\n");
			fprintf(fp," %9.6f %9.6f %9.6f %12.6f\n",TM_ROTMAT[0],TM_ROTMAT[1],TM_ROTMAT[2],TM_ROTMAT[9]);
			fprintf(fp," %9.6f %9.6f %9.6f %12.6f\n",TM_ROTMAT[3],TM_ROTMAT[4],TM_ROTMAT[5],TM_ROTMAT[10]);
			fprintf(fp," %9.6f %9.6f %9.6f %12.6f\n",TM_ROTMAT[6],TM_ROTMAT[7],TM_ROTMAT[8],TM_ROTMAT[11]);
			fclose(fp);
		}
		//[6]matrix_out
		//vector <double> wei_matrix;
		//Calc_Local_Matrix(TM_AMI1,TM_CLE1,TM_AMI2,TM_CLE2,TM_MOLN1,TM_MOLN2,wei_matrix);
		//vector <double> ret_matrix;
		//tm_align.TM_Align_Get_Matrix_MatchWei(TM_MOL1,TM_MOL2,TM_MOLN1,TM_MOLN2,norm_d0,wei_matrix,ret_matrix);
		//sprintf(www_nam,"%s.matrix",outnam.c_str());
		//ofstream fout(www_nam);
		//for(int i=0;i<TM_MOLN1;i++)for(int j=0;j<TM_MOLN2;j++)fout << ret_matrix[i*TM_MOLN2+j] << " ";
		//fout.close();

		//[7]script_out
		if(Out_Script>0) //-> output script
		{
			//generate AFP
			int *AFP=new int[(TM_MOLN1+TM_MOLN2)*4];
			Alignment_To_AFP(alignment_out,TM_MOLN1,TM_MOLN2,AFP);
			//output script
			sprintf(www_nam,"%s.scp",outnam.c_str());
			fp=fopen(www_nam,"wb");
			if(fp==0)
			{
				fprintf(stderr,"ERROR: file %s can't be opened. \n",www_nam);
			}
			else
			{
				if(Out_Script==1)Output_JMol_Script(fp,AFP);
				if(Out_Script==2)Output_RasMol_Script(fp,AFP);
				if(Out_Script==3)Output_PyMol_Script(fp,AFP);
				fclose(fp);
			}
			delete [] AFP;
		}

	}
	else
	{
		kabsch.rot_mol(TM_MOL1,TM_MOL_TMP,TM_MOLN1,TM_ROTMAT);
		FASTA_Output_More(ws_output_tot,nam1,nam2,TM_AMI1,TM_AMI2,TM_CLE1,TM_CLE2,TM_MOL_TMP,TM_MOL2,norm_d0,alignment_out,0);
	}

	//output detailed score for each aligned position
	if(detail_file!="")
	{
		fp=fopen(detail_file.c_str(),"wb");
		if(fp==0)
		{
			fprintf(stderr,"ERROR: file %s can't be opened. \n",detail_file.c_str());
		}
		else
		{
			if(detail_full==0)
			{
				vector <double> distance;
				vector <double> tmsco;
				Calc_Global_Score(TM_MOL1,TM_MOL2,TM_MOLN1,TM_MOLN2,TM_ROTMAT,ali2,norm_d0,distance,tmsco);
				Output_Detailed(fp,TM_PDB1,TM_PDB2,nam1_con,nam2_con,alignment,ali1_seq_pdb,ali2_seq_pdb,
					blos_pos,cles_pos,seqid_pos,match_wei,distance,tmsco);
			}
			else
			{
				vector <vector <string> > recname;
				vector <vector <double> > distance;
				vector <vector <double> > tmsco;
				Calc_Global_Score_FULL(kabsch,TM_PDB1,TM_PDB2,TM_MOLN1,TM_MOLN2,TM_ROTMAT,ali2,norm_d0,recname,distance,tmsco);
				Output_Detailed_FULL(fp,TM_PDB1,TM_PDB2,nam1_con,nam2_con,alignment,ali1_seq_pdb,ali2_seq_pdb,
					blos_pos,cles_pos,seqid_pos,match_wei,recname,distance,tmsco);
			}
			fclose(fp);
		}
	}

	//output score
	//[note]: name according to input files
//	printf("%s %s %4d %4d -> %6d %6d %9.2f -> %4d %7.3f %7.3f -> %6.3f %6.3f %6.3f -> %4d %4d %4d\n",
//		nam1.c_str(),nam2.c_str(),TM_MOLN1,TM_MOLN2,blos_out,cles_out,wms,lali,rms,tms,maxsub,gdt1,gdt2,gapo_out,gape_out,seqid_out);
//	printf("%s %s %4d %4d -> %6d %6d %9.2f -> %4d %7.3f %7.3f -> %6.3f %6.3f %6.3f\n",
//		nam1.c_str(),nam2.c_str(),TM_MOLN1,TM_MOLN2,blos_out,cles_out,wms,lali,rms,tms,maxsub,gdt1,gdt2);

	//--- print out ----//
	char ws_command[300000];
	string out_str="";
	if(Out_Screen==0)  //simplest screen-out
	{
		sprintf(ws_command,"%s %s %4d %4d -> %6d %6d %9.2f -> %4d %7.3f %7.3f -> %6.3f %6.3f %6.3f -> %5d %5d %6.3f %.3f\n",
			nam1.c_str(),nam2.c_str(),TM_MOLN1,TM_MOLN2,blos_out,cles_out,wms,lali,rms,tms,maxsub,gdt1,gdt2,seqid_out,norm_len,Distance_Cutoff,ugdt);
		out_str=ws_command;
	}
	else //-> screen out new
	{
		string fin_tmp_str="";
		//--- header out ---//
		sprintf(ws_command,"#------------------------------------------------------#\n");
		fin_tmp_str+=ws_command;
		sprintf(ws_command,"#                    DeepScore                         #\n");
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
		sprintf(ws_command,"1st input protein: %s    length= %4d \n",nam1.c_str(),TM_MOLN1);
		fin_tmp_str+=ws_command;
		sprintf(ws_command,"2nd input protein: %s    length= %4d \n",nam2.c_str(),TM_MOLN2);
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
			blos_out,cles_out,wms,seqid_out,lali,rms,tms,maxsub,gdt1,gdt2,ugdt);
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
		sprintf(ws_command,"%s",ws_output_tot.c_str());
		fin_tmp_str+=ws_command;
		//-- final return --//
		out_str=fin_tmp_str;
	}

	//screen out
	printf("%s",out_str.c_str());


	//delete
	delete [] TM_MOL1;
	delete [] TM_MCB1;
	delete [] TM_PDB1;
	delete [] TM_PDB_TMP;
	delete [] TM_AMI1;
	delete [] TM_CLE1;
	delete [] TM_MOL2;
	delete [] TM_MCB2;
	delete [] TM_PDB2;
	delete [] TM_AMI2;
	delete [] TM_CLE2;
	delete [] ali1;
	delete [] ali2;
	//delete temp
	delete [] ws_reco;
	delete [] tmp_mcb;
	delete [] tmp_ami;
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


//------- usage ------//
void Usage(void)
{
/*
	fprintf(stderr,"DeepScore V1.02 [Dec-30-2013] \n");
	fprintf(stderr,"---------------------- REFERENCE ----------------------\n");
	fprintf(stderr,"Sheng Wang, Jianzhu Ma, Jian Peng and Jinbo Xu.\n");
	fprintf(stderr,"   PROTEIN STRUCTURE ALIGNMENT BEYOND SPATIAL PROXIMITY\n");
	fprintf(stderr,"                Scientific Reports, 3, 1448, (2013) \n");
	fprintf(stderr,"====================== USAGE: =========================\n");
	fprintf(stderr,"./DeepScore <-t protein_A> <-q protein_B> [-o out_name] \n");
	fprintf(stderr,"   [-a alignment] [-x range_A] [-y range_B] [-n normal_len]\n");
	fprintf(stderr,"Example: ./DeepScore -t 1col.pdb -q 1cpc.pdb -o testout \n");
	fprintf(stderr,"    -a 1col_1cpc.fasta -x A:5-121 -y L:1-132  \n");
	fprintf(stderr,"[range_A or range_B] -> the optional residue range. \n"); 
	fprintf(stderr,"if not specified, all residues in the 1st chain are used.\n");
	fprintf(stderr,"[alignment] -> the alignment file between input proteins.\n");
	fprintf(stderr,"if not specified, the sequence alignment will be used. \n");
	fprintf(stderr,"[out_name] -> if specified, then output evaluation files. \n");
	fprintf(stderr,"    otherwise, only screenout evaluation scores. \n");
	fprintf(stderr,"[normal_len] -> score normalized by the following length,\n");
	fprintf(stderr,"   -2: second, -1: first, [0]: min, or other given length\n");
	fprintf(stderr,"---------------------- SCREENOUT ----------------------\n");
	fprintf(stderr,"name1 name2 len1 len2 -> BLOSUM CLESUM DeepScore -> \n");
	fprintf(stderr,"         LALI RMSDval TMscore -> MAXSUB GDT_TS GDT_HA\n");
	fprintf(stderr,"-------------------------------------------------------\n");
*/

	fprintf(stderr,"DeepScore v1.12 [Aug-11-2018] \n");
	fprintf(stderr,"Sheng Wang, Jianzhu Ma, Jian Peng and Jinbo Xu.\n");
	fprintf(stderr,"   PROTEIN STRUCTURE ALIGNMENT BEYOND SPATIAL PROXIMITY\n");
	fprintf(stderr,"                Scientific Reports, 3, 1448, (2013) \n\n");
	fprintf(stderr,"Usage: \n");
	fprintf(stderr,"./DeepScore protein_1 protein_2 [-x range_1] [-y range_2] [-a alignment] [-o out_name] \n");
	fprintf(stderr,"     [-d/D detail_file] [-s script_option] [-P screenout] [-n normalize_len] [-C distance_cut] [-T tm_type] \n\n");
	fprintf(stderr,"Required input: \n");
	fprintf(stderr," protein_1:             The 1st input protein file in PDB format. \n");
	fprintf(stderr," protein_2:             The 2nd input protein file in PDB format. \n\n");
	fprintf(stderr,"Options: \n\n");
	fprintf(stderr,"-x range_1:             The residue range for the 1st input protein, (e.g., A:1-150) \n");
	fprintf(stderr,"-y range_2:             The residue range for the 2nd input protein, \n");
	fprintf(stderr,"                        If not specified, all residues in the first chain will be used. See README for more details. \n\n");
	fprintf(stderr,"-a alignment:           Specify an alignment in FASTA format from which to evaluate it. \n");
	fprintf(stderr,"                        If not specified, the sequence alignment between the two input proteins will be used. \n");
	fprintf(stderr,"                        If set to -1, then will align the two proteins according to the residue number. \n\n");
	fprintf(stderr,"-o out_name:            Specify an output file name for the evaluation. If not specified, \n");
	fprintf(stderr,"                        screen output the evaluation scores. \n\n");
	fprintf(stderr,"-d/D detail_file:       Specify an output file to show the detailed evaluation information, \n");
	fprintf(stderr,"                        for each aligned position with SeqID, BLOSUM, CLESUM, DeepScore, RMSD, and TMscore.\n");
	fprintf(stderr,"                        If use option -D, then will output the assessment of all heavy atoms. \n\n");
	fprintf(stderr,"-s script_option:      [0], do not output script files. (Set as default)\n");
	fprintf(stderr,"                        1,  output script for JMol, if -o is specified.\n");
	fprintf(stderr,"                        2,  output script for RasMol, if -o is specified.\n");
	fprintf(stderr,"                        3,  output script for PyMol, if -o is specified.\n\n");
	fprintf(stderr,"-P screenout:           0,  simple screenout evaluation scores. \n");
	fprintf(stderr,"                       [1], detailed screenout evaluation scores. (Set as default) \n\n");
	fprintf(stderr,"-n normalize_len:       Specify a normalization length for the calculation of TMscore,MAXSUB,GDT_TS/HA. In particular,\n");
	fprintf(stderr,"                       [0], the minimal length of the 1st and 2nd input protein. (Set as default)\n");
	fprintf(stderr,"                       -1,  the length of the 1st input protein; -3, (max_resnum-min_resnum+1) of the 1st model. \n");
	fprintf(stderr,"                       -2,  the length of the 2nd input protein; -4, (max_resnum-min_resnum+1) of the 2nd model. \n\n");
	fprintf(stderr,"-C distance_cut:        Specify a distance cutoff to remove residue pairs whose distance exceeds the threshold. \n");
	fprintf(stderr,"                       [0], keep all residue pairs. (Set as default) \n");
	fprintf(stderr,"                       -1,  automatically assign a distance cutoff value according to d0 in TMscore \n\n");
	fprintf(stderr,"-T tm_type:             Specify the type of TMscore calculation. \n");
	fprintf(stderr,"                        0,  use the alignment directly to superimpose the input structures. \n");
	fprintf(stderr,"                       [1], use typical approach to calculate the TMscore. (Set as default) \n");
	fprintf(stderr,"                        2,  use simplifed approach to calculate the TMscore. \n\n");
	fprintf(stderr,"Simple screenout description (please refer to README file for more details):\n");
	fprintf(stderr,"   name1 name2 len1 len2 -> BLOSUM CLESUM DeepScore -> LALI RMSDval TMscore -> MAXSUB GDT_TS GDT_HA -> SeqID nLen dCut \n");

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
		if(retv==0) //-> new style: i.e., ./DeepScore pdb1 pdb2
		{
			NEWorOLD=1;
		}
		else        //-> old style: i.e., ./DeepScore -t pdb1 -q pdb2
		{
			NEWorOLD=0;
		}
	}

	//---- DeepAlign_Score -------//
	{
		if(argc<3)
		{
			Usage();
			exit(-1);
		}
		//basic input
		string name1="";
		string name2="";
		string ali_file="";
		int ali_res=0;
		string range1="_";
		string range2="_";
		string out_file="";
		string detail_file="";
		int detail_full=0;
		int Out_Script=0; // no script out
		int Out_Screen=1; // detailed screen out
		int Normalize=0;
		double Distance_Cutoff=0;
		int TM_type=1;

		//---- process argument ----//
		extern char* optarg;
		int c=0;
		if(NEWorOLD==0) //old style
		{
			while((c=getopt(argc,argv,"t:q:o:O:d:D:s:a:x:y:P:n:C:T:"))!=EOF)
			{
				switch(c) 
				{
					//---- pairsie job ---//
					case 't':
						name1 = optarg;
						break;
					case 'q':
						name2 = optarg;
						break;
					case 'o':
						out_file = optarg;
						break;
					case 'O':
						Out_PDB = atoi(optarg);
						break;
					case 'd':
						detail_file = optarg;
						detail_full = 0;
						break;
					case 'D':
						detail_file = optarg;
						detail_full = 1;
						break;
					case 's':
						Out_Script = atoi(optarg);;
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
					case 'P':
						Out_Screen = atoi(optarg);
						break;
					case 'n':
						Normalize = atoi(optarg);
						break;
					case 'C':
						Distance_Cutoff = atof(optarg);
						break;
					case 'T':
						TM_type = atoi(optarg);
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
			name1 = argv[1];
			name2 = argv[2];
			while((c=getopt(argc,argv,"o:d:D:s:a:x:y:P:n:C:T:"))!=EOF)
			{
				switch(c) 
				{
					//---- pairsie job ---//
					case 'o':
						out_file = optarg;
						break;
					case 'd':
						detail_file = optarg;
						detail_full = 0;
						break;
					case 'D':
						detail_file = optarg;
						detail_full = 1;
						break;
					case 's':
						Out_Script = atoi(optarg);;
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
					case 'P':
						Out_Screen = atoi(optarg);
						break;
					case 'n':
						Normalize = atoi(optarg);
						break;
					case 'C':
						Distance_Cutoff = atof(optarg);
						break;
					case 'T':
						TM_type = atoi(optarg);
						break;

					//----- default ----//
					default:
						Usage();
						exit(-1);
				}
			}
		}

		//----- check argument ----//
		if(name1=="")
		{
			fprintf(stderr,"ERROR: FILE1 NAME NULL \n");
			exit(-1);
		}
		if(name2=="")
		{
			fprintf(stderr,"ERROR: FILE2 NAME NULL \n");
			exit(-1);
		}
		if(Normalize<-4)
		{
			fprintf(stderr,"ERROR: normal_len must be interger >= -4 \n");
			exit(-1);
		}
		if(Out_Script<0 || Out_Script>3)
		{
			fprintf(stderr,"ERROR: script_option must be 0,1,2 or 3 \n");
			exit(-1);
		}
		if(ali_file=="-1")
		{
			ali_file="";
			ali_res=1;
		}

		//---- main process ----//
		int SIMPLY_LOAD=0;
		Main_Process(name1,range1,name2,range2,ali_file,ali_res,out_file,detail_file,detail_full,
			Out_Script,Normalize,Distance_Cutoff,Out_Screen,SIMPLY_LOAD,TM_type);
		exit(0);
	}
}
