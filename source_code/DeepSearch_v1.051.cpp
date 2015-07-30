#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <cstring>
#include <algorithm>
#include <time.h>
#include <math.h>
#include <iomanip>
#include "getopt.h"

//--- using multiple thread ----//__130630__//
#include <pthread.h>
//--- using multiple thread ----//__130630__//over

using namespace std;


//------- usage -----//
void Usage()
{
	//-------- header part -------//
	cerr << "DeepSearch v1.135 (pThread version) [Feb-02-2015] \n\n";
	cerr << "Search a protein database for a query protein structure. \n";
	cerr << "-------------------------------------------------------- \n";
	cerr << "Sheng Wang, Jianzhu Ma, Jian Peng and Jinbo Xu.\n";
	cerr << "    PROTEIN STRUCTURE ALIGNMENT BEYOND SPATIAL PROXIMITY \n";
	cerr << "                     Scientific Reports, 3, 1448, (2013) \n";
	cerr << "-------------------------------------------------------- \n\n";

	//-------- option part ------//
	cerr << "Usage: \n";
	cerr << "./DeepSearch -a NP -q query_file [-l data_list] [-d data_root] [-o out_file] [-t tmsco] [-p pval] [-n topN]\n\n";
//	cerr << "                    [-t tmsco] [-p pval] [-n topN] [-s score_func] [-c]\n\n";
	cerr << "Options: \n\n";
	cerr << "-a NP:                 Number of processors.\n\n";
	cerr << "-q query_file:         Query protein file.  \n\n";
	cerr << "-l data_list:          The list of protein database [default = databases/bc40_list].\n\n";
	cerr << "-d data_root:          The folder containing the database files (i.e., .pdb files) [default = databases/pdb_BC100/].\n\n";                
	cerr << "-o out_file:           The file containing a brief summary of the searching results [default = query_name.rank ].\n\n";
	cerr << "-t tmsco:              Apply TMscore cutoff during searching process [default = 0.35]. \n\n";
	cerr << "-p pval:               Keep the results for top proteins according to P-value cutoff [default = 0.001]. \n\n";
	cerr << "-n topN:               Keep the results for top topN proteins [default = 100].\n\n";
//	cerr << "-s score_func:         1:dist-score, 2:vect-score, 4:local-score. [default score_func is 7, i.e., using all]. \n\n";
//	cerr << "-c :                   If specified, then the final template structure will be cut according to the alignment. \n\n";
	cerr << "        Example: \n";
	cerr << "            To search the protein database in which any two proteins share <=70% seq id, please run the following command:\n";
	cerr << "                ./DeepSearch -a 12 -q 1pazA.pdb -l databases/bc70_list -d databases/pdb_BC100/ -t 0.35 -n 200 \n";
	cerr << "            The topN proteins or those proteins with p-value smaller than the given cutoff will be saved to 1pazA.rank,\n";
	cerr << "                which also contains detailed alignments of the target to these top proteins.\n";
	cerr << "            The alignments to individual templates are saved in tmp/1pazA/RESULT/.\n\n";
}


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

//----------upper_case-----------//
void toUpperCase(char *buffer)
{
	for(int i=0;i<(int)strlen(buffer);i++)
	if(buffer[i]>=97 && buffer[i]<=122) buffer[i]-=32;
}
void toUpperCase(string &buffer)
{
	for(int i=0;i<(int)buffer.length();i++)
	if(buffer[i]>=97 && buffer[i]<=122) buffer[i]-=32;
}
//----------lower_case-----------//
void toLowerCase(char *buffer)
{
	for(int i=0;i<(int)strlen(buffer);i++)
	if(buffer[i]>=65 && buffer[i]<=90) buffer[i]+=32;
}
void toLowerCase(string &buffer)
{
	for(int i=0;i<(int)buffer.length();i++)
	if(buffer[i]>=65 && buffer[i]<=90) buffer[i]+=32;
}
//------- int2string -----//
string convertInt(int number)
{
	stringstream ss;  //create a stringstream
	ss << number;     //add number to the stream
	return ss.str();  //return a string with the contents of the stream
}
//------- double2string ----//
string convertDouble(double number,int precision)
{
	stringstream ss;  //create a stringstream
	ss << setprecision(precision) << number;     //add number to the stream
	return ss.str();  //return a string with the contents of the stream
}


//------- shuffle and cut -------//
int Shuffle_And_Cut_List(string &list,string &outnam,int num)
{
	ifstream fin;
	string buf,temp;
	fin.open(list.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"list file not found [%s] !!!\n",list.c_str());
		return -1;
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
	sprintf(command,"%s.%d",outnam.c_str(),iter);
	name=command;
	fp=fopen(name.c_str(),"wb");
	for(i=0;i<totnum;i++)
	{
		if(count==bignum && iter<num-1)
		{
			fclose(fp);
			iter++;
			sprintf(command,"%s.%d",outnam.c_str(),iter);
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
//		fprintf(stderr,"bad here !! totnum=[%d], iter[%d]!=num[%d] \n",totnum,iter,num);
//		return -1;
		for(i=iter+1;i<=num;i++)
		{
			sprintf(command,"%s.%d",outnam.c_str(),i);
			name=command;
			fp=fopen(name.c_str(),"wb");
			fclose(fp);
		}
	}
	return 1;
}



//-------------------- fasta related -----------------//
int ReadToFile_FASTA(string &fn,vector<pair<int, int> > &alignment,
	string &nam1_content,string &nam2_content,string &nam1_full,string &nam2_full,
	string &nam1,string &nam2,int &idrate)
{
	int i;
	int cur1=0;
	int cur2=0;
	int len;
	int len1,len2;
	int lali,match;
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
	idrate=0;
	lali=0;
	match=0;
	for(i=0;i<len;i++)
	{
		if(tmp[i]!='-' && seq[i]!='-') //match
		{
			nam1_content.push_back(tmp[i]);
			nam2_content.push_back(seq[i]);
			cur1++;
			cur2++;
			alignment.push_back(pair<int,int>(cur1,cur2));
			//calculate sequence idenrity
			if(tmp[i]==seq[i])match++;
			lali++;			
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
	idrate=match;
	return 1; //success

badend:
	fprintf(stderr,"alignment file format bad [%s] !!!\n",fn.c_str());
	return -1;
}

string detect_alignment_range(string path,pair<int,int> & template_range,pair<int,int> & query_range,vector<int> & aaa,int &idrate)
{
	vector<pair<int,int> > alignment;string nam1_content,nam2_content,nam1_full,nam2_full,nam1,nam2;
	ReadToFile_FASTA(path,alignment,nam1_content,nam2_content,nam1_full,nam2_full,nam1,nam2,idrate);
	//---- init check length ----//
	if(nam1_full.length()!=nam2_full.length() || nam1_full.length()!= alignment.size())
	{
		cerr<<" alignment wrong!"<<endl;
		exit(-1);
	}	
	
	//-------- get start and end pos -------//
	int align_start, align_end;
	int temp_start,temp_end,query_start,query_end;
	align_start=0;
	align_end=0;
	temp_start=0;
	temp_end=0;
	query_start=0;
	query_end=0;
	//---- get start pos ----//
	for(int i=0;i<alignment.size();i++)
	{
		if(alignment[i].first>0 && alignment[i].second>0)
		{
			align_start = i;
			temp_start = alignment[i].first;
			query_start = alignment[i].second;
			break;
		}
	}
	//---- get end pos ----//
	for(int i=alignment.size()-1;i>=0;i--)
	{
		if(alignment[i].first>0 && alignment[i].second>0)
		{
			align_end = i;
			temp_end = alignment[i].first;
			query_end = alignment[i].second;
			break;
		}
	}
	
	//---- get gap number ---//
	int sum = 0;int align_sum=0;int seq_gap=0;int temp_gap=0;int mid_gap_flag=0;
	for(int i=0;i<nam1_full.length();i++)
	{
		if(nam1_full[i]=='-' && i>=align_start && i<= align_end)
			temp_gap++;
		if(nam2_full[i]=='-' && i>=align_start && i<= align_end)
			seq_gap++;
		if(nam1_full[i]!='-' && nam2_full[i]!='-')
			align_sum ++;
	}

	//---- assign result ----//
	template_range.first = temp_start;
	template_range.second = temp_end;
	query_range.first = query_start;
	query_range.second = query_end;
	aaa[0] = align_sum;
	aaa[1] = seq_gap;
	aaa[2] = temp_gap;
	aaa[4] = nam1_content.length();
	return nam2_content;
}

//------ get list -----//
int WS_Get_Name_List(string &file,vector <string> &nam_list)
{
	ifstream fin;
	string buf;
	fin.open(file.c_str(), ios::in);
	if(fin.fail()!=0)return 0;
	int count=0;
	nam_list.clear();
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		nam_list.push_back(buf);
		count++;
	}
	return count;
}

//-------------------- Zcore Evalue related -----------------//
void Read_Mean_Vari(string &infile,double &mean,double &vari)
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
	mean=atof(buf.c_str());
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
	vari=atof(buf.c_str());
}
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
int WS_Get_Score_List(string &file,vector <double> &score_list)
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
int WS_Get_Score_List(string &file,vector <int> &score_list)
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
		int score=atoi(buf.c_str());
		score_list.push_back(score);
		count++;
	}
	return count;
}
int WS_Get_Cutoff_Given_EVD(vector <double> &score, double miu, double beta, double thres)
{
	int i;
	int size=(int)score.size();
	for(i=0;i<size;i++)
	{
		double ws_pvalue_rank=Pvalue(score[i],1/beta,miu);
		if(ws_pvalue_rank>thres)return i;
	}
}
void WS_Calculate_Pvalue(vector <double> &score,double miu, double beta, vector <double> &output)
{
	int i;
	int size=(int)score.size();
	output.clear();
	for(i=0;i<size;i++)
	{
		double ws_pvalue_rank=Pvalue(score[i],1/beta,miu);
		output.push_back(ws_pvalue_rank);
	}
}

//----------------------- get file length related ----------------//
//-> list length
int Get_List_Len(string &file)
{
	ifstream fin;
	string buf,temp;
	fin.open(file.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"file not found [%s] !!!\n",file.c_str());
		return -1;
	}
	int count=0;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		count++;
	}
	return count;
}
//-> 1pazA len
int Get_Head_Nam(string &list,vector <string> & output)
{
	//read
	ifstream fin;
	string buf,temp;
	fin.open(list.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"list %s not found!!\n",list.c_str());
		return -1;
	}
	output.clear();
	int count=0;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		istringstream www(buf);
		www>>temp;
		output.push_back(temp);
		count++;
	
	}
	return count;
}
//-> 1pdbA len
int Get_Head_Len(string &list)
{
	//read
	ifstream fin;
	string buf;
	fin.open(list.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"list %s not found!!\n",list.c_str());
		return -1;
	}
	if(!(fin>>buf))return -1;
	if(!(fin>>buf))return -1;
	return atoi(buf.c_str());
}
//-> read PDB len
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
		return -1;
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
		temp=buf.substr(0,3);
		if(temp=="TER"||temp=="END")break;
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


//=========== pthread command ============//
//-> global variable
string SEQ_NAME;
string SEQ_FILE;
string CAL_ROOT;
string PDB_ROOT;
string OUT_ROOT;
double MEAN;
double VARI;
double TMSCO_CUTOFF;
int SCORE_FUNC;


//-> for zscore and evalue
void *DeepSearch_Pthread_Part_I(void *arg)
{
	int retv;
	int proc_id = (long) arg;
	char wscommand[30000];
	string deepalign="DeepAlign -P 0 -u 0 "; //-> simple screenout and fastest calculation
	if(proc_id==0)deepalign+=" -v ";   //head node, then display
	//get maxsize
	int maxsize=3000;
	//for zscore
	sprintf(wscommand,"./%s -f tmp/%s_deepalign_refer_list.%d -r %s -q %s -m 3 -g 0 -h 1 -c 0 -p 0 -b tmp/%s_deepalign_zscore_rank.%d -e %d -w tmp/%s.tmpout",
		deepalign.c_str(),SEQ_NAME.c_str(),proc_id,CAL_ROOT.c_str(),SEQ_FILE.c_str(),SEQ_NAME.c_str(),proc_id,maxsize,SEQ_NAME.c_str());
	retv=system(wscommand);
	//for evalue
	sprintf(wscommand,"./%s -f tmp/%s_deepalign_refer_list.%d -r %s -q %s -j 0 -m 0 -p 0 -w tmp/%s_deepalign_evalue_rank.%d -e %d -s %d ",
		deepalign.c_str(),SEQ_NAME.c_str(),proc_id,CAL_ROOT.c_str(),SEQ_FILE.c_str(),SEQ_NAME.c_str(),proc_id,maxsize,SCORE_FUNC);
	retv=system(wscommand);
}

//-> for main search
void *DeepSearch_Pthread_Part_II(void *arg)
{
	int retv;
	int proc_id = (long) arg;
	char wscommand[30000];
	string deepalign="DeepAlign -P 0 -u 0 "; //-> simple screenout and fastest calculation
	if(proc_id==0)deepalign+=" -v ";   //head node, then display
	//get maxsize
	int maxsize=3000;
	//process
	sprintf(wscommand,"./%s -f tmp/%s_deepalign_pdb_list.%d -r %s -q %s -m 2 -g %lf -h %lf -c %lf -p 0 -w tmp/%s_deepalign_pdb_rank.%d -e %d -s %d ",
		deepalign.c_str(),SEQ_NAME.c_str(),proc_id,PDB_ROOT.c_str(),SEQ_FILE.c_str(),MEAN,VARI,TMSCO_CUTOFF,SEQ_NAME.c_str(),proc_id,maxsize,SCORE_FUNC);
	retv=system(wscommand);
}

//-> for final search
void *DeepSearch_Pthread_Part_III(void *arg)
{
	int retv;
	int proc_id = (long) arg;
	char wscommand[30000];
	string deepalign="DeepAlign -P 0 -u 1 ";    //-> simple screenout and normal calculation
	if(proc_id==0)deepalign+=" -v ";   //head node, then display
	string output_dir = "tmp/" + SEQ_NAME + "/" + OUT_ROOT + "/";
	//get maxsize
	int maxsize=3000;
	//process
	sprintf(wscommand,"./%s -f tmp/%s_deepalign_pdb_list.%d -r %s -q %s -m 0 -p 1 -d %s -w tmp/%s_deepalign_pdb_rank.%d -e %d -s %d ",
		deepalign.c_str(),SEQ_NAME.c_str(),proc_id,PDB_ROOT.c_str(),SEQ_FILE.c_str(),output_dir.c_str(),SEQ_NAME.c_str(),proc_id,maxsize,SCORE_FUNC);
	retv=system(wscommand);
}

//-> cut pdb
void *DeepSearch_Pthread_Cut_PDB(void *arg)
{
	int retv;
	int proc_id = (long) arg;
	char wscommand[30000];
	sprintf(wscommand,"tmp/%s_deepalign_pdb_list.%d",SEQ_NAME.c_str(),proc_id);
	string infile=wscommand;
	vector <string> proc_list;
	int size=Get_Head_Nam(infile,proc_list);
	for(int i=0;i<size;i++)
	{
		sprintf(wscommand,"util/Domain_Proc tmp/%s/%s/%s-%s.fasta tmp/%s-%s.fasta_cut 0 0 0 1",
			SEQ_NAME.c_str(),OUT_ROOT.c_str(),proc_list[i].c_str(),SEQ_NAME.c_str(),proc_list[i].c_str(),SEQ_NAME.c_str());
		retv=system(wscommand);
		sprintf(wscommand,"util/PDB_File_Cut %s/%s.pdb tmp/%s-%s.fasta_cut_seq1 tmp/%s_tmp_pdb/%s.pdb 0",
			PDB_ROOT.c_str(),proc_list[i].c_str(),proc_list[i].c_str(),SEQ_NAME.c_str(),SEQ_NAME.c_str(),proc_list[i].c_str());
		retv=system(wscommand);
	}
}

//-> re-align process
void *DeepSearch_Pthread_Re_Align(void *arg)
{
	int retv;
	int proc_id = (long) arg;
	char wscommand[30000];
	string deepalign="DeepAlign -P 0 -u 1 "; //-> simple screenout and normal calculation
	if(proc_id==0)deepalign+=" -v ";   //head node, then display
	string output_dir = "tmp/" + SEQ_NAME + "/" + OUT_ROOT + "/";
	//get maxsize
	int maxsize=3000;
	//process
	sprintf(wscommand,"./%s -f tmp/%s_deepalign_pdb_list__.%d -r tmp/%s_tmp_pdb/ -q %s -m 0 -p 1 -i 0 -d %s -w tmp/%s_deepalign_pdb_rank.%d -e %d -s %d ",
		deepalign.c_str(),SEQ_NAME.c_str(),proc_id,SEQ_NAME.c_str(),SEQ_FILE.c_str(),output_dir.c_str(),SEQ_NAME.c_str(),proc_id,maxsize,SCORE_FUNC);
	retv=system(wscommand);
}

//-> clear temp
void *DeepSearch_Pthread_Clean(void *arg)
{
	int retv;
	int proc_id = (long) arg;
	char wscommand[30000];
	sprintf(wscommand,"tmp/%s_deepalign_pdb_list.%d",SEQ_NAME.c_str(),proc_id);
	string infile=wscommand;
	vector <string> proc_list;
	int size=Get_Head_Nam(infile,proc_list);
	for(int i=0;i<size;i++)
	{
		sprintf(wscommand,"rm -f tmp/%s-%s.fasta_cut*",proc_list[i].c_str(),SEQ_NAME.c_str());
		retv=system(wscommand);
		sprintf(wscommand,"rm -f tmp/%s_tmp_pdb/%s.pdb",SEQ_NAME.c_str(),proc_list[i].c_str());
		retv=system(wscommand);
	}
}


//========================= main ==================//
int main(int argc, char ** argv) 
{
int num_procs=1;

	//-- help --//
	if (argc < 2) 
	{
		Usage();
		return 0;
	}

	//-> file related
	string seq_file;
	string pdb_list = "databases/bc40_list";
	string pdb_root = "databases/pdb_BC100/";
	string out_file;
	//-> parameter related
	double tmsco_cutoff=0.35;   //-> 0.35
	double pvalue_cutoff=0.001; //-> 0.001
	int N=100;                  //-> minimal output number
	int M=1000;                 //-> maximal output number
	//-> others
	int score_func=7;           //-> default: 7 (using all for score_func)
	int cal_flag = 1;           //-> default: 1 (calculate p-value)
	double ws_alpha_value=1.0;  //-> default: 1.0 (mpirun related list divide)
	int MAXIMAL_PDB_SIZE=3000;  //-> default: 3000
	int CUT_TEMPLATE=0;         //-> default: we don't cut template

	//-- parse arguments --//
	extern char* optarg;
	extern int optind;
	int c=0;
	while((c=getopt(argc,argv,"a:q:l:d:o:t:p:n:s:z:x:c"))!=EOF)
	{
		switch(c) 
		{
		case 'a':
			num_procs=atoi(optarg);
			break;
		//-> file related
		case 'q':
			seq_file=optarg;
			break;
		case 'l':
			pdb_list=optarg;
			break;
		case 'd':
			pdb_root=optarg;
			break;
		case 'o':
			out_file=optarg;
			break;
		//-> parameter related
		case 't':
			tmsco_cutoff = atof(optarg);
			break;
		case 'p':
			pvalue_cutoff = atof(optarg);
			break;
		case 'n':
			N=atoi(optarg);
			break;
		case 's':
			score_func=atoi(optarg);
			break;
		//-> maximal size 
		case 'z':
			M=atoi(optarg);   //-> maximal output size
			break;
		case 'x':
			MAXIMAL_PDB_SIZE=atoi(optarg); //-> maximal PDB size
			break;
		//-> cut template
		case 'c':
			CUT_TEMPLATE=1;
			break;
		//-> default
		default:
			Usage();
			return 0;
		}
	}

	//-- input check --//
	if(seq_file.length()==0){
		cerr<<"Error: Query protein file missing!"<<endl;
		return 0;
	}
	if(pdb_list.length()==0){
		cerr<<"Warning: No protein datebase specified, will use default datebase(BC40)"<<endl;
		pdb_list = "databases/bc40_list";
		pdb_root = "databases/pdb_BC100/";
	}
	if(pdb_root.length()==0){
		cerr<<"Warning: No protein files root specified, will use defaut databases/pdb_BC100/"<<endl;;
		pdb_root = "databases/pdb_BC100/";
	}
	//-- list check --//
	ifstream temp_in(pdb_list.c_str());
	if(!temp_in.is_open()){
		cerr << "Error: Cannot open pdb_list: " << pdb_list << endl;
		return 0;
	}
	//-- out check --//
	string seq_name;
	getBaseName(seq_file,seq_name,'/','.');
	if(out_file.length()==0){
		out_file = seq_name + ".rank";
	}
	int tgt_file_len=Get_PDB_File_Len(seq_file);
	if(tgt_file_len==-1)
	{
		cerr << "Error: Query pdb file not found or format bad ! " << seq_file << endl;
		return 0;
	}
	string OUTROOT="RESULT";


	//------------ pThread allocate ------------//__130630__//
	pthread_t *threads;
	threads = new pthread_t[num_procs];
	SEQ_NAME=seq_name;
	SEQ_FILE=seq_file;
	PDB_ROOT=pdb_root;
	OUT_ROOT=OUTROOT;
	TMSCO_CUTOFF=tmsco_cutoff;
	SCORE_FUNC=score_func;
	//------------ pThread allocate ------------//__130630__//over
	

	//======================= PART I: REFERENCE_SCAN =====================//
	int retv;
	int TOTN=Get_List_Len(pdb_list);
	char wscommand[30000];
	//[1] make temporar directory
//	if(proc_id==0)
	{
		sprintf(wscommand,"mkdir -p tmp/%s_tmp_pdb",seq_name.c_str());
		retv=system(wscommand);
		sprintf(wscommand,"mkdir -p tmp/%s/%s",seq_name.c_str(),OUTROOT.c_str());
		retv=system(wscommand);
	}

	//[2] calculate reference_list
	string infile="";
	int maxsize=3000;
	double mean=0;
	double vari=1;
	double miu=0;
	double beta=1;
	string cal_list="databases/reference_pdb_list";
	string cal_root="databases/CAL_PDB";
	CAL_ROOT=cal_root;
	//[2-1] process reference_list with length
//	if(proc_id==0)
	{
		sprintf(wscommand,"rm -f tmp/%s_deepalign_refer_list*",seq_name.c_str());
		retv=system(wscommand);
		sprintf(wscommand,"rm -f tmp/%s_deepalign_zscore_rank*",seq_name.c_str());
		retv=system(wscommand);
		sprintf(wscommand,"rm -f tmp/%s_deepalign_evalue_rank*",seq_name.c_str());
		retv=system(wscommand);
		string out_nam="tmp/"+seq_name+"_deepalign_refer_list";
		retv=Shuffle_And_Cut_List(cal_list,out_nam,num_procs);
	}
	
	//[2-2] mpirun reference_list
	{

		//------------ pThread create ------------//
		for(int i=0; i < num_procs; i++ )
		{
			int rc = pthread_create(&threads[i], NULL, &DeepSearch_Pthread_Part_I, (void *)i );
			if (rc)
			{
				fprintf(stderr,"Error:unable to create thread at Part_I, with rc %d \n",rc);
				exit(-1);
			}
		}
		for(int i=0; i < num_procs; i++ )
		{
			int rc = pthread_join(threads[i], NULL);
			if (rc)
			{
				fprintf(stderr,"Error:unable to join thread at Part_I, with rc %d \n",rc);
				exit(-1);
			}
		}
		//------------ pThread create ------------//over
	}

	//[2-3] calculate zscore and evalue for reference_list
//	if(proc_id==0)
	{
		//-> clear "wsout"
		sprintf(wscommand,"rm -f tmp/%s.tmpout",seq_name.c_str());
		retv=system(wscommand);

		//-------- process temporary refer_list --------//
		//-> delete temporary refer_list
		sprintf(wscommand,"rm -f tmp/%s_deepalign_refer_list*",seq_name.c_str());
		retv=system(wscommand);
		//-> process zscore_rank
		sprintf(wscommand,"cat tmp/%s_deepalign_zscore_rank.* > tmp/%s_deepalign_zscore_rank",seq_name.c_str(),seq_name.c_str());
		retv=system(wscommand);
		sprintf(wscommand,"rm -f tmp/%s_deepalign_zscore_rank.*",seq_name.c_str());
		retv=system(wscommand);
		sprintf(wscommand,"sort -n -r -k6 tmp/%s_deepalign_zscore_rank > tmp/%s/%s.rank_zscore",seq_name.c_str(),seq_name.c_str(),seq_name.c_str());
		retv=system(wscommand);
		sprintf(wscommand,"rm -f tmp/%s_deepalign_zscore_rank",seq_name.c_str());
		retv=system(wscommand);
		//-> process evalue_rank
		sprintf(wscommand,"cat tmp/%s_deepalign_evalue_rank.* > tmp/%s_deepalign_evalue_rank",seq_name.c_str(),seq_name.c_str());
		retv=system(wscommand);
		sprintf(wscommand,"rm -f tmp/%s_deepalign_evalue_rank.*",seq_name.c_str());
		retv=system(wscommand);
		sprintf(wscommand,"sort -n -r -k8 tmp/%s_deepalign_evalue_rank > tmp/%s/%s.rank_evalue",seq_name.c_str(),seq_name.c_str(),seq_name.c_str());
		retv=system(wscommand);
		sprintf(wscommand,"rm -f tmp/%s_deepalign_evalue_rank",seq_name.c_str());
		retv=system(wscommand);
		
		//-------- get mean/vari and miu/beta ---------//
		//-> get mean/vari
		sprintf(wscommand,"awk '{print $6}' tmp/%s/%s.rank_zscore | tail -n+4 > tmp/%s.rank_zscore_val",
			seq_name.c_str(),seq_name.c_str(),seq_name.c_str());
		retv=system(wscommand);
		sprintf(wscommand,"util/Stat_List tmp/%s.rank_zscore_val > tmp/%s.rank_zscore_reso",
			seq_name.c_str(),seq_name.c_str());
		retv=system(wscommand);
		infile="tmp/"+seq_name+".rank_zscore_reso";
		Read_Mean_Vari(infile,mean,vari);
		sprintf(wscommand,"rm -f tmp/%s.rank_zscore_*",seq_name.c_str());
		retv=system(wscommand);
		//-> check mean/vari
		if(fabs(mean)<0.01 && fabs(vari-1)<0.01)
		{
			mean=0;
			vari=0.5; //-> this is for DeepAlign_Search
		}
		MEAN=mean;
		VARI=vari;
		//-> get miu/beta
		sprintf(wscommand,"awk '{print $8}' tmp/%s/%s.rank_evalue | tail -n+4 > tmp/%s.rank_evalue_val",
			seq_name.c_str(),seq_name.c_str(),seq_name.c_str());
		retv=system(wscommand);
		sprintf(wscommand,"util/Fitting_EVD tmp/%s.rank_evalue_val > tmp/%s.rank_evalue_reso",
			seq_name.c_str(),seq_name.c_str());
		retv=system(wscommand);
		infile="tmp/"+seq_name+".rank_evalue_reso";
		Read_Miu_Beta(infile,miu,beta);
		sprintf(wscommand,"rm -f tmp/%s.rank_evalue_*",seq_name.c_str());
		retv=system(wscommand);
	}

/*
#ifdef _MPI
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&mean, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&vari, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&miu, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&beta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
#endif
*/

	//======================= PART II: PDB_LIST_SCAN  =====================//
	//[1] process pdb_list with length
//	if(proc_id==0)
	{
		sprintf(wscommand,"rm -f tmp/%s_deepalign_pdb_list*",seq_name.c_str());
		retv=system(wscommand);
		sprintf(wscommand,"rm -f tmp/%s_deepalign_pdb_rank*",seq_name.c_str());
		retv=system(wscommand);
		string out_nam="tmp/"+seq_name+"_deepalign_pdb_list";
		retv=Shuffle_And_Cut_List(pdb_list,out_nam,num_procs);
	}

	//[2] mpirun pdb_list
	{
		//-> mainprocess
	
		//------------ pThread create ------------//_
		for(int i=0; i < num_procs; i++ )
		{
			int rc = pthread_create(&threads[i], NULL, &DeepSearch_Pthread_Part_II, (void *)i );
			if (rc)
			{
				fprintf(stderr,"Error:unable to create thread at Part_II, with rc %d \n",rc);
				exit(-1);
			}
		}
		for(int i=0; i < num_procs; i++ )
		{
			int rc = pthread_join(threads[i], NULL);
			if (rc)
			{
				fprintf(stderr,"Error:unable to join thread at Part_II, with rc %d \n",rc);
				exit(-1);
			}
		}
		//------------ pThread create ------------//over

	}

	//[3] result collect by Evalue
//	if(proc_id==0)
	{
		sprintf(wscommand,"rm -f tmp/%s_deepalign_pdb_list*",seq_name.c_str());
		retv=system(wscommand);
		sprintf(wscommand,"cat tmp/%s_deepalign_pdb_rank.* > tmp/%s_deepalign_pdb_rank",seq_name.c_str(),seq_name.c_str());
		retv=system(wscommand);
		sprintf(wscommand,"rm -f tmp/%s_deepalign_pdb_rank.*",seq_name.c_str());
		retv=system(wscommand);
		sprintf(wscommand,"sort -n -r -k8 tmp/%s_deepalign_pdb_rank > tmp/%s/%s.rank__",seq_name.c_str(),seq_name.c_str(),seq_name.c_str());
		retv=system(wscommand);
		sprintf(wscommand,"rm -f tmp/%s_deepalign_pdb_rank",seq_name.c_str());
		retv=system(wscommand);
	}

	//========= determine topN by p-value =========//__130401__//
	if( cal_flag==1 )
	{
//		if(proc_id==0)
		{
			//calculate pvalue-based cutoff
			sprintf(wscommand,"awk '{print $8}' tmp/%s/%s.rank__ > tmp/%s.rank__val",
				seq_name.c_str(),seq_name.c_str(),seq_name.c_str());
			retv=system(wscommand);
			string score_file="tmp/"+seq_name+".rank__val";
			vector <double> pvalue_score;
			WS_Get_Score_List(score_file,pvalue_score);
			//determine N
			int retN=WS_Get_Cutoff_Given_EVD(pvalue_score,miu,beta,pvalue_cutoff);
			int newN;
			int oriN=N;             //default 100
			int maxN=M<TOTN?M:TOTN; //default 1000
			if(retN<oriN)newN=oriN;
			else if(retN>maxN)newN=maxN;
			else newN=retN;
			N=newN;
			//delete temporary files
			sprintf(wscommand,"rm -f tmp/%s.rank__val",seq_name.c_str());
			retv=system(wscommand);
		}
	}


	//======================= PART III: PDB_LIST realign  =====================//
	//[1] sort list
//	if(proc_id==0)
	{
		sprintf(wscommand,"head -n %d tmp/%s/%s.rank__ | awk '{print $1}' > tmp/%s_deepalign_pdb_list_",
			N,seq_name.c_str(),seq_name.c_str(),seq_name.c_str());
		retv=system(wscommand);
		string in_nam="tmp/"+seq_name+"_deepalign_pdb_list_";
		string out_nam="tmp/"+seq_name+"_deepalign_pdb_list";
		retv=Shuffle_And_Cut_List(in_nam,out_nam,num_procs);
	}

	//[2] mpirun pdb_list
	{

		//------------ pThread create ------------//
		for(int i=0; i < num_procs; i++ )
		{
			int rc = pthread_create(&threads[i], NULL, &DeepSearch_Pthread_Part_III, (void *)i );
			if (rc)
			{
				fprintf(stderr,"Error:unable to create thread at Part_III, with rc %d \n",rc);
				exit(-1);
			}
		}
		for(int i=0; i < num_procs; i++ )
		{
			int rc = pthread_join(threads[i], NULL);
			if (rc)
			{
				fprintf(stderr,"Error:unable to join thread at Part_III, with rc %d \n",rc);
				exit(-1);
			}
		}
		//------------ pThread create ------------//over

	}

	//[3] cut template
	if(CUT_TEMPLATE==1)
	{
		//-> (1) cut pdb
		//------------ pThread create ------------//
		for(int i=0; i < num_procs; i++ )
		{
			int rc = pthread_create(&threads[i], NULL, &DeepSearch_Pthread_Cut_PDB, (void *)i );
			if (rc)
			{
				fprintf(stderr,"Error:unable to create thread at Cut_PDB, with rc %d \n",rc);
				exit(-1);
			}
		}
		for(int i=0; i < num_procs; i++ )
		{
			int rc = pthread_join(threads[i], NULL);
			if (rc)
			{
				fprintf(stderr,"Error:unable to join thread at Cut_PDB, with rc %d \n",rc);
				exit(-1);
			}
		}
		//------------ pThread create ------------//over
		
		//-> (2) get length 
//		if(proc_id==0)
		{
			sprintf(wscommand,"head -n %d tmp/%s/%s.rank__ | awk '{print $1}' > tmp/%s_deepalign_pdb_list_",
				N,seq_name.c_str(),seq_name.c_str(),seq_name.c_str());
			retv=system(wscommand);
			string in_nam="tmp/"+seq_name+"_deepalign_pdb_list_";
			string out_nam="tmp/"+seq_name+"_deepalign_pdb_list__";
			retv=Shuffle_And_Cut_List(in_nam,out_nam,num_procs);
		}

		//-> (3) re-align
		//------------ pThread create ------------//
		for(int i=0; i < num_procs; i++ )
		{
			int rc = pthread_create(&threads[i], NULL, &DeepSearch_Pthread_Re_Align, (void *)i );
			if (rc)
			{
				fprintf(stderr,"Error:unable to create thread at Re_Align, with rc %d \n",rc);
				exit(-1);
			}
		}
		for(int i=0; i < num_procs; i++ )
		{
			int rc = pthread_join(threads[i], NULL);
			if (rc)
			{
				fprintf(stderr,"Error:unable to join thread at Re_Align, with rc %d \n",rc);
				exit(-1);
			}
		}
		//------------ pThread create ------------//over
		
		//-> (4) clear temporary
		//------------ pThread create ------------//
		for(int i=0; i < num_procs; i++ )
		{
			int rc = pthread_create(&threads[i], NULL, &DeepSearch_Pthread_Clean, (void *)i );
			if (rc)
			{
				fprintf(stderr,"Error:unable to create thread at Clean, with rc %d \n",rc);
				exit(-1);
			}
		}
		for(int i=0; i < num_procs; i++ )
		{
			int rc = pthread_join(threads[i], NULL);
			if (rc)
			{
				fprintf(stderr,"Error:unable to join thread at Clean, with rc %d \n",rc);
				exit(-1);
			}
		}
		//------------ pThread create ------------//over

	}

	//[4] clear all and re-rank
//	if(proc_id==0)
	{
		sprintf(wscommand,"rm -f tmp/%s_deepalign_pdb_list*",seq_name.c_str());
		retv=system(wscommand);
		sprintf(wscommand,"cat tmp/%s_deepalign_pdb_rank.* > tmp/%s_deepalign_pdb_rank",seq_name.c_str(),seq_name.c_str());
		retv=system(wscommand);
		sprintf(wscommand,"rm -f tmp/%s_deepalign_pdb_rank.*",seq_name.c_str());
		retv=system(wscommand);
		sprintf(wscommand,"sort -n -r -k8 tmp/%s_deepalign_pdb_rank > tmp/%s/%s.rank_",seq_name.c_str(),seq_name.c_str(),seq_name.c_str());
		retv=system(wscommand);
		sprintf(wscommand,"rm -f tmp/%s_deepalign_pdb_rank",seq_name.c_str());
		retv=system(wscommand);
	}

	//[4] output p-value and other score
//	if(proc_id==0)
	{
		string score_file;
		vector <string> temp_list;
		vector <int> temp_length;
		vector <double> origin_score;
		vector <double> pvalue_score;
		vector <double> tmscore;
		vector <double> RMSD;
		vector <int> LALI;
		//get temp list
		sprintf(wscommand,"awk '{print $1}' tmp/%s/%s.rank_ > tmp/%s.rank_val",
			seq_name.c_str(),seq_name.c_str(),seq_name.c_str());
		retv=system(wscommand);
		score_file="tmp/"+seq_name+".rank_val";
		int retv_size=WS_Get_Name_List(score_file,temp_list);
		if(retv_size<N)N=retv_size;
		sprintf(wscommand,"rm -f tmp/%s.rank_val",seq_name.c_str());
		retv=system(wscommand);
		//get temp length
		sprintf(wscommand,"awk '{print $3}' tmp/%s/%s.rank_ > tmp/%s.rank_val",
			seq_name.c_str(),seq_name.c_str(),seq_name.c_str());
		retv=system(wscommand);
		score_file="tmp/"+seq_name+".rank_val";
		WS_Get_Score_List(score_file,temp_length);
		sprintf(wscommand,"rm -f tmp/%s.rank_val",seq_name.c_str());
		retv=system(wscommand);
		//get original score
		sprintf(wscommand,"awk '{print $8}' tmp/%s/%s.rank_ > tmp/%s.rank_val",
			seq_name.c_str(),seq_name.c_str(),seq_name.c_str());
		retv=system(wscommand);
		score_file="tmp/"+seq_name+".rank_val";
		WS_Get_Score_List(score_file,origin_score);
		WS_Calculate_Pvalue(origin_score,miu,beta,pvalue_score);
		sprintf(wscommand,"rm -f tmp/%s.rank_val",seq_name.c_str());
		retv=system(wscommand);
		//get tmscore
		sprintf(wscommand,"awk '{print $12}' tmp/%s/%s.rank_ > tmp/%s.rank_val",
			seq_name.c_str(),seq_name.c_str(),seq_name.c_str());
		retv=system(wscommand);
		score_file="tmp/"+seq_name+".rank_val";
		WS_Get_Score_List(score_file,tmscore);
		sprintf(wscommand,"rm -f tmp/%s.rank_val",seq_name.c_str());
		retv=system(wscommand);
		//get RMSD
		sprintf(wscommand,"awk '{print $11}' tmp/%s/%s.rank_ > tmp/%s.rank_val",
			seq_name.c_str(),seq_name.c_str(),seq_name.c_str());
		retv=system(wscommand);
		score_file="tmp/"+seq_name+".rank_val";
		WS_Get_Score_List(score_file,RMSD);
		sprintf(wscommand,"rm -f tmp/%s.rank_val",seq_name.c_str());
		retv=system(wscommand);
		//get LALI
		sprintf(wscommand,"awk '{print $10}' tmp/%s/%s.rank_ > tmp/%s.rank_val",
			seq_name.c_str(),seq_name.c_str(),seq_name.c_str());
		retv=system(wscommand);
		score_file="tmp/"+seq_name+".rank_val";
		WS_Get_Score_List(score_file,LALI);
		sprintf(wscommand,"rm -f tmp/%s.rank_val",seq_name.c_str());
		retv=system(wscommand);
	
		//[2-3] final output
		{
			string out_file_tmp=out_file+"_simp";
			FILE *fp=fopen(out_file_tmp.c_str(),"wb");
			ofstream output(out_file.c_str());
			//[2-3-1] output command
			output << "# Generated by Command: " << endl;
			output << "# mpirun -np "<<num_procs<<" ";
			for(int i=0;i<argc;i++)output<<argv[i]<<" ";
			output << endl;
			//[2-3-2] output reference
			output << "# Reference:" << endl;
			output << "# Sheng Wang, Jianzhu Ma, Jian Peng and Jinbo Xu." << endl;
			output << "#     PROTEIN STRUCTURE ALIGNMENT BEYOND SPATIAL PROXIMITY " << endl;
			output << "#                      Scientific Reports, 3, 1448, (2013) " << endl;
			output << "# ----------------------------------------------------------" <<endl;
			//[2-3-3] output header
			output << "# Query Name = " << seq_name << endl;
			string path_top = "tmp/"+ seq_name + "/" + OUTROOT + "/" + temp_list[0] + "-" + seq_name + ".fasta";
			pair<int,int> template_range2,query_range2;vector<int> bbb(4,0);
			int idrate;
			string query_protein_sequence = detect_alignment_range(path_top,template_range2,query_range2,bbb,idrate);
			output << "# Query Sequence = " << query_protein_sequence<<endl;
			output << "# Query Length = " << query_protein_sequence.length()<<endl;
			output << "# Searched Templates = "<< TOTN <<endl;
			time_t* tp=new(time_t);
			*tp=time(NULL);
			output << "# Date "<<ctime(tp);
			output << "# ---------------- Pvalue parameter: miu = "<<miu<<" beta = "<<beta<<" ----------------" <<endl<<endl;
			
			//[2-3-4] output alignment
			output <<setiosflags(ios::left)<<setw(5)<<"No"
						 <<setiosflags(ios::left)<<setw(8)<<"protein"
						 <<setiosflags(ios::left)<<setw(8)<<"tLength"
						 <<setiosflags(ios::left)<<setw(12)<<"Pvalue"
						 <<setiosflags(ios::left)<<setw(12)<<"Score"
						 <<setiosflags(ios::left)<<setw(12)<<"TMsco"
						 <<setiosflags(ios::left)<<setw(12)<<"RMSD"
						 <<setiosflags(ios::left)<<setw(12)<<"qRange"
						 <<setiosflags(ios::left)<<setw(12)<<"tRange"
						 <<setiosflags(ios::left)<<setw(8)<<"Lali"
						 <<setiosflags(ios::left)<<setw(8)<<"#tGaps"
						 <<setiosflags(ios::left)<<setw(8)<<"#qGaps"
						 <<setiosflags(ios::left)<<setw(8)<<"#seqID";
			output <<endl;
			for(int i=0;i<N;i++)
			{
				//init check
				if(pvalue_score[i]>pvalue_cutoff)continue;
				//printf
				pair<int,int> template_range,query_range;
				string path = "tmp/"+ seq_name +"/" + OUTROOT +"/" + temp_list[i] + "-" + seq_name + ".fasta";
				vector<int> aaa(4,0);
				detect_alignment_range(path,template_range,query_range,aaa,idrate);
				char no[256];
				sprintf(no,"%d",i+1);
				string query_scale = convertInt(query_range.first) + "-" + convertInt(query_range.second);
				string template_scale = convertInt(template_range.first) + "-" + convertInt(template_range.second); 
				output <<setiosflags(ios::left)<<setw(5)<<no
				       <<setiosflags(ios::left)<<setw(8)<<temp_list[i]
				       <<setiosflags(ios::left)<<setw(8)<<temp_length[i]
				       <<setiosflags(ios::left)<<setw(12)<<setprecision(4)<<pvalue_score[i]
				       <<setiosflags(ios::left)<<setw(12)<<setprecision(4)<<origin_score[i]
				       <<setiosflags(ios::left)<<setw(12)<<setprecision(4)<<tmscore[i]
				       <<setiosflags(ios::left)<<setw(12)<<setprecision(4)<<RMSD[i]
				       <<setiosflags(ios::left)<<setw(12)<<query_scale
				       <<setiosflags(ios::left)<<setw(12)<<template_scale
				       <<setiosflags(ios::left)<<setw(8)<<LALI[i]
				       <<setiosflags(ios::left)<<setw(8)<<aaa[2]
				       <<setiosflags(ios::left)<<setw(8)<<aaa[1]
				       <<setiosflags(ios::left)<<setw(8)<<idrate;
				output << endl;
				//fprintf out_tmp
				string pvalu_str=convertDouble(pvalue_score[i],4);
				string origi_str=convertDouble(origin_score[i],4);
				string tmsco_str=convertDouble(tmscore[i],4);
				string rmsd_str=convertDouble(RMSD[i],4);
				fprintf(fp,"%s %d %s %s %s %s %d %d\n",temp_list[i].c_str(),temp_length[i],
					pvalu_str.c_str(),origi_str.c_str(),tmsco_str.c_str(),rmsd_str.c_str(),LALI[i],idrate);
			}
			output << endl;
			fclose(fp);

			//--- final output alignment ----//
			for(int i=0;i<N;i++)
			{
				//init check
				if(pvalue_score[i]>pvalue_cutoff)continue;
				//printf
				string align_detail_file = "tmp/"+ seq_name +"/" + OUTROOT + "/" + temp_list[i] + "-" + seq_name + ".local";
				output<<"No "<<i+1<<endl<<">"<<temp_list[i]<<endl;
				ifstream fin;
				fin.open(align_detail_file.c_str(), ios::in);
				string buf;
				for(;;)
				{
					if(!getline(fin,buf,'\n'))break;
					output<<buf<<endl;
				}
				fin.close();
				fin.clear();
			}
			
			//[2-3-5] output done
			output<<"Done."<<endl;	
			output.close();
		}
	}

	//---------- delete temporary directory ----------//
//	if(proc_id==0)
	{
		sprintf(wscommand,"rmdir tmp/%s_tmp_pdb",seq_name.c_str());
		retv=system(wscommand);
	}

	//----- pthread terminate -----//
//	pthread_exit(NULL);
	delete [] threads;
	//----- pthread terminate -----//over


	return 0;
} // end of main function

