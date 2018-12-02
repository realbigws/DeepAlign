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


//---------- input a string, output a vector -----//
int String_To_Vector(string &input,vector <string> &output)
{
	istringstream www(input);
	output.clear();
	int count=0;
	string value;
	for(;;)
	{
		if(! (www>>value) )break;
		output.push_back(value);
		count++;
	}
	return count;
}
int String_To_Vector(string &input,vector <string> &output, char separator)
{
	istringstream www(input);
	output.clear();
	int count=0;
	string value;
	for(;;)
	{
		if(!getline(www,value,separator))break;
		output.push_back(value);
		count++;
	}
	return count;
}


//------- read EVD parameter: miu and beta ---------//
//-> example
/*
miu = %lf , beta = %lf
*/
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
		cerr<<path<<" alignment wrong!"<<endl;
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
	aaa[3] = nam1_content.length();
	return nam2_content;
}


//--------- load RankSimp file --------------//
//-> format
/*
#---- screen output format -----#
#   name1 name2 len1 len2 -> BLOSUM CLESUM DeepScore -> LALI RMSDval TMscore -> MAXSUB GDT_TS GDT_HA -> SeqID nLen dCut
#       1     2    3    4 (5)     6      7         8 (9)  10      11      12 (13)   14     15     16 (17)  18   19   20
*/
int Load_RankSimp_File(string &filename, vector <string> &temp_rank, vector <vector <double> > &temp_score)
{
	ifstream fin;
	string wbuf,temp;
	fin.open(filename.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"%s not found!\n",filename.c_str());
		exit(-1);
	}
	temp_rank.clear();
	temp_score.clear();
	int count=0;
	int col_size=21;
	for(;;)
	{
		if(!getline(fin,wbuf,'\n'))break;
		vector <string> str_out;
		int retv=String_To_Vector(wbuf,str_out);
		if(retv!=col_size)
		{
			fprintf(stderr,"col_size %d not equal to %d for %s \n",
				retv,col_size,wbuf.c_str());
			exit(-1);
		}
		temp_rank.push_back(str_out[0]);
		vector <double> score(4,0);
		score[0]=atof(str_out[7].c_str());       //-> origin_score[i]
		score[1]=atof(str_out[11].c_str());      //-> tmscore[i]
		score[2]=atof(str_out[10].c_str());      //-> RMSD[i]
		score[3]=atoi(str_out[9].c_str());       //-> LALI[i]
		temp_score.push_back(score);
		count++;
	}
	//return
	return count;
}


//----------- calculate Pvalue related -----------//
double Pvalue(double x, double lamda,double miu)
{
	double h = lamda * (x - miu);
	return (h > 10) ? exp(-h) : double(1.0) - exp(-exp(-h));
}

//--------------------- generate a detailed rank file ----------------//
void DeepSearch_RankFile_Generate(string &output_name, string &seq_name, string &align_root, 
	vector <string> &temp_rank, vector <vector <double> > &temp_score, double &pvalue_cutoff,
	string &command, double &miu, double &beta, string &TOTN, int TopN)
{
	int N=TopN;
	if(N<=0)return;

	//[2-3] final output
	{
		ofstream output(output_name.c_str());
		//[2-3-1] output command
		output << "# Generated by Command: ";
		output << "#   "<<command;    //-> need to provide the command //__180826__//
		output<<endl<<endl;

		//[2-3-2] output reference
		output << "# Reference:" << endl;
		output << "# Sheng Wang, Jianzhu Ma, Jian Peng and Jinbo Xu." << endl;
		output << "#     PROTEIN STRUCTURE ALIGNMENT BEYOND SPATIAL PROXIMITY " << endl;
		output << "#                      Scientific Reports, 3, 1448, (2013) " << endl;
		output << "# ----------------------------------------------------------" <<endl;

		//[2-3-3] output header
		output << "# Query Name = " << seq_name << endl;
		string path_top = align_root + "/" + temp_rank[0] + "-" + seq_name + ".fasta";
		pair<int,int> template_range2,query_range2;vector<int> aaa(4,0);int idrate;
		string query_protein_sequence = detect_alignment_range(path_top,template_range2,query_range2,aaa,idrate);
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
			//init pvalue check
			double pvalue=0;
			if(pvalue_cutoff>0)
			{
				pvalue=Pvalue(temp_score[i][0],1/beta,miu);
				if(pvalue>pvalue_cutoff)continue;
			}
			//output detailed rank
			{
				pair<int,int> template_range,query_range;
				string path = align_root +"/" + temp_rank[i] + "-" + seq_name + ".fasta";
				detect_alignment_range(path,template_range,query_range,aaa,idrate);
				char no[256];
				sprintf(no,"%d",i+1);
				string query_scale = convertInt(query_range.first) + "-" + convertInt(query_range.second);
				string template_scale = convertInt(template_range.first) + "-" + convertInt(template_range.second);
				output <<setiosflags(ios::left)<<setw(5)<<no
				       <<setiosflags(ios::left)<<setw(8)<<temp_rank[i]
				       <<setiosflags(ios::left)<<setw(8)<<aaa[3]                                //-> tLength
				       <<setiosflags(ios::left)<<setw(12)<<setprecision(4)<<pvalue              //-> pvalue_score[i]
				       <<setiosflags(ios::left)<<setw(12)<<setprecision(4)<<temp_score[i][0]    //-> origin_score[i]
				       <<setiosflags(ios::left)<<setw(12)<<setprecision(4)<<temp_score[i][1]    //-> tmscore[i]
				       <<setiosflags(ios::left)<<setw(12)<<setprecision(4)<<temp_score[i][2]    //-> RMSD[i]
				       <<setiosflags(ios::left)<<setw(12)<<query_scale
				       <<setiosflags(ios::left)<<setw(12)<<template_scale
				       <<setiosflags(ios::left)<<setw(8)<<aaa[0]                                //-> LALI
				       <<setiosflags(ios::left)<<setw(8)<<aaa[2]                                //-> #tGaps
				       <<setiosflags(ios::left)<<setw(8)<<aaa[1]                                //-> #qGaps
				       <<setiosflags(ios::left)<<setw(8)<<idrate;                               //-> #seqID
				output << endl;
			}
			//output simplied rank
			{
				string pvalu_str=convertDouble(pvalue,4);
				string origi_str=convertDouble(temp_score[i][0],4);
				string tmsco_str=convertDouble(temp_score[i][1],4);
				string rmsd_str=convertDouble(temp_score[i][2],4);
				printf("%s\t%d\t%s\t%s\t%s\t%s\t%d\t%d\n",temp_rank[i].c_str(),aaa[3],
					pvalu_str.c_str(),origi_str.c_str(),tmsco_str.c_str(),rmsd_str.c_str(),aaa[0],idrate);
			}
		}
		output << endl;

		//--- final output alignment ----//
		for(int i=0;i<N;i++)
		{
			//init pvalue check
			double pvalue=0;
			if(pvalue_cutoff>0)
			{
				pvalue=Pvalue(temp_score[i][0],1/beta,miu);
				if(pvalue>pvalue_cutoff)continue;
			}
			//output alignment
			{
				string align_detail_file = align_root + "/" + temp_rank[i] + "-" + seq_name + ".local";
				output<<"No "<<i+1<<endl<<">"<<temp_rank[i]<<endl;
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
		}

		//[2-3-5] output done
		output<<"Done."<<endl;	
		output.close();
	}
}


//--------- main ----------//
int main(int argc,char **argv)
{
	//---- Generate DeepSearch Rank ----//
	{
		if(argc<10)
		{
			fprintf(stderr,"DeepSearch_Rank <query_nam> <rank_sco> <ali_root> <evd_param> <out_name> \n");
			fprintf(stderr,"                <command> <pval_cut> <totn> <topn> \n");
			exit(-1);
		}
		string query_nam=argv[1];
		string rank_sco=argv[2];
		string ali_root=argv[3];
		string evd_param=argv[4];
		string out_name=argv[5];
		string command=argv[6];
		double pval_cut=atof(argv[7]);
		string totn=argv[8];
		int topn=atoi(argv[9]);
		//---- load rank_pval data ----//
		vector <string> temp_rank;
		vector <vector <double> > temp_score;
		int retv=Load_RankSimp_File(rank_sco,temp_rank,temp_score);
		//---- load evd param ----//
		double miu=0;
		double beta=1;
		if(pval_cut>0)Read_Miu_Beta(evd_param,miu,beta);
		//---- determine topn -----//
		int TopN=topn<retv?topn:retv;
		//---- generate detailed rank file ----//
		DeepSearch_RankFile_Generate(out_name,query_nam,ali_root,temp_rank,temp_score,pval_cut,command,miu,beta,totn,TopN);
		//exit
		exit(0);
	}
}


