#include <string>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <map>
#include <iomanip>
#include <time.h>
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
void Fasta_Output(string &nam1_full,string &nam2_full,int *wsrec,int size,
		string &nam1,string &nam2,FILE *fp,int ONEorTWO=0)
{
	int i;
	//output
	if(ONEorTWO==0) //cut both
	{
		fprintf(fp,">%s\n",nam1.c_str());
		for(i=0;i<size;i++)
		{
			if(wsrec[i]==1)fprintf(fp,"%c",nam1_full[i]);
		}
		fprintf(fp,"\n");
		fprintf(fp,">%s\n",nam2.c_str());
		for(i=0;i<size;i++)
		{
			if(wsrec[i]==1)fprintf(fp,"%c",nam2_full[i]);
		}
		fprintf(fp,"\n");
	}
	if(ONEorTWO==1) //cut first
	{
		fprintf(fp,">%s\n",nam1.c_str());
		for(i=0;i<size;i++)
		{
			if(wsrec[i]==1 || nam2_full[i]!='-')fprintf(fp,"%c",nam1_full[i]);
		}
		fprintf(fp,"\n");
		fprintf(fp,">%s\n",nam2.c_str());
		for(i=0;i<size;i++)
		{
			if(wsrec[i]==1 || nam2_full[i]!='-')fprintf(fp,"%c",nam2_full[i]);
		}
		fprintf(fp,"\n");
	}
	if(ONEorTWO==2) //cut second
	{
		fprintf(fp,">%s\n",nam1.c_str());
		for(i=0;i<size;i++)
		{
			if(wsrec[i]==1 || nam1_full[i]!='-')fprintf(fp,"%c",nam1_full[i]);
		}
		fprintf(fp,"\n");
		fprintf(fp,">%s\n",nam2.c_str());
		for(i=0;i<size;i++)
		{
			if(wsrec[i]==1 || nam1_full[i]!='-')fprintf(fp,"%c",nam2_full[i]);
		}
		fprintf(fp,"\n");
	}
}


//-------- pick alignment -----//
void WS_Get_Alignment(vector<pair<int, int> > &in,int *wsrec,int *wali1,int *wali2)
{
	int i,j;
	int pos;
	int ii,jj;
	int size=(int)in.size();
	for(i=0;i<size;i++)wsrec[i]=0; //default
	for(i=0;i<size;i++)
	{
		ii=in[i].first;
		jj=in[i].second;	
		if(ii>0 && jj>0)wsrec[i]=1;
		else
		{
			if(ii>0)
			{
				if(wali1[ii-1]!=-2)wsrec[i]=1;
			}
			if(jj>0)
			{
				if(wali2[jj-1]!=-2)wsrec[i]=1;
			}
		}
	}
}
void WS_Pick_Alignment(vector<pair<int, int> > &in, int *wali1, int *wali2 ,
	int moln1,int moln2,char *ami1,char *ami2,
	vector <string> &out1,vector <string> &out2,
	vector <int> &start1,vector <int> &start2,
	int &seq_start1,int &seq_start2,
	int cover,int thres,int SIMPorNOT1,int SIMPorNOT2)
{
	int i;
	int ii,jj;
	int *ali1=new int[moln1];
	int *ali2=new int[moln2];
	for(i=0;i<moln1;i++)ali1[i]=-2;
	for(i=0;i<moln2;i++)ali2[i]=-2;
	for(i=0;i<moln1;i++)wali1[i]=-2;
	for(i=0;i<moln2;i++)wali2[i]=-2;
	int size=(int)in.size();
	for(i=0;i<size;i++)
	{
		ii=in[i].first;
		jj=in[i].second;
		if(ii>0 && jj>0)
		{
			ali1[ii-1]=jj-1;
			ali2[jj-1]=ii-1;
		}
	}
	
	//--- extend ---//
	int j;
	int pos;

//seq1

//if(SIMPorNOT1==0)   // complex
{
	seq_start1=-1;
	for(i=0;i<moln1;i++)
	{
		if(ali1[i]!=-2)
		{
			if(seq_start1<0)
			{
				seq_start1=i-cover;
				if(seq_start1<0)seq_start1=0;
			}
			for(j=0;j<2*cover+1;j++)
			{
				pos=i-cover+j;
				if(pos>=0 && pos<moln1)
				{
					if(wali1[pos]==-1)continue;
					if(ali1[pos]==-2)
					{
						wali1[pos]=-1;
					}
					else
					{
						wali1[pos]=ali1[pos];
					}
				}
			}
		}
	}
}
/*
else   // simple
{
	int start=-1;
	for(i=0;i<moln1;i++)
	{
		if(ali1[i]!=-2)
		{
			start=i-cover;
			break;
		}
	}
	int end=moln1;
	for(i=moln1-1;i>=0;i--)
	{
		if(ali1[i]!=-2)
		{
			end=i+cover;
			break;
		}
	}
	if(start<0)start=0;
	if(end>=moln1)end=moln1-1;
	for(i=start;i<=end;i++)
	{
		wali1[i]=ali1[i];
		if(wali1[i]==-2)wali1[i]=-1;
	}
	seq_start1=start;
}
*/

//seq2
/*
if(SIMPorNOT2==0)   // complex
{
	seq_start2=-1;
	for(i=0;i<moln2;i++)
	{
		if(ali2[i]!=-2)
		{
			if(seq_start2<0)
			{
				seq_start2=i-cover;
				if(seq_start2<0)seq_start2=0;
			}
			for(j=0;j<2*cover+1;j++)
			{
				pos=i-cover+j;
				if(pos>=0 && pos<moln2)
				{
					if(wali2[pos]==-1)continue;
					if(ali2[pos]==-2)
					{
						wali2[pos]=-1;
					}
					else
					{
						wali2[pos]=ali2[pos];
					}
				}
			}
		}
	}
}
else   // simple
*/
{
	int start=-1;
	for(i=0;i<moln2;i++)
	{
		if(ali2[i]!=-2)
		{
			start=i-cover;
			break;
		}
	}
	int end=moln2;
	for(i=moln2-1;i>=0;i--)
	{
		if(ali2[i]!=-2)
		{
			end=i+cover;
			break;
		}
	}
	if(start<15)start=0;           //-> modified by WS //__130430__//
	if(end>=moln2-15)end=moln2-1;  //-> modified by WS //__130430__//
	for(i=start;i<=end;i++)
	{
		wali2[i]=ali2[i];
		if(wali2[i]==-2)wali2[i]=-1;
	}
	seq_start2=start;
}

	//--- extract ---//
	int first;
	int start,len;
	// rec1 ---//
	first=1;
	start=0;
	len=0;
	vector <pair<int, int> > record1;
	int wsmax1=0;
	for(i=0;i<moln1;i++)
	{
		if(wali1[i]==-2)
		{
			if(first==1)
			{
				first=0;
				start=i;
				len=0;
			}
			len++;
		}
		else
		{
			if(first==0)
			{
				record1.push_back(pair<int,int>(start,len));
				if(len>wsmax1)wsmax1=len;
			}
			first=1;
		}
	}
	if(first==0)
	{
		record1.push_back(pair<int,int>(start,len));
		if(len>wsmax1)wsmax1=len;
	}
	// rec2 ---//
	first=1;
	start=0;
	len=0;
	vector <pair<int, int> > record2;
	int wsmax2=0;
	for(i=0;i<moln2;i++)
	{
		if(wali2[i]==-2)
		{
			if(first==1)
			{
				first=0;
				start=i;
				len=0;
			}
			len++;
		}
		else
		{
			if(first==0)
			{
				record2.push_back(pair<int,int>(start,len));
				if(len>wsmax2)wsmax2=len;
			}
			first=1;
		}
	}
	if(first==0)
	{
		record2.push_back(pair<int,int>(start,len));
		if(len>wsmax2)wsmax2=len;
	}



	//--- output ---//
	char *tmp1=new char[moln1+1];
	char *tmp2=new char[moln2+1];
	out1.clear();
	size=record1.size();
	for(i=0;i<size;i++)
	{
		start=record1[i].first;
		len=record1[i].second;
		if(len>thres)
		{
			strncpy(tmp1,ami1+start,len);
			tmp1[len]='\0';
			string temp=tmp1;
			out1.push_back(temp);
			start1.push_back(start);
		}
		else
		{
//			for(j=0;j<len;j++)wali1[start+j]=-1;
		}
	}
	out2.clear();
	size=record2.size();
	for(i=0;i<size;i++)
	{
		start=record2[i].first;
		len=record2[i].second;
		if(len>thres)
		{
			strncpy(tmp2,ami2+start,len);
			tmp2[len]='\0';
			string temp=tmp2;
			out2.push_back(temp);
			start2.push_back(start);
		}
		else
		{
//			for(j=0;j<len;j++)wali2[start+j]=-1;
		}
	}
	//--- final delete ---//
	delete [] tmp1;
	delete [] tmp2;
	delete [] ali1;
	delete [] ali2;
}

//------- main proc -------//
void WS_Main_Proc(string &fasta_file,int cover,int thres,string &out_file,int dom_num=1,int ONEorTWO=0,
	int SIMP1=0,int SIMP2=0)
{
	//---- load alignment ----//
	vector<pair<int, int> > alignment;
	string nam1_content,nam2_content;
	string nam1_full_,nam2_full_;
	string nam1,nam2;
	ReadToFile_FASTA(fasta_file,alignment,nam1_content,nam2_content,nam1_full_,nam2_full_,nam1,nam2);
	//--- kill double gap ----//
	string nam1_full,nam2_full;
	int i,size;
	int lali=0;
	size=(int)nam1_full_.length();
	for(i=0;i<size;i++)
	{
		if(nam1_full_[i]=='-' && nam2_full_[i]=='-')continue;
		nam1_full.push_back(nam1_full_[i]);
		nam2_full.push_back(nam2_full_[i]);
		lali++;
	}

	//--- get ini length --//
	int moln1=(int)nam1_content.length();
	int moln2=(int)nam2_content.length();
	char *ami1=new char[moln1+1];
	char *ami2=new char[moln2+2];
	strcpy(ami1,nam1_content.c_str());
	strcpy(ami2,nam2_content.c_str());
	//---- proc domain -----//
	vector <string> out1;
	vector <string> out2;
	vector <int> start1;
	vector <int> start2;
	int seq_start1,seq_start2;
	int *wali1=new int[moln1];
	int *wali2=new int[moln2];
	WS_Pick_Alignment(alignment,wali1,wali2,moln1,moln2,ami1,ami2,out1,out2,start1,start2,seq_start1,seq_start2,cover,thres,SIMP1,SIMP2);

	//--- get alignment --//
	int totnum=(int)alignment.size();
	int *wsrec=new int[totnum];
	WS_Get_Alignment(alignment,wsrec,wali1,wali2);

	//--- output ----//
	FILE *fp;
	FILE *fq;
	string wsout;
	wsout=out_file+"_list1";
	fp=fopen(wsout.c_str(),"wb");
	wsout=out_file+"_start1";
	fq=fopen(wsout.c_str(),"wb");
	size=(int)out1.size();
	for(i=0;i<size;i++)
	{
		fprintf(fp,">%s_%d\n",nam1.c_str(),dom_num+i);
		fprintf(fp,"%s\n",out1[i].c_str());
		fprintf(fq,"%d\n",start1[i]);
	}
	fclose(fp);
	fclose(fq);
	wsout=out_file+"_list2";
	fp=fopen(wsout.c_str(),"wb");
	wsout=out_file+"_start2";
	fq=fopen(wsout.c_str(),"wb");
	size=(int)out2.size();
	for(i=0;i<size;i++)
	{
		fprintf(fp,">%s_%d\n",nam2.c_str(),dom_num+i);
		fprintf(fp,"%s\n",out2[i].c_str());
		fprintf(fq,"%d\n",start2[i]);
	}
	fclose(fp);
	fclose(fq);
	
	//--- output more ---//
	fp=fopen(out_file.c_str(),"wb");
	Fasta_Output(nam1_full,nam2_full,wsrec,totnum,nam1,nam2,fp,ONEorTWO);
	fclose(fp);

	//--- output single sequence ---//
	string wsfasta1,wsfasta2;
	if(ONEorTWO==0)  //cut both
	{
		int strcount;
		size=totnum;
		wsfasta1="";
		strcount=0;
		for(i=0;i<size;i++)
		{
			if(nam1_full[i]!='-' && wsrec[i]==1)
			{
				wsfasta1=wsfasta1+nam1_full[i];
				strcount++;
			}
		}
		wsfasta1[strcount]='\0';
		wsfasta2="";
		strcount=0;
		for(i=0;i<size;i++)
		{
			if(nam2_full[i]!='-' && wsrec[i]==1)
			{
				wsfasta2=wsfasta2+nam2_full[i];
				strcount++;
			}
		}
		wsfasta2[strcount]='\0';
	}
	if(ONEorTWO==1)  //cut first
	{
		int strcount;
		size=totnum;
		wsfasta1="";
		strcount=0;
		for(i=0;i<size;i++)
		{
			if(nam1_full[i]!='-' && wsrec[i]==1)
			{
				wsfasta1=wsfasta1+nam1_full[i];
				strcount++;
			}
		}
		wsfasta1[strcount]='\0';
		wsfasta2=nam2_content;
	}
	if(ONEorTWO==2)  //cut second
	{
		int strcount;
		size=totnum;
		wsfasta2="";
		strcount=0;
		for(i=0;i<size;i++)
		{
			if(nam2_full[i]!='-' && wsrec[i]==1)
			{
				wsfasta2=wsfasta2+nam2_full[i];
				strcount++;
			}
		}
		wsfasta2[strcount]='\0';
		wsfasta1=nam1_content;
	}

	//-- output --//
	wsout=out_file+"_seq1";
	fp=fopen(wsout.c_str(),"wb");
	wsout=out_file+"_seqstart1";
	fq=fopen(wsout.c_str(),"wb");
//	fprintf(fp,">%s_seq\n",nam1.c_str());
	fprintf(fp,"%s\n",wsfasta1.c_str());
	fprintf(fq,"%d\n",seq_start1);
	fclose(fp);
	fclose(fq);
	wsout=out_file+"_seq2";
	fp=fopen(wsout.c_str(),"wb");
	wsout=out_file+"_seqstart2";
	fq=fopen(wsout.c_str(),"wb");
//	fprintf(fp,">%s_seq\n",nam2.c_str());
	fprintf(fp,"%s\n",wsfasta2.c_str());
	fprintf(fq,"%d\n",seq_start2);
	fclose(fp);
	fclose(fq);


	//--- delete --//
	delete [] ami1;
	delete [] ami2;
	delete [] wali1;
	delete [] wali2;
	delete [] wsrec;
}



//-------- main --------//
//-> we default second is sequence, and domain should be coutinous
int main(int argc,char **argv)
{
	//---- domain proc -------//
	{
		if(argc<7)
		{
			printf("Version: 1.41 \n");
			printf("Domain_Proc <fasta_file> <out_file> <cover> <thres> <dom_num> <ONEorTWO> \n");
			printf("cover should be set to 5-10 \n");
			printf("thres should be set to 200 \n");
			printf("dom_num should be set to 1 \n");
			printf("[ONEorTWO]: 0 for both, 1 for first, 2 for second \n");
			exit(-1);
		}
		string fasta_file=argv[1];
		string out_file=argv[2];
		int cover=atoi(argv[3]);
		int thres=atoi(argv[4]);
		int dom_num=atoi(argv[5]);
		int ONEorTWO=atoi(argv[6]);
		int SIMP1=0;    //-> template use complex domain cut
		int SIMP2=1;    //-> sequence use simple domain cut
		WS_Main_Proc(fasta_file,cover,thres,out_file,dom_num,ONEorTWO,SIMP1,SIMP2);
		exit(0);
	}
}
