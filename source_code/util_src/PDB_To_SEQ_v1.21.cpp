#include <string>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
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
//--------- PDB_To_SEQ ----------//
int PDB_To_SEQ(string &pdb,string &ami,char chain='_')
{
	//--- list for mapping ---//
	map<string, int > ws_mapping;
	map<string, int>::iterator iter;
	ws_mapping.clear();
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
	int count=0;
	char cur_chain;
	char rel_chain;
	int first=1;
	ami="";
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
		if(temp!="ATOM"&&temp!="HETA"&&temp!="MISS")continue;
		//check CA
		temp=buf.substr(13,2);
		if(temp!="CA"&&temp!="  ")continue;
		//chain
		cur_chain=buf[21];
		if(chain!='_')
		{
			if(first==1)
			{
				if(cur_chain!=chain)continue;
				else
				{
					first=0;
					rel_chain=cur_chain;
				}
			}
			else
			{
				if(cur_chain!=rel_chain)break;
			}
		}
		else
		{
			if(first==1)
			{
				first=0;
				rel_chain=cur_chain;
			}
			else
			{
				if(cur_chain!=rel_chain)break;
			}
		}
		//record name
		name=buf.substr(21,6);
		iter = ws_mapping.find(name);
		if(iter != ws_mapping.end())continue;
		count++;
		ws_mapping.insert(map < string, int >::value_type(name, count));
		//final
		temp=buf.substr(17,3);
		char c=WWW_Three2One_III(temp.c_str());
		ami.push_back(c);
	}
	return count;
}

//-------- main ---------//
int main(int argc,char **argv)
{
	//------ PDB_To_SEQ -------//
	{
		if(argc<3)
		{
			fprintf(stderr,"Version: 1.21\n");
			fprintf(stderr,"PDB_To_SEQ <pdb_file> <seq_file>\n");
			exit(-1);
		}
		string pdb_file=argv[1];
		string seq_file=argv[2];
		string ami;
		int ret_val=PDB_To_SEQ(pdb_file,ami);
		if(ret_val<0)
		{
			fprintf(stderr,"file %s bad!\n",pdb_file.c_str());
			exit(-1);
		}
		string nam;
		getBaseName(pdb_file,nam,'/','.');
		FILE *fp=fopen(seq_file.c_str(),"wb");
		fprintf(fp,">%s\n",nam.c_str());
		fprintf(fp,"%s\n",ami.c_str());
		fclose(fp);
		exit(0);
	}
}

