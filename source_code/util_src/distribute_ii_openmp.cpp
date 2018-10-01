#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#ifdef _OMP
#include <omp.h>
#endif
using namespace std;


//---- get list length ----//
int WS_List_Proc(string &list,vector <string> &proc_rec)
{
	//load
	ifstream fin;
	string buf,temp;
	fin.open(list.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"list not found [%s] !!!\n",list.c_str());
		return -1;
	}
	//process
	int count=0;
	proc_rec.clear();
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		proc_rec.push_back(buf);
		count++;
	}
	return count;
}

//-------- main -------//
int main(int argc, char ** argv) 
{
	// the command line must be:
	if(argc < 2)
	{
		cerr << "./distribute <exe_list_file> \n";
		return 0;
	}
	string exe_list_file=argv[1];

	//load list
	vector <string> exe_list_proc;
	int totnum=WS_List_Proc(exe_list_file,exe_list_proc);
	if(totnum<=0)
	{
		return 0;
	}

	//proc list
	int total_num=0;
	#pragma omp parallel for schedule(dynamic)
	for(int i=0;i<totnum;i++)
	{
		int retv=system(exe_list_proc[i].c_str());
		total_num++;
		fprintf(stderr,"%d/%d\r",total_num,totnum);
	}

	//finalize
	return 0;
}
