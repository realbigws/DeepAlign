#include "PDB_Chain_Fold.h"
#include "Confo_Back.h"
#include "Confo_Beta.h"
#include "Confo_Lett.h"
#include "Acc_Surface.h"
#include "Mol_File.h"
#include "Mol_Out.h"
#include "getopt.h"
using namespace std;


//-----------------------------------------------------------------------------------------------------------//
//---- print_help_msg => print help message
void print_help_msg(void) 
{
	cout << "========================================================|" << endl;
	cout << "PDB_Tool  (version 4.7) [2014.05.30]                    |" << endl;
	cout << "          Transform original .PDB file to compact form  |" << endl;
	cout << "Usage:   ./PDB_Tool <-i input> <-r range> <-o output>   |" << endl;
	cout << "Or,      ./PDB_Tool <-L list>                           |" << endl;
	cout << "--------------------------------------------------------|" << endl;
	cout << "For <-i input> input type,                              |" << endl;
	cout << "-i input  : input original PDB file, e.g., 1col.pdb     |" << endl;
	cout << "-r range  : input residue range, e.g., A:1-51           |" << endl;
	cout << "-o output : output PDB file, e.g., 1colA.pdb            |" << endl;
	cout << "--------------------------------------------------------|" << endl;
	cout << "For <-L list> input type,                               |" << endl;
	cout << "The 1st line in <list> specifies the input directory.   |" << endl;
	cout << "The 2nd line in <list> specifies the output directory.  |" << endl;
	cout << "The other lines in <list> indicate [input range output] |" << endl;
	cout << "          e.g.,  1col.pdb A:1-51 1colA.pdb              |" << endl;
	cout << "========================================================|" << endl;
	cout << "                Primary options                         |" << endl;
	cout << "--------------------------------------------------------|" << endl;
	cout << "-M mode : Specify the atoms for output.                 |" << endl;
	cout << "          -2,  output CA+CB                             |" << endl;
	cout << "          -1,  output CA                                |" << endl;
	cout << "           0,  output backbone+CB                       |" << endl;
	cout << "          [1], output full atoms (Set as default)       |" << endl;
	cout << "--------------------------------------------------------|" << endl;
	cout << "-N num  : Specify the numbering type for output.        |" << endl;
	cout << "          -1,  full sequential, for atom and residue    |" << endl;
	cout << "           0,  part sequential, for atom only           |" << endl;
	cout << "          [1], PDB numbering (Set as default)           |" << endl;
	cout << "--------------------------------------------------------|" << endl;
	cout << "                Additional options                      |" << endl;
	cout << "--------------------------------------------------------|" << endl;
	cout << "The following arguments for both input types            |" << endl;
	cout << "        -C 1 to consider non-CA atoms                   |" << endl;
	cout << "        -R 1 to reconstruct missing CB and backbone     |" << endl;
	cout << "The following arguments only for <-L list> input type   |" << endl;
	cout << "        -F 1 to output AMI,CLE,SSE, and ACC file        |" << endl;
	cout << "        -G 1 to output three log files                  |" << endl;
	cout << "========================================================|" << endl;
	exit(-1);
}
//---- WebServer's default input ----//
string INPUT_NAM="";
string INPUT_RAN="_";
string INPUT_OUT="";
string INPUT_LIST="";
int LIST_OR_SINGLE=0; //default: single
int INPUT_MODE=1; //main (default:1)
int INPUT_TYPE=1; //main (default:1)
int INPUT_GLYS=1; //vice (default:1)
int INPUT_NOCA=0; //vice (default:0)
int INPUT_RECO=0; //vice (default:0)
int INPUT_FIFI=0; //vice (default:0)
int INPUT_LOGF=0; //vice (default:0)
int WARN_OUT=1;   //vice (default:1)

//-----------------------------------------------------------------------------------------------------------//
//---- parameter editor ----//
static option long_options[] =
{
	{"input",   no_argument,       NULL, 'i'},
	{"range",   no_argument,       NULL, 'r'},
	{"output",  no_argument,       NULL, 'o'},
	{"List",    no_argument,       NULL, 'L'},
	{"Mode",    no_argument,       NULL, 'M'},
	{"Type",    no_argument,       NULL, 'N'},
	{"Glys",    no_argument,       NULL, 'T'},
	{"Noca",    no_argument,       NULL, 'C'},
	{"Reco",    no_argument,       NULL, 'R'},
	{"Fifi",    no_argument,       NULL, 'F'},
	{"Logf",    no_argument,       NULL, 'G'},
	{"Warn",    no_argument,       NULL, 'W'},
	{0, 0, 0, 0}
};
//-----------------------------------------------------------------------------------------------------------//
//---- process_args => process input parameter args
void process_args(int argc,char** argv) 
{
	string buf;
	int opt;
	if(1==argc)print_help_msg();    
	while(true) 
	{
		int option_index=0;
		opt=getopt_long(argc,argv,"i:r:o:L:M:N:T:C:R:F:G:W:",
			   long_options,&option_index);
		if (opt==-1)break;	
		switch(opt) 
		{
			case 'i':
				INPUT_NAM=optarg;
				break;
			case 'r':
				INPUT_RAN=optarg;
				break;
			case 'o':
				INPUT_OUT=optarg;
				break;
			case 'L':
				INPUT_LIST=optarg;
				LIST_OR_SINGLE++;
				break;
			case 'M':
				INPUT_MODE=atoi(optarg);
				break;
			case 'N':
				INPUT_TYPE=atoi(optarg);
				break;
			case 'T':
				INPUT_GLYS=atoi(optarg);
				break;
			case 'C':
				INPUT_NOCA=atoi(optarg);
				break;
			case 'R':
				INPUT_RECO=atoi(optarg);
				break;
			case 'F':
				INPUT_FIFI=atoi(optarg);
				break;
			case 'G':
				INPUT_LOGF=atoi(optarg);
				break;
			case 'W':
				WARN_OUT=atoi(optarg);
				break;
			default:
				exit(-1);
		}
	}
}

//------------- Get_PDB_File_Len -----------//
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

//==================== WS_PDB_Back_Process ===============// (process list)
//[list_style]
//first_line: input_dir
//second_line: output_dir
//others:
//input range output-> e.g., 1col.pdb A:1-51 1colA.res
//[note] -> should only contain two  ' '
int WS_PDB_Back_Process(string &list,int OutType,int OutMode,int OutGlys,int OutNoca,int OutReco,int OutFifi,int OutLogf,int OutWarn)
{
	//class
	Mol_File mol_input;
	Mol_Out mol_output;
	Confo_Beta confo_beta;
	Confo_Back confo_back;
	Confo_Lett confo_lett;
	Acc_Surface acc_surface;
	//init
	ifstream fin;
	string buf;
	string input,chain,output;
	string file,temp;
	string path,outa;
	char c;
	int wwscount=0;
	int i;
	int len;
	int pos[10];
	int wscount;
	int ret_val;
	int totlen=30000;
	int moln;
	PDB_Residue *pdb=new PDB_Residue[totlen];
	FILE *fp=0;
	FILE *fq=0;
	if(OutLogf==1)
	{
		fp=fopen("ws_real_list","wb");
		fq=fopen("ws_valid_list","wb");
		//judge
		if(fp==0 || fq==0)
		{
			fprintf(stderr,"ERROR: logfile can't be opened. \n");
			OutLogf=0;
		}
	}
	string TER="TER                                                                             ";


	//recon_related
	XYZ *mol=new XYZ[totlen];   //CA
	XYZ *mcb=new XYZ[totlen];   //CB
	XYZ **mbb;                //BackBone (N,CA,C,O,CB)
	XYZ **mcc;                //BackBone (N,CA,C,O,CB,...,)
	NewArray2D(&mbb,totlen,5);
	NewArray2D(&mcc,totlen,15);
	int *mcc_side=new int[totlen];
	char *ami=new char[totlen+1];
	char *cle=new char[totlen+1];
	char *acc=new char[totlen+1];


	//open
	fin.open(list.c_str(),ios::in);
	if(fin.fail()!=0)
	{
		printf("Process_List Not Found!![%s]\n",list.c_str());
		fin.close();
		fin.clear();
		return -1;
	}
	if(OutLogf==1)
	{
		mol_input.flog=fopen("ws_pdb_log","wb");
		mol_input.LOGOUT=1;
	}
	mol_input.CbBACK=1;
	if(OutNoca==1)mol_input.CaONLY=0; //consider Non-CA atoms !!//__110408__//
	//macro
	mol_input.OUTP_MODE=OutType;
	mol_input.PROC_MODE=OutMode;
	mol_input.GLYC_MODE=OutGlys;
	//memory limit
	mol_input.WARNING_out=OutWarn;
	mol_input.MEMORY_LIMIT=totlen;


	//---just_temp---//[get path]
	{
		if(!getline(fin,buf,'\n'))
		{
			printf("NO_INPUT_DIRECTORY[%s]!!!\n",buf.c_str());
			exit(-1);
		}
		path=buf;
		if(!getline(fin,buf,'\n'))
		{
			printf("NO_OUTPUT_DIRECTORY[%s]!!!\n",buf.c_str());
			exit(-1);
		}
		outa=buf;
	}
	//process
	int zero=-1;
	int count=0;
	for(;;)
	{
		//read
		if(!getline(fin,buf,'\n'))break;
		getline_end(buf,0x0D);
		//check
		len=(int)buf.length();
		wscount=0;
		for(i=0;i<len;i++)
		{
			if(buf[i]==' ')
			{
				pos[wscount]=i;
				wscount++;
			}
		}
		if(wscount!=2)
		{
			printf("LIST_Format ERROR!!![%s]\n",buf.c_str());
			exit(-2);
		}
		//calc
		input=buf.substr(0,pos[0]);
		chain=buf.substr(pos[0]+1,pos[1]-pos[0]-1);
		c=chain[0];
		output=buf.substr(pos[1]+1,len-pos[1]-1);
		file=path+'/'+input;
		temp="0 "+file+' '+chain;


		//pre_process
		ret_val=mol_input.XYZ_Tranform(temp,moln,0,0,0,0,0,0);
		if(ret_val<=0)
		{
			if(ret_val!=-12345)
			{
				if(OutLogf==1)fprintf(fq,"%s @ %5d\n",buf.c_str(),zero);
				continue;
			}
		}
		//check memory
		if(ret_val==-12345)
		{
			//add memory
			if(moln>totlen)
			{
				totlen=moln;
				//delete
				delete [] pdb;
				delete [] mol;
				delete [] mcb;
				delete [] ami;
				delete [] cle;
				delete [] acc;
				DeleteArray2D(&mbb,totlen);
				DeleteArray2D(&mcc,totlen);
				delete [] mcc_side;
				//create
				pdb=new PDB_Residue[totlen];
				mol=new XYZ[totlen];
				mcb=new XYZ[totlen];
				ami=new char[totlen+1];
				cle=new char[totlen+1];
				acc=new char[totlen+1];
				NewArray2D(&mbb,totlen,5);
				NewArray2D(&mcc,totlen,15);
				mcc_side=new int[totlen];
				//memory limit
				mol_input.MEMORY_LIMIT=totlen;
			}
		}

		//reload
		int PRE_LOAD_=mol_input.PRE_LOAD;
		int WARNING_out_=mol_input.WARNING_out;
		mol_input.PRE_LOAD=1;
		mol_input.WARNING_out=0;
		mol_input.XYZ_Tranform(temp,moln,0,mol,ami,0,0,pdb);
		mol_input.PRE_LOAD=PRE_LOAD_;
		mol_input.WARNING_out=WARNING_out_;


		//check
		{
			int i,j;
			int correct;
			int iret;
			correct=1; //default:OK
			for(i=0;i<moln;i++)
			{
				iret=pdb[i].PDB_residue_backbone_check(4);
				if(iret!=1)
				{
					correct=0;
					break;
				}
				iret=pdb[i].PDB_residue_CB_check();
				if(iret!=1)
				{
					correct=0;
					break;
				}
			}

			//judge
			if(correct==1)
			{
				wwscount++;
				if(OutLogf==1)
				{
					fprintf(fp,"%s\n",output.c_str());
					fprintf(fq,"%s 1 %5d\n",buf.c_str(),moln);
				}
			}
			else
			{
				wwscount++;
				if(OutLogf==1)
				{
					fprintf(fp,"%s\n",output.c_str());
					fprintf(fq,"%s 0 %5d\n",buf.c_str(),moln);
				}
				if(OutReco==1)
				{
					//[1]recon
					confo_lett.btb_ori(0,0,0,moln,mol,cle);
					confo_back.Recon_Back_WS_Main(mol,cle,moln,mbb);      //given CA, recon BackBone (N,CA,C,O,CB)
					confo_beta.WS_Recon_Beta_21(mol,mcb,moln,ami,cle);    //given CA, recon CB
					//[2]assign
					for(i=0;i<moln;i++)
					{
						//CA missing!!
						if(pdb[i].get_backbone_part_index(1)==0)
						{
							//backbone (N,CA,C,O)
							for(j=0;j<4;j++)pdb[i].set_backbone_atom(j,mbb[i][j]);
							//sidechain (CB)
							if(i==0||i==moln-1)pdb[i].set_sidechain_atom(0,mbb[i][4]);
							else pdb[i].set_sidechain_atom(0,mcb[i]);
						}
						else
						{
							//backbone (N,CA,C,O)
							for(j=0;j<4;j++)
							{
								if(pdb[i].get_backbone_part_index(j)==0)
								{
									pdb[i].set_backbone_atom(j,mbb[i][j]);
								}
							}
							//sidechain (CB)
							if(pdb[i].get_sidechain_part_index(0)==0)
							{
								if(i==0||i==moln-1)pdb[i].set_sidechain_atom(0,mbb[i][4]);
								else pdb[i].set_sidechain_atom(0,mcb[i]);
							}
						}
					}
				}
			}

			//output
			FILE *fpdb;
			file="";
			file=file+outa+"/"+output+".pdb";
			fpdb=fopen(file.c_str(),"wb");
			if(fpdb==0)
			{
				fprintf(stderr,"ERROR: file %s can't be opened. \n",file.c_str());
			}
			else
			{
				mol_output.Output_PDB_III(fpdb,moln,pdb,c,OutType,OutMode,OutGlys);
//				fprintf(fpdb,"%s\n",TER.c_str());
				fclose(fpdb);
			}

			//output others
			if(OutFifi==1)
			{
				//init
				FILE *fpp;
				PDB_Chain_Fold chain_fold;
				chain_fold.initialize_simple(moln,' ');
				for(i=0;i<moln;i++)chain_fold.set_residue(i,pdb[i]);
				int retv;
				//calc
				retv=chain_fold.calculate_CLE();
				if(retv!=0)
				{
					printf("[%s]CLE_BAD!!!\n",buf.c_str());
					continue;
				}
				retv=chain_fold.calculate_SSE();
				if(retv!=0)
				{
					printf("[%s]SSE_BAD!!!\n",buf.c_str());
					continue;
				}
				//output AMI,CLE,SSE
				string AMI,CLE,SSE;
				AMI=chain_fold.get_sequence();
				CLE=chain_fold.get_CLE();
				SSE=chain_fold.get_SSE();
				file="";
				file=file+outa+"/"+output+".ami";
				fpp=fopen(file.c_str(),"wb");
				if(fpp==0)
				{
					fprintf(stderr,"ERROR: file %s can't be opened. \n",file.c_str());
				}
				else
				{
					fprintf(fpp,">%s\n",output.c_str());
					fprintf(fpp,"%s\n",AMI.c_str());
					fclose(fpp);
				}
				file="";
				file=file+outa+"/"+output+".cle";
				fpp=fopen(file.c_str(),"wb");
				if(fpp==0)
				{
					fprintf(stderr,"ERROR: file %s can't be opened. \n",file.c_str());
				}
				else
				{
					fprintf(fpp,">%s\n",output.c_str());
					fprintf(fpp,"%s\n",CLE.c_str());
					fclose(fpp);
				}
				file="";
				file=file+outa+"/"+output+".sse";
				fpp=fopen(file.c_str(),"wb");
				if(fpp==0)
				{
					fprintf(stderr,"ERROR: file %s can't be opened. \n",file.c_str());
				}
				else
				{
					fprintf(fpp,">%s\n",output.c_str());
					fprintf(fpp,"%s\n",SSE.c_str());
					fclose(fpp);
				}
			
				//-> output ACC_Code 
				for(i=0;i<moln;i++)pdb[i].get_XYZ_array(mcc[i],mcc_side[i]);
				acc_surface.AC_Calc_SolvAcc(mcc,ami,moln,acc,mcc_side);
				file="";
				file=file+outa+"/"+output+".acc";
				fpp=fopen(file.c_str(),"wb");
				if(fpp==0)
				{
					fprintf(stderr,"ERROR: file %s can't be opened. \n",file.c_str());
				}
				else
				{
					fprintf(fpp,">%s\n",output.c_str());
					fprintf(fpp,"%s\n",acc);
					fclose(fpp);
				}
			}
		}

		//printf
		count++;
//		printf("now[%d]\r",count);
	}
	//terminal
	delete [] pdb;
	delete [] ami;
	delete [] cle;
	delete [] acc;
	delete [] mol;
	delete [] mcb;
	DeleteArray2D(&mbb,totlen);
	DeleteArray2D(&mcc,totlen);
	delete [] mcc_side;
	fin.close();
	fin.clear();
	return wwscount;
}

//===================== single_process ===================//__110710__//
//[single_style]
void WS_PDB_Back_Process_Single(string &input,string &range,string &output,
	int OutType,int OutMode,int OutGlys,int OutNoca,int OutReco,int OutFifi,int OutLogf,int OutWarn)
{
	//class
	Mol_File mol_input;
	Mol_Out mol_output;
	Confo_Beta confo_beta;
	Confo_Back confo_back;
	Confo_Lett confo_lett;
	mol_input.MODRES=1;
	//init
	int ret_val;
	int totlen=Get_PDB_File_Len(input);
	if(totlen<=0)
	{
		fprintf(stderr,"pdbfile %s length error!!\n",input.c_str());
		exit(-1);
	}
	int moln;
	PDB_Residue *pdb=new PDB_Residue[totlen];
	string TER="TER                                                                             ";

	//recon_related
	XYZ *mol=new XYZ[totlen];   //CA
	XYZ *mcb=new XYZ[totlen];   //CB
	XYZ **mbb;                  //BackBone (N,CA,C,O,CB)
	NewArray2D(&mbb,totlen,5);
	char *ami=new char[totlen+1];
	char *cle=new char[totlen+1];

	//open
	mol_input.CbBACK=1;
	if(OutNoca==1)mol_input.CaONLY=0; //consider Non-CA atoms !!//__110408__//
	//macro
	mol_input.OUTP_MODE=OutType;
	mol_input.PROC_MODE=OutMode;
	mol_input.GLYC_MODE=OutGlys;
	//memory limit
	mol_input.WARNING_out=OutWarn;
	mol_input.MEMORY_LIMIT=totlen;

	//process
//	int zero=-1;
//	for(;;)
	{
		//pre_process
		ret_val=mol_input.XYZ_Input(input,range,0,moln,0,0,0,0,0);
		if(ret_val<=0)
		{
			if(ret_val!=-12345)goto end;
		}
		if(ret_val!=1)goto end;
		//check memory
		if(ret_val==-12345)
		{
			//add memory
			if(moln>totlen)
			{
				totlen=moln;
				//delete
				delete [] pdb;
				delete [] mol;
				delete [] mcb;
				delete [] ami;
				delete [] cle;
				DeleteArray2D(&mbb,totlen);
				//create
				pdb=new PDB_Residue[totlen];
				mol=new XYZ[totlen];
				mcb=new XYZ[totlen];
				ami=new char[totlen+1];
				cle=new char[totlen+1];
				NewArray2D(&mbb,totlen,5);
				//memory limit
				mol_input.MEMORY_LIMIT=totlen;
			}
		}
		//reload
		int PRE_LOAD_=mol_input.PRE_LOAD;
		int WARNING_out_=mol_input.WARNING_out;
		mol_input.PRE_LOAD=1;
		mol_input.WARNING_out=0;
		mol_input.XYZ_Input(input,range,0,moln,mol,ami,0,0,pdb);
		mol_input.PRE_LOAD=PRE_LOAD_;
		mol_input.WARNING_out=WARNING_out_;

		//check
		{
			int i,j;
			int correct;
			int iret;
			correct=1; //default:OK
			for(i=0;i<moln;i++)
			{
				iret=pdb[i].PDB_residue_backbone_check(4);
				if(iret!=1)
				{
					correct=0;
					break;
				}
				iret=pdb[i].PDB_residue_CB_check();
				if(iret!=1)
				{
					correct=0;
					break;
				}
			}
			confo_lett.btb_ori(0,0,0,moln,mol,cle);

			//judge
			if(correct!=1)
			{
				if(OutReco==1)
				{
					//[1]recon
					confo_back.Recon_Back_WS_Main(mol,cle,moln,mbb);      //given CA, recon BackBone (N,CA,C,O,CB)
					confo_beta.WS_Recon_Beta_21(mol,mcb,moln,ami,cle);    //given CA, recon CB
					//[2]assign
					for(i=0;i<moln;i++)
					{
						//CA missing!!
						if(pdb[i].get_backbone_part_index(1)==0)
						{
							//backbone (N,CA,C,O)
							for(j=0;j<4;j++)pdb[i].set_backbone_atom(j,mbb[i][j]);
							//sidechain (CB)
							if(i==0||i==moln-1)pdb[i].set_sidechain_atom(0,mbb[i][4]);
							else pdb[i].set_sidechain_atom(0,mcb[i]);
						}
						else
						{
							//backbone (N,CA,C,O)
							for(j=0;j<4;j++)
							{
								if(pdb[i].get_backbone_part_index(j)==0)
								{
									pdb[i].set_backbone_atom(j,mbb[i][j]);
								}
							}
							//sidechain (CB)
							if(pdb[i].get_sidechain_part_index(0)==0)
							{
								if(i==0||i==moln-1)pdb[i].set_sidechain_atom(0,mbb[i][4]);
								else pdb[i].set_sidechain_atom(0,mcb[i]);
							}
						}
					}
				}
			}

			//output
			FILE *fpdb;
			fpdb=fopen(output.c_str(),"wb");
			if(fpdb==0)
			{
				fprintf(stderr,"ERROR: file %s can't be opened. \n",output.c_str());
			}
			else
			{
				mol_output.Output_PDB_III(fpdb,moln,pdb,'_',OutType,OutMode,OutGlys);
//				fprintf(fpdb,"%s\n",TER.c_str());
				fclose(fpdb);
			}
		}
	}
end:
	//terminal
	delete [] pdb;
	delete [] ami;
	delete [] cle;
	delete [] mol;
	delete [] mcb;
	DeleteArray2D(&mbb,totlen);
}

//============== main ===============//
int main(int argc, char** argv)
{
	//---- BackBone ----//(process list or single) PDB_Backbone//
	{
		process_args(argc,argv);
		if(LIST_OR_SINGLE>0)WS_PDB_Back_Process(INPUT_LIST,INPUT_TYPE,INPUT_MODE,INPUT_GLYS,INPUT_NOCA,INPUT_RECO,INPUT_FIFI,INPUT_LOGF,WARN_OUT);
		else WS_PDB_Back_Process_Single(INPUT_NAM,INPUT_RAN,INPUT_OUT,INPUT_TYPE,INPUT_MODE,INPUT_GLYS,INPUT_NOCA,INPUT_RECO,INPUT_FIFI,INPUT_LOGF,WARN_OUT);
		exit(0);
	}
}
