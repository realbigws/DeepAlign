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
	cout << "PDB_Tool  (version 4.72) [2015.07.20]                   |" << endl;
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
	cout << "        -F 1 for AMI,CLE,SSE; 2 for ACC; 4 for FEAT     |" << endl;
	cout << "           these output files could be combined         |" << endl;
	cout << "The following arguments only for <-L list> input type   |" << endl;
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

//------ output protein feature files -------//2015_02_20//
//--> output AMI_SSE_CLE
void Output_Protein_Features_AMI_SSE_CLE(
	string &outroot,string &outname,int moln,PDB_Residue *pdb)
{
	//init
	int i;
	FILE *fpp;
	string file;
	PDB_Chain_Fold chain_fold;
	chain_fold.initialize_simple(moln,' ');
	for(i=0;i<moln;i++)chain_fold.set_residue(i,pdb[i]);
	int retv;
	//calc
	retv=chain_fold.calculate_CLE();
	if(retv!=0)
	{
		fprintf(stderr,"[%s]CLE_BAD!!!\n",outname.c_str());
		exit(-1);
	}
	retv=chain_fold.calculate_SSE();
	if(retv!=0)
	{
		fprintf(stderr,"[%s]SSE_BAD!!!\n",outname.c_str());
		exit(-1);
	}
	//------ output AMI,CLE,SSE ------//
	string AMI,CLE,SSE;
	AMI=chain_fold.get_sequence();
	CLE=chain_fold.get_CLE();
	SSE=chain_fold.get_SSE();
	//output AMI
	file="";
	file=file+outroot+"/"+outname+".ami";
	fpp=fopen(file.c_str(),"wb");
	if(fpp==0)
	{
		fprintf(stderr,"ERROR: file %s can't be opened. \n",file.c_str());
	}
	else
	{
		fprintf(fpp,">%s\n",outname.c_str());
		fprintf(fpp,"%s\n",AMI.c_str());
		fclose(fpp);
	}
	//output CLE
	file="";
	file=file+outroot+"/"+outname+".cle";
	fpp=fopen(file.c_str(),"wb");
	if(fpp==0)
	{
		fprintf(stderr,"ERROR: file %s can't be opened. \n",file.c_str());
	}
	else
	{
		fprintf(fpp,">%s\n",outname.c_str());
		fprintf(fpp,"%s\n",CLE.c_str());
		fclose(fpp);
	}
	//output SSE
	file="";
	file=file+outroot+"/"+outname+".sse";
	fpp=fopen(file.c_str(),"wb");
	if(fpp==0)
	{
		fprintf(stderr,"ERROR: file %s can't be opened. \n",file.c_str());
	}
	else
	{
		fprintf(fpp,">%s\n",outname.c_str());
		fprintf(fpp,"%s\n",SSE.c_str());
		fclose(fpp);
	}
}

//--> output ACC and ACC_Value
void Output_Protein_Features_ACC(
	string &outroot,string &outname,Acc_Surface &acc_surface,
	int moln,PDB_Residue *pdb,XYZ **mcc,int *mcc_side,char *ami,char *acc)
{
	//init
	int i;
	FILE *fpp;
	string file;
	PDB_Chain_Fold chain_fold;
	chain_fold.initialize_simple(moln,' ');
	for(i=0;i<moln;i++)chain_fold.set_residue(i,pdb[i]);
	//------ output ACC_Code and ACC_Value ------//
	for(i=0;i<moln;i++)pdb[i].get_XYZ_array(mcc[i],mcc_side[i]);
	acc_surface.AC_Calc_SolvAcc(mcc,ami,moln,acc,mcc_side);
	//output ACC_Code
	file="";
	file=file+outroot+"/"+outname+".acc";
	fpp=fopen(file.c_str(),"wb");
	if(fpp==0)
	{
		fprintf(stderr,"ERROR: file %s can't be opened. \n",file.c_str());
	}
	else
	{
		fprintf(fpp,">%s\n",outname.c_str());
		fprintf(fpp,"%s\n",acc);
		fclose(fpp);
	}
	//output ACC_Value
	file="";
	file=file+outroot+"/"+outname+".acc_value";
	fpp=fopen(file.c_str(),"wb");
	if(fpp==0)
	{
		fprintf(stderr,"ERROR: file %s can't be opened. \n",file.c_str());
	}
	else
	{
		for(i=0;i<moln;i++)fprintf(fpp,"%d\n",acc_surface.AC_normal[i]);
		fclose(fpp);
	}
}

//------ output protein feature in one file -------//2015_02_20//
// file format should be consistent with TPL file
/*
//////////// Features
  Num Res  Missing   SSE    CLE   ACC   pACC  CNa CNb   Xca       Yca       Zca       Xcb       Ycb       Zcb
   1   E      0       L      1     2     87   2   1     19.400     4.600    31.600    17.889     4.542    31.872
   2   N      0       E      5     2     48   4   3     20.500     7.600    33.700    21.326     8.630    32.887
   3   I      0       E      5     1     18   3   5     18.900     9.100    36.800    18.496     8.208    37.931
   4   E      0       E      5     2     53   3   2     19.900    12.700    37.700    19.630    13.781    36.734
   5   V      0       E      5     0      0   6   9     20.200    13.600    41.400    20.814    12.502    42.276
   6   H      0       E      5     1     32   6   3     20.700    17.200    42.800    19.590    18.263    42.565
   7   M      0       E      5     0      0   9   8     22.700    17.900    45.900    24.187    17.551    45.967
   8   L      0       E      5     1     12   5   4     21.100    20.900    47.700    19.558    20.886    47.424
   9   N      0       E      5     1     33   4   3     21.400    23.000    50.900    22.004    24.410    50.855
*/
//-------- ACC<->Int -------//
int ACC_To_Int(char c)
{
	switch(c)
	{
		case 'B': return 0;
		case 'M': return 1;
		case 'E': return 2;
		default: return 1;
	}
}
//----- output protein features -----//
void Output_Protein_Features(
	string &outroot,string &outname,Acc_Surface &acc_surface,XYZ *mol,XYZ *mcb,
	int moln,PDB_Residue *pdb,XYZ **mcc,int *mcc_side,char *ami,char *acc)
{
	//init
	int i,j;
	int retv;
	PDB_Chain_Fold chain_fold;
	chain_fold.initialize_simple(moln,' ');
	for(i=0;i<moln;i++)chain_fold.set_residue(i,pdb[i]);

	//------ calculate AMI,CLE,SSE ------//
	retv=chain_fold.calculate_CLE();
	if(retv!=0)
	{
		fprintf(stderr,"[%s]CLE_BAD!!!\n",outname.c_str());
		exit(-1);
	}
	retv=chain_fold.calculate_SSE();
	if(retv!=0)
	{
		fprintf(stderr,"[%s]SSE_BAD!!!\n",outname.c_str());
		exit(-1);
	}
	string AMI,CLE,SSE;
	AMI=chain_fold.get_sequence();
	CLE=chain_fold.get_CLE();
	SSE=chain_fold.get_SSE();

	//------ calculate ACC_Code and ACC_Value ------//
	for(i=0;i<moln;i++)pdb[i].get_XYZ_array(mcc[i],mcc_side[i]);
	acc_surface.AC_Calc_SolvAcc(mcc,ami,moln,acc,mcc_side);

	//------ calculate contact number for CA and CB ---//
	double DIST_CUTOFF = 64; //-> 8.0A (note that in older version of TPL, the CA/CB contact cutoff is 7.0A
	vector <int> cn_ca(moln,0);
	vector <int> cn_cb(moln,0);
	for(i=0;i<moln;i++)
	{
		for(j=0;j<moln;j++)
		{
			double dist_ca=mol[i].distance_square(mol[j]);
			double dist_cb=mcb[i].distance_square(mcb[j]);
			//-> normal condition
			if( abs(i-j)>=4 )
			{
				if(dist_ca <= DIST_CUTOFF)cn_ca[i]++;
				if(dist_cb <= DIST_CUTOFF)cn_cb[i]++;
			}
			//-> chain broken
			else if( abs(i-j)>=1 )
			{
				//check PDB_residue
				string str1,str2;
				pdb[i].get_PDB_residue_number(str1);
				pdb[j].get_PDB_residue_number(str2);
				str1=str1.substr(1,4);
				str2=str2.substr(1,4);
				int pos1=atoi(str1.c_str());
				int pos2=atoi(str2.c_str());
				if( abs(pos1-pos2)>=4 )
				{
					if(dist_ca <= DIST_CUTOFF)cn_ca[i]++;
					if(dist_cb <= DIST_CUTOFF)cn_cb[i]++;
				}
			}
		}
	}

/*
	//------ calculate core region ------//
	int MinHelixLen = 4;
	int MinBetaLen = 3;
	int MinContact = 1;
	vector <int> Core(moln,1);
	for(i=0;i<moln;i++)
	{
		//calculate core
		Core[i] = 1;
		if(SSE[i]!='H' && SSE[i]!='E') continue;
		Core[i] = 2;
		int sslen = 1;
		for(j=i-1;j>0;j--)
		{
			if(SSE[j]!=SSE[i]) break;
			sslen++;
		}
		for(j=i+1;j<moln;j++)
		{
			if(SSE[j]!=SSE[i]) break;
			sslen++;
		}
		if(SSE[i]=='H' && sslen >= MinHelixLen && cn_ca[i]>MinContact)Core[i] = 5;
		if(SSE[i]=='E' && sslen >= MinBetaLen  && cn_ca[i]>MinContact)Core[i] = 5;
	}
*/

	//------ output feature files -------//
	string file=outroot+"/"+outname+".feature";
	FILE *fpp=fopen(file.c_str(),"wb");
	if(fpp==0)
	{
		fprintf(stderr,"ERROR: file %s can't be opened. \n",file.c_str());
	}
	else
	{
		fprintf(fpp,"  Num Res  Missing   SSE    CLE   ACC   pACC  CNa CNb   Xca       Yca       Zca       Xcb       Ycb       Zcb\n");
		for(i=0;i<moln;i++)
		{
			char c=ACC_To_Int(acc[i]);
			fprintf(fpp,"%4d   %c      0       %c      %c     %1d    %3d  %2d  %2d   %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f\n",
				i+1,AMI[i],SSE[i],CLE[i],c,acc_surface.AC_normal[i],cn_ca[i],cn_cb[i],mol[i].X,mol[i].Y,mol[i].Z,mcb[i].X,mcb[i].Y,mcb[i].Z);
		}
	}
	fclose(fpp);
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
	XYZ **mbb;                  //BackBone (N,CA,C,O,CB)
	XYZ **mcc;                  //BackBone (N,CA,C,O,CB,...,)
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
		fprintf(stderr,"Process_List Not Found!![%s]\n",list.c_str());
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
			fprintf(stderr,"NO_INPUT_DIRECTORY[%s]!!!\n",buf.c_str());
			exit(-1);
		}
		path=buf;
		if(!getline(fin,buf,'\n'))
		{
			fprintf(stderr,"NO_OUTPUT_DIRECTORY[%s]!!!\n",buf.c_str());
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
			fprintf(stderr,"LIST_Format ERROR!!![%s]\n",buf.c_str());
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

			//get CB
			for(i=0;i<moln;i++)pdb[i].get_sidechain_atom( "CB ",mcb[i] );

			//output
			FILE *fpdb;
			file=outa+"/"+output+".pdb";
			fpdb=fopen(file.c_str(),"wb");
			if(fpdb==0)
			{
				fprintf(stderr,"ERROR: file %s can't be opened. \n",file.c_str());
			}
			else
			{
				mol_output.Output_PDB_III(fpdb,moln,pdb,c,OutType,OutMode,OutGlys);
				fclose(fpdb);
			}

			//output others
			if(OutFifi!=0)
			{
				if(OutFifi==1 || OutFifi==3 || OutFifi==5 || OutFifi==7)
					Output_Protein_Features_AMI_SSE_CLE(outa,output,moln,pdb);
				if(OutFifi==2 || OutFifi==3 || OutFifi==6 || OutFifi==7)
					Output_Protein_Features_ACC(outa,output,acc_surface,moln,pdb,mcc,mcc_side,ami,acc);
				if(OutFifi==4 || OutFifi==5 || OutFifi==6 || OutFifi==7)
					Output_Protein_Features(outa,output,acc_surface,mol,mcb,moln,pdb,mcc,mcc_side,ami,acc);
			}
		}

		//printf
		count++;
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
	Acc_Surface acc_surface;
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
	XYZ **mcc;                  //BackBone (N,CA,C,O,CB,...,)
	NewArray2D(&mbb,totlen,5);
	NewArray2D(&mcc,totlen,15);
	int *mcc_side=new int[totlen];
	char *ami=new char[totlen+1];
	char *cle=new char[totlen+1];
	char *acc=new char[totlen+1];

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

			//judge
			if(correct!=1)
			{
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

			//get CB
			for(i=0;i<moln;i++) pdb[i].get_sidechain_atom( "CB ",mcb[i] );

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
				fclose(fpdb);
			}

			//output others
			if(OutFifi!=0)
			{
				//-> get name and root
				string outroot,outname;
				getBaseName(output,outname,'/','.');
				getRootName(output,outroot,'/');
				//-> output files
				if(OutFifi==1 || OutFifi==3 || OutFifi==5 || OutFifi==7)
					Output_Protein_Features_AMI_SSE_CLE(outroot,outname,moln,pdb);
				if(OutFifi==2 || OutFifi==3 || OutFifi==6 || OutFifi==7)
					Output_Protein_Features_ACC(outroot,outname,acc_surface,moln,pdb,mcc,mcc_side,ami,acc);
				if(OutFifi==4 || OutFifi==5 || OutFifi==6 || OutFifi==7)
					Output_Protein_Features(outroot,outname,acc_surface,mol,mcb,moln,pdb,mcc,mcc_side,ami,acc);
			}
		}
	}
end:
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
}

//============== main ===============//
int main(int argc, char** argv)
{
	//---- BackBone ----//(process list or single) PDB_Backbone//
	{
		process_args(argc,argv);
		//-> if set with -F, then automatically set -R
		if(INPUT_FIFI!=0)INPUT_RECO=1;
		//-> list or single
		if(LIST_OR_SINGLE>0)WS_PDB_Back_Process(INPUT_LIST,INPUT_TYPE,INPUT_MODE,INPUT_GLYS,INPUT_NOCA,INPUT_RECO,INPUT_FIFI,INPUT_LOGF,WARN_OUT);
		else WS_PDB_Back_Process_Single(INPUT_NAM,INPUT_RAN,INPUT_OUT,INPUT_TYPE,INPUT_MODE,INPUT_GLYS,INPUT_NOCA,INPUT_RECO,INPUT_FIFI,INPUT_LOGF,WARN_OUT);
		exit(0);
	}
}
