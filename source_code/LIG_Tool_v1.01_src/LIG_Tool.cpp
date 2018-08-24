#include <map>
#include "PDB_Chain_Fold.h"
#include "Confo_Back.h"
#include "Confo_Beta.h"
#include "Confo_Lett.h"
#include "Acc_Surface.h"
#include "Mol_File.h"
#include "Mol_Out.h"
#include "getopt.h"
#include "Ligand_Utility.h"
using namespace std;


//-----------------------------------------------------------------------------------------------------------//
//---- print_help_msg => print help message
void print_help_msg(void) 
{
	cout << "========================================================|" << endl;
	cout << "LIG_Tool  (version 1.01) [2015.02.20]                   |" << endl;
	cout << "    Extract ligands and chains from official PDB file   |" << endl;
	cout << "Usage:   ./LIG_Tool <-i input_pdb> [-o out_name]        |" << endl;
	cout << "         [-p chain_out_root] [-q ligand_out_root]       |" << endl;
	cout << "         [-n length_cut] [-d distance_cut] [-l log]     |" << endl;
	cout << "--------------------------------------------------------|" << endl;
	cout << "-i input_pdb : input original PDB file, e.g., 1col.pdb  |" << endl;
	cout << "-o out_name  : output file name. [default: input_name]  |" << endl;
	cout << "-p chain_out_root : chain output root [default: ./ ]    |" << endl;
	cout << "-q ligand_out_root: ligand output root [default: ./ ]   |" << endl;
	cout << "-n length_cut : chain length cutoff [default: 25 ]      |" << endl;
	cout << "-d distance_cut : ligand distance cutoff [default: 4.0] |" << endl;
	cout << "-l log : output log files (1) or not (0) [default: 0]   |" << endl;
	cout << "========================================================|" << endl;
	exit(-1);
}
//---- WebServer's default input ----//
string INPUT_FILE="";
string OUTPUT_NAME="";
string CHAIN_OUTROOT="./";  //chain output root
string LIGAND_OUTROOT="./"; //ligand output root
int LENGTH_CUT=25;          //chain length cutoff
double DISTANCE_CUT=4.0;    //ligand distance cutoff
int LOG_OR_NOT=0; //default: don't output log

//-----------------------------------------------------------------------------------------------------------//
//---- parameter editor ----//
static option long_options[] =
{
	{"input",   required_argument, NULL, 'i'},
	{"output",  no_argument,       NULL, 'o'},
	{"chain",   no_argument,       NULL, 'p'},
	{"ligand",  no_argument,       NULL, 'q'},
	{"length",  no_argument,       NULL, 'n'},
	{"dist",    no_argument,       NULL, 'd'},
	{"log",     no_argument,       NULL, 'l'},
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
		opt=getopt_long(argc,argv,"i:o:p:q:n:d:l:",
			   long_options,&option_index);
		if (opt==-1)break;	
		switch(opt) 
		{
			case 'i':
				INPUT_FILE=optarg;
				break;
			case 'o':
				OUTPUT_NAME=optarg;
				break;
			case 'p':
				CHAIN_OUTROOT=optarg;
				break;
			case 'q':
				LIGAND_OUTROOT=optarg;
				break;
			case 'n':
				LENGTH_CUT=atoi(optarg);
				break;
			case 'd':
				DISTANCE_CUT=atof(optarg);
				break;
			case 'l':
				LOG_OR_NOT=atoi(optarg);
				break;
			default:
				exit(-1);
		}
	}
}

//====== data structure ======//
struct Ligand_Struc
{
	string lig_name;
	int lig_moln;
	vector <XYZ> lig_xyz;
	vector <string> lig_data;
};

//--------- Ligand_Mapping --------//
int WS_Ligand_Trans(char code) 
{ 
	switch(code) 
	{ 
		case 'A': return 0; 
		case 'B': return 1; 
		case 'C': return 2; 
		case 'D': return 3; 
		case 'E': return 4; 
		case 'F': return 5; 
		case 'G': return 6; 
		case 'H': return 7; 
		case 'I': return 8; 
		case 'J': return 9; 
		case 'K': return 10; 
		case 'L': return 11; 
		case 'M': return 12; 
		case 'N': return 13; 
		case 'O': return 14; 
		case 'P': return 15; 
		case 'Q': return 16; 
		case 'R': return 17; 
		case 'S': return 18; 
		case 'T': return 19; 
		case 'U': return 20; 
		case 'V': return 21; 
		case 'W': return 22; 
		case 'X': return 23; 
		case 'Y': return 24; 
		case 'Z': return 25; 
		case '0': return 26; 
		case '1': return 27; 
		case '2': return 28; 
		case '3': return 29; 
		case '4': return 30; 
		case '5': return 31; 
		case '6': return 32; 
		case '7': return 33; 
		case '8': return 34; 
		case '9': return 35;
		case ' ': return 36;
		default:return -1; 
	} 
}

//================== All_Ligands_Process =============//__110530__//
//-> ligand_string_to_xyz
int Ligand_String_To_XYZ(vector <string> &input, vector <XYZ> &output)
{
	int i;
	int size=(int)input.size();
	output.clear();
	string buf,temp;
	for(i=0;i<size;i++)
	{
		buf=input[i];
		temp=buf.substr(30,8);
		XYZ xyz;
		xyz.X=atof(temp.c_str());
		temp=buf.substr(38,8);
		xyz.Y=atof(temp.c_str());
		temp=buf.substr(46,8);
		xyz.Z=atof(temp.c_str());
		output.push_back(xyz);
	}
	return size;
}
//--------- extract ligand from PDB ------------//
int PDB_Extract_Ligand(string &pdb,vector <Ligand_Struc> &output)
{
	//---- get MODRES ----//
	Mol_File mol_input;
	map <string,string > modres_map;
	{
		const string sss=pdb;
		mol_input.Process_MODRES_Mapping(sss,modres_map);
	}
	//---- get MODRES ----//over

	output.clear();
	//--- list for mapping ---//
	map<string, int > ws_mapping;
	map<string, int>::iterator iter;
	ws_mapping.clear();
	//--- data for mapping ---//
	vector <int> ws_mapp_int;
	vector <string> ws_mapp_nam;
	vector <vector < vector <string> > > ws_mapp_string;
	ws_mapp_int.clear();
	ws_mapp_nam.clear();
	ws_mapp_string.clear();
	string lig_dummy;

	//--- list for mapping ---//
	ifstream fin;
	string buf,temp;
	string ws_nam;
	getBaseName(pdb,ws_nam,'/','.');
	//read
	fin.open(pdb.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"list %s not found!!\n",pdb.c_str());
		return -1;
	}
	int i;
	int ii;
	int len;
	int first=1;
	int tot_lig=0;
	int key;
	char ws3[4];
	char chain;
	vector <string> ws_record;
	ws_record.clear();
	int prev=-1;
	string prev_chain;
	string current_chain;
	char ws_prev[4];
	//processs
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		len=(int)buf.length();
		if(len<4)continue;
		temp=buf.substr(0,4);
		if(temp=="END " || temp=="ENDM")break; //this might be modified in the future
		if(temp=="TER ") //output "ter"
		{
			//output
			if(first==0)
			{
				first=1;
				//record
				vector <string> ws_record_tmp;
				for(i=0;i<(int)ws_record.size();i++)
				{
					if(ws_record[i][77]!='H')
					{
						ws_record_tmp.push_back(ws_record[i]);
					}
				}
				//output
				lig_dummy="";
				lig_dummy=lig_dummy+ws_prev;
				iter = ws_mapping.find(lig_dummy);
				if(iter != ws_mapping.end()) //already exist ligand
				{
					key=ws_mapping[lig_dummy];
					ws_mapp_int[key-1]++;
					ws_mapp_string[key-1].push_back(ws_record_tmp);
				}
				else                         //newly appear ligand
				{
					tot_lig++;
					ws_mapping.insert(map < string, int >::value_type(lig_dummy, tot_lig));
					string wsss=lig_dummy;
					for(ii=0;ii<(int)wsss.length();ii++)if(wsss[ii]==' ')wsss[ii]='_';
					ws_mapp_nam.push_back(wsss);
					ws_mapp_int.push_back(1);
					vector < vector <string> > dummy_vec;
					ws_mapp_string.push_back(dummy_vec);
					ws_mapp_string[tot_lig-1].push_back(ws_record_tmp);
				}
			}
		}
		//get atom
		if(temp!="ATOM" && temp!="HETA")continue;
		if(len<20)continue;
		temp=buf.substr(17,3);
		if(temp[0]==' ')//transfer
		{
			for(ii=0;ii<3;ii++)ws3[ii]=' ';
			ws3[ii]='\0';
			if(temp[1]==' ')
			{
				ws3[0]=temp[2];
			}
			else
			{
				ws3[0]=temp[1];
				ws3[1]=temp[2];
			}
		}
		else strcpy(ws3,temp.c_str());

		//MODRES MAP
		{
			string temp_ori=ws3;
			string temp_mod=temp_ori;
			mol_input.MODRES_Map(temp_ori,temp_mod,modres_map);
			for(ii=0;ii<3;ii++)ws3[ii]=temp_mod[ii];
		}
		//MODRES MAP over

		//record
		chain=buf[21];
		current_chain=buf.substr(21,6);
		int result=0;
		for(i=0;i<3;i++)result+=(WS_Ligand_Trans(ws3[i]))*(int)pow(37.0,1.0*i);
		int retv=WS_Ligand_Num_Code(result);
		if(retv!=-1)
		{
			if(first==1)
			{
				first=0;
				prev=retv;
				prev_chain=current_chain;
				strcpy(ws_prev,ws3);
				ws_record.clear();
				ws_record.push_back(buf);
			}
			else
			{
				if(retv==prev && current_chain==prev_chain)ws_record.push_back(buf);
				else
				{
					//record
					vector <string> ws_record_tmp;
					for(i=0;i<(int)ws_record.size();i++)
					{
						if(ws_record[i][77]!='H')
						{
							ws_record_tmp.push_back(ws_record[i]);
						}
					}
					//output
					lig_dummy="";
					lig_dummy=lig_dummy+ws_prev;
					iter = ws_mapping.find(lig_dummy);
					if(iter != ws_mapping.end()) //already exist ligand
					{
						key=ws_mapping[lig_dummy];
						ws_mapp_int[key-1]++;
						ws_mapp_string[key-1].push_back(ws_record_tmp);
					}
					else                         //newly appear ligand
					{
						tot_lig++;
						ws_mapping.insert(map < string, int >::value_type(lig_dummy, tot_lig));
						string wsss=lig_dummy;
						for(ii=0;ii<(int)wsss.length();ii++)if(wsss[ii]==' ')wsss[ii]='_';
						ws_mapp_nam.push_back(wsss);
						ws_mapp_int.push_back(1);
						vector < vector <string> > dummy_vec;
						ws_mapp_string.push_back(dummy_vec);
						ws_mapp_string[tot_lig-1].push_back(ws_record_tmp);
					}
					//continue
					prev=retv;
					prev_chain=current_chain;
					strcpy(ws_prev,ws3);
					ws_record.clear();
					ws_record.push_back(buf);
				}
			}
		}
		else
		{
			if(first==0)
			{
				first=1;
				//record
				vector <string> ws_record_tmp;
				for(i=0;i<(int)ws_record.size();i++)
				{
					if(ws_record[i][77]!='H')
					{
						ws_record_tmp.push_back(ws_record[i]);
					}
				}
				//output
				lig_dummy="";
				lig_dummy=lig_dummy+ws_prev;
				iter = ws_mapping.find(lig_dummy);
				if(iter != ws_mapping.end()) //already exist ligand
				{
					key=ws_mapping[lig_dummy];
					ws_mapp_int[key-1]++;
					ws_mapp_string[key-1].push_back(ws_record_tmp);
				}
				else                         //newly appear ligand
				{
					tot_lig++;
					ws_mapping.insert(map < string, int >::value_type(lig_dummy, tot_lig));
					string wsss=lig_dummy;
					for(ii=0;ii<(int)wsss.length();ii++)if(wsss[ii]==' ')wsss[ii]='_';
					ws_mapp_nam.push_back(wsss);
					ws_mapp_int.push_back(1);
					vector < vector <string> > dummy_vec;
					ws_mapp_string.push_back(dummy_vec);
					ws_mapp_string[tot_lig-1].push_back(ws_record_tmp);
				}
			}
		}
	}
	//final
	if(first==0)
	{
		first=1;
		//record
		vector <string> ws_record_tmp;
		for(i=0;i<(int)ws_record.size();i++)
		{
			if(ws_record[i][77]!='H')
			{
				ws_record_tmp.push_back(ws_record[i]);
			}
		}
		//output
		lig_dummy="";
		lig_dummy=lig_dummy+ws_prev;
		iter = ws_mapping.find(lig_dummy);
		if(iter != ws_mapping.end()) //already exist ligand
		{
			key=ws_mapping[lig_dummy];
			ws_mapp_int[key-1]++;
			ws_mapp_string[key-1].push_back(ws_record_tmp);
		}
		else                         //newly appear ligand
		{
			tot_lig++;
			ws_mapping.insert(map < string, int >::value_type(lig_dummy, tot_lig));
			string wsss=lig_dummy;
			for(ii=0;ii<(int)wsss.length();ii++)if(wsss[ii]==' ')wsss[ii]='_';
			ws_mapp_nam.push_back(wsss);
			ws_mapp_int.push_back(1);
			vector < vector <string> > dummy_vec;
			ws_mapp_string.push_back(dummy_vec);
			ws_mapp_string[tot_lig-1].push_back(ws_record_tmp);
		}
	}

	//-> final collect
	int j,k;
	char command[5];
	int tot_num=0;
	for(i=0;i<tot_lig;i++)
	{
		for(j=0;j<ws_mapp_int[i];j++)
		{
			//get name
			sprintf(command,"%-4d",tot_num);
			command[4]='\0';
			for(k=0;k<(int)strlen(command);k++)if(command[k]==' ')command[k]='_';
			string name_tmp=command;
			string name=ws_mapp_nam[i]+"_"+name_tmp;
			//assign
			vector <XYZ> xyz_tmp;
			int moln=Ligand_String_To_XYZ(ws_mapp_string[i][j],xyz_tmp);
			Ligand_Struc ligand_struc;
			ligand_struc.lig_name=name;
			ligand_struc.lig_moln=moln;
			ligand_struc.lig_xyz=xyz_tmp;
			ligand_struc.lig_data=ws_mapp_string[i][j];
			//incremental
			output.push_back(ligand_struc);
			tot_num++;
		}
	}
	//-> final return
	return tot_num;
}


//================== All_Chain_Process =============//__110530__//
//-> Compare ligand and chain 
int Compare_Ligand_and_Chain(vector <XYZ> &protein, vector <XYZ> &ligand,double r_cut)
{
	int i,j;
	int ll=(int)ligand.size();
	int pl=(int)protein.size();
	//distance check
	double dist2;
	double thres2=r_cut*r_cut;
	for(i=0;i<ll;i++)
	{
		for(j=0;j<pl;j++)
		{
			dist2=ligand[i].distance_square(protein[j]);
			if(dist2<thres2)return 1; //success !! (within radius)
		}
	}
	//final return
	return 0; //failed !! (not within radius)
}


//-> Compare ligand and chain complex
double Residue_Ligand_Distance(PDB_Residue & PDB, vector <XYZ> &ligand)
{
	int i,k;
	int number;
	int totnum;
	XYZ xyz;
	double dist2;
	double minval=99999;
	number=PDB.get_sidechain_totnum();
	totnum=(int)ligand.size();
	for(k=0;k<number;k++)
	{
		//get_sidechain
		if(PDB.get_sidechain_part_index(k)==0)continue;
		PDB.get_sidechain_atom(k,xyz);
		//calculate distance
		for(i=0;i<totnum;i++)
		{
			dist2=ligand[i].distance_square(xyz);
			if(dist2<minval)minval=dist2;
		}
	}
	//return
	return minval;
}
void Kill_Space(string &in,string &out)
{
	int i;
	out="";
	for(i=0;i<(int)in.length();i++)
	{
		if(in[i]!=' ')out+=in[i];
	}
}
int Compare_Ligand_and_Chain_Complex(vector <PDB_Residue> &protein, vector <XYZ> &ligand,double r_cut, 
	vector <int> &pos_rec, vector <char> &cha_rec, vector <string> &ind_rec, vector <double> &min_rec)
{
	int i;
	int pl=(int)protein.size();
	//distance check
	double dist2;
	double thres2=r_cut*r_cut;
	pos_rec.clear();
	cha_rec.clear();
	ind_rec.clear();
	min_rec.clear();
	int count=0;
	for(i=0;i<pl;i++)
	{
		dist2=Residue_Ligand_Distance(protein[i], ligand);
		if(dist2<thres2)
		{
			pos_rec.push_back(i);
			char c=protein[i].get_AA();
			string ind_;
			protein[i].get_PDB_residue_number(ind_);
			string ind;
			ind=ind_.substr(1,5);
			Kill_Space(ind,ind_);
			cha_rec.push_back(c);
			ind_rec.push_back(ind_);
			min_rec.push_back(sqrt(dist2));
			count++;
		}
	}
	//return
	return count;
}


//-> reconstruct missing heavy-atom
int Reconstruct_Missing(vector <PDB_Residue> &pdb)
{
	//init check
	int i,j;
	int moln=(int)pdb.size();
	for(i=0;i<moln;i++)if(pdb[i].get_backbone_part_index(1)==0)return -1; //no CA !!
	//data
	Confo_Beta confo_beta(moln);
	Confo_Back confo_back(moln);
	Confo_Lett confo_lett;
	XYZ *mol=new XYZ[moln];      //CA
	char *ami=new char[moln+1];  //ami
	//assign
	char amino;
	for(i=0;i<moln;i++)
	{
		pdb[i].get_backbone_atom(1,mol[i]);
		amino=pdb[i].get_AA();
		ami[i]=amino;
	}
	ami[moln]='\0';
	//check
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
	if(correct==1)
	{
		delete [] mol;
		delete [] ami;
		return 0; //no modification
	}
	//reconstruct
	//[0]data
	char *cle=new char[moln+1];
	XYZ *mcb=new XYZ[moln];     //CB
	XYZ **mbb;                  //BackBone (N,CA,C,O,CB)
	NewArray2D(&mbb,moln,5);
	//[1]recon
	confo_lett.btb_ori(0,0,0,moln,mol,cle);
	cle[moln]='\0';
	confo_back.Recon_Back_WS_Main(mol,cle,moln,mbb);      //given CA, recon BackBone (N,CA,C,O,CB)
	confo_beta.WS_Recon_Beta_21(mol,mcb,moln,ami,cle);    //given CA, recon CB
	//[2]assign
	for(i=0;i<moln;i++)
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
	//final return
	delete [] mol;
	delete [] ami;
	delete [] mcb;
	delete [] cle;
	DeleteArray2D(&mbb,moln);
	return 1; //reconstruct
}

//-> PDB_Residue_To_XYZ
int PDB_Residue_To_XYZ(PDB_Residue &PDB, vector <XYZ> &output)
{
	int k;
	int number;
	XYZ xyz;
	int count=0;
	//backbone_out
	number=PDB.get_backbone_totnum();  //this value should be 4
	for(k=0;k<number;k++)
	{
		//get_backbone
		if(PDB.get_backbone_part_index(k)==0)continue;
		PDB.get_backbone_atom(k,xyz);
		output.push_back(xyz);
		count++;
	}
	char amino=PDB.get_AA();
	if(amino=='G')return count;
	//sidechain_out
	number=PDB.get_sidechain_totnum();
	for(k=0;k<number;k++)
	{
		//get_sidechain
		if(PDB.get_sidechain_part_index(k)==0)continue;
		PDB.get_sidechain_atom(k,xyz);
		output.push_back(xyz);
		count++;
	}
	//return
	return count;
}
int PDB_Residue_To_XYZ_Total(vector <PDB_Residue> &pdb, vector <XYZ> &output)
{
	int i;
	int count=0;
	int retv=0;
	int moln=(int)pdb.size();
	output.clear();
	for(i=0;i<moln;i++)
	{
		retv=PDB_Residue_To_XYZ(pdb[i],output);
		count+=retv;
	}
	return count;
}

//-> Extract all chains for ligands
//[note]: len_cut -> for pdb_chain
//        r_cut   -> for ligand
int PDB_Ligand_All_Process(string &file,string &out_name,
	string &pdb_out_dir,string &ligand_out_dir,int len_cut,double r_cut)
{
	//class
	Mol_File mol_input;
	Mol_Out mol_output;
	mol_input.MODRES=1;
	//data
	vector <PDB_Residue> pdb;
	vector <XYZ> mol;
	vector <XYZ> total_xyz;
	string ami;
	vector <string> ind;
	PDB_Chain pdb_chain;
	PDB_Residue residue;
	int retv;
	//chain load
	vector <PDB_Chain> chains;
	retv=mol_input.PDB_read_pdb_file(file,chains);
	if(retv<0)return 0; //failed
	//ligand load	
	vector <Ligand_Struc> ligands;
	retv=PDB_Extract_Ligand(file,ligands);
	if(retv<0)return 0; //failed
	//chain-ligand process
	int i,j,k;
	int moln;
	string chain;
	string cur_nam;
	string resnum;
	int chain_size;
	chain_size=(int)chains.size();
	int lig_size;
	lig_size=(int)ligands.size();
	//fprintf data
	FILE *fp;
	int has_ligand;
	string outname;
	//log
	FILE *f1,*f2,*f3,*f4,*f5;
	string log_file;
	f1=0;
	f2=0;
	f3=0;
	f4=0;
	f5=0;
	//real process
	int ligand_log_first=1;
	int chain_log_first=1;
	for(i=0;i<chain_size;i++)
	{
		//read
		pdb_chain=chains[i];
		moln=pdb_chain.get_length();
		chain=pdb_chain.get_chain_id();
		//length check
		if(moln<len_cut)continue;
		//assign
		pdb.resize(moln);
		mol.resize(moln);
		ami.resize(moln);
		ind.resize(moln);
		for(k=0;k<moln;k++)
		{
			pdb_chain.get_residue(k,residue);
			pdb[k]=residue;
			residue.get_backbone_atom(1,mol[k]);
			ami[k]=residue.get_AA();
			residue.get_PDB_residue_number(resnum);
			resnum=resnum.substr(1,5);
			ind[k]=resnum;
		}
		Reconstruct_Missing(pdb);
		PDB_Residue_To_XYZ_Total(pdb,total_xyz); 
		//check length and ligand
		has_ligand=0;
		for(j=0;j<lig_size;j++)
		{
			retv=Compare_Ligand_and_Chain(total_xyz,ligands[j].lig_xyz,r_cut);
			if(retv==0)continue;
			//get detailed binding site
			vector <int> pos_rec;
			vector <char> cha_rec;
			vector <string> ind_rec;
			vector <double> min_rec;
			int count=Compare_Ligand_and_Chain_Complex(pdb,ligands[j].lig_xyz,r_cut,pos_rec,cha_rec,ind_rec,min_rec);
			//output ligand
			cur_nam=out_name+chain+"_"+ligands[j].lig_name;
			outname=ligand_out_dir+"/"+cur_nam+".pdb";
			fp=fopen(outname.c_str(),"wb");
			for(k=0;k<(int)ligands[j].lig_data.size();k++)fprintf(fp,"%s\n",ligands[j].lig_data[k].c_str());
			fclose(fp);
			//output ligand_log
			if(LOG_OR_NOT==1)
			{
				if(ligand_log_first==1)
				{
					ligand_log_first=0;
					log_file=out_name+".ligand_log";
					f2=fopen(log_file.c_str(),"wb");
					fprintf(f2,">%s\n",out_name.c_str());
					log_file=out_name+".ligand_chain";
					f3=fopen(log_file.c_str(),"wb");
				}
				fprintf(f2,"%s -> ",cur_nam.c_str());
				for(k=0;k<count;k++)fprintf(f2,"%d|%s|%c|%3.1f ",pos_rec[k],ind_rec[k].c_str(),cha_rec[k],min_rec[k]);
				fprintf(f2,"\n");
			}
			has_ligand++;
		}
		if(LOG_OR_NOT==1)
		{
			if(has_ligand>0)fprintf(f3,"%s%s\n",out_name.c_str(),chain.c_str());
		}
		//output pdb
		cur_nam=out_name+chain;
		outname=pdb_out_dir+"/"+cur_nam+".pdb";
		fp=fopen(outname.c_str(),"wb");
		mol_output.Output_PDB_III(fp,moln,pdb,chain[0]);
		fclose(fp);
		//output pdb_log
		if(LOG_OR_NOT==1)
		{
			if(chain_log_first==1)
			{
				chain_log_first=0;
				log_file=out_name+".chain_log";
				f1=fopen(log_file.c_str(),"wb");
				log_file=out_name+".atom_seq";
				f4=fopen(log_file.c_str(),"wb");
				log_file=out_name+".atom_ind";
				f5=fopen(log_file.c_str(),"wb");
			}
			fprintf(f1,"%s %4d\n",cur_nam.c_str(),moln);
			fprintf(f4,">%s\n",cur_nam.c_str());
			fprintf(f4,"%s\n",ami.c_str());
			fprintf(f5,">%s\n",cur_nam.c_str());
			for(k=0;k<moln;k++)fprintf(f5,"%s",ind[k].c_str());
			fprintf(f5,"\n");
		}
	}
	//final
	if(LOG_OR_NOT==1)
	{
		if(ligand_log_first==0)
		{
			fclose(f2);
			fclose(f3);
		}
		if(chain_log_first==0)
		{
			fclose(f1);
			fclose(f4);
			fclose(f5);
		}
	}
	return 1; //success
}

//============== main ===============//
int main(int argc, char** argv)
{
	//---- PDB_Ligand_Process ---//process all chains and ligands
	{
		process_args(argc,argv);
		//judge
		if(INPUT_FILE=="")
		{
			fprintf(stderr,"-i INPUT_FILE should not be blank \n");
			exit(-1);
		}
		if(OUTPUT_NAME=="")
		{
			getBaseName(INPUT_FILE,OUTPUT_NAME,'/','.');
		}
		if(LENGTH_CUT<0)
		{
			fprintf(stderr,"-n LENGTH_CUT %d should be positive \n",LENGTH_CUT);
			exit(-1);
		}
		if(DISTANCE_CUT<0)
		{
			fprintf(stderr,"-d DISTANCE_CUT %lf should be positive \n",DISTANCE_CUT);
			exit(-1);
		}
		if(LOG_OR_NOT<0 || LOG_OR_NOT>1)
		{
			fprintf(stderr,"-l LOG_OR_NOT %d should be 0 or 1 \n",LOG_OR_NOT);
			exit(-1);
		}
		//assign
		string input_pdb=INPUT_FILE;
		string out_name=OUTPUT_NAME;
		string pdb_out_root=CHAIN_OUTROOT;
		string ligand_out_root=LIGAND_OUTROOT;
		//create output
		char command[30000];
		sprintf(command,"mkdir -p %s",pdb_out_root.c_str());
		system(command);
		sprintf(command,"mkdir -p %s",ligand_out_root.c_str());
		system(command);
		//process
		int len_cut=LENGTH_CUT;
		double distance_cut=DISTANCE_CUT;
		PDB_Ligand_All_Process(input_pdb,out_name,pdb_out_root,ligand_out_root,len_cut,distance_cut);
		exit(0);
	}
}
