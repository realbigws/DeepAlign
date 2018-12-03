#include <time.h>
#include <sstream>
#include <map>
#include "getopt.h"
#include "Fast_Sort.h"
#include "Mol_Load.h"
#include "Mol_Out.h"
#include "BLOSEN_CLE.h"
#include "MultiAlign_Out.h"


//-------- utility ------// -> already defined in Computation_Utility.h
/*
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
*/


//-> [1] for output PDB in linear co-ordinate, otherwise not output
int Out_PDB;          //-> output PDB


//--- macros -----//
int WS_OUTPUT_CA;
int WS_OUTPUT_TYPE;
int WS_REFINEMENT;
int HSFB_LEN_ORI;
double MSA_THRES;
//------------- WS_Neo_BLOMAPS -------------//__110105__//
Mol_Load mol_input;
Mol_Out mol_output;
BLOSEN_CLE blomaps;
MultiAlign_Out multi_out;
int *WS_LEN;
XYZ **WS_MOL;
XYZ **WS_MOL_TMP;
char **WS_AMI;
char **WS_NAM;
int **WS_ALI;
double **WS_ROT;
int WS_MAXLEN;
double *WS_RMSD_REC;  //-> for RMSD value output
int *WS_POS_REC;      //-> for RMSD value output
PDB_Residue **WS_ORIN;  //original input
PDB_Residue **WS_ORIN_TMP;
char **WS_IND;  //original index

//-------- init process -----//
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
void WS_Init_Struc(int maxnum,int maxlen)
{
	//init
	WS_LEN=new int[maxnum];
	NewArray2D(&WS_MOL,maxnum,maxlen);
	NewArray2D(&WS_MOL_TMP,maxnum,maxlen);
	NewArray2D(&WS_AMI,maxnum,maxlen);
	NewArray2D(&WS_NAM,maxnum,maxlen);
	NewArray2D(&WS_ALI,maxnum,maxnum*maxlen);
	NewArray2D(&WS_ROT,maxnum,12);
	WS_RMSD_REC=new double[maxnum*maxlen];
	WS_POS_REC=new int[maxnum*maxlen];
}
int WS_Return_Total(string &list)
{
	ifstream fin;
	string buf;
	//load
	fin.open(list.c_str(),ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"no such list![%s]\n",list.c_str());
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
int WS_Input_List(string &list,int &totnum,int limit_len)
{
	ifstream fin;
	string buf,temp;
	string name;
	string range;
	string file;
	//create
	char **ori_ali;
	NewArray2D(&ori_ali,100,30000);
	//load
	fin.open(list.c_str(),ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"no such list![%s]\n",list.c_str());
		return -1;
	}
	//read
	int i;
	int moln;
	int retv;
	int ret_val;
	int maxlen=0;
	if(WS_OUTPUT_TYPE>0)
	{
		WS_ORIN=new PDB_Residue*[totnum]; //original PDB_Residue
		WS_ORIN_TMP=new PDB_Residue*[totnum];
		WS_IND=new char*[totnum]; //original index
	}
	int cur=0;
	for(i=0;i<totnum;i++)
	{
		if(!getline(fin,buf,'\n'))break;
		retv=String_Process_Line_Str(buf,' ',ori_ali);
		if(retv<2)
		{
			file=ori_ali[0];
			range="_";
		}
		else
		{
			file=ori_ali[0];
			range=ori_ali[1];
		}
		//input CA
		ret_val=Get_PDB_File_Len(file);
		if(ret_val<=0 || ret_val>limit_len)
		{
			fprintf(stderr,"%s -> pdb bad or length overange ! %d \n",file.c_str(),ret_val);
			continue;
		}
		ret_val=mol_input.XYZ_Input(file,range,0,WS_LEN[cur],0,0,0,0);
		if(ret_val<1)
		{
			//---- error specific report ---//
			switch(ret_val)
			{
				case -999:
				{
					fprintf(stderr,"%s -> file not found!\n",buf.c_str());
					break;
				}
				case -333:
				{
					fprintf(stderr,"%s -> overrange problem has been found!\n",buf.c_str());
					break;
				}
				case -888:
				{
					fprintf(stderr,"%s -> the PDB file contains format errors!\n",buf.c_str());
					break;
				}
				case -100:
				{
					fprintf(stderr,"%s -> the resdiue range format is error!\n",buf.c_str());
					break;
				}
				case -987:
				{
					fprintf(stderr,"%s -> the starting position should less than the ending position!\n",buf.c_str());
					break;
				}
				case -432:
				{
					fprintf(stderr,"%s -> the chain identification not found!\n",buf.c_str());
					break;
				}
				case -234:
				{
					fprintf(stderr,"%s -> the residue numbers are out of range!\n",buf.c_str());
					break;
				}
				default:
				{
					fprintf(stderr,"%s -> FILE LOAD ERROR!\n",buf.c_str());
					break;
				}
			}
			exit(-1);
		}
		else if(ret_val==2)
		{
			fprintf(stderr,"%s -> the ending position is over range!!\n",buf.c_str());
			exit(-1);
		}
		if(WS_LEN[cur]<HSFB_LEN_ORI)
		{
			fprintf(stderr,"%s -> the length %d of this protein is less than minimal %d!\n",
				buf.c_str(),WS_LEN[cur],HSFB_LEN_ORI);
			continue;
		}

		//record maximal
		if(WS_LEN[cur]>maxlen)maxlen=WS_LEN[cur];
		strcpy(WS_NAM[cur],buf.c_str());
		//input original
		if(WS_OUTPUT_TYPE>0)
		{
			WS_ORIN[cur]=new PDB_Residue[WS_LEN[cur]];
			WS_ORIN_TMP[cur]=new PDB_Residue[WS_LEN[cur]];
			WS_IND[cur]=new char[WS_LEN[cur]*6];
			int PRE_LOAD_=mol_input.PRE_LOAD;
			int WARNING_out_=mol_input.WARNING_out;
			mol_input.PRE_LOAD=1;
			mol_input.WARNING_out=0;
			mol_input.XYZ_Input(file,range,0,moln,WS_MOL[cur],WS_AMI[cur],0,WS_IND[cur],WS_ORIN[cur]);
			mol_input.PRE_LOAD=PRE_LOAD_;
			mol_input.WARNING_out=WARNING_out_;
		}
		cur++;
		//printf
//		printf("load[%3d]\r",cur);
	}
	totnum=cur;
	return maxlen;
}
void WS_BLOMAPS_Init(int maxnum,int maxlen)
{
	blomaps.Multi_Ali_Create(maxnum,maxlen);
	blomaps.BLOSEN_Create(maxnum,maxlen);
	blomaps.BLOSEN_CLE_Create(maxnum,maxlen);
	blomaps.MultiAlign_Cent_Create(maxnum,maxlen);
	multi_out.BC_Output_Init(maxnum,maxlen);
}
void Transfer_Alignment(int totnum,int totlen,int **ali,vector < vector < int > > &out)
{
	int j,k;
	int tag;
	int count;
	for(j=0;j<totnum;j++)
	{
		count=1;
		for(k=0;k<totlen;k++)
		{
			tag=ali[j][k];
			if(tag==0) //gap
			{
				out[j][k]=-1*count;
			}
			else //capital
			{
				out[j][k]=count;
				count++;
			}
		}
	}
}
void Output_Alignment(FILE *fp,char **in,char **nam,int totnum,int totlen,int **ali)
{
	int j,k;
	int tag;
	int count;
	char c;
	for(j=0;j<totnum;j++)
	{
		count=0;
		fprintf(fp,">%s\n",nam[j]);
		for(k=0;k<totlen;k++)
		{
			tag=ali[j][k];
			if(tag==0) //gap
			{
				fprintf(fp,"-");
			}
			else //capital
			{
				c=in[j][count];
				if(c=='Z')fprintf(fp,".");
				else fprintf(fp,"%c",c);
				count++;
			}
		}
		fprintf(fp,"\n");
	}
}
void Output_Rotmat(FILE *fp,char **nam,int totnum,double **rotmat)
{
	int j;
	for(j=0;j<totnum;j++)
	{
		fprintf(fp,">%s\n",nam[j]);
		//output
		fprintf(fp,"%9.6f ",rotmat[j][0]);
		fprintf(fp,"%9.6f ",rotmat[j][1]);
		fprintf(fp,"%9.6f ",rotmat[j][2]);
		fprintf(fp,"%12.6f\n",rotmat[j][9]);
		fprintf(fp,"%9.6f ",rotmat[j][3]);
		fprintf(fp,"%9.6f ",rotmat[j][4]);
		fprintf(fp,"%9.6f ",rotmat[j][5]);
		fprintf(fp,"%12.6f\n",rotmat[j][10]);
		fprintf(fp,"%9.6f ",rotmat[j][6]);
		fprintf(fp,"%9.6f ",rotmat[j][7]);
		fprintf(fp,"%9.6f ",rotmat[j][8]);
		fprintf(fp,"%12.6f\n",rotmat[j][11]);
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

//-------- chain mapping --------//
int Dom_To_Int(char pos)
{
	if(pos>='A' && pos<='Z')return pos-'A';
	if(pos>='a' && pos<='z')return pos-'a'+26;
	if(pos>='1' && pos<='9')return pos-'1'+52;
	return 61;
}
char Int_To_Dom(int pos)
{
	if(pos>=0 && pos<=25)return pos+'A';
	if(pos>=26 && pos<=51)return pos-26+'a';
	if(pos>=52 && pos<=60)return pos-52+'1';
	return '0';
}

//------- Int Char Mapping ---------//
int Char_To_Int(char pos)
{
	if(pos>='0' && pos<='9')return pos-'0';
	if(pos>='A' && pos<='Z')return pos-'A'+10;
	if(pos>='a' && pos<='z')return pos-'a'+36;
	return 0;
}
char Int_To_Char(int pos)
{
	if(pos>=0 && pos<=9)return pos+'0';
	if(pos>=10 && pos<=35)return pos-10+'A';
	if(pos>=36 && pos<=61)return pos-36+'a';
	return '0';
}




//============= main_process =================//
void WS_NEO_BLOMAPS(string &list,int totnum,string &outnam,string &outroot)
{
	Kabsch kabsch;
	//parameter
	int HSFB_LEN=12;
	int PIVOT_MOL_CUT=totnum;
	int ANCHOR_HSFB_CUT=1;
	int COLOR1_HSFB_CUT=5;
	int COLOR2_HSFB_CUT=5;
	int ZoomIn_TopK=5;
	int CONSEN_STAGE=0;
	blomaps.BLOSEN_CLE_Parameter(HSFB_LEN,PIVOT_MOL_CUT,ANCHOR_HSFB_CUT,COLOR1_HSFB_CUT,COLOR2_HSFB_CUT,ZoomIn_TopK,CONSEN_STAGE);
	blomaps.TM_Score_Type=0;  //-> optimize TMscore instead of DeepScore
//	blomaps.PRITNForNOT=1;  //printf temporary information

	//----- the following parameter is for Large_Fragment_Refinement -------//__110201__//
	if(WS_REFINEMENT==1)
	{
		blomaps.CONSEN_STAGE=3; //DO consensus optimization
		blomaps.TMorLEN=0;      //using TMscore as optimization [0]
		blomaps.RECorNOT=0;     //DON'T record initial [0]
		blomaps.ROTorNOT=0;     //DON'T do final consensus [0]
		blomaps.MERGE_DIST=3.0; //column merge only [3.0]
	}

	//process
//	clock_t Time_Start=clock();
	double sumtms;
	sumtms=blomaps.BLOMAPS_Partial(totnum,WS_LEN,WS_MOL,WS_AMI,WS_ALI,WS_ROT);
//	clock_t Time_End=clock();
	//sort result
	int ret_totnum=(int)blomaps.FIN_SCO_REC.size();
	double *temp_score=new double[ret_totnum];
	int *temp_index=new int[ret_totnum];
	{
		for(int k=0;k<ret_totnum;k++)temp_score[k]=(blomaps.FIN_SCO_REC[k]+0.01)*(blomaps.FIN_CORE_LEN[k]+1);
		Fast_Sort <double> fast_sort;
		fast_sort.fast_sort_1(temp_score,temp_index,ret_totnum);
	}
	//---- MSA similarity check for multi solutions ----//
	vector <int> indec_rec;
	indec_rec.push_back(0);
	double best_sco=temp_score[temp_index[0]];
	if(WS_OUTPUT_TYPE>=2)
	{
		//get position
		vector <vector <vector <int> > > fin_ali_rec;
		for(int k=0;k<ret_totnum;k++)
		{
			WS_MAXLEN=blomaps.FIN_BEST_LEN[k];
			for(int ii=0;ii<totnum;ii++)for(int jj=0;jj<WS_MAXLEN;jj++)WS_ALI[ii][jj]=blomaps.FIN_ALI_REC[k][ii][jj];
			vector <vector <int > > fin_ali_tmp (totnum, vector <int >(WS_MAXLEN));
			Transfer_Alignment(totnum,WS_MAXLEN,WS_ALI,fin_ali_tmp);
			fin_ali_rec.push_back(fin_ali_tmp);
		}
		//check overlap
		for(int k=1;k<ret_totnum;k++)
		{
			//init check
			if(temp_score[temp_index[k]]<best_sco*0.8)break;    //-> we use 0.8 cutoff here !! 
			if(blomaps.FIN_SCO_REC[temp_index[k]]<0.35)break;   //-> we show multi solutions for TMscore > 0.35 !!
			//calculate
			int correct=1;
			for(int l=0;l<k;l++)
			{
				//omit failed one
				if(temp_score[temp_index[l]]<0)continue;
				//calculate overlap_Sco
				double overlap_sco=blomaps.Compare_Two_MSA(fin_ali_rec[temp_index[k]],fin_ali_rec[temp_index[l]],
					blomaps.FIN_BEST_LEN[temp_index[k]],blomaps.FIN_BEST_LEN[temp_index[l]],totnum);
				if( 1.0*overlap_sco/40 > MSA_THRES)
				{
					correct=0;
					break;
				}
			}
			//check
			if(correct==1)indec_rec.push_back(k);
		}
	}
	//multi solution
	int rel_size=(int)indec_rec.size();
	int rel_num=rel_size<10?rel_size:10;
	for(int k=rel_num-1;k>=0;k--)
	{
		if(WS_OUTPUT_TYPE<2)if(k!=0)continue;
		int index=temp_index[indec_rec[k]];
		//superimpose CA structures
		WS_MAXLEN=blomaps.FIN_BEST_LEN[index];
		for(int ii=0;ii<totnum;ii++)for(int jj=0;jj<12;jj++)WS_ROT[ii][jj]=blomaps.FIN_ROT_REC[index][ii][jj];
		for(int ii=0;ii<totnum;ii++)for(int jj=0;jj<WS_MAXLEN;jj++)WS_ALI[ii][jj]=blomaps.FIN_ALI_REC[index][ii][jj];		
		for(int i=0;i<totnum;i++)blomaps.rot_mol(WS_MOL[i],WS_MOL_TMP[i],WS_LEN[i],WS_ROT[i]);
		//output multi solution
		if(WS_OUTPUT_TYPE>0)
		{
			FILE *fp;
			char www_nam[300000];
			if(WS_OUTPUT_TYPE>=2)sprintf(www_nam,"%s/%s.%d",outroot.c_str(),outnam.c_str(),k);
			else sprintf(www_nam,"%s/%s",outroot.c_str(),outnam.c_str());
			string rel_outnam=www_nam;
			//[0] superimpose original structure
			for(int i=0;i<totnum;i++)WS_Superimpose_FullAtom(kabsch,WS_ORIN[i],WS_LEN[i],WS_ORIN_TMP[i],WS_ROT[i]); //original PDB
			//[1] output script
			string f2,f3;
			if(WS_OUTPUT_CA>=1)
			{
				f3=rel_outnam+".scp";
				multi_out.BC_Output_All(f3,WS_MOL_TMP,WS_AMI,WS_LEN,totnum,WS_MAXLEN,WS_ALI,WS_OUTPUT_CA);
			}
			//[2] output alignment
			string f1=rel_outnam+".ali";
			fp=fopen(f1.c_str(),"wb");
			if(fp==0)
			{
				fprintf(stderr,"ERROR: file %s can't be opened. \n",f1.c_str());
			}
			else
			{
				Output_Alignment(fp,WS_AMI,WS_NAM,totnum,WS_MAXLEN,WS_ALI);
				fclose(fp);
			}
			//[3] output structure 
			string TER="TER                                                                             ";
			string END="END                                                                             ";
			f2=rel_outnam+".pdb";
			fp=fopen(f2.c_str(),"wb");
			if(fp==0)
			{
				fprintf(stderr,"ERROR: file %s can't be opened. \n",f2.c_str());
			}
			else
			{
				for(int i=0;i<totnum;i++)
				{
					char chain=Int_To_Dom(i);
					mol_output.Output_PDB_III(fp,WS_LEN[i],WS_ORIN_TMP[i],chain,1);
					fprintf(fp,"%s\n",TER.c_str());
				}
				fprintf(fp,"%s\n",END.c_str());
				//--> output for linear chain
				if(Out_PDB==1)
				{
					// output file
					string f2_=rel_outnam+".pdb_linear";
					FILE *fq=fopen(f2_.c_str(),"wb");
					if(fq==0)
					{
						fprintf(stderr,"ERROR: file %s can't be opened. \n",f2_.c_str());
					}
					// output
					for(int i=0;i<totnum;i++)
					{
						char chain=Int_To_Dom(i);
						mol_output.Output_PDB_III(fq,WS_LEN[i],WS_ORIN_TMP[i],chain,-1);
						fprintf(fq,"%s\n",TER.c_str());
					}
					fprintf(fq,"%s\n",END.c_str());
					fclose(fq);
				}
			}
			fclose(fp);
			//[4] output rotmat
			string f4=rel_outnam+".rmt";
			fp=fopen(f4.c_str(),"wb");
			if(fp==0)
			{
				fprintf(stderr,"ERROR: file %s can't be opened. \n",f4.c_str());
			}
			else
			{
				Output_Rotmat(fp,WS_NAM,totnum,WS_ROT);
				fclose(fp);
			}
			//[5] output RMSD val
			blomaps.BC_Return_Conserved_Core_Pos(totnum,WS_MAXLEN,WS_ALI,WS_POS_REC);
			blomaps.BC_Return_Conserved_RMSD_Pos(totnum,WS_MAXLEN,WS_ALI,WS_MOL_TMP,WS_RMSD_REC);
			string f5=rel_outnam+".rms";
			fp=fopen(f5.c_str(),"wb");
			if(fp==0)
			{
				fprintf(stderr,"ERROR: file %s can't be opened. \n",f5.c_str());
			}
			else
			{
				fprintf(fp,">column_conservation\n");
				for(int i=0;i<WS_MAXLEN;i++)
				{
					char numc;
					if(WS_POS_REC[i]<totnum)numc=Int_To_Char(WS_POS_REC[i]);
					else numc='|';
					fprintf(fp,"%c",numc);
				}
				fprintf(fp,"\n");
				fprintf(fp,">column_RMSD\n");
				for(int i=0;i<WS_MAXLEN;i++)
				{
					double rmsd=WS_RMSD_REC[i];
					int rmsdi;
					if(rmsd<0)rmsdi=0;
					else if(rmsd>9)rmsdi=9;
					else rmsdi=(int)rmsd;
					fprintf(fp,"%d",rmsdi);
				}
				fprintf(fp,"\n");
				fclose(fp);
			}
			//[6] output score
			double totsco=blomaps.BC_Calc_SumTM_Given_Mol(WS_MOL_TMP,WS_LEN,totnum,WS_MAXLEN,WS_ALI);
			int corelen=blomaps.BC_Return_Conserved_Core(totnum,WS_MAXLEN,WS_ALI);
			double totrms=blomaps.BC_Return_Conserved_RMSD(totnum,WS_MAXLEN,WS_ALI,WS_MOL_TMP);
			string f6=rel_outnam+".sco";
			fp=fopen(f6.c_str(),"wb");
			if(fp==0)
			{
				fprintf(stderr,"ERROR: file %s can't be opened. \n",f6.c_str());
			}
			else
			{
				fprintf(fp,"CORE_LEN = %4d \n",corelen);
				fprintf(fp,"RMSD = %7.3f \n",totrms);
				fprintf(fp,"Ave_TMscore = %6.4f \n",totsco);
				fclose(fp);
			}
		}
		//only printout the best
		if(k==0) 
		{
			double totsco=blomaps.BC_Calc_SumTM_Given_Mol(WS_MOL_TMP,WS_LEN,totnum,WS_MAXLEN,WS_ALI);
			int corelen=blomaps.BC_Return_Conserved_Core(totnum,WS_MAXLEN,WS_ALI);
			double totrms=blomaps.BC_Return_Conserved_RMSD(totnum,WS_MAXLEN,WS_ALI,WS_MOL_TMP);
			printf("input:%s, #structures:%d\n",list.c_str(),totnum);
			printf("Results:\n");
			printf("Objective function value is %lf\n",(corelen+1)*(totsco+0.01));
			printf("CORE_LEN=%4d, RMSD=%7.3f, Ave_TMscore=%6.4f\n",corelen,totrms,totsco);
		}
	}
	//time consume
//	double time_consume=1.0*(Time_End-Time_Start)/CLOCKS_PER_SEC;
//	printf("%lf\n",time_consume);
}

//------- usage ------//
void Usage(void)
{
/*
	fprintf(stderr,"3DCOMB V1.18 [Dec-30-2013] \n");
	fprintf(stderr,"----------------------- REFERENCE ---------------------------------\n");
	fprintf(stderr,"Sheng Wang, Jian Peng and Jinbo Xu.\n");
	fprintf(stderr,"Alignment of distantly-related protein structures: algorithm, \n");
	fprintf(stderr,"                 bound and implications to homology modeling. \n");
	fprintf(stderr,"                         Bioinformatics. 2011 Sep 15;27(18):2537-45 \n");
	fprintf(stderr,"======================= USAGE: =====================================\n");  
	fprintf(stderr,"./3DCOMB <-i input_list> [-o out_name] [-p out_option] [-s] [-r] \n");
//	fprintf(stderr,"3DCOMB <-i input_list> [-o out_name] [-d out_root]  \n"); 
//	fprintf(stderr,"       [-p out_option] [-m msa_cutoff] [-s] [-r] \n");
	fprintf(stderr,"----------------------- OPTIONS: -----------------------------------\n");
	fprintf(stderr,"<input_list>: is the input file listing the structures to be aligned.\n");
	fprintf(stderr,"   Each line in <input_list> specifies one input structure.\n");
	fprintf(stderr,"[out_name]: is the output name. Default is the basename of input_list.\n");
//	fprintf(stderr,"[out_root]: is the output directory. Default is the current directory.\n");
	fprintf(stderr,"[out_option]: 0, don't output detail information; \n");
	fprintf(stderr,"   [1], single solution; 2, multiple solutions. \n");
//	fprintf(stderr,"[msa_cutoff]: MSA similarity measure between 0 to 1 for multi-solution.\n");
//	fprintf(stderr,"   Solutions above the given value will be eliminated. [Default=0.75] \n");
	fprintf(stderr,"[-s]: is an optional argument for JMOL script output.\n");
	fprintf(stderr,"[-r]: is an optional argument for more iteration refinement.\n");
	fprintf(stderr,"--------------------------------------------------------------------\n");
*/

	fprintf(stderr,"3DCOMB v1.18 \n\n");
	fprintf(stderr,"Sheng Wang, Jian Peng and Jinbo Xu.\n");
	fprintf(stderr,"ALIGNMENT OF DISTANTLY-RELATED PROTEIN STRUCTURES: ALGORITHM, \n");
	fprintf(stderr,"                    BOUND AND IMPLICATIONS TO HOMOLOGY MODELING. \n");
	fprintf(stderr,"                         Bioinformatics. 2011 Sep 15;27(18):2537-45 \n");
	fprintf(stderr,"Usage: \n");
	fprintf(stderr,"./3DCOMB -i input_list [-o out_name] [-p out_option] [-s script_option] [-r] \n\n");
	fprintf(stderr,"Options:\n\n");
	fprintf(stderr,"-i input_list:          The input file listing the structures to be aligned. \n");
	fprintf(stderr,"                        Each line in <input_list> specifies one input structure.\n\n");
	fprintf(stderr,"-o out_name:            The output name for the alignment files, if specified. \n");
	fprintf(stderr,"                        Otherwise, the output name will be the basename of <input_list>. \n\n");
	fprintf(stderr,"-p out_option:          0,  do not output the alignment files. \n");
	fprintf(stderr,"                       [1], output alignment files for single solution. (Set as default) \n");
	fprintf(stderr,"                        2,  output alignment files for multiple solutions. \n\n");
	fprintf(stderr,"-s script_option:      [0], do not output script files. (Set as default) \n");
	fprintf(stderr,"                        1,  output script for JMol. \n");
	fprintf(stderr,"                        2,  output script for RasMol.   \n\n");
	fprintf(stderr,"-r:                     An optional argument for more iteration refinement. \n\n");

}

//---------- main -----------//
int main(int argc,char **argv)
{
	Out_PDB=0;     //-> default, not output linear PDB

	//---- macros -----//
	WS_OUTPUT_CA=0;   //no script output
	WS_REFINEMENT=0;  //no refinement
	WS_OUTPUT_TYPE=1; //single solution
	HSFB_LEN_ORI=12;  //HSFB length
	MSA_THRES=0.75;   //MSA comparison for multi solution
	mol_input.CbBACK=1;
	
	//-------- 3DCOMB -----------//__101110__//
	{
		if(argc<3)
		{
			Usage();
			exit(-1);
		}
		//basic input
		string input_name="";
		string output_name="";
		string output_root=".";

		//---- process argument ----//
		extern char* optarg;
		int c=0;
		while((c=getopt(argc,argv,"i:o:O:d:p:m:s:r"))!=EOF)
		{
			switch(c) 
			{
				//---- pairsie job ---//
				case 'i':
					input_name = optarg;
					break;
				case 'o':
					output_name = optarg;
					break;
				case 'O':
					//-> [1] for output PDB in linear co-ordinate, otherwise not output
					Out_PDB = atoi(optarg);
					break;
				case 'd':
					output_root = optarg;
					break;
				case 'p':
					WS_OUTPUT_TYPE = atoi(optarg);
					break;
				case 'm':
					MSA_THRES = atof(optarg);
					break;
				case 's':
					WS_OUTPUT_CA = atoi(optarg);
					break;
				case 'r':
					WS_REFINEMENT = 1;
					break;
				//----- default ----//
				default:
					Usage();
					exit(-1);
			}
		}

		//check
		if(output_name=="")getBaseName(input_name,output_name,'/','.');
		if(WS_OUTPUT_TYPE<0 || WS_OUTPUT_TYPE>3)
		{
			fprintf(stderr,"out_option %d should be 0,1 or 2 \n",WS_OUTPUT_TYPE);
			exit(-1);
		}
		if(WS_OUTPUT_CA<0 || WS_OUTPUT_CA>2)
		{
			fprintf(stderr,"script_option %d should be 0,1 or 2 \n",WS_OUTPUT_CA);
			exit(-1);
		}
		if(MSA_THRES<0 || MSA_THRES>1)
		{
			fprintf(stderr,"msa_cutoff %lf should between 0 to 1 \n",MSA_THRES);
			exit(-1);
		}
		int input_num=WS_Return_Total(input_name);
		if(input_num<2)
		{
			fprintf(stderr,"3DCOMB input file failed: too few number %d\n",input_num);
			exit(-1);
		}
		WS_Init_Struc(input_num,PROT_MAX_NUM*2);
		int max_len=WS_Input_List(input_name,input_num,PROT_MAX_NUM*2);
		if(max_len<0)
		{
			fprintf(stderr,"3DCOMB input file failed!\n");
			exit(-1);
		}
		if(input_num<2)
		{
			fprintf(stderr,"3DCOMB input number %d failed!\n",input_num);
			exit(-1);
		}

		//process
		WS_BLOMAPS_Init(input_num,max_len);
		WS_NEO_BLOMAPS(input_name,input_num,output_name,output_root);
		exit(0);
	}
}
