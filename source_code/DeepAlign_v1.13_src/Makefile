CXX = g++
CCFLAGS = -O3 -Wall -ffast-math 
OBJ = Utility.o Computation_Utility.o Confo_Lett.o getopt.o DeepAlign.o Mol_Load.o Mol_Ori.o Mol_Out.o PDB_Utility.o Backbone_Sidechain.o PDB_Residue.o PDB_Chain.o PDB_File.o PDB_To_CLE.o XYZ.o Ali_Ori.o Ali_AFP.o Ali_Ali3.o Kabsch.o Bioinfo_Code.o CLEPAPS_Ori.o CLEPAPS_Out.o CLEFAPS.o CLEFAPS_CLE.o CLEFAPS_Main.o Dynamic_Programming.o TM_score.o TM_align.o TM_Align_Main.o Envo_Align.o Confo_Beta.o

make: $(OBJ)
	${CXX} $(CCFLAGS)  -o DeepAlign $(OBJ)
	rm -f $(OBJ)

Utility.o  : Utility.h
	${CXX} ${CCFLAGS} -c Utility.cpp
Computation_Utility.o  : Computation_Utility.h
	${CXX} ${CCFLAGS} -c Computation_Utility.cpp
Confo_Lett.o  : Confo_Lett.h
	${CXX} ${CCFLAGS} -c Confo_Lett.cpp
getopt.o  : getopt.h
	${CXX} ${CCFLAGS} -c getopt.cpp
DeepAlign.o  : reference_nam.h reference_ami.h reference_cle.h reference_mca.h reference_mcb.h
	${CXX} ${CCFLAGS} -c DeepAlign.cpp
Mol_Load.o  : Mol_Load.h
	${CXX} ${CCFLAGS} -c Mol_Load.cpp
Mol_Ori.o  : Mol_Ori.h
	${CXX} ${CCFLAGS} -c Mol_Ori.cpp
Mol_Out.o  : Mol_Out.h
	${CXX} ${CCFLAGS} -c Mol_Out.cpp
PDB_Utility.o  : PDB_Utility.h
	${CXX} ${CCFLAGS} -c PDB_Utility.cpp
Backbone_Sidechain.o  : Backbone_Sidechain.h
	${CXX} ${CCFLAGS} -c Backbone_Sidechain.cpp
PDB_Residue.o  : PDB_Residue.h
	${CXX} ${CCFLAGS} -c PDB_Residue.cpp
PDB_Chain.o  : PDB_Chain.h
	${CXX} ${CCFLAGS} -c PDB_Chain.cpp
PDB_File.o  : PDB_File.h
	${CXX} ${CCFLAGS} -c PDB_File.cpp
PDB_To_CLE.o  : PDB_To_CLE.h
	${CXX} ${CCFLAGS} -c PDB_To_CLE.cpp
XYZ.o  : XYZ.h
	${CXX} ${CCFLAGS} -c XYZ.cpp
Ali_Ori.o  : Ali_Ori.h
	${CXX} ${CCFLAGS} -c Ali_Ori.cpp
Ali_AFP.o  : Ali_AFP.h Fast_Sort.h
	${CXX} ${CCFLAGS} -c Ali_AFP.cpp
Ali_Ali3.o  : Ali_Ali3.h
	${CXX} ${CCFLAGS} -c Ali_Ali3.cpp	
Kabsch.o  : Kabsch.h
	${CXX} ${CCFLAGS} -c Kabsch.cpp
Bioinfo_Code.o  : Bioinfo_Code.h Align_Utility.h
	${CXX} ${CCFLAGS} -c Bioinfo_Code.cpp
CLEPAPS_Ori.o  : CLEPAPS_Ori.h Align_Utility.h
	${CXX} ${CCFLAGS} -c CLEPAPS_Ori.cpp
CLEPAPS_Out.o  : CLEPAPS_Out.h
	${CXX} ${CCFLAGS} -c CLEPAPS_Out.cpp
CLEFAPS.o  : CLEFAPS.h
	${CXX} ${CCFLAGS} -c CLEFAPS.cpp
CLEFAPS_CLE.o  : CLEFAPS_CLE.h
	${CXX} ${CCFLAGS} -c CLEFAPS_CLE.cpp
CLEFAPS_Main.o  : CLEFAPS_Main.h
	${CXX} ${CCFLAGS} -c CLEFAPS_Main.cpp	
Dynamic_Programming.o  : Dynamic_Programming.h
	${CXX} ${CCFLAGS} -c Dynamic_Programming.cpp
TM_score.o  : TM_score.h
	${CXX} ${CCFLAGS} -c TM_score.cpp
TM_align.o  : TM_align.h
	${CXX} ${CCFLAGS} -c TM_align.cpp
TM_Align_Main.o  : TM_Align_Main.h
	${CXX} ${CCFLAGS} -c TM_Align_Main.cpp
Envo_Align.o  : Envo_Align.h
	${CXX} ${CCFLAGS} -c Envo_Align.cpp
Confo_Beta.o  : Confo_Beta.h
	${CXX} ${CCFLAGS} -c Confo_Beta.cpp

clean:
	rm -f DeepAlign $(OBJ)
