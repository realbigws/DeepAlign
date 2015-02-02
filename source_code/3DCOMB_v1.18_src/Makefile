CXX = g++
CCFLAGS = -O3 -Wall -ffast-math -static
OBJ = Ali_AFP.o Ali_Ali3.o Ali_Ori.o Bioinfo_Code.o BLOMAPS_Ori.o BLOSEN.o BLOSEN_CLE.o Computation_Utility.o Confo_Lett.o Kabsch.o getopt.o 3DCOMB.o Mol_Load.o Mol_Ori.o Mol_Out.o MultiAlign_Cent.o MultiAlign_Out.o Utility.o PDB_Utility.o Backbone_Sidechain.o PDB_Residue.o PDB_Chain.o PDB_File.o PDB_To_CLE.o TM_align.o TM_score.o XYZ.o ZoomIn_Align.o Dynamic_Programming.o

make: $(OBJ)
	${CXX} $(CCFLAGS)  -o 3DCOMB $(OBJ)
	rm -f $(OBJ)

Ali_AFP.o  : Ali_AFP.h Fast_Sort.h
	${CXX} ${CCFLAGS} -c Ali_AFP.cpp
Ali_Ali3.o  : Ali_Ali3.h
	${CXX} ${CCFLAGS} -c Ali_Ali3.cpp
Ali_Ori.o  : Ali_Ori.h
	${CXX} ${CCFLAGS} -c Ali_Ori.cpp
Bioinfo_Code.o  : Bioinfo_Code.h
	${CXX} ${CCFLAGS} -c Bioinfo_Code.cpp
BLOMAPS_Ori.o  : BLOMAPS_Ori.h
	${CXX} ${CCFLAGS} -c BLOMAPS_Ori.cpp
BLOSEN.o  : BLOSEN.h
	${CXX} ${CCFLAGS} -c BLOSEN.cpp
BLOSEN_CLE.o  : BLOSEN_CLE.h
	${CXX} ${CCFLAGS} -c BLOSEN_CLE.cpp
Computation_Utility.o  : Computation_Utility.h
	${CXX} ${CCFLAGS} -c Computation_Utility.cpp
Confo_Lett.o  : Confo_Lett.h
	${CXX} ${CCFLAGS} -c Confo_Lett.cpp
Kabsch.o  : Kabsch.h
	${CXX} ${CCFLAGS} -c Kabsch.cpp
getopt.o  : getopt.h
	${CXX} ${CCFLAGS} -c getopt.cpp
3DCOMB.o  : Fast_Sort.h
	${CXX} ${CCFLAGS} -c 3DCOMB.cpp
Mol_Load.o  : Mol_Load.h
	${CXX} ${CCFLAGS} -c Mol_Load.cpp
Mol_Ori.o  : Mol_Ori.h
	${CXX} ${CCFLAGS} -c Mol_Ori.cpp
Mol_Out.o  : Mol_Out.h
	${CXX} ${CCFLAGS} -c Mol_Out.cpp
MultiAlign_Cent.o  : MultiAlign_Cent.h
	${CXX} ${CCFLAGS} -c MultiAlign_Cent.cpp
MultiAlign_Out.o  : MultiAlign_Out.h
	${CXX} ${CCFLAGS} -c MultiAlign_Out.cpp
Utility.o  : Utility.h
	${CXX} ${CCFLAGS} -c Utility.cpp
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
TM_align.o  : TM_align.h
	${CXX} ${CCFLAGS} -c TM_align.cpp
TM_score.o  : TM_score.h
	${CXX} ${CCFLAGS} -c TM_score.cpp
XYZ.o  : XYZ.h
	${CXX} ${CCFLAGS} -c XYZ.cpp
ZoomIn_Align.o  : ZoomIn_Align.h
	${CXX} ${CCFLAGS} -c ZoomIn_Align.cpp
Dynamic_Programming.o  : Dynamic_Programming.h
	${CXX} ${CCFLAGS} -c Dynamic_Programming.cpp

clean:
	rm -f 3DCOMB $(OBJ)
