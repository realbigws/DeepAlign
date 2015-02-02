CXX = g++
CCFLAGS = -O3 -Wall -ffast-math -static 
OBJ = Utility.o getopt.o Computation_Utility.o DeepScore.o Mol_Load.o Mol_Ori.o Mol_Out.o PDB_Utility.o Backbone_Sidechain.o PDB_Residue.o PDB_Chain.o PDB_File.o XYZ.o Kabsch.o TM_score.o TM_align.o Confo_Beta.o Confo_Lett.o Dynamic_Programming.o PDB_To_CLE.o

make: $(OBJ)
	${CXX} $(CCFLAGS)  -o DeepScore $(OBJ)
	rm -f $(OBJ)

Utility.o  : Utility.h
	${CXX} ${CCFLAGS} -c Utility.cpp
getopt.o : getopt.h
	${CXX} ${CCFLAGS} -c getopt.cpp
Computation_Utility.o  : Computation_Utility.h
	${CXX} ${CCFLAGS} -c Computation_Utility.cpp
DeepScore.o  : 
	${CXX} ${CCFLAGS} -c DeepScore.cpp
Confo_Beta.o  : Confo_Beta.h
	${CXX} ${CCFLAGS} -c Confo_Beta.cpp
Confo_Lett.o  : Confo_Lett.h
	${CXX} ${CCFLAGS} -c Confo_Lett.cpp
PDB_To_CLE.o  : PDB_To_CLE.h
	${CXX} ${CCFLAGS} -c PDB_To_CLE.cpp
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
XYZ.o  : XYZ.h
	${CXX} ${CCFLAGS} -c XYZ.cpp
Kabsch.o  : Kabsch.h
	${CXX} ${CCFLAGS} -c Kabsch.cpp
TM_score.o  : TM_score.h
	${CXX} ${CCFLAGS} -c TM_score.cpp
TM_align.o  : TM_align.h
	${CXX} ${CCFLAGS} -c TM_align.cpp
Dynamic_Programming.o  : Dynamic_Programming.h
	${CXX} ${CCFLAGS} -c Dynamic_Programming.cpp

clean:
	rm -f DeepScore $(OBJ)
