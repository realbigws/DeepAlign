

# make the subdir first
CXX = g++
#FLAG = -O3 -static 
FLAG = -O3
SHELL :=/bin/bash

DeepScoreDir = ./DeepScore_src
DeepAlignDir = ./DeepAlign_src
3DCOMBDIR = ./3DCOMB_src
UTILS = ./util_src/

DeepScoreEXE = /DeepScore
DeepAlignEXE = /DeepAlign
3DCOMBEXE = /3DCOMB


all: DeepScore DeepAlign 3DCOMB PDB_Rotmol util


DeepScore:
	cd $(DeepScoreDir)&&make
	mv $(DeepScoreDir)$(DeepScoreEXE) ../

DeepAlign:
	cd $(DeepAlignDir)&&make
	mv $(DeepAlignDir)$(DeepAlignEXE) ../

3DCOMB:
	cd $(3DCOMBDIR)&&make
	mv $(3DCOMBDIR)$(3DCOMBEXE) ../

PDB_Rotmol:
	$(CXX) PDB_Rotmol_src/*.cpp $(FLAG) -o PDB_Rotmol
	mv PDB_Rotmol ../util

util:
	$(CXX) $(UTILS)Stat_List.cpp $(FLAG) -o Stat_List
	$(CXX) $(UTILS)Fitting_EVD.cpp $(FLAG) -o Fitting_EVD
	$(CXX) $(UTILS)Domain_Proc.cpp $(FLAG) -o Domain_Proc
	$(CXX) $(UTILS)PDB_File_Cut.cpp $(FLAG) -o PDB_File_Cut
	$(CXX) $(UTILS)PDB_To_SEQ.cpp $(FLAG) -o PDB_To_SEQ
	$(CXX) $(UTILS)TopK_by_EVD.cpp $(FLAG) -o TopK_by_EVD
	$(CXX) $(UTILS)List_Div_Shuf.cpp $(FLAG) -o List_Div_Shuf
	$(CXX) $(UTILS)DeepSearch_Rank.cpp $(FLAG) -o DeepSearch_Rank
	$(CXX) $(UTILS)PDB_Replace_XYZ.cpp $(FLAG) -o PDB_Replace_XYZ
	$(CXX) $(UTILS)PDB_Atom_Start.cpp $(FLAG) -o PDB_Atom_Start
	$(CXX) $(UTILS)PDB_Resi_Start.cpp $(FLAG) -o PDB_Resi_Start
	$(CXX) $(UTILS)PDB_Resi_Trans.cpp $(FLAG) -o PDB_Resi_Trans
	$(CXX) $(UTILS)PDB_Resi_Copy.cpp $(FLAG) -o PDB_Resi_Copy
	$(CXX) $(UTILS)PDB_Add_Chain.cpp $(FLAG) -o PDB_Add_Chain
	$(CXX) $(UTILS)Alignment_Mapping.cpp $(FLAG) -o Alignment_Mapping
	mkdir -p ../util
	mv Stat_List Fitting_EVD Domain_Proc PDB_File_Cut PDB_To_SEQ TopK_by_EVD List_Div_Shuf DeepSearch_Rank ../util
	mv PDB_Replace_XYZ PDB_Atom_Start PDB_Resi_Start PDB_Resi_Trans PDB_Resi_Copy PDB_Add_Chain Alignment_Mapping ../util/

# clean all the temp files
clean:
