

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
	$(CXX) $(UTILS)TopK_by_EVD.cpp $(FLAG) -o TopK_by_EVD
	$(CXX) $(UTILS)List_Div_Shuf.cpp $(FLAG) -o List_Div_Shuf
	$(CXX) $(UTILS)DeepSearch_Rank.cpp $(FLAG) -o DeepSearch_Rank
	$(CXX) $(UTILS)PDB_Replace_XYZ.cpp $(FLAG) -o PDB_Replace_XYZ
	mkdir -p ../util
	mv Stat_List Fitting_EVD Domain_Proc PDB_File_Cut TopK_by_EVD List_Div_Shuf DeepSearch_Rank PDB_Replace_XYZ ../util/

# clean all the temp files
clean:
