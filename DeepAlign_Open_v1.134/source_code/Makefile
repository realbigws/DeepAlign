

# make the subdir first
CXX = g++
FLAG = -O3
SHELL :=/bin/bash

PDBToolDir = ./PDB_Tool_src
DeepScoreDir = ./DeepScore_src
DeepAlignDir = ./DeepAlign_src
DCOMBDIR = ./3DCOMB_src
PDBToolEXE = /PDB_Tool
DeepScoreEXE = /DeepScore
DeepAlignEXE = /DeepAlign
DCOMBEXE = /3DCOMB
UTILS = ./util_src/

all:
	cd $(PDBToolDir)&&make
	mv ${PDBToolDir}$(PDBToolEXE) ../
	cd $(DeepScoreDir)&&make
	mv $(DeepScoreDir)$(DeepScoreEXE) ../
	cd $(DeepAlignDir)&&make
	mv $(DeepAlignDir)$(DeepAlignEXE) ../
	cd $(DCOMBDIR)&&make
	mv $(DCOMBDIR)$(DCOMBEXE) ../
	
	$(CXX) DeepSearch.cpp $(FLAG) -lpthread -o DeepSearch
	mv DeepSearch ../
	
	$(CXX) $(UTILS)Stat_List.cpp $(FLAG) -o Stat_List
	$(CXX) $(UTILS)Fitting_EVD.cpp $(FLAG) -o Fitting_EVD
	$(CXX) $(UTILS)Domain_Proc.cpp $(FLAG) -o Domain_Proc
	$(CXX) $(UTILS)PDB_File_Cut.cpp $(FLAG) -o PDB_File_Cut
	mv Stat_List Fitting_EVD Domain_Proc PDB_File_Cut ../util/


PDB_Tool:
	cd $(PDBToolDir)&&make
	mv ${PDBToolDir}$(PDBToolEXE) ../

DeepScore:
	cd $(DeepScoreDir)&&make
	mv $(DeepScoreDir)$(DeepScoreEXE) ../

DeepAlign:
	cd $(DeepAlignDir)&&make
	mv $(DeepAlignDir)$(DeepAlignEXE) ../

3DCOMB:
	cd $(DCOMBDIR)&&make
	mv $(DCOMBDIR)$(DCOMBEXE) ../

DeepSearch:
	$(CXX) DeepSearch.cpp $(FLAG) -lpthread -o DeepSearch
	mv DeepSearch ../

util:
	$(CXX) $(UTILS)Stat_List.cpp $(FLAG) -o Stat_List
	$(CXX) $(UTILS)Fitting_EVD.cpp $(FLAG) -o Fitting_EVD
	$(CXX) $(UTILS)Domain_Proc.cpp $(FLAG) -o Domain_Proc
	$(CXX) $(UTILS)PDB_File_Cut.cpp $(FLAG) -o PDB_File_Cut
	mv Stat_List Fitting_EVD Domain_Proc PDB_File_Cut ../util/



# clean all the temp files
clean:
