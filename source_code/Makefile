

# make the subdir first
CXX = g++
#FLAG = -O3 -static 
FLAG = -O3
SHELL :=/bin/bash

PDBToolDir = ./PDB_Tool_src
LIGToolDir = ./LIG_Tool_src
DeepScoreDir = ./DeepScore_src
DeepAlignDir = ./DeepAlign_src
3DCOMBDIR = ./3DCOMB_src

PDBToolEXE = /PDB_Tool
LIGToolEXE = /LIG_Tool
DeepScoreEXE = /DeepScore
DeepAlignEXE = /DeepAlign
3DCOMBEXE = /3DCOMB
UTILS = ./util_src/


all: PDB_Tool DeepScore DeepAlign 3DCOMB DeepSearch util


PDB_Tool:
	cd $(PDBToolDir)&&make
	mv ${PDBToolDir}$(PDBToolEXE) ../

LIG_Tool:
	cd $(LIGToolDir)&&make
	mv ${LIGToolDir}$(LIGToolEXE) ../

DeepScore:
	cd $(DeepScoreDir)&&make
	mv $(DeepScoreDir)$(DeepScoreEXE) ../

DeepAlign:
	cd $(DeepAlignDir)&&make
	mv $(DeepAlignDir)$(DeepAlignEXE) ../

3DCOMB:
	cd $(3DCOMBDIR)&&make
	mv $(3DCOMBDIR)$(3DCOMBEXE) ../

DeepSearch:
	$(CXX) DeepSearch.cpp getopt.cpp $(FLAG) -lpthread -o DeepSearch
	mv DeepSearch ../

util:
	$(CXX) $(UTILS)Stat_List.cpp $(FLAG) -o Stat_List
	$(CXX) $(UTILS)Fitting_EVD.cpp $(FLAG) -o Fitting_EVD
	$(CXX) $(UTILS)Domain_Proc.cpp $(FLAG) -o Domain_Proc
	$(CXX) $(UTILS)PDB_File_Cut.cpp $(FLAG) -o PDB_File_Cut
	$(CXX) $(UTILS)distribute_ii_openmp.cpp $(FLAG) -D_OMP -fopenmp -o distribute_ii_openmp
	mkdir -p ../util
	mv Stat_List Fitting_EVD Domain_Proc PDB_File_Cut distribute_ii_openmp ../util/



# clean all the temp files
clean:
