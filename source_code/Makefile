

# make the subdir first
CXX = g++
#FLAG = -O3 -static
FLAG = -O3
SHELL :=/bin/bash

LIGToolDir = ./LIG_Tool_src
LIGToolEXE = /LIG_Tool
UTILS = ./util_src/

all: LIG_Tool util

LIG_Tool:
	cd $(LIGToolDir)&&make
	mv ${LIGToolDir}$(LIGToolEXE) ../

util:
	$(CXX) $(UTILS)Chemical_Component_Dictionary.cpp $(FLAG) -o Chemical_Component_Dictionary
	mkdir -p ../util
	mv Chemical_Component_Dictionary ../util/

# clean all the temp files
clean:
