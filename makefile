#Source directory
SDIR := source

#Source file extension
SEXT := cxx

#Object directory
ODIR := object

$(shell mkdir -p $(SDIR); mkdir -p $(ODIR);) 

#Source files
SFILES = $(wildcard $(SDIR)/*.$(SEXT))

#Object files
OFILES = $(patsubst $(SDIR)/%.$(SEXT),$(ODIR)/%.o,$(SFILES))

#Executable file name 
EXE = obddfact

#C++ compiler
CXX = g++

#C++ COMPILER_FLAGS 
CXXFLAGS = 

#Include Paths
IPATHS = -Iheader \
	 -Ie:/Libs/buddy22/include

#Library Paths
LPATHS = -Le:/Libs/buddy22/lib

#Linker Flags
LFLAGS = -lbdd -lgmpxx -lgmp -lpsapi

.PHONY: all
all: $(OFILES)
	$(CXX) -g $(OFILES) $(LPATHS) $(LFLAGS) -o $(EXE)

$(ODIR)/%.o: $(SDIR)/%.$(SEXT)
	$(CXX) -c $< $(CXXFLAGS) $(IPATHS) -o $@

.PHONY: debug release
debug: CXXFLAGS += -Og -g -Wall
release: CXXFLAGS += -O3 -flto
debug release: all

.PHONY: clean
clean: 
	rm $(ODIR)/*.o
