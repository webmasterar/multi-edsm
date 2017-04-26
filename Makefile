
CC=     g++

# -O3
CFLAGS= -Wall -g -D_USE_64 -msse4.2 -funroll-loops -fomit-frame-pointer

LFLAGS= -std=c++11 -DNDEBUG -lz -lm -lpthread -I . \
        -I ./libsdsl/include/ \
		-L ./libsdsl/lib/ -lsdsl -ldivsufsort -ldivsufsort64 -Wl,-rpath=$(PWD)/libsdsl/lib \
        -I ./vcflib/tabixpp/ -I ./vcflib/tabixpp/htslib/ -I ./vcflib/smithwaterman/ -I ./vcflib/multichoose/ -I ./vcflib/filevercmp/ -I ./vcflib/src/ \
        -L ./vcflib/ -L ./vcflib/tabixpp/htslib/ -lvcflib -lhts -Wl,-rpath=$(PWD)/vcflib/ -Wl,-rpath=$(PWD)/vcflib/tabixpp/htslib/

EXE=    multiedsm

SRC=    main.cpp MultiEDSM.cpp UnrestrictedMultiShiftAnd.cpp

#
# No need to edit below this line
#

.SUFFIXES:
.SUFFIXES: .cpp .o

OBJ=    $(SRC:.cpp=.o)

.cpp.o:
	$(CC) $(CFLAGS) -c $(LFLAGS) $<

all:    $(EXE)

$(EXE): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(LFLAGS)

$(OBJ): $(MF) $(HD)

clean:
	rm -f $(OBJ) $(EXE) *~

clean-all:
	rm -f $(OBJ) $(EXE) *~
	rm -r libsdsl
	rm -r sdsl-lite
	rm -rf vcflib
