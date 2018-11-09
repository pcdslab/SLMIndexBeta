BUILD = RELEASE
CXX = g++
MSTOOLKIT = MSToolkit
CometSearch = CometSearch
LIBDIVSUFSORT = CometSearch/libdivsufsort

ifeq (${BUILD}, DEBUG)
 BUILD_FLAGS = -O0 -g -Wno-unused-result
else
 BUILD_FLAGS = -O3 -g -Wno-unused-result
endif

export BUILD_FLAGS

override CXXFLAGS += ${BUILD_FLAGS} -Wall  -Wextra -Wno-char-subscripts -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D__LINUX__ -I$(MSTOOLKIT)/include -I$(CometSearch)
OUTPUT_DIR = output
EXECNAME = comet.exe
OBJS = Comet.o
DEPS = CometSearch/CometData.h CometSearch/CometDataInternal.h CometSearch/CometPreprocess.h CometSearch/CometWriteOut.h CometSearch/CometWriteSqt.h CometSearch/OSSpecificThreading.h CometSearch/CometMassSpecUtils.h CometSearch/CometSearch.h CometSearch/CometWritePepXML.h CometSearch/CometWriteTxt.h CometSearch/Threading.h CometSearch/CometPostAnalysis.h CometSearch/SLMIndex.h CometSearch/CometSearchManager.h CometSearch/CometWritePercolator.h CometSearch/Common.h CometSearch/ThreadPool.h CometSearch/CometMassSpecUtils.cpp CometSearch/CometSearch.cpp CometSearch/CometWritePepXML.cpp CometSearch/CometWriteTxt.cpp CometSearch/CometPostAnalysis.cpp CometSearch/SLMIndex.cpp CometSearch/CometSearchManager.cpp CometSearch/CometWritePercolator.cpp CometSearch/Threading.cpp CometSearch/CometPreprocess.cpp CometSearch/CometWriteOut.cpp CometSearch/CometWriteSqt.cpp

LIBPATHS = -L$(MSTOOLKIT) -L$(CometSearch)
LIBS = -lcometsearch -lmstoolkitlite -lm -lpthread -ldivsufsort
ifdef MSYSTEM
   LIBS += -lws2_32
endif

comet.exe: $(OBJS)
	cd $(MSTOOLKIT) ; make lite ; cd ../CometSearch ; make
	mkdir -p ${OUTPUT_DIR}/bin
	${CXX} -fopenmp $(CXXFLAGS) -Wl,-Map=${OUTPUT_DIR}/bin/${EXECNAME}.map ${OUTPUT_DIR}/objs/$(OBJS) $(LIBPATHS) $(LIBS) -o ${OUTPUT_DIR}/bin/${EXECNAME} 

Comet.o: Comet.cpp $(DEPS)
	mkdir -p ${OUTPUT_DIR}/objs
	${CXX} ${CXXFLAGS} Comet.cpp -c -o ${OUTPUT_DIR}/objs/$@

clean:
	rm -rf *.o ./${OUTPUT_DIR}
	cd $(MSTOOLKIT) ; make clean ; cd ../CometSearch ; make clean
