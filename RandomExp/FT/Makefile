# Makefile
##### openmp #####
ifeq ($(shell uname),Darwin)
	CXX           = gcc-9
	LINK          = gcc-9
	LFLAGS = -Wall
	LIBS          = $(SUBLIBS) /usr/local/Cellar/gcc/9.1.0/lib/gcc/9/libgomp.1.dylib -lpthread -lc++ -lstdc++ 
else
	CXX           = g++
	LINK          = g++
	LFLAGS        = -Wl,-O1
	LIBS          = $(SUBLIBS) /usr/lib/x86_64-linux-gnu/libgomp.so.1 -lpthread
endif
CXXFLAGS = -pipe -fopenmp -O2 -std=gnu++11 -Wall -W -D_REENTRANT -fPIC #-D_OPENMP


# suffixes.
.SUFFIXES: .c .cpp .o 

# macros.
CC = g++
CFLAGS = -O2 -Wall -std=c++11
RM = rm -f

# cleaning rules.
clean:
	$(RM) *.o *~ \#* ../*~

distclean:
	$(RM) *.o *~ \#* ../*~ TREE2 PARALLEL *temp*

OBJS = mylib.o mt19937ar.o cpu_time.o
HEADERS = common.h mylib.h define.h

TREE2: main.o common.o $(OBJS) $(HEADERS) 
	$(CC) -o TREE2 main.o common.o $(OBJS) $(CFLAGS)

PARALLEL: main_omp.o common_omp.o $(OBJS) $(HEADERS) 
	$(LINK) $(LFLAGS) -o PARALLEL  main_omp.o common_omp.o $(OBJS) $(LIBS)

main_omp.o: main.cpp
	$(CXX) -c $(CXXFLAGS) -o main_omp.o main.cpp

common_omp.o: common.cpp
	$(CXX) -c $(CXXFLAGS) -o common_omp.o common.cpp


.cpp.o:
	$(CC) $(CFLAGS) -c $<

