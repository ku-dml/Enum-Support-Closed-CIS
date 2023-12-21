# Makefile

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
	$(RM) *.o *~ \#* ../*~ A1* STATS TRANS* *temp*

OBJS = mylib.o mt19937ar.o cpu_time.o
HEADERS = common.h mylib.h define.h
MAIN = main_alg1.o main_alg1b.o 

all: $(MAIN) common.o $(OBJS) $(HEADERS)
	make A1
	make A1B
	make STATS

A1: main_alg1.o common.o $(OBJS) $(HEADERS) 
	$(CC) -o A1 main_alg1.o common.o $(OBJS) $(CFLAGS)

A1B: main_alg1b.o common.o $(OBJS) $(HEADERS) 
	$(CC) -o A1B main_alg1b.o common.o $(OBJS) $(CFLAGS)

STATS: main_stats.o common.o $(OBJS) $(HEADERS)
	$(CC) -o STATS main_stats.o common.o $(OBJS) $(CFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) -c $<

