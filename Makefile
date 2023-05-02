# compiler
CXX = g++
CXXFLAGS = -O2 -Wall -std=c++11

# src directories
DIR_BOLEY_ET_AL = ./Boley_et_al
DIR_PVALUE = ./statistic
DIR_STAT = ./statistic/stat
DIR_TEST = ./test

# target directories
TARGET = target
TARGET_DEPS = ${TARGET}/deps

# src file
ifeq ($(shell uname),Linux)
	SRC := ${shell find . -mindepth 1 -maxdepth 1 -name "*.cpp"}
	SRC_BOLEY_ET_AL := ${shell find ${DIR_BOLEY_ET_AL} -mindepth 1 -maxdepth 1 -name "*.cpp"}
	SRC_STAT := ${shell find ${DIR_STAT} -mindepth 1 -maxdepth 1 -name "*.cpp"}
	SRC_PVALUE := ${shell find ${DIR_PVALUE} -mindepth 1 -maxdepth 1 -name "*.cpp"}
	SRC_TEST := ${shell find ${DIR_TEST} -mindepth 1 -maxdepth 1 -name "*.cpp"}
else
	SRC := ${shell find . -depth 1 -name "*.cpp"}
	SRC_BOLEY_ET_AL := ${shell find ${DIR_BOLEY_ET_AL} -depth 1 -name "*.cpp"}
	SRC_STAT := ${shell find ${DIR_STAT} -depth 1 -name "*.cpp"}
	SRC_PVALUE := ${shell find ${DIR_PVALUE} -depth 1 -name "*.cpp"}
	SRC_TEST := ${shell find ${DIR_TEST} -depth 1 -name "*.cpp"}
endif

SRCS = ${SRC} ${SRC_BOLEY_ET_AL} ${SRC_STAT} ${SRC_PVALUE} ${SRC_TEST}
SRCS_LIB = ${SRC} ${SRC_STAT} ${SRC_PVALUE}

# objects files
OBJS = ${foreach src,${SRCS},${patsubst %.cpp,${TARGET_DEPS}/%.o,${notdir ${src}}}}
OBJS_LIB = ${foreach src,${SRCS_LIB},${patsubst %.cpp,${TARGET_DEPS}/%.o,${notdir ${src}}}}
OBJS_TEST = ${foreach src,${SRC_TEST},${patsubst %.cpp,${TARGET_DEPS}/%.o,${notdir ${src}}}}
OBJS_BOLEY_ET_AL = ${foreach src,${SRC_BOLEY_ET_AL},${patsubst %.cpp,${TARGET_DEPS}/%.o,${notdir ${src}}}}

.PHONY: clean test run

run: ${TARGET}/BoleyEtAl
	./$<

build: ${TARGET}/BoleyEtAl

# run test
test: ${TARGET}/PValueTest
	./$<

${TARGET}/BoleyEtAl: ${OBJS_LIB} ${OBJS_BOLEY_ET_AL}
	${CXX} ${CXXFLAGS} -o $@ $^

# make ELF for test
${TARGET}/PValueTest: ${OBJS_LIB} ${OBJS_TEST}
	${CXX} ${CXXFLAGS} -o $@ $^

# compile each src file to object file when file updated
define template
${TARGET_DEPS}/${patsubst %.cpp,%.o,${notdir $1}}: $1 ${TARGET_DEPS}
	if [ ! -e ${TARGET_DEPS}/${patsubst %.cpp,%.o,${notdir $1}} ] \
	|| [ ${TARGET_DEPS}/${patsubst %.cpp,%.o,${notdir $1}} -ot ${1} ]; then \
		${CXX} -c ${CXXFLAGS} -o ${TARGET_DEPS}/${patsubst %.cpp,%.o,${notdir $1}} ${1}; \
	fi
endef

${foreach src,${SRCS},${eval ${call template,${src}}}}

# remove binaries
clean:
	rm -rf ${TARGET}

# make target directory
${TARGET_DEPS}:
	mkdir -p $@
