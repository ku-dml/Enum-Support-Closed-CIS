# compiler
CXX = g++
CXXFLAGS = -O2 -Wall -std=c++11

# src directories
DIR_CLOST = ./clost
DIR_PVALUE = ./statistic
DIR_STAT = ./statistic/stat

# target directories
TARGET = target
TARGET_DEPS = ${TARGET}/deps

# src file
SRC := ${shell find . -depth 1 -name "*.cpp"}
SRC_CLOST := ${shell find ${DIR_CLOST} -depth 1 -name "*.cpp"}
SRC_STAT := ${shell find ${DIR_STAT} -depth 1 -name "*.cpp"}
SRC_PVALUE := ${shell find ${DIR_PVALUE} -depth 1 -name "*.cpp"}
SRCS = ${SRC} ${SRC_CLOST} ${SRC_STAT} ${SRC_PVALUE}

# objects files
OBJS = ${foreach src,${SRCS},${patsubst %.cpp,${TARGET_DEPS}/%.o,${notdir ${src}}}}

.PHONY: clean test

# run test
test_stat: ${TARGET}/PValue
	./$<

# make ELF
${TARGET}/PValue: ${OBJS}
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
