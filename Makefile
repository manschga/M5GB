COMPILER=GCC_

SOURCES = main.cpp other_fct.cpp lab_poly.cpp M5GB.cpp poly.cpp signature.cpp 

OBJECTS = $(SOURCES:.cpp=.o)

DEPS = $(SOURCES:.cpp=.h)

PROGRAM	= main.${COMPILER}

PROGRAM_PROFILE	= main_profile.${COMPILER}

CXXFLAGS += -g -O3 #-pg

CXX = g++

LINKER = ${CXX}

ARCHIVE = elements_archive.tar

REFERENCE_COMPILE = compile.log

REFERENCE_FILE = run.log

CXXFLAGS_PROFILE = -g -pg

PROFILE_FILE = profile_data.txt


default: ${PROGRAM}

${PROGRAM}: ${OBJECTS}
	${LINKER} -o $@ $^ ${CXXFLAGS}

before:
	./touch_cpp.sh


%.o: %.cpp
	${LINKER} -c -o $@ $< ${CXXFLAGS}


log:
	make clean_all
	make default 2>&1 | tee ${REFERENCE_COMPILE}
	./${PROGRAM} 2>&1 | tee ${REFERENCE_FILE}

clean:
	@rm -f ${PROGRAM} ${OBJECTS}

clean_all: clean
	@rm -f *~ *.bak *.out *.aux *.orig *.log
	@rm -rf html bin obj

run:
	make clean
	make
	./${PROGRAM}

delete:
	make clean_all
	@rm -f *.cpp *.h *.${COMPILER} *.o *.m *.txt *.sh Doxyfile
	#@rm -r html latex

testing:
	make delete
	tar -xf ${ARCHIVE}
	make default 2>&1 | tee compile_.log
	./${PROGRAM} 2>&1 | tee run_.log
	./diff_check.sh run_.log ${REFERENCE_FILE}
	#make delete

profile: ${PROGRAM}
	valgrind --tool=callgrind ./$^

mem: ${PROGRAM}
	valgrind -v --leak-check=yes --tool=memcheck --undef-value-errors=yes --track-origins=yes --log-file=$^.addr.out --show-reachable=yes ./$^

${PROGRAM_PROFILE}: ${OBJECTS}
	${LINKER} -o $@ $^ ${CXXFLAGS_PROFILE}

gprof: ${PROGRAM_PROFILE}
	./${PROGRAM_PROFILE}
	gprof -f isDivisible ${PROGRAM_PROFILE} gmon.out > ${PROFILE_FILE}
	cat ${PROFILE_FILE}

doc:
	doxygen Doxyfile

copy:
	echo "copy folder"
	scp -r ${PWD} scicomp_19@143.50.47.149:Test_Hauke

sync:
	rsync -a --delete ${PWD} scicomp_19@143.50.47.149:Test_Hauke

push:
	git add ${PWD}
	git commit -m "."
	git push

pop:
	git pull



