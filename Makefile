CXXFLAGS := -std=c++11 -O2 -g -fPIC

all:	libSXrdClasses.so analX

%.o : %.cxx %.h
	g++ ${CXXFLAGS} -c -o $@ `root-config --cflags` $<

SXrdClasses_Dict.h SXrdClasses_Dict.cxx: SXrdClasses.h SXrdClasses_LinkDef.h
	rootcint -f SXrdClasses_Dict.cxx -c -p $^

libSXrdClasses.so: SXrdClasses.o SXrdClasses_Dict.o
	g++ ${CXXFLAGS} -shared -o $@ `root-config --cflags` $^

count_stuff: count_stuff.cxx deep_dump.cxx libSXrdClasses.so
	g++ `root-config --cflags --libs` -Wl,-rpath=. -o count_stuff $^

wisc_anal: wisc_anal.cxx libSXrdClasses.so
	g++ ${CXXFLAGS} -o $@ -Wl,-rpath=. `root-config --cflags --libs` $^

ucsd_anal: ucsd_anal.cxx deep_dump.cxx libSXrdClasses.so
	g++ ${CXXFLAGS} -o $@ -Wl,-rpath=. `root-config --cflags --libs` $^

ANALH := $(wildcard Anal*.h)   $(wildcard AnFi*.h)   $(wildcard AnEx*.h)
ANALS := $(wildcard Anal*.cxx) $(wildcard AnFi*.cxx) $(wildcard AnEx*.cxx)
ANALO := $(ANALS:%.cxx=%.o)

# Poor man's dependency faker: rebuild all on any header change
${ANALO}: ${ANALH}

analX: ${ANALO} libSXrdClasses.so
	g++ ${CXXFLAGS} -o $@ -Wl,-rpath=. `root-config --cflags --libs` $^

clean:
	rm -f *.o *rdict.pcm
	rm -f SXrdClasses_Dict.* libSXrdClasses.so
	rm -f wisc_anal ucsd_anal analX count_stuff
