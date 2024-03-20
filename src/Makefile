
CXX      = g++ -std=c++11# -g -pg #-fpermissive #-g # -pg
LIBGAB   = ../lib/libgab/

CXXFLAGS = -Wall -lm -O3 -lz -I../lib/libgab/ -I../lib/libgab//gzstream/ -I../lib/samtools/ -I../lib/htslib/    -c
LDFLAGS  = -lgsl -lcurl -lgslcblas -llzma -lbz2 -lz  -lpthread 
LDLIBS   = ../lib/samtools/bedidx.o ../lib/htslib/libhts.a ../lib/samtools/libst.a ../lib/samtools/libst.a 


all:  ../lib/libgab/libgab.a  ../lib/htslib/libhts.a ../lib/samtools/bedidx.o insize

%.o: %.cpp #../lib/libgab/utils.o  
	${CXX} ${CXXFLAGS} $^ -o $@

insize:	insize.o   ../lib/libgab/libgab.a   ../lib/libgab/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

static:	insize.o   ../lib/libgab/libgab.a   ../lib/libgab/gzstream/libgzstream.a
	${CXX} -static -o $@ $^ $(LDLIBS) $(LDFLAGS)

../lib/htslib/libhts.a:  ../lib/htslib/hts_internal.h
	cd ../lib/htslib/ && git submodule update --init --recursive && autoreconf -i && ./configure  --disable-libcurl && make && cd ../..

../lib/htslib/hts_internal.h:
	rm -rfv ../lib/htslib/
	git clone --recursive https://github.com/samtools/htslib.git ../lib/htslib/

../lib/samtools/bedidx.o:  ../lib/samtools/sam.h
	cd ../lib/samtools/  && autoreconf -i && ./configure  --disable-libcurl && make && cd ../..

../lib/samtools/sam.h:
	rm -rf ../lib/samtools/
	git clone --recursive https://github.com/samtools/samtools.git ../lib/samtools/

../lib/libgab/libgab.a:  ../lib/libgab/libgab.h
	cd ../lib/libgab/ &&  make libgab.a && make -C gzstream/ && cd ../..

../lib/libgab/libgab.h:
	rm -rfv ../lib/libgab/
	mkdir ../lib/libgab/
	git clone --depth 1 https://github.com/grenaud/libgab.git ../lib/libgab/


clean :
	rm -f *.o bam2prof diffprof


