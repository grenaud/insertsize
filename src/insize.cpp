#include <iostream>
#include <vector>
#include <cstring>
#include <cstring>
#include <set>
#include <ctype.h>
#include <stdlib.h>
#include <sys/mman.h>


#include <sys/types.h>

#include <sys/stat.h>
#include <fcntl.h>
#include "utils.h"

extern "C" {
    //#include "tabix.h"
    //#include "bam.h"
#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include "bam.h"

#include "samtools.h"
#include "sam_opts.h"
#include "bedidx.h"

}


#define bam_is_reverse(b)     (((b)->core.flag&BAM_FREVERSE)    != 0)
#define bam_is_unmapped(b)    (((b)->core.flag&BAM_FUNMAP)      != 0)
#define bam_is_paired(b)      (((b)->core.flag&BAM_FPAIRED)     != 0)
#define bam_is_propaired(b)   (((b)->core.flag&BAM_FPROPER_PAIR) != 0)
#define bam_is_read1(b)       (((b)->core.flag&BAM_FREAD1)      != 0)

#define bam_is_qcfailed(b)    (((b)->core.flag&BAM_FQCFAIL)     != 0)
#define bam_is_rmdup(b)       (((b)->core.flag&BAM_FDUP)        != 0)
#define bam_is_sec(b)         (((b)->core.flag&BAM_FSECONDARY)        != 0)
#define bam_is_supp(b)        (((b)->core.flag&BAM_FSUPPLEMENTARY)    != 0)

#define bam_is_failed(b)      ( bam_is_qcfailed(b) || bam_is_rmdup(b) || bam_is_sec(b) || bam_is_supp(b) )

#define bam_mqual(b)          ((b)->core.qual)
#define bam_isize(b)          ((b)->core.isize)
#define bam_lqseq(b)          ((b)->core.l_qseq)

using namespace std;



int main (int argc, char *argv[]) {

    bool onlyMapped=false;
    bool onlyPP=false;

    string usage=string(""+string(argv[0])+" <options>  [in BAM file]"+
			"\nThis program reads a BAM file and produces the insert sizes\n"+
			"\n"+

			// "\nreads and the puts the rest into another bam file.\n"+
			// "\nTip: if you do not need one of them, use /dev/null as your output\n"+

			"\n\n\tOther options:\n"+
			"\t\t"+"-m\t\t\tRequire the reads to be mapped (Default: "+booleanAsString( onlyMapped )+")\n"+
			"\t\t"+"-p\t\t\tRequire the reads to be properly paired (Default: "+booleanAsString( onlyPP )+")\n"+
			"\n");

    if(argc == 1 ||
       (argc == 2 && (string(argv[0]) == "-h" || string(argv[0]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }

    
    for(int i=1;i<(argc-1);i++){ //all but the last 3 args


        if(string(argv[i]) == "-m"  ){
            onlyMapped=true;
            continue;
        }

        if(string(argv[i]) == "-p"  ){
            onlyPP=true;
            continue;
        }

	cerr<<"Error: unknown option "<<string(argv[i])<<endl;
	return 1;
    }

    
    string bamfiletopen = string( argv[ argc-1 ] );

    samFile  *fp;
    bam1_t    *b;
    bam_hdr_t *h;
    
    fp = sam_open_format(bamfiletopen.c_str(), "r", NULL); 
    if(fp == NULL){
	cerr << "Could not open input BAM file"<< bamfiletopen << endl;
	return 1;
    }

    h = sam_hdr_read(fp);
    if(h == NULL){
	cerr<<"Could not read header for "<<bamfiletopen<<endl;
	return 1;
    }
    b = bam_init1();
    while(sam_read1(fp, h, b) >= 0){
	if(bam_is_failed(b) )          continue;
	
	if(bam_is_unmapped(b) ){
	    //if the read is unmapped and we only consider mapped reads
	    if(onlyMapped){  continue; }
	}else{//read is mapped	    
	    //we accept
	}

	//if(b->core.l_qseq < minLength) continue;
	bool ispaired    = bam_is_paired(b);
	bool isfirstpair = bam_is_read1(b);
	
	if(ispaired){	    
	    if( isfirstpair   ){
		if( !bam_is_propaired(b) ){
		    if(onlyPP){  continue; }//if not properly paired, skip		    
		}

		int32_t isize = bam_isize(b);
		if( isize == 0) continue; //from different chromosomes
		if( isize > 0)
		    cout<<isize<<endl;
		else
		    cout<<-1.0*isize<<endl;
	    }else{
		//ignore
	    }
	}else{
	    cout<<bam_lqseq(b)<<endl;
	}

	
    }
    
    bam_destroy1(b);
    sam_close(fp);
    
   
    return 0;
}

