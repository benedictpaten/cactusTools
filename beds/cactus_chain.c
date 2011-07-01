#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>

#include "cactus.h"
#include "avl.h"
#include "commonC.h"
#include "hashTableC.h"
//#include "cactus_addReferenceSeq.h"

/*
 *Feb 7 2011: modify to include the referene sequence in the flower
 *Jan 19 2011: Correct to only put chain-segments in the same bed-record if there are links between them
 *(any two adjacent segments in the record must belong to two different blocks). 
 *Start a new bed record if hit an end with a self-edge
 *
 */

char *formatSequenceHeader(Sequence *sequence) {
    const char *sequenceHeader = sequence_getHeader(sequence);
    if(strlen(sequenceHeader) > 0) {
        char *cA = st_malloc(sizeof(char) *(1 + strlen(sequenceHeader)));
        sscanf(sequenceHeader, "%s", cA);
        return cA;
    }
    else {
        return cactusMisc_nameToString(sequence_getName(sequence));
    }
}

bool isStubCap(Cap *cap){
    /*
     *Return true if cap is of a stubEnd, otherwise return false
     */
    assert(cap != NULL);
    End *end = cap_getEnd(cap);
    assert(end != NULL);
    return (end_isStubEnd(end)) ? true : false;
}

bool isLinked(End *end1, End *end2){
    //Return true if there is a link between end1 and end2, otherwise return false
    Link *link = group_getLink(end_getGroup(end1));
    if (link != NULL){
        End *link3end = link_get3End(link);
        End *link5end = link_get5End(link);
        end1 = end_getPositiveOrientation(end1);
        end2 = end_getPositiveOrientation(end2);
        if ( (link3end == end1 && link5end == end2) || (link3end == end2 && link5end == end1) ){
            return true;
        }
    }
    return false;
}

struct List *flower_getThreadStarts(Flower *flower, char *name){
    /*
     *Get 3' end Stubs of the sequence by its name
     */
    Cap *cap;
    struct List *startCaps = constructEmptyList(0, free); 
    Flower_CapIterator *capIterator = flower_getCapIterator(flower);
    while((cap= flower_getNextCap(capIterator)) != NULL){
        if(isStubCap(cap)){//dead end or inherited end
            Sequence *sequence = cap_getSequence(cap);
            if(sequence == NULL){continue;}
            char *sequenceHeader = formatSequenceHeader(sequence);
            if(strcmp(sequenceHeader, name) == 0){//cap matched with name
                if(!cap_getStrand(cap)){//if cap is on negative strand - reverse it
                    cap = cap_getReverse(cap);
                }
                if(cap_getSide(cap)){ continue; }//if cap is the 5' end, ignore
                listAppend(startCaps, cap);
            }
            free(sequenceHeader);
        }
    }
    flower_destructCapIterator(capIterator);
    return startCaps;
}

void moveCapToNextBlock(Cap **cap){
    /*Move cap to the next block that its segment is adjacency to*/
    assert((*cap) != NULL);
    st_logInfo("cap %s, %d\t", cactusMisc_nameToString(cap_getName(*cap)), cap_getCoordinate(*cap));
    if(isStubCap(*cap)){
        st_logInfo("STUB\n");
        return;
    }
    Cap *otherCap = cap_getOtherSegmentCap(*cap);
    st_logInfo("other cap %s, %d\t", cactusMisc_nameToString(cap_getName(otherCap)), cap_getCoordinate(otherCap));
    *cap = cap_getAdjacency(otherCap);
    assert(*cap != NULL);
    st_logInfo("moved-cap %s, %d\n", cactusMisc_nameToString(cap_getName(*cap)), cap_getCoordinate(*cap));
    
    /*Cap *adjCap = cap_getAdjacency(*cap);
    assert(adjCap != NULL);
    if(cap_getEnd(*cap) == cap_getEnd(adjCap)){//DOUBLE CHECK ... self connected end
        *cap = adjCap;
    }else{
        *cap = cap_getAdjacency(cap_getOtherSegmentCap(*cap));
    //}*/
}

Cap *end_getCapBySeqName(End *end, char *name){
    Cap *cap;
    End_InstanceIterator *it = end_getInstanceIterator(end);
    while( (cap = end_getNext(it)) != NULL ){
        Sequence *sequence = cap_getSequence(cap);
        if(sequence == NULL){continue;}
        char *sequenceHeader = formatSequenceHeader(sequence);
        st_logInfo("%s\t%d\n", sequenceHeader, cap_getCoordinate(cap));
        if(strstr(sequenceHeader, name) != NULL){//cap matched with name
            break;
        }
    }
    end_destructInstanceIterator(it);
    return cap;
}

int32_t cap_getCoor(Cap *cap){
    int32_t coor = cap_getCoordinate(cap);
    Sequence *seq = cap_getSequence(cap); 
    assert(seq != NULL);
    
    if (cap_getStrand(cap)){//pos strand
        coor -= sequence_getStart(seq);
    }else{
        //coor =  sequence_getLength(seq) - 1 - coor;
        coor = (sequence_getStart(seq) + sequence_getLength(seq) - 1) - coor;
    }
    return coor;
}

int32_t getReverseCoor(int32_t coor, int32_t length){
    return length -1 - coor;
}

int32_t cap_getSeqSize(Cap *cap){
    Sequence *seq = cap_getSequence(cap);
    assert(seq != NULL);
    return sequence_getLength(seq);
}

char *cap_getChr(Cap *cap){
    char *chr;
    Sequence *sequence = cap_getSequence(cap);
    if(sequence == NULL){return NULL;}
    char *seqname = formatSequenceHeader(sequence);
    char sep[] = ".";
    
    strtok(stString_copy(seqname), sep); //query e.g "panTro2"
    chr = strtok(NULL, sep);
    return chr;
}
    
int32_t convertCoor(int32_t coor, int32_t chrsize, int32_t offset, int32_t fragmentsize, char strand){
    if(strand == '+'){
        coor = coor + offset; 
    }else{
        coor = coor + chrsize - (offset + fragmentsize);
    }
    return coor;
}

//void block_getCHAIN(Block *block, FILE *fileHandle, char *query, char *target, char *qchr, int32_t start, int32_t qchrsize, int level) {
int block_getCHAIN(Block *block, FILE *fileHandle, char *query, char *target, char *qchr, int32_t start, int32_t qchrsize, int32_t chainid) {
    /*
     */
    Block_InstanceIterator *instanceIterator = block_getInstanceIterator(block);
    Segment *segment;
    while((segment = block_getNext(instanceIterator)) != NULL) {
        Sequence *sequence = segment_getSequence(segment);
        if(sequence != NULL) {
            char *sequenceHeader = formatSequenceHeader(sequence);
            if(strcmp(sequenceHeader, query) == 0){
                Cap *tcap5 = end_getCapBySeqName(cap_getEnd(segment_get5Cap(segment)), target); 
                if(!cap_getStrand(tcap5)){//if reference segment is on the neg strand, reverse it
                    tcap5 = cap_getOtherSegmentCap( cap_getReverse(tcap5) );
                    segment = segment_getReverse(segment);
                }
                Cap *cap5 = segment_get5Cap(segment);
                Cap *cap3 = segment_get3Cap(segment);
                //int s = cap_getCoordinate(cap5) + start - 2;
                //int e = cap_getCoordinate(cap3) + start + 1 - 2;
                int32_t qstart = cap_getCoor(cap5);
                int32_t qend = cap_getCoor(cap3) + 1;
                int32_t tstart = cap_getCoor(tcap5);
                int32_t tend = cap_getCoor(cap_getOtherSegmentCap(tcap5)) + 1;
                int32_t tsize = cap_getSeqSize(tcap5);
                int32_t qsize = cap_getSeqSize(cap5);
                char *tchr = cap_getChr(tcap5);
                char tstrand = cap_getStrand(tcap5) ? '+' : '-';
                char qstrand = cap_getStrand(cap5) ? '+' : '-';


                qstart = convertCoor(qstart, qchrsize, start, qsize, qstrand);
                qend = convertCoor(qend, qchrsize, start, qsize, qstrand);
                //fprintf(fileHandle, "%s %d %d %s.%d %d %s %d %d %d %d %d %d\n", chr, s, e, "NA", level, 0, ".", s, e, 0, 1, e-s,0);
                //fprintf(fileHandle, "chain %d %s %d %c %d %d %s %d %c %d %d %d\n", 0, tchr, tsize, tstrand, tstart, tend, qchr, qsize, qstrand, qstart, qend, level);
                fprintf(fileHandle, "chain %d %s %d %c %d %d %s %d %c %d %d %d\n", 0, tchr, tsize, tstrand, tstart, tend, qchr, qchrsize, qstrand, qstart, qend, chainid);
                chainid ++;
                fprintf(fileHandle, "%d\n\n", segment_getLength(segment));
                //fprintf(fileHandle, "%s %d %d %d %d %s %d %d %d %d %d %d\n", chr, s, e, level, 0, "+", s, e, 0, 1, e-s,0);
            }
            free(sequenceHeader);
        }
    }
    block_destructInstanceIterator(instanceIterator);
    return chainid;
}


//void chain_getCHAINs(Chain *chain, Cap *qcap, FILE *fileHandle, char *query, char *target, char *qchr, int32_t start, int32_t qchrsize, int level) {
int32_t chain_getCHAINs(Chain *chain, Cap *qcap, FILE *fileHandle, char *query, char *target, char *qchr, int32_t start, int32_t qchrsize, int32_t chainid) {
    /*
     */
    st_logInfo("chain_getCHAINs\n");
    char qstrand = '.';
    char tstrand = '.';
    char *tchr = "";
    //char *chainName = cactusMisc_nameToString(chain_getName(chain));
    int32_t qstart = 0;
    int32_t tstart = 0;
    int32_t qend = -1;
    int32_t tend = -1;
    int32_t qsize = 0;
    int32_t tsize = 0;
    int32_t dq = 0;
    int32_t dt = 0;
    int32_t size = 0;
    struct IntList *blockSizes = constructEmptyIntList(0);
    struct IntList *dtList = constructEmptyIntList(0);
    struct IntList *dqList = constructEmptyIntList(0);
    Block *block;
    End *currEnd;
    End *prevOtherEnd = NULL;
    //End *nextEnd = NULL;
    Cap *othercap;
    Cap *tcap;
    //Cap *tprevCap = NULL;
    //Cap *qprevCap = NULL;

    if(isStubCap(qcap)){ //startCap of the current flower's thread
        qcap = cap_getAdjacency(qcap);
    }//else: current cap is somewhere in the middle of the thread
    
    if(qcap == NULL){ return chainid; }
    while(!isStubCap(qcap)){
        othercap = cap_getOtherSegmentCap(qcap);
        currEnd = cap_getEnd(qcap);
        block = end_getBlock(currEnd);
        Chain *currchain = block_getChain(block);
        st_logInfo("\n\nqcap: %s\t%d\n", query, cap_getCoordinate(qcap));
        tcap = end_getCapBySeqName(currEnd, target);
        if(tcap == NULL){ continue; } //

        if (chain == currchain){
            if( prevOtherEnd == NULL || isLinked(currEnd, prevOtherEnd) ){
                size = cap_getCoordinate(othercap) - cap_getCoordinate(qcap) + 1;
                intListAppend(blockSizes, size);
                
                //if(qprevCap == NULL){
                if(qend == -1){
                    qstrand = cap_getStrand(qcap) ? '+' : '-';
                    tstrand = cap_getStrand(tcap) ? '+' : '-';
                    tstart = cap_getCoor(tcap);
                    qstart = cap_getCoor(qcap);
                    tsize = cap_getSeqSize(tcap);
                    qsize = cap_getSeqSize(qcap);
                    tchr = cap_getChr(tcap);
                }else{
                    assert( (qstrand == '+' && cap_getStrand(qcap)) || (qstrand == '-' && !cap_getStrand(qcap)) );
                    assert( (tstrand == '+' && cap_getStrand(tcap)) || (tstrand == '-' && !cap_getStrand(tcap)) );
                    dt = abs(cap_getCoor(tcap) - tend);
                    dq = abs(cap_getCoor(qcap) - qend);
                    intListAppend(dtList, dt);
                    intListAppend(dqList, dq);
                }
                qend = cap_getCoor(othercap) + 1;
                tend = cap_getCoor(cap_getOtherSegmentCap(tcap)) + 1;
            }else{
                chainid = chain_getCHAINs(chain, qcap, fileHandle, query, target, qchr, start, qchrsize, chainid);
                break;
            }
        }
 
        //st_logInfo("cap %s, %d\t", cactusMisc_nameToString(cap_getName(cap)), cap_getCoordinate(cap));
        prevOtherEnd = end_getOtherBlockEnd(currEnd);
        moveCapToNextBlock(&qcap);
        //st_logInfo("movedTo %s, %d\n", cactusMisc_nameToString(cap_getName(cap)), cap_getCoordinate(cap));
    }
    
    if(blockSizes->length == 0){
       //fprintf(fileHandle, "No block found for chain %s\n", chainName);
        return chainid;
    }

    //Convert query coordinates to assembly coordinates:
    qstart = convertCoor(qstart, qchrsize, start, qsize, qstrand);    
    qend = convertCoor(qend, qchrsize, start, qsize, qstrand);    

    int score = 0;
    //chromStart = chromStart + start -1;
    //chromEnd = chromEnd + start -1;
    //fprintf(fileHandle, "chain %d %s %d %c %d %d %s %d %c %d %d %s\n", score, tchr, tsize, tstrand, tstart, tend, qchr, qchrsize, qstrand, qstart, qend, cactusMisc_nameToString(flower_getName(chain_getFlower(chain))) );

    if(tstrand == '+'){
        fprintf(fileHandle, "chain %d %s %d %c %d %d %s %d %c %d %d %d\n", score, tchr, tsize, tstrand, tstart, tend, qchr, qchrsize, qstrand, qstart, qend, chainid);
        //Print blockSizes
        for(int i=0; i< dtList->length; i++){
            fprintf(fileHandle, "%d %d %d\n", blockSizes->list[i], dtList->list[i], dqList->list[i]);
        }
        fprintf(fileHandle, "%d\n", blockSizes->list[blockSizes->length - 1]);

    }else{ //reverse the chain so that tstrand is '+'
        tstrand = '+';
        qstrand = qstrand == '+' ? '-' : '+';
        int32_t temp = qend;
        qend = getReverseCoor(qstart, qchrsize) + 1;
        qstart = getReverseCoor(temp, qchrsize) + 1; 
        temp = tend;
        tend = getReverseCoor(tstart, tsize) + 1; 
        tstart = getReverseCoor(temp, tsize) + 1;
        fprintf(fileHandle, "chain %d %s %d %c %d %d %s %d %c %d %d %d\n", score, tchr, tsize, tstrand, tstart, tend, qchr, qchrsize, qstrand, qstart, qend, chainid);

        for(int i= dtList->length - 1; i >= 0; i--){
            fprintf(fileHandle, "%d %d %d\n", blockSizes->list[i + 1], dtList->list[i], dqList->list[i]);
        }
        fprintf(fileHandle, "%d\n", blockSizes->list[0]);
    }

    fprintf(fileHandle, "\n");
    chainid ++;
    return chainid;
}


int32_t getSeqLength(Flower *flower, char *header){
    Flower_SequenceIterator *it = flower_getSequenceIterator(flower);
    Sequence *sequence;
    while((sequence = flower_getNextSequence(it)) != NULL){
        char *sequenceHeader = formatSequenceHeader(sequence);
        if(strcmp(sequenceHeader, header) == 0){
            flower_destructSequenceIterator(it);
            return sequence_getLength(sequence);
        }
    }
    flower_destructSequenceIterator(it);
    return 0;
}



//void getBEDs(Flower *flower, FILE *fileHandle, char *species, int level){
int getCHAINs(Flower *flower, FILE *fileHandle, char *query, char *target, int32_t chainid){
    char *qchr;
    //char *tchr;
    char *tok;
    int32_t qstart = 0;
    int32_t qchrsize = 0;
    //int32_t tstart = 0;
    //int start = 1;
    char sep[] = ".";

    //get chrom and start for the query
    assert(query != NULL);
    strtok(stString_copy(query), sep); //query e.g "panTro2"
    qchr = strtok(NULL, sep);
    if(qchr == NULL){
        qchr = "";
    }else{
        tok = strtok(NULL, sep);//chromsize
        if(tok != NULL){
            sscanf(tok, "%d", &qchrsize);
            sscanf(strtok(NULL, sep), "%d", &qstart);
        }
    }

    //sscanf(strtok(NULL, sep), "%d", &len);
    //sscanf(strtok(NULL, sep), "%d", &strand);
    //fprintf(stderr, "chrom *%s*, chromsize %d, start %d\n", chr, chrsize, start);

    //Get beds for all chain of current level:
    struct List *startCaps;
    startCaps = flower_getThreadStarts(flower, query);
    st_logInfo("Number of start Caps at flower %s is %d\n",
                cactusMisc_nameToString(flower_getName(flower)), startCaps->length);
    for(int i=0; i< startCaps->length; i++){
        Cap *startcap = startCaps->list[i];
        Flower_ChainIterator *chainIterator = flower_getChainIterator(flower);
        Chain *chain;
        while((chain = flower_getNextChain(chainIterator)) != NULL){
            chainid = chain_getCHAINs(chain, startcap, fileHandle, query, target, qchr, qstart, qchrsize, chainid);
        
        }
        flower_destructChainIterator(chainIterator);
    }

    //Get beds for non-trivial chains
    Flower_BlockIterator *blockIterator = flower_getBlockIterator(flower);
    Block *block;
    while((block = flower_getNextBlock(blockIterator)) != NULL){
        if(block_getChain(block) == NULL){//non-trivial chain
            chainid = block_getCHAIN(block, fileHandle, query, target, qchr, qstart, qchrsize, chainid);
        }
    }
    flower_destructBlockIterator(blockIterator);

    //Call child flowers recursively.
    Flower_GroupIterator *groupIterator = flower_getGroupIterator(flower);
    Group *group;
    //level ++;
    //fprintf(fileHandle, "Go to level %d\n", level);
    while((group = flower_getNextGroup(groupIterator)) != NULL) {
        Flower *nestedFlower = group_getNestedFlower(group);
        if(nestedFlower != NULL) {
            //fprintf(fileHandle, "level %d\n", level);
            chainid = getCHAINs(group_getNestedFlower(group), fileHandle, query, target, chainid); //recursive call.
        }
    }
    flower_destructGroupIterator(groupIterator);
    return chainid;
}

struct List *getSequences(Flower *flower, char *name){
   //get names of all the sequences in 'flower' that have their names start with 'name'
   Sequence *sequence;
   struct List *seqs = constructEmptyList(0, free);
   Flower_SequenceIterator * seqIterator = flower_getSequenceIterator(flower);
   while((sequence = flower_getNextSequence(seqIterator)) != NULL){
      Event *event = sequence_getEvent(sequence);
      //char *sequenceHeader = formatSequenceHeader(sequence);
      char *sequenceHeader = stString_copy(event_getHeader(event));
      //if 'sequenceHeader' starts with 'name'
      //if(strstr(sequenceHeader, name) == sequenceHeader){
      if(strstr(sequenceHeader, name) != NULL){
         listAppend(seqs, sequenceHeader);
      }
      //free(sequenceHeader);
   }
   flower_destructSequenceIterator(seqIterator);
   return seqs;
}

void usage() {
    fprintf(stderr, "cactus_bedGenerator, version 0.2\n");
    fprintf(stderr, "Prints to output file all segments of the target sequence that are in blocks that contain both query & target\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr, "-q --query\n");
    fprintf(stderr, "-t --target: If target = 'reference' then will build the reference sequence named 'reference'\n");
    fprintf(stderr, "-c --cactusDisk : The cactus database conf string\n");
    fprintf(stderr, "-d --flowerName : The name of the flower (the key in the database)\n");
    fprintf(stderr, "-e --outputFile : The file to write the BEDs in.\n");
    fprintf(stderr, "-h --help : Print this help screen\n");
}

int main(int argc, char *argv[]) {
    CactusDisk *cactusDisk;
    Flower *flower;

    /*
     * Arguments/options
     */
    char * logLevelString = NULL;
    char * cactusDiskDatabaseString = NULL;
    char * flowerName = NULL;
    char * outputFile = NULL;
    char * query = NULL;
    char * target = NULL;

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    while(1) {
        static struct option long_options[] = {
            { "logLevel", required_argument, 0, 'a' },
            //{ "species", required_argument, 0, 'b' },
            { "query", required_argument, 0, 'q' },
            { "target", required_argument, 0, 't' },
            { "cactusDisk", required_argument, 0, 'c' },
            { "flowerName", required_argument, 0, 'd' },
            { "outputFile", required_argument, 0, 'e' },
            { "help", no_argument, 0, 'h' },
            { 0, 0, 0, 0 }
        };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:q:t:c:d:e:h", long_options, &option_index);

        if(key == -1) {
            break;
        }

        switch(key) {
            case 'a':
                logLevelString = stString_copy(optarg);
                break;
            case 'q':
                query = stString_copy(optarg);
                break;
            case 't':
                target = stString_copy(optarg);
                break;
            case 'c':
                cactusDiskDatabaseString = stString_copy(optarg);
                break;
            case 'd':
                flowerName = stString_copy(optarg);
                break;
            case 'e':
                outputFile = stString_copy(optarg);
                break;
            case 'h':
                usage();
                return 0;
            default:
                usage();
                return 1;
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    // (0) Check the inputs.
    ///////////////////////////////////////////////////////////////////////////

    assert(cactusDiskDatabaseString != NULL);
    assert(flowerName != NULL);
    assert(outputFile != NULL);
    assert(query != NULL);
    assert(target != NULL);

    //////////////////////////////////////////////
    //Set up logging
    //////////////////////////////////////////////

    if(logLevelString != NULL && strcmp(logLevelString, "INFO") == 0) {
        st_setLogLevel(ST_LOGGING_INFO);
    }
    if(logLevelString != NULL && strcmp(logLevelString, "DEBUG") == 0) {
        st_setLogLevel(ST_LOGGING_DEBUG);
    }

    //////////////////////////////////////////////
    //Log (some of) the inputs
    //////////////////////////////////////////////

    st_logInfo("Flower name : %s\n", flowerName);
    st_logInfo("Output BED file : %s\n", outputFile);

    //////////////////////////////////////////////
    //Load the database
    //////////////////////////////////////////////

    stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(cactusDiskDatabaseString);
    cactusDisk = cactusDisk_construct(kvDatabaseConf, 0);
    st_logInfo("Set up the flower disk\n");

    ///////////////////////////////////////////////////////////////////////////
    // Parse the basic reconstruction problem
    ///////////////////////////////////////////////////////////////////////////

    flower = cactusDisk_getFlower(cactusDisk, cactusMisc_stringToName(flowerName));
    st_logInfo("Parsed the top level flower of the cactus tree to check\n");

    ///////////////////////////////////////////////////////////////////////////
    // Recursive check the flowers.
    ///////////////////////////////////////////////////////////////////////////

    int64_t startTime = time(NULL);
    FILE *fileHandle = fopen(outputFile, "w");
    //fprintf(fileHandle, "track name=%s\n", species);
    /*if(strstr(target, "reference") != NULL){
        //flower_addReferenceSequence(flower, cactusDisk, target);
        flower_addReferenceSequence(cactusDisk_getFlower(cactusDisk, cactusMisc_stringToName("0")), cactusDisk, target);
    }*/
    //addReferenceSequenceTest(flower);
    char *seq;
    struct List *seqs = getSequences(flower, query);    
    for(int i = 0; i < seqs->length; i++){//each sequence matches the query (each contig)
        seq = seqs->list[i];
        st_logInfo("Getting beds for sequence \"%s\"\n", seq);
        //getChains(flower, fileHandle, seq, 0);
        getCHAINs(flower, fileHandle, seq, target, 0);
    }
    fclose(fileHandle);
    st_logInfo("Got the beds in %i seconds/\n", time(NULL) - startTime);

    ///////////////////////////////////////////////////////////////////////////
    // Clean up.
    ///////////////////////////////////////////////////////////////////////////

    cactusDisk_destruct(cactusDisk);
    stKVDatabaseConf_destruct(kvDatabaseConf);

    return 0;
}
