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
#include "sonLibSortedSet.h"
//#include "cactus_addReferenceSeq.h"

/*
 *Feb 7 2011: modify to include the referene sequence in the flower
 *Jan 19 2011: Correct to only put chain-segments in the same bed-record if there are links between them
 *(any two adjacent segments in the record must belong to two different blocks). 
 *Start a new bed record if hit an end with a self-edge
 *
 */

char *getCoor(char *header, int32_t  *start, int32_t *chrsize);
int segmentCmp(const void *s1,const void *s2);
//int segmentCmp(Segment *s1,Segment *s2);

struct Thread{
    char *header;
    char *chr;
    int32_t start;
    stSortedSet *segments;
};

struct Thread *constructThread( char *header ){
    struct Thread *thread = st_malloc( sizeof(struct Thread) );
    thread->header = stString_copy( header );
    int32_t chrsize = 0;
    thread->chr = stString_copy(getCoor( header, &(thread->start), &chrsize ));
    thread->segments = stSortedSet_construct3(segmentCmp, NULL);
    return thread;
}

void destructThread( struct Thread *thread ){
    stSortedSet_destruct(thread->segments);
    free( thread->header );
    free( thread->chr );
    free( thread);
}

int segmentCmp(const void *sm1, const void *sm2){
    //make sure segment is on the + strand
    /*if( !segment_getStrand(s1) ){ s1 = segment_getReverse( s1 ); }
    if( !segment_getStrand(s2) ){ s2 = segment_getReverse( s2 ); }*/
    Segment *s1 = (Segment *) sm1;
    Segment *s2 = (Segment *) sm2;
    assert( segment_getStrand(s1) );
    assert( segment_getStrand(s2) );
    int32_t start1 = segment_getStart(s1);
    int32_t start2 = segment_getStart(s2);
    if(start1 == start2) {
        return 0;
    }else if( start1 < start2){
        return -1;
    }else{
        return 1;
    }
}

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

Segment *block_getSegmentByEventHeader(Block *block, char *header){
    Segment *segment;
    Block_InstanceIterator *it = block_getInstanceIterator(block);
    while( (segment = block_getNext(it)) != NULL ){
        Sequence *seq = segment_getSequence( segment );
        if(seq != NULL){
            char *eventHeader = stString_copy( event_getHeader( sequence_getEvent(seq) ) );
            if( strcmp(eventHeader, header) == 0 ){
                break;
            }
            free(eventHeader);
        }
    }
    block_destructInstanceIterator( it );
    return segment;
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

char *getCoor(char *header, int32_t  *start, int32_t *chrsize){
    char *chr;
    char *tok;
    *start = 0;
    char sep[] = ".";

    assert(header != NULL);
    strtok( stString_copy(header), sep );//species e.g "hg18"
    chr = strtok(NULL, sep);
    if(chr == NULL){
        chr = "";
    }else{
        tok = strtok(NULL, sep);
        if(tok != NULL){
            sscanf( tok, "%d", chrsize );//chromsize
            sscanf( strtok(NULL,sep), "%d", start);
        }
    }
    return chr;
}

int block_getCHAIN(Block *block, FILE *fileHandle, char *query, char *target, int32_t chainid) {
    /*
     */
    Block_InstanceIterator *instanceIterator = block_getInstanceIterator(block);
    Segment *segment;
    while((segment = block_getNext(instanceIterator)) != NULL) {
        Sequence *sequence = segment_getSequence(segment);
        if(sequence != NULL) {
            char *sequenceHeader = formatSequenceHeader(sequence);
            Event *event = sequence_getEvent( sequence );
            char *eventHeader = stString_copy( event_getHeader( event ) ); 
            if(strcmp(eventHeader, query) == 0){
                Cap *tcap5 = end_getCapBySeqName(cap_getEnd(segment_get5Cap(segment)), target); 
                if(tcap5 == NULL){continue;}
                if(!cap_getStrand(tcap5)){//if reference segment is on the neg strand, reverse it
                    tcap5 = cap_getOtherSegmentCap( cap_getReverse(tcap5) );
                    segment = segment_getReverse(segment);
                }
                Cap *cap5 = segment_get5Cap(segment);
                Cap *cap3 = segment_get3Cap(segment);
                int32_t qstart = cap_getCoor(cap5);
                int32_t qend = cap_getCoor(cap3) + 1;
                int32_t tstart = cap_getCoor(tcap5);
                int32_t tend = cap_getCoor(cap_getOtherSegmentCap(tcap5)) + 1;
                int32_t tsize = cap_getSeqSize(tcap5);
                int32_t qsize = cap_getSeqSize(cap5);
                char *tchr = cap_getChr(tcap5);
                char tstrand = cap_getStrand(tcap5) ? '+' : '-';
                char qstrand = cap_getStrand(cap5) ? '+' : '-';
               
                int32_t qchrsize = qsize; 
                int32_t start = 0;
                //char *qchr = getCoor( sequenceHeader, &start, &qchrsize );
                char *qchr;
                if( strstr(sequenceHeader, "NODE") != NULL ){//HACK for Velvet's contig headers
                    strtok(stString_copy(sequenceHeader), "_");
                    qchr = strtok(NULL, "_");
                }else{
                    qchr = getCoor( sequenceHeader, &start, &qchrsize );
                }

                qstart = convertCoor(qstart, qchrsize, start, qsize, qstrand);
                qend = convertCoor(qend, qchrsize, start, qsize, qstrand);
                fprintf(fileHandle, "chain %d %s %d %c %d %d %s %d %c %d %d %d\n", 0, tchr, tsize, tstrand, tstart, tend, qchr, qchrsize, qstrand, qstart, qend, chainid);
                chainid ++;
                fprintf(fileHandle, "%d\n\n", segment_getLength(segment));
            }
            free(sequenceHeader);
        }
    }
    block_destructInstanceIterator(instanceIterator);
    return chainid;
}

void removeSubSortedSet( stSortedSet *segments, Segment *segment){
    st_logInfo("RemoveSubSortedSet...\n");
    Segment *sm;
    while( ( sm = stSortedSet_searchLessThan(segments, segment) ) != NULL ){
        stSortedSet_remove(segments, sm);
    }
    return;
}

int32_t printThread( struct Thread *thread, FILE *fileHandle, char *query, char *target, int32_t chainid ){
    char qstrand = '.';
    char tstrand = '.';
    char *tchr = "";
    char *qchr = "";
    int32_t start = 0;
    int32_t qchrsize;
    int32_t qstart = 0;
    int32_t tstart = 0;
    int32_t qend = -1;
    int32_t tend = -1;
    int32_t qsize = 0;
    int32_t tsize = 0;
    int32_t dq = 0;
    int32_t dt = 0;
    struct IntList *blockSizes = constructEmptyIntList(0);
    struct IntList *dtList = constructEmptyIntList(0);
    struct IntList *dqList = constructEmptyIntList(0);
    End *prevOtherEnd = NULL;

    stSortedSetIterator *it = stSortedSet_getIterator(thread->segments);
    Segment *sm;
    Cap *qcap;
    Cap *tcap;

    while( (sm = stSortedSet_getNext(it)) != NULL ){
        qcap = segment_get5Cap(sm);
        Cap *qothercap = segment_get3Cap(sm);
        tcap = end_getCapBySeqName( cap_getEnd(qcap), target ); 
        if( prevOtherEnd == NULL || isLinked( prevOtherEnd, cap_getEnd(qcap) ) ){
            intListAppend(blockSizes, segment_getLength(sm));
            if(qend == -1){
                qstrand = cap_getStrand(qcap) ? '+' : '-';
                tstrand = cap_getStrand(tcap) ? '+' : '-';
                qstart = cap_getCoor(qcap);
                tstart = cap_getCoor(tcap);
                qsize = cap_getSeqSize(qcap);
                tsize = cap_getSeqSize(tcap);
                tchr = cap_getChr(tcap);
            }else{
                assert( (qstrand == '+' && cap_getStrand(qcap)) || (qstrand == '-' && !cap_getStrand(qcap)) );
                assert( (tstrand == '+' && cap_getStrand(tcap)) || (tstrand == '-' && !cap_getStrand(tcap)) );
                dt = abs(cap_getCoor(tcap) - tend);
                dq = abs(cap_getCoor(qcap) - qend);
                intListAppend(dtList, dt);
                intListAppend(dqList, dq);
            }
            qend = cap_getCoor(qothercap) + 1;
            tend = cap_getCoor(cap_getOtherSegmentCap(tcap)) + 1;
        }else{
            break;
        }
        prevOtherEnd = cap_getEnd(qothercap);
    }
    stSortedSet_destructIterator(it);
    
    //PRINTING:
    if(blockSizes->length == 0){
       //fprintf(fileHandle, "No block found for chain %s\n", chainName);
        return chainid;
    }
            
    //get query chromosome
    qchrsize = qsize; 
    if( strstr(thread->header, "NODE") != NULL ){//HACK for Velvet's contig headers
        strtok(stString_copy(thread->header), "_");
        qchr = strtok(NULL, "_");
    }else if( strstr(thread->header, "Gi") != NULL &&  strstr(thread->header, "Ref") != NULL ){ // HACK ecoli 
        strtok(stString_copy(thread->header), "_");
        qchr = strtok(NULL, ".");
    }else{
        qchr = getCoor( thread->header, &start, &qchrsize );
    }

    //Convert query coordinates to assembly coordinates:
    qstart = convertCoor(qstart, qchrsize, start, qsize, qstrand);
    qend = convertCoor(qend, qchrsize, start, qsize, qstrand);

    int score = 0;
    
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

    destructIntList(blockSizes);
    destructIntList(dqList);
    destructIntList(dtList);
    
    //Picks up where the chain records split, continue to print the rest of the chain:
    if (sm != NULL){
        //Start new bed record. Remove all the previous Segments first:
        removeSubSortedSet(thread->segments, sm);
        chainid = printThread( thread, fileHandle, query, target, chainid );
    }

    return chainid;
}

void addSegmentToThread( struct Thread *thread, Segment *segment ){
    if( !segment_getStrand(segment) ){//make sure segment is on the + strand
        segment = segment_getReverse( segment );
    }
    stSortedSet_insert(thread->segments, segment);
    return;
}

void addSegments(struct List *threads, Block *block, char *query, char *target){
    Segment *tsegment = block_getSegmentByEventHeader(block, target);
    if( tsegment == NULL ){ return; }

    Segment *segment;
    struct Thread *thread = NULL;
    Block_InstanceIterator *it = block_getInstanceIterator(block);
    while( (segment = block_getNext(it)) != NULL ){
        Sequence *seq = segment_getSequence( segment );
        if(seq != NULL){
            char *seqHeader = formatSequenceHeader( seq );
            char *eventHeader = stString_copy( event_getHeader( sequence_getEvent(seq) ) );
            if( strcmp(eventHeader, query) == 0 ){
                int32_t i;
                for(i = 0; i < threads->length; i++){
                    thread = threads->list[i];
                    if( strcmp(thread->header, seqHeader) == 0 ){
                        break;
                    }
                }
                if( i == threads->length ){
                    thread = constructThread( seqHeader );
                    listAppend(threads, thread);
                }
                addSegmentToThread( thread, segment );
            }
            free(seqHeader);
            free(eventHeader);
        }
    }
    block_destructInstanceIterator( it );
}

int32_t chain_getCHAINs(Chain *chain, FILE *fileHandle, char *query, char *target, int32_t chainid) {
    struct List *threads = constructEmptyList(0, free);
    int32_t numBlocks;
    Block **blocks = chain_getBlockChain( chain, &numBlocks); 
    for(int32_t i=0; i< numBlocks; i++){
        addSegments( threads, blocks[i], query, target );
    }
    free(blocks);

    //Print the chains:
    for( int32_t i = 0; i< threads->length; i++ ){
        struct Thread *thread = threads->list[i];
        chainid = printThread( thread, fileHandle, query, target, chainid );
        destructThread(thread);
    }

    free(threads->list);
    free(threads);
    //destructList(threads);
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

int getCHAINs(Flower *flower, FILE *fileHandle, char *query, char *target, int32_t chainid){
    Chain *chain;
    int32_t startTime;
    Flower_ChainIterator *chainIt = flower_getChainIterator( flower );

    //Get chains of current level
    while( (chain = flower_getNextChain(chainIt) ) != NULL ){
        startTime = time(NULL);
        chainid = chain_getCHAINs( chain, fileHandle, query, target, chainid );
        st_logInfo("chain_getCHAINs in %i seconds/\n", time(NULL) - startTime);
    }
    flower_destructChainIterator( chainIt );

    //Get beds for non-trivial chains
    startTime = time(NULL);
    Flower_BlockIterator *blockIterator = flower_getBlockIterator(flower);
    Block *block;
    while((block = flower_getNextBlock(blockIterator)) != NULL){
        if(block_getChain(block) == NULL){//non-trivial chain
            chainid = block_getCHAIN(block, fileHandle, query, target, chainid);
        }
    }
    flower_destructBlockIterator(blockIterator);
    st_logInfo("(block_getChain)s in %i s/\n", time(NULL) - startTime );

    //Call child flowers recursively.
    Flower_GroupIterator *groupIterator = flower_getGroupIterator(flower);
    Group *group;
    while((group = flower_getNextGroup(groupIterator)) != NULL) {
        Flower *nestedFlower = group_getNestedFlower(group);
        if(nestedFlower != NULL) {
            chainid = getCHAINs(group_getNestedFlower(group), fileHandle, query, target, chainid); //recursive call.
        }
    }
    flower_destructGroupIterator(groupIterator);
    return chainid;
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
    st_setLogLevelFromString(logLevelString);

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
    getCHAINs(flower, fileHandle, query, target, 0);
    fclose(fileHandle);
    st_logInfo("Got the chains in %i seconds/\n", time(NULL) - startTime);

    ///////////////////////////////////////////////////////////////////////////
    // Clean up.
    ///////////////////////////////////////////////////////////////////////////

    cactusDisk_destruct(cactusDisk);
    stKVDatabaseConf_destruct(kvDatabaseConf);

    return 0;
}
