/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

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

/*
 *Jul 4 2011: rewrite for better efficiency
 *Feb 7 2011: modify to include the referene sequence in the flower
 *Jan 19 2011: Correct to only put chain-segments in the same bed-record if there are links between them
 *(any two adjacent segments in the record must belong to two different blocks). 
 *Start a new bed record if hit an end with a self-edge
 *
 *Dec 09 2010: Modify code to adapt the Benedict's List data-structure
 *Sep 08 2010: nknguyen@soe.ucsc.edu
 *Print to outputFile bed record of each chain of the interested/inputed species.
 */
char *getCoor(char *header, int *start);
int segmentCmp(const void *s1,const void *s2);
//int segmentCmp(Segment *s1,Segment *s2);

struct Thread{
    char *header;
    char *chr;
    char *chainName;
    int32_t start;
    int32_t segmentNumber;
    stSortedSet *segments;
};

struct Thread *constructThread( char *header, char *chainName ){
    struct Thread *thread = st_malloc( sizeof(struct Thread) );
    thread->header = stString_copy( header );
    thread->chr = stString_copy(getCoor( header, &(thread->start) ));
    thread->chainName = stString_copy( chainName );
    thread->segments = stSortedSet_construct3(segmentCmp, NULL);
    return thread;
}

void destructThread( struct Thread *thread ){
    stSortedSet_destruct(thread->segments);
    free( thread->header );
    free( thread->chr );
    free( thread->chainName );
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

char *formatSequenceHeader(Sequence *sequence) {
    const char *sequenceHeader = sequence_getHeader(sequence);
    if (strlen(sequenceHeader) > 0) {
        char *cA = st_malloc(sizeof(char) * (1 + strlen(sequenceHeader)));
        sscanf(sequenceHeader, "%s", cA);
        return cA;
    } else {
        return cactusMisc_nameToString(sequence_getName(sequence));
    }
}

Sequence *getSequenceMatchesHeader(Flower *flower, char *header){
    //Returns the first Sequence whose name matches 'header'
    Flower_SequenceIterator *it = flower_getSequenceIterator(flower);
    Sequence *sequence = NULL;
    while((sequence = flower_getNextSequence(it)) != NULL){
        char *sequenceHeader = formatSequenceHeader(sequence);
        if(strstr(sequenceHeader, header) != NULL){
            free(sequenceHeader);
            break;
        }
        free(sequenceHeader);
    }
    flower_destructSequenceIterator(it);
    return sequence;
}

char *getCoor(char *header, int *start){
    char *chr;
    char *tok;
    *start = 0;
    char sep[] = ".";

    assert(header != NULL);
    strtok(stString_copy(header), sep); //species e.g "hg18"
    chr = strtok(NULL, sep);
    if(chr == NULL){
        chr = "";
    }else{
        tok = strtok(NULL, sep);//chromsize
        if(tok != NULL){
            sscanf(strtok(NULL, sep), "%d", start);
        }
    }
    return chr;
}

void block_getBED(Block *block, FILE *fileHandle, char *species, int level) {
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
            if(strcmp(eventHeader, species) == 0){
                int start;
                char *chr = getCoor( sequenceHeader, &start );
                if(!segment_getStrand(segment)){//if segment on the neg strand, reverse it
                    segment = segment_getReverse(segment);
                }
                Cap *cap5 = segment_get5Cap(segment);
                Cap *cap3 = segment_get3Cap(segment);
                //int s = cap_getCoordinate(cap5) + start - 2;
                //int e = cap_getCoordinate(cap3) + start + 1 - 2;
                int s = cap_getCoordinate(cap5) + start - 1;
                int e = cap_getCoordinate(cap3) + start + 1 - 1;
                fprintf(fileHandle, "%s %d %d %s.%d %d %s %d %d %d %d %d %d\n", chr, s, e, "NA", level, 0, ".", s, e, 0, 1, e-s,0);
            }
            free(sequenceHeader);
        }
    }
    block_destructInstanceIterator(instanceIterator);
}

void removeSubSortedSet( stSortedSet *segments, Segment *segment){
    st_logInfo("RemoveSubSortedSet...\n");
    Segment *sm;
    while( ( sm = stSortedSet_searchLessThan(segments, segment) ) != NULL ){
        stSortedSet_remove(segments, sm);
    }
    return;
}

void printThread( struct Thread *thread, FILE *fileHandle, int32_t level ){
    struct IntList *blockStarts = constructEmptyIntList(0);
    struct IntList *blockSizes = constructEmptyIntList(0);

    End *prevOtherEnd = NULL;

    stSortedSetIterator *it = stSortedSet_getIterator(thread->segments);
    Segment *sm;
    
    while( (sm = stSortedSet_getNext(it)) != NULL ){
        if( prevOtherEnd == NULL || isLinked( prevOtherEnd, cap_getEnd(segment_get5Cap(sm)) ) ){
            intListAppend(blockStarts, segment_getStart(sm));
            intListAppend(blockSizes, segment_getLength(sm));
        }else{
            break;
        }
        prevOtherEnd = cap_getEnd(segment_get3Cap(sm));
    }
    stSortedSet_destructIterator(it);

    if (sm != NULL){
        //Start new bed record. Remove all the previous Segments first:
        removeSubSortedSet(thread->segments, sm); 
        fprintf(fileHandle, "SPLIT TO A NEW BED RECORD!\n");
        printThread( thread, fileHandle, level );
    }

    int32_t blockCount = blockStarts->length;
    if( blockCount == 0 ){
        return;
    }
    //st_logInfo("\tblockCount = %d\n", blockCount);
    
    int32_t chromStart = blockStarts->list[0] + thread->start -1;
    int32_t chromEnd = blockStarts->list[blockCount -1] + blockSizes->list[blockCount -1] + thread->start -1;
    
    fprintf(fileHandle, "%s %d %d %s.%d %d %s %d %d %s %d ", thread->chr, chromStart, chromEnd,
                         thread->chainName, level, 0, ".", chromStart, chromEnd, "0", blockCount);
    
    //Print blockSizes
    int32_t i;
    for(i=0; i< blockCount; i++){
        fprintf(fileHandle, "%d,", blockSizes->list[i] );
    }
    destructIntList(blockSizes);
    fprintf(fileHandle, " ");
    for(i=0; i< blockCount; i++){
        fprintf(fileHandle, "%d,", blockStarts->list[i] );
    }
    destructIntList(blockStarts);
    fprintf(fileHandle, "\n");
    return;
}

void addSegmentToThread( struct Thread *thread, Segment *segment ){
    if( !segment_getStrand(segment) ){//make sure segment is on the + strand
        segment = segment_getReverse( segment );
    }
    stSortedSet_insert(thread->segments, segment);
    return;
}

void addSegments( struct List *threads, Block *block, char *header, char *chainName ){
    //st_logInfo("\taddSegments, header %s, chain %s\n", header, chainName);
    Segment *segment;
    int32_t startTime;
    struct Thread *thread = NULL;
    Block_InstanceIterator *it = block_getInstanceIterator(block);
    while( (segment = block_getNext(it)) != NULL ){
        Sequence *seq = segment_getSequence( segment );
        if(seq != NULL){
            char *seqHeader = formatSequenceHeader( seq );
            char *eventHeader = stString_copy( event_getHeader( sequence_getEvent(seq) ) );
            if( strcmp(eventHeader, header) == 0 ){
                //st_logInfo("\t\tFound %s\n", seqHeader);
                int32_t i;
                for(i = 0; i < threads->length; i++){
                    thread = threads->list[i];
                    if( strcmp(thread->header, seqHeader) == 0 ){
                        break;
                    }
                }
                if( i == threads->length ){
                    thread = constructThread( seqHeader, chainName );
                    listAppend(threads, thread);
                }
                startTime = time(NULL);
                addSegmentToThread( thread, segment );
                //st_logInfo("addSegmentToThread in %i seconds/\n", time(NULL) - startTime);
            }
            free(seqHeader);
            free(eventHeader);
        }
    }
    block_destructInstanceIterator( it );
}

void chain_getBEDs(Chain *chain, FILE *fileHandle, char *species, int level) {
    /*
     */
    //st_logInfo("chain_getBEDs: species %s, level %d\n", species, level);
    char *chainName = cactusMisc_nameToString(chain_getName(chain));
    //st_logInfo("\tchainName: %s\n", chainName);
    
    //Get all the threads with eventHeader 'species'
    int32_t startTime = time(NULL);
    struct List * threads = constructEmptyList(0, free);
    int32_t numBlocks;
    Block **blocks = chain_getBlockChain( chain, &numBlocks );
    for( int32_t i = 0; i< numBlocks; i++ ){
        Block * block = blocks[i];
        addSegments( threads, block, species, chainName );
    }
    free(blocks);
    //st_logInfo("Adding all segments for the chain in %i seconds/\n", time(NULL) - startTime);
    //st_logInfo("\t%d blocks; %d threads\n", numBlocks, threads->length);

    //Print the beds:
    startTime = time(NULL);
    for( int32_t i = 0; i< threads->length; i++ ){
        struct Thread *thread = threads->list[i];
        printThread( thread, fileHandle, level );
        destructThread(thread);
    }
    //st_logInfo("printThread in %i seconds/\n", time(NULL) - startTime);

    free(threads->list);
    free(threads);
    //destructList(threads);
    return;
}

void getBEDs(Flower *flower, FILE *fileHandle, char *species, int level){
    //st_logInfo("getBEDs, species %s\n", species);
    Chain *chain;
    int32_t startTime;
    Flower_ChainIterator *chainIt = flower_getChainIterator( flower );
    
    //Get beds for chains at current level
    while( (chain = flower_getNextChain(chainIt)) != NULL ){
        startTime = time(NULL);
        chain_getBEDs(chain, fileHandle, species, level);
        //st_logInfo("chain_getBEDs in %i seconds/\n", time(NULL) - startTime);
    }
    flower_destructChainIterator( chainIt );

    //Get beds for non-trivial chains:
    startTime = time(NULL);
    Flower_BlockIterator *blockIt = flower_getBlockIterator( flower );
    Block *block;
    while( (block = flower_getNextBlock(blockIt) ) != NULL ){
        if( block_getChain(block) == NULL ){//non-trivial chain
            block_getBED(block, fileHandle, species, level);
        }
    }
    flower_destructBlockIterator( blockIt );
    //st_logInfo("(block_getChain)s in %i seconds/\n", time(NULL) - startTime);

    //Call child flowers recursively
    Flower_GroupIterator *groupIt = flower_getGroupIterator( flower );
    Group *group;
    level ++;
    if (level > 5){ return; } // ONLY PRINT OUT THE FIRST 5 LEVEL BEDs
    while( (group = flower_getNextGroup(groupIt) ) != NULL ){
        Flower *nestedFlower = group_getNestedFlower( group );
        if( nestedFlower != NULL ){
            getBEDs(nestedFlower, fileHandle, species, level);
        }
    }
    flower_destructGroupIterator( groupIt );

}

void usage() {
    fprintf(stderr, "cactus_bedGenerator, version 0.2\n");
    fprintf(stderr, "Prints to output file all segments of the target sequence that are in blocks that contain both query & target\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr, "-b --species: species (sequence) of interest.If species = 'reference' then will build the reference sequence named 'reference'\n");
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
    char * species = NULL;

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    while(1) {
        static struct option long_options[] = {
            { "logLevel", required_argument, 0, 'a' },
            { "species", required_argument, 0, 'b' },
            { "cactusDisk", required_argument, 0, 'c' },
            { "flowerName", required_argument, 0, 'd' },
            { "outputFile", required_argument, 0, 'e' },
            { "help", no_argument, 0, 'h' },
            { 0, 0, 0, 0 }
        };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:b:c:d:e:h", long_options, &option_index);

        if(key == -1) {
            break;
        }

        switch(key) {
            case 'a':
                logLevelString = stString_copy(optarg);
                break;
            case 'b':
                species = stString_copy(optarg);
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
    assert(species != NULL);

    //////////////////////////////////////////////
    //Set up logging
    //////////////////////////////////////////////

    st_setLogLevelFromString(logLevelString);

    //////////////////////////////////////////////
    //Log (some of) the inputs
    //////////////////////////////////////////////

    st_logInfo("Flower name : %s\n", flowerName);
    st_logInfo("Output BED file : %s\n", outputFile);
    st_logInfo("Species: %s\n", species);

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
    
    if(strstr(species, "reference") != NULL && getSequenceMatchesHeader(flower, "reference") == NULL){
        fprintf(stderr, "No reference sequence found in cactusDisk\n");
        exit(EXIT_FAILURE);
    }
   
    getBEDs(flower, fileHandle, species, 0);
    fclose(fileHandle);
    st_logInfo("Got the beds in %i seconds/\n", time(NULL) - startTime);

    ///////////////////////////////////////////////////////////////////////////
    // Clean up.
    ///////////////////////////////////////////////////////////////////////////

    cactusDisk_destruct(cactusDisk);
    stKVDatabaseConf_destruct(kvDatabaseConf);

    return 0;
}
