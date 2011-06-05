/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>

#include "cactus.h"

/*
 * Sep 27 2010: nknguyen@soe.ucsc.edu (based on the reference-parts in Benedict's cactus_MAFGenerator.c)
 * Adding the reference sequence into cactus structure
 */

typedef struct _refseq {
    int32_t length;
    char *header;
    int32_t index;
    char *string;
} ReferenceSequence;

End *getPseudoAdjacentEnd(End *end) {
    assert(end != NULL);
    PseudoAdjacency *pseudoAdjacency = end_getPseudoAdjacency(end);
    End *adjend;
    assert(pseudoAdjacency != NULL);

    assert(
            pseudoAdjacency_get3End(pseudoAdjacency) == end
                    || pseudoAdjacency_get5End(pseudoAdjacency) == end);
    //assert( pseudoAdjacency_get3End(pseudoAdjacency) != pseudoAdjacency_get5End(pseudoAdjacency) );
    if (pseudoAdjacency_get3End(pseudoAdjacency) == end) {
        adjend = pseudoAdjacency_get5End(pseudoAdjacency);
    } else {
        adjend = pseudoAdjacency_get3End(pseudoAdjacency);
    }
    return adjend;
}

static int32_t pseudoChromosome_getLength(End *end) {
    /*
     *Return the total number of bases of the pseudochromosome (and of all 
     *pseudochromosomes at lower flowers within this pseudochrom) that 'end' belongs to, 
     */
    int32_t len = 0;
    assert(end_isStubEnd(end));
    Group *group;
    End *inheritedEnd;

    group = end_getGroup(end);
    if (!group_isLeaf(group)) {//has lower level
        inheritedEnd = flower_getEnd(group_getNestedFlower(group),
                end_getName(end));
        len += pseudoChromosome_getLength(inheritedEnd);
    }

    end = getPseudoAdjacentEnd(end);
    while (end_isBlockEnd(end)) {
        Block *block = end_getBlock(end);
        //if (block_getInstanceNumber(block) > 0) {
        len += block_getLength(block);
        //}

        end = end_getOtherBlockEnd(end);
        group = end_getGroup(end);
        if (!group_isLeaf(group)) {//has lower level
            inheritedEnd = flower_getEnd(group_getNestedFlower(group),
                    end_getName(end));
            len += pseudoChromosome_getLength(inheritedEnd);
        }
        end = getPseudoAdjacentEnd(end);
    }
    return len;
}

static ReferenceSequence *referenceSequence_construct(Flower *flower,
        char *header, int length) {
    ReferenceSequence *refseq = st_malloc(sizeof(ReferenceSequence));
    refseq->index = 0;
    refseq->header = stString_copy(header);
    refseq->length = length;
    refseq->string = st_malloc(sizeof(char) * (length + 1));
    strcpy(refseq->string, "");
    return refseq;
}

static void referenceSequence_destruct(ReferenceSequence *refseq) {
    free(refseq->string);
    free(refseq->header);
    free(refseq);
}

char *formatSequenceHeader1(Sequence *sequence) {
    const char *sequenceHeader = sequence_getHeader(sequence);
    if (strlen(sequenceHeader) > 0) {
        char *cA = st_malloc(sizeof(char) * (1 + strlen(sequenceHeader)));
        sscanf(sequenceHeader, "%s", cA);
        return cA;
    } else {
        return cactusMisc_nameToString(sequence_getName(sequence));
    }
}

char *getConsensusString(Block *block) {
    //st_logInfo("\nGetting consensus string, blockOrientation: %d...:\n", block_getOrientation(block));
    Block_InstanceIterator *it = block_getInstanceIterator(block);
    Segment *segment;
    struct List *list = constructEmptyList(0, free);
    while ((segment = block_getNext(it)) != NULL) {
        if (segment_getSequence(segment) != NULL) {
            listAppend(list, segment_getString(segment));
            //st_logInfo("%s\t%d\t%s\n", event_getHeader(segment_getEvent(segment)), segment_getStrand(segment), segment_getString(segment));
        }
    }
    block_destructInstanceIterator(it);
    assert(list->length > 0);
    if(list->length == 1){
        char *cA = stString_copy(list->list[0]);
        destructList(list);
        return cA;
    }
    //st_logInfo("Number of segments: %d\n", list->length);

    char *consensusStr = stString_copy(list->list[0]);
    //st_logInfo("Initialize consensusStr: %s\n", consensusStr);
    char bases[] = {'A', 'C', 'G', 'T', 'N'};
    int numchar = 5;
    for(int32_t i=0; i < strlen(consensusStr); i++){//each base position
        int32_t freq[] = {0, 0, 0, 0, 0};
        for(int32_t j=0; j < list->length; j++){//each segment in the block
            char *str = list->list[j];
            char base = toupper( *(str + i) );
            int k;
            for (k = 0; k < numchar; k++){
                if( base == bases[k] ){
                    freq[k] += 1;
                    break;
                }
            }
            assert(k < numchar);
        }
        int32_t max = 0;
        //st_logInfo("Pos %d, base freq: ", i);
        for(int k = 0; k < numchar; k++){
            //st_logInfo("%c: %d, ", bases[k], freq[k]);
            if (max < freq[k]){
                max = freq[k];
                *(consensusStr + i) = bases[k];
            }
        }
        //st_logInfo("\n");
    }
    //st_logInfo("ConsensusStr: %s\n\n", consensusStr);
    destructList(list);
    return consensusStr;
}

Sequence *getSequenceByHeader(Flower *flower, char *header) {
    Flower_SequenceIterator *it = flower_getSequenceIterator(flower);
    Sequence *sequence;
    while ((sequence = flower_getNextSequence(it)) != NULL) {
        char *sequenceHeader = formatSequenceHeader1(sequence);
        if (strcmp(sequenceHeader, header) == 0) {
            flower_destructSequenceIterator(it);
            free(sequenceHeader);
            return sequence;
        }
        free(sequenceHeader);
    }
    flower_destructSequenceIterator(it);
    return NULL;
}

/*
 *Return: the first cap in input 'end' that belongs to sequence 'header'
 *        NULL if not found
 */
Cap *end_getCapByHeader(End *end, char *header) {
    End_InstanceIterator *it = end_getInstanceIterator(end);
    Cap *cap;
    while ((cap = end_getNext(it)) != NULL) {
        Sequence *sequence = cap_getSequence(cap);
        assert(cap_getSide(cap) == end_getSide(end));
        if (sequence != NULL) {
            char *sequenceHeader = formatSequenceHeader1(sequence);
            if (strcmp(sequenceHeader, header) == 0) {
                end_destructInstanceIterator(it);
                free(sequenceHeader);
                return cap;
            }
            free(sequenceHeader);
        }
    }
    end_destructInstanceIterator(it);
    return NULL;
}

void addAdj(End *end, End *adjEnd, char *header) {
    /*
     *Add adjacency between caps of 'header' sequence in 'end' and 'adjEnd'
     */
    Cap *cap = end_getCapByHeader(end, header);
    Cap *cap2 = end_getCapByHeader(adjEnd, header);
    assert(cap != NULL);
    assert(cap2 != NULL);
    cap_makeAdjacent(cap, cap2);
    //cap_check(cap);
    //cap_check(cap2);
}

Segment *addReferenceSegmentToBlock(End *end, ReferenceSequence *refseq) {
    /*
     */
    End *pseudoAdjEnd = getPseudoAdjacentEnd(end);
    Cap *adjcap = end_getCapByHeader(pseudoAdjEnd, refseq->header);
    assert(adjcap);
    bool strand = end_getSide(pseudoAdjEnd) != end_getSide(end) ? cap_getStrand(adjcap) : cap_getStrand(cap_getReverse(adjcap));
    Block *block = end_getBlock(end);

    Sequence *sequence = getSequenceByHeader(block_getFlower(block),
            refseq->header);
    assert(sequence != NULL);

    //Adding segment to block
    Segment *segment = segment_construct2(block, refseq->index, strand, sequence);
    refseq->index += block_getLength(block);

    /*//DEBUG
    Block_InstanceIterator *it = block_getInstanceIterator(block);
    Segment *sm;
    st_logInfo("\nBLOCK (segment added with strand %d)\n", strand);
    while ((sm = block_getNext(it)) != NULL) {
        if (segment_getSequence(segment) != NULL) {
            st_logInfo("%s\t%d\t%s\n", event_getHeader(segment_getEvent(sm)), segment_getStrand(sm), segment_getString(sm));
        }
    }
    block_destructInstanceIterator(it);
    //END DEBUG 
    */

    //Adding adj
    addAdj(end, pseudoAdjEnd, refseq->header);
    
    return segment;
}

bool block_metaSequence(End *end, ReferenceSequence *refseq, bool prevStrand) {
    /*
     *Adding sequence of end's block to MetaSequence
     */
    End *pseudoAdjEnd = getPseudoAdjacentEnd(end);
    bool strand = end_getSide(end) != end_getSide(pseudoAdjEnd) ? prevStrand : (!prevStrand) ;

    /*fprintf(stdout, "PrevOrientation: %d, currOrientation: %d\n", end_getOrientation(pseudoAdjEnd), end_getOrientation(end));
    
    if(end_isBlockEnd(pseudoAdjEnd) && end_isBlockEnd(end)){
    fprintf(stdout, "PrevOtherEnd: %s - PrevEnd: %s; CurrEnd: %s - currOtherEnd: %s\n", cactusMisc_nameToString(end_getName(end_getOtherBlockEnd(pseudoAdjEnd))), 
                      cactusMisc_nameToString(end_getName(pseudoAdjEnd)), cactusMisc_nameToString(end_getName(end)), cactusMisc_nameToString(end_getName(end_getOtherBlockEnd(end))) );
    }
    fprintf(stdout, "PrevSide: %d, currSide: %d; PrevStrand %d, CurrStrand: %d\n", end_getSide(pseudoAdjEnd), end_getSide(end),  prevStrand, strand);
    */
    Block *block = strand ? end_getBlock(end) : end_getBlock(end_getReverse(end));
    //Block *block = end_getSide(end) != end_getSide(pseudoAdjEnd) ? end_getBlock(end) : end_getBlock(end_getReverse(end));
    //Block *block = end_getBlock(end);
    if (block_getInstanceNumber(block) > 0) {
        char *instanceString = getConsensusString(block);
        refseq->string = strcat(refseq->string, instanceString);
        free(instanceString);
    } else {//if block doesn't have any instance, add 'N'
        assert(block_getLength(block) == 1);
        refseq->string = strcat(refseq->string, "N");
    }
    return strand;
}

Cap *copyRefCapToLowerFlowers(Cap *cap) {
    assert(cap != NULL);
    End *end = cap_getEnd(cap);
    Group *group = end_getGroup(end);
    Flower *nestedflower = group_getNestedFlower(group);
    if (nestedflower != NULL) { //has lower level
        End *inheritedEnd = flower_getEnd(nestedflower, end_getName(end));
        if (end_getSide(end) != end_getSide(inheritedEnd)) {//make sure end & inheritedEnd are the same
            inheritedEnd = end_getReverse(inheritedEnd);
        }
        cap = cap_copyConstruct(inheritedEnd, cap);
        copyRefCapToLowerFlowers(cap);
    }
    return cap;
}

void addStubAdjacency(End *end, char *header){
    End *pseudoAdjEnd = getPseudoAdjacentEnd(end);
    addAdj(end, pseudoAdjEnd, header);
    Group *group = end_getGroup(end);
    Flower *nestedflower = group_getNestedFlower(group);
    if(nestedflower != NULL){
        End *inheritedEnd = flower_getEnd(nestedflower, end_getName(end));
        addStubAdjacency(inheritedEnd, header);
    }
    return;
}

void reference_walkDown(End *end, ReferenceSequence *refseq, struct IntList *prevStrands);

void reference_walkUp(End *end, ReferenceSequence *refseq, struct IntList *prevStrands) {
    assert(end != NULL);
    if (end_isBlockEnd(end)) {
        if (strlen(refseq->string) == refseq->length) {
            //if(block_getInstanceNumber(block) > 0){
            Segment *segment = addReferenceSegmentToBlock(end, refseq);
            segment_check(segment);
            copyRefCapToLowerFlowers(segment_get5Cap(segment));
            copyRefCapToLowerFlowers(segment_get3Cap(segment));
            //}
        } else {
            bool prevStrand = prevStrands->list[prevStrands->length -1];
            prevStrands->list[prevStrands->length - 1] = block_metaSequence(end, refseq, prevStrand);
        }
        reference_walkDown(end_getOtherBlockEnd(end), refseq, prevStrands);
    } else {
        assert(end_isAttached(end));
        Group *parentGroup = flower_getParentGroup(end_getFlower(end));
        if (parentGroup != NULL) {
            End *parentEnd = group_getEnd(parentGroup, end_getName(end));
            if (end_getSide(parentEnd) != end_getSide(end)){
                parentEnd = end_getReverse(parentEnd);
            }
            if( prevStrands ){
                prevStrands->length -= 1;
                //listRemove(prevStrands, prevStrands->list[prevStrands->length - 1]);
            }
            reference_walkUp(parentEnd, refseq, prevStrands);
            if(refseq->index > 0){
                addAdj(end, getPseudoAdjacentEnd(end), refseq->header);
            }
        } else { //We reached the end of a pseudo-chromosome!
            assert(
                    pseudoChromosome_get3End(
                            pseudoAdjacency_getPseudoChromosome(
                                    end_getPseudoAdjacency(end))) == end);
            //adding last Stub (5')
            if (refseq->index > 0) {
                Sequence *sequence = getSequenceByHeader(end_getFlower(end),
                        refseq->header);
                End *pseudoAdjEnd = getPseudoAdjacentEnd(end);
                if( end_getSide(pseudoAdjEnd) == end_getSide(end) ){
                    end = end_getReverse(end);
                }
                Cap *endCap = cap_construct2(end, refseq->index, true, sequence);
                copyRefCapToLowerFlowers(endCap);

                addStubAdjacency(end, refseq->header);
            }
        }
    }
}

void reference_walkDown(End *end, ReferenceSequence *refseq, struct IntList *prevStrands) {
    assert(end != NULL);
    //assert(end_isAttached(end));
    Group *group = end_getGroup(end);
    if (group_isLeaf(group)) { //Walk across
        end = getPseudoAdjacentEnd(end);
        //Now walk up
        reference_walkUp(end, refseq, prevStrands);
    } else { //Walk down
        End *inheritedEnd = flower_getEnd(group_getNestedFlower(group), end_getName(end));
        if (end_getSide(end) != end_getSide(inheritedEnd)){
            inheritedEnd = end_getReverse(inheritedEnd);
        }

        if(prevStrands){
            bool strand = prevStrands->list[prevStrands->length -1];
            intListAppend(prevStrands, strand);
        }
        reference_walkDown(inheritedEnd, refseq, prevStrands);
    }
}

MetaSequence *constructReferenceMetaSequence(End *end, CactusDisk *cactusDisk,
        ReferenceSequence *refseq, Event *event) {
    /*
     *Traverse pseudochromosome (of 'end') and its lower levels to get the reference MetaSequence
     */
    st_logInfo("Getting reference MetaSequence...\n");
    Name eventName = event_getName(event);
    MetaSequence *metaSequence;
    int32_t start = 1;

    struct IntList * strands = constructEmptyIntList(0);
    bool firststrand = true;
    intListAppend(strands, firststrand);
    reference_walkDown(end, refseq, strands);
    destructIntList(strands);
    
    assert(strlen(refseq->string) == refseq->length);
    metaSequence = metaSequence_construct(start, strlen(refseq->string),
            refseq->string, refseq->header, eventName, cactusDisk);
    return metaSequence;
}

void constructReferenceSequence(MetaSequence *metaSequence, Flower *flower, Name name) {
    /*
     *Attach MetaSequence to cactus flowers
     */
    assert(flower != NULL);
    Reference *ref = flower_getReference(flower);

    //add reference sequence to current flower
    PseudoChromosome *pc = reference_getPseudoChromosome(ref, name);
    if( pc != NULL ){
        sequence_construct(metaSequence, flower);
    }

    //recursively add reference sequence to lower-level flowers
    Flower_GroupIterator *it = flower_getGroupIterator(flower);
    Group *group;
    while ((group = flower_getNextGroup(it)) != NULL) {
        Flower *nestedFlower = group_getNestedFlower(group);
        if (nestedFlower != NULL) {
            constructReferenceSequence(metaSequence, nestedFlower, name);
        }
    }
    flower_destructGroupIterator(it);
}

char *getChromName(char *name, int num) {
    /*
     *Add the order of the pseudochromosome to the reference sequence name. E.g 'reference.chr1'
     */
    char chrom[5];//assume there are less than 999 chroms! 
    sprintf(chrom, "%d", num);
    char *chromName = st_malloc(
            sizeof(char) * (strlen(name) + strlen(chrom) + strlen(".chr") + 1));
    strcpy(chromName, name);
    strcat(chromName, ".chr");
    strcat(chromName, chrom);
    return chromName;
}

Event *eventTree_getEventByHeader(EventTree *eventTree, char *header){
    EventTree_Iterator *it = eventTree_getIterator(eventTree);
    Event *event;
    while((event = eventTree_getNext(it)) != NULL){
        if( strcmp(event_getHeader(event), header) == 0 ){
            eventTree_destructIterator(it);
            return event;
        }
    }
    eventTree_destructIterator(it);
    return NULL;
}

void addReferenceEvent2(Flower *flower, char *header, Name eventName){
    EventTree *eventTree = flower_getEventTree(flower);
    Event *rootEvent = eventTree_getRootEvent(eventTree);
    assert(eventName);
    Event *event = event_construct2(eventName, header, INT32_MAX, rootEvent, event_getChild(rootEvent, 0), eventTree);

    //Recursively addRefernceEvent for lower flowers
    Flower_GroupIterator *it = flower_getGroupIterator(flower);
    Group *group;
    while( (group = flower_getNextGroup(it)) != NULL ){
        Flower *nestedFlower = group_getNestedFlower(group);
        if(nestedFlower){
            addReferenceEvent2(nestedFlower, header, event_getName(event));
        }
    }
    flower_destructGroupIterator(it);
    return;
}

Event *addReferenceEvent(Flower *flower, char *header){
    EventTree *eventTree = flower_getEventTree(flower);
    Event *rootEvent = eventTree_getRootEvent(eventTree);
    Event *event = event_construct4(header, INT32_MAX, rootEvent, event_getChild(rootEvent, 0), eventTree);
    
    //Recursively addRefernceEvent for lower flowers
    Flower_GroupIterator *it = flower_getGroupIterator(flower);
    Group *group;
    while( (group = flower_getNextGroup(it)) != NULL ){
        Flower *nestedFlower = group_getNestedFlower(group);
        if(nestedFlower){
            addReferenceEvent2(nestedFlower, header, event_getName(event));
        }
    }
    flower_destructGroupIterator(it);
    return event;
}

/*
 *Adding the reference sequence of each pseudoChrom to cactus structure.
 *'header' is going to be set as the reference's name
 */
Flower *flower_addReferenceSequence(Flower *flower, CactusDisk *cactusDisk,
        char *header) {
    int64_t startTime = time(NULL);
    //Return if event with 'header' has already existed. Otherwise, add event 'header'
    EventTree *eventTree = flower_getEventTree(flower);
    st_logInfo("flower_getEventTree:\t%d seconds\n", time(NULL) - startTime);

    Event *event;

    startTime = time(NULL);
    event = eventTree_getEventByHeader(eventTree, header);
    st_logInfo("eventTree_getEventByHeader:\t%d seconds\n", time(NULL) - startTime);

    if(event) {
        return flower;
    }else{
        startTime = time(NULL);
        event = addReferenceEvent(flower, header);
        st_logInfo("addReferenceEvent:\t%d seconds\n", time(NULL) - startTime);
        //event = eventTree_getEventByHeader(eventTree, header);
        assert(event != NULL);
    }

    startTime = time(NULL);
    Reference *reference = flower_getReference(flower);
    st_logInfo("flower_getReference:\t%d seconds\n", time(NULL) - startTime);
    
    assert(reference != NULL);
    Reference_PseudoChromosomeIterator *it =
            reference_getPseudoChromosomeIterator(reference);
    PseudoChromosome *pseudoChromosome;
    int chromNum = 0;

    while ((pseudoChromosome = reference_getNextPseudoChromosome(it)) != NULL) {//Each pseudoChrom
        chromNum++;
        char *chromHeader = getChromName(header, chromNum);
        End *end = pseudoChromosome_get5End(pseudoChromosome);
        assert(!end_isBlockEnd(end));

        int len = pseudoChromosome_getLength(end);
        if(len == 0){continue;}
        //assert(len != 0);
        ReferenceSequence *refseq = referenceSequence_construct(flower,
                chromHeader, len);

        st_logInfo(
                "\nInitialize refseq: index %d, length %d, header *%s*, string *%s*\n",
                refseq->index, refseq->length, refseq->header, refseq->string);

        //Construct the MetaSequence 
        startTime = time(NULL);
        MetaSequence *metaSequence = constructReferenceMetaSequence(end,
                cactusDisk, refseq, event);
        st_logInfo("constructReferenceMetaSequence:\t%d seconds\n", time(NULL) - startTime);
        /*st_logInfo(
                "Got metasequence: name *%s*, start %d, length %d, header *%s*, event %s\n",
                cactusMisc_nameToString(metaSequence_getName(metaSequence)),
                metaSequence_getStart(metaSequence),
                metaSequence_getLength(metaSequence),
                metaSequence_getHeader(metaSequence),
                event_getHeader(eventTree_getEvent(eventTree, metaSequence_getEventName(metaSequence))));
        */

        //st_logInfo("\nConstructing reference sequence...\n");
        startTime = time(NULL);
        constructReferenceSequence(metaSequence, flower, pseudoChromosome_getName(pseudoChromosome) );
        st_logInfo("constructReferenceSequence:\t%d seconds\n", time(NULL) - startTime);
        //st_logInfo("Constructed reference sequence successfully.\n");

        //Add startStub (3' end)
        startTime = time(NULL);
        Sequence *sequence = getSequenceByHeader(flower, refseq->header);
        st_logInfo("getSequenceByHeader:\t%d seconds\n", time(NULL) - startTime);

        startTime = time(NULL);
        Cap *startcap = cap_construct2(end, refseq->index, true, sequence);
        st_logInfo("cap_construct2:\t%d seconds\n", time(NULL) - startTime);
        //Cap *startcap = cap_construct(end, event);

        startTime = time(NULL);
        cap_check(startcap);
        st_logInfo("cap_check:\t%d seconds\n", time(NULL) - startTime);
        refseq->index++;

        startTime = time(NULL);
        copyRefCapToLowerFlowers(startcap);
        st_logInfo("copyRefCapToLowerFlowers:\t%d seconds\n", time(NULL) - startTime);

        //adding reference Segments to the blocks and creating inherited caps
        //st_logInfo("Adding reference segments and adjacencies...\n");
        startTime = time(NULL);
        reference_walkDown(end, refseq, NULL);
        st_logInfo("reference_walkDown:\t%d seconds\n", time(NULL) - startTime);
        //st_logInfo("Added reference segments and adjacencies successfully.\n");

        //free memory:
        startTime = time(NULL);
        referenceSequence_destruct(refseq);
        st_logInfo("referenceSequence_destruct:\t%d seconds\n", time(NULL) - startTime);
        free(chromHeader);
    }
    reference_destructPseudoChromosomeIterator(it);

    return flower;
}

static void checkAddedReferenceSequence_checkAdjacency(Cap *cap) {
    End *end = cap_getEnd(cap);
    Cap *adjacentCap = cap_getAdjacency(cap);
    assert(adjacentCap != NULL);
    End *adjacentEnd = cap_getEnd(adjacentCap);
    PseudoAdjacency *pseudoAdjacency = end_getPseudoAdjacency(end);
    assert(pseudoAdjacency != NULL);
    assert(pseudoAdjacency_get5End(pseudoAdjacency) != pseudoAdjacency_get3End(pseudoAdjacency));
    assert(end_getPositiveOrientation(end) == pseudoAdjacency_get5End(pseudoAdjacency) || end_getPositiveOrientation(adjacentEnd) == pseudoAdjacency_get5End(pseudoAdjacency));
    assert(end_getPositiveOrientation(end) == pseudoAdjacency_get3End(pseudoAdjacency) || end_getPositiveOrientation(adjacentEnd) == pseudoAdjacency_get3End(pseudoAdjacency));
}

static void checkAddedReferenceSequence(Flower *flower, const char *referenceEventName) {
    /*
     * Checks an added reference sequence is consistent with the reference genome.
     */
    //Check the blocks of the flower have a copy of the reference sequence
    Block *block;
    Flower_BlockIterator *blockIt = flower_getBlockIterator(flower);
    while((block = flower_getNextBlock(blockIt)) != NULL) {
        Segment *segment;
        Block_InstanceIterator *segmentIt = block_getInstanceIterator(block);
        bool b = 0;
        while((segment = block_getNext(segmentIt)) != NULL) {
            if(strcmp(event_getHeader(segment_getEvent(segment)), referenceEventName) == 0) {
                assert(!b);
                b = 1;
                //Check the adjacency of the caps of the reference segment respect the reference structure
                checkAddedReferenceSequence_checkAdjacency(segment_get5Cap(segment));
                checkAddedReferenceSequence_checkAdjacency(segment_get3Cap(segment));
            }
        }
        block_destructInstanceIterator(segmentIt);
        assert(b);
    }
    flower_destructBlockIterator(blockIt);
    //Recurse on the children.
    Group *group;
    Flower_GroupIterator *groupIt = flower_getGroupIterator(flower);
    while((group = flower_getNextGroup(groupIt))) {
        if(!group_isLeaf(group)) { //Now check the adjacencies are linked
            checkAddedReferenceSequence(group_getNestedFlower(group), referenceEventName);
        }
    }
    flower_destructGroupIterator(groupIt);
}

void usage() {
    fprintf(stderr, "cactus_addReferenceSeq, version 0.0\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr, "-b --referenceEventString : Name of the reference event\n");
    fprintf(
            stderr,
            "-c --cactusDisk : The location of the flower disk directory (the databaseString)\n");
    fprintf(stderr, "-h --help : Print this help screen\n");
}

int main(int argc, char *argv[]) {
    char * logLevelString = NULL;
    char * cactusDiskDatabaseString = NULL;
    char * flowerName = "0";
    char * name = NULL;

    while (1) {
        static struct option long_options[] = { { "logLevel",
                required_argument, 0, 'a' }, { "referenceEventString", required_argument, 0,
                'b' }, { "cactusDisk", required_argument, 0, 'c' },
                { "help", no_argument, 0, 'h' }, { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:b:c:h", long_options,
                &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
            case 'a':
                logLevelString = stString_copy(optarg);
                break;
            case 'b':
                name = stString_copy(optarg);
                break;
            case 'c':
                cactusDiskDatabaseString = stString_copy(optarg);
                break;
                //case 'd':
                //    flowerName = stString_copy(optarg);
                //    break;
            case 'h':
                usage();
                return 0;
            default:
                usage();
                return 1;
        }
    }

    //assert(flowerName != NULL);

    //////////////////////////////////////////////
    //Set up logging
    //////////////////////////////////////////////

    if (logLevelString != NULL && strcmp(logLevelString, "INFO") == 0) {
        st_setLogLevel(ST_LOGGING_INFO);
    }
    if (logLevelString != NULL && strcmp(logLevelString, "DEBUG") == 0) {
        st_setLogLevel(ST_LOGGING_DEBUG);
    }

    //////////////////////////////////////////////
    //Load the database
    //////////////////////////////////////////////

    stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(
            cactusDiskDatabaseString);
    
    CactusDisk *cactusDisk = cactusDisk_construct(kvDatabaseConf, 0);
    st_logInfo("Set up the flower disk\n");

    ///////////////////////////////////////////////////////////////////////////
    // Parse the basic reconstruction problem
    ///////////////////////////////////////////////////////////////////////////

    int64_t starttime = time(NULL);
    Flower *flower = cactusDisk_getFlower(cactusDisk,
            cactusMisc_stringToName(flowerName));
    st_logInfo("cactusDisk_getFlower:\t%d seconds\n", time(NULL) - starttime);
    st_logInfo("Parsed the top level flower of the cactus tree to check\n");

    ///////////////////////////////////////////////////////////////////////////
    // Add the sequence
    ///////////////////////////////////////////////////////////////////////////

    int64_t startTime = time(NULL);
    flower = flower_addReferenceSequence(flower, cactusDisk, name);

    //test(flower);
    st_logInfo("Added the reference sequence in %i seconds/\n",
            time(NULL) - startTime);

    ///////////////////////////////////////////////////////////////////////////
    // Check the added sequence
    ///////////////////////////////////////////////////////////////////////////

    checkAddedReferenceSequence(flower, name);
    flower_checkRecursive(flower);

    ///////////////////////////////////////////////////////////////////////////
    // Update the disk
    ///////////////////////////////////////////////////////////////////////////

    //write to Disk:
    startTime = time(NULL);
    cactusDisk_write(cactusDisk);
    st_logInfo("cactusDisk_write:\t%d seconds\n", time(NULL) - startTime);

    ///////////////////////////////////////////////////////////////////////////
    // Clean up.
    ///////////////////////////////////////////////////////////////////////////

    cactusDisk_destruct(cactusDisk);
    stKVDatabaseConf_destruct(kvDatabaseConf);
    free(name);
    if(logLevelString != NULL) {
        free(logLevelString);
    }
    free(cactusDiskDatabaseString);

    //while(1);

    return 0;
}
