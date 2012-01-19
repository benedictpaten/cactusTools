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
#include <ctype.h>

#include "cactus.h"
#include "avl.h"
#include "commonC.h"
#include "hashTableC.h"
#include "cactusTraversal.h"


/*
 * Library for generating mafs from cactus.
 */

static char *formatSequenceHeader(Sequence *sequence) {
    const char *sequenceHeader = sequence_getHeader(sequence);
    if (strlen(sequenceHeader) > 0) {
        char *cA = st_malloc(sizeof(char) * (1 + strlen(sequenceHeader)));
        sscanf(sequenceHeader, "%s", cA);
        return cA;
    } else {
        return cactusMisc_nameToString(sequence_getName(sequence));
    }
}

static char *getSegmentStringShowingOnlySubstitutionsWithRespectToTheReference(Segment *segment) {
    char *string = segment_getString(segment);
    assert(string != NULL);
    Block *block = segment_getBlock(segment);
    Block_InstanceIterator *it = block_getInstanceIterator(block);
    Segment *segment2;
    while((segment2 = block_getNext(it)) != NULL) {
        if(segment2 != segment && strcmp(cactusMisc_getDefaultReferenceEventHeader(), event_getHeader(segment_getEvent(segment2))) == 0) {
            assert(segment != segment_getReverse(segment2));
            char *string2 = segment_getString(segment2);
            assert(string2 != NULL);
            assert(strlen(string) == strlen(string2));
            for(int32_t i=0; i<strlen(string); i++) {
                if(toupper(string[i]) == toupper(string2[i])) {
                    string[i] = '*';
                }
            }
            free(string2);
            break;
        }
    }
    block_destructInstanceIterator(it);
    return string;
}

static void getMAFBlockP2(Segment *segment, FILE *fileHandle, char *(*getString)(Segment *segment)) {
    assert(segment != NULL);
    Sequence *sequence = segment_getSequence(segment);
    if (sequence != NULL) {
        char *sequenceHeader = formatSequenceHeader(sequence);
        int32_t start;
        if (segment_getStrand(segment)) {
            start = segment_getStart(segment) - sequence_getStart(sequence);
        } else { //start with respect to the start of the reverse complement sequence
            start = (sequence_getStart(sequence) + sequence_getLength(sequence)
                    - 1) - segment_getStart(segment);
        }
        int32_t length = segment_getLength(segment);
        char *strand = segment_getStrand(segment) ? "+" : "-";
        int32_t sequenceLength = sequence_getLength(sequence);
        char *instanceString = getString(segment); //segment_getString(segment);
        fprintf(fileHandle, "s\t%s\t%i\t%i\t%s\t%i\t%s\n", sequenceHeader,
                start, length, strand, sequenceLength, instanceString);
        free(instanceString);
        free(sequenceHeader);
    }
}

static void getMAFBlockP(Segment *segment, FILE *fileHandle, char *(*getString)(Segment *segment)) {
    int32_t i;
    for (i = 0; i < segment_getChildNumber(segment); i++) {
        getMAFBlockP(segment_getChild(segment, i), fileHandle, getString);
    }
    getMAFBlockP2(segment, fileHandle, getString);
}

static int32_t getNumberOnPositiveStrand(Block *block) {
    Block_InstanceIterator *it = block_getInstanceIterator(block);
    Segment *segment;
    int32_t i = 0;
    while ((segment = block_getNext(it)) != NULL) {
        if (segment_getChildNumber(segment) == 0) {
            if (segment_getStrand(segment)) {
                i++;
            }
        }
    }
    block_destructInstanceIterator(it);
    return i;
}

static void getMAFBlock2(Block *block, FILE *fileHandle, char *(*getString)(Segment *segment)) {
    //void getMAFBlock(Block *block, FILE *fileHandle, ReferenceSequence *referenceSequence) {
    /*
     * Outputs a MAF representation of the block to the given file handle.
     */
    //Correct the orientation..
    if (getNumberOnPositiveStrand(block) == 0) {
        block = block_getReverse(block);
    }
    if (block_getInstanceNumber(block) > 0) {
        //Add in the header
        if (block_getRootInstance(block) != NULL) {
            /* Get newick tree string with internal labels and no unary events */
            char *newickTreeString = block_makeNewickString(block, 1, 0);
            assert(newickTreeString != NULL);
            fprintf(fileHandle, "a score=%i tree='%s'\n",
                    block_getLength(block) * block_getInstanceNumber(block),
                    newickTreeString);
            free(newickTreeString);
        } else {
            fprintf(fileHandle, "a score=%i\n",
                    block_getLength(block) * block_getInstanceNumber(block));
        }
        //Now for the reference segment
        /*if (referenceSequence != NULL) {
         char *instanceString = getConsensusString(block);
         fprintf(fileHandle, "s\t%s\t%i\t%i\t%s\t%i\t%s\n",
         referenceSequence->header, referenceSequence->index,
         block_getLength(block), "+", referenceSequence->length,
         instanceString);
         free(instanceString);
         referenceSequence->index += block_getLength(block);
         }*/
        //Now add the blocks in
        if (block_getRootInstance(block) != NULL) {
            assert(block_getRootInstance(block) != NULL);
            getMAFBlockP(block_getRootInstance(block), fileHandle, getString);
            fprintf(fileHandle, "\n");
        } else {
            Block_InstanceIterator *iterator = block_getInstanceIterator(block);
            Segment *segment;
            while ((segment = block_getNext(iterator)) != NULL) {
                getMAFBlockP2(segment, fileHandle, getString);
            }
            block_destructInstanceIterator(iterator);
            fprintf(fileHandle, "\n");
        }
    }
}

void getMAFBlock(Block *block, FILE *fileHandle) {
    getMAFBlock2(block, fileHandle, segment_getString);
}

void getMAFBlockShowingOnlySubstitutionsWithRespectToTheReference(Block *block, FILE *fileHandle) {
    getMAFBlock2(block, fileHandle, getSegmentStringShowingOnlySubstitutionsWithRespectToTheReference);
}

void(*prepMafBlockFn)(Block *, FILE *);

void prepMafBlock(stList *caps, FILE *fileHandle) {
    Cap *cap = stList_get(caps, 0);
    assert(cap_getSide(cap));
    if(cap_getSegment(cap)) {
        prepMafBlockFn(segment_getBlock(cap_getSegment(cap)), fileHandle);
    }
}

void getMAFsReferenceOrdered2(const char *referenceEventString, Flower *flower,
        FILE *fileHandle, void(*getMafBlockFn)(Block *, FILE *)) {
    /*
     * Outputs MAF representations of all the block in the flower and its descendants, ordered
     * according to the reference ordering.
     */
    Event *referenceEvent = eventTree_getEventByHeader(flower_getEventTree(flower), referenceEventString);
    End *end;
    Flower_EndIterator *endIt = flower_getEndIterator(flower);
    prepMafBlockFn = getMafBlockFn;
    while ((end = flower_getNextEnd(endIt)) != NULL) {
        if (end_isStubEnd(end) && end_isAttached(end)) {
            Cap *cap = getCapForReferenceEvent(end, event_getName(referenceEvent)); //The cap in the reference
            assert(cap != NULL);
            assert(cap_getSequence(cap) != NULL);
            cap = cap_getStrand(cap) ? cap : cap_getReverse(cap);
            if(!cap_getSide(cap)) {
                traverseCapsInSequenceOrderFrom3PrimeCap(cap, fileHandle, NULL, (void (*)(stList *, void *))prepMafBlock);
            }
        }
    }
    flower_destructEndIterator(endIt);
}

void getMAFsReferenceOrdered(Flower *flower,
        FILE *fileHandle, void(*getMafBlockFn)(Block *, FILE *)) {
    getMAFsReferenceOrdered2(cactusMisc_getDefaultReferenceEventHeader(), flower, fileHandle, getMafBlockFn);
}

void getMAFs(Flower *flower, FILE *fileHandle,
        void(*getMafBlock)(Block *, FILE *)) {
    /*
     * Outputs MAF representations of all the block sin the flower and its descendants.
     */

    //Make MAF blocks for each block
    Flower_BlockIterator *blockIterator = flower_getBlockIterator(flower);
    Block *block;
    while ((block = flower_getNextBlock(blockIterator)) != NULL) {
        getMafBlock(block, fileHandle);
        //getMAFBlock(block, fileHandle, NULL);
    }
    flower_destructBlockIterator(blockIterator);

    //Call child flowers recursively.
    Flower_GroupIterator *groupIterator = flower_getGroupIterator(flower);
    Group *group;
    while ((group = flower_getNextGroup(groupIterator)) != NULL) {
        if (!group_isLeaf(group)) {
            getMAFs(group_getNestedFlower(group), fileHandle, getMafBlock); //recursive call.
        }
    }
    flower_destructGroupIterator(groupIterator);
}

void makeMAFHeader(Flower *flower, FILE *fileHandle) {
    fprintf(fileHandle, "##maf version=1 scoring=N/A\n");
    char *cA = eventTree_makeNewickString(flower_getEventTree(flower));
    fprintf(fileHandle, "# cactus %s\n\n", cA);
    free(cA);
}
