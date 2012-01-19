/*
 * cactusTraversal.c
 *
 *  Created on: 19 Jan 2012
 *      Author: benedictpaten
 */

/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <stdio.h>
#include <ctype.h>
#include "cactus.h"
#include "sonLib.h"

////////////////////////////////////
////////////////////////////////////
//Functions to traverse the caps of a sequence
////////////////////////////////////
////////////////////////////////////

static Cap *getCapUp(Cap *cap) {
    /*
     * Gets the highest level version of a cap.
     */
    while (1) {
        assert(cap != NULL);
        cap = cap_getSide(cap) ? cap : cap_getReverse(cap);
        if (end_isBlockEnd(cap_getEnd(cap))) {
            return cap;
        }
        Group *parentGroup = flower_getParentGroup(
                end_getFlower(cap_getEnd(cap)));
        if (parentGroup == NULL) {
            return cap;
        }
        cap = flower_getCap(group_getFlower(parentGroup), cap_getName(cap));
    }
}

static stList *getCapsDown(Cap *cap, bool side) {
    /*
     * Gets a list of a cap in progressively lower level problems, such that
     * the first member of the return list will be the given cap, and the last member
     * will be the lowest level version of the cap.
     */
    stList *caps = stList_construct();
    assert(end_isAttached(cap_getEnd(cap)) || end_isBlockEnd(cap_getEnd(cap)));
    if (end_isStubEnd(cap_getEnd(cap))) {
        assert(flower_getParentGroup(end_getFlower(cap_getEnd(cap))) == NULL);
    }
    while (1) {
        stList_append(caps,
                cap_getSide(cap) == side ? cap : cap_getReverse(cap));
        assert(end_getGroup(cap_getEnd(cap)) != NULL);
        Flower *nestedFlower = group_getNestedFlower(
                end_getGroup(cap_getEnd(cap)));
        if (nestedFlower != NULL) {
            assert(
                    flower_getEnd(nestedFlower, end_getName(cap_getEnd(cap)))
                            != NULL);
            cap = flower_getCap(nestedFlower, cap_getName(cap));
            assert(cap != NULL);
        } else {
            break;
        }
    }
    return caps;
}

void traverseCapsInSequenceOrderFrom3PrimeCap(Cap *cap, void *extraArg,
        void(*_3PrimeFn)(stList *caps, void *extraArg),
        void(*_5PrimeFn)(stList *caps, void *extraArg)) {
    assert(end_isStubEnd(cap_getEnd(cap)));
    assert(end_isAttached(cap_getEnd(cap)));
    assert(flower_getParentGroup(end_getFlower(cap_getEnd(cap))) == NULL);
    while (1) {
        //Call 3' function
        stList *caps = getCapsDown(cap, 0);
        if (_3PrimeFn != NULL) {
            _3PrimeFn(caps, extraArg);
        }
        //Get the adjacent 5 prime cap
        assert(group_isLeaf(end_getGroup(cap_getEnd(stList_peek(caps)))));
        cap = getCapUp(cap_getAdjacency(stList_peek(caps)));
        stList_destruct(caps);
        //Now call 5' function
        caps = getCapsDown(cap, 1);
        if (_5PrimeFn != NULL) {
            _5PrimeFn(caps, extraArg);
        }
        stList_destruct(caps);
        if (cap_getSegment(cap) != NULL) { //Get the opposite 3 prime cap.
            cap = cap_getOtherSegmentCap(cap);
            assert(cap != NULL);
        } else {
            assert(end_isStubEnd(cap_getEnd(cap)));
            assert(end_isAttached(cap_getEnd(cap)));
            assert(
                    flower_getParentGroup(end_getFlower(cap_getEnd(cap)))
                            == NULL);
            break;
        }
    }
}


Cap *getCapForReferenceEvent(End *end, Name referenceEventName) {
    /*
     * Get the cap for a given event.
     */
    End_InstanceIterator *it = end_getInstanceIterator(end);
    Cap *cap;
    while ((cap = end_getNext(it)) != NULL) {
        if (event_getName(cap_getEvent(cap)) == referenceEventName) {
            end_destructInstanceIterator(it);
            return cap;
        }
    }
    end_destructInstanceIterator(it);
    assert(0);
    return NULL;
}
