/*
 * cactusTraversal.h
 *
 *  Created on: 19 Jan 2012
 *      Author: benedictpaten
 */

#ifndef CACTUS_TRAVERSAL_H_
#define CACTUS_TRAVERSAL_H_

void traverseCapsInSequenceOrderFrom3PrimeCap(Cap *cap, void *extraArg,
        void(*_3PrimeFn)(stList *caps, void *extraArg),
        void(*_5PrimeFn)(stList *caps, void *extraArg));

Cap *getCapForReferenceEvent(End *end, Name referenceEventName);

#endif /* CACTUS_TRAVERSAL_H_ */
