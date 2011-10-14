/*
 * mafs.h
 *
 *  Created on: 24 Feb 2011
 *      Author: benedictpaten
 */

#ifndef MAFS_H_
#define MAFS_H_

void getMAFBlock(Block *block, FILE *fileHandle);

void getMAFBlockShowingOnlySubstitutionsWithRespectToTheReference(Block *block, FILE *fileHandle);

void getMAFsReferenceOrdered2(const char *referenceEventName, Flower *flower,
        FILE *fileHandle, void(*getMafBlockFn)(Block *, FILE *));

void getMAFsReferenceOrdered(Flower *flower,
        FILE *fileHandle, void(*getMafBlockFn)(Block *, FILE *));

void getMAFs(Flower *flower, FILE *fileHandle, void (*getMafBlock)(Block *, FILE *));

void makeMAFHeader(Flower *flower, FILE *fileHandle);

#endif /* MAFS_H_ */
