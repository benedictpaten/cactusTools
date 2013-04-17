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
#include "cactusMafs.h"
#include "cactusUtils.h"
//#include "cactus_addReferenceSeq.h"

static void usage() {
    fprintf(stderr, "cactus_mafGenerator, version 0.2\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr,
            "-c --cactusDisk : The location of the flower disk directory\n");
    fprintf(stderr,
            "-d --flowerName : The name of the flower (the key in the database)\n");
    fprintf(stderr, "-e --outputFile : The file to write the MAFs in.\n");
    fprintf(
            stderr,
            "-g --referenceEventString : String identifying the reference event. This option will include a reference sequence in the blocks.\n");
    fprintf(stderr, "-h --help : Print this help screen\n");
    fprintf(stderr,
            "-i --showOnlySubstitutionsWithRespectToTheReference : Display only substitutions with respect to the reference.\n");
}

int main(int argc, char *argv[]) {
    /*
     * Arguments/options
     */
    char * logLevelString = NULL;
    char * cactusDiskDatabaseString = NULL;
    char * flowerName = NULL;
    char * outputFile = NULL;
    char *referenceEventString = (char *)cactusMisc_getDefaultReferenceEventHeader();
    bool showOnlySubstitutionsWithRespectToTheReference = 0;

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel",
                required_argument, 0, 'a' }, { "cactusDisk", required_argument,
                0, 'c' }, { "flowerName", required_argument, 0, 'd' }, {
                "outputFile", required_argument, 0, 'e' }, {
                "referenceEventString", optional_argument, 0, 'g' }, { "help",
                no_argument, 0, 'h' },
                { "showOnlySubstitutionsWithRespectToTheReference", no_argument, 0, 'i' }, { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:c:d:e:g:hi", long_options,
                &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
            case 'a':
                logLevelString = stString_copy(optarg);
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
            case 'g':
                referenceEventString = stString_copy(optarg);
                break;
            case 'h':
                usage();
                return 0;
            case 'i':
                showOnlySubstitutionsWithRespectToTheReference = 1;
                break;
            default:
                usage();
                return 1;
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    // (0) Check the inputs.
    ///////////////////////////////////////////////////////////////////////////

    assert(flowerName != NULL);
    assert(outputFile != NULL);

    //////////////////////////////////////////////
    //Set up logging
    //////////////////////////////////////////////

    st_setLogLevelFromString(logLevelString);

    //////////////////////////////////////////////
    //Log (some of) the inputs
    //////////////////////////////////////////////

    st_logInfo("Flower name : %s\n", flowerName);
    st_logInfo("Output MAF file : %s\n", outputFile);

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

    Flower *flower = cactusDisk_getFlower(cactusDisk, cactusMisc_stringToName(
            flowerName));
    st_logInfo("Parsed the top level flower of the cactus tree to check\n");

    ///////////////////////////////////////////////////////////////////////////
    // Recursive check the flowers.
    ///////////////////////////////////////////////////////////////////////////

    int64_t startTime = time(NULL);
    FILE *fileHandle = fopen(outputFile, "w");
    makeMAFHeader(flower, fileHandle);

    if(eventTree_getEventByHeader(flower_getEventTree(flower), referenceEventString) == NULL) {
        st_logInfo("No reference event found, so not ordering by reference\n", referenceEventString);
        getMAFs(flower, fileHandle, getMAFBlock);
    }
    else {
        st_logInfo("Ordering by reference by string %s\n", referenceEventString);
        if(showOnlySubstitutionsWithRespectToTheReference) {
            getMAFsReferenceOrdered2(referenceEventString, flower, fileHandle, getMAFBlockShowingOnlySubstitutionsWithRespectToTheReference);
        } else {
            getMAFsReferenceOrdered2(referenceEventString, flower, fileHandle, getMAFBlock);
        }
    }
    fclose(fileHandle);
    st_logInfo("Got the mafs in %" PRIi64 " seconds/\n", time(NULL) - startTime);

    ///////////////////////////////////////////////////////////////////////////
    // Clean up.
    ///////////////////////////////////////////////////////////////////////////

    cactusDisk_destruct(cactusDisk);
    stKVDatabaseConf_destruct(kvDatabaseConf);

    return 0;
}
