#!/usr/bin/env python

#nknguyen@soe.ucsc.edu
#Generate the config "experiment.xml" file needed to run cactus 
#09/21/2010

def writeExperimentXml(options):
    f = open(options.outfile, 'w')
    speciesLi = options.species.split()
    species = ""
    if options.sequenceDir:
        for spc in speciesLi:
            species += options.sequenceDir + "/" + spc + " "
        species = species.strip()
    else:
        species = options.species
        
    f.write("<cactus_workflow_experiment config=\"%s\" sequences=\"%s\" species_tree=\"%s\">\n" % (options.config, species, options.tree))
    f.write("<cactus_disk>\n")
    #f.write("<st_kv_database_conf type=\"tokyo_cabinet\"><tokyo_cabinet database_dir=\"%s\" /></st_kv_database_conf>\n" % options.cactusDisk)
    f.write("%s\n" %options.databaseString)
    f.write("</cactus_disk>\n")
    f.write("</cactus_workflow_experiment>\n")
    f.close()
    return

def main():
    import os
    import sys
    import re
    from optparse import OptionParser
    #Adding Options
    usage = "usage: %prog [options]\n"
    parser = OptionParser(usage=usage)
    parser.add_option("-o", "--output", dest="outfile", help="Output file", default="./experiment.xml")
    parser.add_option("-c", "--databaseString", dest="databaseString", help="cactus database string")
    #parser.add_option("-c", "--cactusDisk", dest="cactusDisk", help="Location of cactusDisk", default="./cactusDisk")
    parser.add_option("-s", "--species", dest="species", help="List of species (required argument).")
    parser.add_option("-d", "--sequenceDir", dest="sequenceDir", help="Directory where the sequences are")
    parser.add_option("-t", "--tree", dest="tree", help="Newick tree of 'species' (required argument).")
    parser.add_option("-b", "--config", dest="config",help="cactus_workflow_config (required argument), values: 'default', 'highSensitivity','highSpecificity', or path to cactus_workflow_config.xml file")
    (options, args) = parser.parse_args()
    if not options.species or not options.tree or not options.config:
        sys.stderr.write(usage)
    else:
        writeExperimentXml(options)

if __name__ == "__main__":
    main()



