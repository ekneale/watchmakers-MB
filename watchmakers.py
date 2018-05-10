#!/usr/bin/python

from watchmakers.load import *
from watchmakers.io_operations import *
from watchmakers.analysis import *
from watchmakers.sensitivity import *
from watchmakers.data import *
from watchmakers.analysis_pmt import *

######################## Start of main function ###########################

if __name__ == "__main__":

    print docstring
    print arguments
#     print defaultValues
    print ""

    if arguments['-D']:
        deleteAllWorkDirectories()

    if arguments['-m']:
        generateMacros(int(arguments['-N']),int(arguments['-e']))

    if arguments['-j']:
        generateJobs(int(arguments['-N']),arguments)

    if arguments['-n']:
        extractNtuple(arguments)

    if arguments['--extractNtup']:
        extractNtupleALL(arguments)

    if arguments['-M']:
        mergeNtupleFiles(arguments)

    if arguments['-A']:
        g,h = {},{}
        g,h = runAnalysisProcess(arguments["-f"],g,h)
        writeResultsToFile(arguments["-o"],g,h)

    if  arguments["--efficiency"]:
        extractHistogramWitCorrectRate()

    if arguments['--pmtanalysis']:
	extractPMTHistograms()

    if arguments['--pmtanalysis2']:
	extractPMTSinglesHistograms()

    if arguments['--pmtanalysis3']:
	extractIBDHistograms()

    if arguments['--pmtanalysis4']:
	extractPMTn9Histograms()
	
    if arguments['--pmtanalysis5']:
	extractPMTSinglesFidCutHistograms()

    if arguments['--pmtanalysis6']:
	extractPMTAccidentalsFidCutHistograms()

    if arguments['--sensitivity']:
        sensitivityMapNew()

    if arguments['--customJob']:
        customJob(arguments)

    if arguments['--pass1Trigger']:
        performPass1(arguments)

    if arguments['--pass2Trigger']:
        performPass2(arguments)

    if arguments['--fileDict']:
        createFileDictionary(arguments)


######################## Waba Luba Dub Dub!! ###########################
