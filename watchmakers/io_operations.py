from load import *



def testCreateDirectory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
    else:
        rmtree(directory)
        os.makedirs(directory)

def deleteDirectory(directory):
    if os.path.exists(directory):
        rmtree(directory)
#


def macroGenerator(percentage,isotope,location,runs,events):

    covPCT = {'10pct':0.1,'15pct':0.15,'20pct':0.2,\
    '25pct':0.25,'30pct':0.30,'35pct':0.35,'40pct':0.40}
    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)


    #Part of the macro that is the same for all jobs
    dir = os.getcwd()
#    print arguments
#    print arguments["--depth"]

    depth = float(arguments["--depth"])

    header = '''
/glg4debug/glg4param omit_muon_processes  0.0
/glg4debug/glg4param omit_hadronic_processes  0.0

/rat/db/set DETECTOR experiment "Watchman"
/rat/db/set DETECTOR detector_factory "Watchman"
/rat/db/set WATCHMAN_PARAMS photocathode_coverage %4.2f
%s

/run/initialize

# BEGIN EVENT LOOP
/rat/proc lesssimpledaq
/rat/proc fitbonsai
/rat/proc fitcentroid
/rat/proc fitpath
/rat/proc count
/rat/procset update 1000

# Use IO.default_output_filename
/rat/proclast outroot
/rat/procset file "%s/root_files%s/%s/%s/watchman_%s_%s_%s_%d.root"
#END EVENT LOOP

''' %(covPCT[percentage],additionalMacOpt,dir,additionalMacStr,isotope,percentage,isotope,percentage,location,runs)


    #Part of macro that varies with the various conditions
    if location == 'PMT':
        line1 = '''
/generator/add decaychain %s:regexfill
/generator/pos/set inner_pmts[0-9]+

/run/beamOn %d''' %(isotope,events*5)
    elif location == 'FV':
        line1 = '''
/generator/add decaychain %s:fill:poisson
/generator/pos/set  0 0 0
/generator/rate/set 6.43

/run/beamOn %d''' %(isotope,events*5)
    elif location == 'FN':
        line1 = '''
/generator/add combo fastneutron:regexfill
/generator/pos/set rock_[0-9]+
/generator/vtx/set 0 0 0
/generator/fastneutron/depth %f
/generator/fastneutron/enthresh 10.0
/generator/fastneutron/sidewalls 1.0

/run/beamOn %d'''%(depth,events)
    elif location == 'FNFairport':
        line1 = '''
/generator/add combo fastneutron:regexfill
/generator/pos/set rock_[0-9]+
/generator/vtx/set 0 0 0
/generator/fastneutron/depth 1434.0
/generator/fastneutron/enthresh 10.0
/generator/fastneutron/sidewalls 1.0

/run/beamOn %d'''%(events)
    elif location == 'FNBoulby':
        line1 = '''
/generator/add combo fastneutron:regexfill
/generator/pos/set rock_[0-9]+
/generator/vtx/set 0 0 0
/generator/fastneutron/depth 2805.
/generator/fastneutron/enthresh 10.0
/generator/fastneutron/sidewalls 1.0

/run/beamOn %d'''%(events)
    elif    location=='I':
        line1 = '''
/generator/add combo ibd:fill
/generator/vtx/set  1 0 0
/generator/pos/set 0 0 0

/run/beamOn %d'''%(events)
    elif location == 'S':
        line1 ='''
/generator/add combo spectrum:fill
/generator/vtx/set e+ %s
/generator/pos/set 0 0 0

/run/beamOn %d'''%(isotope,events)
    elif location =='N':
        line1 = '''
/generator/add combo gun2:fill
/generator/vtx/set %s  0 0 0 0 0.001 0.20
/generator/pos/set 0 0 0

/run/beamOn %d'''%(isotope,events)
    elif location == 'RN':
        AZ = isotope
        A =  int(int(AZ)/1000)
        Z = int(AZ) - A*1000
        line1 = '''
/generator/add combo isotope:fill
/generator/pos/set 0 0 0
/generator/vtx/set GenericIon 0 0 0
/generator/isotope/A %s.0
/generator/isotope/Z %s.0
/generator/isotope/E 0.0

/run/beamOn %d''' %(A,Z,events)
    else:
        line1 = 'A'
        print location
    return header+line1



def jobString(percentage,j,runs,models,arguments):
#    directory = "/p/lscratche/adg/Watchboy/simplifiedData/rp_sim/wm"
    directory   = os.getcwd()
    softDir     = "/usr/gapps/adg/geant4/rat_pac_and_dependency"
    ratDir      = os.environ['RATROOT']
    rootDir     = os.environ['ROOTSYS']
    g4Dir       =  os.environ['G4INSTALL']
    watchmakersDir = os.environ['WATCHENV']
    try:
        sheffield   = os.environ['SHEFFIELD']
        # print 'Running on sheffield cluster'
    except:
        # print 'Not running on sheffield cluster'
        sheffield =  0

    software    = "%s/bin/rat" %(ratDir)
    d,iso,loc,coverage,coveragePCT = loadSimulationParameters()

    ele =  d["%s"%(iso[j])]
    location = loc[j]

    goodness     = float(arguments['-g'])
    case = int(arguments['-j'])

    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)


    N            = int(arguments["-N"])
    rate         = float(arguments["-r"])
    timemask     = float(arguments['-t'])*1000.0
    distancemask = float(arguments['-d'])
    goodness     = float(arguments['-g'])
    dirGoodness  = float(arguments['-G'])
    minNHIT      = float(arguments['-T'])
    fIn          = arguments["-f"]
    fidV         = float(arguments["--fv"])
    pmtV         = float(arguments["--psup"])
    tankV        = float(arguments["--tankDis"])
    outF         = arguments["--ntupleout"]
    superNova    = arguments["--supernovaFormat"]

    if sheffield:

        line1 = """#!/bin/sh
    #MSUB -N WM_%s_%s_%d_%s    #name of job
    #MSUB -A adg         # sets bank account
    #MSUB -l nodes=1:ppn=1,walltime=23:59:59,partition=borax  # uses 1 node
    #MSUB -q pbatch         #pool
    #MSUB -o %s/log_case%s%s/wmpc_%s_%s_%d.log
    #MSUB -e %s/log_case%s%s/wmpc_%s_%s_%d.err
    #MSUB -d %s  # directory to run from
    #MSUB -V
    #MSUB                     # no more psub commands

    source %s/bin/thisroot.sh
    source /usr/local/gcc49/setup.sh
    source /usr/local/geant4/setup.sh 10.3
    source %s/geant4make.sh
    source %s/env.sh
    source %s/env_wm.sh
    export G4NEUTRONHP_USE_ONLY_PHOTONEVAPORATION=1
    export SHEFFIELD=1\n
    """ %(percentage,location,runs,additionalMacStr,\
    directory,case,additionalMacStr,percentage,location,runs,\
    directory,case,additionalMacStr,percentage,location,runs,\
    directory,\
    rootDir,g4Dir,ratDir,watchmakersDir)

    else:
        line1 = """#!/bin/sh
    #MSUB -N WM_%s_%s_%d_%s    #name of job
    #MSUB -A adg         # sets bank account
    #MSUB -l nodes=1:ppn=1,walltime=23:59:59,partition=borax  # uses 1 node
    #MSUB -q pbatch         #pool
    #MSUB -o %s/log_case%s%s/wmpc_%s_%s_%d.log
    #MSUB -e %s/log_case%s%s/wmpc_%s_%s_%d.err
    #MSUB -d %s  # directory to run from
    #MSUB -V
    #MSUB                     # no more psub commands

    source %s/bin/thisroot.sh
    source %s/../../../bin/geant4.sh
    source %s/geant4make.sh
    source %s/env.sh
    source %s/env_wm.sh
    export G4NEUTRONHP_USE_ONLY_PHOTONEVAPORATION=1\n
    """ %(percentage,location,runs,additionalMacStr,\
    directory,case,additionalMacStr,percentage,location,runs,\
    directory,case,additionalMacStr,percentage,location,runs,\
    directory,\
    rootDir,g4Dir,g4Dir,ratDir,watchmakersDir)

    for mods in models:
        if location == "FN":
            line1 += "export PHYSLIST=%s\n" %(mods)
        if case == 1 or case == 2 or case == 4:
            _log = "log_case%s%s/%s/%s/rat.%s_%s_%s_%d.log" %(case,additionalMacStr,mods,percentage,percentage,mods,location,runs)
            _mac = "%s/macro%s/%s/%s/run%s_%s_%d.mac" %(directory,additionalMacStr,mods,percentage,mods,location,runs)
            line1 += "%s -l %s %s\n" %(software,_log,_mac)
        if case >= 3:
            fileN = "root_files%s/%s/%s/watchman_%s_%s_%s_%d.root" %(additionalString,mods,percentage,mods,percentage,location,runs)
            if additionalString != "":
                fileNO = "ntuple_root_files%s/%s/%s/watchman_%s_%s_%s%s_%d.root" %(additionalString,mods,percentage,mods,percentage,location,additionalString,runs)
                if sheffield:
                    line1 += "root -b -l -q %s/watchmakers/\'goldenFileExtractor.C(\"%s\",\"%s\",%f,%f,%f,%f,%f,%f,%f,%f,%f\")\'\n" %(watchmakersDir,fileN,fileNO,minNHIT,goodness,dirGoodness,timemask,\
                                    rate,distancemask,fidV,pmtV,tankV)
                    ###int goldenFileExtractor(const char *file, const char *outfile = "null", double nhit_min =3., double goodness_min = 0.1, double goodness_dir = 0.1, double timeWindow_ns = 100000, double rate = 10.0, double maxDistance = 2.0, double fidBound = 5.4, double pmtBound = 6.4, double tankBound = 8.0000) {
                else:
                    line1 += "watch -n %s -f %s --ntupleout %s\n" %(additionalCommands,fileN,fileNO)
    return line1,case



def generateMacros(N,e):
    d,iso,loc,coverage,coveragePCT = loadSimulationParameters()
    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)
    print additionalMacOpt
#    N = int(arguments['-N'])
    print N,e
    ##Clean or create macro directories
    for j in range(len(iso)):
        for ii in d["%s"%(iso[int(j)])]:
            for idx,cover in enumerate(coverage):
                dir = "macro%s/%s/%s" %(additionalMacStr,ii,cover)
#                print dir
                testCreateDirectory(dir)
#    for idx,cover in enumerate(coverage):
#        dir = "%s" %(cover)
#        testCreateDirectory(dir)

    for j in range(len(iso)):
        for ii in d["%s"%(iso[int(j)])]:
            for idx,cover in enumerate(coverage):
                for val in range(N):
                    line = macroGenerator(cover,ii,loc[j],val,e )
                    dir = "macro%s/%s/%s" %(additionalMacStr,ii,cover)

                    outfile = open("%s/run%s_%s_%d.mac" %(dir,ii,\
                    loc[j],val),"wb")
                    outfile.writelines(line)
                    outfile.close
    return 0



def generateJobs(N,arguments):
    d,iso,loc,coverage,coveragePCT = loadSimulationParameters()
    case = arguments["-j"]
    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)

    try:
        sheffield   = os.environ['SHEFFIELD']
        # print 'Running on sheffield cluster'
    except:
        # print 'Not running on sheffield cluster'
        sheffield =  0

    '''Find wheter the jobs folder exist: if not create, if yes clean and recreate'''


    directory = 'jobs_case%s%s'%(case,additionalMacStr)
    if not os.path.exists(directory):
        os.makedirs(directory)
    else:
        rmtree(directory)
        os.makedirs(directory)

    for ii in loc:
        for idx,cover in enumerate(coverage):
            directory = "jobs_case%s%s/%s/%s" %(case,additionalMacStr,ii,cover)
            if not os.path.exists(directory):
                os.makedirs(directory)
            else:
                rmtree(directory)
                os.makedirs(directory)

    '''Find wheter the jobs folder exist: if no create, if yes clean and recreate'''
    directory = 'log_case%s%s'%(case,additionalMacStr)
    if not os.path.exists(directory):
        os.makedirs(directory)
    else:
        rmtree(directory)
        os.makedirs(directory)

    for j in range(len(iso)):
        for ii in d["%s"%(iso[int(j)])]:
            for idx,cover in enumerate(coverage):
                directory = "root_files%s/%s/%s" %(additionalMacStr,ii,cover)
                if not os.path.exists(directory):
                    os.makedirs(directory)

    #
    # for j in range(len(iso)):
    #     for ii in d["%s"%(iso[int(j)])]:
    #         for idx,cover in enumerate(coverage):
    #             directory = "ntuple_root_files/%s/%s" %(ii,cover)
    #             if not os.path.exists(directory):
    #                 os.makedirs(directory)
    #

    for j in range(len(iso)):
        for ii in d["%s"%(iso[int(j)])]:
            for idx,cover in enumerate(coverage):
                directory = "ntuple_root_files%s/%s/%s" %(additionalString,ii,cover)
                if not os.path.exists(directory):
                    os.makedirs(directory)

    for j in range(len(iso)):
        for ii in d["%s"%(iso[int(j)])]:
            for idx,cover in enumerate(coverage):
                directory = "log_case%s%s/%s/%s" %(case,additionalMacStr,ii,cover)
                if not os.path.exists(directory):
                    os.makedirs(directory)


    '''Make sure that the softlink are correct for Bonsai input'''

    ratDir      = os.environ['RATROOT']

    src = ratDir+'/fit_param.dat'
    dst = os.getcwd()+'/fit_param.dat'
    if not os.path.exists(dst):
        os.symlink(src,dst)

    src = ratDir+'/like.bin'
    dst = os.getcwd()+'/like.bin'
    if not os.path.exists(dst):
        os.symlink(src,dst)

    job = 'jobs_case%s%s'%(case,additionalMacStr)

    job_list = '''#!/bin/sh
'''

    for j in range(len(iso)):
        for idx,cover in enumerate(coverage):
            models  = d["%s" %(iso[j])]
            for index in range(N):
                line,case = jobString(cover,j,index,models,arguments)
                stringFile = "%s/%s/%s/jobs%s_%s_%s_%d_case%d.sh" %(job,loc[j],cover,cover,\
                                                            "%s"%(iso[int(j)]),loc[j],index,case)
		if sheffield:
                    job_list+= 'condor_qsub -l nodes=1:ppn=1 ' + stringFile + '\n'
                    outfile = open(stringFile,"wb")
                    outfile.writelines(line)
                else:
		    if index == 0:
                        job_list+= '(msub ' + stringFile +') || ./'+ stringFile + '\n'
                    outfile = open(stringFile,"wb")
                    outfile.writelines(line)
                    
                    if index < N-1:
                        stringFile1 = "(msub %s/%s/%s/jobs%s_%s_%s_%d_case%d.sh || ./%s/%s/%s/jobs%s_%s_%s_%d_case%d.sh)" %(job,loc[j],cover,cover,\
                                                                                                 "%s"%(iso[int(j)]),loc[j],index+1,case,job,loc[j],cover,cover,\
                                                                                                 "%s"%(iso[int(j)]),loc[j],index+1,case)
                    outfile.writelines(stringFile1)
                outfile.close
                os.chmod(stringFile,S_IRWXU)


    outfile = open('sub_jobs_case%s%s'%(case,additionalMacStr),"wb")
    outfile.writelines(job_list)
    outfile.close
    os.chmod('sub_jobs_case%s%s'%(case,additionalMacStr),S_IRWXG)
    os.chmod('sub_jobs_case%s%s'%(case,additionalMacStr),S_IRWXU)
    return 0


def customJob(arguments):

    directory   = os.getcwd()
    softDir     = "/usr/gapps/adg/geant4/rat_pac_and_dependency"
    ratDir      = os.environ['RATROOT']
    rootDir     = os.environ['ROOTSYS']
    g4Dir       =  os.environ['G4INSTALL']
    watchmakersDir = os.environ['WATCHENV']
    sheffield   = os.environ['SHEFFIELD']

    if sheffield:
        print 'Running on sheffield cluster'
    software    = "%s/bin/rat" %(ratDir)


    ''' Custom job for photocoverage 02-2017 analyis'''

    d,iso,loc,coverage,coveragePCT = loadSimulationParameters()

    cnt = 0

    outF2 = open('sub_cust_job',"wb")
    outF2.writelines('#!/bin/sh\n')


    for j in range(len(iso)):
        for ii in d["%s"%(iso[int(j)])]:
            line1 = """#!/bin/sh
#MSUB -N WM_%s_%s    #name of job
#MSUB -A adg         # sets bank account
#MSUB -l nodes=1:ppn=1,walltime=23:59:59,partition=borax  # uses 1 node
#MSUB -q pbatch         #pool
#MSUB -o %s/log/CJwmpc_%s_%s.log
#MSUB -e %s/log/CJwmpc_%s_%s.err
#MSUB -d %s  # directory to run from
#MSUB -V
#MSUB                     # no more psub commands

source %s/bin/thisroot.sh
source %s/../../../bin/geant4.sh
source %s/geant4make.sh
source %s/env.sh
source %s/env_wm.sh
export G4NEUTRONHP_USE_ONLY_PHOTONEVAPORATION=1\n
""" %(ii,loc[j],\
                      directory,ii,loc[j],\
                      directory,ii,loc[j],\
                      directory,\
                      rootDir,g4Dir,g4Dir,ratDir,watchmakersDir)
#            line2 = "watch --extractNtup -N 600 -P %s -L %s \n" %(ii,loc[j])
#            line3 = "watch --extractNtup -N 600 -P %s -L %s --fv 4.302\n" %(ii,loc[j])
#            line4 = "watch --extractNtup -N 600 -P %s -L %s --fv 4.302 -g 0.65\n" %(ii,loc[j])
#            line5 = "watch --extractNtup -N 600 -P %s -L %s --fv 4.302 -g 0.65 -T 12\n" %(ii,loc[j])
#            line6 = "watch --extractNtup -N 600 -P %s -L %s --fv 4.302 -T 12\n" %(ii,loc[j])
#            line7 = "watch --extractNtup -N 600 -P %s -L %s -g 0.65\n" %(ii,loc[j])
#            line8 = "watch --extractNtup -N 600 -P %s -L %s -g 0.65 -T 12\n" %(ii,loc[j])
#            line9 = "watch --extractNtup -N 600 -P %s -L %s -T 12\n" %(ii,loc[j])

            line2 = "watch --extractNtup -N 600 -P %s -L %s \n" %(ii,loc[j])
            line3 = "watch --extractNtup -N 600 -P %s -L %s -g 0.65\n" %(ii,loc[j])
            line4 = "watch --extractNtup -N 600 -P %s -L %s -g 0.65 -T 8\n" %(ii,loc[j])
            line5 = "watch --extractNtup -N 600 -P %s -L %s -T 8\n" %(ii,loc[j])
            line6 = "watch --extractNtup -N 600 -P %s -L %s -g 0.65 -T 10\n" %(ii,loc[j])
            line7 = "watch --extractNtup -N 600 -P %s -L %s -T 10\n" %(ii,loc[j])
            line8 = "watch --extractNtup -N 600 -P %s -L %s -g 0.65 -T 12\n" %(ii,loc[j])
            line9 = "watch --extractNtup -N 600 -P %s -L %s -T 12\n" %(ii,loc[j])
            line10 = "watch --extractNtup -N 600 -P %s -L %s -g 0.65 -T 14\n" %(ii,loc[j])
            line11 = "watch --extractNtup -N 600 -P %s -L %s -T 14\n" %(ii,loc[j])


            outfile = open('jobs/sub_jobs__%d_2'%(cnt),"wb")
            outF2.writelines('msub jobs/sub_jobs__%d_2\n'%(cnt))
            outfile.writelines(line1)
            outfile.writelines(line2)
            outfile.close

            outfile = open('jobs/sub_jobs__%d_3'%(cnt),"wb")
            outF2.writelines('msub jobs/sub_jobs__%d_3\n'%(cnt))
            outfile.writelines(line1)
            outfile.writelines(line3)
            outfile.close

            outfile = open('jobs/sub_jobs__%d_4'%(cnt),"wb")
            outF2.writelines('msub jobs/sub_jobs__%d_4\n'%(cnt))
            outfile.writelines(line1)
            outfile.writelines(line4)
            outfile.close

            outfile = open('jobs/sub_jobs__%d_5'%(cnt),"wb")
            outF2.writelines('msub jobs/sub_jobs__%d_5\n'%(cnt))
            outfile.writelines(line1)
            outfile.writelines(line5)
            outfile.close

            outfile = open('jobs/sub_jobs__%d_6'%(cnt),"wb")
            outF2.writelines('msub jobs/sub_jobs__%d_6\n'%(cnt))
            outfile.writelines(line1)
            outfile.writelines(line6)
            outfile.close

            outfile = open('jobs/sub_jobs__%d_7'%(cnt),"wb")
            outF2.writelines('msub jobs/sub_jobs__%d_7\n'%(cnt))
            outfile.writelines(line1)
            outfile.writelines(line7)
            outfile.close

            outfile = open('jobs/sub_jobs__%d_8'%(cnt),"wb")
            outF2.writelines('msub jobs/sub_jobs__%d_8\n'%(cnt))
            outfile.writelines(line1)
            outfile.writelines(line8)
            outfile.close

            outfile = open('jobs/sub_jobs__%d_9'%(cnt),"wb")
            outF2.writelines('msub jobs/sub_jobs__%d_9\n'%(cnt))
            outfile.writelines(line1)
            outfile.writelines(line9)
            outfile.close

            outfile = open('jobs/sub_jobs__%d_10'%(cnt),"wb")
            outF2.writelines('msub jobs/sub_jobs__%d_10\n'%(cnt))
            outfile.writelines(line1)
            outfile.writelines(line10)
            outfile.close

            outfile = open('jobs/sub_jobs__%d_11'%(cnt),"wb")
            outF2.writelines('msub jobs/sub_jobs__%d_11\n'%(cnt))
            outfile.writelines(line1)
            outfile.writelines(line11)
            outfile.close

            cnt+=1
    outF2.close





def deleteAllWorkDirectories():
    d,iso,loc,coverage,coveragePCT = loadSimulationParameters()

    dir = "log"
    if os.path.exists(dir):
        rmtree(dir)

    dir = "jobs"
    if os.path.exists(dir):
        rmtree(dir)

    for idx,cover in enumerate(coverage):
        dir = "macro"
        if os.path.exists(dir):
            rmtree(dir)

    if os.path.exists('fit_param.dat'):
        os.remove('fit_param.dat')

    if os.path.exists('like.bin'):
        os.remove('like.bin')

    if os.path.exists('sub_jobs'):
        os.remove('sub_jobs')



def mergeFiles():
    # Read external requirements
    #arguments = docopt.docopt(docstring)
    # Load internal requirements
    d,iso,loc,coverage,coveragePCT = loadSimulationParameters()
    trees = {}
    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)
    pathFinal = "root_files%s/merged_ntuple_watchman" % (additionalMacStr)
    for j in range(len(iso)):
        for ii in d["%s"%(iso[int(j)])]:
            for idx,cover in enumerate(coverage):
                t_name  = "data_%s_%s_%s"%(ii,cover,loc[j])
                trees[t_name] = TChain("data")

                s = "ntuple_root_files%s/%s/%s/watchman_%s_%s_%s_*.root" %(additionalMacStr,ii,cover,ii,cover,loc[j])
                sw = "%s_%s_%s_%s.root"%(pathFinal,ii,cover,loc[j])

                print "Writing ", sw,"from",s
                trees[t_name].Add(s)
                print "Number of entries ",trees[t_name].GetEntries()
                trees[t_name].Merge(sw)
                del trees[t_name]
    del trees
    return 0

def testEnabledCondition(arguments):

    additionalString      = ""
    additionalCommands    = ""

    additionalMacStr = ""
    additionalMacOpt = ""

# Commands required for root_file

    if (arguments['--detectMedia']):
        additionalMacOpt += "/rat/db/set GEO[detector] material \"%s\"\n" %(arguments['--detectMedia'])
        additionalMacStr += "_detectorMedia_%s" %(arguments['--detectMedia'])
        additionalString += "_detectorMedia_%s" %(arguments['--detectMedia'])

    if (arguments['--collectionEff']):
        additionalMacOpt += "/rat/db/set GEO[inner_pmts] efficiency_correction %f\n" %(float(arguments['--collectionEff']))
        additionalMacStr += "_collectionEfficiency_%f" %(float(arguments['--collectionEff']))
        additionalString += "_collectionEfficiency_%f" %(float(arguments['--collectionEff']))



#    if (arguments['--fidThick'] and arguments['--shieldThick'] and arguments['--halfHeight'] and arguments['--tankRadius']):
#        additionalMacOpt += "/rat/db/set GEO[fiducial] r_max %f\n" %(((float(arguments['--tankRadius']))-1.5875)-((float(arguments['--shieldThick']))+(float(arguments['--fidThick']))))
#        additionalMacOpt += "/rat/db/set GEO[fiducial] size_z %f\n" %(((float(arguments['--halfHeight']))-1.5875)-((float(arguments['--shieldThick']))+(float(arguments['--fidThick']))))

    if (arguments['--pmtModel']):
        additionalMacOpt += "/rat/db/set GEO[inner_pmts] pmt_model \"%s\"\n" %((arguments['--pmtModel']))
        additionalMacStr += "_pmtModel_%s" %((arguments['--pmtModel']))
        additionalString += "_pmtModel_%s" %((arguments['--pmtModel']))

    if (arguments['--photocath'] and not arguments['--pmtModel']):
        additionalMacOpt += "/rat/db/set PMT[r7081pe]  photocathode_surface \"photocathode_%s\"\n" %((arguments['--photocath']))
        additionalMacStr += "_photocathode_%s" %((arguments['--photocath']))
        additionalString += "_photocathode_%s" %((arguments['--photocath']))

    if (arguments['--photocath'] and arguments['--pmtModel']):
        additionalMacOpt += "/rat/db/set PMT[%s]  photocathode_surface \"photocathode_%s\"\n" %(arguments['--pmtModel'],arguments['--photocath'])
        additionalMacStr += "_photocathode_%s" %((arguments['--photocath']))
        additionalString += "_photocathode_%s" %((arguments['--photocath']))

    baseValue = 7
    #Analysis strings, usually shows up in ntuple processing
#    print defaultValues[baseValue+1]
    if float(arguments['-r'])          != defaultValues[baseValue+1]:
        additionalString += "_rate_%f" %(float(arguments['-r']))
        additionalCommands += " -r %f " %(float(arguments['-r']))

    if float(arguments['-d'])          != defaultValues[baseValue+2]:
        additionalString += "_deltaR_%f" %(float(arguments['-d']))
        additionalCommands += " -d %f" %(float(arguments['-d']))

    if float(arguments['-t'])          != defaultValues[baseValue+3]:
        additionalString += "_deltaT_%f" %(float(arguments['-t']))
        additionalCommands +=  " -t %f" %(float(arguments['-t']))

    if float(arguments['-T'])            != (defaultValues[baseValue+4]):
        additionalString += "_n9Min_%d" %(int(arguments['-T']))
        additionalCommands += " -T %d" %(int(arguments['-T']))

    if float(arguments['-g'])          != defaultValues[baseValue+5]:
        additionalString += "_posGood_%f" %(float(arguments['-g']))
        additionalCommands += " -g %f" %(float(arguments['-g']))

    if float(arguments['-G'])          != defaultValues[baseValue+6]:
        additionalString += "_dirGood_%f" %(float(arguments['-G']))
        additionalCommands += " -G %f" %(float(arguments['-G']))

    if float(arguments['--fv'])        !=  defaultValues[baseValue+7]:
        additionalString += "_FVboundary_%f" %(float(arguments['--fv']))
        additionalCommands +=  "--fv %f" %(float(arguments['--fv']))

    if float(arguments['--psup'])      != defaultValues[baseValue+8]:
        additionalString += "_PMTboundary_%f" %(float(arguments['--psup']))
        additionalCommands += "--psup %f" %(float(arguments['--psup']))

    if float(arguments['--tankDis'])   != defaultValues[baseValue+9]:
        additionalString += "_Tankboundary_%f" %(float(arguments['--tankDis']))
        additionalCommands +=" --tankDis %f" %(float(arguments['--tankDis']))


    if float(arguments['--tankRadius']) != defaultValues[baseValue+10]:
        additionalMacOpt += "/rat/db/set GEO[tank] r_max %f\n" %(float(arguments['--tankRadius']))
        additionalMacOpt += "/rat/db/set GEO[detector] r_max %f\n" %(float(arguments['--tankRadius'])-1.5875)
        additionalMacOpt += "/rat/db/set GEO[shield] detector_size_d %f\n" %(float(arguments['--tankRadius'])*2)
        additionalMacStr += "_tankRadius_%f" %(float(arguments['--tankRadius']))
        additionalString += "_tankRadius_%f" %(float(arguments['--tankRadius']))

    if float(arguments['--halfHeight'])!= defaultValues[baseValue+11]:
        additionalMacOpt += "/rat/db/set GEO[tank] size_z %f\n" %(float(arguments['--halfHeight']))
        additionalMacOpt += "/rat/db/set GEO[shield] detector_size_z %f\n" %(float(arguments['--halfHeight'])*2)
        additionalMacOpt += "/rat/db/set GEO[detector] size_z %f\n" %(float(arguments['--halfHeight'])-1.5875)
        additionalMacOpt += "/rat/db/set GEO[cables] size_z %f\n" %(float(arguments['--halfHeight'])-1.5875)
        additionalMacStr += "_halfHeight_%f" %(float(arguments['--halfHeight']))
        additionalString += "_halfHeight_%f" %(float(arguments['--halfHeight']))

    if float(arguments['--shieldThick'])!= defaultValues[baseValue+12]:
        additionalMacOpt += "/rat/db/set GEO[shield] shield_thickness %f\n" %(float(arguments['--shieldThick']))
        additionalMacStr += "_shieldThickness_%f" %(float(arguments['--shieldThick']))
        additionalString += "_shieldThickness_%f" %(float(arguments['--shieldThick']))

    if float(arguments['--steelThick'])!= defaultValues[baseValue+13]:
        additionalMacOpt += "/rat/db/set GEO[shield] steel_thickness %f\n" %(float(arguments['--steelThick']))
        additionalMacStr += "_steelThickness_%f" %(float(arguments['--steelThick']))
        additionalString += "_steelThickness_%f" %(float(arguments['--steelThick']))

    if float(arguments['--fidThick'])!= defaultValues[baseValue+14]:
        additionalString += "_fidThickness_%f" %(float(arguments['--fidThick']))
        additionalCommands +=" --fidThick %f" %(float(arguments['--fidThick']))



    if int(arguments['--supernovaFormat']):
        additionalString += "_supernovaFormat"
        additionalCommands +=" --supernovaFormat "

    if additionalString == "":
        additionalString = "_default"



    return  additionalString,additionalCommands,additionalMacStr,additionalMacOpt





def writeResultsToFile(s,g,h):
    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)
    _str = "ntuple_root_files%s/%s" %(additionalString,s)
    d,iso,loc,coverage,coveragePCT = loadSimulationParameters()
    f_root = TFile(_str,"recreate")
    for gE in g:
        g["%s"%(gE)].Write()
    for hE in h:
        h["%s"%(hE)].Write()
    f_root.Close()


def mergeNtupleFiles(arguments):
    # Read external requirements
    #arguments = docopt.docopt(docstring)
    # Load internal requirements
    d,iso,loc,coverage,coveragePCT = loadSimulationParameters()
    trees = {}

    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)


    pathFinal = "ntuple_root_files%s/merged_ntuple_watchman" %(additionalString)

    if arguments["-P"] and arguments["-L"]:
        ii      = arguments["-P"]
        locj    = arguments["-L"]
        for idx,cover in enumerate(coverage):
            t_name  = "data_%s_%s_%s"%(ii,cover,locj)
            try:
                trees[t_name] = TChain("data")

                s = "ntuple_root_files%s/%s/%s/watchman_%s_%s_%s_*.root" %(additionalString,ii,cover,ii,cover,locj)
                sw = "%s_%s_%s_%s.root"%(pathFinal,ii,cover,locj)

                print "Writing ", sw,"from",s
                trees[t_name].Add(s)
                print "Number of entries ",trees[t_name].GetEntries()
                trees[t_name].Merge(sw)
                del trees[t_name]
            except:
                print 'Error for %s' %(t_name)

    if (arguments["-P"] and not arguments["-L"]) or (arguments["-L"] and not arguments["-P"]):
        print "arguments -L and -P must be used at the same time, for now"


    if (not arguments["-P"] and not arguments["-L"]):
        for j in range(len(iso)):
            for ii in d["%s"%(iso[int(j)])]:
                for idx,cover in enumerate(coverage):
                    t_name  = "data_%s_%s_%s"%(ii,cover,loc[j])
                    try:
                        trees[t_name] = TChain("data")

                        s = "ntuple_root_files%s/%s/%s/watchman_%s_%s_%s_*.root" %(additionalString,ii,cover,ii,cover,loc[j])
                        sw = "%s_%s_%s_%s.root"%(pathFinal,ii,cover,loc[j])

                        print "Writing ", sw,"from",s
                        trees[t_name].Add(s)
                        print "Number of entries ",trees[t_name].GetEntries()
                        trees[t_name].Merge(sw)
                        del trees[t_name]
                    except:
                        print 'Error for %s' %(t_name)


    del trees
    return 0





def extractNtuple(arguments):

    N            = int(arguments["-N"])
    rate         = float(arguments["-r"])
    timemask     = float(arguments['-t'])*1000.0
    distancemask = float(arguments['-d'])
    goodness     = float(arguments['-g'])
    dirGoodness  = float(arguments['-G'])
    minNHIT      = float(arguments['-T'])
    fIn          = arguments["-f"]
    fidV         = float(arguments["--fv"])
    pmtV         = float(arguments["--psup"])
    tankV        = float(arguments["--tankDis"])
    outF         = arguments["--ntupleout"]
    superNova    = arguments["--supernovaFormat"]

    print file
    d,iso,loc,coverage,coveragePCT = loadSimulationParameters()
    if not superNova:
        try:
            goldenFileExtractor(fIn,outF,minNHIT,goodness,dirGoodness,timemask,\
                            rate,distancemask,fidV,pmtV,tankV)
        except:
            print "Error.."
    else:
        try:
            supernovaAnalysis(fIn,outF)
        except:
            print "Error.."

def extractNtupleALL(arguments):
    d,iso,loc,coverage,coveragePCT = loadSimulationParameters()

    N            = int(arguments["-N"])

    rate         = float(arguments["-r"])
    timemask     = float(arguments['-t'])*1000.0
    distancemask = float(arguments['-d'])
    goodness     = float(arguments['-g'])
    dirGoodness  = float(arguments['-G'])
    minNHIT      = float(arguments['-T'])
    fidV         = float(arguments["--fv"])
    pmtV         = float(arguments["--psup"])
    tankV        = float(arguments["--tankDis"])

    superNova    = arguments["--supernovaFormat"]

    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)


    for j in range(len(iso)):
        for ii in d["%s"%(iso[int(j)])]:
            for idx,cover in enumerate(coverage):
                directory = "ntuple_root_files%s/%s/%s" %(additionalString,ii,cover)
                if not os.path.exists(directory):
                    os.makedirs(directory)



#    for j in range(len(iso)):
#            for ii in d["%s"%(iso[int(j)])]:
#                for idx,cover in enumerate(coverage):
#                    for run in range(N):
#                        fIn =  "root_files/%s/%s/watchman_%s_%s_%s_%d.root" %(ii,cover,ii,cover,loc[j],run)
#                        fOut = "ntuple_root_files%s/%s/%s/watchman_%s_%s_%s_%d.root" %(additionalString,ii,cover,ii,cover,loc[j],run)
#                        print fIn, " -> ", fOut
#                        goldenFileExtractor(fIn,minNHIT,goodness,dirGoodness,timemask,\
#                                                rate,distancemask,fidV,pmtV,tankV,fOut)



    if arguments["-P"] and arguments["-L"]:
        ii      = arguments["-P"]
        locj    = arguments["-L"]
        for idx,cover in enumerate(coverage):
            for run in range(N):
                fIn =  "root_files%s/%s/%s/watchman_%s_%s_%s_%d.root" %(additionalMacStr,ii,cover,ii,cover,locj,run)
                fOut = "ntuple_root_files%s/%s/%s/watchman_%s_%s_%s_%d.root" %(additionalString,ii,cover,ii,cover,locj,run)
                if os.path.isfile(fIn) and not os.path.isfile(fOut):
                    print fIn, " -> ", fOut
                    if not superNova:
                        try:
                            goldenFileExtractor(fIn,fOut,minNHIT,goodness,dirGoodness,timemask,\
                                            rate,distancemask,fidV,pmtV,tankV)
                        except:
                            print "Error.."
                    else:
                        try:
                            supernovaAnalysis(fIn,fOut)
                        except:
                            print "Error.."


    if (arguments["-P"] and not arguments["-L"]) or (arguments["-L"] and not arguments["-P"]):
        print "arguments -L and -P must be used at the same time, for now"


    if (not arguments["-P"] and not arguments["-L"]):
        for j in range(len(iso)):
            for ii in d["%s"%(iso[int(j)])]:
                for idx,cover in enumerate(coverage):
                    for run in range(N):
                        fIn =  "root_files%s/%s/%s/watchman_%s_%s_%s_%d.root" %(additionalMacStr,ii,cover,ii,cover,loc[j],run)
                        fOut = "ntuple_root_files%s/%s/%s/watchman_%s_%s_%s_%d.root" %(additionalString,ii,cover,ii,cover,loc[j],run)
                        if os.path.isfile(fIn) and not os.path.isfile(fOut):
                            print fIn, " -> ", fOut
                            if not superNova:
                                try:
                                    goldenFileExtractor(fIn,fOut,minNHIT,goodness,dirGoodness,timemask,\
                                                    rate,distancemask,fidV,pmtV,tankV)
                                except:
                                    print "Error.."
                            else:
                                try:
                                    supernovaAnalysis(fIn,fOut)
                                except:
                                    print "Error.."
