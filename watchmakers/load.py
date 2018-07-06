
import os.path

from stat import S_IRWXG,S_IRWXU
from shutil import rmtree
import warnings

import numpy as np
from numpy import sqrt
from numpy import array as npa
from numpy import power,absolute,logical_and,column_stack,zeros,empty,append,\
sqrt,absolute,recarray


from math import pow,exp,log10,pi

try:
    # from root_numpy import root2rec,array2tree,array2root,tree2array
    from rootpy.plotting import Canvas,Hist,Hist2D,Graph
    from rootpy.plotting.style import set_style
    from rootpy.io import root_open
    #from rootpy.interactive import wait
    #    set_style('ATLAS')

    warnings.simplefilter("ignore")
except:
    print "Could not load in root_numpy or rootpy, they are required to run this module."

defaultValues  = [1,3,2500,2805.,'merged_ntuple_watchman',\
'merged_ntuple_watchman','null', 'processed_watchman.root',\
10.,2.0, 100.0, 9, 0.65,0.1,8026.35,8026.35,1600.0,6.35,1000.,\
'day','boulby', 1.0, 0.043, 0.133,50.,10.,32.,0.002]

docstring = """
    Usage: watchmakers.py [options]

    Arguments:

    Options:
    -D                  Delete all current photocoverage directory.
    --newVers           Major revision to Watchmakers. By default off for old results
    --force             Forcing the recreation of the root_file,bonsai_root_file and log folders
    --noRoot            Allows to generate scripts without loading in any ROOT module
    -j=<jobType>        Create submision scripts (1,2,4:rat-pac files|case >3 ntuplefiles) [default %d]
                        >3 option will generate a nutple_root_files_flags folder for results
    -m                  Also generate macro files
    -N=<N>              Number of MC script that were run [Default: %d]
    -e=<runBeamEntry>   Number of entries per macro (U/Th event x5) [Default: %d]
    --depth=<depthD>    Depth of detector (for fast neutron spectra) [Default: %f]

    --evalRate          Evaluate the rates and combine to efficiency histrograms
    --findRate          Find rate as a function of input flags.
    -n                  generate ntuple from single rat-pac root files
    --extractNtup       generate ntuple from all rat-pac root files
    -f=<ifile>          Input file [Default: %s]
    --ft=<ifile2>        Input file type [Default: %s]
    --ntupleout=<outN>  Name of ntuple out [Default: %s]
    -o=<outputfile>     Efficiency output file [Default: %s]
    --supernovaFormat   Record supernova files instead of golden files
    --pass1Trigger      Process rat-pac files with pass1 triggering
    --pass2Trigger      Process pass1 files with pass2 conditions
    --pass2Flag=<p2f>   Use pass2 files in analysis [Default: 1]
    --fileDict          Create a dicionary of exisiting files
    --PDFs              Extract PDFs with series of cuts

    -r=<rate>           rate of accidentals in hz [Default: %f]
    -d=<distance>       Maximal distance between two events (m) [Default: %f]
    -t=<time>           Maximal time between two events (micro) [Default: %f]
    -T=<tubes>          Minimal number of tubes hit [Default: %d]
    --minN9=<_MPE>      Minimal number of photoelectron [Default: 9.]
    -g=<goodness>       Bonsai position goodness parameter [Default: %f]
    -G=<Goodness>       Bonsai direction goodness parameter [Default: %f]
    --RNRedux=<_RNR>    Reduction due to time/spatial veto arround (>0.90) [Default: 0.9]

    -P=<proc>           Pick a single physics process to analyis/merge (used for ntup)
    -L=<loc>            Pick a single physics location to analyis/merge (used for ntup)
    -C=<cov>            Pick a single coverage

    --customJob         Custom job for photocoverage 02-2017

    --tankRadius=<TR>   Total radius of tank (mm) [Default: %f]
    --halfHeight=<HH>   Half height of tank (mm) [Default: %f]
    --shieldThick=<ST>  Steel->PMT distance (mm) [Default: %f]
    --steelThick=<StT>  Steel Thickness (mm)     [Default: %f]
    --fidThick=<fT>     Fiducial volume-> PMT Thickness (mm) [Default: %f]

    -M                  Merge result files from trial ntuples

    --efficiency        Read merged files and perform analysis
    -A                  Read merged files and perform analysis

    -R                  Read analyzed result and evaluate sensitivity
    --sensitivity       Read analyzed results and evaluate sensitivity
    --sensIBD           Read, analyze IBD with different cuts on prompt and delayed

    --timeScale=<_ts>   Integration period (sec,day,month,year) [Default: %s]
    --site=<_site>      Site of the experiment (boulby,fairport) [Default: %s]
    --OnOff=<_OOratio>  Ratio of reactor on to reactor off [Default: %d]
    --cores=<_cores>    Number of cores to discover [Default: 1]
    --bkgdSys=<_BSys>   Systematic background percentage [Default: 0.2]
    --40MWth            Option to sensitivity to do case scenarios
    --40MWthSig=<_SL>   Sigma discovery [Default: 3.0]

    --U238_PPM=<_Uppm>  Concentration of U-238 in glass [Default: %f]
    --Th232_PPM=<_Thp>  Concentration of Th-232 in glass [Default: %f]
    --K_PPM=<_K>        Concentration of K-40 in glass [Default: 16.0]
    --U238_Gd=<_U238Gd>    Activity of U238 in Gd sample mBq/kg [Default: %f ]
    --Th232_Gd=<_ThGd>     Activity of Th232 in Gd sample mBq/kg [Default: %f]
    --U235_Gd=<_U235Gd>    Activity of U235 in Gd sample mBq/kg [Default: %f]
    --Rn222=<_Rn>       Radon activity in water SK 2x10^-3 Bq/m^3 [Default: %f]

    --detectMedia=<_dM>  Detector media (doped_water,...)
    --collectionEff=<CE> Collection efficiency (e.g.: 0.85,0.67,0.475)

    --pmtModel=<_PMTM>   PMT Model (r7081pe for 10inch or r11780_hqe for 12inch)
    --photocath =<_PC>  PMT photocathode (R7081HQE)

    """ % (defaultValues[0],defaultValues[1],defaultValues[2],defaultValues[3],defaultValues[4],\
           defaultValues[5],defaultValues[6],defaultValues[7],defaultValues[8],\
           defaultValues[9],defaultValues[10],defaultValues[11],defaultValues[12],\
           defaultValues[13],defaultValues[14],defaultValues[15],defaultValues[16],\
           defaultValues[17],defaultValues[18],defaultValues[19],defaultValues[20],\
           defaultValues[21],defaultValues[22],defaultValues[23],defaultValues[24],\
           defaultValues[25],defaultValues[26],defaultValues[27])

try:
    import docopt
    arguments = docopt.docopt(docstring)
    print 'using docopt as the user control interface'
except ImportError:
    print 'docopt is not a recognized module, it is required to run this module'


if arguments['--noRoot']:
    print 'Not loading any ROOT modules. usefull for generating files on oslic'

else:
    from ROOT import TRandom3
    from ROOT import TChain,TGraph,TGraphErrors,gSystem,gROOT,TH1D,TH2D,TFile,TCanvas,TF1
    from ROOT import THStack,Double
    from ROOT import kRed,kBlue,kGreen,kCyan,kOrange

    from ROOT import kOrange as kO,kBlue as kB,kGreen as kG
    from ROOT import kMagenta as kM,kAzure as kA,kRed as kR
    from ROOT import TCanvas,TLine, TLatex
    from ROOT import gStyle,gPad,TPaletteAxis

    gSystem.Load("$RATROOT/lib/libRATEvent")
    gSystem.AddIncludePath(" -I$RATROOT/include")

    gROOT.LoadMacro("$WATCHENV/watchmakers/goldenFileExtractor.C")
    from ROOT import goldenFileExtractor

    gROOT.LoadMacro("$WATCHENV/watchmakers/pass1Trigger.C")
    from ROOT import pass1Trigger

    gStyle.SetOptStat(1112211)



# This is deprecated
# gROOT.LoadMacro("$WATCHENV/watchmakers/supernovaAnalysis.C")
# from ROOT import supernovaAnalysis

def loadSimulationParameters():
    #Chain and subsequent isotopes
    d = {}
    # Water and PMT contamination
    # d['CHAIN_238U_NA'] =['238U','234Pa','214Pb','214Bi','210Bi','210Tl']
    # d['CHAIN_232Th_NA'] = ['232Th','228Ac','212Pb','212Bi','208Tl']
    # d['CHAIN_222Rn_NA'] = ['222Rn','214Pb','214Bi','210Bi','210Tl']

    d['CHAIN_238U_NA'] =['234Pa','214Pb','214Bi','210Bi','210Tl']
    d['CHAIN_232Th_NA'] = ['228Ac','212Pb','212Bi','208Tl']
    d['40K_NA']         = ['40K']
    d['CHAIN_222Rn_NA'] = ['214Pb','214Bi','210Bi','210Tl']
    d['TANK_ACTIVITY'] = ['60Co','137Cs']


    # Radioisotope that should have beta-Neutron modes, (beta only generated)
    #    A = ['16','17','18','17','18','8','9','11']
    #    Z = ['6','6','6','7','7','2','3','3']
    #Reduced selection
    A,Z = ['9','11'],['3','3']
    ZA = A
    for i in range(len(A)):
        ZA[i] = str(int(A[i])*1000 +int(Z[i]))
    d['A_Z'] =  ZA

    #Oscillated spectrum at Boulby and IMB site
    d['ibd'] = ['boulby']
    d['N'] = ['neutron']
    d['IBD'] = ['IBD']
    # Fast neutron contamination
    d['FN'] = ['QGSP_BERT_EMV','QGSP_BERT_EMX','QGSP_BERT','QGSP_BIC','QBBC',\
    'QBBC_EMZ','FTFP_BERT','QGSP_FTFP_BERT']
    #
    #List of all physics process clumped together
    iso = ['CHAIN_238U_NA','CHAIN_232Th_NA','CHAIN_222Rn_NA',\
    'CHAIN_222Rn_NA','A_Z','ibd','FN','N','IBD']
    loc = ['PMT','PMT','PMT','FV','RN','S','FN','N','I']


    #Photocoverage selected
    coverage = ['10pct','15pct','20pct','25pct','30pct','35pct','40pct']
    coveragePCT = {'10pct':9.86037,'15pct':14.887,'20pct':19.4453,\
    '25pct':24.994,'30pct':28.8925,'35pct':34.3254,'40pct':39.1385}

    if arguments['-C']:
        coverage = [arguments['-C']]

    return d, iso,loc,coverage,coveragePCT

def loadSimulationParametersNew():
    #Chain and subsequent isotopes
    d = {}

    d['CHAIN_238U_NA'] =['234Pa','214Pb','214Bi','210Bi','210Tl']
    d['CHAIN_232Th_NA'] = ['228Ac','212Pb','212Bi','208Tl']
    d['CHAIN_235U_NA'] = ['231Th','223Fr','211Pb','211Bi','207Tl']
    d['40K_NA']         = ['40K']
    d['CHAIN_222Rn_NA'] = ['214Pb','214Bi','210Bi','210Tl']
    d['TANK_ACTIVITY'] = ['60Co','137Cs']

    d['ibd_p'] = ['promptPositron']
    d['ibd_n'] = ['delayedNeutron']
    d['pn_ibd']   = ['promptDelayedPair']
    # Fast neutron contamination
    d['FN'] = ['QGSP_BERT_EMV','QGSP_BERT_EMX','QGSP_BERT','QGSP_BIC','QBBC',\
    'QBBC_EMZ','FTFP_BERT','QGSP_FTFP_BERT']
    A,Z = ['9','11'],['3','3']
    ZA = A
    for i in range(len(A)):
        ZA[i] = str(int(A[i])*1000 +int(Z[i]))
    d['A_Z'] =  ZA

    process = {'40K_NA':['WaterVolume','PMT','CONCRETE','GUNITE','ROCK'],\
    'CHAIN_238U_NA':['PMT','CONCRETE','GUNITE','ROCK','GD'],\
    'CHAIN_232Th_NA':['PMT','CONCRETE','GUNITE','ROCK','GD'],\
    'CHAIN_235U_NA':['GD'],\
    'CHAIN_222Rn_NA':['WaterVolume'],\
    'TANK_ACTIVITY':['TANK'],\
    'FN':['ROCK'],\
    'A_Z':['WaterVolume'],\
    'ibd_p':['WaterVolume'],\
    'ibd_n':['WaterVolume'],\
    'pn_ibd':['WaterVolume']}

    #Photocoverage selected
    coverage = ['15pct','20pct','25pct','30pct','35pct','SuperK','WatchmanSphere']

    if arguments['-C']:
        coverage = [arguments['-C']]

    return d,process,coverage



def loadPMTInfo():
    import subprocess
    from io_operations import testEnabledCondition
    conditions = testEnabledCondition(arguments)
    cond = conditions[2]
    cmd = ["""grep 'generated PMTs' log_case*%s/boulby/**pct/rat.**pct_boulby_S_0.log"""%(cond)]
    a =  subprocess.check_output(cmd,shell=True)
    b = a.splitlines()
    c = []
    for _b in b:
        c.append(float(_b.split()[3]))

    cmd = ["""grep 'actual photocathode coverage' log_case*%s/boulby/**pct/rat.**pct_boulby_S_0.log"""%(cond)]
    a =  subprocess.check_output(cmd,shell=True)
    b = a.splitlines()
    d = []
    for _b in b:
        d.append(float(_b.split()[4])*100.)

    return c,d


def loadAnalysisParameters(timeScale='day'):

    pmt = loadPMTInfo()

    # Default units are in sec. Conversion factor are below
    timeSec     = 1.0/365./24./3600.

    # Number of free proton
    if timeScale == 'sec':
        timeS   = 1.0
    if timeScale == 'day':
        timeS   = 24.0*3600.
    if timeScale == 'month':
        timeS   = 365.0/12.*24.0*3600.
    if timeScale == 'year':
        timeS   = 365.0*24.0*3600.

    #PMT mass in kilograms
    mass = 1.4 # from Hamamatsu tech details

    ### This had been changed by M.B. from Tamzin implementation
    fidRadius = float(arguments['--tankRadius'])-float(arguments['--steelThick'])-float(arguments['--shieldThick'])-float(arguments['--fidThick'])
    fidHeight = float(arguments['--halfHeight'])-float(arguments['--steelThick'])-float(arguments['--shieldThick'])-float(arguments['--fidThick'])

    tankRadius  = float(arguments["--tankRadius"])-float(arguments['--steelThick'])
    tankHeight  = float(arguments["--halfHeight"])-float(arguments['--steelThick'])

    fidVolume  = pi*pow(fidRadius/1000.,2)*(2.*fidHeight/1000.)
    tankVolume = pi*pow(tankRadius/1000.,2)*(2.*tankHeight/1000.)
    FVkTonRatio = fidVolume/tankVolume
    #print "Change in load.py"

    #Evaluate FV to total detector volume ratio
    nKiloTons   = tankVolume/1000.
    FreeProtons = 0.668559
    TNU         = FreeProtons* nKiloTons *timeSec
    #FVkTonRatio = pow(float(arguments['--fv']),3)/pow(float(arguments['--tankDis']),3)

    #Fast neutrons conversion
    #Rock mass
    # Original Estimate
   # volumeR         = (2.*22.5*23.8*1.0+2.*17*23.8*1.0+2.*22.5*17.*1.0)
    volumeR         = pi*pow(14,2)*28) - pi*pow(13,2)*25.5) # 'tube'-shaped rock layer surrounding cavern: 1m thick, 1.5m thick on top
#    volumeR         = power(22.,3)-power(20.,3)# Rock cavern e.g. (22m x 22m x 22m) - (20m x 20m x 20m)
    density         = 2.39 #from McGrath
    rockMass        = volumeR*power(100.,3)*density


    #Mass of rock evalyated
    avgMuon         = npa([180.,264.])
    avgMuonNC       = power(avgMuon,0.849)
    avgNFluxMag     = 1e-6
    muonRate        = npa([7.06e-7,4.09e-8]) # mu/cm2/s
    tenMeVRatio     = npa([7.51/34.1,1.11/4.86])
    fastNeutrons    = rockMass*avgMuonNC*avgNFluxMag*muonRate*tenMeVRatio

    avgRNYieldRC    = power(avgMuon,0.73)
    skRNRate        = 0.5e-7 # 1/mu/g cm2
    avgMuonSK       = power(219.,0.73)
    skMuFlux        = 1.58e-7 #mu/cm2/sec
    radionuclideRate= (skRNRate*avgRNYieldRC/avgMuonSK)*muonRate*nKiloTons*1e9


    # boulbyIBDRate   = 1120.8*.4/.6 *TNU #//924.48*TNU Taken from website, average corrected
    boulbyIBDRate   = 800.*TNU #new values from geoneutrinos.org
    fairportIBDRate = 7583.*TNU

    inta        = ['si','so','eo','ei']

    dAct        = {}

    #Add the U-238 chain
    M_U238      = 3.953e-25
    Lambda_U238 = 4.916e-18
    PPM_U238    = float(arguments["--U238_PPM"])
    ActivityU238= Lambda_U238*PPM_U238/M_U238/1e6
#    _proc       = ['238U','234Pa','214Pb','214Bi','210Bi','210Tl']
#    _loca       = ['PMT','PMT',  'PMT',  'PMT',  'PMT',  'PMT']
#    acc         = ['chain','acc',  'acc',  'acc',  'acc',  'acc']
#    _br         = [1.0,1.0,     1.0,    1.0,   1.0 ,   0.002]
#    _site        = ['','',      '',     '',     '',     '']
#Changed for Tamzin, as we do not use 238U chain, but it's component
    _proc       = ['234Pa','214Pb','214Bi','210Bi','210Tl']
    _loca       = ['PMT',  'PMT',  'PMT',  'PMT',  'PMT']
    acc         = ['acc',  'acc',  'acc',  'acc',  'acc']
    _br         = [1.0,     1.0,    1.0,   1.0 ,   0.002]
    _site        = ['',      '',     '',     '',     '']
    proc        = _proc
    loca        = _loca
    br          = _br
    site        = _site
    #    decayCnst   = [2.9e-5,  4.31e-4,  5.81e-4,   1.601e-6 , 0.00909]
    arr         = empty(5)
    arr[:]      = ActivityU238
    for index,ele in enumerate(_proc):
        dAct["%s_%s"%(ele,_loca[index])] = ActivityU238*_br[index]*timeS

    Activity    = arr


    #Add the Th-232 chain
    M_Th232      = 3.853145e-25 #kg
    Lambda_Th232 = 1.57e-18 #1/s
    PPM_Th232    = float(arguments["--Th232_PPM"])
    ActivityTh232 = Lambda_Th232*PPM_Th232/M_Th232/1e6
    #    print ActivityU238,ActivityTh232
#Changed for Tamzin, as we do not use 238U chain, but it's component
#    _proc        =['232Th','228Ac','212Pb','212Bi','208Tl']
#    _loca        =['PMT'  ,'PMT',   'PMT', 'PMT',  'PMT'  ]
#    acc          +=['chain'  ,'acc',   'acc', 'acc',  'acc'  ]
#    _br          = [1.0,     1.0,    1.0,   1.0 ,   1.0]
    #    decayCnst   += [1.57e-18,3.3e-5,1.8096e-5, 1.908e-4, 0.003784]
    _site        = ['',      '',     '',     '',     '']
    _proc        =['228Ac','212Pb','212Bi','208Tl']
    _loca        =['PMT',   'PMT', 'PMT',  'PMT'  ]
    acc          +=['acc',   'acc', 'acc',  'acc'  ]
    _br          = [ 1.0,    1.0,   1.0 ,   1.0]
    #    decayCnst   += [1.57e-18,3.3e-5,1.8096e-5, 1.908e-4, 0.003784]
    _site        = ['',     '',     '',     '']

    arr         = empty(4)
    arr[:]      = ActivityTh232
    Activity    = append(   Activity,arr)

    proc        += _proc
    loca        += _loca
    br          += _br
    site        +=_site
    for index,ele in enumerate(_proc):
        dAct["%s_%s"%(ele,_loca[index])] = ActivityTh232*_br[index]*timeS



    #Add the Rn-222 chain
    # N_Rn222     = 2e-3 # Bq/m3
    ActivityRn222     = float(arguments["--Rn222"])*nKiloTons*1e3 #required tons, not ktons
#    _proc       =['222Rn','214Pb','214Bi','210Bi','210Tl']
#    _loca       =['FV', 'FV',   'FV',   'FV',   'FV']
#    acc         +=['chain','acc',  'acc',  'acc',   'acc']
#    _br         = [1.0, 1.0,   1.0,   1.0,     0.002]
    #    decayCnst   += [ 4.31e-4,  5.81e-4,   1.601e-6 , 0.00909]
#    _site        = ['', '',     '',     '',     '']

    _proc       =['214Pb','214Bi','210Bi','210Tl']
    _loca       =['FV',   'FV',   'FV',   'FV']
    acc         +=['acc',  'acc',  'acc',   'acc']
    _br         = [1.0,   1.0,   1.0,     0.002]
    #    decayCnst   += [ 4.31e-4,  5.81e-4,   1.601e-6 , 0.00909]
    _site        = ['', '',     '',     '',     '']

    arr = empty(4)
    arr[:]      = ActivityRn222
    Activity    = append(   Activity,arr)
    proc        += _proc
    loca        += _loca
    br          += _br
    site        += _site
    for index,ele in enumerate(_proc):
        dAct["%s_%s"%(ele,_loca[index])] = ActivityRn222*_br[index]*timeS



    #Add the neutrino signal
    _proc        =['imb','imb','boulby','boulby']
    _loca        =['S','S',     'S',  'S']
    acc         +=['di', 'corr', 'di', 'corr']
    _br          = [1.0,  1.0, 1.0 , 1.0]
    site        += [''   ,'' , '', '']
    arr         = npa([fairportIBDRate,fairportIBDRate,boulbyIBDRate,boulbyIBDRate])
    Activity    = append(    Activity,arr)

    proc        += _proc
    loca        += _loca
    br          += _br
    for index,ele in enumerate(_proc):
        dAct["%s_%s"%(ele,_loca[index])] = arr[index]*timeS

    # Add the neutron
    _proc        =['neutron','neutron']
    _loca        =['N',     'N']
    acc         +=['corr',  'corr']
    _br          = [1.0,   1.0]
    arr         = npa([fairportIBDRate,boulbyIBDRate])
    #    print "Neutrino activity ",arr*timeS/nKiloTons
    Activity    = append(    Activity,arr)
    _site       = [ '','boulby']
    site        += _site
    proc        += _proc
    loca        += _loca
    br          += _br
    for index,ele in enumerate(_proc):
        dAct["%s_%s%s"%(ele,_loca[index],_site[index])] = arr[index]*timeS

    # add a fast neutron at Fairport
    _proc        = ['QGSP_BERT_EMV','QGSP_BERT_EMX','QGSP_BERT','QGSP_BIC',\
    'QBBC','QBBC_EMZ','FTFP_BERT','QGSP_FTFP_BERT']
    _loca        =  ['FN','FN','FN','FN','FN','FN','FN','FN']
    acc         +=  ['corr','corr','corr','corr','corr','corr','corr','corr']
    _br          =  [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0]
    _site        =  ['','','','','','','','']
    arr = empty(8)
    arr[:]      = fastNeutrons[0]
    Activity    = append(Activity,arr)
    proc        += _proc
    loca        += _loca
    br          += _br
    site        += _site
    for index,ele in enumerate(_proc):
        dAct["%s_%s%s"%(ele,_loca[index],_site[index])] = arr[index]*timeS


    # add a fast neutron at Boulby
    _proc        = ['QGSP_BERT_EMV','QGSP_BERT_EMX','QGSP_BERT','QGSP_BIC',\
    'QBBC','QBBC_EMZ','FTFP_BERT','QGSP_FTFP_BERT']
    _loca        =  ['FN','FN','FN','FN','FN','FN','FN','FN']
    acc         +=  ['corr','corr','corr','corr','corr','corr','corr','corr']
    _br          =  [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0]
    arr = empty(8)
    arr[:]      = fastNeutrons[1]
    Activity    = append(Activity,arr)
    _site        = ['boulby','boulby','boulby','boulby','boulby','boulby',\
    'boulby','boulby']
    proc        += _proc
    loca        += _loca
    br          += _br
    site        += _site
    for index,ele in enumerate(_proc):
        dAct["%s_%s%s"%(ele,_loca[index],_site[index])] = arr[index]*timeS


    # Read in the different radionuclide
    _proc        =  ['9003', '11003']
    _loca        =  ['RN','RN']
    acc         +=  ['di','di']
    #normalised to 9Li from SK
    arr         = npa([1.9,0.01])/1.9
    arr         *= radionuclideRate[0]
    Activity    = append(Activity,arr)
    _br         =  [0.495,0.927]
    _site       = ['','']

    proc        += _proc
    loca        += _loca
    br          += _br
    site        += _site
    for index,ele in enumerate(_proc):
        dAct["%s_%s%s"%(ele,_loca[index],_site[index])] = arr[index]*timeS

    _proc        =  ['9003', '11003']
    _loca        =  ['RN','RN']
    acc         +=  ['di','di']
    #normalised to 9Li from SK
    arr         = npa([1.9,0.01])/1.9
    arr         *= radionuclideRate[1]
    Activity    = append(Activity,arr)
    _br         =  [0.495,0.927]
    _site       = ['boulby','boulby']

    proc        += _proc
    loca        += _loca
    br          += _br
    site        += _site
    for index,ele in enumerate(_proc):
        dAct["%s_%s%s"%(ele,_loca[index],_site[index])] = arr[index]*timeS

    _proc        =['IBD','IBD']
    _loca        =['I',     'I']
    acc         +=['corr',  'corr']
    _br          = [1.0,   1.0]
    arr         = npa([fairportIBDRate,boulbyIBDRate])
#    print "Neutrino activity ",arr*timeS/nKiloTons
    Activity    = append(    Activity,arr)
    _site       = [ '','boulby']
    site        += _site
    proc        += _proc
    loca        += _loca
    br          += _br
    for index,ele in enumerate(_proc):
        dAct["%s_%s%s"%(ele,_loca[index],_site[index])] = arr[index]*timeS


    coveNumber    = {'10pct':pmt[0][0],   '15pct':pmt[0][1], '20pct':pmt[0][2],  \
    '25pct':pmt[0][3],  '30pct':pmt[0][4], '35pct':pmt[0][5],  '40pct':pmt[0][6]}
    covePCT       = {'10pct':pmt[1][0], '15pct':pmt[1][1],'20pct':pmt[1][2],\
    '25pct':pmt[1][3],'30pct':pmt[1][4],'35pct':pmt[1][5],'40pct':pmt[1][6]}

    pctTubes   = {"%s"%(pmt[1][0]):pmt[0][0],"%s"%(pmt[1][1]):pmt[0][1],\
    "%s"%(pmt[1][2]):pmt[0][2],"%s"%(pmt[1][3]):pmt[0][3],\
    "%s"%(pmt[1][4]):pmt[0][4],"%s"%(pmt[1][5]):pmt[0][5],\
    "%s"%(pmt[1][6]):pmt[0][6]}

    pct = npa([ float(pmt[1][0]),float(pmt[1][1]),float(pmt[1][2]),\
    float(pmt[1][3]),float(pmt[1][4]),float(pmt[1][5]),\
    float(pmt[1][6]) ])

    return inta,proc,loca,acc,arr,Activity,br,site,timeS,\
    boulbyIBDRate*FVkTonRatio,mass,dAct,coveNumber,covePCT,pctTubes,pct



def loadActivity():

    d,process,coverage = loadSimulationParametersNew()

    ##Evaluate the total mass of PMT glass in kg
    mass, diameter = 1.4, 10.0/0.039 #inch_per_mm # from Hamamatsu tech details
    areaPerPMT = pi*diameter*diameter/4.
    pmtRadius = float(arguments['--tankRadius'])-float(arguments['--steelThick'])-float(arguments['--shieldThick'])
    pmtHeight = float(arguments['--halfHeight'])-float(arguments['--steelThick'])-float(arguments['--shieldThick'])
    psupArea = (2*pmtHeight)*2*pi*pmtRadius + 2.*(pi*pmtRadius**2)
    numPMTs = psupArea/areaPerPMT
    cPMTs = [float(s.strip('pct'))/100.*numPMTs for s in coverage]
    mPMTs = [s*mass for s in cPMTs]

    # print "Num pmts", cPMTs
    # print "Total Mass",mPMTs

    M_U238,Lambda_U238,Abund_U238 = 3.953e-25,4.916e-18,0.992745
    PPM_U238    = float(arguments["--U238_PPM"])
    ActivityU238= Lambda_U238*PPM_U238/M_U238/1e6
    mPMTsU238 = [s*ActivityU238 for s in mPMTs]
    print 'U238',mPMTsU238, ', PPM:',PPM_U238

    M_Th232,Lambda_Th232,Abund_Th232 = 3.853145e-25, 1.57e-18,1.0
    PPM_Th232    = float(arguments["--Th232_PPM"])
    ActivityTh232= Lambda_Th232*PPM_Th232/M_Th232/1e6
    mPMTsTh232 = [s*ActivityTh232 for s in mPMTs]
    print 'Th232',mPMTsTh232, ', PPM:',PPM_Th232

    M_K,Lambda_K,Abund_K = 6.636286e-26,1.842e-18,0.00117
    PPM_K    = float(arguments["--K_PPM"])
    ActivityK= Lambda_K*PPM_K/M_K/1e6
    mPMTsK = [s*ActivityK for s in mPMTs]
    print 'K',mPMTsK, ', PPM:',PPM_K

    print

    return mPMTs




def loadPMTActivity():

    d,process,coverage = loadSimulationParametersNew()

    ##Evaluate the total mass of PMT glass in kg
    mass, diameter = 1.4, 10.0/0.039 #inch_per_mm # from Hamamatsu tech details
    areaPerPMT = pi*diameter*diameter/4.
    pmtRadius = float(arguments['--tankRadius'])-float(arguments['--steelThick'])-float(arguments['--shieldThick'])
    pmtHeight = float(arguments['--halfHeight'])-float(arguments['--steelThick'])-float(arguments['--shieldThick'])
    psupArea = (2*pmtHeight)*2*pi*pmtRadius + 2.*(pi*pmtRadius**2)
    numPMTs = psupArea/areaPerPMT
    cPMTs = [float(s.strip('pct'))/100.*numPMTs for s in coverage]
    mPMTs = [s*mass for s in cPMTs]

    print "Num pmts", cPMTs
    print "Total Mass",mPMTs

    M_U238,Lambda_U238,Abund_U238 = 3.953e-25,4.916e-18,0.992745
    PPM_U238    = float(arguments["--U238_PPM"])
    ActivityU238= Lambda_U238*PPM_U238/M_U238/1e6
    mPMTsU238 = [s*ActivityU238 for s in mPMTs]
    print 'U238',mPMTsU238, ', PPM:',PPM_U238,'activity per PMT:', ActivityU238*mass,'Bq per PMT per isotope in chain'

    M_Th232,Lambda_Th232,Abund_Th232 = 3.853145e-25, 1.57e-18,1.0
    PPM_Th232    = float(arguments["--Th232_PPM"])
    ActivityTh232= Lambda_Th232*PPM_Th232/M_Th232/1e6
    mPMTsTh232 = [s*ActivityTh232 for s in mPMTs]
    print 'Th232',mPMTsTh232, ', PPM:',PPM_Th232, 'activity per PMT:', ActivityTh232*mass,'Bq per PMT per isotope in chain'

    M_K,Lambda_K,Abund_K = 6.636286e-26,1.842e-18,0.00117
    PPM_K    = float(arguments["--K_PPM"])
    ActivityK= Lambda_K*PPM_K/M_K/1e6*Abund_K
    mPMTsK = [s*ActivityK for s in mPMTs]
    print 'K',mPMTsK, ', PPM:',PPM_K, 'activity per PMT:', ActivityK*mass,'Bq per PMT per isotope in chain'

    print

    return mPMTs,mPMTsU238,mPMTsTh232,mPMTsK



def loadTankActivity():                                         ##added by Leah: activity from steel in tank
    ##MASS OF STEEL USED IN KG -- assuming use of steel grade 304
    density = 8000                                              ##kg/m^3
    r1 = (float(arguments["--tankRadius"]))/1000.               ##outer radius inc steel thick in m
    h1 = 2 * (float(arguments["--halfHeight"]))/1000.           ##outer height inc steel thick in m
    V1 = pi * h1 * r1**2                                        ##outer volume inc steel thick in m^3

    r2 = r1 - (float(arguments["--steelThick"]))/1000.         ##inner radius in m
    h2 = h1 - 2*(float(arguments["--steelThick"]))/1000.       ##inner height in m
    V2 = pi * h2 * r2**2                                        ##inner volume in m^3

    tankvol = V1 - V2                                           ##hollow cylinder m^3
    tankmass = tankvol*density

    print "Total steel mass",tankmass,"kg"

    ##STAINLESS STEEL CONTAINS 60CO, 137CS
    act_60co = 19e-3                                           ##mean value from "Measurements of extremely low radioactivity in stainless steel" paper Bq per kg
    tankact_60co = tankmass * act_60co                          ##activity in Bq
    print "60Co in tank steel:\n activity per kilogram:", act_60co,"Bq/kg,\n total activity:",tankact_60co,"Bq"

    act_137cs = 0.77e-3                                         ##mean value from "Measurements of extremely low radioactivity in stainless steel" paper Bq per kg
    tankact_137cs = tankmass * act_137cs                        ##activity in Bq
    print "137Cs in tank steel:\n activity per kilogram:", act_137cs,"Bq/kg,\n total activity:",tankact_137cs,"Bq"

    return tankmass,tankact_60co,tankact_137cs


def loadConcreteActivity():                                     ##added by Leah: activity from concrete
    ##MASS OF CONCRETE USED IN KG -- assuming normal-weight concrete (NWC)
    density = 2300                                              ##kg/m^3
    thickness = 0.5                                             ## slab thickness in metres - alter for desired value
    concvol = 25.5*(pi*pow(13.,2)-pi*pow(12.5,2)) + 0.5*pi*(pow(13.,2)) #0.5m thick, 25m high concrete 'tube': outer diameter 26m, 
									#inner diameter 25m, plus 0.5m base of diameter 26m
    concmass = concvol * density

    print "Total concrete slab mass",concmass,"kg"

##CONCRETE CONTAINS 238U, 232TH, 40K

    act_238u = 61                                              ##mean UK value from "Natural radioactivity in building materials in the European Union" paper Bq per kg
    concact_238u = concmass * act_238u                         ##activity in Bq
    print "238U in concrete slab:\n activity per kilogram:", act_238u,"Bq/kg,\n total activity:",concact_238u,"Bq"

    act_232th = 30                                              ##mean UK value from "Natural radioactivity in building materials in the European Union" paper Bq per kg
    concact_232th = concmass * act_232th                        ##activity in Bq
    print "232Th in concrete slab:\n activity per kilogram:", act_232th,"Bq/kg,\n total activity:",concact_232th,"Bq"

    act_40k = 493                                               ##mean UK value from "Natural radioactivity in building materials in the European Union" paper Bq per kg
    concact_40k = concmass * act_40k                            ##activity in Bq
    print "40K in concrete slab:\n activity per kilogram:", act_40k,"Bq/kg,\n total activity:",concact_40k,"Bq"

    return concmass,concact_238u,concact_232th,concact_40k
###gravel under concrete? 24in


def loadShotcreteActivity():
    ##MASS OF SHOTCRETE USED IN KG -- assuming same shotcrete as used in snolab (Shotcrete Application in SNOLAB paper)
    density = 2400                                              ##kg/m^3 -- approx from "Shotcrete in Tunnel Construction - Introduction to the basic technology of sprayed concrete"
    thickness = 0.1                                             ##thickness in metres
    r1 = 25                                                     ##outer radius in m
    h1 = 25                                                     ##outer height in m
    V1 = pi * h1 * r1**2                                        ##outer volume in m^3

    r2 = r1 - thickness                                         ##inner radius in m
    h2 = h1 - 2*thickness                                       ##inner height in m
    V2 = pi * h2 * r2**2                                        ##inner volume in m^3

    shotvol = V1 - V2                                           ##hollow cylinder m^3
    shotmass = shotvol*density

    print "Total shotcrete mass",shotmass,"kg"

    mass_238u = 3.953e-5
    lambda_238u = 4.916e-18
    abund_238u = 0.992745
    ppm_238u = 2.6
    act_238u =( (lambda_238u * ppm_238u)/mass_238u/1e6)*shotmass
    print "238U in shotcrete coating:\n ppm:", ppm_238u,"\n total activity:",act_238u,"Bq"

    mass_232th = 3.853145e-25
    lambda_232th = 1.57e-18
    abund_232th = 1.0
    ppm_232th = 14
    act_232th = ((lambda_232th * ppm_232th)/mass_232th/1e6)*shotmass
    print "232Th in shotcrete coating:\n ppm:", ppm_232th,"\n total activity:",act_232th,"Bq"

    mass_40k = 6.636286e-26
    lambda_40k = 1.842e-18
    abund_40k = 0.00117
    ppm_40k = 16000           ##1.6%
    act_40k = ((lambda_40k * ppm_40k)/mass_40k/1e6)*shotmass
    print "40K in shotcrete coating:\n ppm:", ppm_40k,"\n total activity:",act_40k,"Bq"

    return shotmass,act_238u,act_232th,act_40k



def loadRockActivity():
    ##NB MOVED TO ROCK SALT CAVERN RATHER THAN ROCK...
    ##5M THICKNESS
    ##MASS OF SALT IN WALL IN KG -- assuming pure rock salt walls (Overview of the European Underground Facilities paper)
    density = 2165                                              ##kg/m^3 -- approx from "Physical Properties Data for Rock Salt"
    thickness = 5                                               ##thickness in metres
    r1 = 13 + thickness                                         ##outer radius in m
    h1 = 25.5 + 2*thickness                                     ##outer height in m
    V1 = pi * h1 * r1**2                                        ##outer volume in m^3

    r2 = 13                                                     ##inner radius in m
    h2 = 25.5                                                   ##inner height in m
    V2 = pi * h2 * r2**2                                        ##inner volume in m^3

    rockvol = V1 - V2                                           ##hollow cylinder m^3 outside 25m cylindrical cavern 
								##plus 0.5m concrete layer on walls and floor 
    rockmass = rockvol*density

    print "Total relevant rock mass",rockmass,"kg"

    mass_238u = 3.953e-5
    lambda_238u = 4.916e-18
    abund_238u = 0.992745
    ppm_238u = 0.067
    act_238u =( (lambda_238u * ppm_238u)/mass_238u/1e6)*rockmass
    print "238U in shotcrete coating:\n ppm:", ppm_238u,"\n total activity:",act_238u,"Bq"

    mass_232th = 3.853145e-25
    lambda_232th = 1.57e-18
    abund_232th = 1.0
    ppm_232th = 0.125
    act_232th = ((lambda_232th * ppm_232th)/mass_232th/1e6)*rockmass
    print "232Th in shotcrete coating:\n ppm:", ppm_232th,"\n total activity:",act_232th,"Bq"

    mass_40k = 6.636286e-26
    lambda_40k = 1.842e-18
    abund_40k = 0.00117
    ppm_40k = 1130
    act_40k = ((lambda_40k * ppm_40k)/mass_40k/1e6)*rockmass
    print "40K in shotcrete coating:\n ppm:", ppm_40k,"\n total activity:",act_40k,"Bq"

    return rockmass,act_238u,act_232th,act_40k


def loadGdActivity():

    d,process,coverage = loadSimulationParamatersNew()

    GdU238    = float(arguments["--U238_Gd"]) / 1000. * nKiloTons * 1e6 * 0.002 # bq/kg * kg of water * Gd(SO4)3 concentration
    GdTh232   = float(arguments["--Th232_Gd"])/ 1000. * nKiloTons * 1e6 * 0.002 # bq/kg * kg of water * Gd(SO4)3 concentration
    GdU235    = float(arguments["--U235_Gd"]) / 1000. * nKiloTons * 1e6 * 0.002  #bq/kg * kg of water * Gd(SO4)3 concentration


    return GdU238,GdTh232,GdU235





