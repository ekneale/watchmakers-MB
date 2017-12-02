from ROOT import TRandom3
from ROOT import TChain,TGraph,TGraphErrors,gSystem,gROOT,TH1D,TH2D,TFile,TCanvas,TF1
from ROOT import THStack,Double
from ROOT import kRed,kBlue,kGreen,kCyan,kOrange

from ROOT import kOrange as kO,kBlue as kB,kGreen as kG
from ROOT import kMagenta as kM,kAzure as kA,kRed as kR
from ROOT import TCanvas,TLine, TLatex
from numpy import sqrt
from ROOT import gStyle,gPad,TPaletteAxis
import os.path
from stat import S_IRWXG,S_IRWXU
from shutil import rmtree
import warnings

import numpy as np
from numpy import array as npa
from numpy import power,absolute,logical_and,column_stack,zeros,empty,append,\
sqrt,absolute,recarray


from math import pow,exp,log10,pi
gStyle.SetOptStat(1112211)

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
'day','boulby', 1.0, 0.043, 0.133]

docstring = """
    Usage: watchmakers.py [options]

    Arguments:

    Options:
    -D                  Delete all current photocoverage directory.

    -j=<jobType>        Create submision scripts (1,2,4:rat-pac files|case >3 ntuplefiles) [default %d]
                        >3 option will generate a nutple_root_files_flags folder for results
    -m                  Also generate macro files
    -N=<N>              Number of MC script that were run [Default: %d]
    -e=<runBeamEntry>   Number of entries per macro (U/Th event x5) [Default: %d]
    --depth=<depthD>    Depth of detector (for fast neutron spectra) [Default: %f]

    -n                  generate ntuple from single rat-pac root files
    --extractNtup       generate ntuple from all rat-pac root files
    -f=<ifile>          Input file [Default: %s]
    --ft=<ifile2>        Input file type [Default: %s]
    --ntupleout=<outN>  Name of ntuple out [Default: %s]
    -o=<outputfile>     Efficiency output file [Default: %s]
    --supernovaFormat   Record supernova files instead of golden files

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
    -C=<cov>            Pick a single coverage to analyse

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

    --timeScale=<_ts>   Integration period (sec,day,month,year) [Default: %s]
    --site=<_site>      Site of the experiment (boulby,fairport) [Default: %s]
    --OnOff=<_OOratio>  Ratio of reactor on to reactor off [Default: %d]
    --cores=<_cores>    Number of cores to discover [Default: 1]
    --bkgdSys=<_BSys>   Systematic background percentage [Default: 0.2]
    --40MWth            Option to sensitivity to do case scenarios
    --40MWthSig=<_SL>   Sigma discovery [Default: 3.0]

    --U238_PPM=<_Uppm>  Concentration of U-238 in glass [Default: %f]
    --Th232_PPM=<_Thp>  Concentration of Th-232 in glass [Default: %f]
    --Rn222=<_Rn>       Radon activity in water [Default: 6.4]

    --detectMedia=<_dM>  Detector media (doped_water,...)
    --collectionEff=<CE> Collection efficiency (e.g.: 0.85,0.67,0.475)

    --pmtModel=<_PMTM>   PMT Model (r7081pe for 10inch or r11780_hqe for 12inch)
    --photocath =<_PC>  PMT photocathode (R7081HQE)

    """ % (defaultValues[0],defaultValues[1],defaultValues[2],defaultValues[3],defaultValues[4],\
           defaultValues[5],defaultValues[6],defaultValues[7],defaultValues[8],\
           defaultValues[9],defaultValues[10],defaultValues[11],defaultValues[12],\
           defaultValues[13],defaultValues[14],defaultValues[15],defaultValues[16],\
           defaultValues[17],defaultValues[18],defaultValues[19],defaultValues[20],\
           defaultValues[21],defaultValues[22],defaultValues[23])

try:
    import docopt
    arguments = docopt.docopt(docstring)
    print 'using docopt as the user control interface'
except ImportError:
    print 'docopt is not a recognized module, it is required to run this module'


gSystem.Load("$RATROOT/lib/libRATEvent")
gSystem.AddIncludePath(" -I$RATROOT/include")


gROOT.LoadMacro("$WATCHENV/watchmakers/goldenFileExtractor.C")
from ROOT import goldenFileExtractor

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
    d['CHAIN_222Rn_NA'] = ['214Pb','214Bi','210Bi','210Tl']

    # Radioisotope that should have beta-Neutron modes, (beta only generated)
    #    A = ['16','17','18','17','18','8','9','11']
    #    Z = ['6','6','6','7','7','2','3','3']
    #Reduced selection
    A = ['9','11']
    Z = ['3','3']
    ZA = A
    for i in range(len(A)):
        ZA[i] = str(int(A[i])*1000 +int(Z[i]))
    d['A_Z'] =  ZA
    #Oscillated spectrum at Boulby and IMB site
    d['ibd'] = ['boulby','imb']
    d['N'] = ['neutron']
    d['IBD'] = ['IBD']
    # Fast neutron contamination
    d['FN'] = ['QGSP_BERT_EMV','QGSP_BERT_EMX','QGSP_BERT','QGSP_BIC','QBBC',\
    'QBBC_EMZ','FTFP_BERT','QGSP_FTFP_BERT']
    #
    #List of all physics process clumped together
    iso = ['CHAIN_238U_NA','CHAIN_232Th_NA','CHAIN_222Rn_NA','A_Z','ibd',\
    'FN','N','IBD']
    loc = ['PMT','PMT','FV','RN','S','FN','N','I']
    #Photocoverage selected
    coverage = ['10pct','15pct','20pct','25pct','30pct','35pct','40pct']
    coveragePCT = {'10pct':9.86037,'15pct':14.887,'20pct':19.4453,\
    '25pct':24.994,'30pct':28.8925,'35pct':34.3254,'40pct':39.1385}

    return d, iso,loc,coverage,coveragePCT

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

    # print '\nExtracting from log files number of PMTs and photocoverage'
    # print cmd
    # print 'Number of PMTs      ', c
    # print 'Actual photocoverage', d
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

    #Evaluate FV to total detector volume ratio
    nKiloTons   = 3.22
    FreeProtons = 0.6065
    TNU         = FreeProtons* nKiloTons *timeSec
    #FVkTonRatio = pow(float(arguments['--fv']),3)/pow(float(arguments['--tankDis']),3)

    ### This had been changed by M.B. from Tamzin implementation
    fidRadius = float(arguments['--tankRadius'])-float(arguments['--steelThick'])-float(arguments['--shieldThick'])-float(arguments['--fidThick'])
    fidHeight = float(arguments['--halfHeight'])-float(arguments['--steelThick'])-float(arguments['--shieldThick'])-float(arguments['--fidThick'])

    tankRadius  = float(arguments["--tankRadius"])-float(arguments['--steelThick'])
    tankHeight  = float(arguments["--halfHeight"])-float(arguments['--steelThick'])

    fidVolume  = pow(fidRadius,2)*(2.*fidHeight)
    tankVolume = pow(tankRadius,2)*(2.*tankHeight)
    FVkTonRatio = fidVolume/tankVolume
    #print "Change in load.py"

    #Fast neutrons conversion
    #Rock mass
    # Original Estimate
    # volumeR         = (2.*22.5*23.8*1.0+2.*17*23.8*1.0+2.*22.5*17.*1.0)
    volumeR         = power(2*float(arguments["--tankRadius"]/1000.)+6,2)*(2*float(arguments["--halfHeight"]/1000.)+6) - power(2*float(arguments["--tankRadius"]/1000.)+4,2)*(2*float(arguments["--halfHeight"]/1000.)+4) # Rock cavern e.g. (22m x 22m x 22m) - (20m x 20m x 20m)
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


    boulbyIBDRate   = 1120.8*.4/.6 *TNU #//924.48*TNU Taken from website, average corrected
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
    N_Rn222     = 2e-3 # Bq/m3
    ActivityRn222     = float(arguments["--Rn222"])
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
    arr[:]      = 6.4
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



    #    # Read in the different radionuclide
    #    proc        +=  ['16006','17006','18006','17007','18007','8002','9003',\
    #    '11003']
    #    loca        +=  ['RN','RN','RN','RN','RN','RN','RN','RN']
    #    acc         +=  ['di','di','di','di','di','di','di','di']
    #    #normalised to 9Li from SK
    #    arr         = npa([0.02,0.001,0.001,0.59*0.002,4e-6,0.23,1.9,0.01])/1.9
    #    arr         *= radionuclideRate[0]
    #    Activity    = append(Activity,arr)
    #    br         +=  [.988,1.0,1.0,0.951,0.143,0.16,0.495,0.927]
    #    site        += ['','','','','','','','']
    #

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

    #    # Read in the different radionuclide
    #    proc        +=  ['16006','17006','18006','17007','18007','8002','9003',\
    #    '11003']
    #    loca        +=  ['RN','RN','RN','RN','RN','RN','RN','RN']
    #    acc         +=  ['di','di','di','di','di','di','di','di']
    #    #normalised to 9Li from SK
    #    arr         = npa([ 0.02,0.001,0.001,0.59*0.002,4e-6,0.23,1.9,0.01])/1.9
    #    arr         *= radionuclideRate[1]
    #    Activity    = append(Activity,arr)
    #    br         +=  [.988,1.0,1.0,0.951,0.143,0.16,0.495,0.927]
    #    site        += ['boulby','boulby','boulby','boulby','boulby','boulby',\
    #    'boulby','boulby']



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

    #
    # return inta,proc,loca,acc,arr,Activity,br,site,timeS,\
    # boulbyIBDRate*FVkTonRatio,mass,dAct,coveNumber,covePCT
