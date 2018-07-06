#import watchmakers as PR

from watchmakers.load import *
from watchmakers.analysis import *
from io_operations import testEnabledCondition
if arguments['--noRoot']:
    print 'Not loading the neutrinoOscillation module'
else:
    import watchmakers.NeutrinoOscillation as nuOsc

from decimal import *
setcontext(ExtendedContext)

from numpy import max
t = arguments['--timeScale']

fidRadius = float(arguments['--tankRadius'])-float(arguments['--steelThick'])-float(arguments['--shieldThick'])-float(arguments['--fidThick'])
fidHeight = float(arguments['--halfHeight'])-float(arguments['--steelThick'])-float(arguments['--shieldThick'])-float(arguments['--fidThick'])

pmtRadius = float(arguments['--tankRadius'])-float(arguments['--steelThick'])-float(arguments['--shieldThick'])
pmtHeight = float(arguments['--halfHeight'])-float(arguments['--steelThick'])-float(arguments['--shieldThick'])

detectorRadius  = float(arguments['--tankRadius'])-float(arguments['--steelThick'])
detectorHeight  = float(arguments['--halfHeight'])-float(arguments['--steelThick'])

def drange(start, stop, step):
    rii= start
    while rii<stop:
        yield rii
        rii+=step


def accidentalCanvas(Graphs):
    # Not currently in use
    proc    = []
    loca    = []
    type    = []
    color   = []
    lineS   = []
    acc     = []

    proc    += ['234Pa','214Pb','214Bi','210Bi','210Tl']
    loca    += ['PMT','PMT','PMT','PMT','PMT']
    type    += ['si','si','si','si','si']
    color   += [kO+0,  kO+1, kO+2, kO+3, kO-2]

    proc    +=['232Th','228Ac','212Pb','212Bi','208Tl']
    loca    +=['PMT','PMT','PMT','PMT','PMT']
    type    += ['si','si','si','si','si']
    color   += [kM+0,  kM-6, kM+2, kM-3, kM-7]

    proc    +=['214Pb','214Bi','210Bi','210Tl']
    loca    +=['FV','FV','FV','FV']
    type    += ['ei','ei','ei','ei']
    color   += [kB+0,  kB+1, kB+2, kB+3]

    proc    +=['214Pb','214Bi','210Bi','210Tl']
    loca    +=['FV','FV','FV','FV']
    type    += ['si','si','si','si']
    color   += [kG+0,  kG+1, kG+2, kG+3]

    proc    +=['QGSP_BERT_EMV','QGSP_BERT_EMX','QGSP_BERT','QGSP_BIC',\
    'QBBC','QBBC_EMZ','FTFP_BERT','QGSP_FTFP_BERT']
    loca    +=['FN','FN','FN','FN','FN','FN','FN','FN']
    type    +=['si','si','si','si','si','si','si','si']
    color   += [kO+2,kO+4,kO+2,kO+4,kO+2,kO+4,kO+2,kO+4]


def siteCanvas(Graphs,site,cut,hist):
    # Historical, not currently in use
    proc        = []
    loca        = []
    type        = []
    color       = []
    lineS       = []
    acc         = []
    scale       = []

    #fast neutrons
    proc        += ['QGSP_BERT_EMV','QGSP_BERT_EMX','QGSP_BERT','QGSP_BIC',\
    'QBBC','QBBC_EMZ','FTFP_BERT','QGSP_FTFP_BERT']
    _t          = 'FN%s' % (site)
    loca        += [_t,_t, _t,_t,_t,_t,_t,_t]
    type        += ['si','si','si','si','si','si','si','si']
    acc         += ['corr','corr','corr','corr','corr','corr','corr', 'corr']
    color       += [kO+6,kO+6,kO+6,kO+6,kO+6,kO+6,kO+6,kO+6]
    lineS       += [2,2,2,2,2,2,2,2]
    scale       += [1./8.,1./8.,1./8.,1./8.,1./8.,1./8.,1./8.,1./8.]

    #ibds
    if site == 'boulby':
        proc    += ['boulby','boulby','neutron']
    else:
        proc    += ['imb','imb','neutron']
    loca        += ['S','S','N%s'%(site)]
    type        += ['ei','ei','ei']
    acc         += ['di','corr','corr']
    color       += [kA+0,kA-0,kA-0]
    lineS       += [1,2,2]
    scale       += [-1.0,0.0,0.0]

    proc        += ['9003','11003']
    _t          =  'RN%s' % (site)
    loca        += [_t,_t]
    type        += ['ei','ei']
    acc         += ['di','di']
    color       += [kG+3,kG+3]


    lineS       += [1,1,1,1,1,1,1,1]
    scale       += [1.,1.,1.,1.,1.,1.,1.,1.]


    for PC in range(10,41):
        S = 0.
        B = Graphs['all_cut_accidental'].Eval(PC)

        for i,arr in enumerate(proc):
            _str = 'scaled_%s_%s_%s_1_abs_%s' %(type[i],proc[i],loca[i],acc[i])
            if scale[i] > 0.:
                B+=Graphs[_str].Eval(PC)*scale[i]
            elif scale[i] < 0:
                S+=Graphs[_str].Eval(PC)
                #        print cut,S,B,S/B,S/sqrt(B+S)
        hist.Fill(PC,cut,S/sqrt(B+S))

    C1 = TCanvas('C1','%s'%(t),1200,800)

    Graphs['all_cut_accidental'].Draw('AC')
    Graphs['all_cut_accidental'].SetLineColor(2)
    Graphs['all_cut_accidental'].SetLineWidth(3)
    Graphs['all_cut_accidental'].SetTitle('doubles rate (all cuts)')
    Graphs['all_cut_accidental'].GetYaxis().SetTitle('event rate [%s^{-1}]'%(t))
    Graphs['all_cut_accidental'].GetXaxis().SetTitle('photo-coverage [%]')
    if site == 'boulby':
        Graphs['all_cut_accidental'].GetYaxis().SetRangeUser(0,70/31.)
    else:
        Graphs['all_cut_accidental'].GetYaxis().SetRangeUser(0,400/31.)

    for i,arr in enumerate(proc):
        _str = 'scaled_%s_%s_%s_1_abs_%s' %(type[i],proc[i],loca[i],acc[i])
        Graphs[_str].Draw('same')
        Graphs[_str].SetLineColor(color[i])
        Graphs[_str].SetLineStyle(lineS[i])
        if loca[i] == 'PMT':
            _descrpt  = '%s %s glass'%(proc[i],loca[i])
            Graphs[_str].SetTitle(_descrpt)
        elif loca[i] == 'FV':
            _descrpt  = '%s %s in water (%s)'%(proc[i],loca[i],type[i])
            Graphs[_str].SetTitle(_descrpt)
        else:
            if proc[i]=='imb':
                Graphs[_str].SetTitle('Prompt-positron at Perry')
            if proc[i]=='boulby':
                Graphs[_str].SetTitle('Prompt-positron at Boulby')
        Graphs[_str].SetMarkerStyle(1)

    t1 = TLatex()
    if site == 'boulby':
        line = TLine(7, 4.,42, 4)
        t1.DrawLatex( 10,65., 'pe > %d.0'%(cut))
    else:
        line = TLine(7, 40,42, 40)
        t1.DrawLatex( 10,363., 'pe > %d.0'%(cut))

    line.SetLineColor(4)
    line.Draw('same')
    C1.SetGridy()
    C1.Update()

    #    C1.SaveAs('gif/pc_%s_%d.gif'%(site,cut))
    C1.SaveAs('gif/pc_%s.gif+'%(site))
    return hist



def sensitivityMapNew():

    # Function to find the optimal signal to background as a function

    detectorMedium      = 1
    detectorMass        = 1000.
    reactorPower        = 0.04
    reactorStandoff     = 25.0
    experiment = nuOsc.NeutrinoOscillation(detectorMedium,detectorMass,reactorPower,reactorStandoff)
    preOsc,afterOsc = experiment.FindRate()
    print preOsc,afterOsc
    site = arguments["--site"]
    if site != 'boulby':
        site = ''
    # Need to fix this for future running

    OnOffRatio = float(arguments["--OnOff"])
    print site,OnOffRatio

    cores = int(arguments["--cores"])

    if arguments["--RNRedux"]:
        rnRedux = float(arguments["--RNRedux"])
        if rnRedux>1:
            print "Value of reduction of radionuclide greater than 1, setting to 0"
            rnRedux = 0.0
    else:
        rnRedux = 0.0

    if t == 'sec':
        timeAdjustment = 24*3600.

    if t == 'day':
        timeAdjustment = 1.0

    if t == 'month':
        timeAdjustment = 1./31.

    if t == 'year':
        timeAdjustment = 1./365.

    maxTime = 14400.*timeAdjustment

    proc        = []
    loca        = []
    type        = []
    color       = []
    lineS       = []
    acc         = []
    scale       = []

    parameters  = loadAnalysisParameters(t)
    rates       = parameters[11]
    FVkTonRatio = (pow(fidRadius,2)*fidHeight)/(pow(detectorRadius,2)*detectorHeight)
    boulbyRate,imbRate = rates["boulby_S"]*FVkTonRatio,rates["imb_S"]*FVkTonRatio
    print 'rates:',imbRate,boulbyRate, ' per ', t, '(',rates,FVkTonRatio,')'

    #fast neutrons
    proc        += ['QGSP_BERT_EMV','QGSP_BERT_EMX','QGSP_BERT','QGSP_BIC',\
    'QBBC','QBBC_EMZ','FTFP_BERT','QGSP_FTFP_BERT']
    _t          = 'FN%s' % (site)
    loca        += [_t,_t, _t,_t,_t,_t,_t,_t]
    type        += ['si','si','si','si','si','si','si','si']
    acc         += ['corr','corr','corr','corr','corr','corr','corr', 'corr']
    color       += [kO+6,kO+6,kO+6,kO+6,kO+6,kO+6,kO+6,kO+6]
    lineS       += [2,2,2,2,2,2,2,2]
    scale       += [1./8.,1./8.,1./8.,1./8.,1./8.,1./8.,1./8.,1./8.]
    #radionuclides
    proc        += ['9003','11003']
    _t          =  'RN%s' % (site)
    loca        += [_t,_t]
    type        += ['ei','ei']
    acc         += ['di','di']
    color       += [kG+3,kG+3]
    lineS       += [1,1]
    scale       += [1.,1.]
    #ibds
    if site == 'boulby':
        proc    += ['boulby','boulby','neutron']
    else:
        proc    += ['imb','imb','neutron']
    loca        += ['S','S','N%s'%(site)]
    type        += ['ei','ei','ei']
    acc         += ['di','corr','corr']
    color       += [kA+0,kA-0,kA-0]
    lineS       += [1,2,2]
    scale       += [-1.0,0.0,0.0]

    c1 = TCanvas('c1','c1',1618,1000)
    c1.SetRightMargin(0.53)
    c1.Divide(2,2)
    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)

    if additionalString == "":
        additionalString = "_default"

    if site == 'boulby':
        location = 'Boulby '
    else:
        location = 'Fairport '

    hist = TH2D('hist','3#sigma discovery phase space -  %s '%(location),31,9.5,40.5,30,9.5,39.5)
    hist.SetXTitle('photocoverage [%]')
    hist.SetYTitle('photoelectron threhsold cut [p.e.]')
    hist.SetZTitle('off-time [%s] for 3 #sigma discovery'%(t))
    hist.GetZaxis().SetTitleOffset(-.55);
    hist.GetZaxis().SetTitleColor(1);
    hist.GetZaxis().CenterTitle();

    gStyle.SetOptStat(0)
    gStyle.SetPalette(55)

    fileIN = 'processed_watchman_March%s.root' %(additionalString)
    fileIN = 'ntuple_root_files%s/processed_watchman.root' %(additionalString)

    print 'Will read in', fileIN
    print "%s %s %10s %10s %10s %10s %10s %10s %10s %10s %10s" %('PC','cut','S','eff','BAC','BRN','BFN','BRS','B_tot','SSBB','T3SIGMA')
    _graph = Graph()
    _graph.SetName('_graph')
    _graph.SetTitle('off-time [%s] required till a 3#sigma discovery'%(t))

    _graphS = Graph()
    _graphS.SetName('_graphS')
    _graphS.SetTitle('IBDs from %d-core at full power'%(cores))

    _graphB = Graph()
    _graphB.SetName('_graphB')
    _graphB.SetTitle('total background rate after optimization')

    _graphBFN = Graph()
    _graphBFN.SetLineColor(kO+6)
    _graphBFN.SetName('_graphBFN')
    _graphBFN.SetTitle('total fast neutron component')

    _graphBRN = Graph()
    _graphBRN.SetLineColor(kG+3)
    _graphBRN.SetName('_graphBFN')
    _graphBRN.SetTitle('total radionuclide component')

    _graphBAC = Graph()
    _graphBAC.SetLineColor(2)
    _graphBAC.SetName('_graphBFN')
    _graphBAC.SetTitle('total accidental component')

    _graphBRS = Graph()
    _graphBRS.SetLineColor(4)
    _graphBRS.SetName('_graphBRS')
    _graphBRS.SetTitle('other reactor antineutrinos')

    _graphSIZE = Graph()
    _graphSIZE.SetLineColor(4)
    _graphSIZE.SetName('_graphBRS')
    _graphSIZE.SetTitle('3-month discovery of 40 MWth at 25 km')

    graphDict = {}
    for pcVal in range(10,41,5):
        graphDict["%dpct"%(pcVal)] = Graph()
        graphDict["%dpct"%(pcVal)].SetLineColor(4)
        graphDict["%dpct"%(pcVal)].SetName('_graph%dpct'%(pcVal))
        graphDict["%dpct"%(pcVal)].SetTitle('40 MWth standoff Mass relation - signal and background scaled from %dpct values'%(pcVal))
        graphDict["%dpct_sys10"%(pcVal)] = Graph()
        graphDict["%dpct_sys10"%(pcVal)].SetLineColor(4)
        graphDict["%dpct_sys10"%(pcVal)].SetName('_graph%dpct_sys10'%(pcVal))
        graphDict["%dpct_sys10"%(pcVal)].SetTitle('%dpct values - sys 10 percent'%(pcVal))
        graphDict["%dpct_sys20"%(pcVal)] = Graph()
        graphDict["%dpct_sys20"%(pcVal)].SetLineColor(4)
        graphDict["%dpct_sys20"%(pcVal)].SetName('_graph%dpct_sys20'%(pcVal))
        graphDict["%dpct_sys20"%(pcVal)].SetTitle('%dpct values - sys 20 percent'%(pcVal))


    _graph.SetLineWidth(2)
    _graphS.SetLineWidth(2)
    _graphS.SetLineColor(4)
    _graphB.SetLineWidth(2)
    _graphBAC.SetLineWidth(2)
    _graphBFN.SetLineWidth(2)
    _graphBRN.SetLineWidth(2)
    _graphBRS.SetLineWidth(2)
    _graphBAC.SetLineStyle(2)
    _graphBFN.SetLineStyle(2)
    _graphBRN.SetLineStyle(2)
    _graphBRS.SetLineStyle(2)

    minY    =   0

    gCnt = 0
    for PC in drange(10.,41.,5.0):
        #        print 'processing photocoverage ',PC,'% :',
        S_tmp           = 0.
        B_tmp           = 0.0
        BFN_tmp           = 0.0
        BRN_tmp           = 0.0
        BAC_tmp           = 0.0
        BRS_tmp           = 0.0
        SoverS_B        = 0.
        SoverS_B_sys20  = 0.
        T3SIGMA_MIN     = 16000000
        cut_tmp     = 0.

        c1.cd(1)
        hist.Draw('colz')
        gPad.SetTicks()

        c1.cd(2)
        _graph.Draw('AL')
        _graph.GetXaxis().SetTitle('photocoverage [%]')
        _graph.GetYaxis().SetTitle('time [%s]'%(t))
        gPad.SetTicks()
        gPad.SetGridy()

        c1.cd(3)
        _graphS.Draw('AL')
        y = _graphS.GetYaxis().GetXmax()
        _graphS.GetYaxis().SetRangeUser(0.,_graphS.GetYaxis().GetXmax()*2.)
        _graphS.GetXaxis().SetTitle('photocoverage [%]')
        _graphS.GetYaxis().SetTitle('event rate [%s^{-1}]'%(t))
        c1.cd(3).BuildLegend(0.42, 0.84, 0.85, 0.88)
        gPad.SetGridy()
        gPad.SetTicks()

        c1.cd(4)
        _graphB.Draw('AL')
        _graphB.GetYaxis().SetRangeUser(0., _graphB.GetYaxis().GetXmax()*2.)
        _graphB.GetXaxis().SetTitle('photocoverage [%]')
        _graphB.GetYaxis().SetTitle('event rate [%s^{-1}]'%(t))
        _graphBAC.Draw('same')
        _graphBFN.Draw('same')
        _graphBRN.Draw('same')
        c1.cd(4).BuildLegend(0.42, 0.72, 0.85, 0.88)
        gPad.SetGridy()
        gPad.SetTicks()

        c1.Update()

        for cut in range(4,40):
            Graphs =obtainAbsoluteEfficiency(fileIN,timeScale=t,cut=cut)

            S       = 0.  # Signal
            BFN     = 0.  # Fast Neutron backgrounds
            BRN     = 0.  # Radionuclides backgrounds
            BAC     = 0.  # Accidental backgrounds
            BRS     = 0.  # Reactor Signal backgrounds
            _str    = 'all_cut_accidental'
            BAC = Graphs[_str].Eval(PC)
            for i,arr in enumerate(proc):
                _str = 'scaled_%s_%s_%s_1_abs_%s' %(type[i],proc[i],loca[i],acc[i])
                if scale[i] > 0.:
                    if loca[i] == 'FN%s'%(site):
			#print BFN, _str
			#print Graphs
                        BFN += Graphs[_str].Eval(PC)*scale[i]
                        # print 'BFN',BFN,Graphs[_str].Eval(PC),scale[i]
                    elif loca[i] == 'RN%s'%(site):
                        BRN += Graphs[_str].Eval(PC)*scale[i]*(1.0-rnRedux)
                        # print 'BRN',BRN,Graphs[_str].Eval(PC),scale[i]
                elif scale[i] < 0:

                    S =Graphs[_str].Eval(PC)

                    #Include 29% other reactor
                    if site == 'boulby' and cores ==1 :
                        BRS    = S*1.29 # 29.% other reactor
                        B      = BAC+BRN+BFN+BRS
                        SSBB   = S/sqrt(B + (S+B)/OnOffRatio)
                        T3SIGMA = 9.*(B +(S+B)/OnOffRatio)/S/S
                        nuRate = boulbyRate
                    elif site == 'boulby' and cores ==2 :
                        S       *= 2.                 #2-cores
                        BRS     = S * 0.29
                        B       = BAC+BRN+BFN+BRS
                        SSBB    = S/sqrt(S+B)
                        T3SIGMA = 9.*(S+B)/S/S #1.29% other reactor
                        nuRate = boulbyRate

                    else:
                        BRS     = S * 0.05
                        B       = BAC+BRN+BFN+BRS
                        SSBB   = S/sqrt(B + (S+B)/OnOffRatio)
                        T3SIGMA = 9.*(B + (S+B)/OnOffRatio)/S/S#5% other reactor
                        nuRate = imbRate
                    # For debugging purposes
                    # print "%3d %3d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f" %(PC,cut,S,S/nuRate,BAC,BRN,BFN,BRS,B,SSBB,T3SIGMA)

                    if  T3SIGMA < maxTime:
                        print "%3d %3d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f" %(PC,cut,S,S/nuRate,BAC,BRN,BFN,BRS,B,SSBB,T3SIGMA)
                        # print "%3d %3d %10.4e %10.4f %10.4e %10.4e %10.4e %10.4e %10.4e" %(PC,cut,S,S/nuRate,BAC,BRN,BFN,BRS,B)
                        hist.Fill(PC,cut,T3SIGMA)
                    if  T3SIGMA >= maxTime:
                        hist.Fill(PC,cut,maxTime)
                    if T3SIGMA <   T3SIGMA_MIN:
                        S_tmp         = S
                        B_tmp         = B
                        BFN_tmp       = BFN
                        BRN_tmp       = BRN
                        BAC_tmp       = BAC
                        BRS_tmp       = BRS
                        T3SIGMA_MIN   = T3SIGMA
                        _signal_Eff   = S/nuRate

        _graph.SetPoint(gCnt,PC,T3SIGMA_MIN)
        _graphS.SetPoint(gCnt,PC,S_tmp)
        _graphB.SetPoint(gCnt,PC,B_tmp)
        _graphBFN.SetPoint(gCnt,PC,BFN_tmp)
        _graphBRN.SetPoint(gCnt,PC,BRN_tmp)
        _graphBAC.SetPoint(gCnt,PC,BAC_tmp)
        _graphBRS.SetPoint(gCnt,PC,BRS_tmp)

        # _sOff,_T3sOff = findEquivalentStandoff(S_tmp,BAC_tmp+BFN_tmp,BRN_tmp,_signal_Eff)

        if arguments["--40MWth"]:
            _size,_T3 = findScaledVolume(S_tmp,BAC_tmp+BFN_tmp,BRN_tmp,1.0)
            graphCnt = 0
            sigmaLvl = float(arguments["--40MWthSig"])
            if 1==1:#PC == 40:
                for _detectSize in drange(1.,50.,1.):
                    _sOff,_sRate,_sRatePostEff = findEquivalentStandoff(S_tmp,\
                    BAC_tmp+BFN_tmp,BRN_tmp,_signal_Eff,0.,6,_detectSize,sigmaLvl)
                    graphDict["%dpct"%(PC)].SetPoint(graphCnt,_sOff,_detectSize)
                    #Do the 10% systematics
                    _sOff,_sRate,_sRatePostEff = findEquivalentStandoff(S_tmp,\
                    BAC_tmp+BFN_tmp,BRN_tmp,_signal_Eff,0.05,6,_detectSize,sigmaLvl)
                    graphDict["%dpct_sys10"%(PC)].SetPoint(graphCnt,_sOff,_detectSize)
                    #Do the 20% systematics
                    _sOff,_sRate,_sRatePostEff = findEquivalentStandoff(S_tmp,\
                    BAC_tmp+BFN_tmp,BRN_tmp,_signal_Eff,0.1,6,_detectSize,sigmaLvl)
                    graphDict["%dpct_sys20"%(PC)].SetPoint(graphCnt,_sOff,_detectSize)
                    graphCnt+=1


            _graphSIZE.SetPoint(gCnt,PC,_size)
        # print PC,_size,_T3,_sOff,_T3sOff

        gCnt+=1

    c1.cd(1)
    hist.Draw('colz')
    gPad.SetTicks()

    c1.cd(2)
    _graph.Draw('AL')
    _graph.GetXaxis().SetTitle('photocoverage [%]')
    _graph.GetYaxis().SetTitle('time [%s]'%(t))
    gPad.SetTicks()
    gPad.SetGridy()


    c1.cd(3)
    _graphS.Draw('AL')
    y = _graphS.GetYaxis().GetXmax()
    _graphS.GetYaxis().SetRangeUser(0.,_graphS.GetYaxis().GetXmax()*2.)
    _graphS.GetXaxis().SetTitle('photocoverage [%]')
    _graphS.GetYaxis().SetTitle('event rate [%s^{-1}]'%(t))
    c1.cd(3).BuildLegend(0.42, 0.84, 0.85, 0.88)
    gPad.SetGridy()
    gPad.SetTicks()

    c1.cd(4)
    _graphS.Draw('AL')
    if site == 'boulby':
        _graphS.GetYaxis().SetRangeUser(0., _graphS.GetYaxis().GetXmax()*3.)
    else:
        _graphS.GetYaxis().SetRangeUser(0., _graphS.GetYaxis().GetXmax()*2.)
    _graphS.GetXaxis().SetTitle('photocoverage [%]')
    _graphS.GetYaxis().SetTitle('event rate [%s^{-1}]'%(t))
    _graphB.Draw('same')
    _graphBAC.Draw('same')
    _graphBFN.Draw('same')
    _graphBRN.Draw('same')
    _graphBRS.Draw('same')
    c1.cd(4).BuildLegend(0.42, 0.72, 0.85, 0.88)
    gPad.SetGridy()
    gPad.SetTicks()

    c3 = TCanvas('c3','c3',1618,1000)
    c3.Divide(1,1)
    c3.cd(1)
    _graphS.Draw('AL')
    if site == 'boulby':
        _graphS.GetYaxis().SetRangeUser(0., _graphS.GetYaxis().GetXmax()*3.)
    else:
        _graphS.GetYaxis().SetRangeUser(0., _graphS.GetYaxis().GetXmax()*2.)
    _graphS.GetXaxis().SetTitle('photocoverage [%]')
    _graphS.GetYaxis().SetTitle('event rate [%s^{-1}]'%(t))
    _graphB.Draw('same')
    _graphBAC.Draw('same')
    _graphBFN.Draw('same')
    _graphBRN.Draw('same')
    _graphBRS.Draw('same')
    c3.cd(1).BuildLegend(0.42, 0.72, 0.85, 0.88)
    gPad.SetGridy()
    gPad.SetTicks()
    c3.Update()


    c2 = TCanvas('c2','c2',1618,1000)
    # c2.SetRightMargin(0.53)
    c2.Divide(1,1)
    c2.cd(1)
    _graph.Draw('AL')
    _graph.GetXaxis().SetTitle('photocoverage [%]')
    _graph.GetYaxis().SetTitle('time [%s]'%(t))
    gPad.SetTicks()
    gPad.SetGridy()

    c4 = TCanvas('c4','c4',1618,1000)
    c4.Divide(1,1)
    c4.cd(1)
    _graphSIZE.Draw('AL')
    _graphSIZE.GetXaxis().SetTitle('photocoverage [%]')
    _graphSIZE.GetYaxis().SetTitle('Detector FV size [kton]')
    c4.Update()
    if arguments["--40MWth"]:
        _pc = ['10','15','20','25','30','35','40']
        _pc = ['40']
        _Canvas = {}
        for pc in _pc:
            _Canvas["c%s"%(pc)] = TCanvas("c%s"%(pc),"c%s"%(pc),1618,1000)
            _Canvas["c%s"%(pc)].Divide(1,1)
            graphDict["%spct"%(pc)].Draw('ALP')
            graphDict["%spct_sys10"%(pc)].Draw('LPsame')
            graphDict["%spct_sys20"%(pc)].Draw('LPsame')
            _Canvas["c%s"%(pc)].Update()
            _Canvas["c%s"%(pc)].SaveAs('results/Canvas%s.gif'%(pc))
            _Canvas["c%s"%(pc)].SaveAs('results/Canvas%s.C'%(pc))


    analysisString =""
    analysisString += arguments["--site"]
    analysisString += arguments["--cores"]
    analysisString += arguments["--RNRedux"]
    analysisString += arguments["--timeScale"]
    analysisString += arguments["--OnOff"]

    hist.SaveAs('results/SensitivityMap%s%s.C'%(additionalString,analysisString))
    _graph.SaveAs('results/SensitivityGraph%s%s.C'%(additionalString,analysisString))
    _graphS.SaveAs('results/SensitivityGraphS%s%s.C'%(additionalString,analysisString))
    _graphB.SaveAs('results/SensitivityGraphB%s%s.C'%(additionalString,analysisString))

    c1.SaveAs('results/SensitivityCanvas%s%s.C'%(additionalString,analysisString))
    c1.SaveAs('results/SensitivityCanvas%s%s.gif'%(additionalString,analysisString))

    c2.SaveAs('results/SensitivityTime%s%s.C'%(additionalString,analysisString))
    c2.SaveAs('results/SensitivityTime%s%s.gif'%(additionalString,analysisString))

    c3.SaveAs('results/SensitivityRate%s%s.C'%(additionalString,analysisString))
    c3.SaveAs('results/SensitivityRate%s%s.gif'%(additionalString,analysisString))

    c4.SaveAs('results/SensitivitySize%s%s.C'%(additionalString,analysisString))
    c4.SaveAs('results/SensitivitySize%s%s.gif'%(additionalString,analysisString))




def sensitivityMapPass2():

    # Function to find the optimal signal to background as a function



    site = arguments["--site"]

    # Need to fix this for future running

    OnOffRatio = float(arguments["--OnOff"])
    print site,'with on-off ratio of ',OnOffRatio

    cores = int(arguments["--cores"])

    if arguments["--RNRedux"]:
        rnRedux = float(arguments["--RNRedux"])
        if rnRedux>1:
            print "Value of reduction of radionuclide greater than 1, setting to 0"
            rnRedux = 0.0
    else:
        rnRedux = 0.0

    if t == 'sec':
        timeAdjustment = 24*3600.
    if t == 'day':
        timeAdjustment = 1.0
    if t == 'month':
        timeAdjustment = 1./31.
    if t == 'year':
        timeAdjustment = 1./365.
    maxTime = 14400.*timeAdjustment

    print '\nEvaluation based on geoneutrinos.org'
    parameters  = loadAnalysisParameters(t)
    rates       = parameters[11]
    sizeDetc    = 2.*pi*pow(fidRadius/1000.,2)*fidHeight/1000./1000.
    sizeTank    = 2.*pi*pow(detectorRadius/1000.,2)*detectorHeight/1000./1000.
    FVkTonRatio = (pow(fidRadius,2)*fidHeight)/(pow(detectorRadius,2)*detectorHeight)
    boulbyRate,imbRate = rates["boulby_S"]*FVkTonRatio,rates["imb_S"]*FVkTonRatio
    print ' boulby rates: %4.2e per %s per %4.2f kton; [r: %4.2f m; z: %4.2f m]'\
    %(boulbyRate,t,sizeDetc,fidRadius/1000.,fidHeight/1000.)
    #fast neutrons
    # print 'Debug',rates["boulby_S"],FVkTonRatio
    print '\nEvaluation not based on geoneutrinos.org'
    detectorMedium,detectorMass,reactorPower,reactorStandoff = 1,sizeDetc*1000.,1.5,24.98
    experiment = nuOsc.NeutrinoOscillation(detectorMedium,detectorMass,reactorPower,reactorStandoff)
    preOsc,afterOsc = experiment.FindRate()
    print ' Neutrino rate pre osc: %4.2f; neutrino rate post osc: %4.2f at %4.2f GWth, at %4.2f km, for %4.2f kton' %(preOsc,afterOsc,reactorPower,reactorStandoff,detectorMass/1000.)
    detectorMedium,detectorMass,reactorPower,reactorStandoff = 1,sizeDetc*1000.,1.575,24.98
    experiment = nuOsc.NeutrinoOscillation(detectorMedium,detectorMass,reactorPower,reactorStandoff)
    preOsc,afterOsc = experiment.FindRate()
    print ' Neutrino rate pre osc: %4.2f; neutrino rate post osc: %4.2f at %4.2f GWth, at %4.2f km, for %4.2f kton' %(preOsc,afterOsc,reactorPower,reactorStandoff,detectorMass/1000.)
    print ''

    proc,loca,type,color,lineS,acc,scale   = [],[],[],[],[],[],[]

    proc        += ['QGSP_BERT_EMV','QGSP_BERT_EMX','QGSP_BERT','QGSP_BIC',\
    'QBBC','QBBC_EMZ','FTFP_BERT','QGSP_FTFP_BERT']
    _t          = 'FN%s' % (site)
    loca        += [_t,_t, _t,_t,_t,_t,_t,_t]
    type        += ['si','si','si','si','si','si','si','si']
    acc         += ['corr','corr','corr','corr','corr','corr','corr', 'corr']
    color       += [kO+6,kO+6,kO+6,kO+6,kO+6,kO+6,kO+6,kO+6]
    lineS       += [2,2,2,2,2,2,2,2]
    scale       += [1./8.,1./8.,1./8.,1./8.,1./8.,1./8.,1./8.,1./8.]

    #radionuclides
    proc        += ['9003','11003']
    _t          =  'RN%s' % (site)
    loca        += [_t,_t]
    type        += ['ei','ei']
    acc         += ['di','di']
    color       += [kG+3,kG+3]
    lineS       += [1,1]
    scale       += [1.,1.]
    #ibds
    if site == 'boulby':
        proc    += ['boulby','boulby','neutron']
    else:
        proc    += ['imb','imb','neutron']
    loca        += ['S','S','N%s'%(site)]
    type        += ['ei','ei','ei']
    acc         += ['di','corr','corr']
    color       += [kA+0,kA-0,kA-0]
    lineS       += [1,2,2]
    scale       += [-1.0,0.0,0.0]

    c1 = TCanvas('c1','c1',1618,1000)
    c1.SetRightMargin(0.53)
    c1.Divide(2,2)

    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)
    if additionalString == "":
        additionalString = "_default"
    location = 'Boulby '

    hist = TH2D('hist','3#sigma discovery phase space -  %s '%(location),31,9.5,40.5,30,9.5,39.5)
    hist.SetXTitle('photocoverage [%]')
    hist.SetYTitle('photoelectron threhsold cut [p.e.]')
    hist.SetZTitle('off-time [%s] for 3 #sigma discovery'%(t))
    hist.GetZaxis().SetTitleOffset(-.55);
    hist.GetZaxis().SetTitleColor(1);
    hist.GetZaxis().CenterTitle();
    gStyle.SetOptStat(0)
    gStyle.SetPalette(55)



    fileIN = 'pass2_root_files%s/processed_watchman.root' %(additionalString)
    h = {}

    binR,rangeRmin,rangeRmax = 31,0.45,3.55
    binwidthR = (rangeRmax-rangeRmin)/binR
    binN,rangeNmin,rangeNmax = 48,7.5,55.5
    binwidthN = (rangeNmax-rangeNmin)/binN

    _cov = arguments['-C']

    proc = ['214Bi_PMT','208Tl_PMT','210Tl_PMT']
    proc = ['214Bi_PMT','208Tl_PMT']
    proc = ['214Bi_WV','208Tl_PMT','214Bi_PMT']
    _proc = 'Sum'
    h['hist%s'%(_proc)] = TH2D('hist%s'%(_proc),'%s Rate of events -  %s '%(_proc,location),binR,rangeRmin,rangeRmax,binN,rangeNmin,rangeNmax)
    h['hist%s'%(_proc)].SetXTitle('distance from wall [m]')
    h['hist%s'%(_proc)].SetYTitle('n9 cut')
    h['hist%s'%(_proc)].SetZTitle('rate per %s'%(t))
    h['hist%s'%(_proc)].GetZaxis().SetTitleOffset(-.55);
    h['hist%s'%(_proc)].GetZaxis().SetTitleColor(1);
    h['hist%s'%(_proc)].GetZaxis().CenterTitle();
    for _proc in proc:
        print '\nEvaluating process ',_proc
        h['hist%s'%(_proc)] = TH2D('hist%s'%(_proc),'%s Rate of events -  %s '%(_proc,location),binR,rangeRmin,rangeRmax,binN,rangeNmin,rangeNmax)
        h['hist%s'%(_proc)].SetXTitle('distance from wall [m]')
        h['hist%s'%(_proc)].SetYTitle('n9 cut')
        h['hist%s'%(_proc)].SetZTitle('rate per %s'%(t))
        h['hist%s'%(_proc)].GetZaxis().SetTitleOffset(-.55);
        h['hist%s'%(_proc)].GetZaxis().SetTitleColor(1);
        h['hist%s'%(_proc)].GetZaxis().CenterTitle();
        gStyle.SetOptStat(0)
        gStyle.SetPalette(55)
        for _d in drange(rangeRmin+binwidthR/2.,rangeRmax,binwidthR):
            _evts,eff,rateHz,minR,tot = obtainNeutronLike(_cov,_proc,_distance2pmt=_d,_n9=8)
            if rateHz == 0:
                rateHz = minR
            rate = rateHz*24.*3600./timeAdjustment
            h['hist%s'%(_proc)].Fill(_d,rangeNmin+binwidthN/2.0,rate)
            print '\n',_d,eff,rateHz*24.*3600./timeAdjustment,
            for _n in range(int(rangeNmin+binwidthN/2.0+1),int(rangeNmax)):
                _evts,eff,rateHz,minR,tot = obtainNeutronLike(_cov,_proc,_distance2pmt=_d,_n9=_n)
                if rateHz == 0:
                    rateHz = minR
                print rateHz*24.*3600./timeAdjustment,
                rate = rateHz*24.*3600./timeAdjustment
                h['hist%s'%(_proc)].Fill(_d,_n,rate)
        print ''
        h['hist%s'%('Sum')].Add(h['hist%s'%(_proc)],1)
        h['hist%s'%(_proc)].SaveAs('h%s%s.C'%(_proc,_cov))
        h['hist%s'%(_proc)].SaveAs('h%s%s.gif'%(_proc,_cov))
    h['hist%s'%('Sum')].SaveAs('h%s%s.C'%('Sum',_cov))

    procS = ['boulby','neutron']

    for _proc in procS:
        print '\nEvaluating process ',_proc
        h['hist%s'%(_proc)] = TH2D('hist%s'%(_proc),'%s Rate of events -  %s '%(_proc,location),binR,rangeRmin,rangeRmax,binN,rangeNmin,rangeNmax)
        h['hist%s'%(_proc)].SetXTitle('distance from wall [m]')
        h['hist%s'%(_proc)].SetYTitle('n9 cut')
        h['hist%s'%(_proc)].SetZTitle('rate per %s'%(t))
        h['hist%s'%(_proc)].GetZaxis().SetTitleOffset(-.55);
        h['hist%s'%(_proc)].GetZaxis().SetTitleColor(1);
        h['hist%s'%(_proc)].GetZaxis().CenterTitle();

        h['eff%s'%(_proc)] = TH2D('eff%s'%(_proc),'%s Rate of events -  %s '%(_proc,location),binR,rangeRmin,rangeRmax,binN,rangeNmin,rangeNmax)
        h['eff%s'%(_proc)].SetXTitle('distance from wall [m]')
        h['eff%s'%(_proc)].SetYTitle('n9 cut')
        h['eff%s'%(_proc)].SetZTitle('efficiency')
        h['eff%s'%(_proc)].GetZaxis().SetTitleOffset(-.55);
        h['eff%s'%(_proc)].GetZaxis().SetTitleColor(1);
        h['eff%s'%(_proc)].GetZaxis().CenterTitle();

        gStyle.SetOptStat(0)
        gStyle.SetPalette(55)
        for _d in drange(rangeRmin+binwidthR/2.,rangeRmax,binwidthR):
            _evts,eff,rateHz,minR,tot = obtainNeutronLike(_cov,_proc,_distance2pmt=_d,_n9=8,_dist=2.0)
            if rateHz == 0:
                rateHz = minR
                eff = 1./tot
            rate = rateHz*24.*3600./timeAdjustment
            sizeFV    =  2.*pi*pow((pmtRadius/1000.-_d),2)*(pmtHeight/1000.-_d)/1000.
            h['hist%s'%(_proc)].Fill(_d,8,rate)
            h['eff%s'%(_proc)].Fill(_d,8,eff*sizeTank/sizeFV)

            print '\n',_d,eff,rateHz*24.*3600./timeAdjustment,
            for _n in range(9,int(rangeNmax)):
                _evts,eff,rateHz,minR,tot = obtainNeutronLike(_cov,_proc,_distance2pmt=_d,_n9=_n,_dist=2.0)
                if rateHz == 0:
                    rateHz = minR
                    eff = 1./tot
                # print rateHz*24.*3600./timeAdjustment,
                sizeFV    = 2.*pi*pow((pmtRadius/1000.-_d),2)*(pmtHeight/1000.-_d)/1000.
                print eff*sizeTank/sizeFV,
                rate = rateHz*24.*3600./timeAdjustment
                h['hist%s'%(_proc)].Fill(_d,_n,rate)
                # print pmtRadius/1000.
                h['eff%s'%(_proc)].Fill(_d,_n,eff*sizeTank/sizeFV)
        print ''
        h['hist%s'%(_proc)].SaveAs('h%s.C'%(_proc))
        h['eff%s'%(_proc)].SaveAs('eff%s.C'%(_proc))

    timeAcc = 0.0001*86400.

    offsets_n9 = [0,2,4,6,8,10,12,14,16]  ## bin numbers
    offsets_dtw = [0]       ## bin numbers

    for offset in offsets_n9:
        for fv_offset in offsets_dtw:
            if fv_offset < 0:
                _proc = '_%d_neg%d_%s'%(offset,-fv_offset,_cov)
            else:
                _proc = '_%d_%d_%s'%(offset,fv_offset,_cov)

            h['S%s'%(_proc)] = TH2D('S%s'%(_proc),'%s Rate of events -  %s '%(_proc,location),binR,rangeRmin,rangeRmax,binN,rangeNmin,rangeNmax)
            h['S%s'%(_proc)].SetXTitle('distance from wall [m]')
            h['S%s'%(_proc)].SetYTitle('n9 cut')
            h['S%s'%(_proc)].SetZTitle('rate per day')
            h['S%s'%(_proc)].GetZaxis().SetTitleOffset(-.55);
            h['S%s'%(_proc)].GetZaxis().SetTitleColor(1);
            h['S%s'%(_proc)].GetZaxis().CenterTitle();

            h['B%s'%(_proc)] = TH2D('B%s'%(_proc),'%s Rate of events -  %s '%(_proc,location),binR,rangeRmin,rangeRmax,binN,rangeNmin,rangeNmax)
            h['B%s'%(_proc)].SetXTitle('distance from wall [m]')
            h['B%s'%(_proc)].SetYTitle('n9 cut')
            h['B%s'%(_proc)].SetZTitle('rate per day')
            h['B%s'%(_proc)].GetZaxis().SetTitleOffset(-.55);
            h['B%s'%(_proc)].GetZaxis().SetTitleColor(1);
            h['B%s'%(_proc)].GetZaxis().CenterTitle();

            h['SoverB%s'%(_proc)] = TH2D('SoverB%s'%(_proc),'%s Rate of events -  %s '%(_proc,location),binR,rangeRmin,rangeRmax,binN,rangeNmin,rangeNmax)
            h['SoverB%s'%(_proc)].SetXTitle('distance from wall [m]')
            h['SoverB%s'%(_proc)].SetYTitle('n9 cut')
            h['SoverB%s'%(_proc)].SetZTitle('rate per day')
            h['SoverB%s'%(_proc)].GetZaxis().SetTitleOffset(-.55);
            h['SoverB%s'%(_proc)].GetZaxis().SetTitleColor(1);
            h['SoverB%s'%(_proc)].GetZaxis().CenterTitle();

            for _d in range(binR-fv_offset-1):
                for _n in range(binN-offset-1):
                    _db=_d+1
                    _nb=_n+1
                    _p_d  = h['eff%s'%('boulby')].GetXaxis().GetBinCenter(_db)
                    _p_n9 = h['eff%s'%('boulby')].GetYaxis().GetBinCenter(_nb)
                    _n_d  = h['eff%s'%('neutron')].GetXaxis().GetBinCenter(_db+fv_offset)
                    _n_n9 = h['eff%s'%('neutron')].GetYaxis().GetBinCenter(_nb+offset)
                    _p_v  = h['eff%s'%('boulby')].GetBinContent(_db,_nb)
                    _n_v  = h['eff%s'%('neutron')].GetBinContent(_db+fv_offset,_nb+offset)
                    _rate_v  = h['hist%s'%('neutron')].GetBinContent(_db+fv_offset,_nb+offset)
                    print "Positron/neutron: Wall distance (%4.1f,%4.1f), n9 cut (%d,%d), efficiency (%4.3f,%4.3f): combined eff/rate : %4.3f per day"\
                    %(_p_d,_n_d,_p_n9,_n_n9,_p_v,_n_v,_rate_v*_p_v*86400.)
                    _signal = _rate_v*_p_v*86400.

                    _p_d  = h['hist%s'%('Sum')].GetXaxis().GetBinCenter(_db)
                    _p_n9 = h['hist%s'%('Sum')].GetYaxis().GetBinCenter(_nb)
                    _n_d  = h['hist%s'%('Sum')].GetXaxis().GetBinCenter(_db+fv_offset)
                    _n_n9 = h['hist%s'%('Sum')].GetYaxis().GetBinCenter(_nb+offset)
                    _p_v  = h['hist%s'%('Sum')].GetBinContent(_db,_nb)
                    _n_v  = h['hist%s'%('Sum')].GetBinContent(_db+fv_offset,_nb+offset)

                    print "Accidental       : Wall distance (%4.1f,%4.1f), n9 cut (%d,%d), rate (%4.3f,%4.3f): combined rate : %4.3f per day"\
                    %(_p_d,_n_d,_p_n9,_n_n9,_p_v,_n_v,_p_v*_n_v*timeAcc)
                    _background = _p_v*_n_v*timeAcc*0.05

                    h['S%s'%(_proc)].SetBinContent(_db,_nb,_signal)
                    h['B%s'%(_proc)].SetBinContent(_db,_nb,_background)
                    h['SoverB%s'%(_proc)].SetBinContent(_db,_nb,_signal/sqrt(_signal+_background))


            h['S%s'%(_proc)].SaveAs("ibdSignal%s.C"%(_proc))
            h['B%s'%(_proc)].SaveAs("ibdBackground%s.C"%(_proc))
            h['SoverB%s'%(_proc)].SaveAs("ibdSignalOverBackground%s.C"%(_proc))
            print _proc,h['SoverB%s'%(_proc)].GetMaximum()



def sensitivityMapPass2New():

    site = arguments["--site"]

    # Need to fix this for future running

    OnOffRatio = float(arguments["--OnOff"])
    print site,'with on-off ratio of ',OnOffRatio

    cores = int(arguments["--cores"])

    if arguments["--RNRedux"]:
        rnRedux = float(arguments["--RNRedux"])
        if rnRedux>1:
            print "Value of reduction of radionuclide greater than 1, setting to 0"
            rnRedux = 0.0
    else:
        rnRedux = 0.0

    if t == 'sec':
        timeAdjustment = 24*3600.
    if t == 'day':
        timeAdjustment = 1.0
    if t == 'month':
        timeAdjustment = 1./31.
    if t == 'year':
        timeAdjustment = 1./365.
    maxTime = 14400.*timeAdjustment

    print '\nEvaluation based on geoneutrinos.org'
    #parameters  = loadAnalysisParametersNew(t)
    rates       = {}#parameters[11]
    rates["boulby_S"]=1.0
    rates["imb_S"]=1.0
    print 'Wrong rates for now'
    sizeDetc    = 2.*pi*pow(fidRadius/1000.,2)*fidHeight/1000./1000.
    sizeTank    = 2.*pi*pow(detectorRadius/1000.,2)*detectorHeight/1000./1000.
    FVkTonRatio = (pow(fidRadius,2)*fidHeight)/(pow(detectorRadius,2)*detectorHeight)
    boulbyRate,imbRate = rates["boulby_S"]*FVkTonRatio,rates["imb_S"]*FVkTonRatio
    print ' boulby rates: %4.2e per %s per %4.2f kton; [r: %4.2f m; z: %4.2f m]'\
    %(boulbyRate,t,sizeDetc,fidRadius/1000.,fidHeight/1000.)
    #fast neutrons
    # print 'Debug',rates["boulby_S"],FVkTonRatio
    print '\nEvaluation not based on geoneutrinos.org'
    detectorMedium,detectorMass,reactorPower,reactorStandoff = 1,sizeDetc*1000.,1.5,24.98
    experiment = nuOsc.NeutrinoOscillation(detectorMedium,detectorMass,reactorPower,reactorStandoff)
    preOsc,afterOsc = experiment.FindRate()
    print ' Neutrino rate pre osc: %4.2f; neutrino rate post osc: %4.2f at %4.2f GWth, at %4.2f km, for %4.2f kton' %(preOsc,afterOsc,reactorPower,reactorStandoff,detectorMass/1000.)
    detectorMedium,detectorMass,reactorPower,reactorStandoff = 1,sizeDetc*1000.,1.575,24.98
    experiment = nuOsc.NeutrinoOscillation(detectorMedium,detectorMass,reactorPower,reactorStandoff)
    preOsc,afterOsc = experiment.FindRate()
    print ' Neutrino rate pre osc: %4.2f; neutrino rate post osc: %4.2f at %4.2f GWth, at %4.2f km, for %4.2f kton' %(preOsc,afterOsc,reactorPower,reactorStandoff,detectorMass/1000.)
    print ''

    proc,loca,type,color,lineS,acc,scale   = [],[],[],[],[],[],[]

    proc        += ['QGSP_BERT_EMV','QGSP_BERT_EMX','QGSP_BERT','QGSP_BIC',\
    'QBBC','QBBC_EMZ','FTFP_BERT','QGSP_FTFP_BERT']
    _t          = 'FN%s' % (site)
    loca        += [_t,_t, _t,_t,_t,_t,_t,_t]
    type        += ['si','si','si','si','si','si','si','si']
    acc         += ['corr','corr','corr','corr','corr','corr','corr', 'corr']
    color       += [kO+6,kO+6,kO+6,kO+6,kO+6,kO+6,kO+6,kO+6]
    lineS       += [2,2,2,2,2,2,2,2]
    scale       += [1./8.,1./8.,1./8.,1./8.,1./8.,1./8.,1./8.,1./8.]

    #radionuclides
    proc        += ['9003','11003']
    _t          =  'RN%s' % (site)
    loca        += [_t,_t]
    type        += ['ei','ei']
    acc         += ['di','di']
    color       += [kG+3,kG+3]
    lineS       += [1,1]
    scale       += [1.,1.]
    #ibds
    if site == 'boulby':
        proc    += ['boulby','boulby','neutron']
    else:
        proc    += ['imb','imb','neutron']
    loca        += ['S','S','N%s'%(site)]
    type        += ['ei','ei','ei']
    acc         += ['di','corr','corr']
    color       += [kA+0,kA-0,kA-0]
    lineS       += [1,2,2]
    scale       += [-1.0,0.0,0.0]

    c1 = TCanvas('c1','c1',1618,1000)
    c1.SetRightMargin(0.53)
    c1.Divide(2,2)

    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)
    if additionalString == "":
        additionalString = "_default"
    location = 'Boulby '

    hist = TH2D('hist','3#sigma discovery phase space -  %s '%(location),31,9.5,40.5,30,9.5,39.5)
    hist.SetXTitle('photocoverage [%]')
    hist.SetYTitle('photoelectron threhsold cut [p.e.]')
    hist.SetZTitle('off-time [%s] for 3 #sigma discovery'%(t))
    hist.GetZaxis().SetTitleOffset(-.55);
    hist.GetZaxis().SetTitleColor(1);
    hist.GetZaxis().CenterTitle();
    gStyle.SetOptStat(0)
    gStyle.SetPalette(55)


    h = {}

    binR,rangeRmin,rangeRmax = 31,0.45,3.55
    binwidthR = (rangeRmax-rangeRmin)/binR
    binN,rangeNmin,rangeNmax = 48,7.5,55.5
    binwidthN = (rangeNmax-rangeNmin)/binN

    _cov = arguments['-C']

    d,proc,coverage = loadSimulationParametersNew()

    minAchieve = 0

    for _p in proc:
        for _loc in proc[_p]:
            for idx,_cover in enumerate(coverage):
                for _element in d[_p]:
                    _tag = "%s_%s_%s_%s"%(_cover,_loc,_element,_p)
                    _file = "bonsai_root_files%s/%s/merged_%s_%s_%s.root"%(additionalMacStr,_cover,_loc,_element,_p)
                    print _tag
                    obtainEventEfficiency(_cover,_file,_tag)

                    print ''

def readEfficiencyHistogram():

    hist = {}

    d,proc,coverage = loadSimulationParametersNew()
    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)

    print 'Loading in all .C file...'
    for _p in proc:
        for _loc in proc[_p]:
            for idx,_cover in enumerate(coverage):
                for _element in d[_p]:
                    _tag = "%s_%s_%s_%s"%(_cover,_loc,_element,_p)
                    _hist = 'hist'+_tag
                    _dir = "bonsai_root_files%s/%s/"%(additionalMacStr,_cover)
                    _file =  _dir+'hist'+_tag+'.C'
                    gROOT.ProcessLine('.x %s'%(_file))
                    hist[_tag] =  TH2D()
                    gROOT.GetObject(_hist,hist[_tag])

    print 'done processing all .C files. What is in the directory:'
    gROOT.ProcessLine('.ls')

    print '\nLoading PMT activity:'
    mPMTs,mPMTsU238,mPMTsTh232,mPMTsK40 = loadPMTActivity()
    print 'done.'

    tankRadius  = float(arguments["--tankRadius"])-float(arguments['--steelThick'])
    tankHeight  = float(arguments["--halfHeight"])-float(arguments['--steelThick'])
    nKiloTons = pi*pow(tankRadius/1000.,2)*(2.*tankHeight/1000.)
    rRn222 = float(arguments["--Rn222"])*nKiloTons
    print '\nLoaded Rn-222 activity of ',rRn222,'Bq per water volume, assumed a rate of %4.3e Bq/m^3'%(float(arguments["--Rn222"]))


    FreeProtons = 0.668559
    TNU         = FreeProtons* nKiloTons /1000./365./24./3600.
    boulbyIBDRate   = 800.*TNU
    print '\nLoaded an IBD rate of ',boulbyIBDRate,' events per water volume per second, assumed a rate of %4.3e TNU'%(boulbyIBDRate/TNU)

    innerRad = 12.5 #meters
    outerRad = 13.5 #meters
    rockMass = (pi*pow(outerRad,2)*(2.*outerRad)-pi*pow(innerRad,2)*(2.*innerRad))*power(100.,3)*2.39
    # volumeR         = power(22.,3)-power(20.,3)# Rock cavern e.g. (22m x 22m x 22m) - (20m x 20m x 20m)
    # density         = 2.39 #from McGrath
    # rockMass        = volumeR*power(100.,3)*density
    #Mass of rock evalyated
    avgMuon         = npa([180.,264.])
    avgMuonNC       = power(avgMuon,0.849)
    avgNFluxMag     = 1e-6
    muonRate        = npa([7.06e-7,4.09e-8]) # mu/cm2/s
    tenMeVRatio     = npa([7.51/34.1,1.11/4.86])
    fastNeutrons    = rockMass*avgMuonNC*avgNFluxMag*muonRate*tenMeVRatio
    FN_boulby = fastNeutrons[1]

    avgRNYieldRC    = power(avgMuon,0.73)
    skRNRate        = 0.5e-7 # 1/mu/g cm2
    avgMuonSK       = power(219.,0.73)
    skMuFlux        = 1.58e-7 #mu/cm2/sec
    radionuclideRate= (skRNRate*avgRNYieldRC/avgMuonSK)*muonRate*nKiloTons*1e9
    RN_boulby        = radionuclideRate[1]
    print '\nLoaded mass of rock %e g. Fast Neutron Yield %e per sec; radionuclide yield %e per sec'%(rockMass,FN_boulby,RN_boulby)




    print '\n What are the maximum efficiency/rate found in each histogram:'
    _sing = 0.0
    lineU238PMT,lineTh232PMT,lineKPMT = '','',''
    lineFNROCK = ''
    lineU238GUN,lineTh232GUN,lineKGUN = '','',''
    lineU238ROCK,lineTh232ROCK,lineKROCK = '','',''
    lineU238CONC,lineTh232CONC,lineKCONC = '','',''
    lineRn222WaterVolume = ''
    lineU238GD,lineU235GD,lineTH232GD = '','',''
    linePromptWaterVolume,lineDelayedWaterVolume,linePromptDelayedWaterVolume = '','',''
    lineELSE = ''

    firstGo = 1
    for _t in hist:
        if firstGo:
            h = hist[_t].Clone()
            h.SetZTitle('singles rate (Hz)')
            h.SetTitle('Singles rate')
            h.SetName('hSinglesRate')
            h.Reset()
            hn = hist[_t].Clone()
            hn.SetZTitle('neutron rate (Hz)')
            hn.SetTitle('rate')
            hn.SetName('hNeutronRate')
            hn.Reset()
            hee = hist[_t].Clone()
            hee.SetZTitle('positron efficiency')
            hee.SetTitle('efficiency')
            hee.SetName('hPositronEfficiency')
            hee.Reset()
            hne = hist[_t].Clone()
            hne.SetZTitle('neutron efficiency')
            hne.SetTitle('efficiency')
            hne.SetName('hNeutronEfficiency')
            hne.Reset()

            firstGo =0
        if 'PMT' in _t and 'CHAIN_238U_NA' in _t:
            if '210Tl' in _t:
                _sing+=hist[_t].GetMaximum()*mPMTsU238[0]*0.002
                h.Add(hist[_t],mPMTsU238[0]*0.002)
                lineU238PMT += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsU238[0]*0.002)
            else:
                _sing+=hist[_t].GetMaximum()*mPMTsU238[0]
                h.Add(hist[_t],mPMTsU238[0])
                lineU238PMT+= "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsU238[0])
        elif 'PMT' in _t and 'CHAIN_232Th_NA' in _t:
            _sing+=hist[_t].GetMaximum()*mPMTsTh232[0]
            h.Add(hist[_t],mPMTsTh232[0])
            lineTh232PMT += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsTh232[0])
        elif 'PMT' in _t and '40K_NA' in _t:
            _sing+=hist[_t].GetMaximum()*mPMTsK40[0]
            h.Add(hist[_t],mPMTsK40[0])
            lineKPMT += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsK40[0])

        elif 'GUNITE' in _t and 'CHAIN_238U_NA' in _t:
            if '210Tl' in _t:
                _sing+=hist[_t].GetMaximum()*mPMTsU238[0]*0.002
                lineU238GUN += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsU238[0]*0.002)
            else:
                _sing+=hist[_t].GetMaximum()*mPMTsU238[0]
                lineU238GUN+= "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsU238[0])
        elif 'GUNITE' in _t and 'CHAIN_232Th_NA' in _t:
            _sing+=hist[_t].GetMaximum()*mPMTsTh232[0]
            lineTh232GUN += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsTh232[0])
        elif 'GUNITE' in _t and '40K_NA' in _t:
            _sing+=hist[_t].GetMaximum()*mPMTsK40[0]
            lineKGUN += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsK40[0])

        elif 'ROCK' in _t and 'CHAIN_238U_NA' in _t:
            if '210Tl' in _t:
                _sing+=hist[_t].GetMaximum()*mPMTsU238[0]*0.002
                lineU238ROCK += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsU238[0]*0.002)
            else:
                _sing+=hist[_t].GetMaximum()*mPMTsU238[0]
                lineU238ROCK+= "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsU238[0])
        elif 'ROCK' in _t and 'CHAIN_232Th_NA' in _t:
            _sing+=hist[_t].GetMaximum()*mPMTsTh232[0]
            lineTh232ROCK += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsTh232[0])
        elif 'ROCK' in _t and '40K_NA' in _t:
            _sing+=hist[_t].GetMaximum()*mPMTsK40[0]
            lineKROCK += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsK40[0])
        elif 'ROCK' in _t and '_FN' in _t:
            lineFNROCK += "%50s %e %15.10f (%15.10f per day)\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*FN_boulby,hist[_t].GetMaximum()*FN_boulby*3600.*24.)
        elif 'CONC' in _t and 'CHAIN_238U_NA' in _t:
            if '210Tl' in _t:
                _sing+=hist[_t].GetMaximum()*mPMTsU238[0]*0.002
                lineU238CONC += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsU238[0]*0.002)
            else:
                _sing+=hist[_t].GetMaximum()*mPMTsU238[0]
                lineU238CONC+= "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsU238[0])
        elif 'CONC' in _t and 'CHAIN_232Th_NA' in _t:
            _sing+=hist[_t].GetMaximum()*mPMTsTh232[0]
            lineTh232CONC += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsTh232[0])
        elif 'CONC' in _t and '40K_NA' in _t:
            _sing+=hist[_t].GetMaximum()*mPMTsK40[0]
            lineKCONC += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsK40[0])

	elif 'GD' in _t and 'CHAIN_238U_NA' in _t:
            if '210Tl' in _t:
                _sing+=hist[_t].GetMaximum()*GdU238_l*0.002
                h.Add(hist[_t],GdU238_l*0.002)
                lineU238GD += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*GdU238_l*0.002)
            elif '234Pa' in _t:
                _sing+=hist[_t].GetMaximum()*GdU238
                h.Add(hist[_t],GdU238)
                lineU238GD += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*GdU238)
            else:
                _sing+=hist[_t].GetMaximum()*GdU238_l
                h.Add(hist[_t],GdU238_l)
                lineU238GD += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*GdU238_l)
        elif 'GD' in _t and 'CHAIN_232Th_NA' in _t:
	    if '228Ac' in _t:
            	_sing+=hist[_t].GetMaximum()*GdTh232
                h.Add(hist[_t],GdTh232)
            	lineTH232GD += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*GdTh232)
	    else:
            	_sing+=hist[_t].GetMaximum()*GdTh232_l
                h.Add(hist[_t],GdTh232_l)
            	lineTH232GD += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*GdTh232_l)
        elif 'GD' in _t and 'CHAIN_235U_NA' in _t:
	    if '231Th' in _t:
            	_sing+=hist[_t].GetMaximum()*GdU235
                h.Add(hist[_t],GdU235)
            	lineU235GD += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*GdU235)
	    else:
            	_sing+=hist[_t].GetMaximum()*GdU235_l
                h.Add(hist[_t],GdU235_l)
            	lineU235GD += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*GdU235_l)

        elif 'WaterVolume' in _t and 'CHAIN_222Rn_NA' in _t:
            if '210Tl' in _t:
                _sing+=hist[_t].GetMaximum()*rRn222*0.002
                h.Add(hist[_t],rRn222*0.002)
                lineRn222WaterVolume += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*rRn222*0.002)
            else:
                _sing+=hist[_t].GetMaximum()*rRn222
                h.Add(hist[_t],rRn222)
                lineRn222WaterVolume+= "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*rRn222)

        elif 'WaterVolume' in _t and 'promptPositron' in _t:
            _day=hist[_t].GetMaximum()*boulbyIBDRate*3600.*24
            hee.Add(hist[_t],1.)#Only add efficiency for positron
            linePromptWaterVolume += "%50s %e %15.10f per sec (%15.10f per day)\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*boulbyIBDRate,_day)

        elif 'WaterVolume' in _t and 'delayedNeutron' in _t:
            _day=hist[_t].GetMaximum()*boulbyIBDRate*3600.*24
            hn.Add(hist[_t],boulbyIBDRate)
            hne.Add(hist[_t],1.)
            linePromptWaterVolume += "%50s %e %15.10f per sec (%15.10f per day)\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*boulbyIBDRate,_day)

        elif 'WaterVolume' in _t and 'promptDelayedPair' in _t:
            _day=hist[_t].GetMaximum()*boulbyIBDRate*3600.*24
            linePromptWaterVolume += "%50s %e %15.10f per sec (%15.10f per day)\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*boulbyIBDRate,_day)

        else:
            lineELSE += "%50s %e\n"%(_t,hist[_t].GetMaximum())


                    # hist =
                    # print _hist.GetMaximum()
    print ''
    print 'PMT U-238 \n', lineU238PMT
    print 'PMT Th-232\n', lineTh232PMT
    print 'PMT K\n', lineKPMT,'\n'
    print 'Water Rn-222\n',lineRn222WaterVolume,'\n'
    print 'Gunite U-238 \n', lineU238GUN
    print 'Gunite Th-232\n', lineTh232GUN
    print 'Gunite K\n', lineKGUN,'\n'
    print 'Concrete U-238 \n', lineU238CONC
    print 'Concrete Th-232\n', lineTh232CONC
    print 'Concrete K\n', lineKCONC,'\n'
    print 'ROCK U-238 \n', lineU238ROCK
    print 'ROCK Th-232\n', lineTh232ROCK
    print 'ROCK K\n', lineKROCK,'\n'
    print 'ROCK Fast Neutron\n',lineFNROCK,'\n'
    print 'Else  \n', lineELSE,'\n'
    print 'Total singles rate:\t\t\t',_sing,'events per sec at minimum buffer distance of 0.5 m\n'
    print 'Signal information'
    print 'Prompt positron Water volume \n', linePromptWaterVolume
    signal = ['WaterVolume_delayedNeutron_ibd_n','WaterVolume_promptPositron_ibd_p']
    _str =  "bonsai_root_files%s/%s/histograms_U238_%4.3fPPM_Th232_%4.3fPPM_K_%4.3fPPM.root"%(additionalMacStr,_cover,float(arguments["--U238_PPM"]),float(arguments["--Th232_PPM"]),float(arguments["--K_PPM"]))


    timeAcc = 0.0001*86400.

    offsets_n9 = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]  ## bin numbers
    offsets_dtw = [0]       ## bin numbers

    binR,rangeRmin,rangeRmax = 31,0.45,3.55
    binwidthR = (rangeRmax-rangeRmin)/binR
    binN,rangeNmin,rangeNmax = 48,7.5,55.5
    binwidthN = (rangeNmax-rangeNmin)/binN

    _cov = arguments['-C']


    pmtRadius = float(arguments['--tankRadius'])-float(arguments['--steelThick'])-float(arguments['--shieldThick'])
    pmtHeight = float(arguments['--halfHeight'])-float(arguments['--steelThick'])-float(arguments['--shieldThick'])
    detectorRadius  = float(arguments['--tankRadius'])-float(arguments['--steelThick'])
    detectorHeight  = float(arguments['--halfHeight'])-float(arguments['--steelThick'])

    _maxSignal,_maxBkgd,_maxSoverB,_maxOffn9,_maxOff_dtw = -1,-1,-1,-1,-1
    _maxSignal2,_maxBkgd2,_maxSoverB2,_maxOffn92,_maxOff_dtw2,_maxOffset2 = -1,-1,-1,-1,-1,-1


    line,_line,_line2= ("",),"",""
    _histograms = {}
    for offset in offsets_n9:
        for fv_offset in offsets_dtw:
            _histograms["sOverB_%d"%(offset)]= h.Clone()
            _histograms["sOverB_%d"%(offset)].SetZTitle('signal/sqrt(signal+background)')
            _histograms["sOverB_%d"%(offset)].SetYTitle('n9 cut on prompt')
            _histograms["sOverB_%d"%(offset)].SetTitle('hSoB - offset %2d'%(offset))
            _histograms["sOverB_%d"%(offset)].SetName('hSoB%d'%(offset))
            _histograms["sOverB_%d"%(offset)].Reset()

            _histograms["signal_%d"%(offset)]= h.Clone()
            _histograms["signal_%d"%(offset)].SetZTitle('evts/day')
            _histograms["signal_%d"%(offset)].SetYTitle('n9 cut on prompt')
            _histograms["signal_%d"%(offset)].SetTitle('Signal - offset %2d'%(offset))
            _histograms["signal_%d"%(offset)].SetName('hSignal%d'%(offset))
            _histograms["signal_%d"%(offset)].Reset()

            _histograms["background_%d"%(offset)]= h.Clone()
            _histograms["background_%d"%(offset)].SetZTitle('evts/day')
            _histograms["background_%d"%(offset)].SetYTitle('n9 cut on prompt')
            _histograms["background_%d"%(offset)].SetTitle('accidentals - offset %2d'%(offset))
            _histograms["background_%d"%(offset)].SetName('hBackground%d'%(offset))
            _histograms["background_%d"%(offset)].Reset()



            _proc = '_%d_%d_%s'%(offset,fv_offset,_cov)

            _maxSignal,_maxBkgd,_maxSoverB,_maxOffn9,_maxOff_dtw = -1,-1,-1,-1,-1
            for _d in range(binR-fv_offset-1):
                for _n in range(binN-offset-1):
                    _db,_nb=_d+1,_n+1
                    _p_d  = hee.GetXaxis().GetBinCenter(_db)
                    _p_n9 = hee.GetYaxis().GetBinCenter(_nb)
                    _n_d  = hne.GetXaxis().GetBinCenter(_db+fv_offset)
                    _n_n9 = hne.GetYaxis().GetBinCenter(_nb+offset)
                    fidRadius,fidHeight = pmtRadius-_p_d*1e3,pmtHeight-_p_d*1e3
                    ratioScaling = (pow(fidRadius,2)*fidHeight)/(pow(detectorRadius,2)*detectorHeight)
                    _p_v  = hee.GetBinContent(_db,_nb)/ratioScaling
                    _n_v  = hne.GetBinContent(_db+fv_offset,_nb+offset)/ratioScaling
                    _rate_v  = hn.GetBinContent(_db+fv_offset,_nb+offset)

                    _signal = _rate_v*_p_v*86400.

                    _p_d_1  = h.GetXaxis().GetBinCenter(_db)
                    _p_n9_1 = h.GetYaxis().GetBinCenter(_nb)
                    _n_d_1  = h.GetXaxis().GetBinCenter(_db+fv_offset)
                    _n_n9_1 = h.GetYaxis().GetBinCenter(_nb+offset)
                    _p_v_1  = h.GetBinContent(_db,_nb)
                    _n_v_1  = h.GetBinContent(_db+fv_offset,_nb+offset)

                    _background = _p_v_1*_n_v_1*timeAcc*0.05
                    sob = _signal/sqrt(_signal+_background)

                    _histograms["sOverB_%d"%(offset)].SetBinContent(_db,_nb,sob)
                    _histograms["signal_%d"%(offset)].SetBinContent(_db,_nb,_signal)
                    _histograms["background_%d"%(offset)].SetBinContent(_db,_nb,_background)

                    if sob >_maxSoverB:
                        _maxSoverB = _signal/sqrt(_signal+_background)
                        _maxSignal = _signal
                        _maxBkgd   = _background
                        _maxOffn9   = _n_n9
                        _maxOff_dtw = _n_d
                        _line = ("Offset:%3d, beta/n: Wall-dist (%4.1f,%4.1f) m, n9 cut (%d,%d), rel. efficiency (%4.2f,%4.2f), neutron rate :(%4.2f per day), combined eff/rate : %4.2f per day;"\
                            %(offset,_p_d,_n_d\
                            ,_p_n9,_n_n9\
                            ,_p_v,_n_v\
                            ,_rate_v*86400.\
                            ,_rate_v*_p_v*86400.),)
                        _line2 =    ("acc. rate (%5.3f,%5.3f): acc. combined rate: %4.3f per day (pre-prox)"\
                         %(_p_v_1,_n_v_1,_p_v_1*_n_v_1*timeAcc),)

                    if sob >_maxSoverB2:
                        _maxSoverB2 = _signal/sqrt(_signal+_background)
                        _maxSignal2 = _signal
                        _maxBkgd2   = _background
                        _maxOffn92   = _n_n9
                        _maxOff_dtw2 = _n_d
                        _maxOffset2 = offset
                        # _line = ("Offset:%3d, beta/n: Wall-dist (%4.1f,%4.1f) m, n9 cut (%d,%d), rel. efficiency (%4.2f,%4.2f), neutron rate :(%4.2f per day), combined eff/rate : %4.2f per day;"\
                        #     %(offset,_p_d,_n_d\
                        #     ,_p_n9,_n_n9\
                        #     ,_p_v,_n_v\
                        #     ,_rate_v*86400.\
                        #     ,_rate_v*_p_v*86400.),)
                        # _line2 =    ("acc. rate (%5.3f,%5.3f): acc. combined rate: %4.3f per day (pre-prox)"\
                        #  %(_p_v_1,_n_v_1,_p_v_1*_n_v_1*timeAcc),)
                        #


            print 'Offset:',str(offset).rjust(3,' '),',Found max S/sqrt(S+B)',_maxSoverB,',(S,B,n9,dtw):(',_maxSignal,_maxBkgd,_maxOffn9,_maxOff_dtw,')'
            line += (_line + _line2,)

    print '\n\nMore info on the maximal sensitivity found:'
    # print line

    for _l in line:
        for i in range(len(_l)):
            print _l[i],
        print ''

    cut_signal = 0.9
    sigma = 4.65 # 3-sigma 95% of the time, according to Owen report
    _S = _maxSignal2*cut_signal
    _B = _maxSignal2*cut_signal*1.15 + _maxBkgd2 # Includes a 15% other reactor component
    OnOffRatio = 1.5

    if _S >0:
        T3SIGMA = sigma**2*(_B +(_S+_B)/OnOffRatio)/_S/_S

    else:
        T3SIGMA = 1e999
    metric = T3SIGMA*OnOffRatio + T3SIGMA
    _res = "%s %4.1f %3d %3d %4.3f %4.3f %4.1f %4.1f" %(_cover,_maxOff_dtw2,_maxOffn92,_maxOffn92-_maxOffset2,_maxSignal2,_maxBkgd2,T3SIGMA,metric)
    _strRes = "results_DTW_%dmm_U238_%4.3fPPM_Th232_%4.3fPPM_K_%4.3fPPM.txt"%(float(arguments['--shieldThick']),float(arguments["--U238_PPM"]),float(arguments["--Th232_PPM"]),float(arguments["--K_PPM"]))
    # _strRes = 'res.txt'
    with open(_strRes,'a') as file:
        file.write(_res+'\n')

    print '\n\nWriting histograms to file',_str
    f_root = TFile(_str,"recreate")
    h.Write()
    hn.Write()
    hne.Write()
    hee.Write()
    for offset in offsets_n9:
        _histograms["sOverB_%d"%(offset)].Write()
        _histograms["signal_%d"%(offset)].Write()
        _histograms["background_%d"%(offset)].Write()
    f_root.Close()

    print '\n\nall done.'



def findRate():


    d,proc,coverage = loadSimulationParametersNew()
    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)

    print '\nLoading TANK activity:'
    tankmass,tankact_60co,tankact_137cs = loadTankActivity()
    print 'done.'

    print '\nLoading CONCRETE activity:'
    concmass,concact_238u,concact_232th,concact_40k = loadConcreteActivity()
    print 'done.'

    print '\nLoading Shotcrete activity:'
    shotmass,act_238u,act_232th,act_40k = loadShotcreteActivity()
    print 'done.'

    print '\nLoading ROCK activity:'
    rockmass,act_238u,act_232th,act_40k = loadRockActivity()
    print 'done.'


    print '\nLoading PMT activity:'
    mPMTs,mPMTsU238,mPMTsTh232,mPMTsK40 = loadPMTActivity()
    print 'done.'

    tankRadius  = float(arguments["--tankRadius"])-float(arguments['--steelThick'])
    tankHeight  = float(arguments["--halfHeight"])-float(arguments['--steelThick'])
    nKiloTons = pi*pow(tankRadius/1000.,2)*(2.*tankHeight/1000.)
    rRn222 = float(arguments["--Rn222"])*nKiloTons
    print '\nLoaded Rn-222 activity of ',rRn222,'Bq per water volume, assumed a rate of %4.3e Bq/m^3'%(float(arguments["--Rn222"]))

    FreeProtons = 0.668559
    TNU         = FreeProtons* nKiloTons /1000./365./24./3600.
    boulbyIBDRate   = 800.*TNU
    print '\nLoaded an IBD rate of ',boulbyIBDRate,' events per water volume per second, assumed a rate of %4.3e TNU'%(boulbyIBDRate/TNU)

    innerRad = 12.5 #meters
    outerRad = 13.5 #meters
    rockMass = (pi*pow(outerRad,2)*(2.*outerRad)-pi*pow(innerRad,2)*(2.*innerRad))*power(100.,3)*2.39
    # volumeR         = power(22.,3)-power(20.,3)# Rock cavern e.g. (22m x 22m x 22m) - (20m x 20m x 20m)
    # density         = 2.39 #from McGrath
    # rockMass        = volumeR*power(100.,3)*density
    #Mass of rock evalyated
    avgMuon         = npa([180.,264.])
    avgMuonNC       = power(avgMuon,0.849)
    avgNFluxMag     = 1e-6
    muonRate        = npa([7.06e-7,4.09e-8]) # mu/cm2/s
    tenMeVRatio     = npa([7.51/34.1,1.11/4.86])
    fastNeutrons    = rockMass*avgMuonNC*avgNFluxMag*muonRate*tenMeVRatio
    FN_boulby = fastNeutrons[1]

    avgRNYieldRC    = power(avgMuon,0.73)
    skRNRate        = 0.5e-7 # 1/mu/g cm2
    avgMuonSK       = power(219.,0.73)
    skMuFlux        = 1.58e-7 #mu/cm2/sec
    radionuclideRate= (skRNRate*avgRNYieldRC/avgMuonSK)*muonRate*nKiloTons*1e9
    RN_boulby        = radionuclideRate[1]
    print '\nLoaded mass of rock %e g. Fast Neutron Yield %e per sec; radionuclide yield %e per sec'%(rockMass,FN_boulby,RN_boulby)
    _strRes = "rate_%dmm_U238_%4.3fPPM_Th232_%4.3fPPM_K_%4.3fPPM_%s.txt"%(float(arguments['--shieldThick']),float(arguments["--U238_PPM"]),float(arguments["--Th232_PPM"]),float(arguments["--K_PPM"]),arguments["-C"])
    _res = "%e %e %d %e %e %e"%(boulbyIBDRate,rRn222,mPMTs[0],mPMTsU238[0],mPMTsTh232[0],mPMTsK40[0])
    with open(_strRes,'a') as file:
        file.write(_res+'\n')

    return boulbyIBDRate,rRn222,mPMTsU238,mPMTsTh232,mPMTsK40


def runSensitivity():
    hBoulby = TH2D('hBoulby','hBoulby',50,0.5,50.5,50,0.5,50.5)
    print 'Boulby:'
    for cut in range(4,40):
        Graphs = obtainAbsoluteEfficiency('processed_data_watchman.root',\
        timeScale=t,cut=cut)
        hBoulby =  siteCanvas(Graphs,'boulby',cut,hBoulby)
    hBoulby.SaveAs('hBoulby.C')


    hFairport = TH2D('hFairport','hFairport',50,0.5,50.5,50,0.5,50.5)
    print 'Fairport:'
    for cut in range(4,40):
        Graphs = obtainAbsoluteEfficiency('processed_data_watchman.root',\
        timeScale=t,cut=cut)
        hFairport =  siteCanvas(Graphs,'',cut,hFairport)
    hFairport.SaveAs('hFairport.C')

    return 0

def findScaledVolume(S,B_Surface,B_Volume,OnOffRatio = 1.0):
    B_Volume += S*0.02 # Reactor background for site of interest 15TNU/747TNU
    S *= 40./1500. # Transform Boubly rate to 40MWth rate

    for size in drange(1.,500.,0.5):
        _rad    = pow(size*1000./pi/2.,1./3.)
        _surf   = pow(_rad,2)/pow(fidRadius,2)
        _volume = pow(_rad,3)/(pow(fidRadius,2)*fidHeight)
        _S      = S * _volume
        _B      = B_Surface * _surf + B_Volume * _volume
        SSBB   = _S/sqrt(_B + (_S + _B)/OnOffRatio)

        T3SIGMA = 9.*(_B +(_S+_B)/OnOffRatio)/_S/_S
        # except:
        #     print 'Could not'
        # print "%3d %10.2f %10.2f %10.2f %10.2f"%(size, SSBB, T3SIGMA,_S,_B)
        if T3SIGMA < 3.00001:
            break
    return size,T3SIGMA



def findEquivalentStandoff(S,B_Surface,B_Volume,_signal_Eff,_sigma_b = 0.1,\
t = 3.0,detectorMass=1.,N_sigma=3):
    B_Volume += S*0.02 # Reactor background for site of interest 15TNU/747TNU

    # First, find rate on observed signal needed
    _rad    = pow(detectorMass*1000./pi/2.,1./3.)
    _surf   = pow(_rad,2)/pow(fidRadius,2)
    _volume = pow(_rad,3)/(pow(fidRadius,2)*fidHeight)
    _b = B_Surface*_surf + B_Volume*_volume
    _sigma_b*=_b

    # Evaluate the required signal to obtain N_sigma over t
    C_square = N_sigma*N_sigma
    _s = C_square/(2.*t)*(1+sqrt(1. + 4.*_b*t/C_square+\
    4.*_sigma_b*_sigma_b*t*t/C_square))

    if _s > _sigma_b:
        _t = C_square*(_s+_b)/(_s*_s - _sigma_b*_sigma_b*C_square)
    else:
        print 'Background systematics greater than signal'
        return -1,-1,-1
    print 'd,r,s,v: ',detectorMass,_rad,_surf,_volume,_b,_sigma_b,_s,t,_t,

    # Define the medium of the detector (1:water) and power of test reactor
    detectorMedium,reactorPower      = 1,0.04 #0.04
    # Cycle through the standoff and stop when signal is observed
    experiment = nuOsc.NeutrinoOscillation(detectorMedium,\
    1.*1000.,reactorPower,25)
    preOsc,afterOsc = experiment.FindRate()
    print preOsc,afterOsc,
    for _standoff in drange(0.05,50,0.05):
            experiment = nuOsc.NeutrinoOscillation(detectorMedium,\
            detectorMass*1000.,reactorPower,_standoff)
            preOsc,afterOsc = experiment.FindRate()
            preOsc,afterOsc = preOsc*365./12.,afterOsc*365./12.
            if afterOsc*_signal_Eff < _s:
                break

    print 'result:',_standoff,afterOsc*_signal_Eff
    return _standoff, afterOsc,afterOsc*_signal_Eff
