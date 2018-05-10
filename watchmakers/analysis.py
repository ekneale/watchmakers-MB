from load import *
from io_operations import testEnabledCondition,writeResultsToFile

def fillHistograms(inFilePrefix,a1,t1,h,cover,ii,locj,covPCT):
    #   Obtain logarithmic binnings
    nbins, xbins, ybins = logx_logy_array()

    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)

    #fiducialVolume = float(arguments["--fv"])
    #pmtDist         = float(arguments["--psup"])
    # pmtDist = (float(arguments["--halfHeight"]))-(float(arguments["--shieldThick"]))
    pmtDist  = float(arguments["--tankRadius"])-float(arguments['--steelThick'])-float(arguments["--shieldThick"])
    pmtDistZ = float(arguments["--halfHeight"])-float(arguments['--steelThick'])-float(arguments["--shieldThick"])
    fiducialHalfHeight = float(arguments['--halfHeight'])-float(arguments['--steelThick'])-float(arguments['--shieldThick'])-float(arguments['--fidThick'])
    fiducialRadius     = float(arguments['--tankRadius'])-float(arguments['--steelThick'])-float(arguments['--shieldThick'])-float(arguments['--fidThick'])
    print "First change in analysis.py"



    #    print "Fiducial volume is ", fiducialVolume
    #   Read-in file
    try:
        s =  "ntuple_root_files%s/%s_%s_%s_%s.root"%(additionalString,inFilePrefix,ii,cover,locj)
        print "Reading in ",s
        #        data = root2array(s)
        t           = root2rec(s)

        #Apply some analysis
        r           = npa(t.reco_r<fiducialRadius,dtype=bool)
        z           = npa(absolute(t.reco_z)<fiducialHalfHeight,dtype=bool)

        #isFV        = logical_and(r,z,dtype=bool)
        isFV        = npa(t.FV==1,dtype=bool)
        notFV       = npa(t.FV!=1,dtype=bool)
        isFV_t      = npa(t.FV_truth==1,dtype=bool)
        notFV_t     = npa(t.FV_truth!=1,dtype=bool)

        iCandidate  = npa(t.candidate==1,dtype=bool)
        totalEvtGen = npa(t.all_ev==t.all_ev_tot,dtype=bool)

        tot         = float(sum(totalEvtGen))
        totD        = float(len(isFV))

        si          = logical_and(isFV,notFV_t,dtype=bool)
        so          = logical_and(notFV,isFV_t,dtype=bool)
        ei          = logical_and(isFV,isFV_t,dtype=bool)
        eo          = logical_and(notFV,notFV_t,dtype=bool)

        ssi         = sum(si)
        sso         = sum(so)
        sei         = sum(ei)
        seo         = sum(eo)
        print 'MC events generated',tot,', number of MC recorded',len(t.FV)
        print ssi,sso,sei,seo,', absolute: ', ssi/tot,sso/tot,sei/tot,seo/tot, \
        'relative',ssi/totD,sso/totD,sei/totD,seo/totD

        #Define set of histograms to fill out
        for subE in range(5):
            if subE ==0:
                subEv0   = npa(t.candidate==1,dtype=bool)
            else:
                subEv0   = npa(t.detected_ev==subE,dtype=bool)
            suma = sum(subEv0)
            if suma:
                print 'Sub event ',subE,' : ',suma
                totRel = float(suma)
                mask = logical_and(subEv0,si,dtype=bool)
                if sum(mask):
                    type = "si"
                    print subE,type,totRel,'(',sum(mask),')'

                    string = "%s_%s_%s_%d_abs"%(type,ii,locj,subE)
                    if string in h:
                        point = h[string].GetN()
                        h[string].SetPoint(point,covPCT,sum(mask)/tot)
                    else:
                        h[string] = Graph()
                        h[string].SetName(string)
                        h[string].SetPoint(0,covPCT,sum(mask)/tot)
                    string = "%s_%s_%s_%d_rel"%(type,ii,locj,subE)
                    if string in h:
                        point = h[string].GetN()
                        h[string].SetPoint(point,covPCT,sum(mask)/totRel)
                    else:
                        h[string] = Graph()
                        h[string].SetName(string)
                        h[string].SetPoint(0,covPCT,sum(mask)/totRel)

                    anaVar      = t.pe[mask,...]
                    s_dl        = "%s_pe_%s_%s_%s_%d"%(type,cover,ii,locj,subE)
                    h[s_dl]     = Hist(5000,0,500,name=s_dl,title=s_dl)
                    h[s_dl].SetXTitle('photoelectrons')
                    h[s_dl].SetYTitle('counts [a.u.]')
                    h[s_dl].fill_array(anaVar)

                    anaVar      = t.nhit[mask,...]
                    s_dl        ="%s_nhit_%s_%s_%s_%d"%(type,cover,ii,locj,subE)
                    h[s_dl]     = Hist(5000,0,500,name=s_dl,title=s_dl)
                    h[s_dl].SetXTitle('photoelectrons')
                    h[s_dl].SetYTitle('counts [a.u.]')
                    h[s_dl].fill_array(anaVar)

                    i_d         = t.inner_dist[mask,...]
                    i_t         = t.inner_time[mask,...]
                    anaVar      = column_stack((i_d,i_t/1e9))
                    s_dl       ="%s_dDdT_%s_%s_%s_%d"%(type,cover,ii,locj,subE)
                    h[s_dl]     = Hist2D(xbins,ybins,name=s_dl,title=s_dl);
                    h[s_dl].SetXTitle('inter-event distance [m]')
                    h[s_dl].SetYTitle('inter-event time [s]')
                    h[s_dl].fill_array(anaVar)
                    r_r         = t.reco_r[mask,...]
                    r_z         = t.reco_z[mask,...]
                    anaVar      = column_stack((power(r_r/pmtDist,2),r_z/pmtDistZ))
                    s_dl       = "%s_dRdZ_%s_%s_%s_%d"%(type,cover,ii,locj,subE)
                    h[s_dl]= Hist2D(100,0,1.6,200,-1.6,1.6,name=s_dl,title=s_dl)
                    h[s_dl].SetXTitle('(r / R_{pmt})^{2} [unitless]')
                    h[s_dl].SetYTitle('(z / Z_{pmt}) [unitless]')
                    h[s_dl].fill_array(anaVar)
                    t_r         = t.true_r[mask,...]
                    t_z         = t.true_z[mask,...]
                    anaVar      = column_stack((power(t_r/pmtDist,2),t_z/pmtDistZ))
                    s_dl       = "%s_tRtZ_%s_%s_%s_%d"%(type,cover,ii,locj,subE)
                    h[s_dl]= Hist2D(100,0,1.6,200,-1.6,1.6,name=s_dl,title=s_dl)
                    h[s_dl].SetXTitle('(r / R_{pmt})^{2} [unitless]')
                    h[s_dl].SetYTitle('(z / Z_{pmt}) [unitless]')
                    h[s_dl].fill_array(anaVar)




                mask = logical_and(subEv0,so)
                if sum(mask):

                    type        = "so"
                    print subE,type,totRel,'(',sum(mask),')'

                    string = "%s_%s_%s_%d_abs"%(type,ii,locj,subE)
                    if string in h:
                        point = h[string].GetN()
                        h[string].SetPoint(point,covPCT,sum(mask)/tot)
                    else:
                        h[string] = Graph()
                        h[string].SetName(string)
                        h[string].SetPoint(0,covPCT,sum(mask)/tot)
                    string = "%s_%s_%s_%d_rel"%(type,ii,locj,subE)
                    if string in h:
                        point = h[string].GetN()
                        h[string].SetPoint(point,covPCT,sum(mask)/totRel)
                    else:
                        h[string] = Graph()
                        h[string].SetName(string)
                        h[string].SetPoint(0,covPCT,sum(mask)/totRel)
                    anaVar      = t.pe[mask,...]
                    s_dl        = "%s_pe_%s_%s_%s_%d"%(type,cover,ii,locj,subE)
                    h[s_dl]     = Hist(5000,0,500,name=s_dl,title=s_dl)
                    h[s_dl].SetXTitle('photoelectrons')
                    h[s_dl].SetYTitle('counts [a.u.]')
                    h[s_dl].fill_array(anaVar)

                    anaVar      = t.nhit[mask,...]
                    s_dl        ="%s_nhit_%s_%s_%s_%d"%(type,cover,ii,locj,subE)
                    h[s_dl]     = Hist(5000,0,500,name=s_dl,title=s_dl)
                    h[s_dl].SetXTitle('photoelectrons')
                    h[s_dl].SetYTitle('counts [a.u.]')
                    h[s_dl].fill_array(anaVar)


                    i_d         = t.inner_dist[mask,...]
                    i_t         = t.inner_time[mask,...]
                    anaVar      = column_stack((i_d,i_t/1e9))
                    s_dl        ="%s_dDdT_%s_%s_%s_%d"%(type,cover,ii,locj,subE)
                    h[s_dl]     = Hist2D(xbins,ybins,name=s_dl,title=s_dl);
                    h[s_dl].SetXTitle('inter-event distance [m]')
                    h[s_dl].SetYTitle('inter-event time [s]')
                    h[s_dl].fill_array(anaVar)
                    r_r         = t.reco_r[mask,...]
                    r_z         = t.reco_z[mask,...]

                    anaVar      = column_stack((power(r_r/pmtDist,2),r_z/pmtDistZ))
                    s_dl        ="%s_dRdZ_%s_%s_%s_%d"%(type,cover,ii,locj,subE)
                    h[s_dl]= Hist2D(100,0,1.6,200,-1.6,1.6,name=s_dl,title=s_dl)
                    h[s_dl].SetXTitle('(r / R_{pmt})^{2} [unitless]')
                    h[s_dl].SetYTitle('(z / Z_{pmt}) [unitless]')
                    h[s_dl].fill_array(anaVar)
                    t_r         = t.true_r[mask,...]
                    t_z         = t.true_z[mask,...]

                    anaVar      = column_stack((power(t_r/pmtDist,2),t_z/pmtDistZ))
                    s_dl       = "%s_tRtZ_%s_%s_%s_%d"%(type,cover,ii,locj,subE)
                    h[s_dl]= Hist2D(100,0,1.6,200,-1.6,1.6,name=s_dl,title=s_dl)
                    h[s_dl].SetXTitle('(r / R_{pmt})^{2} [unitless]')
                    h[s_dl].SetYTitle('(z / Z_{pmt}) [unitless]')
                    h[s_dl].fill_array(anaVar)


                mask = logical_and(subEv0,ei)
                if sum(mask):

                    type        = "ei"
                    print subE,type,totRel,'(',sum(mask),')'


                    string = "%s_%s_%s_%d_abs"%(type,ii,locj,subE)
                    if string in h:
                        point = h[string].GetN()
                        h[string].SetPoint(point,covPCT,sum(mask)/tot)
                    else:
                        h[string] = Graph()
                        h[string].SetName(string)
                        h[string].SetPoint(0,covPCT,sum(mask)/tot)
                    string = "%s_%s_%s_%d_rel"%(type,ii,locj,subE)
                    if string in h:
                        point = h[string].GetN()
                        h[string].SetPoint(point,covPCT,sum(mask)/totRel)
                    else:
                        h[string] = Graph()
                        h[string].SetName(string)
                        h[string].SetPoint(0,covPCT,sum(mask)/totRel)
                    anaVar      = t.pe[mask,...]
                    s_dl        = "%s_pe_%s_%s_%s_%d"%(type,cover,ii,locj,subE)
                    h[s_dl]     = Hist(5000,0,500,name=s_dl,title=s_dl)
                    h[s_dl].SetXTitle('photoelectrons')
                    h[s_dl].SetYTitle('counts [a.u.]')
                    h[s_dl].fill_array(anaVar)


                    anaVar      = t.nhit[mask,...]
                    s_dl= "%s_nhit_%s_%s_%s_%d"%(type,cover,ii,locj,subE)
                    h[s_dl]     = Hist(5000,0,500,name=s_dl,title=s_dl)
                    h[s_dl].SetXTitle('photoelectrons')
                    h[s_dl].SetYTitle('counts [a.u.]')
                    h[s_dl].fill_array(anaVar)

                    i_d         = t.inner_dist[mask,...]
                    i_t         = t.inner_time[mask,...]
                    anaVar      = column_stack((i_d,i_t/1e9))
                    s_dl= "%s_dDdT_%s_%s_%s_%d"%(type,cover,ii,locj,subE)
                    h[s_dl]     = Hist2D(xbins,ybins,name=s_dl,title=s_dl);
                    h[s_dl].SetXTitle('inter-event distance [m]')
                    h[s_dl].SetYTitle('inter-event time [s]')
                    h[s_dl].fill_array(anaVar)
                    r_r         = t.reco_r[mask,...]
                    r_z         = t.reco_z[mask,...]
                    anaVar      = column_stack((power(r_r/pmtDist,2),r_z/pmtDistZ))
                    s_dl= "%s_dRdZ_%s_%s_%s_%d"%(type,cover,ii,locj,subE)
                    h[s_dl]= Hist2D(100,0,1.6,200,-1.6,1.6,name=s_dl,title=s_dl)
                    h[s_dl].SetXTitle('(r / R_{pmt})^{2} [unitless]')
                    h[s_dl].SetYTitle('(z / Z_{pmt}) [unitless]')
                    h[s_dl].fill_array(anaVar)

                    t_r         = t.true_r[mask,...]
                    t_z         = t.true_z[mask,...]

                    anaVar      = column_stack((power(t_r/pmtDist,2),t_z/pmtDistZ))
                    s_dl       = "%s_tRtZ_%s_%s_%s_%d"%(type,cover,ii,locj,subE)
                    h[s_dl]= Hist2D(100,0,1.6,200,-1.6,1.6,name=s_dl,title=s_dl)
                    h[s_dl].SetXTitle('(r / R_{pmt})^{2} [unitless]')
                    h[s_dl].SetYTitle('(z / Z_{pmt}) [unitless]')
                    h[s_dl].fill_array(anaVar)

                mask = logical_and(subEv0,eo)
                if sum(mask):

                    type        = "eo"
                    print subE,type,totRel,'(',sum(mask),')'
                    string = "%s_%s_%s_%d_abs"%(type,ii,locj,subE)
                    if string in h:
                        point = h[string].GetN()
                        h[string].SetPoint(point,covPCT,sum(mask)/tot)
                    else:
                        h[string] = Graph()
                        h[string].SetName(string)
                        h[string].SetPoint(0,covPCT,sum(mask)/tot)
                    string = "%s_%s_%s_%d_rel"%(type,ii,locj,subE)
                    if string in h:
                        point = h[string].GetN()
                        h[string].SetPoint(point,covPCT,sum(mask)/totRel)
                    else:
                        h[string] = Graph()
                        h[string].SetName(string)
                        h[string].SetPoint(0,covPCT,sum(mask)/totRel)

                    anaVar      = t.pe[mask,...]
                    s_dl         = "%s_pe_%s_%s_%s_%d"%(type,cover,ii,locj,subE)
                    h[s_dl]     = Hist(5000,0,500,name=s_dl,title=s_dl)
                    h[s_dl].SetXTitle('photoelectrons')
                    h[s_dl].SetYTitle('counts [a.u.]')
                    h[s_dl].fill_array(anaVar)

                    anaVar      = t.nhit[mask,...]
                    s_dl        ="%s_nhit_%s_%s_%s_%d"%(type,cover,ii,locj,subE)
                    h[s_dl]     = Hist(5000,0,500,name=s_dl,title=s_dl)
                    h[s_dl].SetXTitle('photoelectrons')
                    h[s_dl].SetYTitle('counts [a.u.]')
                    h[s_dl].fill_array(anaVar)

                    i_d         = t.inner_dist[mask,...]
                    i_t         = t.inner_time[mask,...]
                    anaVar      = column_stack((i_d,i_t/1e9))
                    s_dl        ="%s_dDdT_%s_%s_%s_%d"%(type,cover,ii,locj,subE)
                    h[s_dl]= Hist2D(xbins,ybins,name=s_dl,title=s_dl);
                    h[s_dl].SetXTitle('inter-event distance [m]')
                    h[s_dl].SetYTitle('inter-event time [s]')
                    h[s_dl].fill_array(anaVar)
                    r_r         = t.reco_r[mask,...]
                    r_z         = t.reco_z[mask,...]
                    anaVar      = column_stack((power(r_r/pmtDist,2),r_z/pmtDistZ))
                    s_dl       = "%s_dRdZ_%s_%s_%s_%d"%(type,cover,ii,locj,subE)
                    h[s_dl]= Hist2D(100,0,1.6,200,-1.6,1.6,name=s_dl,title=s_dl)
                    h[s_dl].SetXTitle('(r / R_{pmt})^{2} [unitless]')
                    h[s_dl].SetYTitle('(z / Z_{pmt}) [unitless]')
                    h[s_dl].fill_array(anaVar)
                    t_r         = t.true_r[mask,...]
                    t_z         = t.true_z[mask,...]
                    anaVar      = column_stack((power(t_r/pmtDist,2),t_z/pmtDistZ))
                    s_dl       = "%s_tRtZ_%s_%s_%s_%d"%(type,cover,ii,locj,subE)
                    h[s_dl]= Hist2D(100,0,1.6,200,-1.6,1.6,name=s_dl,title=s_dl)
                    h[s_dl].SetXTitle('(r / R_{pmt})^{2} [unitless]')
                    h[s_dl].SetYTitle('(z / Z_{pmt}) [unitless]')
                    h[s_dl].fill_array(anaVar)


    except:
        print "Could not read file ",s

    print ""
    #        f.Close()
    return h

def extractHistogramWitCorrectRate():

    g,h = {},{}

    d,iso,loc,coverage,coveragePCT = loadSimulationParameters()

    # Obtain logarithmic binnings
    nbins, xbins, ybins = logx_logy_array()

    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)

    boolSUPERNOVA_FORMAT    = arguments["--supernovaFormat"]

    #fiducialVolume          = float(arguments["--fv"])
    fiducialHalfHeight= float(arguments['--halfHeight'])- float(arguments['--steelThick'])- float(arguments['--shieldThick'])- float(arguments['--fidThick'])
    fiducialRadius    = float(arguments['--tankRadius'])- float(arguments['--steelThick'])- float(arguments['--shieldThick'])- float(arguments['--fidThick'])

    #pmtDist                 = float(arguments["--psup"])
    pmtDist  = float(arguments["--tankRadius"])-float(arguments['--steelThick'])-float(arguments["--shieldThick"])
    pmtDistZ = float(arguments["--halfHeight"])-float(arguments['--steelThick'])-float(arguments["--shieldThick"])
    print "Values of Fid. Vol. (HalfHeight,Radius) and PMT placement (HalfHeight,Radius) ",fiducialHalfHeight,fiducialRadius,pmtDist,pmtDistZ, ' mm'

    timeScale               = arguments["--timeScale"]
    inFilePrefix            = arguments["--ft"]
    timeCut                 = float(arguments["-t"])*1e3
    distCut                 = float(arguments["-d"])

    parameters  = loadAnalysisParameters(timeScale)
    rates       = parameters[11]
    mass        = parameters[10]
    pc_num      = parameters[12]
    pc_val      = parameters[13]
    timeS       = parameters[8]
    timeRedux     = float(arguments["-t"]) # cut is in microsecond
    timeRedux     *=1e-6/timeS
    print "\nThe raw rate of events per %s are"%(timeScale)
    _rates = sorted(rates,key=rates.__getitem__, reverse=True)
    for _r in _rates:
        if 'PMT' in _r:
            print "%25s %5.3e per %s per PMT"%(_r,rates[_r],timeScale)
        else:
            print "%25s %5.3e per %s"%(_r,rates[_r],timeScale)
    print "The cuts are : time window %4.3e %s" %(timeRedux,timeScale)
    #    print "Fiducial volume is ", fiducialVolume
    #Read-in file

    #inFilePrefix,h,cover,ii,locj,covPCT


    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)
    _str = "ntuple_root_files%s/%s" %(additionalString,arguments["-o"])
    f_root = TFile(_str,"recreate")

    # First find the PE per MeV
    procConsidered = ['boulby','imb']
    locj           = 'S'
    PEMEV    = {}

    string  = "pePerMeV_boulby"
    g[string] = TGraph()
    string  = "nhitPerMeV_boulby"
    g[string] = TGraph()
    string  = "n9PerMeV_boulby"
    g[string] = TGraph()
    cntB    = 0

    string  = "pePerMeV_imb"
    g[string] = TGraph()
    string  = "nhitPerMeV_imb"
    g[string] = TGraph()
    string  = "n9PerMeV_imb"
    g[string] = TGraph()

    cntF    = 0

    for ii in procConsidered:
        for idx,cover in enumerate(coverage):
            covPCT  = coveragePCT[cover]
            s =  "ntuple_root_files%s/%s_%s_%s_%s.root"%(additionalString,inFilePrefix,ii,cover,locj)
            print "\nEvaluating pe/MeV in ",s
            rfile = TFile(s)
            t   = rfile.Get('data')

            recoFVstring    = "(sqrt(pow(posReco.X(),2) + pow(posReco.Y(),2))<%f && sqrt(pow(posReco.Z(),2))<%f)"%(fiducialRadius,fiducialHalfHeight)
            trueFVstring    = "(sqrt(pow(posTruth.X(),2) + pow(posTruth.Y(),2))<%f && sqrt(pow(posTruth.Z(),2))<%f)"%(fiducialRadius,fiducialHalfHeight)
            posGood        = "(pos_goodness>%f)" %(float(arguments["-g"]))

            backgroundNoise = 1.0e3*float(pc_num["%s"%(cover)])*1500.*1e-9

            fPolyFit1  = TF1('fPolyFit1',"[0]+[1]*x",0.0,10.0)
            fPolyFit1.SetParameters(backgroundNoise,15.)
            s_nhitPerMeV_eisi        = "%s_nhitPerMeV_%s_%s_%s_%d"%('si',cover,ii,locj,1)
            if boolSUPERNOVA_FORMAT:
                t.Draw("pe:mc_energy>>h%s(100,0,10,200,0,200)"%(s_nhitPerMeV_eisi),"%s && %s "%(recoFVstring,posGood),"goff")
            if not boolSUPERNOVA_FORMAT:
                t.Draw("nhit:mc_prim_energy>>h%s(100,0,10,200,0,200)"%(s_nhitPerMeV_eisi),"%s && %s "%(recoFVstring,posGood),"goff")

            h1 = t.GetHistogram()
            h[s_nhitPerMeV_eisi] = h1.ProfileX()
            h[s_nhitPerMeV_eisi].Fit("fPolyFit1","MREQ","",2,6.5)
            fitRes1 = h[s_nhitPerMeV_eisi].GetFunction("fPolyFit1")
            print ' nhit results of fit :',fitRes1.GetParameter(0),fitRes1.GetParameter(1)

            if ii == 'boulby':
                g["nhitPerMeV_boulby"].SetPoint(cntB,pc_val["%s"%(cover)],fitRes1.GetParameter(1))
                g["nhitPerMeV_boulby"].SetName("nhitPerMeV_boulby")
                g["nhitPerMeV_boulby"].GetXaxis().SetTitle('PMT coverage')
                g["nhitPerMeV_boulby"].GetYaxis().SetTitle('pe / MeV')

            if ii == 'imb':
                g["nhitPerMeV_imb"].SetPoint(cntB,pc_val["%s"%(cover)],fitRes1.GetParameter(1))
                g["nhitPerMeV_imb"].SetName("nhitPerMeV_imb")
                g["nhitPerMeV_imb"].GetXaxis().SetTitle('PMT coverage')
                g["nhitPerMeV_imb"].GetYaxis().SetTitle('pe / MeV')



            fPolyFit2  = TF1('fPolyFit2',"[0]+[1]*x",0.0,10.0)
            fPolyFit2.SetParameters(backgroundNoise,15.)
            s_pePerMeV_eisi        = "%s_pePerMeV_%s_%s_%s_%d"%('si',cover,ii,locj,1)
            if boolSUPERNOVA_FORMAT:
                t.Draw("pe:mc_energy>>h%s(100,0,10,200,0,200)"%(s_pePerMeV_eisi),"%s && %s "%(recoFVstring,posGood),"goff")
            if not boolSUPERNOVA_FORMAT:
                t.Draw("pe:mc_prim_energy>>h%s(100,0,10,200,0,200)"%(s_pePerMeV_eisi),"%s && %s "%(recoFVstring,posGood),"goff")

            h1 = t.GetHistogram()
            h[s_pePerMeV_eisi] = h1.ProfileX()
            h[s_pePerMeV_eisi].Fit("fPolyFit2","MREQ","",2,6.5)
            fitRes2 = h[s_pePerMeV_eisi].GetFunction("fPolyFit2")
            print ' pe results of fit :',fitRes2.GetParameter(0),fitRes2.GetParameter(1)
            PEMEV['%s'%(cover)] = fitRes2.GetParameter(1)

            if ii == 'boulby':
                g["pePerMeV_boulby"].SetPoint(cntB,pc_val["%s"%(cover)],fitRes2.GetParameter(1))
                g["pePerMeV_boulby"].SetName("pePerMeV_boulby")
                g["pePerMeV_boulby"].GetXaxis().SetTitle('PMT coverage')
                g["pePerMeV_boulby"].GetYaxis().SetTitle('pe / MeV')

            if ii == 'imb':
                g["pePerMeV_imb"].SetPoint(cntF,pc_val["%s"%(cover)],fitRes2.GetParameter(1))
                g["pePerMeV_imb"].SetName("pePerMeV_imb")
                g["pePerMeV_imb"].GetXaxis().SetTitle('PMT coverage')
                g["pePerMeV_imb"].GetYaxis().SetTitle('pe / MeV')


            fPolyFit3  = TF1('fPolyFit3',"[0]+[1]*x",0.0,10.0)
            fPolyFit3.SetParameters(backgroundNoise,15.)
            s_n9PerMeV_eisi        = "%s_n9PerMeV_%s_%s_%s_%d"%('si',cover,ii,locj,1)
            if boolSUPERNOVA_FORMAT:
                t.Draw("n9:mc_energy>>h%s(100,0,10,200,0,200)"%(s_n9PerMeV_eisi),"%s && %s"%(recoFVstring,posGood),"goff")
            if not boolSUPERNOVA_FORMAT:
                t.Draw("n9:mc_prim_energy>>h%s(100,0,10,200,0,200)"%(s_n9PerMeV_eisi),"%s && %s"%(recoFVstring,posGood),"goff")


            h1 = t.GetHistogram()
            h[s_n9PerMeV_eisi] = h1.ProfileX()
            h[s_n9PerMeV_eisi].Fit("fPolyFit3","MREQ","",2,6.5)
            fitRes3 = h[s_n9PerMeV_eisi].GetFunction("fPolyFit3")
            print ' n9 results of fit :',fitRes3.GetParameter(0),fitRes3.GetParameter(1)

            if ii == 'boulby':
                g["n9PerMeV_boulby"].SetPoint(cntB,pc_val["%s"%(cover)],fitRes3.GetParameter(1))
                g["n9PerMeV_boulby"].SetName("n9PerMeV_boulby")
                g["n9PerMeV_boulby"].GetXaxis().SetTitle('PMT coverage')
                g["n9PerMeV_boulby"].GetYaxis().SetTitle('pe / MeV')
                cntB+=1
            if ii == 'imb':
                g["n9PerMeV_imb"].SetPoint(cntB,pc_val["%s"%(cover)],fitRes3.GetParameter(1))
                g["n9PerMeV_imb"].SetName("n9PerMeV_imb")
                g["n9PerMeV_imb"].GetXaxis().SetTitle('PMT coverage')
                g["n9PerMeV_imb"].GetYaxis().SetTitle('pe / MeV')
                cntF+=1

            f_root.cd()
            h[s_nhitPerMeV_eisi].Write()
            h[s_pePerMeV_eisi].Write()
            h[s_n9PerMeV_eisi].Write()



    print PEMEV
    f_root.cd()
    g["pePerMeV_boulby"].Write()
    g["pePerMeV_imb"].Write()
    g["nhitPerMeV_boulby"].Write()
    g["nhitPerMeV_imb"].Write()
    g["n9PerMeV_boulby"].Write()
    g["n9PerMeV_imb"].Write()

    rawTotalRateEISI = {}
#    rawTotalRateEISI["10pct"] = 0.0
#    rawTotalRateEISI["15pct"] = 0.0
    rawTotalRateEISI["20pct"] = 0.0
    rawTotalRateEISI["25pct"] = 0.0
    rawTotalRateEISI["30pct"] = 0.0
#    rawTotalRateEISI["35pct"] = 0.0
#    rawTotalRateEISI["40pct"] = 0.0

    rawTotalRateTot = {}
#    rawTotalRateTot["10pct"] = 0.0
#    rawTotalRateTot["15pct"] = 0.0
    rawTotalRateTot["20pct"] = 0.0
    rawTotalRateTot["25pct"] = 0.0
    rawTotalRateTot["30pct"] = 0.0
#    rawTotalRateTot["35pct"] = 0.0
#    rawTotalRateTot["40pct"] = 0.0


    for j in range(len(iso)):
        for ii in d["%s"%(iso[int(j)])]:
            locj    = loc[j]

            string  = "si_%s_%s_1_abs" %(ii,locj)
            g[string] = TGraph()
            g[string].SetName(string)
            string  = "ei_%s_%s_1_abs" %(ii,locj)
            g[string] = TGraph()
            g[string].SetName(string)
            string  = "so_%s_%s_1_abs" %(ii,locj)
            g[string] = TGraph()
            g[string].SetName(string)
            string  = "eo_%s_%s_1_abs" %(ii,locj)
            g[string] = TGraph()
            g[string].SetName(string)


            string  = "si_%s_%s_singlesRate" %(ii,locj)
            g[string] = TGraph()
            g[string].SetName(string)
            string  = "ei_%s_%s_singlesRate" %(ii,locj)
            g[string] = TGraph()
            g[string].SetName(string)
            string  = "so_%s_%s_singlesRate" %(ii,locj)
            g[string] = TGraph()
            g[string].SetName(string)
            string  = "eo_%s_%s_singlesRate" %(ii,locj)
            g[string] = TGraph()
            g[string].SetName(string)


            string  = "eisi_%s_%s_singlesRate_ProxCut" %(ii,locj)
            g[string] = TGraph()
            g[string].SetName(string)

            string  = "eisi_%s_%s_abs_ProxCut" %(ii,locj)
            g[string] = TGraph()
            g[string].SetName(string)
            cntG    = 0

            string  = "eisi_%s_%s_singlesRate_ProxCut_TimeCut" %(ii,locj)
            g[string] = TGraph()
            g[string].SetName(string)

            for idx,cover in enumerate(coverage):
                covPCT  = coveragePCT[cover]
                # rawTotalRateEISI["%s"%(cover)] = 0.0

                try:
                # if 1==1:
                    #                        branches = 'pe','nhit','n9','delta_time_s', 'detected_ev','detected_ev_tot','all_ev','all_ev_tot',subevents,event_number,candidate,mc_prim_energy,pos_goodness,posReco,reco_r,reco_z,posTruth,true_r,true_z,dir_goodness,dirReco,dirPrimaryMC,FV,GSV,EV,OV,IV,FV_truth,GSV_truth,EV_truth,OV_truth,IV_truth,inner_dist,inner_time,inner_dist_fv,tot_FV,consecutive_FV')'
                    s =  "ntuple_root_files%s/%s_%s_%s_%s.root"%(additionalString,inFilePrefix,ii,cover,locj)
                    print "\nReading in ",s

                    rfile = TFile(s)
                    t   = rfile.Get('data')

                    recoFVstring    = "(sqrt(pow(posReco.X(),2) + pow(posReco.Y(),2))<%f && sqrt(pow(posReco.Z(),2))<%f)"%(fiducialRadius,fiducialHalfHeight)
                    trueFVstring    = "(sqrt(pow(posTruth.X(),2) + pow(posTruth.Y(),2))<%f && sqrt(pow(posTruth.Z(),2))<%f)"%(fiducialRadius,fiducialHalfHeight)
                    posGood        = "(pos_goodness>%f)" %(0.5)
                    n9Good = "(n9>%f)" %(float(arguments["--minN9"]))

                    eiCond = " %s &&  %s && %s && %s" %(recoFVstring,\
                    trueFVstring,posGood,n9Good)
                    siCond = " %s && !%s && %s && %s" %(recoFVstring,\
                    trueFVstring,posGood,n9Good)
                    eoCond = "!%s && !%s && %s && %s" %(recoFVstring,\
                    trueFVstring,posGood,n9Good)
                    soCond = "!%s &&  %s && %s && %s" %(recoFVstring,\
                    trueFVstring,posGood,n9Good)
                    # print eiCond,siCond
                    # Add additional condition for correlated MC
                    if not boolSUPERNOVA_FORMAT:
                        multCond  = "&& detected_ev==1 && detected_ev_tot>=2"
                        multCond  += "&& tot_FV==2"
                        multCondD  = "&& detected_ev==2 && detected_ev_tot>=2"
                        multCondD += "&& tot_FV==2"
                        multCondD += "&& inner_time<%f " %(timeCut)
                        multCondD += "&& inner_dist_fv<%f" %(distCut)
                        totEvtStr  = "all_ev_tot == all_ev"
                    if boolSUPERNOVA_FORMAT:
                        multCond   = "&& sub_ev ==1 && sub_ev_cnt >=2"
                        multCond  += "&& !single"
                        multCondD  = "&& sub_ev==2 && sub_ev>=2"
                        multCondD += "&& !single"
                        multCondD += "&& local_time_ns<%f " %(timeCut)
                        # multCondD += "&& inner_dist_fv<%f" %(distCut) No such value
                        totEvtStr  = "sub_ev_cnt == sub_ev"

                    tt = t.Draw("pe",totEvtStr,"goff")

                    er = float(rates["%s_%s"%(ii,locj)])

                    if locj == 'PMT':
                        print " adjusting rate for pmt %4.3e %s^{-1}, pmt mass %4.2f, number of PMTs %d : %4.3e %s^{-1}"%(er,timeScale,mass,pc_num["%s"%(cover)],er*pc_num["%s"%(cover)]*mass,timeScale)
                        er*=pc_num["%s"%(cover)]*mass
                    print ' Total raw event rate is ',\
                    er, ' per ', timeScale

                    # Evalute pe distribution (used in sensitivity map)
                    # Fast Neutrons (FN) and Inverse Beta Decay (IBD)
                    # process have correlated part simulated
                    if locj == 'FN' or locj == 'I':
                        s_dl    = "%s_pe_%s_%s_%s_%d"%('ei',cover,ii,locj,1)
                        # h[s_dl]     = TH1D(s_dl,s_dl,5000,0,500)
                        drawArg = "n9>>hei(5000,0,500)"
                        ei1 = t.Draw(drawArg,eiCond,"goff")
                        ei  = t.Draw(drawArg,eiCond+multCondD,"goff")
                        eiP = t.Draw(drawArg,eiCond+multCond,"goff")
                        print locj,ei1,eiP,ei
                        h[s_dl] = t.GetHistogram()
                        h[s_dl].SetName(s_dl)
                        h[s_dl].SetXTitle('photoelectrons')
                        h[s_dl].SetYTitle('counts for %d MC events'%(tt))

                        s_dlsi  = "%s_pe_%s_%s_%s_%d"%('si',cover,ii,locj,1)
                        drawArg = "n9>>hsi(5000,0,500)"
                        si1 = t.Draw(drawArg,siCond,"goff")
                        si = t.Draw(drawArg,siCond+multCondD,"goff")
                        siP = t.Draw(drawArg,siCond+multCond,"goff")
                        print locj,si1,siP,si
                        h[s_dlsi] = t.GetHistogram()
                        h[s_dlsi].SetName(s_dlsi)
                        h[s_dlsi].SetXTitle('photoelectrons')
                        h[s_dlsi].SetYTitle('counts for %d MC events'%(tt))
                    else:
                        s_dl    = "%s_pe_%s_%s_%s_%d"%('ei',cover,ii,locj,1)
                        # h[s_dl]     = TH1D(s_dl,s_dl,5000,0,500)
                        drawArg = "n9>>hei(5000,0,500)"
                        ei = t.Draw(drawArg,eiCond,"goff")
                        h[s_dl] = t.GetHistogram()
                        h[s_dl].SetName(s_dl)
                        h[s_dl].SetXTitle('photoelectrons')
                        h[s_dl].SetYTitle('counts for %d MC events'%(tt))

                        s_dlsi  = "%s_pe_%s_%s_%s_%d"%('si',cover,ii,locj,1)
                        drawArg = "n9>>hsi(5000,0,500)"
                        si = t.Draw(drawArg,siCond,"goff")
                        h[s_dlsi] = t.GetHistogram()
                        h[s_dlsi].SetName(s_dlsi)
                        h[s_dlsi].SetXTitle('photoelectrons')
                        h[s_dlsi].SetYTitle('counts for %d MC events'%(tt))
                    # Estimate the reconstructed energy from pe/MeV fits
                    s_dlT     = "%s_Teff_%s_%s_%s_%d"%('ei',cover,ii,locj,1)
                    # h[s_dlT]     = TH1D(s_dlT,s_dlT,5000,0,500)
                    drawArg = "pe/%f>>heiT(2000,0,20)"%(PEMEV['%s'%(cover)])
                    eiT       = t.Draw(drawArg,eiCond,"goff")
                    h[s_dlT]  = t.GetHistogram()
                    h[s_dlT].SetName(s_dlT)
                    h[s_dlT].SetXTitle('T_{eff} [MeV]')
                    h[s_dlT].SetYTitle('counts [%s]^{-1}'%(timeScale))
                    h[s_dlT].Scale(ei/float(tt)*er)

                    s_dlTsi   = "%s_Teff_%s_%s_%s_%d"%('si',cover,ii,locj,1)
                    # h[s_dlTsi]     = TH1D(s_dlTsi,s_dlTsi,5000,0,500)
                    drawArg = "pe/%f>>hsiT(2000,0,20)"%(PEMEV['%s'%(cover)])
                    siT         = t.Draw(drawArg,siCond,"goff")
                    h[s_dlTsi] = t.GetHistogram()
                    h[s_dlTsi].SetName(s_dlTsi)
                    h[s_dlTsi].SetXTitle('T_{eff} [MeV]')
                    h[s_dlTsi].SetYTitle('counts [%s]^{-1}'%(timeScale))
                    h[s_dlTsi].Scale(si/float(tt)*er)

                    s_dlso    = "%s_Teff_%s_%s_%s_%d"%('so',cover,ii,locj,1)
                    drawArg =  "pe/%f>>hso(2000,0,20)"%(PEMEV['%s'%(cover)])
                    so = t.Draw(drawArg,soCond,"goff")
                    h[s_dlso] = t.GetHistogram()
                    h[s_dlso].SetName(s_dlso)
                    h[s_dlso].SetXTitle('photoelectrons')
                    h[s_dlso].SetYTitle('counts for %d MC events'%(tt))

                    s_dleo    = "%s_Teff_%s_%s_%s_%d"%('eo',cover,ii,locj,1)
                    drawArg =  "pe/%f>>heo(2000,0,20)"%(PEMEV['%s'%(cover)])
                    eo = t.Draw(drawArg,eoCond,"goff")
                    h[s_dleo] = t.GetHistogram()
                    h[s_dleo].SetName(s_dleo)
                    h[s_dleo].SetXTitle('photoelectrons')
                    h[s_dleo].SetYTitle('counts for %d MC events'%(tt))

                    print " (%7s %7s %7s %7s ) %8s |            (%9s %9s %9s %9s)" %('ei','si','so','eo','mc events','ei','si','so','eo')
                    print " (%7d %7d %7d %7d ) %8d  | efficiency (%4.3e %4.3e %4.3e %4.3e)" %(ei,si,so,eo, tt,ei/float(tt),si/float(tt),so/float(tt),eo/float(tt))
                    print " (%7d %7d %7d %7d ) %8d  | event rate (%4.3e %4.3e %4.3e %4.3e) per %s" % (ei,si,so,eo, tt,ei/float(tt)*er,si/float(tt)*er,so/float(tt)*er,eo/float(tt)*er,timeScale)
                    if locj == 'PMT' or locj == 'FV':
                        print " (%7d %7d %7d %7d ) %8d  | event rate (%4.3e %4.3e %4.3e %4.3e) per %s (after time redux)" % (ei,si,so,eo, tt,timeRedux*power(ei/float(tt)*er,2),timeRedux*power(si/float(tt)*er,2),timeRedux*power(so/float(tt)*er,2),timeRedux*power(eo/float(tt)*er,2),timeScale)
                        rawTotalRateEISI["%s"%(cover)] += si/float(tt)*er + ei/float(tt)*er
                        rawTotalRateTot["%s"%(cover)] += si/float(tt)*er + ei/float(tt)*er + so/float(tt)*er + eo/float(tt)*er
                    g["si_%s_%s_1_abs" %(ii,locj)].SetPoint(cntG,pc_val["%s"%(cover)],si/float(tt))
                    g["ei_%s_%s_1_abs" %(ii,locj)].SetPoint(cntG,pc_val["%s"%(cover)],ei/float(tt))
                    g["so_%s_%s_1_abs" %(ii,locj)].SetPoint(cntG,pc_val["%s"%(cover)],so/float(tt))
                    g["eo_%s_%s_1_abs" %(ii,locj)].SetPoint(cntG,pc_val["%s"%(cover)],eo/float(tt))
                    g["si_%s_%s_1_abs" %(ii,locj)].GetXaxis().SetTitle('PMT coverage')
                    g["ei_%s_%s_1_abs" %(ii,locj)].GetXaxis().SetTitle('PMT coverage')
                    g["so_%s_%s_1_abs" %(ii,locj)].GetXaxis().SetTitle('PMT coverage')
                    g["eo_%s_%s_1_abs" %(ii,locj)].GetXaxis().SetTitle('PMT coverage')
                    g["si_%s_%s_1_abs" %(ii,locj)].GetYaxis().SetTitle('efficiency')
                    g["ei_%s_%s_1_abs" %(ii,locj)].GetYaxis().SetTitle('efficiency')
                    g["so_%s_%s_1_abs" %(ii,locj)].GetYaxis().SetTitle('efficiency')
                    g["eo_%s_%s_1_abs" %(ii,locj)].GetYaxis().SetTitle('efficiency')

                    g["si_%s_%s_singlesRate" %(ii,locj)].SetPoint(cntG,pc_val["%s"%(cover)],si/float(tt)*er)
                    g["ei_%s_%s_singlesRate" %(ii,locj)].SetPoint(cntG,pc_val["%s"%(cover)],ei/float(tt)*er)
                    g["so_%s_%s_singlesRate" %(ii,locj)].SetPoint(cntG,pc_val["%s"%(cover)],so/float(tt)*er)
                    g["eo_%s_%s_singlesRate" %(ii,locj)].SetPoint(cntG,pc_val["%s"%(cover)],eo/float(tt)*er)
                    g["si_%s_%s_singlesRate" %(ii,locj)].GetXaxis().SetTitle('PMT coverage')
                    g["ei_%s_%s_singlesRate" %(ii,locj)].GetXaxis().SetTitle('PMT coverage')
                    g["so_%s_%s_singlesRate" %(ii,locj)].GetXaxis().SetTitle('PMT coverage')
                    g["eo_%s_%s_singlesRate" %(ii,locj)].GetXaxis().SetTitle('PMT coverage')
                    g["si_%s_%s_singlesRate" %(ii,locj)].GetYaxis().SetTitle('counts [%s]^{-1}'%(timeScale))
                    g["ei_%s_%s_singlesRate" %(ii,locj)].GetYaxis().SetTitle('counts [%s]^{-1}'%(timeScale))
                    g["so_%s_%s_singlesRate" %(ii,locj)].GetYaxis().SetTitle('counts [%s]^{-1}'%(timeScale))
                    g["eo_%s_%s_singlesRate" %(ii,locj)].GetYaxis().SetTitle('counts [%s]^{-1}'%(timeScale))

                    # See below
                    #                        cntG+=1

                    N  = t.Draw("pe:posReco.X():posReco.Y():posReco.Z()"," %s && %s && %s " %(recoFVstring,posGood,n9Good),"goff")
                    pe   = t.GetV1()
                    x    = t.GetV2()
                    y    = t.GetV3()
                    z    = t.GetV4()

                    s_dl1        = "%s_dD_%s_%s_%s_%s"%('eisi',cover,ii,locj,'proxCutEfficiency')
                    h[s_dl1]     = TH1D(s_dl1,s_dl1,5000,0,15)
                    h[s_dl1].SetName(s_dl1)
                    h[s_dl1].SetXTitle('distance')
                    h[s_dl1].SetYTitle('efficiency')

                    s_dl2       ="%s_pep_ped_%s_%s_%s_%s"%('eisi',cover,ii,locj,'proxCutEfficiency')
                    h[s_dl2]     = TH2D(s_dl2,s_dl2,200,0,200,200,0,200)
                    h[s_dl2].SetXTitle('prompt energy [pe]')
                    h[s_dl2].SetYTitle('delayed energy [pe]')
                    h[s_dl2].SetZTitle('efficiency')


                    s_dl1b        = "%s_dD_%s_%s_%s_%s"%('eisi',cover,ii,locj,'proxCutRate')
                    h[s_dl1b]     = TH1D(s_dl1b,s_dl1b,5000,0,15)
                    h[s_dl1b].SetName(s_dl1b)
                    h[s_dl1b].SetXTitle('distance')
                    h[s_dl1b].SetYTitle('event rate [%s]^{-1}'%(timeScale))

                    s_dl2b       ="%s_pep_ped_%s_%s_%s_%s"%('eisi',cover,ii,locj,'proxCutRate')
                    h[s_dl2b]     = TH2D(s_dl2b,s_dl2b,200,0,200,200,0,200)
                    h[s_dl2b].SetXTitle('prompt energy [pe]')
                    h[s_dl2b].SetYTitle('delayed energy [pe]')
                    h[s_dl2b].SetZTitle('event rate [%s]^{-1}'%(timeScale))

                    cntEISI    = 0.
                    cntAcc  = 0.
                    for index in range(N-1):
                        _rad = sqrt(power(x[index]-x[index+1],2)+power(y[index]-y[index+1],2)+power(z[index]-z[index+1],2))/1000.
                        _pe0 = pe[index]
                        _pe1 = pe[index+1]
                        if _rad < float(arguments["-d"]):
                            h[s_dl2].Fill(_pe0,_pe1)
                            h[s_dl1].Fill(_rad)
                            cntEISI+=1
                        h[s_dl2b].Fill(_pe0,_pe1)
                        h[s_dl1b].Fill(_rad)
                        cntAcc+=1.
                    if cntAcc !=0:
                        print 'Cut effieciency ',cntEISI,cntAcc,cntEISI/cntAcc
                    g["eisi_%s_%s_singlesRate_ProxCut" %(ii,locj)].SetPoint(cntG,pc_val["%s"%(cover)],cntEISI/float(tt)*er)
                    g["eisi_%s_%s_singlesRate_ProxCut" %(ii,locj)].GetXaxis().SetTitle('PMT coverage')
                    g["eisi_%s_%s_singlesRate_ProxCut" %(ii,locj)].GetYaxis().SetTitle('counts [%s]^{-1}'%(timeScale))

                    g["eisi_%s_%s_abs_ProxCut" %(ii,locj)].SetPoint(cntG,pc_val["%s"%(cover)],cntEISI/float(tt))
                    g["eisi_%s_%s_abs_ProxCut" %(ii,locj)].GetXaxis().SetTitle('PMT coverage')
                    g["eisi_%s_%s_abs_ProxCut" %(ii,locj)].GetYaxis().SetTitle('efficiency')

                    if locj == 'PMT' or locj == 'FV':
                        print "(------- %7d  ------  ------ ) %8d  | event rate (--------- %4.3e --------- ---------) per %s (after time redux and prox)" % (si, tt,timeRedux*power(cntEISI/float(tt)*er,2),timeScale)
                        print "Total in-fidual raw rate ","%s"%(cover), rawTotalRateEISI["%s"%(cover)]," %s^{-1}"%(timeScale)
                        er*=pc_num["%s"%(cover)]*mass
                        print "Total total raw rate ","%s"%(cover), rawTotalRateTot["%s"%(cover)]," %s^{-1}"%(timeScale)

                    cntG+=1

                    h[s_dl2].Scale(1./tt)
                    h[s_dl1].Scale(1./tt)
                    h[s_dl2b].Scale(1./tt*er)
                    h[s_dl1b].Scale(1./tt*er)

                    f_root.cd()
                    h[s_dl].Write()
                    h[s_dlsi].Write()
                    h[s_dlT].Write()
                    h[s_dlTsi].Write()
                    h[s_dl1].Write()
                    h[s_dl1b].Write()
                    h[s_dl2].Write()
                    h[s_dl2b].Write()
                    if cover == '30pct':
                        g["si_%s_%s_1_abs" %(ii,locj)].Write()
                        g["ei_%s_%s_1_abs" %(ii,locj)].Write()
                        g["so_%s_%s_1_abs" %(ii,locj)].Write()
                        g["eo_%s_%s_1_abs" %(ii,locj)].Write()
                        g["si_%s_%s_singlesRate" %(ii,locj)].Write()
                        g["ei_%s_%s_singlesRate" %(ii,locj)].Write()
                        g["so_%s_%s_singlesRate" %(ii,locj)].Write()
                        g["eo_%s_%s_singlesRate" %(ii,locj)].Write()
                        g["eisi_%s_%s_singlesRate_ProxCut" %(ii,locj)].Write()
                        g["eisi_%s_%s_abs_ProxCut" %(ii,locj)].Write()

                except:
                    print "Could not read file ",s
                    #        writeResultsToFile(arguments["-o"],g,h)
                    #        print h

    print "\n\n\nThe following file has been created for your convenience: ",_str,"\n\n"
    print "Total in-fiducial raw rate ",rawTotalRateEISI
    print "Total raw rate ",rawTotalRateTot
    #        f.Close()


    return h


def logx_logy_array(nbins = 500,xmin = 1e-2,xmax = 30.,ymin = 1e-9,ymax = 1e3):
    #x-axis
    logxmin = log10(xmin)
    logxmax = log10(xmax)
    xbinwidth= (logxmax-logxmin)/float(nbins)
    xbins	= zeros(nbins,dtype=float)
    xbins[0] = xmin
    #y-axis
    logymin = log10(ymin)
    logymax = log10(ymax)
    ybinwidth= (logymax-logymin)/float(nbins)
    ybins	= zeros(nbins,dtype=float)
    ybins[0] = ymin
    for i in range(nbins):
        xbins[i] = xmin + pow(10,logxmin+i*xbinwidth)
        ybins[i] = ymin + pow(10,logymin+i*ybinwidth)
    return nbins,xbins,ybins

def obtainAbsoluteEfficiency(f,timeScale='day',cut = 10.0):

#   covPCT      = {'9.86037':1432., '14.887':2162.,'19.4453':2824.,'24.994':3558.,'28.8925':4196.,'34.3254':4985.,'39.1385':5684.}
#   pct         = npa([9.86037,14.887,19.4453,24.994,28.8925,34.3254,39.1385])
    pctVal      = ['20pct','25pct','30pct']

    f           = TFile(f,'read')
    EG          = {}
    # print 'Before'
    # print covPCT
    # print pct

    inta,proc,loca,acc,arr,Activity,br,site,timeS,boulbyNeutronFV,mass,dAct,cove,covePCT,covPCT,pct = loadAnalysisParameters(timeScale)
    #print boulbyNeutronFV
    #print 'activity',Activity
    #print dAct
    # print 'After'
    # print covPCT
    # print pct

    x,y         = Double(0.),Double(0.)
    boolFirst   = True
    for _inta in inta:
        for ii in range(len(proc)):
            _str1               = "%s_%s_%s_1_abs" %(_inta,proc[ii],loca[ii])
            _scal_str1          = "scaled_%s_%s_%s%s_1_abs_%s"%(_inta,proc[ii],\
            loca[ii],site[ii],acc[ii])
            _strEff             = "eff_%s_%s_%s%s_1_abs" %(_inta,proc[ii],\
            loca[ii],site[ii])
            EG[_strEff]         = Graph()
            cnter = 0
            EG[_strEff].SetPoint(cnter,0.,0.)
            cnter+=1
            for cc,val in enumerate(pctVal):
                s              = "%s_pe_%s_%s_%s_1"%(_inta,val,proc[ii],loca[ii])
                eff = histIntegral(s,f,cut)
                if eff>0.00:
                    EG[_strEff].SetPoint(cnter,pct[cc],eff)
                    cnter+=1
                    #                    if loca[ii]=="RN":
                    # print s,_inta,val,proc[ii],loca[ii],cnter,pct[cc],eff
                    # print s,cnter,eff
	    if 1==1:
 #           try:
 #		print _str1
                EG[_str1] = f.Get(_str1)
                EG[_scal_str1] = EG[_str1].Clone()

                for i in range(EG[_str1].GetN()):
                    EG[_str1].GetPoint(i,x,y)
                    nY      = y*dAct["%s_%s%s"%(proc[ii],loca[ii],site[ii])]*EG[_strEff].Eval(x)
                    _act = "%s_%s"%(proc[ii],loca[ii])
                    # print 'D',_str1,_act,x,y,dAct[_act],nY
                    # nY1     = y*dAct["%s_%s"%(ii,loca[ii])]
                    #                    if _inta == 'ei' and loca[ii] == 'S':
                    # print 'A',_str1,_strEff,y,Activity[ii],timeS,br[ii],EG[_strEff].Eval(x),nY
                    if loca[ii] == 'PMT':
                        nY      *= mass*covPCT["%s"%(x)]
                        EG[_scal_str1].SetPoint(i,x,nY)
                        # print 'B',nY,i,x,y,covPCT["%s"%(x)]
                    else:
                        EG[_scal_str1].SetPoint(i,x,nY)
                        # print 'C',nY,x,y
  #          except:
#                a = 0
#		print 'Did not succeed'

    _scal_acc,_scal_acc_notFV,_scal_acc1,_scal_acc_notFV1 ="scaled_accidental",\
    "scaled_accidental_notFV","cut_accidental","all_cut_accidental"
    EG[_scal_acc],EG[_scal_acc_notFV],EG[_scal_acc1],EG[_scal_acc_notFV1] = \
    Graph(),Graph(),Graph(),Graph()

    for i,p  in enumerate(pct):
        EG[_scal_acc].SetPoint(i ,p ,0.0)
        EG[_scal_acc_notFV].SetPoint(i,p,0.)
        EG[_scal_acc1].SetPoint(i,p,0.)
        EG[_scal_acc_notFV1].SetPoint(i,p,0.)


    for _inta in inta:
        for ii in range(len(proc)):
            _scal_str1          = "scaled_%s_%s_%s%s_1_abs_%s"%(_inta,\
            proc[ii],loca[ii],site[ii],acc[ii])
            for iii,value in enumerate(pct):
                try:
                    if acc[ii] == 'acc' and (_inta == 'si' or _inta == 'ei'):
                        x   = float(value)
                        oY  = EG[_scal_acc].Eval(x)
                        nY  = EG[_scal_str1].Eval(float(value))
                        EG[_scal_acc].SetPoint(iii,x,oY+nY)
                    if acc[ii] == 'di' and (_inta == 'si' or _inta == 'ei'):
                        x   = float(value)
                        aY  = EG['scaled_ei_neutron_Nboulby_1_abs_corr'].Eval(x)
                        oY  = aY/(boulbyNeutronFV*timeS)
                        # print 'neutron',x,aY, oY, boulbyNeutronFV*timeS
                        nY  = EG[_scal_str1].Eval(float(value))
                        EG[_scal_str1].SetPoint(iii,x,oY*nY)
                    if acc[ii] == 'acc' and _inta != 'si' and _inta != 'ei':
                        x   = float(value)
                        oY  = EG[_scal_acc_notFV].Eval(x)
                        nY  = EG[_scal_str1].Eval(float(value))
                        EG[_scal_acc_notFV].SetPoint(iii,x,oY+nY)
                except:
                    a = 1

    #DistanceEff was taken from watchmanDetector.pdf (0.004)
    # New values are estimated from watch -a results
    distanceEff = 0.025
    dT = float(arguments["-t"])*1e-6 # Change microsecond to rate
    for iii,value in enumerate(pct):
        x   = float(value)
        oY  = EG[_scal_acc].Eval(x)
        nY  = 0.0001/timeS*oY*oY
        EG[_scal_acc1].SetPoint(iii,x,nY)
        EG[_scal_acc_notFV1].SetPoint(iii,x,nY*distanceEff)
    return EG

def pickColor(H,_loc,r_c,o_c,b_c,c_c ):
    if _loc=='PMT':
        H.SetLineColor(kOrange+o_c)
        H.SetFillColor(kOrange+o_c)
        H.SetMarkerColor(kOrange+o_c)
        o_c+=1
    if _loc=='RN':
        H.SetLineColor(7)
        H.SetMarkerColor(7)
    if _loc=='FV':
        H.SetLineColor(kRed+r_c)
        H.SetFillColor(kRed+r_c)
        H.SetMarkerColor(kRed+r_c)
        r_c+=1
    if _loc=='S':
        H.SetLineColor(kBlue+b_c)
        H.SetFillColor(kBlue+b_c)
        H.SetMarkerColor(kBlue+b_c)
        b_c+=1
    if _loc=='N':
        H.SetLineColor(kCyan+c_c)
        H.SetMarkerColor(kCyan+c_c)
        c_c+=1
    return r_c,o_c,b_c,c_c, H

def integralCoincidence(R,lowerBound,upperBound):
    low = -exp(-lowerBound) * (1+lowerBound )
    up  = -exp(-upperBound) * (1+upperBound)
    return up - low

def histIntegral(s,f,cut):
    # H = {}
    Hs = f.Get(s)

    if Hs!= None:
        a = Hs.Integral(0,int(cut*10))
        N = Hs.Integral(0,5000)
    else:
        a = 1.0
        N = 1.0

    # print 'histInt',cut,s,Hs,a,N
    if N !=0:
        return (1.0 - a/N )
    else:
        return 0

def runAnalysisProcess(f,g,h):
    d,iso,loc,coverage,coveragePCT = loadSimulationParameters()
    aDT1 =tDT1 = a1 = t1 = 0

    for j in range(len(iso)):
        for ii in d["%s"%(iso[int(j)])]:
            for idx,cover in enumerate(coverage):
                # print  coverage[idx],coveragePCT[cover],ii,loc[j]
                h  = fillHistograms(f,a1,t1,h,cover,ii,loc[j],\
                float(coveragePCT[cover]))
    return g,h
