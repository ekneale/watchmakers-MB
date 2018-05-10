from load import *
from io_operations import testEnabledCondition,writeResultsToFile
from rootpy.interactive import wait
from ROOT import TH1D, TH1F,TH2F,TProfile2D


# get histograms of PMT accidentals and IBD events as a function of fiducial cut and n9 cut
def extractPMTHistograms():   #--pmtanalysis

    d,iso,loc,coverage,coveragePCT = loadSimulationParameters()
    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)
    _str = "ntuple_root_files%s/SBR.root" %(additionalString)
    f_root = TFile(_str,"recreate")
    #get the parameters from the watch arguments
    pmtDist  		= float(arguments["--tankRadius"])-float(arguments['--steelThick'])-float(arguments["--shieldThick"])
    pmtDistZ 		= float(arguments["--halfHeight"])-float(arguments['--steelThick'])-float(arguments["--shieldThick"])

    timeScale           = arguments["--timeScale"]
    inFilePrefix        = arguments["--ft"] 
    timeCut             = float(arguments["-t"])*1e3
    distCut             = float(arguments["-d"])
    
    #import the analysis parameters from load.py
    parameters  	= loadAnalysisParameters(timeScale)
    boulbyRate		= parameters[9]
    rates       	= parameters[11]
    mass        	= parameters[10]
    pc_num      	= parameters[12]
    timeS		= parameters[8]  
 
    cover = "25pct"

    # test
#    hist1 = TProfile2D("hist1","",10,6.4,7.4,32,8.,40.) 
#    hist2 = TProfile2D("hist2","",10,6.4,7.4,32,8.,40.)
#    hist3 = TProfile2D("hist3","",10,6.4,7.4,32,8.,40.)
    hist1 = TProfile2D("hist1","",20,6.,8.0,32,8.,40.) 
    hist2 = TProfile2D("hist2","",20,6.,8.0,32,8.,40.)
    hist3 = TProfile2D("hist3","",20,6.,8.0,32,8.,40.)

    #test
#    fid = np.arange(5.0,8.0,0.5)
#    n9  = np.arange(8.,24.,2.)
    fid = np.arange(6.0,8.1,0.1)
    n9 = np.arange(8.,40.,1.)
    no = len(fid)
    eiPMT={}
    for jj,fidCut in enumerate(fid):
	for i,n9Cut in enumerate(n9):
	    eiPMT["%s_%s"%(fidCut,n9Cut)] = 0 
    eiIBD={}
    s2b={}
    
    #iterate over processes
    for j in range(len(iso)):
	for ii in d["%s"%(iso[int(j)])]:
            locj = loc[j]		
	    isotope = ii
	    s = "ntuple_root_files%s/%s_%s_%s_%s.root"%(additionalString,inFilePrefix,isotope,cover,locj)
	    rfile = TFile(s)
	    t   = rfile.Get('data')
	    er		= float(rates["%s_%s"%(ii,locj)]) # activity (per day)
            
	    # get the data for the PMT backgrounds
	    if locj == 'PMT':

		print "Reading in ",s
		print "Rate is ",er

	        #iterate over fiducial cut
  		for jj,fidCut in enumerate(fid):
		    #iterate over n9 cut
		    for i,n9Cut in enumerate(n9):
 		    # set out the conditions that define a detected event
		        posGood		= "(pos_goodness>%f)" %(float(arguments["-g"]))
		        n9Good		= "(n9>%f)" %(float(arguments["--minN9"]))
		        totEvtStr	= "all_ev_tot == all_ev"
		        recoFIDstring   = "(sqrt(pow(posReco.X(),2) + pow(posReco.Y(),2))/1000. <%f && sqrt(pow(posReco.Z(),2))/1000. <%f)"%(fidCut,fidCut)
		        tt		= t.Draw("pe",totEvtStr,"goff") # total number of events simulated

			# find the pe efficiency for each cut on n9 (i.e. perform the cut on n9)
			drawArg		= "n9>>hei(5000,0,500)"
			ei 		= t.Draw(drawArg,"%s && %s && %s" %(posGood,n9Good,recoFIDstring),"goff") 
			Hs		= t.GetHistogram()
			if Hs != None:
			    pe_cut 	= Hs.Integral(0,int(n9Cut*10))
			    pe_tot	= Hs.Integral(0,5000)
			else:
			    pe_cut	= 1.0
			    pe_tot	= 1.0
			if pe_tot !=0:
			    pe_eff	= 1.0 - pe_cut/pe_tot
			else:
			    pe_eff 	= 0
			# get number of single events in the FV
                        N  		= t.Draw("pe:posReco.X():posReco.Y():posReco.Z()","%s && %s && %s" %(recoFIDstring,posGood,n9Good),"goff")
                        x    		= t.GetV2()
                        y    		= t.GetV3()
                        z    		= t.GetV4()
			# find the number of single events which are within 2m of each other 
		        cnt = 0.
                        for index in range(N-1):
                            dist = sqrt(power(x[index]-x[index+1],2)+power(y[index]-y[index+1],2)+power(z[index]-z[index+1],2))/1000.
                            if dist < float(arguments["-d"]):
                                cnt+=1
			#scale the number of distance-correlated events by the reconstruction efficiency, total radioactivity rate (per day) and cut on n9
		        eiPMT_tmp = cnt/tt*er*pc_num[cover]*mass*pe_eff
			# find the time-correlated rate per day
			eiPMT_tmp *= eiPMT_tmp*0.0001/(24.*3600.)
			#fill the histogram
  		        hist1.Fill(fidCut,n9Cut,eiPMT_tmp)
		        eiPMT["%s_%s"%(fidCut,n9Cut)]	+= eiPMT_tmp
		        print fidCut,n9Cut,eiPMT_tmp,eiPMT["%s_%s"%(fidCut,n9Cut)]

# 	   backgroundNoise = 1.0e3*float(pc_num["%s"%(cover)])*1500.*1e-9 
	    
            # get the data for the IBD signal
	    elif locj == 'S':
	        if ii == 'boulby':
	    	    print 'reading in ',s	
		    # iterate over fiducial cut
   	    	    for jj,fidCut in enumerate(fid):
			for i,n9Cut in enumerate(n9):
 	        	# set out the conditions that define a detected event
	        	    posGood        = "(pos_goodness>%f)" %(float(arguments["-g"]))
	        	    n9Good         = "(n9>%f)" %(float(arguments["--minN9"]))
	        	    totEvtStr      = "all_ev_tot == all_ev"
			    tt             = t.Draw("pe",totEvtStr,"goff") # total number of events simulated
	        	    boulbyIBDRate  = rates["boulby_S"] # boulby rate per day
	    		    recoFIDstring  = "(sqrt(pow(posReco.X(),2) + pow(posReco.Y(),2))/1000. <%f && sqrt(pow(posReco.Z(),2))/1000. <%f)"%(fidCut,fidCut) # event is reconstructed within the inner volume

			    # find the pe efficiency for each cut on n9 (perform the n9 cut)
			    drawArg		= "n9>>hei(5000,0,500)"
		 	    ei 		= t.Draw(drawArg,"%s && %s && %s" %(posGood,n9Good,recoFIDstring),"goff")
			    Hs		= t.GetHistogram()
			    if Hs != None:
				pe_cut 	= Hs.Integral(0,int(n9Cut*10))
			    	pe_tot	= Hs.Integral(0,5000)
			    else:
			        pe_cut	= 1.0
			        pe_tot	= 1.0
			    if pe_tot !=0:
				pe_eff	= 1.0 - pe_cut/pe_tot
			    else: 
				pe_eff	= 0

			    eiIBD_tmp = t.Draw("reco_r","%s && %s && %s"%(recoFIDstring,posGood,n9Good),"goff") #get no of detected events
			    eiIBD_tmp *= 1./tt*boulbyIBDRate*pe_eff #scale by rate, reconstruction efficiency and pe efficiency
			    hist2.Fill(fidCut,n9Cut,eiIBD_tmp)
			    eiIBD["%s_%s"%(fidCut,n9Cut)] = eiIBD_tmp
			    print fidCut,n9Cut,eiIBD_tmp

    for jj,fidCut in enumerate(fid):
	for i,n9Cut in enumerate(n9):
	    s2b["%s_%s"%(fidCut,n9Cut)] = eiIBD["%s_%s"%(fidCut,n9Cut)]/sqrt(eiPMT["%s_%s"%(fidCut,n9Cut)]+eiIBD["%s_%s"%(fidCut,n9Cut)])
	    hist3.Fill(fidCut,n9Cut,s2b["%s_%s"%(fidCut,n9Cut)])
	    
    f_root.cd()
    c1 = TCanvas("c1","")
    c1.SetRightMargin(0.15)
    hist1.Draw("COLZ2")
    hist1.GetXaxis().SetTitle("Fiducial cut (m)")
    hist1.GetXaxis().SetTitleSize(0.045)
    hist1.GetYaxis().SetTitle("n9 cut")
    hist1.GetYaxis().SetTitleSize(0.045)
    hist1.GetZaxis().SetTitle("Events per day")
    hist1.GetZaxis().SetTitleSize(0.045)
    hist1.SetTitle("PMT accidentals)")

    c2 = TCanvas("c2","")
    c2.SetRightMargin(0.15)
    hist2.Draw("COLZ2")
    hist2.GetXaxis().SetTitle("Fiducial cut (m)")
    hist2.GetXaxis().SetTitleSize(0.045)
    hist2.GetYaxis().SetTitle("n9 cut")
    hist2.GetYaxis().SetTitleSize(0.045)
    hist2.GetZaxis().SetTitle("Events per day")
    hist2.GetZaxis().SetTitleSize(0.045)
    hist2.SetTitle("IBD events")
    
    c3 = TCanvas("c3","")
    c3.SetRightMargin(0.15)
    hist3.Draw("COLZ2")
    hist3.GetXaxis().SetTitle("Fiducial cut (m)")
    hist3.GetXaxis().SetTitleSize(0.045)
    hist3.GetYaxis().SetTitle("n9Cut")
    hist3.GetYaxis().SetTitleSize(0.045)
    hist3.GetZaxis().SetTitle("Signal to background ratio")
    hist3.GetZaxis().SetTitleSize(0.045)
    hist3.SetTitle("Signal to PMT background")

    wait(True)
    wait

    return hist1.Write(),hist2.Write() ,hist3.Write() 

# get histograms of pmt singles as a function of radius r showing contributions from the different isotopes
def extractPMTSinglesHistograms(): # --pmtanalysis2   

    d,iso,loc,coverage,coveragePCT = loadSimulationParameters()
    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)
    _str = "ntuple_root_files%s/SBR.root" %(additionalString)
    f_root = TFile(_str,"recreate")
    #get the parameters from the watch arguments
    pmtDist  		= float(arguments["--tankRadius"])-float(arguments['--steelThick'])-float(arguments["--shieldThick"])
    pmtDistZ 		= float(arguments["--halfHeight"])-float(arguments['--steelThick'])-float(arguments["--shieldThick"])

    timeScale           = arguments["--timeScale"]
    inFilePrefix        = arguments["--ft"] 
    timeCut             = float(arguments["-t"])*1e3
    distCut             = float(arguments["-d"])
    
    #import the analysis parameters from load.py
    parameters  	= loadAnalysisParameters(timeScale)
    boulbyRate		= parameters[9]
    rates       	= parameters[11]
    mass        	= parameters[10]
    pc_num      	= parameters[12]
    timeS		= parameters[8]  
 
    cover = "25pct"

    h1,h1new = {},{}
    h2,h2new = {},{}
    ths1 = THStack("ths1","PMT singles")
    ths2 = THStack("ths2","PMT singles (close-up)")
    h4 = TH1D("h4","",850,0.,8.5)

    n={'234Pa':2,'214Pb':46,'214Bi':4,'210Bi':6,'210Tl':5,'228Ac':3,'212Pb':51,'212Bi':7,'208Tl':1}

    #iterate over processes
    for j in range(len(iso)):
	for ii in d["%s"%(iso[int(j)])]:
            locj = loc[j]		
	    isotope = ii
	    s = "ntuple_root_files%s/%s_%s_%s_%s.root"%(additionalString,inFilePrefix,isotope,cover,locj)
	    rfile = TFile(s)
	    t   = rfile.Get('data')
            
	    # set out the conditions that define an event as a relevant background
   	    recoPMTstring	= "(sqrt(pow(posReco.X(),2) + pow(posReco.Y(),2))/1000. <%f && sqrt(pow(posReco.Z(),2))/1000. <%f)"%(pmtDist,pmtDistZ) # event is reconstructed within the inner volume
	    posGood		= "(pos_goodness>%f)" %(float(arguments["-g"]))
	    n9Good		= "(n9>%f)" %(float(arguments["--minN9"]))
	    totEvtStr	    	= "all_ev_tot == all_ev"

	    if locj == "PMT":

		print "Reading in ",s

		tt = t.Draw("pe",totEvtStr,"goff") # number of events detected / total number of events simulated
		er = float(rates["%s_%s"%(ii,locj)]) # activity (Bq kg^-1)
			
		# perform the cut on n9 
		drawArg		= "n9>>hei(5000,0,500)"
		ei 		= t.Draw(drawArg,"%s && %s && %s" %(posGood,n9Good,recoPMTstring),"goff") 
		Hs		= t.GetHistogram()
		if Hs != None:
		    pe_cut 	= Hs.Integral(0,90) # minN9==9 TODO work out why we multiply by 10 here
		    pe_tot	= Hs.Integral(0,5000)
		else:
		    pe_cut	= 1.0
		    pe_tot	= 1.0
		if pe_tot != 0:
		    pe_eff	= 1.0 - pe_cut/pe_tot
		else:
		    pe_eff	= 0
				
	
		# now let's create histograms of single events in the PMTs
		string = "h1_%s"%(ii)
		h1[string] = TH1D("h1_%s"%(ii),"",400,4.,8.)
		h1[string].SetMarkerStyle(1)
		h1[string].SetLineColor(n["%s"%(ii)])
		# get histogram of events at position r
		t.Draw("reco_r>>h1_%s"%(ii),"%s && %s && %s"%(recoPMTstring,posGood,n9Good),"goff") 
		# scale number of events by total activity and efficiency and n9 cut
		h1[string].Scale(1./tt*er*pc_num[cover]*mass*pe_eff)
		# add the  histogram to the stack
		gROOT.cd()
		stringNew = "h1new_%s"%(ii)
		h1new[stringNew] = h1[string].Clone()
		ths1.Add(h1new[stringNew])
		

		string = "h2_%s"%(ii)
		h2[string] = TH1D("h2_%s"%(ii),"",300,4.5,7.5)
		h2[string].SetMarkerStyle(1)
		h2[string].SetLineColor(n["%s"%(ii)])
		h2[string].SetLineWidth(1)
		t.Draw("reco_r>>h2_%s"%(ii),"%s && %s && %s"%(recoPMTstring,posGood,n9Good),"goff") # get histogram of events at position r
		h2[string].Scale(1./tt*er*pc_num['25pct']*mass) # scale number of events by total activity and efficiency
		# add the  histogram to the stack
		gROOT.cd()
		stringNew = "h2new_%s"%(ii)
		h2new[stringNew] = h2[string].Clone()
		ths2.Add(h2new[stringNew])

#		string = "h3_%s"%(ii)
#		h3[string] = TH1D("h3_%s"%(ii),"",850,0.,8.5)
#		h3[string].SetMarkerStyle(1)
##		h3[string].SetFillColor(n["%s"%(ii)])
#		h3[string].SetLineColor(n["%s"%(ii)])
#		t.Draw("reco_r>>h3_%s"%(ii),"%s && %s && %s"%(recoPMTstring,posGood,n9Good),"goff") # get histogram of events at position r
#		h3[string].Scale(1./tt*er*pc_num['25pct']*mass) # scale number of events by total activity and efficiency
#		# add the  histogram to the stack
#		gROOT.cd()
#		stringNew = "h3new_%s"%(ii)
#		h3new[stringNew] = h3[string].Clone()
#		ths3.Add(h3new[stringNew])
	   	

    # Draw the plots of PMT events
    f_root.cd()

    # draw the first stack - histogram over full tank 
    ths1.ls()
    c1 = TCanvas('c1','Single events due to PMT glass radioactivity')
    c1.cd()
    ths1.SetTitle("Single events due to PMT glass radioactivity")
    ths1.Draw("hist")
    ths1.GetXaxis().SetTitle("r(m)")
    ths1.GetXaxis().SetTitleSize(0.045)
    ths1.GetYaxis().SetTitle("Events per day")
    ths1.GetYaxis().SetTitleSize(0.045)
    c1.BuildLegend(0.1,0.6,0.3,0.9,"Decaying isotopes")
    

    # draw the second stack - close-up
    ths2.ls()
    c2 = TCanvas('c2','Single events due to PMT glass radioactivity')
    c2.cd()
    ths2.Draw("hist")
    ths2.SetTitle("Single events due to PMT glass radioactivity")
    ths2.GetXaxis().SetTitle("r(m)")
    ths2.GetXaxis().SetTitleSize(0.045)
    ths2.GetYaxis().SetTitle("Events per day")
    ths2.GetYaxis().SetTitleSize(0.045)
    c2.BuildLegend(0.1,0.6,0.3,0.9,"Decaying isotopes")

    #draw the third stack - log plot
#    ths3.ls()
#    c3 = TCanvas('c3','PMT singles')
#    c3.cd()
#    c3.SetLogy()
#    ths3.SetTitle("PMT singles - R7801 low activity glass")
#    ths3.GetStack().Draw("hist")

    wait() 

    return ths1.Write(),ths2.Write() #,h4.Write() 

# get histograms of correlated ibd events, positrons and neutrons as a function of n9 cut and fiducial cut
def extractIBDHistograms():   #--pmtanalysis3

    d,iso,loc,coverage,coveragePCT = loadSimulationParameters()
    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)
    _str = "ntuple_root_files%s/SBR.root" %(additionalString)
    f_root = TFile(_str,"recreate")
    #get the parameters from the watch arguments
    pmtDist  		= float(arguments["--tankRadius"])-float(arguments['--steelThick'])-float(arguments["--shieldThick"])
    pmtDistZ 		= float(arguments["--halfHeight"])-float(arguments['--steelThick'])-float(arguments["--shieldThick"])

    timeScale           = arguments["--timeScale"]
    inFilePrefix        = arguments["--ft"] 
    timeCut             = float(arguments["-t"])*1e3
    distCut             = float(arguments["-d"])
    
    #import the analysis parameters from load.py
    parameters  	= loadAnalysisParameters(timeScale)
    boulbyRate		= parameters[9]
    rates       	= parameters[11]
    mass        	= parameters[10]
    pc_num      	= parameters[12]
    timeS		= parameters[8]  
 
    cover = "25pct"

    # test
#    hist1 = TProfile2D("hist1","",10,6.4,7.4,32,8.,40.) 
#    hist2 = TProfile2D("hist2","",10,6.4,7.4,32,8.,40.)
#    hist3 = TProfile2D("hist3","",10,6.4,7.4,32,8.,40.)
    hist1 = TProfile2D("hist1","",24,6.,8.4,32,8.,40.) 
    hist2 = TProfile2D("hist2","",24,6.,8.4,32,8.,40.)
    hist3 = TProfile2D("hist3","",24,6.,8.4,32,8.,40.)

    #test
#    fid = np.arange(5.0,8.0,0.5)
#    n9  = np.arange(8.,24.,2.)
    fid = np.arange(6.0,8.5,0.1)
    n9 = np.arange(8.,40.,1.)
    no = len(fid)
    eiPMT={}
    for jj,fidCut in enumerate(fid):
	for i,n9Cut in enumerate(n9):
	    eiPMT["%s_%s"%(fidCut,n9Cut)] = 0 
    eiIBD={}
    s2b={}
    
    #iterate over processes
    for j in range(len(iso)):
	for ii in d["%s"%(iso[int(j)])]:
            locj = loc[j]		
	    isotope = ii
	    s = "ntuple_root_files%s/%s_%s_%s_%s.root"%(additionalString,inFilePrefix,isotope,cover,locj)
	    rfile = TFile(s)
	    t   = rfile.Get('data')
	    er		= float(rates["%s_%s"%(ii,locj)]) # activity (per day)
            
            # get the data for the IBD signal
	    if locj == 'S':
	        if ii == 'boulby':
	    	    print 'reading in ',s
		    print 'fidCut','n9Cut','IBDs','     eff','     pe_eff'	
		    # iterate over fiducial cut
   	    	    for jj,fidCut in enumerate(fid):
			for i,n9Cut in enumerate(n9):
 	        	# set out the conditions that define a detected event
	        	    posGood        = "(pos_goodness>%f)" %(float(arguments["-g"]))
	        	    n9Good         = "(n9>%f)" %(float(arguments["--minN9"]))
	        	    totEvtStr      = "all_ev_tot == all_ev"
			    tt             = t.Draw("pe",totEvtStr,"goff") # total number of events simulated
	        	    boulbyIBDRate  = rates["boulby_S"] # boulby rate per day
	    		    recoFIDstring  = "(sqrt(pow(posReco.X(),2) + pow(posReco.Y(),2))/1000. <%f && sqrt(pow(posReco.Z(),2))/1000. <%f)"%(fidCut,fidCut) # event is reconstructed within the inner volume
			    # find the pe efficiency for each cut on n9 (perform the n9 cut)
			    drawArg		= "n9>>hei(5000,0,500)"
		 	    ei 		= t.Draw(drawArg,"%s && %s && %s" %(posGood,n9Good,recoFIDstring),"goff")
			    eff		= ei/float(tt)
			    Hs		= t.GetHistogram()
			    if Hs != None:
				pe_cut 	= Hs.Integral(0,int(n9Cut*10))
			    	pe_tot	= Hs.Integral(0,5000)
			    else:
			        pe_cut	= 1.0
			        pe_tot	= 1.0
			    if pe_tot !=0:
				pe_eff	= 1.0 - pe_cut/pe_tot
			    else: 
				pe_eff	= 0

			    eiIBD_tmp = t.Draw("reco_r","%s && %s && %s"%(recoFIDstring,posGood,n9Good),"goff") #get no of detected events
			    eiIBD_tmp *= 1./float(tt)*boulbyIBDRate*pe_eff #scale by rate, reconstruction efficiency and pe efficiency
			    hist1.Fill(fidCut,n9Cut,eiIBD_tmp)
			    print fidCut,'    ',n9Cut,eiIBD_tmp,eff,pe_eff

            # get the data for the IBD signal
	    elif locj == 'I':
	    	print 'reading in ',s	
		print 'fidCut','n9Cut','Pos','     eff','     pe_eff'	
		# iterate over fiducial cut
   	    	for jj,fidCut in enumerate(fid):
		    for i,n9Cut in enumerate(n9):
 	            # set out the conditions that define a detected event
	        	posGood        = "(pos_goodness>%f)" %(float(arguments["-g"]))
	        	n9Good         = "(n9>%f)" %(float(arguments["--minN9"]))
	        	totEvtStr      = "all_ev_tot == all_ev"
			tt             = t.Draw("pe",totEvtStr,"goff") # total number of events simulated
	        	boulbyIBDRate  = rates["boulby_S"] # boulby rate per day
	    		recoFIDstring  = "(sqrt(pow(posReco.X(),2) + pow(posReco.Y(),2))/1000. <%f && sqrt(pow(posReco.Z(),2))/1000. <%f)"%(fidCut,fidCut) # event is reconstructed within the inner volume

			# find the pe efficiency for each cut on n9 (perform the n9 cut)
			drawArg		= "n9>>hei(5000,0,500)"
		 	ei 		= t.Draw(drawArg,"%s && %s && %s" %(posGood,n9Good,recoFIDstring),"goff")
			eff		= ei/float(tt)
			Hs		= t.GetHistogram()
			if Hs != None:
			    pe_cut 	= Hs.Integral(0,int(n9Cut*10))
			    pe_tot	= Hs.Integral(0,5000)
			else:
			    pe_cut	= 1.0
			    pe_tot	= 1.0
			if pe_tot !=0:
			    pe_eff	= 1.0 - pe_cut/pe_tot
			else: 
			    pe_eff	= 0

			eiPos_tmp = t.Draw("reco_r","%s && %s && %s"%(recoFIDstring,posGood,n9Good),"goff") #get no of detected events
			eiPos_tmp *= 1./float(tt)*boulbyIBDRate*pe_eff #scale by rate, reconstruction efficiency and pe efficiency
			hist2.Fill(fidCut,n9Cut,eiPos_tmp)
			print fidCut,'    ',n9Cut,eiPos_tmp,eff,pe_eff

            # get the data for the IBD signal
	    elif locj == 'N':
	    	print 'reading in ',s	
		print 'fidCut','n9Cut','N','     eff','      pe_eff'	
		# iterate over fiducial cut
   	    	for jj,fidCut in enumerate(fid):
		    for i,n9Cut in enumerate(n9):
 	            # set out the conditions that define a detected event
	        	posGood        = "(pos_goodness>%f)" %(float(arguments["-g"]))
	        	n9Good         = "(n9>%f)" %(float(arguments["--minN9"]))
	        	totEvtStr      = "all_ev_tot == all_ev"
			tt             = t.Draw("pe",totEvtStr,"goff") # total number of events simulated
	        	boulbyIBDRate  = rates["boulby_S"] # boulby rate per day
	    		recoFIDstring  = "(sqrt(pow(posReco.X(),2) + pow(posReco.Y(),2))/1000. <%f && sqrt(pow(posReco.Z(),2))/1000. <%f)"%(fidCut,fidCut) # event is reconstructed within the inner volume

			# find the pe efficiency for each cut on n9 (perform the n9 cut)
			drawArg		= "n9>>hei(5000,0,500)"
		 	ei 		= t.Draw(drawArg,"%s && %s && %s" %(posGood,n9Good,recoFIDstring),"goff")
			eff		= ei/float(tt)
			Hs		= t.GetHistogram()
			if Hs != None:
			    pe_cut 	= Hs.Integral(0,int(n9Cut*10))
			    pe_tot	= Hs.Integral(0,5000)
			else:
			    pe_cut	= 1.0
			    pe_tot	= 1.0
			if pe_tot !=0:
			    pe_eff	= 1.0 - pe_cut/pe_tot
			else: 
			    pe_eff	= 0

			eiN_tmp = t.Draw("reco_r","%s && %s && %s"%(recoFIDstring,posGood,n9Good),"goff") #get no of detected events
			eiN_tmp *= 1./float(tt)*boulbyIBDRate*pe_eff #scale by rate, reconstruction efficiency and pe efficiency
			hist3.Fill(fidCut,n9Cut,eiN_tmp)
			print fidCut,'    ',n9Cut,eiN_tmp,eff,pe_eff
	    
    f_root.cd()
    c1 = TCanvas("c1","")
    c1.SetRightMargin(0.15)
    hist1.Draw("COLZ2")
    hist1.GetXaxis().SetTitle("Fiducial cut (distance from centre in m)")
    hist1.GetXaxis().SetTitleSize(0.045)
    hist1.GetYaxis().SetTitle("n9 cut")
    hist1.GetYaxis().SetTitleSize(0.045)
    hist1.GetZaxis().SetTitle("Events per day")
    hist1.GetZaxis().SetTitleSize(0.045)
    hist1.SetTitle("IBD events")

    c2 = TCanvas("c2","")
    c2.SetRightMargin(0.15)
    hist2.Draw("COLZ2")
    hist2.GetXaxis().SetTitle("Fiducial cut (distance from centre in m)")
    hist2.GetXaxis().SetTitleSize(0.045)
    hist2.GetYaxis().SetTitle("n9 cut")
    hist2.GetYaxis().SetTitleSize(0.045)
    hist2.GetZaxis().SetTitle("Events per day")
    hist2.GetZaxis().SetTitleSize(0.045)
    hist2.SetTitle("Positron events")
    
    c3 = TCanvas("c3","")
    c3.SetRightMargin(0.15)
    hist3.Draw("COLZ2")
    hist3.GetXaxis().SetTitle("Fiducial cut (distance from centre in m)")
    hist3.GetXaxis().SetTitleSize(0.045)
    hist3.GetYaxis().SetTitle("n9 cut")
    hist3.GetYaxis().SetTitleSize(0.045)
    hist3.GetZaxis().SetTitle("Events per day")
    hist3.GetZaxis().SetTitleSize(0.045)
    hist3.SetTitle("Neutron events")


    wait(True)
    wait

    return hist1.Write(),hist2.Write() ,hist3.Write() 
 

'''
#get histograms of pmt and ibds as a function of n9 only 
extractPMTn9Histograms()# --pmtanalysis4

    d,iso,loc,coverage,coveragePCT = loadSimulationParameters()
    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)
    _str = "ntuple_root_files%s/SBR.root" %(additionalString)
    f_root = TFile(_str,"recreate")
    #get the parameters from the watch arguments
    pmtDist  		= float(arguments["--tankRadius"])-float(arguments['--steelThick'])-float(arguments["--shieldThick"])
    pmtDistZ 		= float(arguments["--halfHeight"])-float(arguments['--steelThick'])-float(arguments["--shieldThick"])

    timeScale           = arguments["--timeScale"]
    inFilePrefix        = arguments["--ft"] 
    timeCut             = float(arguments["-t"])*1e3
    distCut             = float(arguments["-d"])
    
    #import the analysis parameters from load.py
    parameters  	= loadAnalysisParameters(timeScale)
    boulbyRate		= parameters[9]
    rates       	= parameters[11]
    mass        	= parameters[10]
    pc_num      	= parameters[12]
    timeS		= parameters[8]  
 
    cover = "25pct"

    # test
#    hist1 = TProfile2D("hist1","",7,5.,8.5,8,8.,24.) 
#    hist1b = TProfile2D("hist1b","",7,5.,8.5,8,8.,24.) 
#    hist2 = TProfile2D("hist2","",7,5.,8.5,8,8.,24.)
#    hist3 = TProfile2D("hist3","",7,5.,8.5,8,8.,24.)
    histN91 = TProfile2D("histN91","",30,5.0,8.0,32,8.,40.) 
    histN92 = TProfile2D("histN92","",30,5.0,8.0,32,8.,40.)
    histN93 = TProfile2D("histN93","",30,5.0,8.0,32,8.,40.)

    #test
#    fid = np.arange(5.0,8.5,0.5)
#    n9  = np.arange(8.,24.,2.)
    fid = np.arange(5.0,8.1,0.1)
    n9 = np.arange(8.,40.,1.)
    no = len(fid)
    eiPMT={}
    for i,n9Cut in enumerate(n9):
	eiPMT["%s"%(n9Cut)] = 0 
    eiIBD={}
    s2b={}
    
    #iterate over processes
    for j in range(len(iso)):
	for ii in d["%s"%(iso[int(j)])]:
            locj = loc[j]		
	    isotope = ii
	    s = "ntuple_root_files%s/%s_%s_%s_%s.root"%(additionalString,inFilePrefix,isotope,cover,locj)
	    rfile = TFile(s)
	    t   = rfile.Get('data')
	    er		= float(rates["%s_%s"%(ii,locj)]) # activity (per day)
            
	    # get the data for the PMT backgrounds
	    if locj == 'PMT':

		print "Reading in ",s
		print "Rate is ",er

		#iterate over n9 cut
		for i,n9Cut in enumerate(n9):
 		# set out the conditions that define a detected event
		    r = t.reco_r
		    posGood		= "(pos_goodness>%f)" %(float(arguments["-g"]))
		    n9Good		= "(n9>%f)" %(float(arguments["--minN9"]))
		    totEvtStr	= "all_ev_tot == all_ev"
		    recoFIDstring   = "(sqrt(pow(posReco.X(),2) + pow(posReco.Y(),2))/1000. <%f && sqrt(pow(posReco.Z(),2))/1000. <%f)"%()
		    tt		= t.Draw("pe",totEvtStr,"goff") # total number of events simulated

		    # find the pe efficiency for each cut on n9 (i.e. perform the cut on n9)
		    drawArg	= "n9>>hei(5000,0,500)"
		    ei 		= t.Draw(drawArg,"%s && %s && %s" %(posGood,n9Good,recoPMTstring),"goff") 
		    Hs		= t.GetHistogram()
		    if Hs != None:
		        pe_cut 	= Hs.Integral(0,int(n9Cut*10))
		        pe_tot	= Hs.Integral(0,5000)
		    else:
			pe_cut	= 1.0
			pe_tot	= 1.0
		    if pe_tot  !=0:
			pe_eff	= 1.0 - pe_cut/pe_tot
		    else:
			pe_eff 	= 0
	  	    # get number of single events in the detector
                    N  		= t.Draw("pe:posReco.X():posReco.Y():posReco.Z()","%s && %s && %s" %(recoPMTstring,posGood,n9Good),"goff")
                    x    		= t.GetV2()
                    y    		= t.GetV3()
                    z    		= t.GetV4()
		    # find the number of single events which are within 2m of each other
		    cnt = 0.
                    for index in range(N-1):
                        dist = sqrt(power(x[index]-x[index+1],2)+power(y[index]-y[index+1],2)+power(z[index]-z[index+1],2))/1000.
                        if dist < float(arguments["-d"]):
                            cnt+=1
		    eiPMT_tmp = t.Draw("reco_r","%s && %s && %s"%(posGood,n9Good,recoPMTstring),"goff")
		    eiPMT_tmp *= cnt/eiPMT_tmp*tt*er*pc_num[cover]*mass*pe_eff #scale by reconstruction efficiency, total radioactivity rate, and pe efficiency
 		    # find the correlated event rate per day
		    eiPMT_tmp *= eiPMT_tmp* 0.0001/(24.*3600.)
  		    histN91.Fill(n9Cut,eiPMT_tmp)
		    eiPMT["%s"%(n9Cut)]	+= eiPMT_tmp
		    print n9Cut,eiPMT_tmp,eiPMT["%s"%(n9Cut)]

# 	    backgroundNoise = 1.0e3*float(pc_num["%s"%(cover)])*1500.*1e-9 
	    
            # get the data for the IBD signal
	    if locj == 'S':
	        if ii == 'boulby':
	    	    print 'reading in ',s	
		    for i,n9Cut in enumerate(n9):
 	            # set out the conditions that define a detected event
			posGood        = "(pos_goodness>%f)" %(float(arguments["-g"]))
	                n9Good         = "(n9>%f)" %(float(arguments["--minN9"]))
	                totEvtStr      = "all_ev_tot == all_ev"
		        tt             = t.Draw("pe",totEvtStr,"goff") # total number of events simulated
	                boulbyIBDRate  = rates["boulby_S"] # boulby rate per day
	    	        recoPMTstring  = "(sqrt(pow(posReco.X(),2) + pow(posReco.Y(),2))/1000. <%f && sqrt(pow(posReco.Z(),2))/1000. <%f)"%(pmtDist,pmtDistZ) # event is reconstructed within the detector volume

		        # find the pe efficiency for each cut on n9 (perform the n9 cut)
		        drawArg		= "n9>>hei(5000,0,500)"
		        ei 		= t.Draw(drawArg,"%s && %s && %s" %(posGood,n9Good,recoPMTstring),"goff")
		        Hs		= t.GetHistogram()
		        if Hs != None:
		    	    pe_cut 	= Hs.Integral(0,int(n9Cut*10))
			    pe_tot	= Hs.Integral(0,5000)
			else:
			    pe_cut	= 1.0
			    pe_tot	= 1.0
			if pe_tot !=0:
			    pe_eff	= 1.0 - pe_cut/pe_tot
			else: 
			    pe_eff	= 0

			# no detected events and scale by rate, reconstruction efficiency and pe efficiency (n9 cut)
			eiIBD_tmp = t.Draw("n9Cut:reco_r>>histN92","%s && %s && %s"%(recoPMTstring,posGood,n9Good),"goff")/tt*boulbyIBDRate*pe_eff
#			eiIBD_tmp *= 1./tt*boulbyIBDRate*pe_eff #scale by rate, reconstruction efficiency and pe efficiency (n9 cut)
#			histN92.Fill(r,n9Cut,eiIBD_tmp)
			eiIBD["%s"%(n9Cut)] = eiIBD_tmp
			print n9Cut,eiIBD_tmp

#    for i,n9Cut in enumerate(n9):
#	s2b["%s"%(n9Cut)] = eiIBD["%s"%(n9Cut)]/sqrt(eiPMT["%s"%(n9Cut)]+eiIBD["%s"%(n9Cut)])
#	histN93.Fill(reco_r,n9Cut,s2b["%s"%(n9Cut)])
	    
    f_root.cd()
    c1 = TCanvas("c1","")
    histN91.Draw("COLZ2")
    histN91.GetXaxis().SetTitle("r (m)")
    histN91.GetYaxis().SetTitle("n9 cut")
    histN91.SetTitle("PMT accidentals (proximity-correlated)")

    c2 = TCanvas("c2","")
    histN92.Draw("COLZ2")
    histN92.GetXaxis().SetTitle("r (m)")
    histN92.GetYaxis().SetTitle("n9 cut")
    histN92.SetTitle("IBD events")
    
#    c3 = TCanvas("c3","")
#    histN93.Draw("COLZ2")
#    histN93.GetXaxis().SetTitle("r (m)")
#    histN93.GetYaxis().SetTitle("s/sqrt(s+b)")
#    histN93.SetTitle("Signal to PMT background")


    wait(True)
    wait

    return histN91.Write(),histN92.Write() ,histN93.Write() 


# extract the histograms of singles as a function of fiducial cut
def extractPMTSinglesFidCutHistograms():   # --pmtanalysis5

    d,iso,loc,coverage,coveragePCT = loadSimulationParameters()
    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)
    _str = "ntuple_root_files%s/SBR.root" %(additionalString)
    f_root = TFile(_str,"recreate")
    #get the parameters from the watch arguments
    pmtDist  		= float(arguments["--tankRadius"])-float(arguments['--steelThick'])-float(arguments["--shieldThick"])
    pmtDistZ 		= float(arguments["--halfHeight"])-float(arguments['--steelThick'])-float(arguments["--shieldThick"])

    timeScale           = arguments["--timeScale"]
    inFilePrefix        = arguments["--ft"] 
    timeCut             = float(arguments["-t"])*1e3
    distCut             = float(arguments["-d"])
    
    #import the analysis parameters from load.py
    parameters  	= loadAnalysisParameters(timeScale)
    boulbyRate		= parameters[9]
    rates       	= parameters[11]
    mass        	= parameters[10]
    pc_num      	= parameters[12]
    timeS		= parameters[8]  
 
    cover = "25pct"

    histSinglesFid = TH1F("histSingles2","PMT singles",55,3.,8.5)
    histIBDfid = TH1F("histIBDfid","",55,3.,8.5)

    fid = np.arange(3.,8.5,0.1)
    n = len(fid)
    eiPMTvals=[0]*n
    eiIBDvals=[]

    #iterate over processes
    for j in range(len(iso)):
	for ii in d["%s"%(iso[int(j)])]:
            locj = loc[j]		
	    isotope = ii
	    s = "ntuple_root_files%s/%s_%s_%s_%s.root"%(additionalString,inFilePrefix,isotope,cover,locj)
	    rfile = TFile(s)
	    t   = rfile.Get('data')
            

	    if locj == "PMT":

		print "Reading in ",s
		for jj,fidCut in enumerate(fid):
		    posGood		= "(pos_goodness>%f)" %(float(arguments["-g"]))
		    n9Good		= "(n9>%f)" %(float(arguments["--minN9"]))
		    totEvtStr	    	= "all_ev_tot == all_ev"
   	    	    recoFIDstring	= "(sqrt(pow(posReco.X(),2) + pow(posReco.Y(),2))/1000. <%f && sqrt(pow(posReco.Z(),2))/1000. <%f)"%(fidCut,fidCut) # event is reconstructed within the inner volume
		    tt = t.Draw("pe",totEvtStr,"goff") # number of events detected / total number of events simulated
		    er = float(rates["%s_%s"%(ii,locj)]) # activity (Bq kg^-1)
			
		    # find the pe efficiency
		    drawArg	= "n9>>hei(5000,0,500)"
		    ei 		= t.Draw(drawArg,"%s && %s && %s" %(posGood,n9Good,recoFIDstring),"goff") 
		    Hs		= t.GetHistogram()
		    if Hs != None:
			pe_cut 	= Hs.Integral(0,90) # minN9==9.
		    	pe_tot	= Hs.Integral(0,5000)
		    else:
		        pe_cut	= 1.0
		        pe_tot	= 1.0
		    if pe_tot != 0:
			pe_eff	= 1.0 - pe_cut/pe_tot
		    else:
			pe_eff  = 0

		    eiPMT = t.Draw("reco_r","%s && %s && %s"%(recoFIDstring,posGood,n9Good),"goff")
		    eiPMT*= 1./tt*er*pc_num[cover]*mass*pe_eff #scale by reconstruction efficiency, total radioactivity rate, and pe efficiency
  		    histSinglesFid.Fill(fidCut,eiPMT)
		    eiPMTvals[jj] += eiPMT
				
	    elif locj == 'S':
	        if ii == 'boulby':
	    	    print 'reading in ',s	
		    # iterate over fiducial cut
   	    	    for jj,fidCut in enumerate(fid):
 	        	# set out the conditions that define a detected event
	        	posGood        = "(pos_goodness>%f)" %(float(arguments["-g"]))
	        	n9Good         = "(n9>%f)" %(float(arguments["--minN9"]))
	        	totEvtStr      = "all_ev_tot == all_ev"
			tt             = t.Draw("pe",totEvtStr,"goff") # total number of events simulated
	        	boulbyIBDRate  = rates["boulby_S"] # boulby rate per day
	    		recoFIDstring  = "(sqrt(pow(posReco.X(),2) + pow(posReco.Y(),2))/1000. <%f && sqrt(pow(posReco.Z(),2))/1000. <%f)"%(fidCut,fidCut) # event is reconstructed within the inner volume

			# find the pe efficiency
			drawArg		= "n9>>hei(5000,0,500)"
		 	ei 		= t.Draw(drawArg,"%s && %s && %s" %(posGood,n9Good,recoFIDstring),"goff") 
			Hs		= t.GetHistogram()
			if Hs != None:
  			    pe_cut 	= Hs.Integral(0,90)
			    pe_tot	= Hs.Integral(0,5000)
			else:
			    pe_cut	= 1.0
			    pe_tot	= 1.0
			if pe_tot != 0:
			    pe_eff	= 1.0 - pe_cut/pe_tot
			else:
			    pe_eff      = 0

	 	        eiIBD = t.Draw("reco_r","%s && %s && %s"%(recoFIDstring,posGood,n9Good),"goff") #get no of detected events
		        eiIBD *= 1./tt*boulbyIBDRate*pe_eff #scale by rate and reconstruction efficiency
			histIBDfid.Fill(fidCut,eiIBD)
			eiIBDvals.append(eiIBD)

    # now let's draw the plots of PMT events
    f_root.cd()

    # draw the first stack - histogram over full tank 
    c1 = TCanvas('c1','PMT singles')
    c1.cd()
    histSinglesFid.SetTitle("PMT singles")
    histSinglesFid.Draw("hist")
    histSinglesFid.SetFillColor(9)
    histSinglesFid.SetFillStyle(3001)
    histSinglesFid.GetXaxis().SetTitle("Fiducial cut (m)")
    histSinglesFid.GetYaxis().SetTitle("Events")

    # draw the second stack - close-up
    c2 = TCanvas('c2','IBD events')
    c2.cd()
    histIBDfid.Draw("hist")
    histIBDfid.SetFillColor(9)
    histIBDfid.SetFillStyle(3001)
    histIBDfid.SetTitle("IBD events")
    histIBDfid.GetXaxis().SetTitle("Fiducial cut (m)")
    histIBDfid.GetYaxis().SetTitle("Events")

    wait() 

    return histSinglesFid.Write(),histIBDfid.Write() #,h4.Write() 

# get the accidentals and ibd events as a function of fiducial cut
def extractPMTAccidentalsFidCutHistograms():  #pmtanalysis6 

    d,iso,loc,coverage,coveragePCT = loadSimulationParameters()
    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)
    _str = "ntuple_root_files%s/SBR.root" %(additionalString)
    f_root = TFile(_str,"recreate")
    #get the parameters from the watch arguments
    pmtDist  		= float(arguments["--tankRadius"])-float(arguments['--steelThick'])-float(arguments["--shieldThick"])
    pmtDistZ 		= float(arguments["--halfHeight"])-float(arguments['--steelThick'])-float(arguments["--shieldThick"])

    timeScale           = arguments["--timeScale"]
    inFilePrefix        = arguments["--ft"] 
    timeCut             = float(arguments["-t"])*1e3
    distCut             = float(arguments["-d"])
    
    #import the analysis parameters from load.py
    parameters  	= loadAnalysisParameters(timeScale)
    boulbyRate		= parameters[9]
    rates       	= parameters[11]
    mass        	= parameters[10]
    pc_num      	= parameters[12]
    timeS		= parameters[8]  
 
    cover = "25pct"

    hfidcut1 = TH1F("hfidcut1","",55,3.,8.5)
    hfidcut2 = TH1F("hfidcut2","",55,3.,8.5)
    hfidcut3 = TH1F("hfidcut3","",55,3.,8.5)

    fid = np.arange(0,8.5,0.1)
    n   = len(fid)
    eiPMT=[0]*n

    #iterate over processes
    for j in range(len(iso)):
	for ii in d["%s"%(iso[int(j)])]:
            locj = loc[j]		
	    isotope = ii
	    s = "ntuple_root_files%s/%s_%s_%s_%s.root"%(additionalString,inFilePrefix,isotope,cover,locj)
	    rfile = TFile(s)
	    t   = rfile.Get('data')
            
	    # set out the conditions that define an event as a relevant background
	    posGood		= "(pos_goodness>%f)" %(float(arguments["-g"]))
	    n9Good		= "(n9>%f)" %(float(arguments["--minN9"]))
	    totEvtStr	    	= "all_ev_tot == all_ev"

	    if locj == "PMT":

		print "Reading in ",s
		
		for jj,fidCut in enumerate(fid):
   	    	    recoFIDstring	= "(sqrt(pow(posReco.X(),2) + pow(posReco.Y(),2))/1000. <%f && sqrt(pow(posReco.Z(),2))/1000. <%f)"%(fidCut,fidCut) # event is reconstructed within the inner volume
		    tt = t.Draw("pe",totEvtStr,"goff") # number of events detected / total number of events simulated
		    er = float(rates["%s_%s"%(ii,locj)]) # activity (Bq kg^-1)
			
		    #  find the pe efficiency
		    drawArg		= "n9>>hei(5000,0,500)"
		    ei 		= t.Draw(drawArg,"%s && %s && %s" %(posGood,n9Good,recoFIDstring),"goff") 
		    Hs		= t.GetHistogram()
		    if Hs != None:
			pe_cut 	= Hs.Integral(0,90)
		    	pe_tot	= Hs.Integral(0,5000)
		    else: 
		        pe_cut	= 1.0
		        pe_tot	= 1.0
		    if pe_tot  != 0:
			pe_eff	= 1.0 - pe_cut/pe_tot
		    else:
			pe_eff  = 0
				
		    # get number of single events in the FV
                    N  		= t.Draw("pe:posReco.X():posReco.Y():posReco.Z()","%s && %s && %s" %(recoFIDstring,posGood,n9Good),"goff")
                    x    		= t.GetV2()
                    y    		= t.GetV3()
                    z    		= t.GetV4()
		    # find the number of single events which within ?m of each other
		    cnt = 0.
                    for index in range(N-1):
                        _rad = sqrt(power(x[index]-x[index+1],2)+power(y[index]-y[index+1],2)+power(z[index]-z[index+1],2))/1000.
                        if _rad < float(arguments["-d"]):
                            cnt+=1

		    eiPMT_tmp = cnt/tt*er*pc_num[cover]*mass*pe_eff #scale by reconstruction efficiency, total radioactivity rate, and pe efficiency
		    eiPMT_tmp *= eiPMT_tmp*0.0001/(24.*3600.)
		    eiPMT[jj] +=eiPMT_tmp
#			eiPMT = 0.0001/timeS*eiPMT_tmp*eiPMT_tmp
	
		    # first let's create histograms of single events in the PMTs
		    # get histogram of events at position r
		    hfidcut1.Fill(fidCut,eiPMT_tmp) 
		    
	    elif locj == 'S':
		
		if ii == 'boulby':
		    print "Reading in ",s
		    
   	    	    for jj,fidCut in enumerate(fid):
 	        	# set out the conditions that define a detected event
	        	posGood        = "(pos_goodness>%f)" %(float(arguments["-g"]))
	        	n9Good         = "(n9>%f)" %(float(arguments["--minN9"]))
	        	totEvtStr      = "all_ev_tot == all_ev"
			tt             = t.Draw("pe",totEvtStr,"goff") # total number of events simulated
	        	boulbyIBDRate  = rates["boulby_S"] # boulby rate per day
	    		recoFIDstring  = "(sqrt(pow(posReco.X(),2) + pow(posReco.Y(),2))/1000. <%f && sqrt(pow(posReco.Z(),2))/1000. <%f)"%(fidCut,fidCut) # event is reconstructed within the inner volume

			# find the pe efficiency
			drawArg		= "n9>>hei(5000,0,500)"
		 	ei 		= t.Draw(drawArg,"%s && %s && %s" %(posGood,n9Good,recoFIDstring),"goff") 
			Hs		= t.GetHistogram()
			if Hs != None:
  			    pe_cut 	= Hs.Integral(0,90)
			    pe_tot	= Hs.Integral(0,5000)
			else:
			    pe_cut	= 1.0
			    pe_tot	= 1.0
			if pe_tot != 0:
			    pe_eff	= 1.0 - pe_cut/pe_tot
			else:
			    pe_eff      = 0

	 	        eiIBD = t.Draw("reco_r","%s && %s && %s"%(recoFIDstring,posGood,n9Good),"goff") #get no of detected events
		        eiIBD *= 1./tt*boulbyIBDRate*pe_eff #scale by rate and reconstruction efficiency
		    	hfidcut2.Fill(fidCut,eiIBD)
        	    	s2b = eiIBD/sqrt(eiPMT[jj]+eiIBD)
	  		hfidcut3.Fill(fidCut,s2b)

    f_root.cd()

    # draw the pmt accidentals 
    c1 = TCanvas('c1','PMT accidentals')
    c1.cd()
    hfidcut1.Draw("hist")
    hfidcut1.SetTitle('PMT accidentals per day')
    hfidcut1.GetXaxis().SetTitle("Fiducial cut (m from centre)")
    hfidcut1.GetYaxis().SetTitle("Events per day")
    hfidcut1.SetFillColor(9)
    hfidcut1.SetFillStyle(3001)

    #draw the IBD events
    c2 = TCanvas('c2','Boulby IBD events')
    c2.cd()
    hfidcut2.Draw("hist")
    hfidcut2.SetTitle('Boulby IBD events per day')
    hfidcut2.GetXaxis().SetTitle('Fiducial cut (m from centre)')
    hfidcut2.GetYaxis().SetTitle('Events per day')
    hfidcut2.SetFillColor(9)
    hfidcut2.SetFillStyle(3001)
    
    #draw the signal to PMT background rate
    c3 = TCanvas('c3','Signal to PMT background')
    c3.cd()
    hfidcut3.Draw("hist")
    hfidcut3.SetTitle('Signal to PMT background')
    hfidcut3.GetXaxis().SetTitle('Fiducial cut (m from centre)')
    hfidcut3.GetYaxis().SetTitle('s/sqrt(s+b)')
    hfidcut3.SetFillColor(9)
    hfidcut3.SetFillStyle(3001)
    
    wait()

    return hfidcut1.Write(),hfidcut2.Write(),hfidcut3.Write() 
    
'''
   
