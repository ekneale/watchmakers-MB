// M.B
#include <iostream>
#include <iomanip>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TClass.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TVector3.h>
#include <vector>
#include <TRandom3.h>
#include <TMath.h>

TRandom3 randNum;


// The purpose of this code is to do a first pass analysis on the rat-pac
// format data to extract a "triggered" dataset. A timestamp is generated
// that emulate the data; the process is assumed to be poisson. This code
// is applied on each physics process independently, a code to merge all the
// timestamps for different physics process will be applied at a later stage.
//
// This code requires two input that have no defaults, the rate of events, and
// a interaction code [PMT(code 2) 214Bi (Z:83 A:214) = 20830214].
//
// If the event could have triggered either a prompt or delayed event in
// WATCHMAN, the event is recorded and the associated booleen is set to on.
// This DOES NOT imply that the events is a prompt or delayed, just that it
// has the possibility to be one. This step is the reduce the size of the dataset
// by emulating a trigger condition.
//
// Finally, a root file with the PDFs for the different analysis variable is
// generated:
// 1) non-trigger PDF : all events that left some signature in the detector
// 2) prompt PDF: distributions for events that pass the prompt trigger requierements
// 3) delayed PDF: distributions for events that pass the delayed trigger requierements
//


//#include <libRATEvent.h>

int pass1Trigger(const char *file, double rate, int code,int nPMT ,const char *outfile = "null",
double pmtBoundR = 6.42,double pmtBoundZ = 6.42,
double tankBoundR = 8.02635,double tankBoundZ = 8.02635,
double nhit_min_p = 3., double good_pos_p = 0.1, double good_dir_p = 0.1,
double pe_p = 5.5, double n9_p = 5, double n9over_nhit_p = 0.0001,
double nhit_min_d = 3., double good_pos_d = 0.222, double good_dir_d = 0.1,
double pe_d = 28.7, double n9_d = 5, double n9over_nhit_d = 0.0001
) {

  // Define the incoming out outgoing Trees
  TFile *f = new TFile(file);
  TTree *tree = (TTree*) f->Get("T");
  if (tree==0x0){
    return -1;
  }
  TFile *f_out;
  printf("outfile: %s\n",outfile);
  if (TString(outfile) == TString("null")) {
    outfile = Form("pass1_%s",f->GetName());
    f_out = new TFile(Form("pass1_%s",f->GetName()),"Recreate");
  }else{
    f_out = new TFile(outfile,"Recreate");
  }

  RAT::DS::Root *rds = new RAT::DS::Root();
  tree->SetBranchAddress("ds", &rds);
  int nEvents = tree->GetEntries();

  //Define all the analysis parameters
  Double_t totPE = 0.0,goodness,dirGoodness,timeTmp,n9;
  Double_t newX,newY,newZ,dirX,dirY,dirZ;

  Int_t ibd=0,es=0,cc=0,icc=0,nc=0,old_singal,evt;
  Double_t cosTheta,cosThetaSN,cosThetaSNIBD, local_time,local_time_tmp,delta_time,mc_energy;

  int sub_ev=0,cnt_all = 0,cnt_p =0,cnt_d = 0;
  Int_t subevents = 0,cnt,cntLoop;
  Int_t candidate,_candidate;
  Double_t timeDiff,r,z,r_t,z_t;

  TVector3 posTruth,posReco,dirTruth,dirNu,dirReco,dirIBD,pos1,pos2;

  Int_t maybePrompt, maybeDelay;
  Int_t totNHIT = 0,nhits =0,od_hit=0;
  Int_t particleCountMC;
  ULong64_t             timestamp;
  Int_t                 timestamp_ns,timestamp_s;
  Double_t              x,y,z,u,v,w;
  Double_t              mcX,mcY,mcZ,mcU,mcV,mcW,mcP,p2ToB,p2W,closestPMT;
  Double_t              timeLapse;

  TTree *data = new TTree("data","low-energy detector triggered events");

  data->Branch("maybePrompt",&maybePrompt,"maybePrompt/I");
  data->Branch("maybeDelay",&maybeDelay,"maybeDelay/I");
  data->Branch("timestamp",&timestamp,"timestamp/L");
  data->Branch("timestamp_ns",&timestamp_ns,"timestamp_ns/I");
  data->Branch("timestamp_s",&timestamp_s,"timestamp_s/I");

  data->Branch("nhit",&nhits,"nhit/I");
  data->Branch("id_plus_dr_hit",&totNHIT,"id_plus_dr_hit/I");//Inner detector plus dark rate hits
  data->Branch("od_hit",&od_hit,"od_hit/I");//Inner detector plus dark rate hits

  data->Branch("pe",&totPE,"pe/D");
  data->Branch("n9",&n9,"n9/D");
  data->Branch("good_pos",&goodness,"good_pos/D");
  data->Branch("good_dir",&dirGoodness,"good_dir/D");
  data->Branch("x",&newX,"x/D");
  data->Branch("y",&newY,"y/D");
  data->Branch("z",&newZ,"z/D");
  data->Branch("u",&dirX,"u/D");
  data->Branch("v",&dirY,"v/D");
  data->Branch("w",&dirZ,"w/D");
  data->Branch("particleCountMC",&particleCountMC ,"particleCountMC/I");
  data->Branch("mc_energy",&mc_energy,"mc_energy/D");
  data->Branch("mcx",&mcX,"mcx/D");
  data->Branch("mcy",&mcY,"mcy/D");
  data->Branch("mcz",&mcZ,"mcz/D");
  data->Branch("mcu",&mcU,"mcu/D");
  data->Branch("mcv",&mcV,"mcv/D");
  data->Branch("mcw",&mcW,"mcw/D");
  data->Branch("code",&code,"code/I");
  // data->Branch("p2W",&p2W,"p2W/D");//Proximity to PMT wall
  // data->Branch("p2ToB",&p2ToB,"p2ToB/D");//Proximity to Top or Bottom PMT (closest)
  data->Branch("closestPMT",&closestPMT,"closestPMT/D");//Proximity to PMT wall

  TTree *nodata =  new TTree("nodata","low-energy detector untriggered events");
  nodata->Branch("timestamp",&timestamp,"timestamp/L");
  nodata->Branch("timestamp_ns",&timestamp_ns,"timestamp_ns/I");
  nodata->Branch("timestamp_s",&timestamp_s,"timestamp_s/I");
  nodata->Branch("nhit",&nhits,"nhit/I");
  nodata->Branch("id_plus_dr_hit",&totNHIT,"id_plus_dr_hit/I");//Inner detector plus dark rate hits
  nodata->Branch("od_hit",&od_hit,"od_hit/I");//Inner detector plus dark rate hits
  nodata->Branch("pe",&totPE,"pe/D");
  nodata->Branch("n9",&n9,"n9/D");
  nodata->Branch("good_pos",&goodness,"good_pos/D");
  nodata->Branch("good_dir",&dirGoodness,"good_dir/D");
  nodata->Branch("x",&newX,"x/D");
  nodata->Branch("y",&newY,"y/D");
  nodata->Branch("z",&newZ,"z/D");
  nodata->Branch("u",&dirX,"u/D");
  nodata->Branch("v",&dirY,"v/D");
  nodata->Branch("w",&dirZ,"w/D");
  nodata->Branch("particleCountMC",&particleCountMC ,"particleCountMC/I");
  nodata->Branch("mc_energy",&mc_energy,"mc_energy/D");
  nodata->Branch("mcx",&mcX,"mcx/D");
  nodata->Branch("mcy",&mcY,"mcy/D");
  nodata->Branch("mcz",&mcZ,"mcz/D");
  nodata->Branch("mcu",&mcU,"mcu/D");
  nodata->Branch("mcv",&mcV,"mcv/D");
  nodata->Branch("mcw",&mcW,"mcw/D");
  nodata->Branch("code",&code,"code/I");


  TTree *runSummary =  new TTree("runSummary","mc run summary");
  runSummary->Branch("nEvents",&nEvents,"nEvents/I");
  Int_t subEventTally[20] = {};
  runSummary->Branch("subEventTally",subEventTally,"subEventTally[20]/I");
  runSummary->Branch("rateHZ",&rate,"rateHz/D");
  runSummary->Branch("inputFile",file,"inputFile/C");
  runSummary->Branch("outputFile",outfile,"outputFile/C");
  runSummary->Branch("code",&code,"code/I");
  runSummary->Branch("runEndTime",&timestamp,"runEndTime/L");//Record last time stamp of the file
  runSummary->Branch("runEndTime_ns",&timestamp_ns,"runEndTime_ns/I");//Record last time stamp of the file
  runSummary->Branch("runEndTime_s",&timestamp_s,"runEndTime_s/I");//Record last time stamp of the file
  runSummary->Branch("potential_prompts",&cnt_p,"potential_prompts/I");
  runSummary->Branch("potential_delayed",&cnt_d,"potential_delayed/I");
  Double_t eff_p,eff_d;
  runSummary->Branch("eff_prompts",&eff_p,"eff_prompts/D");
  runSummary->Branch("eff_delayed",&eff_d,"eff_delayed/D");
  runSummary->Branch("pmtBoundR",&pmtBoundR,"pmtBoundR/D");
  runSummary->Branch("pmtBoundZ",&pmtBoundZ,"pmtBoundZ/D");
  runSummary->Branch("tankBoundR",&tankBoundR,"tankBoundR/D");
  runSummary->Branch("tankBoundZ",&tankBoundZ,"tankBoundZ/D");
  runSummary->Branch("nPMT",&nPMT,"nPMT/I");
  runSummary->Branch("nhit_min_p",&nhit_min_p,"nhit_min_p/D");
  runSummary->Branch("good_pos_p",&good_pos_p,"good_pos_p/D");
  runSummary->Branch("good_dir_p",&good_dir_p,"good_dir_p/D");
  runSummary->Branch("pe_p",&pe_p,"pe_p/D");
  runSummary->Branch("n9_p",&n9_p,"n9_p/D");
  runSummary->Branch("n9over_nhit_p",&n9over_nhit_p,"n9over_nhit_p/D");
  runSummary->Branch("nhit_min_d",&nhit_min_d,"nhit_min_d/D");
  runSummary->Branch("good_pos_d",&good_pos_d,"good_pos_d/D");
  runSummary->Branch("good_dir_d",&good_dir_d,"good_dir_d/D");
  runSummary->Branch("pe_d",&pe_d,"pe_d/D");
  runSummary->Branch("n9_d",&n9_d,"n9_d/D");
  runSummary->Branch("n9over_nhit_d",&n9over_nhit_d,"n9over_nhit_d/D");

  vector <double> subeventInfo;
  vector<vector <double> > eventInfo;
  RAT::DS::MC *mc;
  RAT::DS::MCParticle *prim;
  //    Int_t particleCountMC;
  for (evt = 0; evt < nEvents; evt++) {

    tree->GetEntry(evt);

    // Obtain MC truth information
    mc                          = rds->GetMC();
    particleCountMC             = mc->GetMCParticleCount();
    //Get the direction of the neutrino. Saved as last particle
    prim                        = mc->GetMCParticle(0);
    mc_energy                   = prim->ke;
    posTruth                    = prim->GetPosition();
    dirNu                       = prim->GetMomentum();

    mcX 			    = prim->GetPosition().X()/1000.;
    mcY                         = prim->GetPosition().Y()/1000.;
    mcZ                         = prim->GetPosition().Z()/1000.;
    mcU                         = prim->GetMomentum().X();
    mcV                         = prim->GetMomentum().Y();
    mcW                         = prim->GetMomentum().Z();
    mcP		            = sqrt(mcU**2+mcV**2+mcW**2);
    mcU			    = mcU/mcP;
    mcV                         = mcV/mcP;
    mcW                         = mcW/mcP;

    timeLapse 		    = findNextTime(rate);
    timestamp                   += timeLapse;
    timestamp_ns                = int((float(timestamp/1.0e9)-int(timestamp/1.0e9))*1.0e9);//For some reason Modulo does not work
    timestamp_s                = int(timestamp/1e9);
    // printf("%llu %10d %10d \n",timestamp,timestamp_ns,timestamp_s);
    //Find out how many subevents:
    subevents                   = rds->GetEVCount();
    subEventTally[subevents]+=1;

    if (subevents==0){
      nhits                    = 0;//ev->Nhits();
      //RAT::DS::BonsaiFit *pb  = ev->GetBonsaiFit();
      goodness                = -1.;//pb->GetGoodness();
      dirGoodness             = -1.;//pb->GetDirGoodness();

      totNHIT                 = -1;//pb->GetIDHit();
      totPE                   = -1.;//pb->GetIDCharge();
      od_hit                  = -1;//pb->GetODHit();
      //dirReco                 = pb->GetDirection();
      n9                      = -1.;//pb->GetN9();
      newX                    = -1.;//posReco.X();
      newY                    = -1.;//posReco.Y();
      newZ                    = -1.;//posReco.Z();
      dirX                    = -1.;//dirReco.X();
      dirY                    = -1.;//dirReco.Y();
      dirZ                    = -1.;//dirReco.Z();
      nodata->Fill();
    }else{
      cnt = 0;

      subeventInfo.resize(0);
      eventInfo.resize(0);

      cnt_all = cnt = 0 ;

      for (int k = 0; k<subevents; k++) {
        RAT::DS::EV *ev         = rds->GetEV(k);
        nhits                    = ev->Nhits();
        RAT::DS::BonsaiFit *pb  = ev->GetBonsaiFit();
        goodness                = pb->GetGoodness();
        dirGoodness             = pb->GetDirGoodness();
        posReco                 = pb->GetPosition();
        totNHIT                 = pb->GetIDHit();
        totPE                   = pb->GetIDCharge();
        od_hit                  = pb->GetODHit();
        dirReco                 = pb->GetDirection();
        n9                      = pb->GetN9();
        newX                    = posReco.X()/1000.;
        newY                    = posReco.Y()/1000.;
        newZ                    = posReco.Z()/1000.;
        dirX                    = dirReco.X();
        dirY                    = dirReco.Y();
        dirZ                    = dirReco.Z();


        timeTmp = ev->GetCalibratedTriggerTime(); // 0 for first subevent, detlta for all others
        timestamp+=timeTmp;
        //printf("%d %d %d %d\n",timestamp,timeTmp,timeLapse,k);
        /*r = sqrt(pow(posReco.X(),2)+ pow(posReco.Y(),2))/1000.;
        z = posReco.Z()/1000.;
        */
        // Currentlt saving n9over_nhit_p/d for pass3
        if ( nhits > nhit_min_p && totPE > pe_p && n9 > n9_p && goodness > good_pos_p && dirGoodness > good_dir_p && float(n9)/float(nhits)> n9over_nhit_p && newX != -99.999) {
          maybePrompt = 1;
          p2W = pmtBoundR-sqrt(newX**2+newY**2);
          p2ToB = pmtBoundZ-sqrt(newZ**2);
          closestPMT = TMath::Min(p2W,p2ToB);
          cnt_p+=1;
        }
        if ( nhits > nhit_min_d && totPE > pe_d && n9 > n9_d && goodness > good_pos_d && dirGoodness > good_dir_d && float(n9)/float(nhits)> n9over_nhit_d && newX != -99.999) {
          maybeDelay  = 1;
          p2W = pmtBoundR-sqrt(newX**2+newY**2);
          p2ToB = pmtBoundZ-sqrt(newZ**2);
          closestPMT = TMath::Min(p2W,p2ToB);
          cnt_d+=1;
        }
        if (!maybePrompt && !maybeDelay){
          nodata->Fill();
        }
        if (maybePrompt || maybeDelay){
          data->Fill();
          maybePrompt = maybeDelay = 0;
        }

      }//for (int k = 0; k<subevents; k++)
    }//else
  }//for (evt = 0; evt < nEvents; evt++)
  data->Write();
  nodata->Write();
  eff_p = float(cnt_p)/float(nEvents);
  eff_d = float(cnt_d)/float(nEvents);

  runSummary->Fill();
  runSummary->Write();
  return 0;

}//int pdfGenerator(..

Double_t findNextTime(Double_t rate){
    return randNum.Exp(1./rate)*1.0e9;
}
