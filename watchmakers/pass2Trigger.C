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

int pass2Trigger(const char *cumulativeFile, const char *addfile,int first = 0) {
  Double_t totPE = 0.0,goodness,dirGoodness,timeTmp,n9;
  Double_t newX,newY,newZ,dirX,dirY,dirZ;

  Int_t ibd=0,es=0,cc=0,icc=0,nc=0,old_singal,evt;
  Double_t cosTheta,cosThetaSN,cosThetaSNIBD, local_time,local_time_tmp,delta_time,mc_energy;

  int sub_ev=0,cnt_all = 0,cnt_p =0,cnt_d = 0;
  Int_t subevents = 0,cnt,cntLoop;
  Int_t candidate,_candidate;
  Double_t timeDiff,r,z,r_t,z_t;

  TVector3 posTruth,posReco,dirTruth,dirNu,dirReco,dirIBD,pos1,pos2;
  Double_t eff_p,eff_d;

  Int_t maybePrompt, maybeDelay,nEvents;
  Int_t totNHIT = 0,nhits =0,od_hit=0;
  Int_t particleCountMC;
  Long64_t             timestamp,runEndTime;
  Int_t                 timestamp_ns,timestamp_s,runEndTime_ns,runEndTime_s,code;
  Double_t              x,y,z,u,v,w,rate;
  Double_t              mcX,mcY,mcZ,mcU,mcV,mcW,mcP,closestPMT;
  Double_t              timeLapse;
  const char *file,outfile;
  int nPMT;
  double pmtBoundR, pmtBoundZ, tankBoundR ,tankBoundZ;
  double nhit_min_p,good_pos_p,good_dir_p,pe_p,n9_p,n9over_nhit_p,nhit_min_d,good_pos_d ,good_dir_d,pe_d,n9_d,n9over_nhit_d;


  Int_t subEventTally[20] = {};


  // Load in the new file (or daughter file)
  TFile *af = new TFile(addfile);
  TTree *dataDaughter = (TTree*) af->Get("data");
  dataDaughter->SetBranchAddress("maybePrompt",&maybePrompt);
  dataDaughter->SetBranchAddress("maybeDelay",&maybeDelay);
  dataDaughter->SetBranchAddress("timestamp",&timestamp);
  dataDaughter->SetBranchAddress("timestamp_ns",&timestamp_ns);
  dataDaughter->SetBranchAddress("timestamp_s",&timestamp_s);
  dataDaughter->SetBranchAddress("nhit",&nhits);
  dataDaughter->SetBranchAddress("id_plus_dr_hit",&totNHIT);//Inner detector plus dark rate hits
  dataDaughter->SetBranchAddress("od_hit",&od_hit);//Inner detector plus dark rate hits
  dataDaughter->SetBranchAddress("pe",&totPE);
  dataDaughter->SetBranchAddress("n9",&n9);
  dataDaughter->SetBranchAddress("good_pos",&goodness);
  dataDaughter->SetBranchAddress("good_dir",&dirGoodness);
  dataDaughter->SetBranchAddress("x",&newX);
  dataDaughter->SetBranchAddress("y",&newY);
  dataDaughter->SetBranchAddress("z",&newZ);
  dataDaughter->SetBranchAddress("u",&dirX);
  dataDaughter->SetBranchAddress("v",&dirY);
  dataDaughter->SetBranchAddress("w",&dirZ);
  dataDaughter->SetBranchAddress("closestPMT",&closestPMT);
  dataDaughter->SetBranchAddress("particleCountMC",&particleCountMC);
  dataDaughter->SetBranchAddress("mc_energy",&mc_energy);
  dataDaughter->SetBranchAddress("mcx",&mcX);
  dataDaughter->SetBranchAddress("mcy",&mcY);
  dataDaughter->SetBranchAddress("mcz",&mcZ);
  dataDaughter->SetBranchAddress("mcu",&mcU);
  dataDaughter->SetBranchAddress("mcv",&mcV);
  dataDaughter->SetBranchAddress("mcw",&mcW);
  dataDaughter->SetBranchAddress("code",&code);

  // TTree *nodataDaughter = (TTree*) af->Get("nodata");
  // nodataDaughter->SetBranchAddress("timestamp",&timestamp);
  // nodataDaughter->SetBranchAddress("timestamp_ns",&timestamp_ns);
  // nodataDaughter->SetBranchAddress("timestamp_s",&timestamp_s);
  // nodataDaughter->SetBranchAddress("nhit",&nhits);
  // nodataDaughter->SetBranchAddress("id_plus_dr_hit",&totNHIT);//Inner detector plus dark rate hits
  // nodataDaughter->SetBranchAddress("od_hit",&od_hit);//Inner detector plus dark rate hits
  // nodataDaughter->SetBranchAddress("pe",&totPE);
  // nodataDaughter->SetBranchAddress("n9",&n9);
  // nodataDaughter->SetBranchAddress("good_pos",&goodness);
  // nodataDaughter->SetBranchAddress("good_dir",&dirGoodness);
  // nodataDaughter->SetBranchAddress("x",&newX);
  // nodataDaughter->SetBranchAddress("y",&newY);
  // nodataDaughter->SetBranchAddress("z",&newZ);
  // nodataDaughter->SetBranchAddress("u",&dirX);
  // nodataDaughter->SetBranchAddress("v",&dirY);
  // nodataDaughter->SetBranchAddress("w",&dirZ);
  // nodataDaughter->SetBranchAddress("particleCountMC",&particleCountMC);
  // nodataDaughter->SetBranchAddress("mc_energy",&mc_energy);
  // nodataDaughter->SetBranchAddress("mcx",&mcX);
  // nodataDaughter->SetBranchAddress("mcy",&mcY);
  // nodataDaughter->SetBranchAddress("mcz",&mcZ);
  // nodataDaughter->SetBranchAddress("mcu",&mcU);
  // nodataDaughter->SetBranchAddress("mcv",&mcV);
  // nodataDaughter->SetBranchAddress("mcw",&mcW);
  // nodataDaughter->SetBranchAddress("code",&code);

  TTree *runSummaryDaughter = (TTree*) af->Get("runSummary");
  runSummaryDaughter->SetBranchAddress("nEvents",&nEvents);
  runSummaryDaughter->SetBranchAddress("subEventTally",subEventTally);//,"subEventTally[20]/I");
  runSummaryDaughter->SetBranchAddress("rateHZ",&rate);//,"rateHz/D");
  runSummaryDaughter->SetBranchAddress("inputFile",file);//,"inputFile/C");
  runSummaryDaughter->SetBranchAddress("outputFile",outfile);//,"outputFile/C");
  runSummaryDaughter->SetBranchAddress("code",&code);//,"code/I");
  runSummaryDaughter->SetBranchAddress("runEndTime",&runEndTime);//,"runEndTime/L");//Record last time stamp of the file
  runSummaryDaughter->SetBranchAddress("runEndTime_ns",&runEndTime_ns);//,"runEndTime_ns/I");//Record last time stamp of the file
  runSummaryDaughter->SetBranchAddress("runEndTime_s",&runEndTime_s);//,"runEndTime_s/I");//Record last time stamp of the file
  runSummaryDaughter->SetBranchAddress("potential_prompts",&cnt_p);//,"potential_prompts/I");
  runSummaryDaughter->SetBranchAddress("potential_delayed",&cnt_d);//,"potential_delayed/I");
  runSummaryDaughter->SetBranchAddress("eff_prompts",&eff_p);//,"eff_prompts/D");
  runSummaryDaughter->SetBranchAddress("eff_delayed",&eff_d);//,"eff_delayed/D");
  // runSummaryDaughter->SetBranchAddress("fidBoundR",&fidBoundR);//,"fidBoundR/D");
  // runSummaryDaughter->SetBranchAddress("fidBoundZ",&fidBoundZ);//,"fidBoundZ/D");
  runSummaryDaughter->SetBranchAddress("pmtBoundR",&pmtBoundR);//,"pmtBoundR/D");
  runSummaryDaughter->SetBranchAddress("pmtBoundZ",&pmtBoundZ);//,"pmtBoundZ/D");
  runSummaryDaughter->SetBranchAddress("tankBoundR",&tankBoundR);//,"tankBoundR/D");
  runSummaryDaughter->SetBranchAddress("tankBoundZ",&tankBoundZ);//,"tankBoundZ/D");
  runSummaryDaughter->SetBranchAddress("nPMT",&nPMT);//,"nPMT/I");
  runSummaryDaughter->SetBranchAddress("nhit_min_p",&nhit_min_p);//,,"nhit_min_p/D");
  runSummaryDaughter->SetBranchAddress("good_pos_p",&good_pos_p);//,,"good_pos_p/D");
  runSummaryDaughter->SetBranchAddress("good_dir_p",&good_dir_p);//,,"good_dir_p/D");
  runSummaryDaughter->SetBranchAddress("pe_p",&pe_p);//,,"pe_p/D");
  runSummaryDaughter->SetBranchAddress("n9_p",&n9_p);//,,"n9_p/D");
  runSummaryDaughter->SetBranchAddress("n9over_nhit_p",&n9over_nhit_p);//,,"n9over_nhit_p/D");
  runSummaryDaughter->SetBranchAddress("nhit_min_d",&nhit_min_d);//,,"nhit_min_d/D");
  runSummaryDaughter->SetBranchAddress("good_pos_d",&good_pos_d);//,,"good_pos_d/D");
  runSummaryDaughter->SetBranchAddress("good_dir_d",&good_dir_d);//,,"good_dir_d/D");
  runSummaryDaughter->SetBranchAddress("pe_d",&pe_d);//,,"pe_d/D");
  runSummaryDaughter->SetBranchAddress("n9_d",&n9_d);//,,"n9_d/D");
  runSummaryDaughter->SetBranchAddress("n9over_nhit_d",&n9over_nhit_d);//,,"n9over_nhit_d/D");

  if (dataDaughter==0x0||runSummaryDaughter==0x0){
    return -1;
  }

  // Does the parent file exist, if not create it
  FileStat_t buf;
  Int_t test =  gSystem->GetPathInfo(cumulativeFile,buf);
  if(test!=0 || first!=0){
    TFile *f = new TFile(cumulativeFile,"update");
    TTree *data = dataDaughter->CloneTree();
    // TTree *nodata = nodataDaughter->CloneTree();
    TTree *runSummary = runSummaryDaughter->CloneTree();
    data->Write("", TObject::kOverwrite);
    // nodata->Write("", TObject::kOverwrite);
    runSummary->Write("", TObject::kOverwrite);
    af->Close();
    return 0;
  }else{
    printf("File exists, extraction trees.\n");
    TFile *f = new TFile(cumulativeFile,"update");
    TTree *dataParent = (TTree*) f->Get("data");
    dataParent->SetBranchAddress("maybePrompt",&maybePrompt);
    dataParent->SetBranchAddress("maybeDelay",&maybeDelay);
    dataParent->SetBranchAddress("timestamp",&timestamp);
    dataParent->SetBranchAddress("timestamp_ns",&timestamp_ns);
    dataParent->SetBranchAddress("timestamp_s",&timestamp_s);
    dataParent->SetBranchAddress("nhit",&nhits);
    dataParent->SetBranchAddress("id_plus_dr_hit",&totNHIT);//Inner detector plus dark rate hits
    dataParent->SetBranchAddress("od_hit",&od_hit);//Inner detector plus dark rate hits
    dataParent->SetBranchAddress("pe",&totPE);
    dataParent->SetBranchAddress("n9",&n9);
    dataParent->SetBranchAddress("good_pos",&goodness);
    dataParent->SetBranchAddress("good_dir",&dirGoodness);
    dataParent->SetBranchAddress("x",&newX);
    dataParent->SetBranchAddress("y",&newY);
    dataParent->SetBranchAddress("z",&newZ);
    dataParent->SetBranchAddress("u",&dirX);
    dataParent->SetBranchAddress("v",&dirY);
    dataParent->SetBranchAddress("w",&dirZ);
    dataParent->SetBranchAddress("particleCountMC",&particleCountMC);
    dataParent->SetBranchAddress("mc_energy",&mc_energy);
    dataParent->SetBranchAddress("mcx",&mcX);
    dataParent->SetBranchAddress("mcy",&mcY);
    dataParent->SetBranchAddress("mcz",&mcZ);
    dataParent->SetBranchAddress("mcu",&mcU);
    dataParent->SetBranchAddress("mcv",&mcV);
    dataParent->SetBranchAddress("mcw",&mcW);
    dataParent->SetBranchAddress("code",&code);
    dataParent->SetBranchAddress("closestPMT",&closestPMT);
    //
    // TTree *nodataParent = (TTree*) f->Get("nodata");
    // nodataParent->SetBranchAddress("timestamp",&timestamp);
    // nodataParent->SetBranchAddress("timestamp_ns",&timestamp_ns);
    // nodataParent->SetBranchAddress("timestamp_s",&timestamp_s);
    // nodataParent->SetBranchAddress("nhit",&nhits);
    // nodataParent->SetBranchAddress("id_plus_dr_hit",&totNHIT);//Inner detector plus dark rate hits
    // nodataParent->SetBranchAddress("od_hit",&od_hit);//Inner detector plus dark rate hits
    // nodataParent->SetBranchAddress("pe",&totPE);
    // nodataParent->SetBranchAddress("n9",&n9);
    // nodataParent->SetBranchAddress("good_pos",&goodness);
    // nodataParent->SetBranchAddress("good_dir",&dirGoodness);
    // nodataParent->SetBranchAddress("x",&newX);
    // nodataParent->SetBranchAddress("y",&newY);
    // nodataParent->SetBranchAddress("z",&newZ);
    // nodataParent->SetBranchAddress("u",&dirX);
    // nodataParent->SetBranchAddress("v",&dirY);
    // nodataParent->SetBranchAddress("w",&dirZ);
    // nodataParent->SetBranchAddress("particleCountMC",&particleCountMC);
    // nodataParent->SetBranchAddress("mc_energy",&mc_energy);
    // nodataParent->SetBranchAddress("mcx",&mcX);
    // nodataParent->SetBranchAddress("mcy",&mcY);
    // nodataParent->SetBranchAddress("mcz",&mcZ);
    // nodataParent->SetBranchAddress("mcu",&mcU);
    // nodataParent->SetBranchAddress("mcv",&mcV);
    // nodataParent->SetBranchAddress("mcw",&mcW);
    // nodataParent->SetBranchAddress("code",&code);
    TTree *runSummaryParent = (TTree*) f->Get("runSummary");
    runSummaryParent->SetBranchAddress("nEvents",&nEvents);
    runSummaryParent->SetBranchAddress("subEventTally",subEventTally);//,"subEventTally[20]/I");
    runSummaryParent->SetBranchAddress("rateHZ",&rate);//,"rateHz/D");
    runSummaryParent->SetBranchAddress("inputFile",file);//,"inputFile/C");
    runSummaryParent->SetBranchAddress("outputFile",outfile);//,"outputFile/C");
    runSummaryParent->SetBranchAddress("code",&code);//,"code/I");
    runSummaryParent->SetBranchAddress("runEndTime",&runEndTime);//,"runEndTime/L");//Record last time stamp of the file
    runSummaryParent->SetBranchAddress("runEndTime_ns",&runEndTime_ns);//,"runEndTime_ns/I");//Record last time stamp of the file
    runSummaryParent->SetBranchAddress("runEndTime_s",&runEndTime_s);//,"runEndTime_s/I");//Record last time stamp of the file
    runSummaryParent->SetBranchAddress("potential_prompts",&cnt_p);//,"potential_prompts/I");
    runSummaryParent->SetBranchAddress("potential_delayed",&cnt_d);//,"potential_delayed/I");
    runSummaryParent->SetBranchAddress("eff_prompts",&eff_p);//,"eff_prompts/D");
    runSummaryParent->SetBranchAddress("eff_delayed",&eff_d);//,"eff_delayed/D");
    // runSummaryParent->SetBranchAddress("fidBoundR",&fidBoundR);//,"fidBoundR/D");
    // runSummaryParent->SetBranchAddress("fidBoundZ",&fidBoundZ);//,"fidBoundZ/D");
    runSummaryParent->SetBranchAddress("pmtBoundR",&pmtBoundR);//,"pmtBoundR/D");
    runSummaryParent->SetBranchAddress("pmtBoundZ",&pmtBoundZ);//,"pmtBoundZ/D");
    runSummaryParent->SetBranchAddress("tankBoundR",&tankBoundR);//,"tankBoundR/D");
    runSummaryParent->SetBranchAddress("tankBoundZ",&tankBoundZ);//,"tankBoundZ/D");
    runSummaryParent->SetBranchAddress("nPMT",&nPMT);//,"nPMT/I");
    runSummaryParent->SetBranchAddress("nhit_min_p",&nhit_min_p);//,,"nhit_min_p/D");
    runSummaryParent->SetBranchAddress("good_pos_p",&good_pos_p);//,,"good_pos_p/D");
    runSummaryParent->SetBranchAddress("good_dir_p",&good_dir_p);//,,"good_dir_p/D");
    runSummaryParent->SetBranchAddress("pe_p",&pe_p);//,,"pe_p/D");
    runSummaryParent->SetBranchAddress("n9_p",&n9_p);//,,"n9_p/D");
    runSummaryParent->SetBranchAddress("n9over_nhit_p",&n9over_nhit_p);//,,"n9over_nhit_p/D");
    runSummaryParent->SetBranchAddress("nhit_min_d",&nhit_min_d);//,,"nhit_min_d/D");
    runSummaryParent->SetBranchAddress("good_pos_d",&good_pos_d);//,,"good_pos_d/D");
    runSummaryParent->SetBranchAddress("good_dir_d",&good_dir_d);//,,"good_dir_d/D");
    runSummaryParent->SetBranchAddress("pe_d",&pe_d);//,,"pe_d/D");
    runSummaryParent->SetBranchAddress("n9_d",&n9_d);//,,"n9_d/D");
    runSummaryParent->SetBranchAddress("n9over_nhit_d",&n9over_nhit_d);//,,"n9over_nhit_d/D");



    if (dataParent==0x0||runSummaryParent==0x0){
      printf("Error in loading pass2 current file.");
      return -1;
    }
  }

  // TFile *f = new TFile(cumulativeFile,"update");
  TTree *runSummary = new TTree("runSummary","mc run summary");
  runSummary->Branch("nEvents",&nEvents,"nEvents/I");
  runSummary->Branch("subEventTally",subEventTally,"subEventTally[20]/I");
  runSummary->Branch("rateHZ",&rate,"rateHz/D");
  runSummary->Branch("inputFile",file,"inputFile/C");
  runSummary->Branch("outputFile",outfile,"outputFile/C");
  runSummary->Branch("code",&code,"code/I");
  runSummary->Branch("runEndTime",&runEndTime,"runEndTime/L");//Record last time stamp of the file
  runSummary->Branch("runEndTime_s",&runEndTime_s,"runEndTime_s/I");//Record last time stamp of the file
  runSummary->Branch("runEndTime_ns",&runEndTime_ns,"runEndTime_ns/I");//Record last time stamp of the file
  runSummary->Branch("potential_prompts",&cnt_p,"potential_prompts/I");
  runSummary->Branch("potential_delayed",&cnt_d,"potential_delayed/I");
  runSummary->Branch("eff_prompts",&eff_p,"eff_prompts/D");
  runSummary->Branch("eff_delayed",&eff_d,"eff_delayed/D");
  // runSummary->Branch("fidBoundR",&fidBoundR,"fidBoundR/D");
  // runSummary->Branch("fidBoundZ",&fidBoundZ,"fidBoundZ/D");
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

  TTree *data = new TTree("data","low-energy detector triggered events");
  data->Branch("maybePrompt",&maybePrompt,"maybePrompt/I");
  data->Branch("maybeDelay",&maybeDelay,"maybeDelay/I");
  data->Branch("timestamp",&timestamp,"timestamp/L");
  data->Branch("timestamp_s",&timestamp_s,"timestamp_s/I");
  data->Branch("timestamp_ns",&timestamp_ns,"timestamp_ns/I");
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
  data->Branch("closestPMT",&closestPMT,"closestPMT/D");
  data->Branch("particleCountMC",&particleCountMC ,"particleCountMC/I");
  data->Branch("mc_energy",&mc_energy,"mc_energy/D");
  data->Branch("mcx",&mcX,"mcx/D");
  data->Branch("mcy",&mcY,"mcy/D");
  data->Branch("mcz",&mcZ,"mcz/D");
  data->Branch("mcu",&mcU,"mcu/D");
  data->Branch("mcv",&mcV,"mcv/D");
  data->Branch("mcw",&mcW,"mcw/D");
  data->Branch("code",&code,"code/I");
  // data->Branch("timestamp_ns",&timestamp_ns,"timestamp_ns/I");
  // data->Branch("timestamp_s",&timestamp_s,"timestamp_s/I");

  // TTree *nodata =  new TTree("nodata","low-energy detector untriggered events");
  // nodata->Branch("timestamp",&timestamp,"timestamp/L");
  // nodata->Branch("timestamp_s",&timestamp_s,"timestamp_s/I");
  // nodata->Branch("timestamp_ns",&timestamp_ns,"timestamp_ns/I");
  // nodata->Branch("nhit",&nhits,"nhit/I");
  // nodata->Branch("id_plus_dr_hit",&totNHIT,"id_plus_dr_hit/I");//Inner detector plus dark rate hits
  // nodata->Branch("od_hit",&od_hit,"od_hit/I");//Inner detector plus dark rate hits
  // nodata->Branch("pe",&totPE,"pe/D");
  // nodata->Branch("n9",&n9,"n9/D");
  // nodata->Branch("good_pos",&goodness,"good_pos/D");
  // nodata->Branch("good_dir",&dirGoodness,"good_dir/D");
  // nodata->Branch("x",&newX,"x/D");
  // nodata->Branch("y",&newY,"y/D");
  // nodata->Branch("z",&newZ,"z/D");
  // nodata->Branch("u",&dirX,"u/D");
  // nodata->Branch("v",&dirY,"v/D");
  // nodata->Branch("w",&dirZ,"w/D");
  // nodata->Branch("particleCountMC",&particleCountMC ,"particleCountMC/I");
  // nodata->Branch("mc_energy",&mc_energy,"mc_energy/D");
  // nodata->Branch("mcx",&mcX,"mcx/D");
  // nodata->Branch("mcy",&mcY,"mcy/D");
  // nodata->Branch("mcz",&mcZ,"mcz/D");
  // nodata->Branch("mcu",&mcU,"mcu/D");
  // nodata->Branch("mcv",&mcV,"mcv/D");
  // nodata->Branch("mcw",&mcW,"mcw/D");
  // nodata->Branch("code",&code,"code/I");

  //    Int_t particleCountMC;

  Long64_t parentEndTime;
  for (int evt = 0; evt < runSummaryParent->GetEntries(); evt++) {
    runSummaryParent->GetEntry(evt);
    parentEndTime= runEndTime;
    runSummary->Fill();
  }
  for (int evt = 0; evt < runSummaryDaughter->GetEntries(); evt++) {
    runSummaryDaughter->GetEntry(evt);
    runEndTime+=parentEndTime;
    runEndTime_ns                = int((float(runEndTime/1.0e9)-int(runEndTime/1.0e9))*1.0e9);//For some reason Modulo does not work
    runEndTime_s                = int(runEndTime/1e9);
    runSummary->Fill();
  }

  for (int evt = 0; evt < dataParent->GetEntries(); evt++) {
    dataParent->GetEntry(evt);
    // printf("P: %d %d %d\n",evt,dataParent->GetEntries(),parentEndTime);

    data->Fill();
 }
  for (int evt = 0; evt < dataDaughter->GetEntries(); evt++) {
    dataDaughter->GetEntry(evt);
    // printf("S: %d %d %d\n",evt,dataDaughter->GetEntries(),parentEndTime);

    timestamp+=parentEndTime;
    timestamp_ns                = int((float(timestamp/1.0e9)-int(timestamp/1.0e9))*1.0e9);//For some reason Modulo does not work
    timestamp_s                = int(timestamp/1e9);
    data->Fill();
 }
 //    Int_t particleCountMC;
 // Choosing not to save the nodata tree
 // for (int evt = 0; evt < nodataParent->GetEntries(); evt++) {
 //   nodataParent->GetEntry(evt);
 //   nodata->Fill();
 // }
 // for (int evt = 0; evt < nodataDaughter->GetEntries(); evt++) {
 //   nodataDaughter->GetEntry(evt);
 //   timestamp+=parentEndTime;
 //   timestamp_ns                = int((float(timestamp/1.0e9)-int(timestamp/1.0e9))*1.0e9);//For some reason Modulo does not work
 //   timestamp_s                = int(timestamp/1e9);
 //   nodata->Fill();
 // }
 //


 f->cd();
 data->Write("", TObject::kOverwrite);
 // nodata->Write("", TObject::kOverwrite);
 runSummary->Write("", TObject::kOverwrite);
 f->Close();
 af->Close();
 return 0;

}//int pdfGenerator(..
