


struct SumTF1 { 
  
SumTF1(const std::vector<TF1 *> & flist) : fFuncList(flist) {}  
  double operator() (const double * x, const double *p) {
    double result = 0;
    for (unsigned int i = 0; i < fFuncList.size(); ++i) 
      result += fFuncList[i]->EvalPar(x,p); 
    return result; 
  } 
  
  std::vector<TF1*> fFuncList; 
  
};


struct EscapeTF1 { 
EscapeTF1(const std::vector<TF1 *> & flist) : fFuncList(flist) {}

  double operator() (const double * x, const double *p) {
    double result = 0;
    double xx = p[2]*(p[1]-x[0]); 
    for (unsigned int i = 0; i < fFuncList.size(); ++i) 
      result += fFuncList[i]->Eval(xx); 
    return  p[0]*result; 
  }
  std::vector<TF1*> fFuncList;   
};


struct GRespTF1 { //
GRespTF1(TF1* fpeak) : fPeak(fpeak) {}
  double dx = 1;

  double operator() (const double * x, const double *p) {
    //dynamically set params -> this should allow for fitting
    //fDR->SetParameter(2,p[0]);
    //for(int i=0; i<fPeak->GetNpar(); i++)
    //  fPeak->SetParameter(i, p[i-1]); 

    // find integration range
    double lo = -6*p[0];
    double hi =  6*p[0];
    
    //integrate
    double result = 0;
    for(double xx = lo; xx<hi; xx+=dx) 
      result += fPeak->Eval(x[0]+xx)*TMath::Gaus(xx,0,p[0],kTRUE); 

    return  result; 
  }
  TF1* fPeak;   
};



/*
Double_t LtzMinusFunc(TF1* fFunc, Double_t *x, Double_t *par){
  return par[0]*TMath::CauchyDist(x[0],par[1],par[2])-fFunc->EvalPar(x,p);
}
*/ 

void plotlines(){

  std::vector<TF1*> fAllPeaks_;


  TTree *t = new TTree("t", "snkas");
  t->ReadFile("snkas_reduced.csv", "energy/D:intensity/D:width/D");
  Double_t energy, intensity, width; 
  t->SetBranchAddress("energy",    &energy); 
  t->SetBranchAddress("intensity", &intensity);
  t->SetBranchAddress("width",     &width);


  // set characteristic intensity
  t->GetEntry(0); double intensity0 = intensity; 

  

  std::vector<TF1*> fKa_;
  std::vector<TF1*> fKb_;
  
  int npts = 50000; 
  float e_lo =  80000;
  float e_hi = 150000;

  
  
  for(int i=0; i<t->GetEntries();i++){
    t->GetEntry(i); 
    if(energy>26000){
      fKb_.push_back(new TF1(Form("fL%d",i),"[0]*TMath::CauchyDist(x,[1],[2])",e_lo,e_hi));
      fKb_.back()->SetParameters(intensity/intensity0, energy, width);
      fKb_.back()->SetNpx(npts); 	
      //fKb_.back()->SetLineColor(kYellow+i);
    }else{
      fKa_.push_back(new TF1(Form("fL%d",i),"[0]*TMath::CauchyDist(x,[1],[2])",e_lo,e_hi));
      fKa_.back()->SetParameters(intensity/intensity0, energy, width);
      fKa_.back()->SetNpx(npts); 	
      //fKa_.back()->SetLineColor(kAzure+i);
    }
  }

  TF1* fSum_kb = new TF1("fSum_kb",SumTF1(fKb_),e_lo,e_hi,0);
  fSum_kb->SetTitle("Tin Calculated Spectrum; Energy [eV]; Intensity [Arb]"); 
  fSum_kb->SetLineColor(kViolet); 
  fSum_kb->SetNpx(npts); 	

  //fAllPeaks_.push_back(fSum_kb); 
  
  TF1* fSum_ka = new TF1("fSum_ka",SumTF1(fKa_),e_lo,e_hi,0);
  fSum_ka->SetTitle("Tin Calculated Spectrum; Energy [eV]; Intensity [Arb]");
  fSum_ka->SetLineColor(kMagenta); 
  fSum_ka->SetNpx(npts); 	

  //fAllPeaks_.push_back(fSum_ka);

  //------------- Prepare Calibraiton Lines -------------//  
  float E122 = 122060.65;  
  float E136 = 136473.56;
  float I122 = 10000000.; 
  float I136 = 10.68/85.60; I136*=I122;  

  std::vector<TF1*> fCal_; 

  fCal_.push_back(new TF1(Form("fCal122"),
			  "[0]*TMath::CauchyDist(x,[1],[2])+.5*[3]*[0]*TMath::ATan(-[4]*(x-[1]))/TMath::Pi()",
			  e_lo,e_hi));
  fCal_.back()->SetParameters(I122, E122, .1, .00001, 5);
  fCal_.back()->SetLineColor(kGreen); 
  fCal_.back()->SetNpx(npts);
  fAllPeaks_.push_back(fCal_.back());

  
  fCal_.push_back(new TF1(Form("fCal136"),
			  "[0]*TMath::CauchyDist(x,[1],[2])+.5*[3]*[0]*TMath::ATan(-[4]*(x-[1]))/TMath::Pi()",
			  e_lo,e_hi));
  fCal_.back()->SetParameters(I136, E136, .1, .00001, 5);
  fCal_.back()->SetLineColor(kGreen+1); 
  fCal_.back()->SetNpx(npts);
  fAllPeaks_.push_back(fCal_.back());

  
  
  //------------- Prepare Escapes -------------//
  
  double iEscape = I122*.01; 
  
  TF1* fSum_Es_122_ka = new TF1("fSum_Es_122_ka", EscapeTF1(fKa_),e_lo,e_hi,3);
  fSum_Es_122_ka->SetParameters(iEscape, E122, 1);
  fSum_Es_122_ka->SetLineColor(kOrange); 
  fSum_Es_122_ka->SetNpx(npts);
  fAllPeaks_.push_back(fSum_Es_122_ka);

  TF1* fSum_Es_122_kb = new TF1("fSum_Es_122_kb", EscapeTF1(fKb_),e_lo,e_hi,3);
  fSum_Es_122_kb->SetParameters(iEscape, E122, 1);
  fSum_Es_122_kb->SetLineColor(kOrange+1); 
  fSum_Es_122_kb->SetNpx(npts); 	
  fAllPeaks_.push_back(fSum_Es_122_kb);

  TF1* fSum_Es_136_ka = new TF1("fSum_Es_136_ka", EscapeTF1(fKa_),e_lo,e_hi,3);
  fSum_Es_136_ka->SetParameters(iEscape*I136/I122, E136, 1);
  fSum_Es_136_ka->SetLineColor(kOrange+2); 
  fSum_Es_136_ka->SetNpx(npts);
  fAllPeaks_.push_back(fSum_Es_136_ka);

  TF1* fSum_Es_136_kb = new TF1("fSum_Es_136_kb", EscapeTF1(fKb_),e_lo,e_hi,3);
  fSum_Es_136_kb->SetParameters(iEscape*I136/I122, E136, 1);
  fSum_Es_136_kb->SetLineColor(kOrange+3); 
  fSum_Es_136_kb->SetNpx(npts); 	
  fAllPeaks_.push_back(fSum_Es_136_kb);

  //------------- Make Sum Spectrum -------------//    
  TF1* fSum = new TF1("fSum",SumTF1(fAllPeaks_),e_lo,e_hi,0);
  fSum->SetTitle("Sum Spectrum; Energy [eV]; Intensity [Arb]");
  fSum->SetLineColor(kRed); 
  fSum->SetNpx(npts); 	


  double sigma = 40;     
  TF1 *fSim = new TF1("fSim", GRespTF1(fSum) , e_lo, e_hi, 1);
  fSim->SetParameter(0,sigma); //sigma  = XX eV
  fSim->SetLineColor(kBlue); 
  fSim->SetNpx(10000);
  
  
  //------------- Draw -------------//  
  TCanvas* c_temp = new TCanvas("c_temp","c_temp", 1600, 900);
  fSum->Draw(); 
  fSim->Draw("SAME");
  fCal_[0]->Draw("SAME"); 
  fCal_[1]->Draw("SAME");
  //fSum_ka->Draw("SAME");
  //fSum_kb->Draw("SAME");

  fSum_Es_122_ka->Draw("SAME"); 
  fSum_Es_122_kb->Draw("SAME"); 
  fSum_Es_136_ka->Draw("SAME"); 
  fSum_Es_136_kb->Draw("SAME"); 
 
  auto legend = new TLegend(.3,0.70,0.5,0.9);
  legend->AddEntry(fSum,    "Sum",    "l"); 
  legend->AddEntry(fSim,    "Sum x Det. Resp.",    "l"); 
  legend->AddEntry(fCal_[0],    "Calibration",    "l"); 
  //legend->AddEntry(fSum_ka,     "Sn-Ka",    "l"); 
  // legend->AddEntry(fSum_kb,     "Sn-Kb",    "l"); 
  //  legend->AddEntry(fSum_Es_122_ka,  "Escapes",    "l"); 
  legend->AddEntry(fSum_Es_122_ka,  "122-ka",    "l");
  legend->AddEntry(fSum_Es_122_ka,  "122-kb",    "l");
  legend->AddEntry(fSum_Es_136_ka,  "136-ka",    "l");
  legend->AddEntry(fSum_Es_136_kb,  "136-kb",    "l");
 
  
  
  legend->Draw();

  gPad->SetLogy(); 
}
