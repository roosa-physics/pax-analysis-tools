struct LineTF1 { 
  
LineTF1(TGraph* g) : gGraph(g) {}  
  double operator() (const double * x, const double *p) {
    double xx = p[2]*(x[0]-p[1]); 
   
    //if out of range return 0
    double min = gGraph->GetPointX(0);
    double max = gGraph->GetPointX(gGraph->GetN()-1); 
    if(xx<min || xx>max) return 0; // assumes x is sorted

    return p[0]*gGraph->Eval(xx); 
  } 
  
  TGraph* gGraph;   
};


struct EscapeTF1 { 
  
EscapeTF1(TGraph* g) : gGraph(g) {}  
  double operator() (const double * x, const double *p) {
    double xx = p[2]*(p[1]-x[0]); 
       
    //if out of range return 0
    double min = gGraph->GetPointX(0);
    double max = gGraph->GetPointX(gGraph->GetN()-1); 
    if(xx<min || xx>max) return 0; // assumes x is sorted

    return p[0]*gGraph->Eval(xx); 
  } 
  
  TGraph* gGraph;   
};



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

struct GRespTF1 { //
GRespTF1(TF1* fpeak) : fPeak(fpeak) {}
  double dx = 1;
  
  double operator() (const double * x, const double *p) {
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

void plotlines_template(){

  std::vector<TF1*> fAllPeaks_;  

  // set characteristic intensity
  double intensity0 = 1567.1100000000; 
  double norm = 1./intensity0;  

  double e_lo = 90e3; 
  double e_hi = 140e3;
  int npts = 10000; 
  
  TGraph* gKa = new TGraph("xray_templates.txt", "%lg %lg %*lg %*lg");  
  TF1* fKa = new TF1("fKa",LineTF1(gKa),e_lo,e_hi,3);
  fKa->SetParameters(norm, 0., 1.); 
  fKa->SetTitle("Tin-Ka Calculated Spectrum; Energy [eV]; Intensity [Arb]"); 
  fKa->SetLineColor(kRed); 
  fKa->SetNpx(npts); 	
  // fAllPeaks_.push_back(fKa);

  
  TGraph* gKb = new TGraph("xray_templates.txt", "%*lg %*lg %lg %lg");  
  double min_kb, max_kb;
  TF1* fKb = new TF1("fKb",LineTF1(gKb),e_lo, e_hi,3);
  fKb->SetParameters(norm, 0., 1.); 
  fKb->SetTitle("Tin-Kb Calculated Spectrum; Energy [eV]; Intensity [Arb]"); 
  fKb->SetLineColor(kViolet); 
  fKb->SetNpx(npts); 	
  // fAllPeaks_.push_back(fKb);
  
  
  
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
  
  double iEscape = I122*.001; 
  double lo = 23500;
  double hi = 26500;
  double e0 = 4000; 
  
  

  
  TF1* fEs_122_ka = new TF1("fEs_122_ka", EscapeTF1(gKa),
			    E122-hi, E122-lo, 3);
  fEs_122_ka->SetParameters(iEscape, E122, 1);
  fEs_122_ka->SetLineColor(kOrange); 
  fEs_122_ka->SetNpx(npts);
  fAllPeaks_.push_back(fEs_122_ka);

  TF1* fEs_122_kb = new TF1("fEs_122_kb", EscapeTF1(gKb),
			    E122-hi-e0, E122-lo-e0, 3);
  fEs_122_kb->SetParameters(iEscape, E122, 1);
  fEs_122_kb->SetLineColor(kOrange+1); 
  fEs_122_kb->SetNpx(npts); 	
  fAllPeaks_.push_back(fEs_122_kb);

  TF1* fEs_136_ka = new TF1("fEs_136_ka", EscapeTF1(gKa),
			    E136-hi, E136-lo, 3);
  fEs_136_ka->SetParameters(iEscape*I136/I122, E136, 1);
  fEs_136_ka->SetLineColor(kOrange+2); 
  fEs_136_ka->SetNpx(npts);
  fAllPeaks_.push_back(fEs_136_ka);

  TF1* fEs_136_kb = new TF1("fEs_136_kb", EscapeTF1(gKb),
			    E136-hi-e0, E136-lo-e0, 3);
  fEs_136_kb->SetParameters(iEscape*I136/I122, E136, 1);
  fEs_136_kb->SetLineColor(kOrange+3); 
  fEs_136_kb->SetNpx(npts); 	
  fAllPeaks_.push_back(fEs_136_kb);
 
  //------------- Make Sum Spectrum -------------//    
  TF1* fSum = new TF1("fSum",SumTF1(fAllPeaks_),e_lo,e_hi,0);
  fSum->SetTitle("Sum Spectrum; Energy [eV]; Intensity [Arb]");
  fSum->SetLineColorAlpha(kCyan, .5); 
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
  fKa->Draw("SAME");
  fKb->Draw("SAME");

  fEs_122_ka->Draw("SAME"); 
  fEs_122_kb->Draw("SAME"); 
  fEs_136_ka->Draw("SAME"); 
  fEs_136_kb->Draw("SAME"); 

  
  auto legend = new TLegend(.1,0.70,0.2,0.9);
  legend->AddEntry(fSum,    "Sum",    "l"); 
  legend->AddEntry(fSim,    "Sum x Det. Resp.",    "l"); 
  legend->AddEntry(fCal_[0],    "Calibration",    "l"); 
  //legend->AddEntry(fKa,     "Sn-Ka",    "l"); 
  //legend->AddEntry(fKb,     "Sn-Kb",    "l"); 
  //  legend->AddEntry(fEs_122_ka,  "Escapes",    "l"); 
  legend->AddEntry(fEs_122_ka,  "122-ka",    "l");
  legend->AddEntry(fEs_122_ka,  "122-kb",    "l");
  legend->AddEntry(fEs_136_ka,  "136-ka",    "l");
  legend->AddEntry(fEs_136_kb,  "136-kb",    "l");
 
  
  
  legend->Draw();
  
  
  gPad->SetLogy(); 
}
