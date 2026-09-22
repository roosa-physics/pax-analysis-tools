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

void maketemplates(){
  std::vector<TF1*> fAllPeaks_;
  
  TTree *t = new TTree("t", "snkas");
  t->ReadFile("sn_xrays.csv", "energy/D:intensity/D:width/D");
  Double_t energy, intensity, width; 
  t->SetBranchAddress("energy",    &energy); 
  t->SetBranchAddress("intensity", &intensity);
  t->SetBranchAddress("width",     &width);

  // set characteristic intensity
  t->GetEntry(0); double intensity0 = intensity; 

  std::vector<TF1*> fKa_;
  std::vector<TF1*> fKb_;
  
  int npts = 50000; 
  float e_lo = 20000;
  float e_hi = 40000;
  
  for(int i=0; i<t->GetEntries();i++){
    t->GetEntry(i); 
    if(energy>26000){
      fKb_.push_back(new TF1(Form("fL%d",i),"[0]*TMath::CauchyDist(x,[1],[2])",e_lo,e_hi));
      fKb_.back()->SetParameters(intensity/intensity0, energy, width);
      fKb_.back()->SetNpx(npts); 	
      //fKb_.back()->SetLineColor(kYellow+i);
    }else if(energy>24000){
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
  
  TF1* fSum_ka = new TF1("fSum_ka",SumTF1(fKa_),e_lo,e_hi,0);
  fSum_ka->SetTitle("Tin Calculated Spectrum; Energy [eV]; Intensity [Arb]");
  fSum_ka->SetLineColor(kMagenta); 
  fSum_ka->SetNpx(npts); 	
  
  double de = .1;
  double lo = 23500;
  double hi = 26500;
  double e0 = 4000; 
  for(double e = lo; e<hi; e+=de)
    printf("%g\t%g\t%g\t%g\n",
	   e,    fSum_ka->Eval(e),
	   e+e0, fSum_kb->Eval(e+e0)); 
    
 

  
}
