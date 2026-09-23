ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "Minimize");

struct EscapeTF1_ka {
EscapeTF1_ka(TGraph2D* g) : gGraph(g) {}  
  double operator() (const double * x, const double *p) {
    float xx = 25269.7516781+(p[1]-x[0])*p[2];
    float sigma  = p[4]*p[2]; 
    return p[0]*gGraph->Interpolate(xx, sigma)+p[3]; 
  } 
  
  TGraph2D* gGraph;   
};

struct EscapeTF1_kb {
EscapeTF1_kb(TGraph2D* g) : gGraph(g) {}  
  double operator() (const double * x, const double *p) {
    float xx = 28482.3202783000+(p[1]-x[0])*p[2];
    float sigma  = p[4]*p[2]; 
    return p[0]*gGraph->Interpolate(xx,sigma)+p[3]; 
  } 
  
  TGraph2D* gGraph;   
};

double extrap(double x, double x1, double y1, double x2, double y2){
  double slope = (y2-y1)/(x2-x1);
  double dx = x-x1;
  return slope*dx+y1; 
}

void PrintPars(TF1* f){
  printf("Parameters: \n");
  for(int i=0; i<f->GetNpar()-1; i++)
    printf("\t\t %g", f->GetParameter(i)); 
  printf("\n"); 
}

void Notify(TFitResultPtr ptr){
  if(ptr->Status()){
    printf("Problem Fitting %s; Status: %d\n", ptr->GetName(), ptr->Status());
    ptr->Print(); 
  }
}

void clean_gr(TGraph* gr, double thr, double min){
  for(int i=gr->GetN(); i>0; i--)
    if(gr->GetPointX(i)<min) gr->RemovePoint(i);
  for(int i=0; i<gr->GetN(); i++){    
    for(int j=i+1; j<gr->GetN(); j++){
      double de = abs(gr->GetPointX(i)-gr->GetPointX(j)); 
      while(de<thr){
	//printf("cutting %g\n", gr->GetPointX(j) ); 
	gr->RemovePoint(j);
	de = abs(gr->GetPointX(i)-gr->GetPointX(j));	
      }
    }
  }	
}	



void calibrate_57Co(TH1F* htemp, TString sName = TString("Cal57Co"), bool vis = false){
  //pPH_tdc      = {"fValue_tdc",  "phu",
  //		  2e3, 6e3,   .5};
  int npts = 10e3; 


  //Gaussian lines
  auto fC122 = new TF1("fC122", "gausn", 0, 12e3); 
  fC122->SetNpx(npts); 
  
  auto fC136 = new TF1("fC136", "gausn", 0, 12e3); 
  fC136->SetNpx(npts); 	

  //ka escapes
  TFile* fKaMap = new TFile("DR_lookup_ka.root");
  TGraph2D* gKa = (TGraph2D*)fKaMap->Get("DR_lookup_ka");
    
  TF1* fEs_122_ka = new TF1("fEs_122_ka", EscapeTF1_ka(gKa),
			    pPH_tdc.rmin, pPH_tdc.rmax,5);  
  fEs_122_ka->SetLineColor(kOrange); 
  fEs_122_ka->SetNpx(npts);
  
  TF1* fEs_136_ka = new TF1("fEs_136_ka", EscapeTF1_ka(gKa),
			    pPH_tdc.rmin, pPH_tdc.rmax ,5);
  fEs_136_ka->SetLineColor(kPink+2); 
  fEs_136_ka->SetNpx(npts);
  
  //kb escapes
  TFile* fKbMap = new TFile("DR_lookup_kb.root");
  TGraph2D* gKb = (TGraph2D*)fKbMap->Get("DR_lookup_kb");
  
  TF1* fEs_122_kb = new TF1("fEs_122_kb", EscapeTF1_kb(gKb),
			    pPH_tdc.rmin, pPH_tdc.rmax ,5); 
  fEs_122_kb->SetLineColor(kOrange+1); 
  fEs_122_kb->SetNpx(npts);
  
  TF1* fEs_136_kb = new TF1("fEs_136_kb", EscapeTF1_kb(gKb),
			    pPH_tdc.rmin, pPH_tdc.rmax ,5);
  fEs_136_kb->SetLineColor(kPink+1); 
  fEs_136_kb->SetNpx(npts);
  
  //Energy Library
  double E_122 = 122060.65; //eV
  double E_136 = 136473.56; //eV
  double E_Ka = 25269.75;  //eV 
  double E_Kb = 28482.32;  //eV
  
  double E_122esKa = E_122-E_Ka;  
  double E_122esKb = E_122-E_Kb;  
  double E_136esKa = E_136-E_Ka;  
  double E_136esKb = E_136-E_Kb;  

  double E_[] = {E_122esKb, E_122esKa, E_136esKb, E_136esKa, E_122, E_136}; 
  
  //Fit parameters
  float guess_sigma = 40.; //eV
  float range = 105; 
  
  // fit ranges for the escapes
  double ka_hi = 1100.; 
  double ka_lo = -800.; 
  
  double kb_hi = 550.; 
  double kb_lo = -1450.; 

  // ---- Fit 122 ---- //
  printf("Fitting 122Co57\n");
    
  //Guess values
  double i0 = 10000; // intensity  
  int   max_bin = htemp->GetMaximumBin(); 
  float max_phv = htemp->GetBinCenter(max_bin);

  double g0    = max_phv/E_122; //calculate gain @ 122 ->  phv / eV
  double ig0   = E_122/max_phv; // inverse ev/phv
  double sigma = guess_sigma*g0; //convert eV to PhV
  
  float lo = max_phv-range*g0; 
  float hi = max_phv+range*g0;

  fC122->SetParameters(i0, max_phv, sigma);
    
  fC122->SetParLimits(1, lo, hi);              // mean
  fC122->SetParLimits(2, sigma/2., sigma*2);   // sigma
    
  //PrintPars(fC122); 
  auto  r122 = htemp->Fit("fC122","LRSQ","SAME", lo, hi); // Fit call <---
  //Notify(r122); //r122->Print(); 

  i0 = r122->Parameter(0);      //update intensity
  max_phv = r122->Parameter(1); //update max intensity est 
  sigma = r122->Parameter(2);   //update sigma est
  g0 = max_phv/E_122;            //update gain
  ig0 = E_122/max_phv;           //update inv gain

  // print nice stats for quick goodness checking
  printf("122Co57: SEM %g eV, ", r122->ParError(1)*ig0); 
  printf("Sigma %g eV\n",        r122->Parameter(2)*ig0); 

  printf("done!\n\n");
  // ---- Finish ---- //

    
  // ---- Estimate Background ---- // 
  //auto hbkg = htemp->ShowBackground(20, "BackDecreasingWindow");
  //hbkg->Draw("SAME"); 
  // ---- Finish ---- //

  // ---- Find Peaks ---- //
  // sigma, option, theshold
  float threshold  = .005; 
  htemp->ShowPeaks(r122->Parameter(2), "", threshold);
 
  auto fList  = htemp->GetListOfFunctions(); 
  auto mPeaks = (TPolyMarker*)fList->FindObject("TPolyMarker"); 
  Int_t nPeaks = mPeaks->GetN();

  Double_t* vPeaks = mPeaks->GetX();
  TGraph* gPhE_guess = new TGraph(nPeaks, vPeaks, vPeaks);
  printf("------------- strt %d peaks found with threshold = %g\n", gPhE_guess->GetN(), threshold);
  // for(int i = 0; i < gPhE_guess->GetN(); i++) printf("%g, ", gPhE_guess->GetPointX(i)); 
  printf("\n"); 

  
  float min = g0*93e3; //eV
  printf("Peak search in = %g\n", min); 
  // if 2 peaks are closer than 2keV drop the lower intensity one
  clean_gr(gPhE_guess, g0*3e3, min); 
  
  //Repeat the peak finding until there are 9 peaks spaced by more than 2keV 
  while(gPhE_guess->GetN()<6){
    threshold -= .00001;
    htemp->ShowPeaks(r122->Parameter(2), "", threshold);
    mPeaks = (TPolyMarker*)fList->FindObject("TPolyMarker"); 
    nPeaks = mPeaks->GetN();
    vPeaks = mPeaks->GetX();
    gPhE_guess = new TGraph(nPeaks, vPeaks, vPeaks);
    
    // if 2 peaks are closer than 2keV drop the lower intensity one
    clean_gr(gPhE_guess, g0*3e3, min);
  }
    
  
  printf("------------- stop %d peaks found with threshold = %g\n", gPhE_guess->GetN(), threshold);  
  gPhE_guess->Sort(); //Order By Energy        
  //  for(int i = 0; i < gPhE_guess->GetN(); i++) printf("%g, ", gPhE_guess->GetPointX(i)); 
    printf("\n"); 


  gPhE_guess->SetPoint(0, gPhE_guess->GetPointX(0), E_122esKb);
  gPhE_guess->SetPoint(1, gPhE_guess->GetPointX(1), E_122esKa);
  gPhE_guess->SetPoint(2, gPhE_guess->GetPointX(2), E_136esKb);
  gPhE_guess->SetPoint(3, gPhE_guess->GetPointX(3), E_136esKa);
  gPhE_guess->SetPoint(4, gPhE_guess->GetPointX(4), E_122);
  gPhE_guess->SetPoint(5, gPhE_guess->GetPointX(5), E_136);

  double slope = (E_136-E_122esKb)/(gPhE_guess->GetPointX(5)-gPhE_guess->GetPointX(0));
  double ph136 = gPhE_guess->GetPointX(5);   

  // check for too non-linear guesses
  // -> wrong peak of doublets
  for(int i=0; i<gPhE_guess->GetN()-2; i++){  
    double ph = gPhE_guess->GetPointX(i);
    double en = gPhE_guess->GetPointY(i);
    
    double dph = ph-ph136; 
    double guess = E_136 + slope * dph;
    double de = en-guess; 
    if(i==3 && de<-150)
      gPhE_guess->SetPointX(i, ph-de*g0);
    if(i==2 && de>500)
      gPhE_guess->SetPointX(i, ph+de*g0);
  
    // printf("%g %g\n", guess, en-guess);
  }


  
  // ---- Finish ---- //
    
  // ---- Fit 136 ---- //
  printf("Fitting 136..."); 
  double v136 = gPhE_guess->GetPointX(5);      
  lo = v136-range*g0; 
  hi = v136+range*g0;

  // Set Fit pars
  fC136->SetParameters(i0, v136, sigma);
   
  fC136->SetParLimits(1, lo, hi);              // mean  
  fC136->SetParLimits(2, sigma/2., sigma*2);   // sigma
    
  //PrintPars(fC136); 
  auto r136 = htemp->Fit("fC136","LRSQ","SAME", lo, hi); // Fit call <--- 
  Notify(r136); 

  // print nice stats for quick goodness checking
  //printf("136Co57: SEM %g eV, ", r136->ParError(1)*ig0); 
  //printf("Sigma %g eV\n",        r136->Parameter(2)*ig0);  
  printf("done!\n"); 
  // ---- Finish ---- //
    
  // ---- Fit escapes ---- //
  printf("Fitting 122 escape ka..."); 
  double v122ea = gPhE_guess->GetPointX(1);  
  lo = v122ea + ka_lo * g0; 
  hi = v122ea + ka_hi * g0; 
    
  fEs_122_ka->SetRange(lo, hi);
  fEs_122_ka->SetNpx(npts);
    
  fEs_122_ka->SetParameters(.6*i0, v122ea, ig0, sigma);
  fEs_122_ka->SetParLimits(0,.2*i0, 1.5*i0); 
  fEs_122_ka->SetParLimits(1, v122ea-25, v122ea+25);
  fEs_122_ka->SetParLimits(2, .5*ig0, 2*ig0);   
  fEs_122_ka->SetParLimits(3, 1, 50);
    
  fEs_122_ka->ReleaseParameter(2);
  fEs_122_ka->ReleaseParameter(4); 
  fEs_122_ka->FixParameter(2, ig0);
  fEs_122_ka->FixParameter(4, sigma);
    
  //PrintPars(fEs_122_ka);  
  auto r122eka = htemp->Fit("fEs_122_ka","LRSQ","GOFF", lo, hi); 
  if(r122eka->Status()>2) //Try 1 more time
    r122eka = htemp->Fit("fEs_122_ka","LRSQ","GOFF", lo, hi); 
    
  Notify(r122eka);  
  //printf("\n SEM %s %g eV\n", fEs_122_ka->GetName(), E_122esKa*r122eka->ParError(1)/r122eka->Parameter(1));  
  printf("done!\n"); 
    
  printf("Fitting 122 escape kb..."); 
  double v122eb=gPhE_guess->GetPointX(0); 
  lo = v122eb + kb_lo * g0;  
  hi = v122eb + kb_hi * g0; 
 
  fEs_122_kb->SetParameters(.6*i0, v122eb, ig0, sigma);
  fEs_122_kb->SetParLimits(0,.1*i0, 2*i0); 
  fEs_122_kb->SetParLimits(1, v122eb-25, v122eb+25);
  fEs_122_kb->SetParLimits(2, .5*ig0, 2*ig0);   
  fEs_122_kb->SetParLimits(3, 1, 50);

  fEs_122_kb->ReleaseParameter(2); 
  fEs_122_kb->ReleaseParameter(4);   
  fEs_122_kb->FixParameter(2, ig0);
  fEs_122_kb->FixParameter(4, sigma);

  //PrintPars(fEs_122_kb);
  auto r122ekb = htemp->Fit("fEs_122_kb","LRSQ","GOFF", lo, hi);  // Fit call <---
  if(r122ekb->Status()>2) //Try 1 more time
    r122ekb = htemp->Fit("fEs_122_kb","LRSQ","GOFF", lo, hi); 
  
  fEs_122_kb->SetRange(lo, hi); // for visualization
   
  //Notify(r122ekb); 
  //printf("SEM %s %g eV\n", fEs_122_kb->GetName(),
  //E_122esKb*r122ekb->ParError(1)/r122ekb->Parameter(1)); 
  printf("done!\n"); 
 
      
  printf("Fitting 136 escape kb..."); 
  double v136eb = gPhE_guess->GetPointX(2); 
  lo = v136eb + kb_lo * g0;  
  hi = v136eb + kb_hi * g0;
    
  fEs_136_kb->SetParameters(.6*r136->Parameter(0), v136eb, ig0, sigma);
  fEs_136_kb->SetParLimits(0,.1*r136->Parameter(0), 2*r136->Parameter(0)); 
  fEs_136_kb->SetParLimits(1, v136eb-45, v136eb+25);
  fEs_136_kb->SetParLimits(2, .5*ig0, 2*ig0);    
  fEs_136_kb->SetParLimits(3, 1, 50);

  fEs_136_kb->ReleaseParameter(2); 
  fEs_136_kb->ReleaseParameter(4);   
  fEs_136_kb->FixParameter(2, ig0);
  fEs_136_kb->FixParameter(4, sigma);

  //PrintPars(fEs_136_kb);
  auto r136ekb = htemp->Fit("fEs_136_kb","LRSQ","GOFF", lo, hi);  // Fit call <---
  if(r136ekb->Status()>2) //Try 1 more time
    r136ekb = htemp->Fit("fEs_136_kb","LRSQ","GOFF", lo, hi);
  //r136ekb->Print(); 

  fEs_136_kb->SetRange(lo,hi); // for visualization
  Notify(r136ekb); 
  //printf("SEM %s %g eV\n", fEs_136_kb->GetName(),
  //E_136esKb*r136ekb->ParError(1)/r136ekb->Parameter(1));   
  printf("done!\n");
  
  
  printf("Fitting 136 escape ka..."); 
  double v136ea = -5 + gPhE_guess->GetPointX(3);
  lo = v136ea + ka_lo * g0;  
  hi = v136ea + ka_hi * g0; 

  fEs_136_ka->SetParameters(.5*r136->Parameter(0), v136ea, ig0, sigma);
  fEs_136_ka->SetParLimits(0,.1*r136->Parameter(0), 2*r136->Parameter(0)); 
  fEs_136_ka->SetParLimits(1, v136ea-25, v136ea+25);
  fEs_136_ka->SetParLimits(2, .5*ig0, 2*ig0);    
  fEs_136_ka->SetParLimits(3, 1, 50);

  fEs_136_ka->ReleaseParameter(2); 
  fEs_136_ka->ReleaseParameter(4);   
  fEs_136_ka->FixParameter(2, ig0);
  fEs_136_ka->FixParameter(4, sigma);     
    
  //PrintPars(fEs_136_ka);   
  auto r136eka = htemp->Fit("fEs_136_ka","LRSQ","GOFF", lo, hi); // Fit call <---
  if(r136eka->Status()>2) //Try 1 more time
    r136eka = htemp->Fit("fEs_136_ka","LRSQ","GOFF", lo, hi);
 
  fEs_136_ka->SetRange(lo,hi); // for visualization

  Notify(r136eka); 
  // printf("SEM %s %g eV\n", fEs_136_ka->GetName(),
  // E_136esKa*r136eka->ParError(1)/r136eka->Parameter(1)); 
  printf("done!\n");
  // ---- Finish ---- //    
  printf("Finished fitting !! \n");   


  TGraphErrors* gResult = new TGraphErrors(5);
  gResult->SetTitle(Form("gCal_%s; O.F. Values; Energy (eV)", sName.Data()));
  gResult->SetLineColor(kBlue+1); 
    
  double err122 = .12; 
  double err136 = .29; 
  double ph, dph; 
    
  //gammas
  gResult->SetPoint(4, r122->Parameter(1), E_122);
  gResult->SetPointError(4, r122->Error(1), err122);
  
  gResult->SetPoint(5, r136->Parameter(1), E_136);
  gResult->SetPointError(5, r136->Error(1), err136);

  //122 ka escape
  ph  = r122eka->Parameter(1); gResult->SetPoint(      1,  ph, E_122esKa);
  dph = r122eka->Error(1);     gResult->SetPointError( 1, dph, err122);

  //122 kb escape
  ph  = r122ekb->Parameter(1); gResult->SetPoint(      0,  ph, E_122esKb);
  dph = r122ekb->Error(1);     gResult->SetPointError( 0, dph, err122);
    
  //136 ka escape
  ph  = r136eka->Parameter(1); gResult->SetPoint(      3,  ph, E_136esKa ); 
  dph = r136eka->Error(1);     gResult->SetPointError( 3, dph, err136);

  //136 kb escape
  ph  = r136ekb->Parameter(1); gResult->SetPoint(      2,  ph, E_136esKb);
  dph = r136ekb->Error(1);     gResult->SetPointError( 2, dph, err136);
    
    
  TFile* fout = new TFile(Form("%s.root",sName.Data()),"RECREATE"); 
  fout->cd();
  gResult->SetName(Form("%s_cal",sName.Data()));
  gResult->Write();
  r122->SetName(Form("%s_fr122",sName.Data())); 
  r122->Write(); 
  // ---- Finish ---- //

  bool visualize = true; 
  if(visualize){
    auto b1 = new TBrowser("b1", "b1", 1600, 900); 
    auto c1 = new TCanvas(Form("%s_c",sName.Data()),Form("c_%s",sName.Data()), 2400, 750); 
    c1->Divide(2); 

    // Energy spectrum with fit outputs in panel 1
    c1->cd(1);   
    htemp->Draw(); 
    htemp->GetXaxis()->SetRangeUser(90e3*g0, 140e3*g0);    
    gPad->SetLogy();
    
    fC122->Draw("SAME"); 
    fC136->Draw("SAME");
    fEs_122_ka->Draw("SAME");  
    fEs_122_kb->Draw("SAME");
    fEs_136_ka->Draw("SAME"); 
    fEs_136_kb->Draw("SAME");
    

    // Guess and results minus the slope in panel 2 
    c1->cd(2);
    auto l_c2 = new TLegend(.35,0.70,0.5,0.9);
   
    TGraph* gSlope = new TGraph(2);
    gSlope->SetPoint(0,
		     gResult->GetPointX(0),
		     gResult->GetPointY(0)); 
    gSlope->SetPoint(1,
		     gResult->GetPointX(gResult->GetN()-1),
		     gResult->GetPointY(gResult->GetN()-1)); 
      
    auto gResult_sub = (TGraphErrors*) gResult->Clone("gResult_sub");
    gResult_sub->SetMarkerStyle(20); 
    for(int i=0; i < gResult_sub->GetN(); i++){
      double xx = gResult->GetPointX(i);
      double yy = gResult->GetPointY(i)-gSlope->Eval(xx);
      gResult_sub->SetPoint(i,xx, yy); 
    }

    auto gPhE_guess_sub = (TGraph*) gPhE_guess->Clone("gPhE_guess_sub");
    gPhE_guess_sub->SetMarkerStyle(20);
    gPhE_guess_sub->SetLineColor(kRed); 
    gPhE_guess_sub->SetTitle("Fit outputs and initial Guesses - Slope"); 
    for(int i=0; i < gPhE_guess_sub->GetN(); i++){
      double xx = gPhE_guess->GetPointX(i);
      double yy = gPhE_guess->GetPointY(i)-gSlope->Eval(xx);
      gPhE_guess_sub->SetPoint(i,xx, yy); 
    }
     
    gResult_sub->Draw("AP"); l_c2->AddEntry(gResult_sub,    "Fit Results", "P");
    gResult_sub->SetTitle("Fit Results and initial guess - Slope"); 
    double x0, x1, x2, x3, y0, y1, y2, y3;
    gResult_sub->ComputeRange(x0, y0, x1, y1);
    gPhE_guess_sub->ComputeRange(x2, y2, x3, y3);    
    double min = (y0<y2)?y0:y2; 
    double max = (y1>y3)?y1:y3;
    gResult_sub->GetYaxis()->SetRangeUser(min-2,max+2); 

    gPhE_guess_sub->Draw("SAME");
    l_c2->AddEntry(gPhE_guess_sub, "Initial Guess", "L");

    l_c2->Draw(); 
    
    c1->cd(1); 
    fout->cd(); 
    c1->Write(); 
    //c1->Close(); 
    printf("\n\n"); 
  }
  fout->Close(); 
}
