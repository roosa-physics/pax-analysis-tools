//void gGausResp(TGraph* g_old, TGraph* g_new, double sigma){  
//    g_new->SetPoint(i, x0, result); 
//  }   
//}


TGraph* gKa = new TGraph("xray_templates.txt", "%lg %lg %*lg %*lg");  
TGraph* gKb = new TGraph("xray_templates.txt", "%*lg %*lg %lg %lg");  


void DR_lookup(){ 
  int  n_sigs = 26; 

  int nKa = gKa->GetN(); 
  TGraph2D* g2Ka = new TGraph2D(nKa*n_sigs);   
  g2Ka->SetName("DR_lookup_ka"); 
  
  int nKb = gKb->GetN();
  TGraph2D* g2Kb = new TGraph2D(nKb*n_sigs);   
  g2Kb->SetName("DR_lookup_kb");
  
  for(int j = 0; j<n_sigs; j++){
    float sigma = 18+3*j;  

    double dx = 1.5;
    double lo = -6*sigma;
    double hi =  6*sigma;
    
    for(int i=0; i<nKa; i++){
      double result = 0;
      double x0 = gKa->GetPointX(i); 
      for(double xx = lo; xx<hi; xx+=dx) 
	result += gKa->Eval(x0+xx)*TMath::Gaus(xx,0,sigma,kTRUE); 
      g2Ka->SetPoint(i+nKa*j, x0, sigma, result); 
    }

    
    for(int i=0; i<nKb; i++){
      double result = 0;
      double x0 = gKb->GetPointX(i); 
      for(double xx = lo; xx<hi; xx+=dx) 
	result += gKb->Eval(x0+xx)*TMath::Gaus(xx,0,sigma,kTRUE); 
      g2Kb->SetPoint(i+nKb*j, x0, sigma, result); 
    }
  }


  g2Ka->SaveAs("DR_lookup_ka.root"); 
  g2Kb->SaveAs("DR_lookup_kb.root"); 

  TCanvas* c1 = new TCanvas();
  c1->Divide(2);
  c1->cd(1); g2Ka->Draw("COL"); 
  c1->cd(2); g2Kb->Draw("COL");

  c1->SaveAs("DRPlot.png"); 
  c1->SaveAs("DRPlot.root");
}
