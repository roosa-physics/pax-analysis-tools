#include <iostream>
#include <fstream>

using namespace std;



// Define a struct that holds the parameter info necessary for plotting
struct par{
  TString name;
  TString units; 
  float rmin;
  float rmax;
  float binsize;  
};
//Make a list of defaults

// define some parameter structs
struct par pEnergy      = {"Energy",  "eV",
			   10e3, 250e3,   15.}; 
struct par pdT_Trig     = {"dT_trig", "ms",
			   1.0,   1.1, .0005}; 
struct par pRiseTime    = {"RiseTime", "Frac. Max",
			   0.55,   1.4, .0002}; 
struct par pPostPDeriv  = {"PostPDeriv", "Phu/Tu",
			   0.0,   6, .2}; 
struct par pfRes_avg    = {"fRes_avg", "Phu",
			   -50,   90, .5}; 
struct par pfRes_rms    = {"fRes_rms", "Phu",
			   1,   80, .1}; 
struct par pPTM         = {"PTrMean", "Phu",
			   4000,   6000, 1}; 
struct par pfRes_lnt    = {"fRes_lnt", "Phu",
			   //			    -200, 200, 1}; 
			   -7000, 10000, 5}; 
struct par pPH_tdc      = {"fValue_tdc",  "phu",
			   2.e3, 6e3, .7};




// Define a function that can take the parameter plotting info and make a 1D plot with lables
TH1F* Plot_1p(TTree* tree, struct par par ,TString title, TString cut_str="", bool overwrite = false){  
  // in case of a temporary hist creation loop with repeat names
  // it is  useful to find the old hist and clear it
  // but for now I will default to appending mode 
  float bins = (par.rmax-par.rmin)/par.binsize; 


  TH1F* htemp; 
  if(gDirectory->FindObject(title)){ // check if a hist of given title exists
    htemp = (TH1F*)gDirectory->Get(title);
    if(overwrite) htemp->Clear(); 
  }
  else //otherwise create the desired hist
    htemp = new TH1F(title,
		     Form("%s; %s (%s); Counts per %g %s",
			  cut_str.Data(),                     // title 
			  par.name.Data(), par.units.Data(),  // x label 
			  par.binsize,     par.units.Data()), // y label
		     bins, par.rmin, par.rmax); 
     
  int counts = tree->Draw(Form("%s>>+%s", par.name.Data(), title.Data()), cut_str.Data(), "");
  printf("Draw call %s %s : %d counts\n", title.Data(), cut_str.Data(), counts); 
  
  return htemp; 
}


// Define a function that can take the parameter plotting info and make a 2D plot with lables
TH2F* Plot_2p(TTree* tree, struct par par_x, struct par par_y, TString title, TString cut_str="", bool overwrite = false){  
  // in case of a temporary hist creation loop with repeat names
  // it is  useful to find the old hist and clear it
  // but for now I will default to appending mode 
  float bins_x = (par_x.rmax-par_x.rmin)/par_x.binsize; 
  float bins_y = (par_y.rmax-par_y.rmin)/par_y.binsize;

  TH2F* htemp; 
  if(gDirectory->FindObject(title)){ // check if a hist of given title exists
    htemp = (TH2F*)gDirectory->Get(title);
    if(overwrite) htemp->Clear(); 
  }
  else //otherwise create the desired hist
    htemp = new TH2F(title,
		     Form("%s; %s (%s); %s (%s)",
			  cut_str.Data(),                         // title 
			  par_x.name.Data(), par_x.units.Data(),  // x label 
			  par_y.name.Data(), par_y.units.Data()), // y label 
		     bins_x, par_x.rmin, par_x.rmax,
		     bins_y, par_y.rmin, par_y.rmax); 
     
  int counts = tree->Draw(Form("%s:%s>>+%s", par_y.name.Data(), par_x.name.Data(), title.Data()), cut_str.Data(), "");
  printf("Draw call %s %s : %d counts\n", title.Data(), cut_str.Data(), counts); 

  
  return htemp; 
}




TTree* Get_TTree(const char* fname){
  auto fFile = TFile::Open(fname);
  return (TTree*)fFile->Get("tTES");
}


void writeTGE(TGraphErrors* g, TString fname = "ascii.txt"){
  std::ofstream file;
  file.open(fname.Data());
  
  for(int i = 0; i < g->GetN(); i++){
    file << g->GetPointX(i) << ", "
	 << g->GetPointY(i) << ", "
	 << g->GetErrorX(i) << ", "
	 << g->GetErrorY(i) << endl;
  }
  file.close();
}


