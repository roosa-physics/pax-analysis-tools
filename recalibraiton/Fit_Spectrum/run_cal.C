#include"Utils.C"
#include"Graphical_Cut_List.C"
#include"calibrate_57Co.C"
#include"filelist.C"


void run_cal(){
  //load the data
  //for(int i = 0; i < sFileList.size(); i++){

  int i = 12; 
  printf("Opening %s\n", sFileList[i]->Data()); 
  TTree* tTES = Get_TTree(sFileList[i]->Data()); 

  auto strip_name = new TString(sFileList[i]->Data()); ;
  strip_name->ReplaceAll("coadd_0708/", ""); 
  strip_name->ReplaceAll("_0708.root", "");
    
  vector< par > pars = {pPH_tdc}; 
  int npars = pars.size(); 
    
  int no_chs = 32; 
  
  for(int ch = 0; ch<no_chs; ch++){

    printf("Make histos and calibrate %s, ch%d\n", strip_name->Data(), ch); 

    //declare and fill the per chan energy histograms  
    TString good_string = Form("isGood==1");  
    good_string += Form("&&Energy>70e3");
    good_string += Form("&&RiseTime>.99");
    good_string += Form("&&dT_trig>50");
    good_string += Form("&&Chan==%d", ch);  
    
    TH1F* htemp = Plot_1p(tTES, pars.at(0),
			  Form("h%s_ch%d", pars.at(0).name.Data(), ch),
			  good_string, true);
    
    if(htemp->GetEntries()==0){printf("0 entries, skipping\n"); continue;}
      
    calibrate_57Co(htemp, Form("escal_%s_ch%d",strip_name->Data(), ch), true); 
  }
 
}
