// plots trigger rates vs. fill
// use trg to control trg=fmsOnly setup
//  - 2 to use all trg setups
//  - 1 pp200_production_fms_2012
//  - 0 pp200_production_2012
// use rCutType to cut by run number
//  - <0 for runs before and equal to rCut
//  - =0 for all runs (ignore rCut)
//  - >0 for runs after rCut

void trgMon(Int_t beam_en=200, Int_t trg=2, Int_t rCutType=0,Int_t rCut=0)
{
  gStyle->SetOptStat(0);


  /* trigger list */
  ///////////////////////////////////////////////////
  const Int_t N_TRIG = 3;
  enum trigger_enum {kLED,kZDC,kVPD};
  const Int_t N_MB = 2;
  enum minbias_enum {kZDCMB,kVPDMB};


  char trigger_str[N_TRIG][16];
  char trigger_ps_str[N_TRIG][16];
  strcpy(trigger_str[kLED],"LED");
  strcpy(trigger_str[kZDC],"ZDC");
  strcpy(trigger_str[kVPD],"VPD");
  for(Int_t i=0; i<N_TRIG; i++) sprintf(trigger_ps_str[i],"%sps",trigger_str[i]);

  char minbias_str[N_MB][16];
  char minbias_ps_str[N_MB][16];
  strcpy(minbias_str[kZDCMB],"ZDC");
  strcpy(minbias_str[kVPDMB],"VPD");
  for(Int_t i=0; i<N_MB; i++) sprintf(minbias_ps_str[i],"%sps",minbias_str[i]);

  Color_t trigger_col[N_TRIG];
  trigger_col[kLED] = kRed;
  trigger_col[kZDC] = kAzure;
  trigger_col[kVPD] = kMagenta+2;

  ///////////////////////////////////////////////////

  // SELECTION OF MB TRIGGER
  // using VPD rather than BBC since for fms-only runs, no BBC triggers
  const Int_t MB_SELECT = kVPDMB;
  


  /* exclusion list */
  ///////////////////////////////////////////////////
  const Int_t N_EXCLUDE = 1;
  Int_t exclusion_list[N_EXCLUDE] = 
  {
    17000000,
  };
  /*
   * 56047 -- outlier on most triggers
   * 62008 -- outlier on most triggers
   * 64076 -- JP2 prescale = 140 (should be 1)
   * 67044 -- outlier on most triggers
   */
  ///////////////////////////////////////////////////
  
  
  
  /* define time and outlier cut */
  ///////////////////////////////////////////////////
  Int_t time_cut=180; // [seconds]
  ///////////////////////////////////////////////////






  /* read trigger table */
  TFile * fout = new TFile("runlog.root","RECREATE");
  char trigger_tree_str[2048];
  strcpy(trigger_tree_str,"index/I:runnum/I:day/I:run/I:type/C:fill/F:energy/F:time/I");
  for(Int_t i=0; i<N_TRIG; i++) 
    sprintf(trigger_tree_str,"%s:%s/D:%s/D",trigger_tree_str,trigger_str[i],trigger_ps_str[i]);
  // !!!
  // do not scan for minbias branches, since minbias branches included in regular triggers branches
  //for(Int_t i=0; i<N_MB; i++)
    //sprintf(trigger_tree_str,"%s:%s/D:%s/D",trigger_tree_str,minbias_str[i],minbias_ps_str[i]);
    //
  printf("%s\n",trigger_tree_str);
  TTree * tr = new TTree("tr","tr");
  tr->ReadFile("trigger.dat",trigger_tree_str);
  tr->Print();


  /* obtain run index numbers */
  Int_t nRuns = tr->GetEntries();
  Int_t index;
  tr->SetBranchAddress("index",&index);
  tr->GetEntry(0);
  Int_t iterl = index;
  tr->GetEntry(nRuns-1);
  Int_t iterh = index;
  Int_t niter = iterh-iterl;
  printf("plotting run index %d-%d, %d total\n",iterl,iterh,niter);
  const Int_t N_RUNS = nRuns;

  /* initialize fail code (counting left to right)
   * bit 1 -- on exclusion list or bad epoch or before first run
   * bit 2 -- mb = 0 
   * bit 3 -- number of events < 50 in any trigger (except LED)
   * bit 4 -- run too short
   * bit 5 -- trigger 0 outside of range
   * bit 6 -- trigger 1 outside of range
   *    .
   *    .
   *    .
   * bit 15 -- trigger 10 outside of range
   * bit 16 -- trigger 11 = LED; not QAing with LED, so this bit is unused
   */
  Long_t failcode[N_RUNS];
  for(Int_t r=0; r<N_RUNS; r++) failcode[r]=0;

  /* initialise bad run array */
  Bool_t badrun[N_RUNS];
  for(Int_t r=0; r<N_RUNS; r++)
  {
    if(r<iterl || r>iterh) 
    {
      badrun[r]=true;
      failcode[r] = failcode[r] | 0x8000;
    }
    else badrun[r]=false;
  };


  /* epoch lines */
  ///////////////////////////////////////////////////
  const Int_t N_EPOCH=2;
  Int_t epoch[N_EPOCH+1][N_TRIG];
  for(Int_t i=0; i<N_TRIG; i++)
  {
    // "epoch n" <--> epoch[n] <= index < epoch[n+1]
    epoch[0][i] = iterl; // first epoch bad since triggers not at physics
    epoch[1][i] = 1;
    epoch[2][i] = iterh;
  };

  /* bad epochs */
  const Int_t N_BAD_EPOCH=1;
  Int_t bad_epoch[N_BAD_EPOCH] = {0};
  /*
   * epoch 0 -- before triggers elevated to physics
   * epoch 2 -- several hot towers causing high rates
   * epoch 10 -- several hot towers causing high rates
   */
 
  /* bold epochs */ // for drawing certain epoch lines in red rather than blue
  const Int_t N_BOLD_EPOCH=1;
  Int_t bold_epoch[N_BOLD_EPOCH] = {0};
  ///////////////////////////////////////////////////

  /* epoch to fit with exponential decay */
  const Int_t N_EXP_FITS=1;
  Int_t exp_fit_list[N_EXP_FITS] = {1};



  /* run QA lines */
  ///////////////////////////////////////////////////
  Double_t qa_cut_high[N_TRIG];
  Double_t qa_cut_low[N_TRIG];
  const Double_t QA_PLOT_MAX = 10;
  const Double_t QA_PLOT_MIN = 0.1;
  ////
  qa_cut_high[kLED] = QA_PLOT_MAX;
  qa_cut_high[kZDC] = QA_PLOT_MAX;
  qa_cut_high[kVPD] = QA_PLOT_MAX;
  ////
  qa_cut_low[kLED] = QA_PLOT_MIN;
  qa_cut_low[kZDC] = QA_PLOT_MIN;
  qa_cut_low[kVPD] = QA_PLOT_MIN;
  ///////////////////////////////////////////////////


  /* print exclusion list */
  char exclusion_list_print[400];
  strcpy(exclusion_list_print,"");
  for(Int_t it=0; it<N_EXCLUDE; it++)
  {
    if(it+1 == N_EXCLUDE) sprintf(exclusion_list_print,"%srunnum==%d",exclusion_list_print,exclusion_list[it]);
    else sprintf(exclusion_list_print,"%srunnum==%d||",exclusion_list_print,exclusion_list[it]);
  };
  printf("%s\n",exclusion_list_print);
  printf("\nOUTLIER LIST:\n");
  tr->SetScanField(0);
  tr->Scan("index:runnum",exclusion_list_print);
  system("touch outliers.dat; rm outliers.dat"); 
  gSystem->RedirectOutput("outliers.dat");
  tr->Scan("index:runnum",exclusion_list_print);
  gSystem->RedirectOutput(0);
  printf("these runs are manually removed from the plots\n\n");



  /* set trigger tree branch addresses */
  Double_t nev[N_TRIG];
  Double_t nev_ps[N_TRIG];
  Double_t mb[N_MB];
  Double_t mb_ps[N_MB];
  Int_t runnum,day,run,time;
  Float_t energy,fill;
  char type[64];
  tr->SetBranchAddress("runnum",&runnum);
  tr->SetBranchAddress("day",&day);
  tr->SetBranchAddress("run",&run);
  tr->SetBranchAddress("fill",&fill);
  tr->SetBranchAddress("energy",&energy);
  tr->SetBranchAddress("time",&time);
  tr->SetBranchAddress("type",type);
  for(Int_t i=0; i<N_TRIG; i++)
  {
    tr->SetBranchAddress(trigger_str[i],&(nev[i]));
    tr->SetBranchAddress(trigger_ps_str[i],&(nev_ps[i]));
  };
  /*
  for(Int_t i=0; i<N_MB; i++)
  {
    tr->SetBranchAddress(minbias_str[i],&(mb[i]));
    tr->SetBranchAddress(minbias_ps_str[i],&(mb_ps[i]));
  };
  */


  /* initialise plots */
  TGraph * ev_gr[N_EPOCH][N_TRIG]; // (nev*nev_ps) / (mb*mb_ps)
  TGraph * nv_gr[N_EPOCH][N_TRIG]; // (nev*nev_ps) / (mb*mb_ps*epochMean) [normalised ev_gr]
  for(Int_t i=0; i<N_TRIG; i++)
  {
    for(Int_t e=0; e<N_EPOCH; e++)
    {
      ev_gr[e][i] = new TGraph();
      nv_gr[e][i] = new TGraph();
    };
  };


  /* build ev_gr plots */
  Bool_t PLOT;
  Double_t value;
  Int_t current_epoch[N_TRIG];
  Int_t gr_i[N_EPOCH][N_TRIG];
  Double_t ev_mgr_max[N_TRIG];
  Int_t ev_mgr_max_index[N_TRIG];
  for(Int_t i=0; i<N_TRIG; i++)
  {
    ev_mgr_max[i]=0;
    for(Int_t e=0; e<N_EPOCH; e++) gr_i[e][i]=0;
  };
  for(Int_t r=0; r<tr->GetEntries(); r++)
  {
    tr->GetEntry(r);
    mb[kZDCMB] = nev[kZDC]; // MB MAP !!!!!!!!!!!!!!!!!!!!!!!+++++++++++++++
    mb[kVPDMB] = nev[kVPD];
    mb_ps[kVPDMB] = nev_ps[kVPD];
    mb_ps[kVPDMB] = nev_ps[kVPD];
    PLOT=true;

    // check if run duration is greater than time_cut
    if(time < time_cut) 
    {
      PLOT=false;
      failcode[index] = failcode[index] | 0x1000;
    };

    // check if run is in exclusion list
    for(Int_t o=0; o<N_EXCLUDE; o++)
    {
      if(runnum == exclusion_list[o]) 
      {
        PLOT=false;
        failcode[index] = failcode[index] | 0x8000;
      };
    };

    // make sure we won't divide by zero
    if(mb[MB_SELECT]*mb_ps[MB_SELECT]<1) 
    {
      PLOT=false;
      failcode[index] = failcode[index] | 0x4000;
    };

    // check to see if we have enough events in each
    // trigger, except for LED
    /*
    for(Int_t i=0; i<N_TRIG; i++)
    {
      if(i!=kLED && nev[i]<50)
      {
        PLOT=false;
        failcode[index] = failcode[index] | 0x2000;
      };
    };
    */


    if(PLOT)
    {
      // determine which epoch we are in
      for(Int_t i=0; i<N_TRIG; i++)
      {
        current_epoch[i] = -1;
        for(Int_t e=0; e<N_EPOCH; e++)
        {
          if(index>=epoch[e][i] && index<=epoch[e+1][i] && epoch[e][i]!=epoch[e+1][i])
            current_epoch[i] = e;
        };
      };

      // check if it's a bad epoch and mark its runs bad;
      // if it's not really really off scale, plot it anyway and turn the point colors to red
      for(Int_t i=0; i<N_TRIG; i++)
      {
        for(Int_t e=0; e<N_EPOCH; e++)
        {
          for(Int_t b=0; b<N_BAD_EPOCH; b++)
          {
            if(bad_epoch[b] == current_epoch[i]) 
            {
              badrun[index]=true;
              failcode[index] = failcode[index] | 0x8000;
            };
          };
        };
      };


      // add point to ev_gr
      for(Int_t i=0; i<N_TRIG; i++)
      {
        if(current_epoch[i]>=0)
        {
          value = (nev[i]*nev_ps[i])/(mb[MB_SELECT]*mb_ps[MB_SELECT]);
          ev_gr[current_epoch[i]][i]->SetPoint(gr_i[current_epoch[i]][i],index,value);
          (gr_i[current_epoch[i]][i])++;
          if(value > ev_mgr_max[i])
          {
            ev_mgr_max[i] = value;
            ev_mgr_max_index[i] = index;
          };
        };
      };
    }
    else badrun[index]=true;
  };
  for(Int_t i=0; i<N_TRIG; i++) 
    printf("%s max=%f index=%d\n",trigger_str[i],ev_mgr_max[i],ev_mgr_max_index[i]);
  

  /* compute means and build nv_gr plots (normalised ev_gr)
   * -OR-
   * fit this epoch to an exponential decay function, if the
   * epoch number is in the array exp_fit_list
   */
  Double_t epoch_mean[N_EPOCH][N_TRIG];
  TF1 * epoch_exp[N_EPOCH][N_TRIG];
  char epoch_exp_n[N_EPOCH][N_TRIG][64];
  char exp_formula[N_EPOCH][N_TRIG][64];
  Double_t xx,yy,denom;
  Bool_t use_exp_fit;
  Double_t ave;
  for(Int_t i=0; i<N_TRIG; i++)
  {
    for(Int_t e=0; e<N_EPOCH; e++)
    {
      sprintf(exp_formula[e][i],"[0]*exp([1]*(%d-x))",epoch[e][i]);
      sprintf(epoch_exp_n[e][i],"epoch_exp_e%d_i%d",e,i);
      epoch_exp[e][i] = new TF1(epoch_exp_n[e][i],exp_formula[e][i],epoch[e][i],epoch[e+1][i]);
      epoch_exp[e][i]->SetLineColor(kRed);
      epoch_exp[e][i]->SetLineWidth(2);


      if(ev_gr[e][i]->GetN())
      {
        use_exp_fit=false;
        for(Int_t ee=0; ee<N_EXP_FITS; ee++)
        {
          if(e==exp_fit_list[ee])
          {
            use_exp_fit=true;

            ave = 0;
            for(Int_t pp=0; pp<10; pp++)
            {
              ev_gr[e][i]->GetPoint(0,xx,yy);
              ave += yy;
            };
            ave /= 10.;
            epoch_exp[e][i]->SetParameter(0,ave);
            epoch_exp[e][i]->SetParLimits(1,0,1);
            ev_gr[e][i]->Fit(epoch_exp[e][i],"","",epoch[e][i],epoch[e+1][i]);
          };
        };

        epoch_mean[e][i] = ev_gr[e][i]->GetMean(2);


        for(Int_t p=0; p<ev_gr[e][i]->GetN(); p++)
        {
          ev_gr[e][i]->GetPoint(p,xx,yy);

          if(use_exp_fit) denom = epoch_exp[e][i]->Eval(xx);
          else denom = epoch_mean[e][i];

          yy /= denom;
          nv_gr[e][i]->SetPoint(p,xx,yy);

          // check to see if nv is out of qa cut ranges (for all triggers except LED)
          /*
          if((yy>qa_cut_high[i] || yy<qa_cut_low[i]) && i!=kLED)
          {
            badrun[(Int_t)xx]=true;
            failcode[(Int_t)xx] = failcode[(Int_t)xx] | (Int_t) pow(2,(N_TRIG-i-1));
          };
          */
        };
      };
    };
  };

  
  /* build multigraphs */
  TMultiGraph * ev_mgr[N_TRIG];
  TMultiGraph * nv_mgr[N_TRIG];
  for(Int_t i=0; i<N_TRIG; i++)
  {
    ev_mgr[i] = new TMultiGraph();
    nv_mgr[i] = new TMultiGraph();
    for(Int_t e=0; e<N_EPOCH; e++)
    {
      if(ev_gr[e][i]->GetN())
      {
        ev_mgr[i]->Add(ev_gr[e][i]);
        nv_mgr[i]->Add(nv_gr[e][i]);
      };
    };
  };


  /* multigraph aesthetics */
  char ev_mgr_t[N_TRIG][64];
  char nv_mgr_t[N_TRIG][64];
  for(Int_t i=0; i<N_TRIG; i++)
  {
    sprintf(ev_mgr_t[i],"%s events",trigger_str[i]);
    sprintf(nv_mgr_t[i],"%s epoch-normalised events",trigger_str[i]);
    ev_mgr[i]->SetTitle(ev_mgr_t[i]);
    nv_mgr[i]->SetTitle(nv_mgr_t[i]);
    for(Int_t e=0; e<N_EPOCH; e++)
    {
      ev_gr[e][i]->SetMarkerStyle(kFullCircle);
      nv_gr[e][i]->SetMarkerStyle(kFullCircle);
      ev_gr[e][i]->SetMarkerSize(0.5);
      nv_gr[e][i]->SetMarkerSize(0.5);
      ev_gr[e][i]->SetMarkerColor(trigger_col[i]);
      nv_gr[e][i]->SetMarkerColor(trigger_col[i]);
      for(Int_t b=0; b<N_BAD_EPOCH; b++)
      {
        if(e==bad_epoch[b])
        {
          ev_gr[e][i]->SetMarkerColor(kRed);
          nv_gr[e][i]->SetMarkerColor(kRed);
        };
      };
    };
  };
  for(Int_t i=0; i<N_TRIG; i++)
  {
    ev_mgr[i]->SetMinimum(0);
    ev_mgr[i]->SetMaximum(ev_mgr_max[i]*1.01);
    nv_mgr[i]->SetMinimum(QA_PLOT_MIN);
    nv_mgr[i]->SetMaximum(QA_PLOT_MAX);
    for(Int_t e=0; e<N_EPOCH; e++)
    {
      //nv_gr[e][i]->GetYaxis()->SetLimits(0.1,10);
    };
  };
  


  /* epoch, mean, and QA line aesthetics */
  TLine * epoch_ev_line[N_EPOCH+1][N_TRIG];
  TLine * epoch_nv_line[N_EPOCH+1][N_TRIG];
  TLine * mean_line[N_EPOCH][N_TRIG];
  TLine * qa_line_high[N_TRIG];
  TLine * qa_line_low[N_TRIG];
  TLine * unity_line = new TLine(iterl,1,iterh,1);
  Double_t mean;
  Bool_t isBold;
  for(Int_t i=0; i<N_TRIG; i++)
  {
    for(Int_t e=0; e<N_EPOCH+1; e++)
    {
      epoch_ev_line[e][i] = new TLine(epoch[e][i],0,epoch[e][i],ev_mgr_max[i]*1.01);
      epoch_nv_line[e][i] = new TLine(epoch[e][i],0.1,epoch[e][i],10);
      isBold=false;
      for(Int_t bb=0; bb<N_BOLD_EPOCH; bb++) { if(e==bold_epoch[bb]) isBold=true;};
      if(isBold)
      {
        epoch_ev_line[e][i]->SetLineColor(kRed);
        epoch_nv_line[e][i]->SetLineColor(kRed);
        epoch_ev_line[e][i]->SetLineWidth(3);
        epoch_nv_line[e][i]->SetLineWidth(3);
      }
      else
      {
        epoch_ev_line[e][i]->SetLineColor(kAzure);
        epoch_nv_line[e][i]->SetLineColor(kAzure);
        epoch_ev_line[e][i]->SetLineWidth(2);
        epoch_nv_line[e][i]->SetLineWidth(2);
      };
    };
    for(Int_t e=0; e<N_EPOCH; e++)
    {
      mean = ev_gr[e][i]->GetMean(2);
      printf("%s epoch %d mean = %f\n",trigger_str[i],e,mean);
      mean_line[e][i] = new TLine(epoch[e][i],mean,epoch[e+1][i],mean);
      mean_line[e][i]->SetLineColor(kBlack);
      mean_line[e][i]->SetLineWidth(2);
    };
    qa_line_high[i] = new TLine(iterl,qa_cut_high[i],iterh,qa_cut_high[i]);
    qa_line_low[i] = new TLine(iterl,qa_cut_low[i],iterh,qa_cut_low[i]);
    qa_line_high[i]->SetLineColor(kOrange-6);
    qa_line_low[i]->SetLineColor(kOrange-6);
    qa_line_high[i]->SetLineWidth(2);
    qa_line_low[i]->SetLineWidth(2);
  };
  unity_line->SetLineColor(kBlack);
  unity_line->SetLineWidth(2);


  /*
   * Canvas Layout:
   *
   * 1.JP0     2.LgBS1     3.SmBS1     4.DiBS
   * 5.JP1     6.LgBS2     7.SmBS2     8.DiJP   
   * 9.JP2     10.LgBS3    11.SmBS3    12.LED
   *
   */
  Int_t padn[N_TRIG];
  padn[kLED] = 1;
  padn[kZDC] = 2;
  padn[kVPD] = 3;

  // draw canvases
  Int_t ncolumns=1; // default 4
  TCanvas * ev_canv = new TCanvas("ev_canv","ev_canv",2800,250*((N_TRIG-1)/ncolumns+3));
  TCanvas * nv_canv = new TCanvas("nv_canv","nv_canv",2800,250*((N_TRIG-1)/ncolumns+3));
  ev_canv->Clear();
  nv_canv->Clear();
  ev_canv->Divide(ncolumns,(N_TRIG-1)/ncolumns+1);
  nv_canv->Divide(ncolumns,(N_TRIG-1)/ncolumns+1);
  for(Int_t i=0; i<N_TRIG; i++)
  {
    ev_canv->cd(padn[i]);
    ev_mgr[i]->Draw("AP");
    for(Int_t e=0; e<N_EPOCH+1; e++) 
    {
      epoch_ev_line[e][i]->Draw();
      if(e<N_EPOCH) 
      {
        use_exp_fit=false;
        for(Int_t ee=0; ee<N_EXP_FITS; ee++)
        {
          if(e==exp_fit_list[ee]) use_exp_fit=true;
        };

        if(use_exp_fit) epoch_exp[e][i]->Draw("same");
        else mean_line[e][i]->Draw();
      };
    };

    nv_canv->cd(padn[i]); 
    nv_canv->GetPad(padn[i])->SetLogy();
    nv_mgr[i]->Draw("AP");
    for(Int_t e=0; e<N_EPOCH+1; e++) epoch_nv_line[e][i]->Draw();
    unity_line->Draw();
    qa_line_high[i]->Draw();
    qa_line_low[i]->Draw();
  };


  /* minbias cross sections */
  Double_t xsec[N_MB]; // cross section in millibarns
  //xsec[kBBCMB] = 26.1; // from vernier scan
  Double_t max_xsec[N_MB];
  //max_xsec[kBBCMB] = 30;
  max_xsec[kZDCMB] = 0.2;
  max_xsec[kVPDMB] = 20;
  TH1D * xsec_dist[N_MB]; //superfluously includes BBC xsec distribution for completion (it's unused)
  char xsec_dist_n[N_MB][32];
  char xsec_proj[N_MB][256];
  char xsec_proj_cut[256];
  sprintf(xsec_proj_cut,"ZDC*VPD>0 && time>%d",time_cut);
  for(Int_t i=0; i<N_MB; i++)
  {
    sprintf(xsec_dist_n[i],"%s_xsec_dist",minbias_str[i]);
    xsec_dist[i] = new TH1D(xsec_dist_n[i],xsec_dist_n[i],100,0,max_xsec[i]);
    sprintf(xsec_proj[i],"%f*%s*%sps",26.1,minbias_str[i],minbias_str[i]);
    tr->Project(xsec_dist_n[i],xsec_proj[i],xsec_proj_cut);
  };
  xsec[kZDCMB] = xsec_dist[kZDCMB]->GetMean();
  xsec[kVPDMB] = xsec_dist[kVPDMB]->GetMean();
  printf("\n");
  for(Int_t i=0; i<N_MB; i++) printf("[x] %s Xsec = %.5f mb\n",minbias_str[i],xsec[i]);
  printf("currently QAing with %s\n\n",minbias_str[MB_SELECT]);

  gSystem->RedirectOutput("mb_select","w");
  printf("%s\n",minbias_str[MB_SELECT]);
  gSystem->RedirectOutput(0);


  /* luminosity tracker */
  TGraph * lum_all[N_MB];
  TGraph * lum_good[N_MB];
  TGraph * lum_prog_all[N_MB];
  TGraph * lum_prog_good[N_MB];
  Int_t lum_all_i[N_MB];
  Int_t lum_good_i[N_MB];
  Double_t intlum_all[N_MB];
  Double_t intlum_good[N_MB];
  TLine * lumi_bold_epoch[N_BOLD_EPOCH][N_MB];
  for(Int_t m=0; m<N_MB; m++)
  {
    lum_all[m] = new TGraph();
    lum_good[m] = new TGraph();
    lum_prog_all[m] = new TGraph();
    lum_prog_good[m] = new TGraph();
    lum_all_i[m]=0;
    lum_good_i[m]=0;
    intlum_all[m]=0;
    intlum_good[m]=0;
  };
  for(Int_t i=0; i<tr->GetEntries(); i++)
  {
    tr->GetEntry(i);
    for(Int_t m=0; m<N_MB; m++)
    {
      lum_all[m]->SetPoint(lum_all_i[m],index,mb[m]*mb_ps[m]/(xsec[m]*pow(10,9)));
        // n.b.: xsec in millibarns --> *10^-3; to report lumi in pb^-1, we
        // multiply by factor 10^-9
      intlum_all[m] += mb[m]*mb_ps[m]/(xsec[m]*pow(10,9));
      lum_prog_all[m]->SetPoint(lum_all_i[m],index,intlum_all[m]);
      lum_all_i[m]++;
      if(badrun[index]==false)
      {
        lum_good[m]->SetPoint(lum_good_i[m],index,mb[m]*mb_ps[m]/(xsec[m]*pow(10,9)));
        intlum_good[m] += mb[m]*mb_ps[m]/(xsec[m]*pow(10,9));
        lum_prog_good[m]->SetPoint(lum_good_i[m],index,intlum_good[m]);
        lum_good_i[m]++;
      };
    };
  };
  TMultiGraph * lum_prog[N_MB];
  char lum_prog_t[N_MB][128];
  char lum_good_t[N_MB][128];
  char lum_all_t[N_MB][128];
  for(Int_t m=0; m<N_MB; m++)
  {
    lum_prog[m] = new TMultiGraph();
    sprintf(lum_all_t[m],"%s Luminosity -- All Runs;run index;L_{int} [pb^{-1}]",minbias_str[m]);
    sprintf(lum_good_t[m],"%s Luminosity -- Good Runs;run index;L_{int} [pb^{-1}]",minbias_str[m]);
    lum_all[m]->SetTitle(lum_all_t[m]);
    lum_good[m]->SetTitle(lum_good_t[m]);
    lum_all[m]->SetMarkerStyle(kFullCircle);
    lum_good[m]->SetMarkerStyle(kFullCircle);
    lum_all[m]->SetMarkerSize(0.5);
    lum_good[m]->SetMarkerSize(0.5);
    lum_all[m]->SetMarkerColor(kBlack);
    lum_good[m]->SetMarkerColor(kRed);
    lum_prog_all[m]->SetLineColor(kBlack);
    lum_prog_good[m]->SetLineColor(kRed);
    lum_prog_all[m]->SetLineWidth(2);
    lum_prog_good[m]->SetLineWidth(2);
    lum_prog[m]->Add(lum_prog_all[m]);
    lum_prog[m]->Add(lum_prog_good[m]);
    sprintf(lum_prog_t[m],
      "%s Integrated Luminosity -- (black=all red=good);run index;L_{int} [pb^{-1}]",minbias_str[m]);
    lum_prog[m]->SetTitle(lum_prog_t[m]);
    for(Int_t bb=0; bb<N_BOLD_EPOCH; bb++)
    {
      lumi_bold_epoch[bb][m] = new TLine(epoch[bold_epoch[bb]][0],0,
                                         epoch[bold_epoch[bb]][0],intlum_all[m]);
      lumi_bold_epoch[bb][m]->SetLineColor(kMagenta);
      lumi_bold_epoch[bb][m]->SetLineWidth(3);
    };
  };


  /* luminosity canvas */
  /*
  TCanvas * lum_canv = new TCanvas("lum_canv","lum_canv",1200,400*N_MB);
  lum_canv->Clear();
  lum_canv->Divide(3,N_MB);
  for(Int_t i=1; i<=(3*N_MB); i++) lum_canv->GetPad(i)->SetGrid(1,1);
  for(Int_t m=0; m<N_MB; m++)
  {
    lum_canv->cd(m*3+1); lum_all[m]->Draw("AP");
    lum_canv->cd(m*3+2); lum_good[m]->Draw("AP");
    lum_canv->cd(m*3+3); lum_prog[m]->Draw("AL");
    for(Int_t bb=0; bb<N_BOLD_EPOCH; bb++) lumi_bold_epoch[bb][m]->Draw();
  };
  */


  /* print data table */
  gROOT->ProcessLine(".! touch goodruns.dat ; rm goodruns.dat");
  gROOT->ProcessLine(".! touch badruns.dat ; rm badruns.dat");
  char printstr[2048];
  for(Int_t i=0; i<tr->GetEntries(); i++)
  {
    tr->GetEntry(i);
    sprintf(printstr,"%d %d %d %d",runnum,index,(Int_t)fill,time);
    for(Int_t j=0; j<N_TRIG; j++) sprintf(printstr,"%s %d %.2f",printstr,nev[j],nev_ps[j]);
    //for(Int_t j=0; j<N_MB; j++) sprintf(printstr,"%s %d %.2f",printstr,mb[j],mb_ps[j]);
    sprintf(printstr,"%s %16b %s",printstr,failcode[index],type);
    if(badrun[index]) gSystem->RedirectOutput("badruns.dat");
    else gSystem->RedirectOutput("goodruns.dat");
    printf("%s\n",printstr);
    gSystem->RedirectOutput(0);
  };


  /* ratio plots */
  const Int_t N_RAT_PLOTS = 3;
  Int_t nn[N_RAT_PLOTS]; // numerator plot
  Int_t dd[N_RAT_PLOTS]; // denominator plot
  Color_t ratcol[N_RAT_PLOTS];
  Int_t rpn[N_RAT_PLOTS]; // pad number (4x6 pads in the canvas.. you gotta draw it to understand)
  Int_t ii=0;
  nn[ii] = kLED;  dd[ii] = kZDC;  rpn[ii] = 1;  ratcol[ii++] = kRed;
  nn[ii] = kLED;  dd[ii] = kVPD;  rpn[ii] = 2;  ratcol[ii++] = kBlue;
  nn[ii] = kZDC;  dd[ii] = kVPD;  rpn[ii] = 3;  ratcol[ii++] = kGreen+2;

  
  /* compute ratios */
  TGraph * rat_gr[N_RAT_PLOTS];
  Int_t rat_gr_i[N_RAT_PLOTS];
  char rat_gr_t[N_RAT_PLOTS][64];
  Double_t yynn[N_RUNS];
  Double_t yydd[N_RUNS];
  for(Int_t i=0; i<N_RAT_PLOTS; i++)
  {
    rat_gr[i] = new TGraph();
    rat_gr_i[i]=0;

    for(Int_t r=0; r<N_RUNS; r++)
    {
      yynn[r] = 0;
      yydd[r] = 0;
    };

    // get numerator points
    for(Int_t e=0; e<N_EPOCH; e++)
    {
      for(Int_t p=0; p<ev_gr[e][nn[i]]->GetN(); p++)
      {
        ev_gr[e][nn[i]]->GetPoint(p,xx,yy);
        yynn[(Int_t)xx] = yy;
      };
      for(Int_t p=0; p<ev_gr[e][dd[i]]->GetN(); p++)
      {
        ev_gr[e][dd[i]]->GetPoint(p,xx,yy);
        yydd[(Int_t)xx] = yy;
      };
    };

    // compute ratios
    for(Int_t r=0; r<N_RUNS; r++)
    {
      if(yynn[r]*yydd[r]>0) 
      {
        rat_gr[i]->SetPoint(rat_gr_i[i],r,yynn[r]/yydd[r]);
        (rat_gr_i[i])++;
      };
    };

    // aesthetics
    sprintf(rat_gr_t[i],"%s / %s trigger rato",trigger_str[nn[i]],trigger_str[dd[i]]);
    rat_gr[i]->SetTitle(rat_gr_t[i]);
    rat_gr[i]->SetMarkerStyle(kFullCircle);
    rat_gr[i]->SetMarkerSize(0.5);
    rat_gr[i]->SetMarkerColor(ratcol[i]);
  };


  /* ratio canvas layout */
  TCanvas * rat_canv = new TCanvas("rat_canv","rat_canv",2800,250*((N_RAT_PLOTS-1)/ncolumns+3));
  rat_canv->Clear();
  rat_canv->Divide(ncolumns,(N_RAT_PLOTS-1)/ncolumns+1);
  for(Int_t i=0; i<N_RAT_PLOTS; i++)
  {
    rat_canv->cd(rpn[i]);
    rat_canv->GetPad(rpn[i])->SetGrid(1,1);
    rat_gr[i]->Draw("AP");
  };



  /* prescale vs. run index plots */
  TGraph * prescale_gr[N_TRIG];
  Int_t prescale_gr_i[N_TRIG];
  char prescale_gr_t[N_TRIG][64];
  for(Int_t i=0; i<N_TRIG; i++)
  {
    prescale_gr[i] = new TGraph();
    prescale_gr_i[i] = 0;
    sprintf(prescale_gr_t[i],"%s prescale factor",trigger_str[i]);
    prescale_gr[i]->SetTitle(prescale_gr_t[i]);
    prescale_gr[i]->SetMarkerStyle(kFullCircle);
    prescale_gr[i]->SetMarkerSize(0.5);
    prescale_gr[i]->SetMarkerColor(trigger_col[i]);
  };

  // set JP1 plot draw limit so we can see the prescale for relevant runs
  //prescale_gr[kJP1]->SetMinimum(0);
  //prescale_gr[kJP1]->SetMaximum(20);


  for(Int_t i=0; i<tr->GetEntries(); i++)
  {
    tr->GetEntry(i);
    if(badrun[index] != 1)
    {
      for(Int_t j=0; j<N_TRIG; j++)
      {
        prescale_gr[j]->SetPoint(prescale_gr_i[j],index,nev_ps[j]);
        (prescale_gr_i[j])++;
      };
    };
  };

  
  /* draw prescale canvas */
  TCanvas * prescale_canv = new TCanvas("prescale_canv","prescale_canv",2800,250*((N_TRIG-1)/ncolumns+3));
  prescale_canv->Clear();
  prescale_canv->Divide(ncolumns,(N_TRIG-1)/ncolumns+1);
  for(Int_t i=0; i<N_TRIG; i++)
  {
    prescale_canv->cd(padn[i]);
    prescale_canv->GetPad(padn[i])->SetGrid(1,1);
    prescale_gr[i]->Draw("AP");
  };
  

  /* write run index table */
  Int_t good_count=0;
  for(Int_t r=0; r<N_RUNS; r++)
  {
    if(badrun[r]==false) good_count++;
  };
  gROOT->ProcessLine(".! touch runindex.txt ; rm runindex.txt");
  gSystem->RedirectOutput("runindex.txt");
  /*
  for(Int_t m=0; m<N_MB; m++)
  {
    printf("%s estimated integrated luminosity = %s*ps/(%.3f mb)\n",minbias_str[m],minbias_str[m],xsec[m]);
    printf("--> for all runs: %f pb<sup>-1</sup>\n",intlum_all[m]);
    printf("--> for good runs: %f pb<sup>-1</sup>\n\n",intlum_good[m]);
  };
  */
  printf("Total number of runs: %d\n",N_RUNS);
  printf("Number of good runs: %d\n",good_count);
  printf("\nRun Index <--> Run Number\n");
  tr->SetScanField(0);
  tr->Scan("index:runnum","");
  gSystem->RedirectOutput(0);
  printf("---> Run Index: runindex.txt\n");


  /* print pngs */
  ev_canv->Print("ev_canv.png","png");
  nv_canv->Print("nv_canv.png","png");
  //lum_canv->Print("lum_canv.png","png");
  rat_canv->Print("rat_canv.png","png");
  prescale_canv->Print("prescale_canv.png","png");

  /* print exponential fits */
  Int_t eee;
  for(Int_t i=0; i<N_TRIG; i++)
  {
    for(Int_t ee=0; ee<N_EXP_FITS; ee++)
    {
      eee = exp_fit_list[ee];
      if(epoch_exp[eee][i]!=NULL)
      {
        printf("%s epoch=%d N0=%f s=%f\n",trigger_str[i],eee,
          epoch_exp[eee][i]->GetParameter(0),
          epoch_exp[eee][i]->GetParameter(1));
      };
    };
  };
      
    


  
  /* write rootfile */
  ev_canv->Write();
  nv_canv->Write();
  //lum_canv->Write();
  rat_canv->Write();
  prescale_canv->Write();
  tr->Write("tr");
  fout->Close();
};
