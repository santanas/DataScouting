{

  gROOT->Reset();

  // ** Run2012BCD - dijet **
  char input_root_file[500] = "/afs/cern.ch/work/s/santanas/Workspace/DiJetSearch2012/Scouting_8TeV_2012/histoRunBCD.root";
  // ** Run2012B - dijet **
  //char input_root_file[500] = "/afs/cern.ch/work/s/santanas/Workspace/DiJetSearch2012/Scouting_8TeV_2012/histoRunB.root";
  // ** Run2012C - dijet **
  //  char input_root_file[500] = "/afs/cern.ch/work/s/santanas/Workspace/DiJetSearch2012/Scouting_8TeV_2012/histoRunC.root";
  // ** Run2012D - dijet **
  //  char input_root_file[500] = "/afs/cern.ch/work/s/santanas/Workspace/DiJetSearch2012/Scouting_8TeV_2012/histoRunD.root";
  char input_directory[500] = "scoutingDiJetVariables";//new Bora

  //======================================================
  // ** OLD **
  // Run2012B - all analyses 
  //char input_root_file[500] = "root://eoscms//eos/cms/store/cmst3/user/santanas/DataScouting/DQM_histograms/DataScouting_V00-01-06_Run2012B_runrange_193752-197044_dijet_alfaT_razor_dijetpairs_trijetpairs.root";
  // Run2012C - all analyses
  //char input_root_file[500] = "root://eoscms//eos/cms/store/cmst3/user/santanas/DataScouting/DQM_histograms/DataScouting_V00-01-06_Run2012C_runrange_197885-203755_dijet_alfaT_razor_dijetpairs_trijetpairs.root";
  // Run2012B+Run2012C - all analyses
  //char input_root_file[500] = "root://eoscms//eos/cms/store/cmst3/user/santanas/DataScouting/DQM_histograms/DataScouting_V00-01-06_Run2012B_Run2012C_runrange_193752-203755_dijet_alfaT_razor_dijetpairs_trijetpairs.root";
  //char input_directory[500] = "DQMData/Run 999999/DataScouting/Run summary/DiJet";//old CVS
  //======================================================


  TFile *file0=TFile::Open( input_root_file );
  TDirectoryFile* DQMData_Merged_Runs_DataScouting_Run_summary_DiJet = (TDirectoryFile*) file0->Get( input_directory );

  file0->ls();

  //--------------------------
  TCanvas c1("c1","c1",800,600);
  TH1D* h1_cutFlow = (TH1D*) DQMData_Merged_Runs_DataScouting_Run_summary_DiJet->Get( "h1_cutFlow;1" );  
  h1_cutFlow->SetStats(0);
  h1_cutFlow->Draw();
  c1.SaveAs("cutFlow_Run2012BCD.png");

  //--------------------------
  TCanvas c2("c2","c2",800,600);
  TH1D* h1_selJets_emEnergyFraction = (TH1D*) DQMData_Merged_Runs_DataScouting_Run_summary_DiJet->Get( "h1_selJets_emEnergyFraction;1" );  
  h1_selJets_emEnergyFraction->SetStats(0);
  h1_selJets_emEnergyFraction->Draw();
  c2.SaveAs("EM_Energy_Fraction_Run2012BCD.png");

  //--------------------------
  TCanvas c3("c3","c3",800,600);
  TH1D* h1_selJets_hadEnergyFraction = (TH1D*) DQMData_Merged_Runs_DataScouting_Run_summary_DiJet->Get( "h1_selJets_hadEnergyFraction;1" );  
  h1_selJets_hadEnergyFraction->SetStats(0);
  h1_selJets_hadEnergyFraction->Draw();
  c3.SaveAs("HAD_Energy_Fraction_Run2012BCD.png");

  //--------------------------
  TCanvas c4("c4","c4",800,600);
  c4.SetLogy();
  TH1D* h1_DphijjWide_finalSel = (TH1D*) DQMData_Merged_Runs_DataScouting_Run_summary_DiJet->Get( "h1_DphijjWide_finalSel;1" );  
  h1_DphijjWide_finalSel->SetStats(0);
  h1_DphijjWide_finalSel->Draw();
  c4.SaveAs("dphi_finalSel_Run2012BCD.png");

  //--------------------------
  TCanvas c5("c5","c5",800,600);
  c5.SetLogz();
  TH2D* h2_metVSmetclean = (TH2D*) DQMData_Merged_Runs_DataScouting_Run_summary_DiJet->Get( "h2_metVSmetclean;1" );  
  h2_metVSmetclean->SetStats(0);
  h2_metVSmetclean->Draw("colz");
  c5.SaveAs("met_vs_metClean_Run2012BCD.png");

  //--------------------------
  TCanvas c6("c6","c6",800,600);
  c6.SetLogy();
  TH1D* h1_HT_inclusive = (TH1D*) DQMData_Merged_Runs_DataScouting_Run_summary_DiJet->Get( "h1_HT_inclusive;1" );  
  TH1D* h1_HT_finalSel = (TH1D*) DQMData_Merged_Runs_DataScouting_Run_summary_DiJet->Get( "h1_HT_finalSel;1" );  
  h1_HT_inclusive->GetYaxis()->SetRangeUser(1,h1_HT_inclusive->GetMaximum()*10);
  h1_HT_inclusive->SetLineColor(kRed);
  h1_HT_inclusive->Draw();
  c6.Update();
  statBox_h1_HT_inclusive = (TPaveStats*)h1_HT_inclusive->GetListOfFunctions()->FindObject("stats");
  statBox_h1_HT_inclusive->SetX1NDC(0.28392);
  statBox_h1_HT_inclusive->SetX2NDC(0.582915);
  statBox_h1_HT_inclusive->SetY1NDC(0.480769);
  statBox_h1_HT_inclusive->SetY2NDC(0.86014);
  statBox_h1_HT_inclusive->SetLineColor(kRed);
  statBox_h1_HT_inclusive->SetTextColor(kRed);
  //
  c6.Modified();
  c6.Update();
  //
  h1_HT_finalSel->SetLineColor(kBlue);
  h1_HT_finalSel->Draw("sames");
  c6.Update();
  statBox_h1_HT_finalSel = (TPaveStats*)h1_HT_finalSel->GetListOfFunctions()->FindObject("stats");
  statBox_h1_HT_finalSel->SetX1NDC(0.59);
  statBox_h1_HT_finalSel->SetX2NDC(0.8889);
  statBox_h1_HT_finalSel->SetY1NDC(0.480769);
  statBox_h1_HT_finalSel->SetY2NDC(0.86014);
  statBox_h1_HT_finalSel->SetLineColor(kBlue);
  statBox_h1_HT_finalSel->SetTextColor(kBlue);
  //
  c6.Modified();
  c6.Update();
  c6.SaveAs("HT_noiseRejection_Run2012BCD.png");

  //--------------------------
  TCanvas c7("c7","c7",800,600);
  c7.SetLogy();
  TH1D* h1_MjjWide_finalSel_WithoutNoiseFilter_varbin = (TH1D*) DQMData_Merged_Runs_DataScouting_Run_summary_DiJet->Get( "h1_MjjWide_finalSel_WithoutNoiseFilter_varbin;1" );  
  TH1D* h1_MjjWide_finalSel_varbin = (TH1D*) DQMData_Merged_Runs_DataScouting_Run_summary_DiJet->Get( "h1_MjjWide_finalSel_varbin;1" );  
  h1_MjjWide_finalSel_WithoutNoiseFilter_varbin->GetYaxis()->SetRangeUser(1,h1_MjjWide_finalSel_WithoutNoiseFilter_varbin->GetMaximum()*10);
  h1_MjjWide_finalSel_WithoutNoiseFilter_varbin->SetLineColor(kRed);
  h1_MjjWide_finalSel_WithoutNoiseFilter_varbin->Draw();
  c7.Update();
  statBox_h1_MjjWide_finalSel_WithoutNoiseFilter_varbin = (TPaveStats*)h1_MjjWide_finalSel_WithoutNoiseFilter_varbin->GetListOfFunctions()->FindObject("stats");
  statBox_h1_MjjWide_finalSel_WithoutNoiseFilter_varbin->SetX1NDC(0.28392);
  statBox_h1_MjjWide_finalSel_WithoutNoiseFilter_varbin->SetX2NDC(0.582915);
  statBox_h1_MjjWide_finalSel_WithoutNoiseFilter_varbin->SetY1NDC(0.480769);
  statBox_h1_MjjWide_finalSel_WithoutNoiseFilter_varbin->SetY2NDC(0.86014);
  statBox_h1_MjjWide_finalSel_WithoutNoiseFilter_varbin->SetLineColor(kRed);
  statBox_h1_MjjWide_finalSel_WithoutNoiseFilter_varbin->SetTextColor(kRed);
  //
  c7.Modified();
  c7.Update();
  //
  h1_MjjWide_finalSel_varbin->SetLineColor(kBlue);
  h1_MjjWide_finalSel_varbin->Draw("sames");
  c7.Update();
  statBox_h1_MjjWide_finalSel_varbin = (TPaveStats*)h1_MjjWide_finalSel_varbin->GetListOfFunctions()->FindObject("stats");
  statBox_h1_MjjWide_finalSel_varbin->SetX1NDC(0.59);
  statBox_h1_MjjWide_finalSel_varbin->SetX2NDC(0.8889);
  statBox_h1_MjjWide_finalSel_varbin->SetY1NDC(0.480769);
  statBox_h1_MjjWide_finalSel_varbin->SetY2NDC(0.86014);
  statBox_h1_MjjWide_finalSel_varbin->SetLineColor(kBlue);
  statBox_h1_MjjWide_finalSel_varbin->SetTextColor(kBlue);
  //
  c7.Modified();
  c7.Update();
  c7.SaveAs("Mjj_noiseRejection_Run2012BCD.png");


}
