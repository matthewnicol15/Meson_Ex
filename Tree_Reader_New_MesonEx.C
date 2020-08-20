#include <cstdlib>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TChain.h>
#include <TCanvas.h>
#include <vector>
#include <TLorentzVector.h>


double Calc_dtfInterDOCA(const TVector3 &locUnitDir1, const TVector3 &locUnitDir2, const TVector3 &locVertex1, const TVector3 &locVertex2, TVector3 &locInterDOCA1, TVector3 &locInterDOCA2){
  //originated from code by Jörn Langheinrich
  //you can use this function to find the DOCA to a fixed point by calling this function with locUnitDir1 and 2 parallel, and the fixed vertex as locVertex2
  double locUnitDot = locUnitDir1*locUnitDir2;
  double locDenominator = locUnitDot*locUnitDot - 1.0; /// scalar product of directions
  double locDistVertToInterDOCA1 = 0.0, locDistVertToInterDOCA2 = 0.0; //distance from vertex to DOCA point


  if(fabs(locDenominator) < 1.0e-15) //parallel
  locDistVertToInterDOCA1 = (locVertex2 - locVertex1)*locUnitDir2/locUnitDot; //the opposite
  else{
    double locA = (locVertex1 - locVertex2)*locUnitDir1;
    double locB = (locVertex1 - locVertex2)*locUnitDir2;
    locDistVertToInterDOCA1 = (locA - locUnitDot*locB)/locDenominator;
    locDistVertToInterDOCA2 = (locUnitDot*locA - locB)/locDenominator;
  }


  locInterDOCA1 = locVertex1 + locDistVertToInterDOCA1*locUnitDir1; //intersection point of DOCA line and track 1
  locInterDOCA2 = locVertex2 + locDistVertToInterDOCA2*locUnitDir2; //intersection point of DOCA line and track 2
  double locDOCA = (locInterDOCA1 - locInterDOCA2).Mag();
  return ((locVertex2.Z() > locVertex1.Z()) ? locDOCA : -1.0*locDOCA);
}

TLorentzVector  Correct_Electron_g(TLorentzVector x){

  Double_t E_new, Px_el, Py_el, Pz_el;
  TLorentzVector el_new;

  E_new = x.E()-0.03689+0.1412*x.E()-0.04316*pow(x.E(),2)+0.007046*pow(x.E(),3)-0.0004055*pow(x.E(),4);

  Px_el = E_new*(x.Px()/x.Rho());
  Py_el = E_new*(x.Py()/x.Rho());
  Pz_el = E_new*(x.Pz()/x.Rho());

  el_new.SetXYZM(Px_el, Py_el, Pz_el, 0.000511);

  return el_new;
}

void Tree_Reader_New_MesonEx(){

  gROOT->ProcessLine(".L ./Loader.C+");
  TFile *f = new TFile("Trees/skim11_Pass1_Tree_090620_10.root");
  TTree *t1 = (TTree*)f->Get("skim11_Tree_090620_10");


  vector<TLorentzVector> *v_p4=0;

  TLorentzVector *readbeam=NULL;
  TLorentzVector *readtarget=NULL;

  vector<TLorentzVector> *v_vertex=0;

  vector<double> *v_path=0;
  Double_t path;

  vector<double> *v_time=0;
  Double_t time;

  vector<double> *v_beta=0;

  Double_t start_time;
  vector<double> *v_energy=0;
  vector<Int_t> *v_charge=0;
  vector<Int_t> *v_PID=0;
  vector<double> *v_chi2PID=0;
  vector<double> *v_status=0;

  Int_t readchargetracks;
  Int_t readprotonno;
  Int_t readneutronno;
  Int_t readphotonno;
  Int_t readantiprotonno;
  Int_t readpositronno;
  Int_t readpipno;
  Int_t readpimno;
  Int_t readelno;
  Int_t readkaonpno;
  Int_t readkaonmno;
  Int_t readothertracks;
  Int_t readeventno;
  Int_t readrunno;
  Int_t readtriggerno;

  t1->SetBranchAddress("eventno",&readeventno);
  t1->SetBranchAddress("runno",&readrunno);
  t1->SetBranchAddress("triggerno",&readtriggerno);
  t1->SetBranchAddress("start_time",&start_time);
  t1->SetBranchAddress("p4",&v_p4);
  t1->SetBranchAddress("vertex",&v_vertex);
  t1->SetBranchAddress("beta",&v_beta);
  t1->SetBranchAddress("status",&v_status);
  t1->SetBranchAddress("energy",&v_energy);
  t1->SetBranchAddress("charge",&v_charge);
  t1->SetBranchAddress("PID",&v_PID);
  t1->SetBranchAddress("chi2PID",&v_chi2PID);
  t1->SetBranchAddress("time",&v_time);
  t1->SetBranchAddress("path",&v_path);
  t1->SetBranchAddress("beam",&readbeam);
  t1->SetBranchAddress("target",&readtarget);
  t1->SetBranchAddress("elno",&readelno);
  t1->SetBranchAddress("positronno",&readpositronno);
  t1->SetBranchAddress("protonno",&readprotonno);
  t1->SetBranchAddress("antiprotonno",&readantiprotonno);
  t1->SetBranchAddress("neutronno",&readneutronno);
  t1->SetBranchAddress("photonno",&readphotonno);
  t1->SetBranchAddress("pipno",&readpipno);
  t1->SetBranchAddress("pimno",&readpimno);
  t1->SetBranchAddress("kaonpno",&readkaonpno);
  t1->SetBranchAddress("kaonmno",&readkaonmno);

  TFile fileOutput1("New_Data/skim11_MesonEx_110620_01.root","recreate");

  //Creating histograms
  // Tests to check how many particles are in each event
  auto* hkaonpno_1=new TH1F("hkaonpno_1","Number of K^{+} (1);Number of K^{+};Counts",9,0,8);
  auto* hkaonmno_1=new TH1F("hkaonmno_1","Number of K^{-} (1);Number of K^{-};Counts",9,0,8);
  auto* hpimno_1=new TH1F("hpimno_1","Number of pi^{-} (1);Number of pi^{-};Counts",9,0,8);
  auto* hpipno_1=new TH1F("hpipno_1","Number of pi^{+} (1);Number of pi^{+};Counts",9,0,8);
  auto* hprotonno_1=new TH1F("hprotonno_1","Number of p (1);Number of p;Counts",9,0,8);
  auto* helectronno_1=new TH1F("helectronno_1","Number of e^{-} (1);Number of e^{-};Counts",9,0,8);

  //PID only cut
  auto* hmiss_mass_all=new TH1F("miss_all","MM^2(e' K^{+} K^{-} p);MM^2(e' K^{+} K^{-} p) [GeV];Counts",400,-1,4);
  auto* hmiss_momentum_all=new TH1F("hmiss_momentum_all","P(B + T - e' - K^{+} - K^{-} - p);P(B + T - e' - K^{+} - K^{-} - p) [GeV];Counts",400,-1,4);
  auto* hinv_K_P_K_M=new TH1F("hinv_K_P_K_M","M(K^{+} K^{-});M(K^{+} K^{-}) [GeV];Counts",400,-1,4);
  auto* hpath_kp=new TH1F("hpath_kp","Path of K^{+};Path of K^{+} ;Counts",500,-1000,1000);
  auto* hvertex_time_kp=new TH1F("hvertex_time_kp","Vertex time K^{+};Vertex time K^{+} (1);Counts",500,-1000,1000);




  //delta beta and momentum cuts on kaons and protons
  auto* hmiss_mass_allc=new TH1F("miss_allc","MM^2(e' K^{+} p #pi^{-}) after cuts on delta beta and momentum of K^{+} and p;MM^2(e' K^{+} p #pi^{-}) [GeV];Counts",400,-1,4);
  auto* hmiss_momentum_allc=new TH1F("hmiss_momentum_allc","P(B + T - e' - K^{+} - p - #pi^{-}) after cuts on delta beta and momentum of K^{+} and p;P(B + T - e' - K^{+} - p - #pi^{-}) [GeV];Counts",400,-1,4);
  auto* hinv_K_P_K_Mc=new TH1F("hinv_K_P_K_Mc","M(K^{+} K^{-});M(K^{+} K^{-}) [GeV];Counts",400,-1,4);
  auto* hvertex_time_kpc=new TH1F("hvertex_time_kpc","Vertex time K^{+};Vertex time K^{+} (1);Counts",500,-1000,1000);



  //Exclusivity cuts as well
  auto* hinv_K_P_K_Mt=new TH1F("hinv_K_P_K_Mt","M(K^{+} K^{-});M(K^{+} K^{-}) [GeV];Counts",400,-1,4);
  auto* hvertex_time_kpt=new TH1F("hvertex_time_kpt","Vertex time K^{+};Vertex time K^{+} (1);Counts",500,-1000,1000);




  vector<TLorentzVector> v_el;
  vector<TLorentzVector> v_pip;
  vector<TLorentzVector> v_pim;
  vector<TLorentzVector> v_pr;
  vector<TLorentzVector> v_kp;
  vector<TLorentzVector> v_km;
  vector<TLorentzVector> v_othertracks;
  vector<TLorentzVector> v_vertex_el;
  vector<TLorentzVector> v_vertex_pip;
  vector<TLorentzVector> v_vertex_pim;
  vector<TLorentzVector> v_vertex_pr;
  vector<TLorentzVector> v_vertex_kp;
  vector<TLorentzVector> v_vertex_km;

  TLorentzVector missall;
  TLorentzVector missallcasc;
  TLorentzVector inv_K_P_K_M;
  TLorentzVector vertex_el;
  TLorentzVector vertex_pip;
  TLorentzVector vertex_pim;
  TLorentzVector vertex_pr;
  TLorentzVector vertex_kp;
  TLorentzVector vertex_km;

  vector<Double_t> v_beta_tof_el;
  Double_t beta_tof_el;
  vector<Double_t> v_P_el;
  Double_t P_el;
  vector<Double_t>v_path_el;
  Double_t path_el;
  vector<Double_t>v_TOF_el;
  Double_t TOF_el;
  vector<Double_t> v_beta_calc_el;
  Double_t beta_calc_el;
  vector<Double_t> v_delta_beta_el;
  Double_t delta_beta_el;
  vector<Double_t> v_vt_el;
  Double_t vt_el;
  vector<Double_t> v_vertex_time_el;
  Double_t vertex_time_el;

  vector<Double_t> v_beta_tof_pip;
  Double_t beta_tof_pip;
  vector<Double_t> v_P_pip;
  Double_t P_pip;
  vector<Double_t>v_path_pip;
  Double_t path_pip;
  vector<Double_t>v_TOF_pip;
  Double_t TOF_pip;
  vector<Double_t> v_beta_calc_pip;
  Double_t beta_calc_pip;
  vector<Double_t> v_delta_beta_pip;
  Double_t delta_beta_pip;
  vector<Double_t> v_vt_pip;
  Double_t vt_pip;
  vector<Double_t> v_vertex_time_pip;
  Double_t vertex_time_pip;

  vector<Double_t> v_beta_tof_pim;
  Double_t beta_tof_pim;
  vector<Double_t> v_P_pim;
  Double_t P_pim;
  vector<Double_t>v_path_pim;
  Double_t path_pim;
  vector<Double_t>v_TOF_pim;
  Double_t TOF_pim;
  vector<Double_t> v_beta_calc_pim;
  Double_t beta_calc_pim;
  vector<Double_t> v_delta_beta_pim;
  Double_t delta_beta_pim;
  vector<Double_t> v_vt_pim;
  Double_t vt_pim;
  vector<Double_t> v_vertex_time_pim;
  Double_t vertex_time_pim;

  vector<Double_t> v_beta_tof_pr;
  Double_t beta_tof_pr;
  vector<Double_t> v_P_pr;
  Double_t P_pr;
  vector<Double_t>v_path_pr;
  Double_t path_pr;
  vector<Double_t>v_TOF_pr;
  Double_t TOF_pr;
  vector<Double_t> v_beta_calc_pr;
  Double_t beta_calc_pr;
  vector<Double_t> v_delta_beta_pr;
  Double_t delta_beta_pr;
  vector<Double_t> v_vt_pr;
  Double_t vt_pr;
  vector<Double_t> v_vertex_time_pr;
  Double_t vertex_time_pr;

  vector<Double_t> v_beta_tof_kp;
  Double_t beta_tof_kp;
  vector<Double_t> v_P_kp;
  Double_t P_kp;
  vector<Double_t>v_path_kp;
  Double_t path_kp;
  vector<Double_t>v_TOF_kp;
  Double_t TOF_kp;
  vector<Double_t> v_beta_calc_kp;
  Double_t beta_calc_kp;
  vector<Double_t> v_delta_beta_kp;
  Double_t delta_beta_kp;
  vector<Double_t> v_vt_kp;
  Double_t vt_kp;
  vector<Double_t> v_time_kp;
  Double_t time_kp;
  vector<Double_t> v_vertex_time_kp;
  Double_t vertex_time_kp;

  vector<Double_t> v_beta_tof_km;
  Double_t beta_tof_km;
  vector<Double_t> v_P_km;
  Double_t P_km;
  vector<Double_t>v_path_km;
  Double_t path_km;
  vector<Double_t>v_TOF_km;
  Double_t TOF_km;
  vector<Double_t> v_beta_calc_km;
  Double_t beta_calc_km;
  vector<Double_t> v_delta_beta_km;
  Double_t delta_beta_km;
  vector<Double_t> v_vt_km;
  Double_t vt_km;
  vector<Double_t> v_vertex_time_km;
  Double_t vertex_time_km;

  vector<Double_t> v_beta_tof_othertracks;
  Double_t beta_tof_othertracks;
  vector<Double_t> v_P_othertracks;
  Double_t P_othertracks;
  vector<Double_t> v_beta_calc_othertracks;
  Double_t beta_calc_othertracks;
  vector<Double_t> v_delta_beta_othertracks;
  Double_t delta_beta_othertracks;
  vector<Double_t> v_vt_othertracks;
  Double_t vt_othertracks;

  Double_t delta_theta, delta_phi_S1_pr_pim, delta_P, Delta_Phi_Unidentified;

  Double_t Delta_Vertex_Time;

  Double_t c=30;

  TLorentzVector el;
  TLorentzVector el_corrected;
  TLorentzVector pip;
  TLorentzVector pim;
  TLorentzVector pr;
  TLorentzVector kp;
  TLorentzVector km;
  TLorentzVector othertracks;


  double doca_S2_kp_1_e_scattered, doca_S2_p_pi_1, doca_S2_p_pi_2;
  TVector3 old_vert_S2_kp_1, old_vert_S2_kp_2, old_vert_S2_e_scattered, old_vert_S2_p, old_vert_S2_pi_1,old_vert_S2_pi_2, old_dir_S2_kp_1, old_dir_S2_kp_2, old_dir_S2_e_scattered, old_dir_S2_p, old_dir_S2_pi_1,old_dir_S2_pi_2, new_vert_S2_kp_1, new_vert_S2_e_scattered, new_vert_S2_p, new_vert_S2_pi_1, new_vert_S2_pi_2;
  TVector3 loc_DOCA_vertex_1, loc_DOCA_vertex_2;
  new_vert_S2_kp_1.SetXYZ(0.0,0.0,0.0);
  new_vert_S2_e_scattered.SetXYZ(0.0,0.0,0.0);
  new_vert_S2_p.SetXYZ(0.0,0.0,0.0);
  new_vert_S2_pi_1.SetXYZ(0.0,0.0,0.0);
  new_vert_S2_pi_2.SetXYZ(0.0,0.0,0.0);

  Int_t Total_Events = 0;
  Int_t Good_Events = 0;
  Int_t empty = 0;

  Long64_t nentries = t1->GetEntries();
  Int_t Percentage = nentries/100;
  for(Long64_t i=0; i<nentries;i++){
    t1->GetEntry(i);
    if (i % Percentage == 0){
      fprintf (stderr, "%lld\r", i/Percentage);
      fflush (stderr);
    }


    Total_Events++;

    v_el.clear();
    v_beta_tof_el.clear();
    v_P_el.clear();
    v_path_el.clear();
    v_TOF_el.clear();
    v_beta_calc_el.clear();
    v_delta_beta_el.clear();
    v_vertex_time_el.clear();
    v_vertex_el.clear();

    v_pip.clear();
    v_beta_tof_pip.clear();
    v_P_pip.clear();
    v_path_pip.clear();
    v_TOF_pip.clear();
    v_beta_calc_pip.clear();
    v_delta_beta_pip.clear();
    v_vertex_pip.clear();
    v_vertex_time_pip.clear();
    v_vertex_pip.clear();


    v_pim.clear();
    v_beta_tof_pim.clear();
    v_P_pim.clear();
    v_path_pim.clear();
    v_TOF_pim.clear();
    v_beta_calc_pim.clear();
    v_delta_beta_pim.clear();
    v_vertex_pim.clear();
    v_vertex_time_pim.clear();
    v_vertex_pim.clear();

    v_pr.clear();
    v_beta_tof_pr.clear();
    v_P_pr.clear();
    v_path_pr.clear();
    v_TOF_pr.clear();
    v_beta_calc_pr.clear();
    v_delta_beta_pr.clear();
    v_vertex_pr.clear();
    v_vertex_time_pr.clear();
    v_vertex_pr.clear();

    v_kp.clear();
    v_beta_tof_kp.clear();
    v_P_kp.clear();
    v_path_kp.clear();
    v_TOF_kp.clear();
    v_beta_calc_kp.clear();
    v_delta_beta_kp.clear();
    v_time_kp.clear();
    v_vertex_kp.clear();
    v_vertex_time_kp.clear();
    v_vertex_kp.clear();


    v_km.clear();
    v_beta_tof_km.clear();
    v_P_km.clear();
    v_path_km.clear();
    v_TOF_km.clear();
    v_beta_calc_km.clear();
    v_delta_beta_km.clear();
    v_vertex_km.clear();
    v_vertex_time_km.clear();
    v_vertex_km.clear();


    v_othertracks.clear();
    v_beta_tof_othertracks.clear();
    v_P_othertracks.clear();
    v_beta_calc_othertracks.clear();
    v_delta_beta_othertracks.clear();

    Int_t Nparticles = v_p4->size();
    if(Nparticles==0)continue;
    for(Int_t j=0; j<Nparticles; j++){

      if(v_PID->at(j)==11){
        el.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),v_p4->at(j).M());
        P_el = sqrt((pow(v_p4->at(j).Px(),2))+(pow(v_p4->at(j).Py(),2))+(pow(v_p4->at(j).Pz(),2)));
        TOF_el = v_time->at(j);
        path_el = v_path->at(j);
        beta_tof_el = v_beta->at(j);
        beta_calc_el = P_el/(sqrt((pow(P_el,2))+(pow(el.M(),2))));
        delta_beta_el = beta_calc_el-beta_tof_el;
        vertex_time_el = TOF_el - path_el / (beta_tof_el*c);
        vertex_el.SetXYZT(v_vertex->at(j).X(), v_vertex->at(j).Y(), v_vertex->at(j).Z(), vertex_time_el);

        v_el.push_back(el);
        v_beta_tof_el.push_back(beta_tof_el);
        v_P_el.push_back(P_el);
        v_path_el.push_back(path_el);
        v_TOF_el.push_back(TOF_el);
        v_beta_calc_el.push_back(beta_calc_el);
        v_delta_beta_el.push_back(delta_beta_el);
        v_vertex_time_el.push_back(vertex_time_el);
        v_vertex_el.push_back(vertex_el);
      }
      else if(v_PID->at(j)==211){
        pip.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),v_p4->at(j).M());
        P_pip = sqrt((pow(v_p4->at(j).Px(),2))+(pow(v_p4->at(j).Py(),2))+(pow(v_p4->at(j).Pz(),2)));
        TOF_pip = v_time->at(j);
        path_pip = v_path->at(j);
        beta_tof_pip = v_beta->at(j);
        beta_calc_pip = P_pip/(sqrt((pow(P_pip,2))+(pow(pip.M(),2))));
        delta_beta_pip = beta_calc_pip-beta_tof_pip;
        vertex_time_pip = TOF_pip - path_pip / (beta_tof_pip*c);
        vertex_pip.SetXYZT(v_vertex->at(j).X(), v_vertex->at(j).Y(), v_vertex->at(j).Z(), vertex_time_pip);

        v_pip.push_back(pip);
        v_beta_tof_pip.push_back(beta_tof_pip);
        v_P_pip.push_back(P_pip);
        v_path_pip.push_back(path_pip);
        v_TOF_pip.push_back(TOF_pip);
        v_beta_calc_pip.push_back(beta_calc_pip);
        v_delta_beta_pip.push_back(delta_beta_pip);
        v_vertex_time_pip.push_back(vertex_time_pip);
        v_vertex_pip.push_back(vertex_pip);
      }
      else if(v_PID->at(j)==-211){
        pim.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),v_p4->at(j).M());
        P_pim = sqrt((pow(v_p4->at(j).Px(),2))+(pow(v_p4->at(j).Py(),2))+(pow(v_p4->at(j).Pz(),2)));
        TOF_pim = v_time->at(j);
        path_pim = v_path->at(j);
        beta_tof_pim = v_beta->at(j);
        beta_calc_pim = P_pim/(sqrt((pow(P_pim,2))+(pow(pim.M(),2))));
        delta_beta_pim = beta_calc_pim-beta_tof_pim;
        vertex_time_pim = TOF_pim - path_pim / (beta_tof_pim*c);
        vertex_pim.SetXYZT(v_vertex->at(j).X(), v_vertex->at(j).Y(), v_vertex->at(j).Z(), vertex_time_pim);

        v_pim.push_back(pim);
        v_beta_tof_pim.push_back(beta_tof_pim);
        v_P_pim.push_back(P_pim);
        v_path_pim.push_back(path_pim);
        v_TOF_pim.push_back(TOF_pim);
        v_beta_calc_pim.push_back(beta_calc_pim);
        v_delta_beta_pim.push_back(delta_beta_pim);
        v_vertex_time_pim.push_back(vertex_time_pim);
        v_vertex_pim.push_back(vertex_pim);
      }
      else if(v_PID->at(j)==2212){
        pr.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),v_p4->at(j).M());
        P_pr = sqrt((pow(v_p4->at(j).Px(),2))+(pow(v_p4->at(j).Py(),2))+(pow(v_p4->at(j).Pz(),2)));
        TOF_pr = v_time->at(j);
        path_pr = v_path->at(j);
        beta_tof_pr = v_beta->at(j);
        beta_calc_pr = P_pr/(sqrt((pow(P_pr,2))+(pow(pr.M(),2))));
        delta_beta_pr = beta_calc_pr-beta_tof_pr;
        vertex_time_pr = TOF_pr - path_pr / (beta_tof_pr*c);
        vertex_pr.SetXYZT(v_vertex->at(j).X(), v_vertex->at(j).Y(), v_vertex->at(j).Z(), vertex_time_pr);

        v_pr.push_back(pr);
        v_beta_tof_pr.push_back(beta_tof_pr);
        v_P_pr.push_back(P_pr);
        v_path_pr.push_back(path_pr);
        v_TOF_pr.push_back(TOF_pr);
        v_beta_calc_pr.push_back(beta_calc_pr);
        v_delta_beta_pr.push_back(delta_beta_pr);
        v_vertex_time_pr.push_back(vertex_time_pr);
        v_vertex_pr.push_back(vertex_pr);
      }
      else if(v_PID->at(j)==321){
        kp.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),v_p4->at(j).M());
        P_kp = sqrt((pow(v_p4->at(j).Px(),2))+(pow(v_p4->at(j).Py(),2))+(pow(v_p4->at(j).Pz(),2)));
        TOF_kp = v_time->at(j);
        path_kp = v_path->at(j);
        beta_tof_kp = v_beta->at(j);
        beta_calc_kp = P_kp/(sqrt((pow(P_kp,2))+(pow(kp.M(),2))));
        delta_beta_kp = beta_calc_kp-beta_tof_kp;
        vertex_time_kp = TOF_kp - path_kp / (beta_tof_kp*c);
        vertex_kp.SetXYZT(v_vertex->at(j).X(), v_vertex->at(j).Y(), v_vertex->at(j).Z(), vertex_time_kp);

        v_kp.push_back(kp);
        v_beta_tof_kp.push_back(beta_tof_kp);
        v_P_kp.push_back(P_kp);
        v_path_kp.push_back(path_kp);
        v_TOF_kp.push_back(TOF_kp);
        v_beta_calc_kp.push_back(beta_calc_kp);
        v_delta_beta_kp.push_back(delta_beta_kp);
        v_vertex_time_kp.push_back(vertex_time_kp);
        v_vertex_kp.push_back(vertex_kp);
      }
      else if(v_PID->at(j)==-321){
        km.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),v_p4->at(j).M());
        P_km = sqrt((pow(v_p4->at(j).Px(),2))+(pow(v_p4->at(j).Py(),2))+(pow(v_p4->at(j).Pz(),2)));
        TOF_km = v_time->at(j);
        path_km = v_path->at(j);
        beta_tof_km = v_beta->at(j);
        beta_calc_km = P_km/(sqrt((pow(P_km,2))+(pow(km.M(),2))));
        delta_beta_km = beta_calc_km-beta_tof_km;
        vertex_time_km = TOF_km - path_km / (beta_tof_km*c);
        vertex_km.SetXYZT(v_vertex->at(j).X(), v_vertex->at(j).Y(), v_vertex->at(j).Z(), vertex_time_km);

        v_km.push_back(km);
        v_beta_tof_km.push_back(beta_tof_km);
        v_P_km.push_back(P_km);
        v_path_km.push_back(path_km);
        v_TOF_km.push_back(TOF_km);
        v_beta_calc_km.push_back(beta_calc_km);
        v_delta_beta_km.push_back(delta_beta_km);
        v_vertex_time_km.push_back(vertex_time_km);
        v_vertex_km.push_back(vertex_km);
      }
      else{
        othertracks.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),v_p4->at(j).M());
        P_othertracks = sqrt((pow(v_p4->at(j).Px(),2))+(pow(v_p4->at(j).Py(),2))+(pow(v_p4->at(j).Pz(),2)));
        beta_tof_othertracks = v_beta->at(j);
        beta_calc_othertracks = P_othertracks/(sqrt((pow(P_othertracks,2))+(pow(othertracks.M(),2))));
        delta_beta_othertracks = beta_calc_othertracks-beta_tof_othertracks;

        v_othertracks.push_back(othertracks);
        v_beta_tof_othertracks.push_back(beta_tof_othertracks);
        v_P_othertracks.push_back(P_othertracks);
        v_beta_calc_othertracks.push_back(beta_calc_othertracks);
        v_delta_beta_othertracks.push_back(delta_beta_othertracks);
      }


    }
    hkaonpno_1->Fill(v_kp.size());
    hkaonmno_1->Fill(v_km.size());
    helectronno_1->Fill(v_el.size());
    hpimno_1->Fill(v_pim.size());
    hprotonno_1->Fill(v_pr.size());

    // if(v_km.size()==1 && v_kp.size()==1){

      // if(readeventno==12101 || readeventno==16263 || readeventno==16299 || readeventno==19714 || readeventno==23224 || readeventno==24508)cout<<"match"<<endl;

    // }

    if(v_el.size()==0)continue;
    el_corrected = Correct_Electron_g(v_el.at(0));


    // // Selecting Events
    if(v_kp.size()==1 && v_km.size()==1 && v_el.size()==1 && v_pr.size()==1){
      // if(i>600)break;
      // cout<<"event no: "<<readeventno<<"uncorrected e' P: "<<v_el.at(0).Rho()<<"uncorrected e' E: "<<v_el.at(0).E()<<endl;

      Good_Events++;
      missall = (TLorentzVector)*readbeam + (TLorentzVector)*readtarget - el_corrected - v_kp.at(0) - v_pr.at(0) - v_km.at(0);
      inv_K_P_K_M = v_kp.at(0) + v_km.at(0);
      // Testing vertex time and path functions
      hpath_kp->Fill(v_path_kp.at(0));
      hvertex_time_kp->Fill(v_vertex_time_kp.at(0));

      hmiss_mass_all->Fill(missall.M2());
      hmiss_momentum_all->Fill(missall.Rho());
      hinv_K_P_K_M->Fill(inv_K_P_K_M.M());

      if(fabs(v_delta_beta_kp.at(0))<0.01 && v_P_kp.at(0)>0.2 && v_P_kp.at(0)<2.6 && fabs(v_delta_beta_pr.at(0))<0.01 && v_P_pr.at(0)>0.2 && v_P_pr.at(0)<4.8){


        hmiss_mass_allc->Fill(missall.M2());
        hmiss_momentum_allc->Fill(missall.Rho());
        hinv_K_P_K_Mc->Fill(inv_K_P_K_M.M());
        hvertex_time_kpc->Fill(v_vertex_time_kp.at(0));


        if(fabs(missall.M2())<0.05 && fabs(missall.Rho())<0.5){
          hinv_K_P_K_Mt->Fill(inv_K_P_K_M.M());
          hvertex_time_kpt->Fill(v_vertex_time_kp.at(0));


        }
      }
    }
    // // Strangeness 2 channels
    // if(v_kp.size()==2 && v_pr.size()==1 && v_el.size()==1 && v_pim.size()==2){
    //   missallcasc = (TLorentzVector)*readbeam + (TLorentzVector)*readtarget - el_corrected - v_kp.at(0) - v_kp.at(1) -v_pr.at(0) - v_pim.at(0) - v_pim.at(1);
    //   miss3 = (TLorentzVector)*readbeam + (TLorentzVector)*readtarget - el_corrected - v_kp.at(0) - v_kp.at(1);
    //   cascade_lambda_1 = v_pr.at(0) + v_pim.at(0);
    //   cascade_lambda_2 = v_pr.at(0) + v_pim.at(1);
    //   cascade = v_pr.at(0) + v_pim.at(0) + v_pim.at(1);
    //   Delta_Vertex_Time = v_vertex_kp.at(0).T() - v_vertex_kp.at(1).T();
    //   hdelta_vertex_time_kp2->Fill(Delta_Vertex_Time);
    //   hvertex_time_kp2->Fill(v_vertex_time_kp.at(0), v_vertex_time_kp.at(1));
    //
    //   old_vert_S2_e_scattered.SetXYZ(v_vertex_el.at(0).X(), v_vertex_el.at(0).Y(),v_vertex_el.at(0).Z());
    //   old_vert_S2_kp_1.SetXYZ(v_vertex_kp.at(0).X(), v_vertex_kp.at(0).Y(),v_vertex_kp.at(0).Z());
    //   old_vert_S2_kp_2.SetXYZ(v_vertex_kp.at(1).X(), v_vertex_kp.at(1).Y(),v_vertex_kp.at(1).Z());
    //   old_vert_S2_p.SetXYZ(v_vertex_pr.at(0).X(), v_vertex_pr.at(0).Y(),v_vertex_pr.at(0).Z());
    //   old_vert_S2_pi_1.SetXYZ(v_vertex_pim.at(0).X(), v_vertex_pim.at(0).Y(),v_vertex_pim.at(0).Z());
    //   old_vert_S2_pi_2.SetXYZ(v_vertex_pim.at(1).X(), v_vertex_pim.at(1).Y(),v_vertex_pim.at(1).Z());
    //   old_dir_S2_e_scattered.SetXYZ(el_corrected.Px()/el_corrected.Rho(), el_corrected.Py()/el_corrected.Rho(), el_corrected.Pz()/el_corrected.Rho());
    //   old_dir_S2_kp_1.SetXYZ(v_kp.at(0).Px()/v_kp.at(0).Rho(), v_kp.at(0).Py()/v_kp.at(0).Rho(), v_kp.at(0).Pz()/v_kp.at(0).Rho());
    //   old_dir_S2_kp_2.SetXYZ(v_kp.at(1).Px()/v_kp.at(1).Rho(), v_kp.at(1).Py()/v_kp.at(1).Rho(), v_kp.at(1).Pz()/v_kp.at(1).Rho());
    //   old_dir_S2_p.SetXYZ(v_pr.at(0).Px()/v_pr.at(0).Rho(), v_pr.at(0).Py()/v_pr.at(0).Rho(), v_pr.at(0).Pz()/v_pr.at(0).Rho());
    //   old_dir_S2_pi_1.SetXYZ(v_pim.at(0).Px()/v_pim.at(0).Rho(), v_pim.at(0).Py()/v_pim.at(0).Rho(), v_pim.at(0).Pz()/v_pim.at(0).Rho());
    //   old_dir_S2_pi_2.SetXYZ(v_pim.at(1).Px()/v_pim.at(1).Rho(), v_pim.at(1).Py()/v_pim.at(1).Rho(), v_pim.at(1).Pz()/v_pim.at(1).Rho());
    //
    //
    //
    //
    //   doca_S2_kp_1_e_scattered=Calc_dtfInterDOCA(old_dir_S2_kp_1, old_dir_S2_e_scattered, old_vert_S2_kp_1, old_vert_S2_e_scattered, new_vert_S2_kp_1, new_vert_S2_e_scattered);
    //   doca_S2_p_pi_1=Calc_dtfInterDOCA(old_dir_S2_p, old_dir_S2_pi_1, old_vert_S2_p, old_vert_S2_pi_1, new_vert_S2_p, new_vert_S2_pi_1);
    //   doca_S2_p_pi_2=Calc_dtfInterDOCA(old_dir_S2_p, old_dir_S2_pi_2, old_vert_S2_p, old_vert_S2_pi_2, new_vert_S2_p, new_vert_S2_pi_2);
    //   //THE FOLLOWING IS POCA
    //   loc_DOCA_vertex_1.SetX(new_vert_S2_kp_1.X()-(new_vert_S2_kp_1.X()-new_vert_S2_e_scattered.X())/2.0);
    //   loc_DOCA_vertex_1.SetY(new_vert_S2_kp_1.Y()-(new_vert_S2_kp_1.Y()-new_vert_S2_e_scattered.Y())/2.0);
    //   loc_DOCA_vertex_1.SetZ(new_vert_S2_kp_1.Z()-(new_vert_S2_kp_1.Z()-new_vert_S2_e_scattered.Z())/2.0);
    //
    //   loc_DOCA_vertex_2.SetX(new_vert_S2_p.X()-(new_vert_S2_p.X()-new_vert_S2_pi_1.X())/2.0);
    //   loc_DOCA_vertex_2.SetY(new_vert_S2_p.Y()-(new_vert_S2_p.Y()-new_vert_S2_pi_1.Y())/2.0);
    //   loc_DOCA_vertex_2.SetZ(new_vert_S2_p.Z()-(new_vert_S2_p.Z()-new_vert_S2_pi_1.Z())/2.0);
    //
    //   hdocappi->Fill(abs(doca_S2_kp_1_e_scattered));
    //   hdocaeK->Fill(abs(doca_S2_p_pi_1));
    //
    //   hmiss_mass_all_casc->Fill(missallcasc.M2());
    //   hmiss_momentum_all_casc->Fill(missallcasc.Rho());
    //   hmiss_3->Fill(miss3.M());
    //   hcascade->Fill(cascade_lambda_1.M(), cascade_lambda_2.M());
    //   hcascade_lambda_1->Fill(cascade_lambda_1.M());
    //   hcascade_lambda_2->Fill(cascade_lambda_2.M());
    //   hmiss_inv_cascade_1->Fill(cascade.M(),miss3.M());
    //   if(fabs(v_delta_beta_kp.at(0))<0.01 && v_P_kp.at(0)>0.5 && v_P_kp.at(0)<2.6 &&
    //   fabs(v_delta_beta_kp.at(1))<0.01 && v_P_kp.at(1)>0.5 && v_P_kp.at(1)<2.6 &&
    //   fabs(v_delta_beta_pr.at(0))<0.01 && v_P_pr.at(0)>0.5 && v_P_pr.at(0)<4.8){
    //
    //     hmiss_momentum_all_cascc->Fill(missallcasc.Rho());
    //     hmiss_mass_all_cascc->Fill(missallcasc.M2());
    //     hmiss_3c->Fill(miss3.M());
    //     hcascade_c->Fill(cascade_lambda_1.M(), cascade_lambda_2.M());
    //     hcascade_lambda_1c->Fill(cascade_lambda_1.M());
    //     hcascade_lambda_2c->Fill(cascade_lambda_2.M());
    //     hmiss_inv_cascade_1c->Fill(cascade.M(),miss3.M());
    //
    //     if(fabs(missallcasc.M2())<0.05){
    //       hmiss_3t->Fill(miss3.M());
    //       hcascade_t->Fill(cascade_lambda_1.M(), cascade_lambda_2.M());
    //       hcascade_lambda_1t->Fill(cascade_lambda_1.M());
    //       hcascade_lambda_2t->Fill(cascade_lambda_2.M());
    //     }
    //   }
    // }
  }
  cout<<"Good Events: "<<Good_Events<<endl;
  cout<<"Total Events: "<<Total_Events<<endl;
  cout<<"Empty Events: "<<empty<<endl;
  fileOutput1.Write();
}
