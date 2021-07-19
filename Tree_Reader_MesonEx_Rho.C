// electron correction only needed for real data not simulation data
// To change code accordingly, search for "sim data" and "real data"
// To change between mass and PID data, search for "Mass data" and "PID data"


#include <TDatabasePDG.h>
#include <cstdlib>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <math.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TH3.h>
#include <TChain.h>
#include <TCanvas.h>
#include <vector>
#include <TLorentzVector.h>


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

void Tree_Reader_MesonEx_Rho(){
  gROOT->ProcessLine(".L /mnt/f/PhD/Macros/Loader.C+");

  TFile *f = new TFile("/mnt/f/PhD/Trees/MesonEx/PID_Trees/Rho_Simulations_RGA_Fall2018_Inbending_3282_PID_Tree_190721_01.root");
  TTree *t1 = (TTree*)f->Get("Rho_Tree_190721_01");


  Int_t Nparticles=0;


  TLorentzVector *readbeam=NULL;
  TLorentzVector *readtarget=NULL;

  vector<TLorentzVector> *v_p4=0;
  vector<TLorentzVector> *v_vertex=0;
  vector<double> *v_path=0;
  vector<double> *v_time=0;
  vector<double> *v_beta=0;
  vector<double> *v_beta_FT=0;
  vector<double> *v_energy=0;
  vector<Int_t> *v_charge=0;
  vector<Double_t> *v_PID=0;
  vector<double> *v_chi2PID=0;
  vector<double> *v_status=0;
  vector<int> *v_region=0;
  vector<double> *v_mass=0;
  vector<int> *v_Pos_position=0;

  Int_t readchargetracks;
  Int_t readothertracks;
  Int_t readeventno;
  Int_t readrunno;
  Int_t readtriggerno;

  t1->SetBranchAddress("eventno",&readeventno);
  t1->SetBranchAddress("runno",&readrunno);
  t1->SetBranchAddress("triggerno",&readtriggerno);
  // t1->SetBranchAddress("start_time",&start_time);
  t1->SetBranchAddress("p4",&v_p4);
  t1->SetBranchAddress("vertex",&v_vertex);
  t1->SetBranchAddress("beta",&v_beta);
  t1->SetBranchAddress("beta_FT",&v_beta_FT);
  t1->SetBranchAddress("status",&v_status);
  t1->SetBranchAddress("energy",&v_energy);
  t1->SetBranchAddress("charge",&v_charge);
  t1->SetBranchAddress("PID",&v_PID);
  t1->SetBranchAddress("chi2PID",&v_chi2PID);
  t1->SetBranchAddress("time",&v_time);
  t1->SetBranchAddress("region",&v_region);
  t1->SetBranchAddress("path",&v_path);
  t1->SetBranchAddress("beam",&readbeam);
  t1->SetBranchAddress("target",&readtarget);
  /* // Mass data
  t1->SetBranchAddress("mass",&v_mass);
  t1->SetBranchAddress("Pos_position",&v_Pos_position);
  */ // Mass data
  TFile fileOutput1("/mnt/f/PhD/Meson_Analysis/Simulations/Rho_Sim_3282_RGA_Fall2018_Inbending_190721_01.root","recreate");

  //Creating histograms
  // photon polarization
  auto* h_photon_polarization=new TH1D("h_photon_polarization","photon polarization",100,0,1);

  // Tests to check how many particles are in each event
  auto* hkaonpno_1=new TH1F("hkaonpno_1","Number of K^{+} (1);Number of K^{+};Counts",9,0,8);
  auto* hkaonmno_1=new TH1F("hkaonmno_1","Number of K^{-} (1);Number of K^{-};Counts",9,0,8);
  auto* hprotonno_1=new TH1F("hprotonno_1","Number of p (1);Number of p;Counts",9,0,8);
  auto* helectronno_1=new TH1F("helectronno_1","Number of e^{-} (1);Number of e^{-};Counts",9,0,8);

  // Testing the positive masses
  auto* h_mass_positive=new TH2D("h_mass_positive","Positive Masses",200,0,2,200,0,2);
  auto* hbeta_kp=new TH2D("hbeta_kp","Beta measured of K^{+}",200,0,11,200,0,1);
  auto* hbeta_p=new TH2D("hbeta_p","Beta measured of p",200,0,11,200,0,1);

  // Detector distribution
  auto* h_Kp_Detectors=new TH1F("h_Kp_Detectors","K^{+} Detector distribution",3,0,3);
  auto* h_Km_Detectors=new TH1F("h_Km_Detectors","K^{-} Detector distribution",3,0,3);
  auto* h_pr_Detectors=new TH1F("h_pr_Detectors","p Detector distribution",3,0,3);
  auto* h_el_Detectors=new TH1F("h_el_Detectors","e^{-} Detector distribution",3,0,3);


  //PID only cut
  auto* hmiss_mass_all=new TH1F("miss_all","MM^2(e' K^{+} K^{-} p);MM^2(e' K^{+} K^{-} p) [GeV];Counts",200,-2,2);
  auto* hmissp=new TH1F("hmissp","MM(e' K^{+} K^{-});MM(e' K^{+} K^{-}) [GeV];Counts",200,0,2);
  auto* hmissKpKm=new TH2D("hmissKpKm","Missing mass for K^{+} against K^{-};MM(e'K^{-} p) [GeV];MM(e'K^{+} p) [GeV]",200,-1,3,200,-1,3);
  auto* hphoton=new TH1F("hphoton","MM(e');MM(e') [GeV];Counts",200,-1,3);
  auto* hmiss_momentum_all=new TH1F("hmiss_momentum_all","P(B + T - e' - K^{+} - K^{-} - p);P(B + T - e' - K^{+} - K^{-} - p) [GeV];Counts",200,-1,2);
  auto* hinvkpkm=new TH1F("hinvkpkm","M(K^{+} K^{-});M(K^{+} K^{-}) [GeV];Counts",200,0,2);
  auto* hmiss_mass_res=new TH1F("hmiss_mass_res","MM(e' p);MM(e' p) [GeV];Counts",200,0,2);
  auto* hpath_kp=new TH1F("hpath_kp","Path of K^{+};Path of K^{+} ;Counts",500,-1000,1000);
  auto* hvertex_time_kp=new TH1F("hvertex_time_kp","Vertex time K^{+};Vertex time K^{+} (1);Counts",500,-1000,1000);

  // delta beta Histograms
  auto* h_delta_beta_kp=new TH2D("h_delta_beta_kp","#Delta#Beta for K^{+}",500,0,11,400,-1,1);
  auto* h_delta_beta_km=new TH2D("h_delta_beta_km","#Delta#Beta for K^{-}",500,0,11,400,-1,1);
  auto* h_delta_beta_pr=new TH2D("h_delta_beta_pr","#Delta#Beta for p",500,0,11,400,-1,1);


  //delta beta and momentum cuts on kaons and protons
  auto* hinvkpkmc=new TH1F("hinvkpkmc","M(K^{+} K^{-});M(K^{+} K^{-}) [GeV];Counts",200,0,2);
  auto* hmiss_mass_allc=new TH1F("miss_allc","MM^2(e' K^{+} p #pi^{-}) after cuts on delta beta and momentum of K^{+} and p;MM^2(e' K^{+} p #pi^{-}) [GeV];Counts",250,-0.05,0.05);
  auto* hmisspc=new TH1F("hmisspc","MM(e' K^{+} K^{-});MM(e' K^{+} K^{-}) [GeV];Counts",200,-1,3);
  auto* hmiss_momentum_allc=new TH1F("hmiss_momentum_allc","P(B + T - e' - K^{+} - p - #pi^{-}) after cuts on delta beta and momentum of K^{+} and p;P(B + T - e' - K^{+} - p - #pi^{-}) [GeV];Counts",200,-1,2);
  auto* hvertex_time_kpc=new TH1F("hvertex_time_kpc","Vertex time K^{+};Vertex time K^{+} (1);Counts",500,-1000,1000);


  //Exclusivity cuts as well
  auto* hinvkpkmt=new TH1F("hinvkpkmt","M(K^{+} K^{-});M(K^{+} K^{-}) [GeV];Counts",200,0,2);
  auto* hvertex_time_kpt=new TH1F("hvertex_time_kpt","Vertex time K^{+};Vertex time K^{+} (1);Counts",500,-1000,1000);
  auto* hmisspt=new TH1F("hmisspt","MM(e' K^{+} K^{-});MM(e' K^{+} K^{-}) [GeV];Counts",200,-1,3);
  auto* hmissKpKmt=new TH2D("hmissKpKmt","Missing mass for K^{+} against K^{-};MM(e'K^{-} p) [GeV];MM(e'K^{+} p) [GeV]",200,-1,3,200,-1,3);

  //Cut around missing mass of proton instead of exclusivity
  auto* hinvkpkm_p=new TH1F("hinvkpkm_p","M(K^{+} K^{-});M(K^{+} K^{-}) [GeV];Counts",200,0,2);
  auto* hmissp_p=new TH1F("hmissp_p","MM(e' K^{+} K^{-});MM(e' K^{+} K^{-}) [GeV];Counts",200,-1,3);

  // Cut around the missing masses of the K+ and K-
  auto* hinvkpkm_all=new TH1F("hinvkpkm_all","M(K^{+} K^{-});M(K^{+} K^{-}) [GeV];Counts",200,0,2);

  // Photon information
  auto* h_photon_energy=new TH1F("h_photon_energy","#gamma energy;Energy [GeV];Counts",500,0,10);
  auto* h_photon_energy_t=new TH1F("h_photon_energy_t","#gamma energy;Energy [GeV];Counts",500,0,10);
  auto* h_photon_energy_t_boost=new TH1F("h_photon_energy_t_boost","#gamma energy;Energy [GeV];Counts",500,0,10);

  // Theta of invariant KK in CoM frame
  auto* h_invkpkm_CoM_Theta=new TH1F("h_invkpkm_CoM_Theta","#theta of resonance in CoM frame;cos(#theta);Counts",200,-1,1);

  // Angle between boosted kaons and parent particles or photon in invariant mass frame
  auto* h_kp_inv_parent_angle=new TH1F("h_kp_inv_parent_angle","Angle between boosted K^{+} and parent particle;Angle [deg];Counts",180,0,180);
  auto* h_km_inv_parent_angle=new TH1F("h_km_inv_parent_angle","Angle between boosted K^{-} and parent particle;Angle [deg];Counts",180,0,180);
  auto* h_kp_inv_photon_angle=new TH1F("h_kp_inv_photon_angle","Angle between boosted K^{+} and parent particle;Angle [deg];Counts",180,0,180);
  auto* h_km_inv_photon_angle=new TH1F("h_km_inv_photon_angle","Angle between boosted K^{-} and parent particle;Angle [deg];Counts",180,0,180);

  // Angle between boosted kaons and parent particles or photon in missing mass frame
  auto* h_kp_miss_parent_angle=new TH1F("h_kp_miss_parent_angle","Angle between boosted K^{+} and parent particle;Angle [deg];Counts",180,0,180);
  auto* h_km_miss_parent_angle=new TH1F("h_km_miss_parent_angle","Angle between boosted K^{-} and parent particle;Angle [deg];Counts",180,0,180);
  auto* h_kp_miss_photon_angle=new TH1F("h_kp_miss_photon_angle","Angle between boosted K^{+} and parent particle;Angle [deg];Counts",180,0,180);
  auto* h_km_miss_photon_angle=new TH1F("h_km_miss_photon_angle","Angle between boosted K^{-} and parent particle;Angle [deg];Counts",180,0,180);
  auto* h_kp_miss_phi=new TH1F("h_kp_miss_phi","Boosted K^{+} #Phi in missing mass frame;Angle [deg];Counts",180,0,180);
  auto* h_km_miss_phi=new TH1F("h_km_miss_phi","Boosted K^{-} #Phi in missing mass frame;Angle [deg];Counts",180,0,180);

  // Looking at difference between measured and reconstructed electron
  auto* hdelta_P=new TH1F("hdelta_P","",200,-10,10);
  auto* hdelta_Theta=new TH1F("hdelta_Theta","",100,-8,8);
  auto* hdelta_Phi=new TH1F("hdelta_Phi","",100,-8,8);
  auto* hdelta_Phi_2=new TH1F("hdelta_Phi_2","",100,-8,8);

  auto* h_phi=new TH1F("h_phi","phi",120,-180,180);
  auto* h_phi1=new TH1F("h_phi1","phi",120,-180,180);
  auto* h_phi2=new TH1F("h_phi2","phi",120,-180,180);
  auto* h_phi3=new TH1F("h_phi3","phi",120,-180,180);
  auto* h_phi4=new TH1F("h_phi4","phi",120,-180,180);


  auto* h_region=new TH3D("h_region","region",3,0,3,3,0,3,3,0,3);

  auto* h_t=new TH1D("h_t","t distribution in lab frame",100,0,10);
  auto* h_t_boost=new TH1D("h_t_boost","t distribution in rho rest frame",100,0,300);

  vector<TLorentzVector> v_el;
  vector<TLorentzVector> v_pr;
  vector<TLorentzVector> v_kp;
  vector<TLorentzVector> v_km, v_pip, v_pim;
  vector<TLorentzVector> v_vertex_el;
  vector<TLorentzVector> v_vertex_pr;
  vector<TLorentzVector> v_vertex_kp;
  vector<TLorentzVector> v_vertex_km,v_vertex_pim, v_vertex_pip;

  TLorentzVector beam, beam_boost_inv;
  TLorentzVector missall;
  TLorentzVector missp;
  TLorentzVector missKp;
  TLorentzVector missKm;
  TLorentzVector pip_boost_rho;
  TLorentzVector misse;
  TLorentzVector CoM_Frame;
  TVector3 X_Axis,Y_Axis;
  Double_t Phi;
  Double_t pip_x_component, pip_y_component;
  TLorentzVector invkpkm, invkpkm_boost_CoM, invrho;
  TLorentzVector miss_mass_res;
  TLorentzVector photon, photon_boost_inv, photon_boost_miss;
  TVector3 inv_KpKm_boost, miss_mass_res_3, photon_3, photon_3_boost_inv, photon_3_boost_miss, beam_3;
  TVector3 el_corrected_3, el_3;
  TLorentzVector vertex_el, vertex_pip, vertex_pim;
  TLorentzVector vertex_pr;
  TLorentzVector vertex_kp;
  TLorentzVector vertex_km;
  Double_t v, Q_2;  // properties of photon for polarization calculations
  Double_t el_theta; // scattering angle of electron for polarization calculations
  Double_t Photon_polarization; // calculated photon polarization

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
  vector<Double_t> v_region_el;
  Double_t region_el;

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
  vector<Double_t> v_region_pip;
  Double_t region_pip;

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
  vector<Double_t> v_region_pim;
  Double_t region_pim;

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
  vector<Double_t> v_region_pr;
  Double_t region_pr;

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
  vector<Double_t> v_region_kp;
  Double_t region_kp;

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
  vector<Double_t> v_region_km;
  Double_t region_km;

  Double_t c=30;
  Double_t t;
  TLorentzVector momentum_exchange; // momentum transfer

  TLorentzVector el, pip, pim;
  TLorentzVector el_corrected, el_corrected_boost_inv;
  TLorentzVector pr;
  TLorentzVector kp, kp_boost_inv, kp_boost_miss;
  TLorentzVector km, km_boost_inv, km_boost_miss;


  double doca_S2_kp_1_e_scattered, doca_S2_p_pi_1, doca_S2_p_pi_2;
  TVector3 old_vert_S2_kp_1, old_vert_S2_kp_2, old_vert_S2_e_scattered, old_vert_S2_p, old_vert_S2_pi_1,old_vert_S2_pi_2, old_dir_S2_kp_1, old_dir_S2_kp_2, old_dir_S2_e_scattered, old_dir_S2_p, old_dir_S2_pi_1,old_dir_S2_pi_2, new_vert_S2_kp_1, new_vert_S2_e_scattered, new_vert_S2_p, new_vert_S2_pi_1, new_vert_S2_pi_2;
  TVector3 loc_DOCA_vertex_1, loc_DOCA_vertex_2;
  new_vert_S2_kp_1.SetXYZ(0.0,0.0,0.0);
  new_vert_S2_e_scattered.SetXYZ(0.0,0.0,0.0);
  new_vert_S2_p.SetXYZ(0.0,0.0,0.0);
  new_vert_S2_pi_1.SetXYZ(0.0,0.0,0.0);
  new_vert_S2_pi_2.SetXYZ(0.0,0.0,0.0);

  Int_t empty = 0;
  Long64_t nentries = t1->GetEntries();
  // Long64_t nentries = 5000000;
  cout<<nentries<<endl;
  Int_t Percentage = nentries/100;
  for(Long64_t i=0; i<nentries;i++){
    t1->GetEntry(i);
    if (i % Percentage == 0){
      fprintf (stderr, "%lld\r", i/Percentage);
      fflush (stderr);
    }
    v_pip.clear();
    v_beta_tof_pip.clear();
    v_P_pip.clear();
    v_path_pip.clear();
    v_TOF_pip.clear();
    v_beta_calc_pip.clear();
    v_delta_beta_pip.clear();
    v_vertex_time_pip.clear();
    v_vertex_pip.clear();
    v_region_pip.clear();

    v_pim.clear();
    v_beta_tof_pim.clear();
    v_P_pim.clear();
    v_path_pim.clear();
    v_TOF_pim.clear();
    v_beta_calc_pim.clear();
    v_delta_beta_pim.clear();
    v_vertex_time_pim.clear();
    v_vertex_pim.clear();
    v_region_pim.clear();

    v_el.clear();
    v_beta_tof_el.clear();
    v_P_el.clear();
    v_path_el.clear();
    v_TOF_el.clear();
    v_beta_calc_el.clear();
    v_delta_beta_el.clear();
    v_vertex_time_el.clear();
    v_vertex_el.clear();
    v_region_el.clear();

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
    v_region_pr.clear();


    v_pip.clear();
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
    v_region_kp.clear();

    v_pim.clear();
    v_beta_tof_km.clear();
    v_P_km.clear();
    v_path_km.clear();
    v_TOF_km.clear();
    v_beta_calc_km.clear();
    v_delta_beta_km.clear();
    v_vertex_km.clear();
    v_vertex_time_km.clear();
    v_vertex_km.clear();
    v_region_km.clear();

    // Get the number of particles in each event
    if(v_p4->size()) Nparticles = v_p4->size();
    // Loop over all the particles in the event
    for(Int_t j=0; j<Nparticles; j++){

      /* // Mass data
      if(v_charge->at(j) < 0){

        if(v_PID->at(j)==11){

          el.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),0.000511);
          P_el = sqrt((pow(v_p4->at(j).Px(),2))+(pow(v_p4->at(j).Py(),2))+(pow(v_p4->at(j).Pz(),2)));
          TOF_el = v_time->at(j);
          path_el = v_path->at(j);
          beta_tof_el = v_beta->at(j);
          beta_calc_el = P_el/(sqrt((pow(P_el,2))+(pow(el.M(),2))));
          delta_beta_el = beta_calc_el-beta_tof_el;
          vertex_time_el = TOF_el - path_el / (beta_tof_el*c);
          vertex_el.SetXYZT(v_vertex->at(j).X(), v_vertex->at(j).Y(), v_vertex->at(j).Z(), vertex_time_el);
          region_el = v_region->at(j);
          h_el_Detectors->Fill(region_el);

          v_el.push_back(el);
          v_beta_tof_el.push_back(beta_tof_el);
          v_P_el.push_back(P_el);
          v_path_el.push_back(path_el);
          v_TOF_el.push_back(TOF_el);
          v_beta_calc_el.push_back(beta_calc_el);
          v_delta_beta_el.push_back(delta_beta_el);
          v_vertex_time_el.push_back(vertex_time_el);
          v_vertex_el.push_back(vertex_el);
          v_region_el.push_back(region_el);
        }
        else {

          km.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),0.493677);
          P_km = sqrt((pow(v_p4->at(j).Px(),2))+(pow(v_p4->at(j).Py(),2))+(pow(v_p4->at(j).Pz(),2)));
          TOF_km = v_time->at(j);
          path_km = v_path->at(j);
          beta_tof_km = v_beta->at(j);
          beta_calc_km = P_km/(sqrt((pow(P_km,2))+(pow(km.M(),2))));
          delta_beta_km = beta_calc_km-beta_tof_km;
          vertex_time_km = TOF_km - path_km / (beta_tof_km*c);
          vertex_km.SetXYZT(v_vertex->at(j).X(), v_vertex->at(j).Y(), v_vertex->at(j).Z(), vertex_time_km);
          region_km = v_region->at(j);
          h_Km_Detectors->Fill(region_km);

          v_pim.push_back(km);
          v_beta_tof_km.push_back(beta_tof_km);
          v_P_km.push_back(P_km);
          v_path_km.push_back(path_km);
          v_TOF_km.push_back(TOF_km);
          v_beta_calc_km.push_back(beta_calc_km);
          v_delta_beta_km.push_back(delta_beta_km);
          v_vertex_time_km.push_back(vertex_time_km);
          v_vertex_km.push_back(vertex_km);
          v_region_km.push_back(region_km);
        }
      }

      // Using pure MC data
      if(v_mcPID->at(j)==11){
        mc_el.SetXYZM(v_mcp4->at(j).Px(),v_mcp4->at(j).Py(),v_mcp4->at(j).Pz(),0.000511);
      }
      else if(v_mcPID->at(j)==321){
        mc_kp.SetXYZM(v_mcp4->at(j).Px(),v_mcp4->at(j).Py(),v_mcp4->at(j).Pz(),0.493677);
      }
      else if(v_mcPID->at(j)==-321){
        mc_km.SetXYZM(v_mcp4->at(j).Px(),v_mcp4->at(j).Py(),v_mcp4->at(j).Pz(),0.493677);
      }
      else if(v_mcPID->at(j)==2212){
        mc_pr.SetXYZM(v_mcp4->at(j).Px(),v_mcp4->at(j).Py(),v_mcp4->at(j).Pz(),0.938272);
      }
    }
    // Checking the masses of positive particles to assign as proton or kaon
    if(v_mass->at(v_Pos_position->at(0)) > v_mass->at(v_Pos_position->at(1))){

      // Setting values for the proton
      pr.SetXYZM(v_p4->at(v_Pos_position->at(0)).Px(),v_p4->at(v_Pos_position->at(0)).Py(),v_p4->at(v_Pos_position->at(0)).Pz(),0.938272);
      P_pr = sqrt((pow(v_p4->at(v_Pos_position->at(0)).Px(),2))+(pow(v_p4->at(v_Pos_position->at(0)).Py(),2))+(pow(v_p4->at(v_Pos_position->at(0)).Pz(),2)));
      TOF_pr = v_time->at(v_Pos_position->at(0));
      path_pr = v_path->at(v_Pos_position->at(0));
      beta_tof_pr = v_beta->at(v_Pos_position->at(0));
      beta_calc_pr = P_pr/(sqrt((pow(P_pr,2))+(pow(pr.M(),2))));
      delta_beta_pr = beta_calc_pr-beta_tof_pr;
      vertex_time_pr = TOF_pr - path_pr / (beta_tof_pr*c);
      vertex_pr.SetXYZT(v_vertex->at(v_Pos_position->at(0)).X(), v_vertex->at(v_Pos_position->at(0)).Y(), v_vertex->at(v_Pos_position->at(0)).Z(), vertex_time_pr);
      region_pr = v_region->at(v_Pos_position->at(0));
      h_pr_Detectors->Fill(region_pr);

      v_pr.push_back(pr);
      v_beta_tof_pr.push_back(beta_tof_pr);
      v_P_pr.push_back(P_pr);
      v_path_pr.push_back(path_pr);
      v_TOF_pr.push_back(TOF_pr);
      v_beta_calc_pr.push_back(beta_calc_pr);
      v_delta_beta_pr.push_back(delta_beta_pr);
      v_vertex_time_pr.push_back(vertex_time_pr);
      v_vertex_pr.push_back(vertex_pr);
      v_region_pr.push_back(region_pr);

      // Setting values for the K+
      kp.SetXYZM(v_p4->at(v_Pos_position->at(1)).Px(),v_p4->at(v_Pos_position->at(1)).Py(),v_p4->at(v_Pos_position->at(1)).Pz(),0.493677);
      P_kp = sqrt((pow(v_p4->at(v_Pos_position->at(1)).Px(),2))+(pow(v_p4->at(v_Pos_position->at(1)).Py(),2))+(pow(v_p4->at(v_Pos_position->at(1)).Pz(),2)));
      TOF_kp = v_time->at(v_Pos_position->at(1));
      path_kp = v_path->at(v_Pos_position->at(1));
      beta_tof_kp = v_beta->at(v_Pos_position->at(1));
      beta_calc_kp = P_kp/(sqrt((pow(P_kp,2))+(pow(kp.M(),2))));
      delta_beta_kp = beta_calc_kp-beta_tof_kp;
      vertex_time_kp = TOF_kp - path_kp / (beta_tof_kp*c);
      vertex_kp.SetXYZT(v_vertex->at(v_Pos_position->at(1)).X(), v_vertex->at(v_Pos_position->at(1)).Y(), v_vertex->at(v_Pos_position->at(1)).Z(), vertex_time_kp);
      region_kp = v_region->at(v_Pos_position->at(1));
      h_Kp_Detectors->Fill(region_kp);

      v_kp.push_back(kp);
      v_beta_tof_kp.push_back(beta_tof_kp);
      v_P_kp.push_back(P_kp);
      v_path_kp.push_back(path_kp);
      v_TOF_kp.push_back(TOF_kp);
      v_beta_calc_kp.push_back(beta_calc_kp);
      v_delta_beta_kp.push_back(delta_beta_kp);
      v_vertex_time_kp.push_back(vertex_time_kp);
      v_vertex_kp.push_back(vertex_kp);
      v_region_kp.push_back(region_kp);
    }

    else if(v_mass->at(v_Pos_position->at(0)) < v_mass->at(v_Pos_position->at(1))) {

      // Setting values for the proton
      pr.SetXYZM(v_p4->at(v_Pos_position->at(1)).Px(),v_p4->at(v_Pos_position->at(1)).Py(),v_p4->at(v_Pos_position->at(1)).Pz(),0.938272);
      P_pr = sqrt((pow(v_p4->at(v_Pos_position->at(1)).Px(),2))+(pow(v_p4->at(v_Pos_position->at(1)).Py(),2))+(pow(v_p4->at(v_Pos_position->at(1)).Pz(),2)));
      TOF_pr = v_time->at(v_Pos_position->at(1));
      path_pr = v_path->at(v_Pos_position->at(1));
      beta_tof_pr = v_beta->at(v_Pos_position->at(1));
      beta_calc_pr = P_pr/(sqrt((pow(P_pr,2))+(pow(pr.M(),2))));
      delta_beta_pr = beta_calc_pr-beta_tof_pr;
      vertex_time_pr = TOF_pr - path_pr / (beta_tof_pr*c);
      vertex_pr.SetXYZT(v_vertex->at(v_Pos_position->at(1)).X(), v_vertex->at(v_Pos_position->at(1)).Y(), v_vertex->at(v_Pos_position->at(1)).Z(), vertex_time_pr);
      region_pr = v_region->at(v_Pos_position->at(1));
      h_pr_Detectors->Fill(region_pr);

      v_pr.push_back(pr);
      v_beta_tof_pr.push_back(beta_tof_pr);
      v_P_pr.push_back(P_pr);
      v_path_pr.push_back(path_pr);
      v_TOF_pr.push_back(TOF_pr);
      v_beta_calc_pr.push_back(beta_calc_pr);
      v_delta_beta_pr.push_back(delta_beta_pr);
      v_vertex_time_pr.push_back(vertex_time_pr);
      v_vertex_pr.push_back(vertex_pr);
      v_region_pr.push_back(region_pr);

      // Setting values for the K+
      kp.SetXYZM(v_p4->at(v_Pos_position->at(0)).Px(),v_p4->at(v_Pos_position->at(0)).Py(),v_p4->at(v_Pos_position->at(0)).Pz(),0.493677);
      P_kp = sqrt((pow(v_p4->at(v_Pos_position->at(0)).Px(),2))+(pow(v_p4->at(v_Pos_position->at(0)).Py(),2))+(pow(v_p4->at(v_Pos_position->at(0)).Pz(),2)));
      TOF_kp = v_time->at(v_Pos_position->at(0));
      path_kp = v_path->at(v_Pos_position->at(0));
      beta_tof_kp = v_beta->at(v_Pos_position->at(0));
      beta_calc_kp = P_kp/(sqrt((pow(P_kp,2))+(pow(kp.M(),2))));
      delta_beta_kp = beta_calc_kp-beta_tof_kp;
      vertex_time_kp = TOF_kp - path_kp / (beta_tof_kp*c);
      vertex_kp.SetXYZT(v_vertex->at(v_Pos_position->at(0)).X(), v_vertex->at(v_Pos_position->at(0)).Y(), v_vertex->at(v_Pos_position->at(0)).Z(), vertex_time_kp);
      region_kp = v_region->at(v_Pos_position->at(0));
      h_Kp_Detectors->Fill(region_kp);

      v_kp.push_back(kp);
      v_beta_tof_kp.push_back(beta_tof_kp);
      v_P_kp.push_back(P_kp);
      v_path_kp.push_back(path_kp);
      v_TOF_kp.push_back(TOF_kp);
      v_beta_calc_kp.push_back(beta_calc_kp);
      v_delta_beta_kp.push_back(delta_beta_kp);
      v_vertex_time_kp.push_back(vertex_time_kp);
      v_vertex_kp.push_back(vertex_kp);
      v_region_kp.push_back(region_kp);
      */ // Mass data

      // /* // PID data
      if(v_PID->at(j)==11){
      el.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),0.000511);
      P_el = el.Rho();
      TOF_el = v_time->at(j);
      path_el = v_path->at(j);
      beta_tof_el = v_beta_FT->at(j);
      beta_calc_el = P_el/(sqrt((pow(P_el,2))+(pow(el.M(),2))));
      delta_beta_el = beta_calc_el-beta_tof_el;
      vertex_time_el = TOF_el - path_el / (beta_tof_el*c);
      vertex_el.SetXYZT(v_vertex->at(j).X(), v_vertex->at(j).Y(), v_vertex->at(j).Z(), vertex_time_el);
      region_el = v_region->at(j);
      h_el_Detectors->Fill(region_el);

      v_el.push_back(el);
      v_beta_tof_el.push_back(beta_tof_el);
      v_P_el.push_back(P_el);
      v_path_el.push_back(path_el);
      v_TOF_el.push_back(TOF_el);
      v_beta_calc_el.push_back(beta_calc_el);
      v_delta_beta_el.push_back(delta_beta_el);
      v_vertex_time_el.push_back(vertex_time_el);
      v_vertex_el.push_back(vertex_el);
      v_region_el.push_back(region_el);
    }
    else if(v_PID->at(j)==-211){
    pim.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),0.139570);
    P_pim = pim.Rho();
    TOF_pim = v_time->at(j);
    path_pim = v_path->at(j);
    beta_tof_pim = v_beta_FT->at(j);
    beta_calc_pim = P_pim/(sqrt((pow(P_pim,2))+(pow(pim.M(),2))));
    delta_beta_pim = beta_calc_pim-beta_tof_pim;
    vertex_time_pim = TOF_pim - path_pim / (beta_tof_pim*c);
    vertex_pim.SetXYZT(v_vertex->at(j).X(), v_vertex->at(j).Y(), v_vertex->at(j).Z(), vertex_time_pim);
    region_pim = v_region->at(j);

    v_pim.push_back(pim);
    v_beta_tof_pim.push_back(beta_tof_pim);
    v_P_pim.push_back(P_pim);
    v_path_pim.push_back(path_pim);
    v_TOF_pim.push_back(TOF_pim);
    v_beta_calc_pim.push_back(beta_calc_pim);
    v_delta_beta_pim.push_back(delta_beta_pim);
    v_vertex_time_pim.push_back(vertex_time_pim);
    v_vertex_pim.push_back(vertex_pim);
    v_region_pim.push_back(region_pim);
  }
  else if(v_PID->at(j)==211){
  pip.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),0.139570);
  P_pip = pip.Rho();
  TOF_pip = v_time->at(j);
  path_pip = v_path->at(j);
  beta_tof_pip = v_beta_FT->at(j);
  beta_calc_pip = P_pip/(sqrt((pow(P_pip,2))+(pow(pip.M(),2))));
  delta_beta_pip = beta_calc_pip-beta_tof_pip;
  vertex_time_pip = TOF_pip - path_pip / (beta_tof_pip*c);
  vertex_pip.SetXYZT(v_vertex->at(j).X(), v_vertex->at(j).Y(), v_vertex->at(j).Z(), vertex_time_pip);
  region_pip = v_region->at(j);

  v_pip.push_back(pip);
  v_beta_tof_pip.push_back(beta_tof_pip);
  v_P_pip.push_back(P_pip);
  v_path_pip.push_back(path_pip);
  v_TOF_pip.push_back(TOF_pip);
  v_beta_calc_pip.push_back(beta_calc_pip);
  v_delta_beta_pip.push_back(delta_beta_pip);
  v_vertex_time_pip.push_back(vertex_time_pip);
  v_vertex_pip.push_back(vertex_pip);
  v_region_pip.push_back(region_pip);
}
    else if(v_PID->at(j)==-321){
    km.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),0.49368);
    P_km = sqrt((pow(v_p4->at(j).Px(),2))+(pow(v_p4->at(j).Py(),2))+(pow(v_p4->at(j).Pz(),2)));
    TOF_km = v_time->at(j);
    path_km = v_path->at(j);
    beta_tof_km = v_beta_FT->at(j);
    beta_calc_km = P_km/(sqrt((pow(P_km,2))+(pow(km.M(),2))));
    delta_beta_km = beta_calc_km-beta_tof_km;
    vertex_time_km = TOF_km - path_km / (beta_tof_km*c);
    vertex_km.SetXYZT(v_vertex->at(j).X(), v_vertex->at(j).Y(), v_vertex->at(j).Z(), vertex_time_km);
    region_km = v_region->at(j);
    h_Km_Detectors->Fill(region_km);

    v_pim.push_back(km);
    v_beta_tof_km.push_back(beta_tof_km);
    v_P_km.push_back(P_km);
    v_path_km.push_back(path_km);
    v_TOF_km.push_back(TOF_km);
    v_beta_calc_km.push_back(beta_calc_km);
    v_delta_beta_km.push_back(delta_beta_km);
    v_vertex_time_km.push_back(vertex_time_km);
    v_vertex_km.push_back(vertex_km);
    v_region_km.push_back(region_km);
  }

  else if(v_PID->at(j)==2212){
  pr.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),0.938);
  P_pr = pr.Rho();
  TOF_pr = v_time->at(j);
  path_pr = v_path->at(j);
  beta_tof_pr = v_beta_FT->at(j);
  beta_calc_pr = P_pr/(sqrt((pow(P_pr,2))+(pow(pr.M(),2))));
  delta_beta_pr = beta_calc_pr-beta_tof_pr;
  vertex_time_pr = TOF_pr - path_pr / (beta_tof_pr*c);
  vertex_pr.SetXYZT(v_vertex->at(j).X(), v_vertex->at(j).Y(), v_vertex->at(j).Z(), vertex_time_pr);
  region_pr = v_region->at(j);
  h_pr_Detectors->Fill(region_pr);

  v_pr.push_back(pr);
  v_beta_tof_pr.push_back(beta_tof_pr);
  v_P_pr.push_back(P_pr);
  v_path_pr.push_back(path_pr);
  v_TOF_pr.push_back(TOF_pr);
  v_beta_calc_pr.push_back(beta_calc_pr);
  v_delta_beta_pr.push_back(delta_beta_pr);
  v_vertex_time_pr.push_back(vertex_time_pr);
  v_vertex_pr.push_back(vertex_pr);
  v_region_pr.push_back(region_pr);
}
else if(v_PID->at(j)==321){
kp.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),0.49368);
P_kp = sqrt((pow(v_p4->at(j).Px(),2))+(pow(v_p4->at(j).Py(),2))+(pow(v_p4->at(j).Pz(),2)));
TOF_kp = v_time->at(j);
path_kp = v_path->at(j);
beta_tof_kp = v_beta_FT->at(j);
beta_calc_kp = P_kp/(sqrt((pow(P_kp,2))+(pow(kp.M(),2))));
delta_beta_kp = beta_calc_kp-beta_tof_kp;
vertex_time_kp = TOF_kp - path_kp / (beta_tof_kp*c);
vertex_kp.SetXYZT(v_vertex->at(j).X(), v_vertex->at(j).Y(), v_vertex->at(j).Z(), vertex_time_kp);
region_kp = v_region->at(j);
h_Kp_Detectors->Fill(region_kp);

v_kp.push_back(kp);
v_beta_tof_kp.push_back(beta_tof_kp);
v_P_kp.push_back(P_kp);
v_path_kp.push_back(path_kp);
v_TOF_kp.push_back(TOF_kp);
v_beta_calc_kp.push_back(beta_calc_kp);
v_delta_beta_kp.push_back(delta_beta_kp);
v_vertex_time_kp.push_back(vertex_time_kp);
v_vertex_kp.push_back(vertex_kp);
v_region_kp.push_back(region_kp);
}
// */ //PID data
}


if(v_el.size()==0)continue;
el_corrected = Correct_Electron_g(v_el.at(0));

// Selecting Events
    if(v_pip.size()==1 && v_pim.size()==1 && v_el.size()==1 && v_pr.size()==1){

      // if(v_region_pip.at(0)==1 && v_region_pim.at(0)==1 && v_region_pr.at(0)==1){

      // Use the non-corrected electron for sim data
      missall = (TLorentzVector)*readbeam + (TLorentzVector)*readtarget - v_el.at(0) - v_pip.at(0) - v_pr.at(0) - v_pim.at(0);
      missp = (TLorentzVector)*readbeam + (TLorentzVector)*readtarget - v_el.at(0) - v_pip.at(0) - v_pim.at(0);
      missKp = (TLorentzVector)*readbeam + (TLorentzVector)*readtarget - v_el.at(0) - v_pr.at(0) - v_pim.at(0);
      missKm = (TLorentzVector)*readbeam + (TLorentzVector)*readtarget - v_el.at(0) - v_pip.at(0) - v_pr.at(0);


      // Use the corrected electron for real data
      // missall = (TLorentzVector)*readbeam + (TLorentzVector)*readtarget - el_corrected - v_pip.at(0) - v_pr.at(0) - v_pim.at(0);
      // missp = (TLorentzVector)*readbeam + (TLorentzVector)*readtarget - el_corrected - v_pip.at(0) - v_pim.at(0);
      // missKp = (TLorentzVector)*readbeam + (TLorentzVector)*readtarget - el_corrected - v_pr.at(0) - v_pim.at(0);
      // missKm = (TLorentzVector)*readbeam + (TLorentzVector)*readtarget - el_corrected - v_pip.at(0) - v_pr.at(0);
      // misse = (TLorentzVector)*readbeam + (TLorentzVector)*readtarget - v_pip.at(0) - v_pr.at(0) - v_pim.at(0);

    beam = (TLorentzVector)*readbeam;
    beam_3 = beam.Vect();
    // Calculate momentum change between proton and target
    momentum_exchange = v_pr.at(0) - (TLorentzVector)*readtarget;
    t = -(momentum_exchange.M2());



    // defining the invariant and missing mass of the resonance
    invrho = v_pip.at(0) + v_pim.at(0);
    miss_mass_res = (TLorentzVector)*readbeam + (TLorentzVector)*readtarget - v_el.at(0) - v_pr.at(0);

    // Getting information on the virtual photon
    // photon = (TLorentzVector)*readbeam - el_corrected; // Use this for real data
    photon = (TLorentzVector)*readbeam - v_el.at(0); // Use this for sim data

    v = photon.E(); // photon energy
    Q_2 = -(photon.M2()); // photon 4 vector squared
    el_theta = el_corrected.Theta(); // scattered electron theta

    // Calculating photon polarization
    Photon_polarization = 1 / (1 + (2 * ((Q_2 + pow(v,2)) / Q_2)) * pow(tan(el_theta / 2),2));

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Boosting all particles into rho rest frame

    inv_KpKm_boost = invrho.BoostVector(); // Get boost vector of rho
    // inv_KpKm_boost = invrho.BoostVector(); // Get boost vector of rho

    //boost the pion in rho rest frame
    pip_boost_rho = v_pip.at(0); // Get TLorentzVector of pi+
    pip_boost_rho.Boost(-inv_KpKm_boost); // Boost pi+ in rho rest frame

    //boost the photon in rho rest frame
    photon_boost_inv = photon;
    photon_boost_inv.Boost(-inv_KpKm_boost); // Boost photon in rho rest frame

    //boost the electron in rho rest frame
    el_corrected_boost_inv = el_corrected;
    el_corrected_boost_inv.Boost(-inv_KpKm_boost); // Boost electron in rho rest frame

    //boost the beam in rho rest frame
    beam_boost_inv = beam;
    beam_boost_inv.Boost(-inv_KpKm_boost); // Boost electron in rho rest frame


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Looking at difference between measured and reconstructed electron
    hdelta_P->Fill(el_corrected.Rho() - misse.Rho());
    hdelta_Theta->Fill(el_corrected.Theta() - misse.Theta());
    hdelta_Phi->Fill(el_corrected.Phi() - misse.Phi());


    // Plotting the missing and invariant mass plots
    hmiss_mass_all->Fill(missall.M2());
    hmiss_momentum_all->Fill(missall.Rho());
    hinvkpkm->Fill(invrho.M());
    hmiss_mass_res->Fill(miss_mass_res.M());
    hmissp->Fill(missp.M());

    // Plotting the delta beta plots
    h_delta_beta_kp->Fill(v_P_pip.at(0),v_delta_beta_pip.at(0));
    h_delta_beta_km->Fill(v_P_pim.at(0),v_delta_beta_pim.at(0));
    h_delta_beta_pr->Fill(v_P_pr.at(0),v_delta_beta_pr.at(0));


    // Getting the 3 vector for scattered electron
    el_corrected_3.SetXYZ(el_corrected.Px()/el_corrected.E(),el_corrected.Py()/el_corrected.E(),el_corrected.Pz()/el_corrected.E());

    // Getting the photon 3 vector
    photon_3.SetXYZ(photon.Px() / photon.E(),photon.Py() / photon.E(), photon.Pz() / photon.E());


    if(fabs(v_delta_beta_pip.at(0))<0.02 && /* v_P_pip.at(0)>0.5 && v_P_kp.at(0)<2.6 &&*/ fabs(v_delta_beta_pim.at(0))<0.02 /*&& v_P_pim.at(0)>0.5 && v_P_km.at(0)<2.6*/ && fabs(v_delta_beta_pr.at(0))<0.02 /*&& v_P_pr.at(0)>0.2 && v_P_pr.at(0)<4.8*/){

      hmiss_mass_allc->Fill(missall.M2());
      hmiss_momentum_allc->Fill(missall.Rho());
      hinvkpkmc->Fill(invrho.M());
      hmisspc->Fill(missp.M());

      // if(fabs(missall.M2())<0.01 && missall.Rho()>0.3)hinvkpkmt->Fill(invkpkm.M());
      if(fabs(missall.M2())<0.01 && missall.Rho()<0.5){
        hinvkpkmt->Fill(invrho.M());
        hmisspt->Fill(missp.M());
        h_photon_energy_t->Fill(photon.E());
        hmissKpKmt->Fill(missKp.M(),missKm.M());
        h_t->Fill(t);


        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Determining phi

        // Defining reference frame
        X_Axis = (el_corrected_boost_inv.Vect()).Cross(beam_boost_inv.Vect()); // scattered e cross with beam
        Y_Axis = (photon_boost_inv.Vect()).Cross(X_Axis); // z-axis cross with x-axis to get y-axis

        // Find xy-components of pi+
        pip_x_component = pip_boost_rho.Vect() * X_Axis;
        pip_y_component = pip_boost_rho.Vect() * Y_Axis;

        // Calculate phi from the xy-components of pi+
        Phi = atan2(pip_x_component,pip_y_component);

        // Checking the sign of the phi angle
        // if(X_Axis.Angle(Y_Prime_Axis)*TMath::RadToDeg()>90) Phi = -Phi;

        // Applying cut on t and photon energy
        // if(photon.E() < 2.4 || photon.E() > 2.6)continue;
        // if(t < 0.1 || t > 0.2)continue;

        h_photon_polarization->Fill(Photon_polarization);

        // if(Photon_polarization < 0.85 || Photon_polarization > 0.95)continue;

        // Plotting the phi angles
        h_phi->Fill(TMath::RadToDeg()*Phi);

        if(invrho.M() > 0.67 && invrho.M() < 0.87) h_phi1->Fill(TMath::RadToDeg()*Phi);
        if(invrho.M() > 1.1 && invrho.M() < 1.5) h_phi2->Fill(TMath::RadToDeg()*Phi);
        if(invrho.M() > 1.5 && invrho.M() < 1.9) h_phi3->Fill(TMath::RadToDeg()*Phi);
        if(invrho.M() > 1.9) h_phi4->Fill(TMath::RadToDeg()*Phi);

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      }
    // }
  }
}
}



cout<<"Total Events: "<<nentries<<endl;
fileOutput1.Write();
}
