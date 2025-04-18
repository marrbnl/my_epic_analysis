//================================================
// Study the kinematics of produced hadron in 
// certain (x, Q2, nu) bin
//
// To-do:
// 1. Exclude/include feed-down?
//================================================
R__LOAD_LIBRARY(libeicsmear);

#include "bins.h"

const double degree = 180./TMath::Pi();

const int verbosity = 1; // use this to control debugging messages

using namespace std;

class hadron_gen
{
  private:
    // evt info
    int hadron_id;

    double x_true;
    double Q2_true;
    double nu_true;

    TLorentzVector hadron_beam;
    TLorentzVector struck_quark;

    TH1D* h1d_parent_id;

    int iQ2bin;
    int ixbin;
    int inubin;

    // kinematics of produced hadron in Q2 and x bin
    TH2D* h2d_hadron_pt_vs_eta_gen_in_Q2_x[Q2bin][xbin];
    TH2D* h2d_hadron_z_vs_eta_gen_in_Q2_x[Q2bin][xbin];
    
    // kinematics of produced hadron in Q2 and nu bin
    TH2D* h2d_hadron_pt_vs_eta_gen_in_Q2_nu[Q2bin][nubin];
    TH2D* h2d_hadron_z_vs_eta_gen_in_Q2_nu[Q2bin][nubin];

    // kinematics of produced hadron in x and nu bin
    TH2D* h2d_hadron_pt_vs_eta_gen_in_x_nu[xbin][nubin];
    TH2D* h2d_hadron_z_vs_eta_gen_in_x_nu[xbin][nubin];

  public:
    hadron_gen(int _hadron_id)
    {
      cout << "Constructing analyzing module to study hadronization of particle with ID " << _hadron_id << endl;
      hadron_id = abs(_hadron_id);

      x_true = 1E1; // unphysical
      Q2_true = 1E-5; // out of range
      nu_true = -9999;

      iQ2bin = -9999;
      ixbin = -9999;
      inubin = -9999;

      h1d_parent_id = new TH1D(Form("h1d_hadron_%d_parent_id",hadron_id),"parent ID",10000,-0.5,9999.5);
      h1d_parent_id->Sumw2();

      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        for (int ix = 0; ix < xbin; ++ix)
        {
          h2d_hadron_pt_vs_eta_gen_in_Q2_x[iQ2][ix] = new TH2D(Form("h2d_hadron_%d_pt_vs_eta_gen_in_Q2_x_%d_%d",hadron_id,iQ2,ix),"D0 pt vs eta",100,0,10,40,-10,10);
          h2d_hadron_pt_vs_eta_gen_in_Q2_x[iQ2][ix]->Sumw2();

          h2d_hadron_z_vs_eta_gen_in_Q2_x[iQ2][ix] = new TH2D(Form("h2d_hadron_%d_z_vs_eta_gen_in_Q2_x_%d_%d",hadron_id,iQ2,ix),"D0 z vs eta",100,0,1,40,-10,10);
          h2d_hadron_z_vs_eta_gen_in_Q2_x[iQ2][ix]->Sumw2();
        }
      }

      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        for (int inu = 0; inu < nubin; ++inu)
        {
          h2d_hadron_pt_vs_eta_gen_in_Q2_nu[iQ2][inu] = new TH2D(Form("h2d_hadron_%d_pt_vs_eta_gen_in_Q2_nu_%d_%d",hadron_id,iQ2,inu),"D0 pt vs eta",100,0,10,40,-10,10);
          h2d_hadron_pt_vs_eta_gen_in_Q2_nu[iQ2][inu]->Sumw2();

          h2d_hadron_z_vs_eta_gen_in_Q2_nu[iQ2][inu] = new TH2D(Form("h2d_hadron_%d_z_vs_eta_gen_in_Q2_nu_%d_%d",hadron_id,iQ2,inu),"D0 z vs eta",100,0,1,40,-10,10);
          h2d_hadron_z_vs_eta_gen_in_Q2_nu[iQ2][inu]->Sumw2();
        }
      }

      for (int ix = 0; ix < xbin; ++ix)
      {
        for (int inu = 0; inu < nubin; ++inu)
        {
          h2d_hadron_pt_vs_eta_gen_in_x_nu[ix][inu] = new TH2D(Form("h2d_hadron_%d_pt_vs_eta_gen_in_x_nu_%d_%d",hadron_id,ix,inu),"D0 pt vs eta",100,0,10,40,-10,10);
          h2d_hadron_pt_vs_eta_gen_in_x_nu[ix][inu]->Sumw2();

          h2d_hadron_z_vs_eta_gen_in_x_nu[ix][inu] = new TH2D(Form("h2d_hadron_%d_z_vs_eta_gen_in_x_nu_%d_%d",hadron_id,ix,inu),"D0 z vs eta",100,0,1,40,-10,10);
          h2d_hadron_z_vs_eta_gen_in_x_nu[ix][inu]->Sumw2();
        }
      }
    }
    virtual ~hadron_gen()
    {
      delete h1d_parent_id;
      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        for (int ix = 0; ix < xbin; ++ix)
        {
          delete h2d_hadron_pt_vs_eta_gen_in_Q2_x[iQ2][ix];
          delete h2d_hadron_z_vs_eta_gen_in_Q2_x[iQ2][ix];
        }
      }
      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        for (int inu = 0; inu < nubin; ++inu)
        {
          delete h2d_hadron_pt_vs_eta_gen_in_Q2_nu[iQ2][inu];
          delete h2d_hadron_z_vs_eta_gen_in_Q2_nu[iQ2][inu];
        }
      }
      for (int ix = 0; ix < xbin; ++ix)
      {
        for (int inu = 0; inu < nubin; ++inu)
        {
          delete h2d_hadron_pt_vs_eta_gen_in_x_nu[ix][inu];
          delete h2d_hadron_z_vs_eta_gen_in_x_nu[ix][inu];
        }
      }
    };

    void Reset()
    {
      x_true = 1E1;
      Q2_true = 1E-5;
      nu_true = -9999;

      iQ2bin = -9999;
      ixbin = -9999;

      h1d_parent_id->Reset("ICESM");
      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        for (int ix = 0; ix < xbin; ++ix)
        {
          h2d_hadron_pt_vs_eta_gen_in_Q2_x[iQ2][ix]->Reset("ICESM");
          h2d_hadron_z_vs_eta_gen_in_Q2_x[iQ2][ix]->Reset("ICESM");
        }
      }
      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        for (int inu = 0; inu < nubin; ++inu)
        {
          h2d_hadron_pt_vs_eta_gen_in_Q2_nu[iQ2][inu]->Reset("ICESM");
          h2d_hadron_z_vs_eta_gen_in_Q2_nu[iQ2][inu]->Reset("ICESM");
        }
      }
      for (int ix = 0; ix < xbin; ++ix)
      {
        for (int inu = 0; inu < nubin; ++inu)
        {
          h2d_hadron_pt_vs_eta_gen_in_x_nu[ix][inu]->Reset("ICESM");
          h2d_hadron_z_vs_eta_gen_in_x_nu[ix][inu]->Reset("ICESM");
        }
      }
    }

    void SetQ2True(double _Q2_true)
    { 
      Q2_true = _Q2_true; 

      iQ2bin = -9999;
      for (int iQ2 = 0; iQ2 < Q2bin-1; ++iQ2)
      { // NB: do not loop the last bin which is inclusive
        if (Q2_true>=Q2_lo[iQ2] && Q2_true<Q2_hi[iQ2]) iQ2bin = iQ2;
      }      
    }

    void SetXTrue(double _x_true)
    { 
      x_true = _x_true; 

      ixbin = -9999;
      for (int ix = 0; ix < xbin-1; ++ix)
      { // NB: do not loop the last bin which is inclusive
        if (x_true>=x_lo[ix] && x_true<x_hi[ix]) ixbin = ix;
      }
    }

    void SetNuTrue(double _nu_true)
    {
      nu_true = _nu_true;

      inubin = -9999;
      for (int inu = 0; inu < nubin-1; ++inu)
      { // NB: do not loop the last bin which is inclusive
        if (nu_true>=nu_lo[inu] && nu_true<nu_hi[inu]) inubin = inu;
      }
    }

    void FillGenKin(erhic::EventMC* py_evt)
    {
      erhic::ParticleMC* proton = py_evt->GetTrack(1);
      if (proton!=NULL)
      {
        assert(abs(proton->Id())!=2212);
        hadron_beam = proton->Get4Vector();
      }
      else return; // if incoming proton not found, skip the whole event

      TVector3 boost_vec = hadron_beam.BoostVector();

      for(int ipart = 0; ipart < py_evt->GetNTracks(); ipart++)
      {
        erhic::ParticleMC* part = py_evt->GetTrack(ipart);

        if (abs(part->Id())!=hadron_id) continue;
        if (part->GetStatus()!=1 && part->GetStatus()!=11 && part->GetStatus()!=2) continue; // 1 -- stable particles, 11 -- decay particles in Pythia, 2 -- decay particles in BeAGLE

        // if (hadron_id==421 && abs(part->GetParentId())==4) continue; // FIX ME: temp change, select D0 from charm

        h1d_parent_id->Fill(abs(part->GetParentId()));

        TLorentzVector hadron_mom4_gen = part->Get4Vector();
        double frag_z = hadron_beam.Dot(hadron_mom4_gen)/(nu_true*hadron_beam.M());

        // if (hadron_mom4_gen.Pt()<0.1) continue;

        if (verbosity>1) std::cout << "Particle id " << hadron_id << " with pt " << hadron_mom4_gen.Pt() << " z " << frag_z << " eta " << hadron_mom4_gen.PseudoRapidity() << std::endl;
        if (fabs(frag_z-part->GetZ())>0.01)
        {
          std::cout << "HUGE ISSUE!!!" << frag_z << " " << part->GetZ() << endl;
          TLorentzVector hadron_mom4_gen_rest = hadron_mom4_gen;
          hadron_mom4_gen_rest.Boost(-boost_vec); 
          cout << "DIAG " << part->GetPt() << " => " << part->GetZ() << " with nu " << nu_true << " hadron beam E " << hadron_beam.E() << " produced hadron " << hadron_mom4_gen_rest.E() << " Eh/nu = " << hadron_mom4_gen_rest.E()/nu_true << endl;
        } 
        // if (nu_true>=50 && nu_true<100 && x_true>0.1 && x_true<1 && part->GetEta()>=1 && part->GetEta()<3 && part->GetZ()>0.82 && part->GetZ()<0.92)
        // {
          
        //   TLorentzVector hadron_mom4_gen_rest = hadron_mom4_gen;
        //   hadron_mom4_gen_rest.Boost(-boost_vec); 
        //   cout << "DIAG " << part->GetPt() << " => " << part->GetZ() << " with nu " << nu_true << " hadron beam E " << hadron_beam.E() << " produced hadron " << hadron_mom4_gen_rest.E() << " Eh/nu = " << hadron_mom4_gen_rest.E()/nu_true << endl;
        // }

        if (verbosity>1) std::cout << "Q2, x " << Q2_true << ", " << x_true << " bin " << iQ2bin << ", " << ixbin << std::endl;

        if (iQ2bin>=0 && ixbin>=0)
        {
          h2d_hadron_pt_vs_eta_gen_in_Q2_x[iQ2bin][ixbin]->Fill(hadron_mom4_gen.Pt(),hadron_mom4_gen.PseudoRapidity());
          h2d_hadron_z_vs_eta_gen_in_Q2_x[iQ2bin][ixbin]->Fill(frag_z,hadron_mom4_gen.PseudoRapidity());

          // Now fill the inclusive bin
          h2d_hadron_pt_vs_eta_gen_in_Q2_x[iQ2bin][xbin-1]->Fill(hadron_mom4_gen.Pt(),hadron_mom4_gen.PseudoRapidity());
          h2d_hadron_z_vs_eta_gen_in_Q2_x[iQ2bin][xbin-1]->Fill(frag_z,hadron_mom4_gen.PseudoRapidity());
          h2d_hadron_pt_vs_eta_gen_in_Q2_x[Q2bin-1][ixbin]->Fill(hadron_mom4_gen.Pt(),hadron_mom4_gen.PseudoRapidity());
          h2d_hadron_z_vs_eta_gen_in_Q2_x[Q2bin-1][ixbin]->Fill(frag_z,hadron_mom4_gen.PseudoRapidity());
          h2d_hadron_pt_vs_eta_gen_in_Q2_x[Q2bin-1][xbin-1]->Fill(hadron_mom4_gen.Pt(),hadron_mom4_gen.PseudoRapidity());
          h2d_hadron_z_vs_eta_gen_in_Q2_x[Q2bin-1][xbin-1]->Fill(frag_z,hadron_mom4_gen.PseudoRapidity());
        } 

        if (iQ2bin>=0 && inubin>=0)
        {
          h2d_hadron_pt_vs_eta_gen_in_Q2_nu[iQ2bin][inubin]->Fill(hadron_mom4_gen.Pt(),hadron_mom4_gen.PseudoRapidity());
          h2d_hadron_z_vs_eta_gen_in_Q2_nu[iQ2bin][inubin]->Fill(frag_z,hadron_mom4_gen.PseudoRapidity());

          // Now fill the inclusive bin
          h2d_hadron_pt_vs_eta_gen_in_Q2_nu[iQ2bin][nubin-1]->Fill(hadron_mom4_gen.Pt(),hadron_mom4_gen.PseudoRapidity());
          h2d_hadron_z_vs_eta_gen_in_Q2_nu[iQ2bin][nubin-1]->Fill(frag_z,hadron_mom4_gen.PseudoRapidity());
          h2d_hadron_pt_vs_eta_gen_in_Q2_nu[Q2bin-1][inubin]->Fill(hadron_mom4_gen.Pt(),hadron_mom4_gen.PseudoRapidity());
          h2d_hadron_z_vs_eta_gen_in_Q2_nu[Q2bin-1][inubin]->Fill(frag_z,hadron_mom4_gen.PseudoRapidity());
          h2d_hadron_pt_vs_eta_gen_in_Q2_nu[Q2bin-1][nubin-1]->Fill(hadron_mom4_gen.Pt(),hadron_mom4_gen.PseudoRapidity());
          h2d_hadron_z_vs_eta_gen_in_Q2_nu[Q2bin-1][nubin-1]->Fill(frag_z,hadron_mom4_gen.PseudoRapidity());
        }

        if (ixbin>=0 && inubin>=0)
        {
          h2d_hadron_pt_vs_eta_gen_in_x_nu[ixbin][inubin]->Fill(hadron_mom4_gen.Pt(),hadron_mom4_gen.PseudoRapidity());
          h2d_hadron_z_vs_eta_gen_in_x_nu[ixbin][inubin]->Fill(frag_z,hadron_mom4_gen.PseudoRapidity());

          // Now fill the inclusive bin
          h2d_hadron_pt_vs_eta_gen_in_x_nu[ixbin][nubin-1]->Fill(hadron_mom4_gen.Pt(),hadron_mom4_gen.PseudoRapidity());
          h2d_hadron_z_vs_eta_gen_in_x_nu[ixbin][nubin-1]->Fill(frag_z,hadron_mom4_gen.PseudoRapidity());
          h2d_hadron_pt_vs_eta_gen_in_x_nu[xbin-1][inubin]->Fill(hadron_mom4_gen.Pt(),hadron_mom4_gen.PseudoRapidity());
          h2d_hadron_z_vs_eta_gen_in_x_nu[xbin-1][inubin]->Fill(frag_z,hadron_mom4_gen.PseudoRapidity());
          h2d_hadron_pt_vs_eta_gen_in_x_nu[xbin-1][nubin-1]->Fill(hadron_mom4_gen.Pt(),hadron_mom4_gen.PseudoRapidity());
          h2d_hadron_z_vs_eta_gen_in_x_nu[xbin-1][nubin-1]->Fill(frag_z,hadron_mom4_gen.PseudoRapidity());
        } 
      }
    }

    void Write()
    {
      h1d_parent_id->Write();
      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        for (int ix = 0; ix < xbin; ++ix)
        {
          h2d_hadron_pt_vs_eta_gen_in_Q2_x[iQ2][ix]->Write();
          h2d_hadron_z_vs_eta_gen_in_Q2_x[iQ2][ix]->Write();
        }
      }
      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        for (int inu = 0; inu < nubin; ++inu)
        {
          h2d_hadron_pt_vs_eta_gen_in_Q2_nu[iQ2][inu]->Write();
          h2d_hadron_z_vs_eta_gen_in_Q2_nu[iQ2][inu]->Write();
        }
      }
      for (int ix = 0; ix < xbin; ++ix)
      {
        for (int inu = 0; inu < nubin; ++inu)
        {
          h2d_hadron_pt_vs_eta_gen_in_x_nu[ix][inu]->Write();
          h2d_hadron_z_vs_eta_gen_in_x_nu[ix][inu]->Write();
        }
      }
    }
};

bool event_w_charm(erhic::EventMC* event, int gen_type)
{
  if (gen_type==0)
  { // Pythia 6
    bool flag_search_group = false; // flag if the group of particles has been searched
    for (int ipart = 0; ipart < event->GetNTracks(); ++ipart)
    {
      erhic::ParticleMC* part = event->GetTrack(ipart);

      if ( part->GetStatus()!=21 && !flag_search_group ) continue;
      if ( part->GetStatus()==21 )
      { // entering into the group of particles with KS=21
        flag_search_group = true;
        if ( part->Id()==4 || part->Id()==-4 ) return true;
      }
      if ( part->GetStatus()!=21 && flag_search_group ) break;
    }
  }
  else
  { // BeAGLE
    bool flag_search_group = false; // flag if the group of particles has been searched
    for (int ipart = event->GetNTracks()-1; ipart >=0; ipart--)
    { // faster looping backwards
      erhic::ParticleMC* part = event->GetTrack(ipart);

      if ( part->GetStatus()!=3 && !flag_search_group ) continue;
      if ( part->GetStatus()==3 )
      { // entering into the group of particles with KS=3
        flag_search_group = true;
        if ( part->Id()==4 || part->Id()==-4 ) return true;
      }
      if ( part->GetStatus()!=3 && flag_search_group ) break;
    }
  }

  return false;
}

bool event_selection(erhic::EventMC* event, int data_type)
{
  if (data_type<2)
  { // EIC event selection (0.01<y<0.85, modifed from Yuxiang's version 0.05<y<0.8)
    if (event->GetY()<0.01 || event->GetY()>0.85) return false;
  }
  else if (data_type==2)
  { // HERMES event selections (from HERMES paper)
    if (event->GetY()>0.85) return false;
    if (event->GetW2()<4) return false;
    if (event->GetNu()<6) return false;
  }
  else
  { // CLAS event selection (to be implemented)

  }

  return true;
}

void ana_hadron_gen(const char* inFile = "eA.root", const char* outFile = "hist.root", int nevt = 0, int data_type = 0, int gen_type = 1)
{
  cout << "Data Type: "; 
  if (data_type==0) cout << "EIC" << endl;
  else if (data_type==1) cout << "HERMES" << endl;
  else cout << "CLAS" << endl;

  cout << "Generator Type: "; 
  if (gen_type==0) cout << "Pythia6" << endl;
  else cout << "BeAGLE" << endl;

  TFile *f = new TFile(inFile);

  //Get EICTree Tree
  TTree *tree = (TTree*)f->Get("EICTree");

  Int_t nEntries = tree->GetEntries();
  cout<<"-------------------------------"<<endl;
  cout<<"Total Number of Events = "<<nEntries<<endl<<endl;

  // Event Class
  erhic::EventMC *event_pythia(NULL);
  erhic::EventBeagle *event_beagle(NULL);

  // Access event Branch
  if (gen_type==0) tree->SetBranchAddress("event",&event_pythia);
  else tree->SetBranchAddress("event",&event_beagle);

  // Inclusive event counter before any event selections
  TH2D* h2d_logx_logQ2 = new TH2D("h2d_logx_logQ2","inclusive event counter",100,-4,0,100,0,3);
  h2d_logx_logQ2->Sumw2();
  TH2D* h2d_logx_lognu = new TH2D("h2d_logx_lognu","inclusive event counter",100,-4,0,100,0,3);
  h2d_logx_lognu->Sumw2();

  // Event counter after event selections
  TH1D* h1d_nevt_in_Q2_x[Q2bin][xbin] = {0};
  TH1D* h1d_nevt_w_charm_in_Q2_x[Q2bin][xbin] = {0};
  for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
  {
    for (int ix = 0; ix < xbin; ++ix)
    {
      h1d_nevt_in_Q2_x[iQ2][ix] = new TH1D(Form("h1d_nevt_in_Q2_x_%d_%d",iQ2,ix),"event counter",1,0.5,1.5);
      h1d_nevt_in_Q2_x[iQ2][ix]->Sumw2(); 

      h1d_nevt_w_charm_in_Q2_x[iQ2][ix] = new TH1D(Form("h1d_nevt_w_charm_in_Q2_x_%d_%d",iQ2,ix),"event counter",1,0.5,1.5);
      h1d_nevt_w_charm_in_Q2_x[iQ2][ix]->Sumw2(); 
    }
  }

  TH1D* h1d_nevt_in_Q2_nu[Q2bin][nubin] = {0};
  TH1D* h1d_nevt_w_charm_in_Q2_nu[Q2bin][nubin] = {0};
  for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
  {
    for (int inu = 0; inu < nubin; ++inu)
    {
      h1d_nevt_in_Q2_nu[iQ2][inu] = new TH1D(Form("h1d_nevt_in_Q2_nu_%d_%d",iQ2,inu),"event counter",1,0.5,1.5);
      h1d_nevt_in_Q2_nu[iQ2][inu]->Sumw2(); 

      h1d_nevt_w_charm_in_Q2_nu[iQ2][inu] = new TH1D(Form("h1d_nevt_w_charm_in_Q2_nu_%d_%d",iQ2,inu),"event counter",1,0.5,1.5);
      h1d_nevt_w_charm_in_Q2_nu[iQ2][inu]->Sumw2(); 
    }
  }

  TH1D* h1d_nevt_in_x_nu[xbin][nubin] = {0};
  TH1D* h1d_nevt_w_charm_in_x_nu[xbin][nubin] = {0};
  for (int ix = 0; ix < xbin; ++ix)
  {
    for (int inu = 0; inu < nubin; ++inu)
    {
      h1d_nevt_in_x_nu[ix][inu] = new TH1D(Form("h1d_nevt_in_x_nu_%d_%d",ix,inu),"event counter",1,0.5,1.5);
      h1d_nevt_in_x_nu[ix][inu]->Sumw2(); 

      h1d_nevt_w_charm_in_x_nu[ix][inu] = new TH1D(Form("h1d_nevt_w_charm_in_x_nu_%d_%d",ix,inu),"event counter",1,0.5,1.5);
      h1d_nevt_w_charm_in_x_nu[ix][inu]->Sumw2(); 
    }
  }

  // Hadron analyzers
  hadron_gen ana_D0(421); // which will also analyze -421
  hadron_gen ana_Lc(4122);
  hadron_gen ana_D(411);
  hadron_gen ana_Jpsi(443);
  hadron_gen ana_B0(511);
  hadron_gen ana_B(521);
  hadron_gen ana_charged_pion(211);
  hadron_gen ana_neutral_pion(111);
  hadron_gen ana_charged_kaon(321);
  hadron_gen ana_proton(2212);
  
  //Loop Over Events
  if (nevt == 0) nevt = nEntries;
  for(Int_t ievt = 0; ievt < nevt; ievt++)
  {    
    if (ievt%10000==0) cout<<"Processing event = "<<ievt<<"/"<<nevt<<endl;

    tree->GetEntry(ievt);

    erhic::EventMC* event = NULL;

    if (gen_type==0) event = event_pythia;
    else event = event_beagle;

    //Write Out Q2
    double Q2 = event->GetQ2(); // Can also do event->QSquared
    double x = event->GetX();
    double nu = event->GetNu();

    h2d_logx_logQ2->Fill(log10(x),log10(Q2));
    h2d_logx_lognu->Fill(log10(x),log10(nu));
    // printf("For Event %d, Q^2 = %.3f GeV^2!\n",ievt,Q2);

    bool flag_event_select = true;
    event_selection(event, data_type);
    
    if (!flag_event_select) continue;

    int iQ2bin = -9999, ixbin = -9999, inubin = -9999;
    for (int iQ2 = 0; iQ2 < Q2bin-1; ++iQ2)
    { // NB: do not loop the last bin which is inclusive
      if (Q2>=Q2_lo[iQ2] && Q2<Q2_hi[iQ2]) iQ2bin = iQ2;
    } 
    for (int ix = 0; ix < xbin-1; ++ix)
    { // NB: do not loop the last bin which is inclusive
      if (x>=x_lo[ix] && x<x_hi[ix]) ixbin = ix;
    } 
    for (int inu = 0; inu < nubin-1; ++inu)
    { // NB: do not loop the last bin which is inclusive
      if (nu>=nu_lo[inu] && nu<nu_hi[inu]) inubin = inu;
    }

    cout << endl << "[i] Event " << ievt << endl;

    for(int ipart = 0; ipart < event->GetNTracks(); ipart++)
      {
	erhic::ParticleMC* part = event->GetTrack(ipart);
	
	//if (abs(part->Id())!=hadron_id) continue;
	//if (part->GetStatus()!=1 && part->GetStatus()!=11 && part->GetStatus()!=2) continue; // 1 -- stable particles, 11 -- decay particles in Pythia, 2 -- decay particles in BeAGLE
	
	//if (fabs(part->Id())==421 && abs(part->GetParentId())==4) continue;

	if (abs(part->Id())==4)
	  {
	    cout << "[i] Found a charm quark" << endl;
	  }
	if (abs(part->Id())==421)
	  {
	    cout << "[i] Found a D0" << endl;
	  }
      }
 
    continue;

    if (iQ2bin>=0 && ixbin>=0)
    {
      h1d_nevt_in_Q2_x[iQ2bin][ixbin]->Fill(1);
      h1d_nevt_in_Q2_x[iQ2bin][xbin-1]->Fill(1);
      h1d_nevt_in_Q2_x[Q2bin-1][ixbin]->Fill(1);
      h1d_nevt_in_Q2_x[Q2bin-1][xbin-1]->Fill(1);

      if ( event_w_charm(event,gen_type) )
      {
        h1d_nevt_w_charm_in_Q2_x[iQ2bin][ixbin]->Fill(1);
        h1d_nevt_w_charm_in_Q2_x[iQ2bin][xbin-1]->Fill(1);
        h1d_nevt_w_charm_in_Q2_x[Q2bin-1][ixbin]->Fill(1);
        h1d_nevt_w_charm_in_Q2_x[Q2bin-1][xbin-1]->Fill(1);
      }
    } 

    if (iQ2bin>=0 && inubin>=0)
    {
      h1d_nevt_in_Q2_nu[iQ2bin][inubin]->Fill(1);
      h1d_nevt_in_Q2_nu[iQ2bin][nubin-1]->Fill(1);
      h1d_nevt_in_Q2_nu[Q2bin-1][inubin]->Fill(1);
      h1d_nevt_in_Q2_nu[Q2bin-1][nubin-1]->Fill(1);

      if ( event_w_charm(event,gen_type) )
      {
        h1d_nevt_w_charm_in_Q2_nu[iQ2bin][inubin]->Fill(1);
        h1d_nevt_w_charm_in_Q2_nu[iQ2bin][nubin-1]->Fill(1);
        h1d_nevt_w_charm_in_Q2_nu[Q2bin-1][inubin]->Fill(1);
        h1d_nevt_w_charm_in_Q2_nu[Q2bin-1][nubin-1]->Fill(1);
      }
    }

    if (ixbin>=0 && inubin>=0)
    {
      h1d_nevt_in_x_nu[ixbin][inubin]->Fill(1);
      h1d_nevt_in_x_nu[ixbin][nubin-1]->Fill(1);
      h1d_nevt_in_x_nu[xbin-1][inubin]->Fill(1);
      h1d_nevt_in_x_nu[xbin-1][nubin-1]->Fill(1);

      if ( event_w_charm(event,gen_type) )
      {
        h1d_nevt_w_charm_in_x_nu[ixbin][inubin]->Fill(1);
        h1d_nevt_w_charm_in_x_nu[ixbin][nubin-1]->Fill(1);
        h1d_nevt_w_charm_in_x_nu[xbin-1][inubin]->Fill(1);
        h1d_nevt_w_charm_in_x_nu[xbin-1][nubin-1]->Fill(1);
      }
    } 

    // for(int ipart = 0; ipart < event->GetNTracks(); ipart++)
    // {
    //   erhic::ParticleMC* part = event->GetTrack(ipart);

    //   if (abs(part->Id())!=321) continue;
    //   if (part->GetStatus()!=1 && part->GetStatus()!=11 && part->GetStatus()!=2) continue; // 1 -- stable particles, 11 -- decay particles in Pythia, 2 -- decay particles in BeAGLE

    //   erhic::ParticleMCeA* part_eA = part->eA;
    //   cout << "part index " << ipart << " pdg id " << part->Id() << " NoBam " << part_eA->NoBam << endl;
    // }
    
    ana_D0.SetQ2True(event->GetQ2());
    ana_D0.SetXTrue(event->GetX());
    ana_D0.SetNuTrue(event->GetNu());
    ana_D0.FillGenKin(event);

    ana_Lc.SetQ2True(event->GetQ2());
    ana_Lc.SetXTrue(event->GetX());
    ana_Lc.SetNuTrue(event->GetNu());
    ana_Lc.FillGenKin(event);

    ana_D.SetQ2True(event->GetQ2());
    ana_D.SetXTrue(event->GetX());
    ana_D.SetNuTrue(event->GetNu());
    ana_D.FillGenKin(event);

    ana_Jpsi.SetQ2True(event->GetQ2());
    ana_Jpsi.SetXTrue(event->GetX());
    ana_Jpsi.SetNuTrue(event->GetNu());
    ana_Jpsi.FillGenKin(event);

    ana_B0.SetQ2True(event->GetQ2());
    ana_B0.SetXTrue(event->GetX());
    ana_B0.SetNuTrue(event->GetNu());
    ana_B0.FillGenKin(event);

    ana_B.SetQ2True(event->GetQ2());
    ana_B.SetXTrue(event->GetX());
    ana_B.SetNuTrue(event->GetNu());
    ana_B.FillGenKin(event);

    ana_charged_pion.SetQ2True(event->GetQ2());
    ana_charged_pion.SetXTrue(event->GetX());
    ana_charged_pion.SetNuTrue(event->GetNu());
    ana_charged_pion.FillGenKin(event);

    ana_neutral_pion.SetQ2True(event->GetQ2());
    ana_neutral_pion.SetXTrue(event->GetX());
    ana_neutral_pion.SetNuTrue(event->GetNu());
    ana_neutral_pion.FillGenKin(event);

    ana_charged_kaon.SetQ2True(event->GetQ2());
    ana_charged_kaon.SetXTrue(event->GetX());
    ana_charged_kaon.SetNuTrue(event->GetNu());
    ana_charged_kaon.FillGenKin(event);

    ana_proton.SetQ2True(event->GetQ2());
    ana_proton.SetXTrue(event->GetX());
    ana_proton.SetNuTrue(event->GetNu());
    ana_proton.FillGenKin(event);
  }

  TFile* fout = new TFile(outFile,"recreate");
  h2d_logx_logQ2->Write();
  h2d_logx_lognu->Write();
  for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
  {
    for (int ix = 0; ix < xbin; ++ix)
    {
      h1d_nevt_in_Q2_x[iQ2][ix]->Write();
      h1d_nevt_w_charm_in_Q2_x[iQ2][ix]->Write();
    }
  }
  for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
  {
    for (int inu = 0; inu < nubin; ++inu)
    {
      h1d_nevt_in_Q2_nu[iQ2][inu]->Write();
      h1d_nevt_w_charm_in_Q2_nu[iQ2][inu]->Write();
    }
  }
  for (int ix = 0; ix < xbin; ++ix)
  {
    for (int inu = 0; inu < nubin; ++inu)
    {
      h1d_nevt_in_x_nu[ix][inu]->Write();
      h1d_nevt_w_charm_in_x_nu[ix][inu]->Write();
    }
  }
  
  ana_D0.Write();
  ana_Lc.Write();
  ana_D.Write();
  ana_Jpsi.Write();
  ana_B0.Write();
  ana_B.Write();
  ana_charged_pion.Write();
  ana_neutral_pion.Write();
  ana_charged_kaon.Write();
  ana_proton.Write();

  fout->Write();
}
