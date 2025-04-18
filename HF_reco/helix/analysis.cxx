#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class PlotFile;
#endif

#ifndef __CINT__
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include "math.h"
#include "string.h"

#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TSystem.h"
#include "TUnixSystem.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TLatex.h"
#endif


#include "StPhysicalHelix.h"
#include "SystemOfUnits.h"
#include "PhysicalConstants.h"

using namespace std;

const double gPionMass = 0.13957;
const double gKaonMass = 0.493677;

const double twoPi = 2.*3.1415927;
const double eMass = 0.000511;

const double bField = -1.7; // Tesla

TVector3 getDcaToVtx(const int index, TVector3 vtx);

TLorentzVector getPairParent(const int index1, const int index2, TVector3 vtx,
			     float &dcaDaughters, float &cosTheta, float &decayLength, float &V0DcaToVtx);

TTreeReaderArray<float> *rcMomPx2;
TTreeReaderArray<float> *rcMomPy2;
TTreeReaderArray<float> *rcMomPz2;
TTreeReaderArray<float> *rcCharge2;

TTreeReaderArray<float> *rcTrkLoca2;
TTreeReaderArray<float> *rcTrkLocb2;
TTreeReaderArray<float> *rcTrkTheta2;
TTreeReaderArray<float> *rcTrkPhi2;

int main(int argc, char **argv)
{
  if(argc!=3 && argc!=1) return 0;

  TString listname;
  TString outname;

  if(argc==1)
    {
      listname  = "test.list";
      outname = "test.root";
    }

  if(argc==3)
    {
      listname = argv[1];
      outname = argv[2];
    }  

  TChain *chain = new TChain("events");

  int nfiles = 0;
  char filename[512];
  ifstream *inputstream = new ifstream;
  inputstream->open(listname.Data());
  if(!inputstream)
    {
      printf("[e] Cannot open file list: %s\n", listname.Data());
    }
  while(inputstream->good())
    {
      inputstream->getline(filename, 512);
      if(inputstream->good())
	{
	  TFile *ftmp = TFile::Open(filename, "read");
	  if(!ftmp||!(ftmp->IsOpen())||!(ftmp->GetNkeys())) 
	    {
	      printf("[e] Could you open file: %s\n", filename);
	    } 
	  else
	    {
	      cout<<"[i] Add "<<nfiles<<"th file: "<<filename<<endl;
	      chain->Add(filename);
	      nfiles++;
	    }
	}
    }
  inputstream->close();
  printf("[i] Read in %d files with %lld events in total\n", nfiles, chain->GetEntries());

  TH1F *hEventStat = new TH1F("hEventStat", "Event statistics", 5, 0, 5);
  hEventStat->GetXaxis()->SetBinLabel(1, "MC events");
  hEventStat->GetXaxis()->SetBinLabel(2, "D0");
  hEventStat->GetXaxis()->SetBinLabel(3, "D0 -> pi+K");
  hEventStat->GetXaxis()->SetBinLabel(4, "Reco D0");

  TH1F *hMcMult = new TH1F("hMcMult", "MC multiplicity (|#eta| < 3.5);N_{MC}", 50, 0, 50);

  TH1F *hMcVtxX = new TH1F("hMcVtxX", "x position of MC vertex;x (mm)", 100, -5.05, 4.95);
  TH1F *hMcVtxY = new TH1F("hMcVtxY", "y position of MC vertex;y (mm)", 500, -5.01, 4.99);
  TH1F *hMcVtxZ = new TH1F("hMcVtxZ", "z position of MC vertex;z (mm)", 400, -200, 200);

  TH2F *hD0DecayVxVy = new TH2F("hD0DecayVxVy", "D^{0} decay vertex to primary vertex;#Deltav_{x} (mm);#Deltav_{y} (mm)", 400, -1-0.0025, 1-0.0025, 400, -1-0.0025, 1-0.0025);
  TH2F *hD0DecayVrVz = new TH2F("hD0DecayVrVz", "D^{0} decay vertex to primary vertex;#Deltav_{z} (mm);#Deltav_{r} (mm)", 100, -2, 2, 100, -0.2, 1.8);

  TH2F *hMCD0PtRap = new TH2F("hMCD0PtRap", "MC D^{0};y;p_{T} (GeV/c)", 20, -5, 5, 100, 0, 10);

  TH2F *hMcPiPtEta = new TH2F("hMcPiPtEta", "MC #pi from D^{0} decay;#eta^{MC};p_{T}^{MC} (GeV/c)", 20, -5, 5, 100, 0, 10);
  TH2F *hMcPiPtEtaReco = new TH2F("hMcPiPtEtaReco", "RC #pi from D^{0} decay;#eta^{MC};p_{T}^{MC} (GeV/c)", 20, -5, 5, 100, 0, 10);

  TH2F *hMcKPtEta = new TH2F("hMcKPtEta", "MC K from D^{0} decay;#eta^{MC};p_{T}^{MC} (GeV/c)", 20, -5, 5, 100, 0, 10);
  TH2F *hMcKPtEtaReco = new TH2F("hMcKPtEtaReco", "RC K from D^{0} decay;#eta^{MC};p_{T}^{MC} (GeV/c)", 20, -5, 5, 100, 0, 10);

  TH1F *hNRecoVtx = new TH1F("hNRecoVtx", "Number of reconstructed vertices;N", 10, 0, 10);

  const char* part_name[3] = {"Pi", "K", "P"};
  const char* part_title[3] = {"#pi", "K", "P"};
  TH3F *hRcSecPartLocaToRCVtx[2];
  TH3F *hRcSecPartLocbToRCVtx[2];
  TH3F *hRcPrimPartLocaToRCVtx[2];
  TH3F *hRcPrimPartLocbToRCVtx[2];
  for(int i=0; i<2; i++)
    {
      hRcSecPartLocaToRCVtx[i] = new TH3F(Form("hRcSec%sLocaToRCVtx",part_name[i]), Form( "DCA_{xy} distribution for D^{0} decayed %s;p_{T} (GeV/c);#eta;DCA_{xy} (mm)", part_title[i]), 100, 0, 10, 20, -5, 5, 100, 0, 1);
      hRcSecPartLocbToRCVtx[i] = new TH3F(Form("hRcSec%sLocbToRCVtx",part_name[i]), Form( "DCA_{z} distribution for D^{0} decayed %s;p_{T} (GeV/c);#eta;DCA_{z} (mm)", part_title[i]), 100, 0, 10, 20, -5, 5, 100, -0.5, 0.5);
      hRcPrimPartLocaToRCVtx[i] = new TH3F(Form("hRcPrim%sLocaToRCVtx",part_name[i]), Form( "DCA_{xy} distribution for primary %s;p_{T} (GeV/c);#eta;DCA_{xy} (mm)", part_title[i]), 100, 0, 10, 20, -5, 5, 100, 0, 1);
      hRcPrimPartLocbToRCVtx[i] = new TH3F(Form("hRcPrim%sLocbToRCVtx",part_name[i]), Form( "DCA_{z} distribution for primary %s;p_{T} (GeV/c);#eta;DCA_{z} (mm)", part_title[i]), 100, 0, 10, 20, -5, 5, 100, -0.5, 0.5);
    }

  const char* axis_name[3] = {"x", "y", "z"};
  const int nDimDca = 4;
  const int nBinsDca[nDimDca] = {50, 20, 500, 50};
  const double minBinDca[nDimDca] = {0, -5, -1+0.002, 0};
  const double maxBinDca[nDimDca] = {5, 5, 1+0.002, 50};
  THnSparseF *hPrimTrkDcaToRCVtx[3][3];
  for(int i=0; i<3; i++)
    {
      for(int j=0; j<3; j++)
	{
	  hPrimTrkDcaToRCVtx[i][j] = new THnSparseF(Form("hPrim%sDca%sToRCVtx",part_name[i],axis_name[j]), Form("DCA_{%s} distribution for primary %s;p_{T} (GeV/c);#eta;DCA_{%s} (mm);N_{MC}",axis_name[j],part_title[i],axis_name[j]), nDimDca, nBinsDca, minBinDca, maxBinDca);
	}
    }

  TH3F *h3PairDca12[2];
  TH3F *h3PairCosTheta[2];
  TH3F *h3PairDca[2];
  TH3F *h3PairDecayLength[2];
  const char* pair_name[2] = {"signal", "bkg"};
  const char* pair_title[2] = {"Signal", "Background"};
  for(int i=0; i<2; i++)
    {
      h3PairDca12[i] = new TH3F(Form("h3PairDca12_%s", pair_name[i]), Form("%s pair DCA_{12};p_{T} (GeV/c);#eta;DCA_{12} (mm)", pair_title[i]), 100, 0, 10, 20, -5, 5, 100, 0, 1);

      h3PairCosTheta[i] = new TH3F(Form("h3PairCosTheta_%s", pair_name[i]), Form("%s pair cos(#theta);p_{T} (GeV/c);#eta;cos(#theta)", pair_title[i]), 100, 0, 10, 20, -5, 5, 100, -1, 1);

      h3PairDca[i] = new TH3F(Form("h3PairDca_%s", pair_name[i]), Form("%s pair DCA;p_{T} (GeV/c);#eta;DCA_{pair} (mm)", pair_title[i]), 100, 0, 10, 20, -5, 5, 100, 0, 1);

      h3PairDecayLength[i] = new TH3F(Form("h3PairDecayLength_%s", pair_name[i]), Form("%s pair decay length;p_{T} (GeV/c);#eta;L (mm)", pair_title[i]), 100, 0, 10, 20, -5, 5, 100, 0, 1);
    }

  // Invariant mass
  const char* cut_name[2] = {"all", "DCA"};
  TH3F *h3InvMass[2][2];
  for(int i=0; i<2; i++)
    {
      for(int j=0; j<2; j++)
	{
	  h3InvMass[i][j] = new TH3F(Form("h3InvMass_%s_%s", pair_name[i], cut_name[j]), "Invariant mass of unlike-sign #piK pairs;p_{T} (GeV/c);y;M_{#piK} (GeV/c^{2})", 100, 0, 10, 20, -5, 5, 100, 1.6, 2.0);
	}
    }

  TTreeReader treereader(chain);
  // MC  
  TTreeReaderArray<int> mcPartGenStatus = {treereader, "MCParticles.generatorStatus"};
  TTreeReaderArray<int> mcPartPdg = {treereader, "MCParticles.PDG"};
  TTreeReaderArray<float> mcPartCharge = {treereader, "MCParticles.charge"};
  TTreeReaderArray<unsigned int> mcPartParent_begin = {treereader, "MCParticles.parents_begin"};
  TTreeReaderArray<unsigned int> mcPartParent_end = {treereader, "MCParticles.parents_end"};
  TTreeReaderArray<int> mcPartParent_index = {treereader, "_MCParticles_parents.index"};
  TTreeReaderArray<unsigned int> mcPartDaughter_begin = {treereader, "MCParticles.daughters_begin"};
  TTreeReaderArray<unsigned int> mcPartDaughter_end = {treereader, "MCParticles.daughters_end"};
  TTreeReaderArray<int> mcPartDaughter_index = {treereader, "_MCParticles_daughters.index"};
  TTreeReaderArray<double> mcPartMass = {treereader, "MCParticles.mass"};
  TTreeReaderArray<double> mcPartVx = {treereader, "MCParticles.vertex.x"};
  TTreeReaderArray<double> mcPartVy = {treereader, "MCParticles.vertex.y"};
  TTreeReaderArray<double> mcPartVz = {treereader, "MCParticles.vertex.z"};
  TTreeReaderArray<double> mcMomPx = {treereader, "MCParticles.momentum.x"};
  TTreeReaderArray<double> mcMomPy = {treereader, "MCParticles.momentum.y"};
  TTreeReaderArray<double> mcMomPz = {treereader, "MCParticles.momentum.z"};
  TTreeReaderArray<double> mcEndPointX = {treereader, "MCParticles.endpoint.x"};
  TTreeReaderArray<double> mcEndPointY = {treereader, "MCParticles.endpoint.y"};
  TTreeReaderArray<double> mcEndPointZ = {treereader, "MCParticles.endpoint.z"};

  TTreeReaderArray<unsigned int> assocChSimID = {treereader, "ReconstructedChargedParticleAssociations.simID"};
  TTreeReaderArray<unsigned int> assocChRecID = {treereader, "ReconstructedChargedParticleAssociations.recID"};
  
  TTreeReaderArray<float> rcMomPx = {treereader, "ReconstructedChargedParticles.momentum.x"};
  TTreeReaderArray<float> rcMomPy = {treereader, "ReconstructedChargedParticles.momentum.y"};
  TTreeReaderArray<float> rcMomPz = {treereader, "ReconstructedChargedParticles.momentum.z"};
  TTreeReaderArray<float> rcPosx = {treereader, "ReconstructedChargedParticles.referencePoint.x"};
  TTreeReaderArray<float> rcPosy = {treereader, "ReconstructedChargedParticles.referencePoint.y"};
  TTreeReaderArray<float> rcPosz = {treereader, "ReconstructedChargedParticles.referencePoint.z"};
  TTreeReaderArray<float> rcCharge = {treereader, "ReconstructedChargedParticles.charge"};
  TTreeReaderArray<int>   rcPdg = {treereader, "ReconstructedChargedParticles.PDG"};

  TTreeReaderArray<float> rcTrkLoca = {treereader, "CentralCKFTrackParameters.loc.a"};
  TTreeReaderArray<float> rcTrkLocb = {treereader, "CentralCKFTrackParameters.loc.b"};
  TTreeReaderArray<float> rcTrkqOverP = {treereader, "CentralCKFTrackParameters.qOverP"};
  TTreeReaderArray<float> rcTrkTheta = {treereader, "CentralCKFTrackParameters.theta"};
  TTreeReaderArray<float> rcTrkPhi = {treereader, "CentralCKFTrackParameters.phi"};

  rcMomPx2 = new TTreeReaderArray<float>{treereader, "ReconstructedChargedParticles.momentum.x"};
  rcMomPy2 = new TTreeReaderArray<float>{treereader, "ReconstructedChargedParticles.momentum.y"};
  rcMomPz2 = new TTreeReaderArray<float>{treereader, "ReconstructedChargedParticles.momentum.z"};
  rcCharge2 = new TTreeReaderArray<float>{treereader, "ReconstructedChargedParticles.charge"};

  rcTrkLoca2 = new TTreeReaderArray<float>{treereader, "CentralCKFTrackParameters.loc.a"};
  rcTrkLocb2 = new TTreeReaderArray<float>{treereader, "CentralCKFTrackParameters.loc.b"};
  rcTrkTheta2 = new TTreeReaderArray<float>{treereader, "CentralCKFTrackParameters.theta"};
  rcTrkPhi2 = new TTreeReaderArray<float>{treereader, "CentralCKFTrackParameters.phi"};

  TTreeReaderArray<float> CTVx = {treereader, "CentralTrackVertices.position.x"};
  TTreeReaderArray<float> CTVy = {treereader, "CentralTrackVertices.position.y"};
  TTreeReaderArray<float> CTVz = {treereader, "CentralTrackVertices.position.z"};
  TTreeReaderArray<int> CTVndf = {treereader, "CentralTrackVertices.ndf"};
  TTreeReaderArray<float> CTVchi2 = {treereader, "CentralTrackVertices.chi2"};
  TTreeReaderArray<float> CTVerr_xx = {treereader, "CentralTrackVertices.positionError.xx"};
  TTreeReaderArray<float> CTVerr_yy = {treereader, "CentralTrackVertices.positionError.yy"};
  TTreeReaderArray<float> CTVerr_zz = {treereader, "CentralTrackVertices.positionError.zz"};

  TTreeReaderArray<int> prim_vtx_index = {treereader, "PrimaryVertices_objIdx.index"};

  TTreeReaderArray<unsigned int> vtxAssocPart_begin = {treereader, "CentralTrackVertices.associatedParticles_begin"};
  TTreeReaderArray<unsigned int> vtxAssocPart_end = {treereader, "CentralTrackVertices.associatedParticles_end"};
  TTreeReaderArray<int> vtxAssocPart_index = {treereader, "_CentralTrackVertices_associatedParticles.index"};
  
  int nevents = 0;
  while(treereader.Next())
    {
      if(nevents%1000==0) printf("\n[i] New event %d\n",nevents);
      //if(nevents==20) break;

      // find MC primary vertex
      int nMCPart = mcPartMass.GetSize();
      TVector3 vertex_mc(-999., -999., -999.);
      for(int imc=0; imc<nMCPart; imc++)
	{
	  if(mcPartGenStatus[imc] == 4 && mcPartPdg[imc] == 11)
	    {
	      vertex_mc.SetXYZ(mcEndPointX[imc], mcEndPointY[imc], mcEndPointZ[imc]);
	      break;
	    }
	}
      hEventStat->Fill(0.5);
      hMcVtxX->Fill(vertex_mc.x());
      hMcVtxY->Fill(vertex_mc.y());
      hMcVtxZ->Fill(vertex_mc.z());

      // get RC primary vertex
      TVector3 vertex_rc(-999., -999., -999.);
      if(prim_vtx_index.GetSize()>0)
	{
	  int rc_vtx_index = prim_vtx_index[0];
	  vertex_rc.SetXYZ(CTVx[rc_vtx_index], CTVy[rc_vtx_index], CTVz[rc_vtx_index]);
	}

      // map MC and RC particles
      int nAssoc = assocChRecID.GetSize();
      map<int, int> assoc_map_to_rc;
      map<int, int> assoc_map_to_mc;
      for(int j=0; j<nAssoc; j++)
	{
	  assoc_map_to_rc[assocChSimID[j]] = assocChRecID[j];
	  assoc_map_to_mc[assocChRecID[j]] = assocChSimID[j];
	}

      for(int j=0; j<nAssoc; j++)
	{
	  if(assoc_map_to_rc.count(assocChSimID[j])>1) cout << "Found mc " << assocChSimID[j] << endl;
	  if(assoc_map_to_mc.count(assocChRecID[j])>1) cout << "Found rc " << assocChRecID[j] << endl;
	}

      // Loop over primary particles
      int nMcPart = 0;
      for(int imc=0; imc<nMCPart; imc++)
	{
	  if(mcPartGenStatus[imc] == 1 && mcPartCharge[imc] != 0)
	    {
	      double dist = sqrt( pow(mcPartVx[imc]-vertex_mc.x(),2) + pow(mcPartVy[imc]-vertex_mc.y(),2) + pow(mcPartVz[imc]-vertex_mc.z(),2));      
	      if(dist < 1e-4)
		{
		  // count charged particles within |eta| < 3.5
		  TVector3 mc_mom(mcMomPx[imc], mcMomPy[imc], mcMomPz[imc]);
		  double mcEta = mc_mom.PseudoRapidity();
		  if(fabs(mcEta) < 3.5) nMcPart++;
		}
	    }
	}
      hMcMult->Fill(nMcPart);
      
      for(int imc=0; imc<nMCPart; imc++)
	{
	  if(mcPartGenStatus[imc] == 1 && mcPartCharge[imc] != 0)
	    {
	      double dist = sqrt( pow(mcPartVx[imc]-vertex_mc.x(),2) + pow(mcPartVy[imc]-vertex_mc.y(),2) + pow(mcPartVz[imc]-vertex_mc.z(),2));      
	      if(dist < 1e-4)
		{		  
		  // check if the MC particle is reconstructed
		  int rc_index = -1;
		  if(assoc_map_to_rc.find(imc) != assoc_map_to_rc.end()) rc_index = assoc_map_to_rc[imc];

		  if(rc_index>=0)
		    {
		      TVector3 dcaToVtx = getDcaToVtx(rc_index, vertex_rc);
		      
		      int ip = -1;
		      if(fabs(mcPartPdg[imc]) == 211) ip = 0;
		      if(fabs(mcPartPdg[imc]) == 321) ip = 1;
		      if(fabs(mcPartPdg[imc]) == 2212) ip = 2;
		      if(ip>=0)
			{
			  TVector3 mom(rcMomPx[rc_index], rcMomPy[rc_index], rcMomPz[rc_index]);
			  if(ip<2)
			    {
			      hRcPrimPartLocaToRCVtx[ip]->Fill(mom.Pt(), mom.Eta(), dcaToVtx.Perp());
			      hRcPrimPartLocbToRCVtx[ip]->Fill(mom.Pt(), mom.Eta(), dcaToVtx.z());
			    }

			  double fill1[] = {mom.Pt(), mom.Eta(), dcaToVtx.x(), nMcPart*1.};
			  double fill2[] = {mom.Pt(), mom.Eta(), dcaToVtx.y(), nMcPart*1.};
			  double fill3[] = {mom.Pt(), mom.Eta(), dcaToVtx.z(), nMcPart*1.};
			  hPrimTrkDcaToRCVtx[ip][0]->Fill(fill1);
			  hPrimTrkDcaToRCVtx[ip][1]->Fill(fill2);
			  hPrimTrkDcaToRCVtx[ip][2]->Fill(fill3);
			}
		    }
		}
            }
	}

      // look for D0
      bool hasD0 = false;
      vector<int> mc_index_D0_pi;
      vector<int> mc_index_D0_k;
      mc_index_D0_pi.clear();
      mc_index_D0_k.clear();
      
      for(int imc=0; imc<nMCPart; imc++)
	{
	  if(fabs(mcPartPdg[imc]) != 421) continue;
	  hEventStat->Fill(1.5);
	  
	  int nDuaghters = mcPartDaughter_end[imc]-mcPartDaughter_begin[imc];
	  if(nDuaghters!=2) continue;

	  // find D0 that decay into pi+K
	  bool is_pik_decay = false;	  
	  int daug_index_1 = mcPartDaughter_index[mcPartDaughter_begin[imc]];
	  int daug_index_2 = mcPartDaughter_index[mcPartDaughter_begin[imc]+1];
	  int daug_pdg_1 = mcPartPdg[daug_index_1];
	  int daug_pdg_2 = mcPartPdg[daug_index_2];
	  if( (fabs(daug_pdg_1)==321 && fabs(daug_pdg_2)==211) || (fabs(daug_pdg_1)==211 && fabs(daug_pdg_2)==321) )
	    {
	      is_pik_decay = true;
	    }
	  if(!is_pik_decay) continue;
	  if(fabs(daug_pdg_1)==211)
	    {
	      mc_index_D0_pi.push_back(daug_index_1);
	      mc_index_D0_k.push_back(daug_index_2);
	    }
	  else
	    {
	      mc_index_D0_pi.push_back(daug_index_2);
	      mc_index_D0_k.push_back(daug_index_1);
	    }
	  hasD0 = true;
	  hEventStat->Fill(2.5);

	  // D0 kinematics
	  TLorentzVector mc_mom_vec;
	  mc_mom_vec.SetXYZM(mcMomPx[imc], mcMomPy[imc], mcMomPz[imc], mcPartMass[imc]);
	  double mcRap = mc_mom_vec.Rapidity();
	  double mcPt = mc_mom_vec.Pt();
	  hMCD0PtRap->Fill(mcRap, mcPt);

	  // decay dauther kinematics
	  for(int ip = 0; ip<2; ip++)
	    {
	      int mc_part_index;
	      if(ip==0) mc_part_index = mc_index_D0_pi[mc_index_D0_pi.size()-1];
	      if(ip==1) mc_part_index = mc_index_D0_k[mc_index_D0_k.size()-1];
	      
	      TLorentzVector mc_part_vec;
	      mc_part_vec.SetXYZM(mcMomPx[mc_part_index], mcMomPy[mc_part_index], mcMomPz[mc_part_index], mcPartMass[mc_part_index]);
	      if(ip==0) hMcPiPtEta->Fill(mc_part_vec.Eta(), mc_part_vec.Pt());
	      if(ip==1) hMcKPtEta->Fill(mc_part_vec.Eta(), mc_part_vec.Pt());
		  
	      int rc_part_index = -1;
	      if(assoc_map_to_rc.find(mc_part_index) != assoc_map_to_rc.end()) rc_part_index = assoc_map_to_rc[mc_part_index];
	      if(rc_part_index>=0)
		{
		  TVector3 dcaToVtx = getDcaToVtx(rc_part_index, vertex_rc);
		  
		  TVector3 mom(rcMomPx[rc_part_index], rcMomPy[rc_part_index], rcMomPz[rc_part_index]);
		  hRcSecPartLocaToRCVtx[ip]->Fill(mom.Pt(), mom.Eta(), dcaToVtx.Pt());
		  hRcSecPartLocbToRCVtx[ip]->Fill(mom.Pt(), mom.Eta(), dcaToVtx.z());
		  
		  //printf("Sec %d: (%2.4f, %2.4f, %2.4f), mcStartPoint = (%2.4f, %2.4f, %2.4f)\n", rc_part_index, pos.x(), pos.y(), pos.z(), mcPartVx[mc_part_index], mcPartVy[mc_part_index], mcPartVz[mc_part_index]);
		}
	    }
	}

      // Get reconstructed pions and kaons
      hNRecoVtx->Fill(CTVx.GetSize());
      const int pid_mode = 1; // 0 - truth; 1 - realistic
      vector<unsigned int> pi_index;
      vector<unsigned int> k_index;
      pi_index.clear();
      k_index.clear();
      for(unsigned int rc_index=0; rc_index<rcMomPx.GetSize(); rc_index++)
	{	  
	  if(pid_mode==0)
	    {
	      int iSimPartID = -1;
	      if(assoc_map_to_mc.find(rc_index) != assoc_map_to_mc.end()) iSimPartID = assoc_map_to_mc[rc_index];
	      if(iSimPartID>=0)
		{
		  if(fabs(mcPartPdg[iSimPartID]) == 211) pi_index.push_back(rc_index);
		  if(fabs(mcPartPdg[iSimPartID]) == 321) k_index.push_back(rc_index);
		}
	    }
	  else if(pid_mode==1)
	    {
	      if(fabs(rcPdg[rc_index]) == 211) pi_index.push_back(rc_index);
	      if(fabs(rcPdg[rc_index]) == 321) k_index.push_back(rc_index);
	    }
	}

      // pair pion and kaon
      for(unsigned int i=0; i<pi_index.size(); i++)
	{
	  TVector3 dcaToVtx = getDcaToVtx(pi_index[i], vertex_rc);
	  for(unsigned int j=0; j<k_index.size(); j++)
	    {
	      TVector3 dcaToVtx2 = getDcaToVtx(k_index[j], vertex_rc);
	      if(rcCharge[pi_index[i]]*rcCharge[k_index[j]]<0)
		{
		  //printf("[i] Check pair (%d, %d)\n", pi_index[i], k_index[j]);
		  // -- only look at unlike-sign pi+k pair
		  bool is_D0_pik = false;
		  int mc_index_pi = -1, mc_index_k = -1;
		  if(assoc_map_to_mc.find(pi_index[i]) != assoc_map_to_mc.end()) mc_index_pi = assoc_map_to_mc[pi_index[i]];
		  if(assoc_map_to_mc.find(k_index[j])  != assoc_map_to_mc.end()) mc_index_k  = assoc_map_to_mc[k_index[j]];

		  for(unsigned int k=0; k<mc_index_D0_pi.size(); k++)
		    {
		      if(mc_index_pi==mc_index_D0_pi[k] && mc_index_k==mc_index_D0_k[k])
			{
			  is_D0_pik = true;
			  break;
			}
		    }

		  float dcaDaughters, cosTheta, decayLength, V0DcaToVtx;
		  TLorentzVector parent = getPairParent(pi_index[i], k_index[j], vertex_rc, dcaDaughters, cosTheta, decayLength, V0DcaToVtx);
				  
		  if(is_D0_pik)
		    {
		      h3PairDca12[0]->Fill(parent.Pt(), parent.Rapidity(), dcaDaughters);
		      h3PairCosTheta[0]->Fill(parent.Pt(), parent.Rapidity(), cosTheta);
		      h3PairDca[0]->Fill(parent.Pt(), parent.Rapidity(), V0DcaToVtx);
		      h3PairDecayLength[0]->Fill(parent.Pt(), parent.Rapidity(), decayLength);
		      //printf("Signal: dca12 = %2.4f, cosTheta = %2.4f, D0dca = %2.4f, decay = %2.4f\n", dcaDaughters, cosTheta, V0DcaToVtx, decayLength);
		      h3InvMass[0][0]->Fill(parent.Pt(), parent.Rapidity(), parent.M());
		    }
		  else
		    {
		      h3PairDca12[1]->Fill(parent.Pt(), parent.Rapidity(), dcaDaughters);
		      h3PairCosTheta[1]->Fill(parent.Pt(), parent.Rapidity(), cosTheta);
		      h3PairDca[1]->Fill(parent.Pt(), parent.Rapidity(), V0DcaToVtx);
		      h3PairDecayLength[1]->Fill(parent.Pt(), parent.Rapidity(), decayLength);
			  
		      //printf("Bkg: dca12 = %2.4f, cosTheta = %2.4f, D0dca = %2.4f, decay = %2.4f\n", dcaDaughters, cosTheta, V0DcaToVtx, decayLength);
		      h3InvMass[1][0]->Fill(parent.Pt(), parent.Rapidity(), parent.M());
		    }

		  if(dcaToVtx.Perp() >= 0.02 && dcaToVtx2.Perp() >= 0.02 &&
		     dcaDaughters < 0.07 && cosTheta > 0.95 && decayLength > 0.05 && V0DcaToVtx < 0.1)
		    {
		      if(is_D0_pik)
			{
			  h3InvMass[0][1]->Fill(parent.Pt(), parent.Rapidity(), parent.M());
			}
		      else
			{
			  h3InvMass[1][1]->Fill(parent.Pt(), parent.Rapidity(), parent.M());
			}
		    }
		}
	    }
	}
	      
      nevents++;
    }

  TFile *outfile = new TFile(outname.Data(), "recreate");

  hEventStat->Write();
  hMcMult->Write();
  hMcVtxX->Write();
  hMcVtxY->Write();
  hMcVtxZ->Write();
  
  hD0DecayVxVy->Write();
  hD0DecayVrVz->Write();
  
  hMCD0PtRap->Write();

  hMcPiPtEta->Write();
  hMcPiPtEtaReco->Write();
  hMcKPtEta->Write();
  hMcKPtEtaReco->Write();
  
  hNRecoVtx->Write();

  for(int ip=0; ip<2; ip++)
    {
      hRcSecPartLocaToRCVtx[ip]->Write();
      hRcSecPartLocbToRCVtx[ip]->Write();
      hRcPrimPartLocaToRCVtx[ip]->Write();
      hRcPrimPartLocbToRCVtx[ip]->Write();
    }

  for(int i=0; i<3; i++)
    {
      for(int j=0; j<3; j++)
	{
	  hPrimTrkDcaToRCVtx[i][j]->Write();
	}
    }

  for(int i=0; i<2; i++)
    {
      h3PairDca12[i]->Write();
      h3PairCosTheta[i]->Write();
      h3PairDca[i]->Write();
      h3PairDecayLength[i]->Write();
    }

  for(int i=0; i<2; i++)
    {
      for(int j=0; j<2; j++)
	{
	  h3InvMass[i][j]->Write();
	}
    }
  
  
  outfile->Close();

}

//======================================
TVector3 getDcaToVtx(const int index, TVector3 vtx)
{
  //printf("check %d: (%2.4f, %2.4f, %2.4f) =? (%2.4f, %2.4f, %2.4f)\n", index, mom.x(), mom.y(), mom.z(), rcMomPx2->At(index), rcMomPy2->At(index), rcMomPz2->At(index));

  TVector3 pos(rcTrkLoca2->At(index) * sin(rcTrkPhi2->At(index)) * -1 * millimeter, rcTrkLoca2->At(index) * cos(rcTrkPhi2->At(index)) * millimeter, rcTrkLocb2->At(index) * millimeter);
  TVector3 mom(rcMomPx2->At(index), rcMomPy2->At(index), rcMomPz2->At(index));
   
  StPhysicalHelix pHelix(mom, pos, bField * tesla, rcCharge2->At(index));

  TVector3 vtx_tmp;
  vtx_tmp.SetXYZ(vtx.x()*millimeter, vtx.y()*millimeter, vtx.z()*millimeter);
  
  pHelix.moveOrigin(pHelix.pathLength(vtx_tmp));
  TVector3 dcaToVtx = pHelix.origin() - vtx_tmp;

  dcaToVtx.SetXYZ(dcaToVtx.x()/millimeter, dcaToVtx.y()/millimeter, dcaToVtx.z()/millimeter);
  
  return dcaToVtx;
}

//======================================
TLorentzVector getPairParent(const int index1, const int index2, TVector3 vtx,
			     float &dcaDaughters, float &cosTheta, float &decayLength, float &V0DcaToVtx)
{
  // -- get helix
  TVector3 pos1(rcTrkLoca2->At(index1) * sin(rcTrkPhi2->At(index1)) * -1 * millimeter, rcTrkLoca2->At(index1) * cos(rcTrkPhi2->At(index1)) * millimeter, rcTrkLocb2->At(index1) * millimeter);
  TVector3 mom1(rcMomPx2->At(index1), rcMomPy2->At(index1), rcMomPz2->At(index1));

  TVector3 pos2(rcTrkLoca2->At(index2) * sin(rcTrkPhi2->At(index2)) * -1 * millimeter, rcTrkLoca2->At(index2) * cos(rcTrkPhi2->At(index2)) * millimeter, rcTrkLocb2->At(index2) * millimeter);
  TVector3 mom2(rcMomPx2->At(index2), rcMomPy2->At(index2), rcMomPz2->At(index2));

  float charge1 = rcCharge2->At(index1);
  float charge2 = rcCharge2->At(index2);
  
  StPhysicalHelix p1Helix(mom1, pos1, bField * tesla, charge1);
  StPhysicalHelix p2Helix(mom2, pos2, bField * tesla, charge2);

  TVector3 vtx_tmp;
  vtx_tmp.SetXYZ(vtx.x()*millimeter, vtx.y()*millimeter, vtx.z()*millimeter);

  // -- use straight lines approximation to get point of DCA of particle1-particle2 pair
  // -- this is used to save time in heavy-ion collisions 
  // -- move origin
  // p1Helix.moveOrigin(p1Helix.pathLength(vtx));
  // p2Helix.moveOrigin(p2Helix.pathLength(vtx));
		  
  // TVector3 const p1Mom = p1Helix.momentum(bField * tesla);
  // TVector3 const p2Mom = p2Helix.momentum(bField * tesla);

  // StPhysicalHelix const p1StraightLine(p1Mom, p1Helix.origin(), 0, charge1);
  // StPhysicalHelix const p2StraightLine(p2Mom, p2Helix.origin(), 0, charge2);
  
  // pair<double, double> const ss = p1StraightLine.pathLengths(p2StraightLine);
  // TVector3 const p1AtDcaToP2 = p1StraightLine.at(ss.first);
  // TVector3 const p2AtDcaToP1 = p2StraightLine.at(ss.second);

  // -- Full calculation 
  pair<double, double> const ss = p1Helix.pathLengths(p2Helix);
  TVector3 const p1AtDcaToP2 = p1Helix.at(ss.first);
  TVector3 const p2AtDcaToP1 = p2Helix.at(ss.second);
  // cout << ss.first << "  " << ss.second << endl;
  // printf("p1AtDcaToP2 origin = (%2.4f, %2.4f, %2.4f)\n", p1AtDcaToP2.x(), p1AtDcaToP2.y(), p1AtDcaToP2.z());
  // printf("p2AtDcaToP1 origin = (%2.4f, %2.4f, %2.4f)\n", p2AtDcaToP1.x(), p2AtDcaToP1.y(), p2AtDcaToP1.z());
  
  // -- calculate DCA of particle1 to particle2 at their DCA
  dcaDaughters = (p1AtDcaToP2 - p2AtDcaToP1).Mag()/millimeter;
	
  // -- calculate Lorentz vector of particle1-particle2 pair
  TVector3 const p1MomAtDca = p1Helix.momentumAt(ss.first,  bField * tesla);
  TVector3 const p2MomAtDca = p2Helix.momentumAt(ss.second, bField * tesla);
  
  TLorentzVector p1FourMom(p1MomAtDca, sqrt(p1MomAtDca.Mag2()+gPionMass*gPionMass));
  TLorentzVector p2FourMom(p2MomAtDca, sqrt(p2MomAtDca.Mag2()+gKaonMass*gKaonMass));
  
  TLorentzVector parent = p1FourMom + p2FourMom;
	
  // -- calculate cosThetaStar
  // TLorentzVector pairFourMomReverse(-parent.Px(), -parent.Py(), -parent.Pz(), parent.E());
  // TLorentzVector p1FourMomStar = p1FourMom;
  // p1FourMomStar.Boost(pairFourMomReverse.Vect());
  // float cosThetaStar = std::cos(p1FourMomStar.Vect().Angle(parent.Vect()));
	
  // -- calculate decay vertex (secondary or tertiary)
  TVector3 decayVertex = (p1AtDcaToP2 + p2AtDcaToP1) * 0.5 ;
	
  // -- calculate pointing angle and decay length with respect to primary vertex
  //    if decay vertex is a tertiary vertex
  //    -> only rough estimate -> needs to be updated after secondary vertex is found
  TVector3 vtxToV0 = decayVertex - vtx_tmp;
  float pointingAngle = vtxToV0.Angle(parent.Vect());
  cosTheta = std::cos(pointingAngle);
  decayLength = vtxToV0.Mag()/millimeter;

  // -- calculate V0 DCA to primary vertex
  V0DcaToVtx = decayLength * std::sin(pointingAngle);
    
  //TVector3 dcaToVtx = getDcaToVtx(parent.Vect(), decayVertex, 0, vtx);
  //V0DcaToVtx = dcaToVtx.Mag();

  return parent;
}

