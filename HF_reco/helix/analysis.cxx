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

TVector3 getDcaToVtx(TVector3 mom, TVector3 pos, float charge, TVector3 vtx);

TLorentzVector getPairParent(TVector3 mom1, TVector3 pos1, float charge1, TVector3 mom2, TVector3 pos2, float charge2, TVector3 vtx,
			     float &dcaDaughters, float &cosTheta, float &decayLength, float &V0DcaToVtx);

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

  const char* part_name[2] = {"Pi", "K"};
  const char* part_title[2] = {"#pi", "K"};
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

  TH3F *h3PairDca12[2];
  TH3F *h3PairCosTheta[2];
  TH3F *h3PairDca[2];
  TH3F *h3PairDecayLength[2];
  const char* pair_name[2] = {"signal", "bkg"};
  const char* pair_title[2] = {"Signal", "Background"};
  for(int i=0; i<2; i++)
    {
      h3PairDca12[i] = new TH3F(Form("h3PairDca12_%s", pair_name[i]), Form("%s pair DCA_{12};p_{T} (GeV/c);#eta;DCA_{12} (mm)", pair_title[i]), 100, 0, 10, 20, -5, 5, 100, 0, 0.2);

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
  TTreeReaderArray<float> mcMomPx = {treereader, "MCParticles.momentum.x"};
  TTreeReaderArray<float> mcMomPy = {treereader, "MCParticles.momentum.y"};
  TTreeReaderArray<float> mcMomPz = {treereader, "MCParticles.momentum.z"};
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

  TTreeReaderArray<float> rcTrkLoca = {treereader, "CentralCKFTrackParameters.loc.a"};
  TTreeReaderArray<float> rcTrkLocb = {treereader, "CentralCKFTrackParameters.loc.b"};
  TTreeReaderArray<float> rcTrkqOverP = {treereader, "CentralCKFTrackParameters.qOverP"};
  TTreeReaderArray<float> rcTrkTheta = {treereader, "CentralCKFTrackParameters.theta"};
  TTreeReaderArray<float> rcTrkPhi = {treereader, "CentralCKFTrackParameters.phi"};

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
      //printf("\n+++++ New Event %d +++++\n", nevents);

      // find MC primary vertex
      int nMCPart = mcPartMass.GetSize();
      TVector3 vertex_mc(-999., -999., -999.);
      for(int imc=0; imc<nMCPart; imc++)
	{
	  if(mcPartGenStatus[imc] == 4 && mcPartPdg[imc] == 11)
	    {
	      vertex_mc.SetXYZ(mcEndPointX[imc], mcEndPointY[imc], mcEndPointZ[imc]);
	      //printf("[i] Primary vertex (x, y, z) = (%2.4f, %2.4f, %2.4f)\n", vertex_mc.x(), vertex_mc.y(), vertex_mc.z());
	      break;
	    }
	}

      // get RC primary vertex
      TVector3 vertex_rc(-999., -999., -999.);
      if(prim_vtx_index.GetSize()>0)
	{
	  int rc_vtx_index = prim_vtx_index[0];
	  vertex_rc.SetXYZ(CTVx[rc_vtx_index], CTVy[rc_vtx_index], CTVz[rc_vtx_index]);
	}

      // map MC and RC particles
      int nAssoc = assocChRecID.GetSize();
      map<int, int> assoc_map;
      for(int j=0; j<nAssoc; j++)
	{
	  assoc_map[assocChSimID[j]] = assocChRecID[j];
	}

      // Loop over primary particles
      for(int imc=0; imc<nMCPart; imc++)
	{
	  if(mcPartGenStatus[imc] == 1 && mcPartCharge[imc] != 0 && (fabs(mcPartPdg[imc]) == 211 || fabs(mcPartPdg[imc]) == 321))
	    {
	      double dist = sqrt( pow(mcPartVx[imc]-vertex_mc.x(),2) + pow(mcPartVy[imc]-vertex_mc.y(),2) + pow(mcPartVz[imc]-vertex_mc.z(),2));      
	      if(dist < 1e-4)
		{
		  // check if the MC particle is reconstructed
		  int rc_index = -1;
		  if(assoc_map.find(imc) != assoc_map.end()) rc_index = assoc_map[imc];

		  if(rc_index>=0)
		    {
		      TVector3 pos(rcPosx[rc_index], rcPosy[rc_index], rcPosz[rc_index]);
		      TVector3 mom(rcMomPx[rc_index], rcMomPy[rc_index], rcMomPz[rc_index]);
		      TVector3 dcaToVtx = getDcaToVtx(mom, pos, rcCharge[rc_index], vertex_rc);
		      
		      int ip = 0;
		      if(fabs(mcPartPdg[imc]) == 321) ip = 1;
		      hRcPrimPartLocaToRCVtx[ip]->Fill(mom.Pt(), mom.Eta(), dcaToVtx.Pt());
		      hRcPrimPartLocbToRCVtx[ip]->Fill(mom.Pt(), mom.Eta(), dcaToVtx.z());
		    }
		}
            }
	}

      // look for D0
      bool hasD0 = false;
      int mc_index_D0_pi = -1;
      int mc_index_D0_k  = -1;
      for(int imc=0; imc<nMCPart; imc++)
	{
	  if(fabs(mcPartPdg[imc]) == 421)
	    {
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
	      hasD0 = true;
	      hEventStat->Fill(2.5);

	      // D0 kinematics
	      TLorentzVector mc_mom_vec;
	      mc_mom_vec.SetXYZM(mcMomPx[imc], mcMomPy[imc], mcMomPz[imc], mcPartMass[imc]);
	      double mcRap = mc_mom_vec.Rapidity();
	      double mcPt = mc_mom_vec.Pt();
	      hMCD0PtRap->Fill(mcRap, mcPt);

	      // decay dauther kinematics
	      mc_index_D0_pi = fabs(daug_pdg_1)==211 ? daug_index_1 : daug_index_2;
	      mc_index_D0_k = fabs(daug_pdg_1)==321 ? daug_index_1 : daug_index_2;
	      for(int ip = 0; ip<2; ip++)
		{
		  int mc_part_index;
		  if(ip==0) mc_part_index = mc_index_D0_pi;
		  if(ip==1) mc_part_index = mc_index_D0_k;

		  TLorentzVector mc_part_vec;
		  mc_part_vec.SetXYZM(mcMomPx[mc_part_index], mcMomPy[mc_part_index], mcMomPz[mc_part_index], mcPartMass[mc_part_index]);
		  if(ip==0) hMcPiPtEta->Fill(mc_part_vec.Eta(), mc_part_vec.Pt());
		  if(ip==1) hMcKPtEta->Fill(mc_part_vec.Eta(), mc_part_vec.Pt());

		  int rc_part_index = -1;
		  if(assoc_map.find(mc_part_index) != assoc_map.end()) rc_part_index = assoc_map[mc_part_index];
		  if(rc_part_index>=0)
		    {
		      TVector3 pos(rcPosx[rc_part_index], rcPosy[rc_part_index], rcPosz[rc_part_index]);
		      TVector3 mom(rcMomPx[rc_part_index], rcMomPy[rc_part_index], rcMomPz[rc_part_index]);
		      TVector3 dcaToVtx = getDcaToVtx(mom, pos, rcCharge[rc_part_index], vertex_rc);

		      hRcSecPartLocaToRCVtx[ip]->Fill(mom.Pt(), mom.Eta(), dcaToVtx.Pt());
		      hRcSecPartLocbToRCVtx[ip]->Fill(mom.Pt(), mom.Eta(), dcaToVtx.z());
		    }
		}
	    }
	}
      hEventStat->Fill(0.5);
      hMcVtxX->Fill(vertex_mc.x());
      hMcVtxY->Fill(vertex_mc.y());
      hMcVtxZ->Fill(vertex_mc.z());

      // Get pair information
      vector<int> pi_index;
      vector<int> k_index;
      pi_index.clear();
      k_index.clear();
      for(int rc_index=0; rc_index<rcMomPx.GetSize(); rc_index++)
	{	  
	  int iSimPartID = -1;
	  for(int j=0; j<nAssoc; j++)
	    {
	      if(assocChRecID[j]==rc_index)
		{
		  iSimPartID = assocChSimID[j];
		  break;
		}
	    }
	  if(iSimPartID<0) continue;
	  if(fabs(mcPartPdg[iSimPartID]) == 211) pi_index.push_back(rc_index);
	  if(fabs(mcPartPdg[iSimPartID]) == 321) k_index.push_back(rc_index);
	}

      // pair pion and kaon
      for(int i=0; i<pi_index.size(); i++)
	{
	  TLorentzVector pi_mom_vec;
	  pi_mom_vec.SetXYZM(rcMomPx[pi_index[i]], rcMomPy[pi_index[i]], rcMomPz[pi_index[i]], gPionMass);
	  for(int j=0; j<k_index.size(); j++)
	    {
	      TLorentzVector k_mom_vec;
	      k_mom_vec.SetXYZM(rcMomPx[k_index[j]], rcMomPy[k_index[j]], rcMomPz[k_index[j]], gKaonMass);
	      if(rcCharge[pi_index[i]]*rcCharge[k_index[j]]<0)
		{

		  // -- only look at signal pi+k pair
		  bool is_D0_pik = true;
		  if(hasD0)
		    {
		      for(int k=0; k<nAssoc; k++)
			{
			  if(assocChRecID[k]==pi_index[i])
			    {
			      if(assocChSimID[k]!=mc_index_D0_pi)
				{
				  is_D0_pik = false;
				  break;
				}
			    }

			  if(assocChRecID[k]==k_index[j])
			    {
			      if(assocChSimID[k]!=mc_index_D0_k)
				{
				  is_D0_pik = false;
				  break;
				}
			    }
			}
		    }

		  /*
		  // -- get helix
		  TVector3 pos1(rcPosx[pi_index[i]], rcPosy[pi_index[i]], rcPosz[pi_index[i]]);
		  TVector3 mom1(rcMomPx[pi_index[i]], rcMomPy[pi_index[i]], rcMomPz[pi_index[i]]);
		  StPhysicalHelix p1Helix(mom1, pos1, bField * tesla, rcCharge[pi_index[i]]);

		  TVector3 pos2(rcPosx[k_index[j]], rcPosy[k_index[j]], rcPosz[k_index[j]]);
		  TVector3 mom2(rcMomPx[k_index[j]], rcMomPy[k_index[j]], rcMomPz[k_index[j]]);
		  StPhysicalHelix p2Helix(mom2, pos2, bField * tesla,  rcCharge[k_index[j]]);
		  
		  // -- use straight lines approximation to get point of DCA of particle1-particle2 pair
		  StPhysicalHelix const p1StraightLine(mom1, pos1, 0, rcCharge[pi_index[i]]);
		  StPhysicalHelix const p2StraightLine(mom2, pos2, 0, rcCharge[k_index[j]]);

		  //printf("p2Mom = (%2.4f, %2.4f, %2.4f)\n",mom2.x(), mom2.y(), mom2.z());
	
		  pair<double, double> const ss = p1StraightLine.pathLengths(p2StraightLine);
		  TVector3 const p1AtDcaToP2 = p1StraightLine.at(ss.first);
		  TVector3 const p2AtDcaToP1 = p2StraightLine.at(ss.second);
		  //cout << ss.first << "  " << ss.second << endl;
	
		  // -- calculate DCA of particle1 to particle2 at their DCA
		  float dcaDaughters = (p1AtDcaToP2 - p2AtDcaToP1).Mag();
	
		  // -- calculate Lorentz vector of particle1-particle2 pair
		  TVector3 const p1MomAtDca = p1Helix.momentumAt(ss.first,  bField * tesla);
		  TVector3 const p2MomAtDca = p2Helix.momentumAt(ss.second, bField * tesla);
	
		  TLorentzVector p1FourMom(p1MomAtDca, sqrt(p1MomAtDca.Mag2()+gPionMass*gPionMass));
		  TLorentzVector p2FourMom(p2MomAtDca, sqrt(p2MomAtDca.Mag2()+gKaonMass*gKaonMass));
	
		  TLorentzVector parent = p1FourMom + p2FourMom;


	
		  // -- calculate cosThetaStar
		  TLorentzVector pairFourMomReverse(-parent.Px(), -parent.Py(), -parent.Pz(), parent.E());
		  TLorentzVector p1FourMomStar = p1FourMom;
		  p1FourMomStar.Boost(pairFourMomReverse.Vect());
		  float cosThetaStar = std::cos(p1FourMomStar.Vect().Angle(parent.Vect()));
	
		  // -- calculate decay vertex (secondary or tertiary)
		  TVector3 decayVertex = (p1AtDcaToP2 + p2AtDcaToP1) * 0.5 ;
	
		  // -- calculate pointing angle and decay length with respect to primary vertex
		  //    if decay vertex is a tertiary vertex
		  //    -> only rough estimate -> needs to be updated after secondary vertex is found
		  TVector3 vtxToV0 = decayVertex - vertex_rc;
		  float pointingAngle = vtxToV0.Angle(parent.Vect());
		  float cosTheta = std::cos(pointingAngle);
		  float decayLength = vtxToV0.Mag();

		  // -- calculate V0 DCA to primary vertex
		  TVector3 dcaToVtx = getDcaToVtx(parent.Vect(), decayVertex, 0, vertex_rc);
		  float V0Dca = dcaToVtx.Mag();
		  */

		  TVector3 pos1(rcPosx[pi_index[i]], rcPosy[pi_index[i]], rcPosz[pi_index[i]]);
		  TVector3 mom1(rcMomPx[pi_index[i]], rcMomPy[pi_index[i]], rcMomPz[pi_index[i]]);

		  TVector3 pos2(rcPosx[k_index[j]], rcPosy[k_index[j]], rcPosz[k_index[j]]);
		  TVector3 mom2(rcMomPx[k_index[j]], rcMomPy[k_index[j]], rcMomPz[k_index[j]]);

		  float dcaDaughters, cosTheta, decayLength, V0DcaToVtx; 
		  TLorentzVector parent = getPairParent(mom1, pos1, rcCharge[pi_index[i]], mom2, pos2, rcCharge[k_index[j]], vertex_rc,
							dcaDaughters, cosTheta, decayLength, V0DcaToVtx);

		  if(hasD0)
		    {
		      if(is_D0_pik)
			{
			  //cout << p1FourMom.M() << "  " << p2FourMom.M() << "  " << parent.M() << endl;
			  h3PairDca12[0]->Fill(parent.Pt(), parent.Rapidity(), dcaDaughters);
			  h3PairCosTheta[0]->Fill(parent.Pt(), parent.Rapidity(), cosTheta);
			  h3PairDca[0]->Fill(parent.Pt(), parent.Rapidity(), V0DcaToVtx);
			  h3PairDecayLength[0]->Fill(parent.Pt(), parent.Rapidity(), decayLength);
			  //printf("Signal: dca12 = %2.4f, cosTheta = %2.4f, D0dca = %2.4f, decay = %2.4f\n", dcaDaughters, cosTheta, V0DcaToVtx, decayLength);
			}
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
		}
	    }
	}



      // loop over reconstructed particles
      vector<int> pi_dca_index;
      vector<int> k_dca_index;
      pi_dca_index.clear();
      k_dca_index.clear();
      for(int rc_index=0; rc_index<rcMomPx.GetSize(); rc_index++)
	{
	  int iSimPartID = -1;
	  for(int j=0; j<nAssoc; j++)
	    {
	      if(assocChRecID[j]==rc_index)
		{
		  iSimPartID = assocChSimID[j];
		  break;
		}
	    }
	  if(iSimPartID<0) continue;
	  if(fabs(mcPartPdg[iSimPartID]) != 211 && fabs(mcPartPdg[iSimPartID]) != 321) continue;

	  TVector3 pos(rcPosx[rc_index], rcPosy[rc_index], rcPosz[rc_index]);
	  TVector3 mom(rcMomPx[rc_index], rcMomPy[rc_index], rcMomPz[rc_index]);
	  TVector3 dcaToVtx = getDcaToVtx(mom, pos, rcCharge[rc_index], vertex_rc);

	  if(dcaToVtx.Pt() > 0.02)
	    {
	      if(fabs(mcPartPdg[iSimPartID]) == 211) pi_dca_index.push_back(rc_index);
	      if(fabs(mcPartPdg[iSimPartID]) == 321) k_dca_index.push_back(rc_index);
	    }
	}
      // pair pion and kaon
      for(int i=0; i<pi_dca_index.size(); i++)
	{
	  for(int j=0; j<k_dca_index.size(); j++)
	    {
	      TVector3 pos1(rcPosx[pi_dca_index[i]], rcPosy[pi_dca_index[i]], rcPosz[pi_dca_index[i]]);
	      TVector3 mom1(rcMomPx[pi_dca_index[i]], rcMomPy[pi_dca_index[i]], rcMomPz[pi_dca_index[i]]);

	      TVector3 pos2(rcPosx[k_dca_index[j]], rcPosy[k_dca_index[j]], rcPosz[k_dca_index[j]]);
	      TVector3 mom2(rcMomPx[k_dca_index[j]], rcMomPy[k_dca_index[j]], rcMomPz[k_dca_index[j]]);

	      float dcaDaughters, cosTheta, decayLength, V0DcaToVtx; 
	      TLorentzVector parent = getPairParent(mom1, pos1, rcCharge[pi_dca_index[i]], mom2, pos2, rcCharge[k_dca_index[j]], vertex_rc,
						    dcaDaughters, cosTheta, decayLength, V0DcaToVtx);

	      if(dcaDaughters < 0.002 && cosTheta > 0.95 && decayLength > 0.05 && V0DcaToVtx < 0.1)
		{
		  if(hasD0)
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
	      
      nevents++;
    }

  TFile *outfile = new TFile(outname.Data(), "recreate");

  hEventStat->Write();
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


TVector3 getDcaToVtx(TVector3 mom, TVector3 pos, float charge, TVector3 vtx)
{
  StPhysicalHelix pHelix(mom, pos, bField * tesla, charge);
  pHelix.moveOrigin(pHelix.pathLength(vtx));
  TVector3 dcaToVtx = pHelix.origin() - vtx;
  return dcaToVtx;
}

TLorentzVector getPairParent(TVector3 mom1, TVector3 pos1, float charge1, TVector3 mom2, TVector3 pos2, float charge2, TVector3 vtx,
			     float &dcaDaughters, float &cosTheta, float &decayLength, float &V0DcaToVtx)
{
  // -- get helix
  StPhysicalHelix p1Helix(mom1, pos1, bField * tesla, charge1);
  StPhysicalHelix p2Helix(mom2, pos2, bField * tesla, charge2);

  // -- move origin
  p1Helix.moveOrigin(p1Helix.pathLength(vtx));
  p2Helix.moveOrigin(p2Helix.pathLength(vtx));
		  
  // -- use straight lines approximation to get point of DCA of particle1-particle2 pair
  TVector3 const p1Mom = p1Helix.momentum(bField * tesla);
  TVector3 const p2Mom = p2Helix.momentum(bField * tesla);
  StPhysicalHelix const p1StraightLine(p1Mom, p1Helix.origin(), 0, charge1);
  StPhysicalHelix const p2StraightLine(p2Mom, p2Helix.origin(), 0, charge2);
  
  pair<double, double> const ss = p1StraightLine.pathLengths(p2StraightLine);
  TVector3 const p1AtDcaToP2 = p1StraightLine.at(ss.first);
  TVector3 const p2AtDcaToP1 = p2StraightLine.at(ss.second);
  //cout << ss.first << "  " << ss.second << endl;
  
  // -- calculate DCA of particle1 to particle2 at their DCA
  dcaDaughters = (p1AtDcaToP2 - p2AtDcaToP1).Mag();
	
  // -- calculate Lorentz vector of particle1-particle2 pair
  TVector3 const p1MomAtDca = p1Helix.momentumAt(ss.first,  bField * tesla);
  TVector3 const p2MomAtDca = p2Helix.momentumAt(ss.second, bField * tesla);
  
  TLorentzVector p1FourMom(p1MomAtDca, sqrt(p1MomAtDca.Mag2()+gPionMass*gPionMass));
  TLorentzVector p2FourMom(p2MomAtDca, sqrt(p2MomAtDca.Mag2()+gKaonMass*gKaonMass));
  
  TLorentzVector parent = p1FourMom + p2FourMom;
	
  // -- calculate cosThetaStar
  TLorentzVector pairFourMomReverse(-parent.Px(), -parent.Py(), -parent.Pz(), parent.E());
  TLorentzVector p1FourMomStar = p1FourMom;
  p1FourMomStar.Boost(pairFourMomReverse.Vect());
  float cosThetaStar = std::cos(p1FourMomStar.Vect().Angle(parent.Vect()));
	
  // -- calculate decay vertex (secondary or tertiary)
  TVector3 decayVertex = (p1AtDcaToP2 + p2AtDcaToP1) * 0.5 ;
	
  // -- calculate pointing angle and decay length with respect to primary vertex
  //    if decay vertex is a tertiary vertex
  //    -> only rough estimate -> needs to be updated after secondary vertex is found
  TVector3 vtxToV0 = decayVertex - vtx;
  float pointingAngle = vtxToV0.Angle(parent.Vect());
  cosTheta = std::cos(pointingAngle);
  decayLength = vtxToV0.Mag();

  // -- calculate V0 DCA to primary vertex
  V0DcaToVtx = decayLength * std::sin(pointingAngle);
    
  //TVector3 dcaToVtx = getDcaToVtx(parent.Vect(), decayVertex, 0, vtx);
  //V0DcaToVtx = dcaToVtx.Mag();

  return parent;
}

