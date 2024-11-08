#include <fstream>
#include <map>
#include "TROOT.h"
#include "TClass.h"
#include "TFile.h"
#include "TLine.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TProfile.h"
#include "TRandom.h"
#include <vector>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TLatex.h"

#include "DD4hep/Detector.h"
#include "DD4hep/DetElement.h"

#include "Acts/Plugins/DD4hep/ConvertDD4hepDetector.hpp"
#include "Acts/Plugins/DD4hep/DD4hepFieldAdapter.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Plugins/DD4hep/ConvertDD4hepDetector.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Vertexing/ImpactPointEstimator.hpp"

using namespace std;

const double gPionMass = 0.13957;
const double gKaonMass = 0.493677;

void ana_d0(TString listname = "file.list", const TString outname = "test.root")
{
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

  TH2F *hD0DecayVxVyReco = new TH2F("hD0DecayVxVyReco", "Reconstructed D^{0} decay vertex to primary vertex;#Deltav_{x} (mm);#Deltav_{y} (mm)", 400, -1-0.0025, 1-0.0025, 400, -1-0.0025, 1-0.0025);
  TH2F *hD0DecayVrVzReco = new TH2F("hD0DecayVrVzReco", "Reconstructed D^{0} decay vertex to primary vertex;#Deltav_{z} (mm);#Deltav_{r} (mm)", 100, -2, 2, 100, -0.2, 1.8);

  TH2F *hMCD0PtRap = new TH2F("hMCD0PtRap", "MC D^{0};y;p_{T} (GeV/c)", 20, -5, 5, 100, 0, 10);
  TH2F *hMCD0PtRapReco = new TH2F("hMCD0PtRapReco", "Reconstructed D^{0} decay vertex;y;p_{T} (GeV/c)", 20, -5, 5, 100, 0, 10);
  TH2F *hD0VtxDist = new TH2F("hD0VtxDist", "Reconstructed D^{0} vertex to true D^{0} vertex;#Deltav_{z} (mm);#Deltav_{r} (mm)", 100, -1, 1, 100, -1, 1);

  TH2F *hMcPiPtEta = new TH2F("hMcPiPtEta", "MC #pi from D^{0} decay;#eta^{MC};p_{T}^{MC} (GeV/c)", 20, -5, 5, 100, 0, 10);
  TH2F *hMcPiPtEtaReco = new TH2F("hMcPiPtEtaReco", "RC #pi from D^{0} decay;#eta^{MC};p_{T}^{MC} (GeV/c)", 20, -5, 5, 100, 0, 10);

  TH2F *hMcKPtEta = new TH2F("hMcKPtEta", "MC K from D^{0} decay;#eta^{MC};p_{T}^{MC} (GeV/c)", 20, -5, 5, 100, 0, 10);
  TH2F *hMcKPtEtaReco = new TH2F("hMcKPtEtaReco", "RC K from D^{0} decay;#eta^{MC};p_{T}^{MC} (GeV/c)", 20, -5, 5, 100, 0, 10);

  TH1F *hNRecoVtx = new TH1F("hNRecoVtx", "Number of reconstructed vertices;N", 10, 0, 10);

  const char* part_name[2] = {"Pi", "K"};
  const char* part_title[2] = {"#pi", "K"};
  TH3F *hRcSecPartLocaToMCVtx[2];
  TH3F *hRcSecPartLocbToMCVtx[2];
  TH3F *hRcPrimPartLocaToMCVtx[2];
  TH3F *hRcPrimPartLocbToMCVtx[2];
  TH3F *hRcSecPartLocaToRCVtx[2];
  TH3F *hRcSecPartLocbToRCVtx[2];
  TH3F *hRcPrimPartLocaToRCVtx[2];
  TH3F *hRcPrimPartLocbToRCVtx[2];
  for(int i=0; i<2; i++)
    {
      hRcSecPartLocaToMCVtx[i] = new TH3F(Form("hRcSec%sLocaToMCVtx",part_name[i]), Form("Loc.a distribution for D^{0} decayed %s;p_{T} (GeV/c);#eta;loc.a (mm)", part_title[i]), 100, 0, 10, 20, -5, 5, 100, 0, 1);
      hRcSecPartLocbToMCVtx[i] = new TH3F(Form("hRcSec%sLocbToMCVtx",part_name[i]), Form( "Loc.b distribution for D^{0} decayed %s;p_{T} (GeV/c);#eta;loc.b (mm)", part_title[i]), 100, 0, 10, 20, -5, 5, 100, -0.5, 0.5);
      hRcPrimPartLocaToMCVtx[i] = new TH3F(Form("hRcPrim%sLocaToMCVtx",part_name[i]), Form( "Loc.a distribution for primary %s;p_{T} (GeV/c);#eta;loc.a (mm)", part_title[i]), 100, 0, 10, 20, -5, 5, 100, 0, 1);
      hRcPrimPartLocbToMCVtx[i] = new TH3F(Form("hRcPrim%sLocbToMCVtx",part_name[i]), Form( "Loc.b distribution for primary %s;p_{T} (GeV/c);#eta;loc.b (mm)", part_title[i]), 100, 0, 10, 20, -5, 5, 100, -0.5, 0.5);

      hRcSecPartLocaToRCVtx[i] = new TH3F(Form("hRcSec%sLocaToRCVtx",part_name[i]), Form( "Loc.a distribution for D^{0} decayed %s;p_{T} (GeV/c);#eta;loc.a (mm)", part_title[i]), 100, 0, 10, 20, -5, 5, 100, 0, 1);
      hRcSecPartLocbToRCVtx[i] = new TH3F(Form("hRcSec%sLocbToRCVtx",part_name[i]), Form( "Loc.b distribution for D^{0} decayed %s;p_{T} (GeV/c);#eta;loc.b (mm)", part_title[i]), 100, 0, 10, 20, -5, 5, 100, -0.5, 0.5);
      hRcPrimPartLocaToRCVtx[i] = new TH3F(Form("hRcPrim%sLocaToRCVtx",part_name[i]), Form( "Loc.a distribution for primary %s;p_{T} (GeV/c);#eta;loc.a (mm)", part_title[i]), 100, 0, 10, 20, -5, 5, 100, 0, 1);
      hRcPrimPartLocbToRCVtx[i] = new TH3F(Form("hRcPrim%sLocbToRCVtx",part_name[i]), Form( "Loc.b distribution for primary %s;p_{T} (GeV/c);#eta;loc.b (mm)", part_title[i]), 100, 0, 10, 20, -5, 5, 100, -0.5, 0.5);
    }

  // Invariant mass
  const char* pair_name[2] = {"all", "DCA"};
  TH3F *h3InvMass[2];
  for(int i=0; i<2; i++)
    {
      h3InvMass[i] = new TH3F(Form("h3InvMass_%s", pair_name[i]), "Invariant mass of unlike-sign #piK pairs;p_{T} (GeV/c);y;M_{#piK} (GeV/c^{2})", 100, 0, 10, 20, -5, 5, 100, 1.6, 2.0);
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

    // Load DD4Hep geometry
  dd4hep::Detector& detector = dd4hep::Detector::getInstance();
  detector.fromCompact("/opt/detector/epic-main/share/epic/epic_craterlake.xml");
  dd4hep::DetElement geometry = detector.world();

  // Convert DD4Hep geometry to tracking geometry
  Acts::GeometryContext trackingGeoCtx;
  auto logger = Acts::getDefaultLogger("DD4hepConversion", Acts::Logging::Level::INFO);
  Acts::BinningType bTypePhi = Acts::equidistant;
  Acts::BinningType bTypeR = Acts::equidistant;
  Acts::BinningType bTypeZ = Acts::equidistant;
  double layerEnvelopeR = Acts::UnitConstants::mm;
  double layerEnvelopeZ = Acts::UnitConstants::mm;
  double defaultLayerThickness = Acts::UnitConstants::fm;
  using Acts::sortDetElementsByID;

  std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry{nullptr};
  trackingGeometry = Acts::convertDD4hepDetector(geometry,*logger,bTypePhi,bTypeR,bTypeZ,layerEnvelopeR,layerEnvelopeZ,defaultLayerThickness,sortDetElementsByID,trackingGeoCtx);

  // Define Perigee surface at which reconstructed track parameters are set
  auto perigee = Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3(0,0,0));

  // Get Magnetic field context
  Acts::MagneticFieldContext fieldctx;
  std::shared_ptr<const Acts::DD4hepFieldAdapter> field_provider = std::make_shared<const Acts::DD4hepFieldAdapter>(detector.field());
  Acts::MagneticFieldProvider::Cache field_cache = field_provider->makeCache(fieldctx);

  // auto lookupResult = field_provider->getField(Acts::Vector3{0, 0, 0}, field_cache);
  // Acts::Vector3 fieldValue = *lookupResult;
  // cout << fieldValue.x() << "  " << fieldValue.y() << "  " << fieldValue.z() << endl;

  // dd4hep::Position pos(0,0,0);
  // auto field = detector.field().magneticField(pos);
  // cout << Acts::UnitConstants::T  << "  " <<  dd4hep::tesla << endl;
  // cout << field.x() << "  " << field.y() << "  " << field.z() << endl;
  // cout << field.x()* (Acts::UnitConstants::T / dd4hep::tesla) << "  " << field.y()* (Acts::UnitConstants::T / dd4hep::tesla) << "  " << field.z()* (Acts::UnitConstants::T / dd4hep::tesla) << endl;
  // Stepper and Propagator
  using Stepper    = Acts::EigenStepper<>;
  using Propagator = Acts::Propagator<Stepper>;

  Stepper stepper(field_provider);
  Propagator propagator(stepper);

  // Create Impact Point Estimator
  Acts::ImpactPointEstimator::Config ImPoEs_cfg(field_provider,std::make_shared<Propagator>(propagator));

  Acts::ImpactPointEstimator::State ImPoEs_state;
  ImPoEs_state.fieldCache = field_cache;

  Acts::ImpactPointEstimator ImPoEs(ImPoEs_cfg);

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
      
      int nAssoc = assocChRecID.GetSize();
      map<int, int> assoc_map;

      for(int j=0; j<nAssoc; j++)
	{
	  assoc_map[assocChSimID[j]] = assocChRecID[j];
	}

      bool hasD0 = false;
      for(int imc=0; imc<nMCPart; imc++)
	{
	  // loop over primary particles
	  if(mcPartGenStatus[imc] == 1 && mcPartCharge[imc] != 0)
	    {
	      double dist = sqrt( pow(mcPartVx[imc]-vertex_mc.x(),2) + pow(mcPartVy[imc]-vertex_mc.y(),2) + pow(mcPartVz[imc]-vertex_mc.z(),2));      
	      if(dist < 1e-4)
		{
		  // check if the MC particle is reconstructed
		  int rc_index = -1;
		  if(assoc_map.find(imc) != assoc_map.end()) rc_index = assoc_map[imc];
		  if(rc_index>-1)
		    {
		      TVector3 rc_vec(rcMomPx[rc_index], rcMomPy[rc_index], rcMomPz[rc_index]);

		      Acts::BoundVector params;
		      params(Acts::eBoundLoc0)   = rcTrkLoca[rc_index];
		      params(Acts::eBoundLoc1)   = rcTrkLocb[rc_index];
		      params(Acts::eBoundPhi)    = rcTrkPhi[rc_index];
		      params(Acts::eBoundTheta)  = rcTrkTheta[rc_index];
		      params(Acts::eBoundQOverP) = rcTrkqOverP[rc_index];
		      params(Acts::eBoundTime)   = 0;

		      //FIXME: Set covariance matrix based on input ROOT file information
		      Acts::BoundSquareMatrix cov = Acts::BoundSquareMatrix::Zero();

		      Acts::Vector3 mc_vtx_pos(vertex_mc.x() * Acts::UnitConstants::mm, vertex_mc.y() * Acts::UnitConstants::mm, vertex_mc.z() * Acts::UnitConstants::mm);

		      Acts::Vector3 rc_vtx_pos(vertex_rc.x() * Acts::UnitConstants::mm, vertex_rc.y() * Acts::UnitConstants::mm, vertex_rc.z() * Acts::UnitConstants::mm);
			  
		      if(fabs(mcPartPdg[imc]) == 211 || fabs(mcPartPdg[imc]) == 321)
			{
			  // Acts::ParticleHypothesis particle_hypothesis;
			  // if(fabs(mcPartPdg[imc]) == 211) particle_hypothesis = Acts::ParticleHypothesis::pion();
			  // if(fabs(mcPartPdg[imc]) == 321) particle_hypothesis = Acts::ParticleHypothesis::kaon();

			  int ip = 0;
			  if(fabs(mcPartPdg[imc]) == 321) ip = 1;

			  Acts::BoundTrackParameters track_parameters(perigee,params,cov,Acts::ParticleHypothesis::pion());
			  if(fabs(mcPartPdg[imc]) == 321) track_parameters = Acts::BoundTrackParameters(perigee,params,cov,Acts::ParticleHypothesis::kaon());

			  //--- Get track parameters at 3D DCA to MC primary vertex ----
			  auto result = ImPoEs.estimate3DImpactParameters(trackingGeoCtx,fieldctx,track_parameters,mc_vtx_pos,ImPoEs_state);
			  if(result.ok())
			    {
			      Acts::BoundTrackParameters trk_boundpar_vtx = result.value();
			      const auto& trk_vtx_params  = trk_boundpar_vtx.parameters();
			      auto trk_vtx_gbl_pos = trk_boundpar_vtx.position(trackingGeoCtx);
			      //cout << "real: " << trk_vtx_params[Acts::eBoundLoc0] << "  " << trk_vtx_params[Acts::eBoundLoc1] << endl;

			      double dca_xy = sqrt( pow(trk_vtx_gbl_pos.x()-mc_vtx_pos.x(),2) + pow(trk_vtx_gbl_pos.y()-mc_vtx_pos.y(),2) );
			      double dca_z = trk_vtx_gbl_pos.z()-mc_vtx_pos.z();
       
			      hRcPrimPartLocaToMCVtx[ip]->Fill(rc_vec.Pt(), rc_vec.Eta(), dca_xy);
			      hRcPrimPartLocbToMCVtx[ip]->Fill(rc_vec.Pt(), rc_vec.Eta(), dca_z);
			      // hRcPrimPartLocaToMCVtx[ip]->Fill(rc_vec.Pt(), rc_vec.Eta(), trk_vtx_params[Acts::eBoundLoc0]);
			      // hRcPrimPartLocbToMCVtx[ip]->Fill(rc_vec.Pt(), rc_vec.Eta(), trk_vtx_params[Acts::eBoundLoc1]);
			    }

			  //--- Get track parameters at 3D DCA to RC primary vertex ----
			  if(prim_vtx_index.GetSize()>0)
			    {
			      auto result = ImPoEs.estimate3DImpactParameters(trackingGeoCtx,fieldctx,track_parameters,rc_vtx_pos,ImPoEs_state);
			      if(result.ok())
				{
				  Acts::BoundTrackParameters trk_boundpar_vtx = result.value();
				  const auto& trk_vtx_params  = trk_boundpar_vtx.parameters();
				  auto trk_vtx_gbl_pos = trk_boundpar_vtx.position(trackingGeoCtx);
				  double dca_xy = sqrt( pow(trk_vtx_gbl_pos.x()-rc_vtx_pos.x(),2) + pow(trk_vtx_gbl_pos.y()-rc_vtx_pos.y(),2) );
				  double dca_z = trk_vtx_gbl_pos.z()-rc_vtx_pos.z();
				  hRcPrimPartLocaToRCVtx[ip]->Fill(rc_vec.Pt(), rc_vec.Eta(), dca_xy);
				  hRcPrimPartLocbToRCVtx[ip]->Fill(rc_vec.Pt(), rc_vec.Eta(), dca_z);
				}
			    }
			}
		    }
		}
            }

	  // D0
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
	      int daug_pi_index = fabs(daug_pdg_1)==211 ? daug_index_1 : daug_index_2;
	      int daug_k_index  = fabs(daug_pdg_1)==321 ? daug_index_1 : daug_index_2;
	      for(int ip = 0; ip<2; ip++)
		{
		  int mc_part_index;
		  if(ip==0) mc_part_index = daug_pi_index;
		  if(ip==1) mc_part_index = daug_k_index;

		  TLorentzVector mc_part_vec;
		  mc_part_vec.SetXYZM(mcMomPx[mc_part_index], mcMomPy[mc_part_index], mcMomPz[mc_part_index], mcPartMass[mc_part_index]);
		  if(ip==0) hMcPiPtEta->Fill(mc_part_vec.Eta(), mc_part_vec.Pt());
		  if(ip==1) hMcKPtEta->Fill(mc_part_vec.Eta(), mc_part_vec.Pt());

		  int rc_part_index = -1;
		  if(assoc_map.find(mc_part_index) != assoc_map.end()) rc_part_index = assoc_map[mc_part_index];

		  if(rc_part_index>=0)
		    {
		      if(ip==0) hMcPiPtEtaReco->Fill(mc_part_vec.Eta(), mc_part_vec.Pt());
		      if(ip==1) hMcKPtEtaReco->Fill(mc_part_vec.Eta(), mc_part_vec.Pt());
		  
		      TVector3 rc_vec(rcMomPx[rc_part_index], rcMomPy[rc_part_index], rcMomPz[rc_part_index]);
		      
		      Acts::BoundVector params;
		      params(Acts::eBoundLoc0)   = rcTrkLoca[rc_part_index];
		      params(Acts::eBoundLoc1)   = rcTrkLocb[rc_part_index];
		      params(Acts::eBoundPhi)    = rcTrkPhi[rc_part_index];
		      params(Acts::eBoundTheta)  = rcTrkTheta[rc_part_index];
		      params(Acts::eBoundQOverP) = rcTrkqOverP[rc_part_index];
		      params(Acts::eBoundTime)   = 0;

		      //FIXME: Set covariance matrix based on input ROOT file information
		      Acts::BoundSquareMatrix cov = Acts::BoundSquareMatrix::Zero();

		      Acts::Vector3 mc_vtx_pos(vertex_mc.x() * Acts::UnitConstants::mm, vertex_mc.y() * Acts::UnitConstants::mm, vertex_mc.z() * Acts::UnitConstants::mm);

		      Acts::Vector3 rc_vtx_pos(vertex_rc.x() * Acts::UnitConstants::mm, vertex_rc.y() * Acts::UnitConstants::mm, vertex_rc.z() * Acts::UnitConstants::mm);

		      
		      Acts::BoundTrackParameters track_parameters(perigee,params,cov,Acts::ParticleHypothesis::pion());
		      if(ip==1) track_parameters = Acts::BoundTrackParameters(perigee,params,cov,Acts::ParticleHypothesis::kaon());

		      //--- Get track parameters at 3D DCA to MC primary vertex ----
		      auto result = ImPoEs.estimate3DImpactParameters(trackingGeoCtx,fieldctx,track_parameters,mc_vtx_pos,ImPoEs_state);
		      if(result.ok())
			{
			  Acts::BoundTrackParameters trk_boundpar_vtx = result.value();
			  const auto& trk_vtx_params  = trk_boundpar_vtx.parameters();
			  auto trk_vtx_gbl_pos = trk_boundpar_vtx.position(trackingGeoCtx);
			  double dca_xy = sqrt( pow(trk_vtx_gbl_pos.x()-mc_vtx_pos.x(),2) + pow(trk_vtx_gbl_pos.y()-mc_vtx_pos.y(),2) );
			  double dca_z = trk_vtx_gbl_pos.z()-mc_vtx_pos.z();
			  hRcSecPartLocaToMCVtx[ip]->Fill(rc_vec.Pt(), rc_vec.Eta(), dca_xy);
			  hRcSecPartLocbToMCVtx[ip]->Fill(rc_vec.Pt(), rc_vec.Eta(), dca_z);
			}

		      
		      //--- Get track parameters at 3D DCA to RC primary vertex ----
		      if(prim_vtx_index.GetSize()>0)
			{
			  auto result = ImPoEs.estimate3DImpactParameters(trackingGeoCtx,fieldctx,track_parameters,rc_vtx_pos,ImPoEs_state);
			  if(result.ok())
			    {
			      Acts::BoundTrackParameters trk_boundpar_vtx = result.value();
			      const auto& trk_vtx_params  = trk_boundpar_vtx.parameters();
			      auto trk_vtx_gbl_pos = trk_boundpar_vtx.position(trackingGeoCtx);
			      double dca_xy = sqrt( pow(trk_vtx_gbl_pos.x()-rc_vtx_pos.x(),2) + pow(trk_vtx_gbl_pos.y()-rc_vtx_pos.y(),2) );
			      double dca_z = trk_vtx_gbl_pos.z()-rc_vtx_pos.z();
			      hRcSecPartLocaToRCVtx[ip]->Fill(rc_vec.Pt(), rc_vec.Eta(), dca_xy);
			      hRcSecPartLocbToRCVtx[ip]->Fill(rc_vec.Pt(), rc_vec.Eta(), dca_z);
			    }
			}
		    }
		}

	      // decay vertex
	      TVector3 mc_vtx_decay(-999.,-999.,-999.);
	      mc_vtx_decay.SetXYZ(mcEndPointX[imc], mcEndPointY[imc], mcEndPointZ[imc]);
	      double mc_decay_dvx = mc_vtx_decay.x()-vertex_mc.x();
	      double mc_decay_dvy = mc_vtx_decay.y()-vertex_mc.y();
	      double mc_decay_dvz = mc_vtx_decay.z()-vertex_mc.z();
	      double mc_decay_dvr = mc_decay_dvx*mc_decay_dvx + mc_decay_dvy*mc_decay_dvy;
	      hD0DecayVxVy->Fill(mc_decay_dvx, mc_decay_dvy);
	      hD0DecayVrVz->Fill(mc_decay_dvz, mc_decay_dvr);

	      //printf("[i] Found D0 decay at (%2.4f, %2.4f, %2.4f)\n", mc_vtx_decay.x(), mc_vtx_decay.y(), mc_vtx_decay.z() );

	      // check if the decay vertex is reconstructed
	      TVector3 rc_vtx_decay(-999, -999, -999);
	      for(unsigned int v=0; v<CTVx.GetSize(); v++)
		{
		  if( vtxAssocPart_end[v]-vtxAssocPart_begin[v] != 2) continue;

		  bool found_d0 = true;
		  for(int irc = vtxAssocPart_begin[v]; irc < vtxAssocPart_end[v]; irc++)
		    {
		      int index = vtxAssocPart_index[irc];
		      int iSimPartID = -1;
		      for(int j=0; j<nAssoc; j++)
			{
			  if(assocChRecID[j]==index)
			    {
			      iSimPartID = assocChSimID[j];
			      break;
			    }
			}

		      if(iSimPartID!=daug_index_1 && iSimPartID!=daug_index_2)
			{
			  found_d0 = false;
			  break;
			}
		    }
		  if(found_d0)
		    {
		      hEventStat->Fill(3.5);
		      hMCD0PtRapReco->Fill(mcRap, mcPt);
		      rc_vtx_decay.SetXYZ(CTVx[v], CTVy[v], CTVz[v]);
		      double rc_decay_dvx = rc_vtx_decay.x() - mc_vtx_decay.x();
		      double rc_decay_dvy = rc_vtx_decay.y() - mc_vtx_decay.y();
		      double rc_decay_dvz = rc_vtx_decay.z() - mc_vtx_decay.z();
		      double rc_decay_dvr = rc_decay_dvx*rc_decay_dvx + rc_decay_dvy*rc_decay_dvy;
		      hD0VtxDist->Fill(rc_decay_dvz, rc_decay_dvr);
		      hD0DecayVxVyReco->Fill(mc_decay_dvx, mc_decay_dvy);
		      hD0DecayVrVzReco->Fill(mc_decay_dvz, mc_decay_dvr);
		      printf("[i] Reco D0 decay at (%2.4f, %2.4f, %2.4f)\n", rc_vtx_decay.x(), rc_vtx_decay.y(), rc_vtx_decay.z() );
		    }
		}

	    }
	}
      hEventStat->Fill(0.5);
      hMcVtxX->Fill(vertex_mc.x());
      hMcVtxY->Fill(vertex_mc.y());
      hMcVtxZ->Fill(vertex_mc.z());

      hNRecoVtx->Fill(CTVx.GetSize());
      if(CTVx.GetSize()==2)
	{
	  for(int v=0; v<CTVx.GetSize(); v++)
	    {
	      // printf("[i] Vertex %d: (x, y, z) = (%2.2f, %2.2f, %2.2f), Ntrk = %d, chi2/ndf = %2.2f, (dx, dy, dz) to MC vtx = (%2.2f, %2.2f, %2.2f)\n",
	      // 	     v+1, CTVx[v], CTVy[v], CTVz[v], vtxAssocPart_end[v]-vtxAssocPart_begin[v], CTVchi2[v]/(CTVndf[v]%100000), CTVx[v]-vertex_mc.x(), CTVy[v]-vertex_mc.y(), CTVz[v]-vertex_mc.z());
	    }
	}

      // loop over reconstructed particles
      const bool select_d0 = true;
      if(!select_d0) hasD0 = true;
      if(hasD0 && prim_vtx_index.GetSize()>0)
	{
	  // find pion and kaon based on true pdg
	  vector<int> pi_index;
	  vector<int> k_index;
	  vector<int> pi_dca_index;
	  vector<int> k_dca_index;
	  pi_index.clear();
	  k_index.clear();
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
	      if(fabs(mcPartPdg[iSimPartID]) == 211 || fabs(mcPartPdg[iSimPartID]) == 321)
		{
		  if(fabs(mcPartPdg[iSimPartID]) == 211) pi_index.push_back(rc_index);
		  if(fabs(mcPartPdg[iSimPartID]) == 321) k_index.push_back(rc_index);
		  
		  Acts::BoundVector params;
		  params(Acts::eBoundLoc0)   = rcTrkLoca[rc_index];
		  params(Acts::eBoundLoc1)   = rcTrkLocb[rc_index];
		  params(Acts::eBoundPhi)    = rcTrkPhi[rc_index];
		  params(Acts::eBoundTheta)  = rcTrkTheta[rc_index];
		  params(Acts::eBoundQOverP) = rcTrkqOverP[rc_index];
		  params(Acts::eBoundTime)   = 0;
		  
		  //FIXME: Set covariance matrix based on input ROOT file information
		  Acts::BoundSquareMatrix cov = Acts::BoundSquareMatrix::Zero();
		  Acts::Vector3 rc_vtx_pos(vertex_rc.x() * Acts::UnitConstants::mm, vertex_rc.y() * Acts::UnitConstants::mm, vertex_rc.z() * Acts::UnitConstants::mm);
		  Acts::BoundTrackParameters track_parameters(perigee,params,cov,Acts::ParticleHypothesis::pion());
		  if(fabs(mcPartPdg[iSimPartID]) == 321) track_parameters = Acts::BoundTrackParameters(perigee,params,cov,Acts::ParticleHypothesis::kaon());

		  auto result = ImPoEs.estimate3DImpactParameters(trackingGeoCtx,fieldctx,track_parameters,rc_vtx_pos,ImPoEs_state);
		  if(result.ok())
		    {
		      Acts::BoundTrackParameters trk_boundpar_vtx = result.value();
		      auto trk_vtx_gbl_pos = trk_boundpar_vtx.position(trackingGeoCtx);
		      double dca_xy = sqrt( pow(trk_vtx_gbl_pos.x()-rc_vtx_pos.x(),2) + pow(trk_vtx_gbl_pos.y()-rc_vtx_pos.y(),2) );
		      double dca_z = trk_vtx_gbl_pos.z()-rc_vtx_pos.z();
		      if(dca_xy>0.04)
			{
			  if(fabs(mcPartPdg[iSimPartID]) == 211) pi_dca_index.push_back(rc_index);
			  if(fabs(mcPartPdg[iSimPartID]) == 321) k_dca_index.push_back(rc_index);
			}
		    }
		}
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
		  if(rcTrkqOverP[pi_index[i]]*rcTrkqOverP[k_index[j]]<0)
		    {
		      TLorentzVector parent = pi_mom_vec + k_mom_vec;
		      h3InvMass[0]->Fill(parent.Pt(), parent.Rapidity(), parent.M());
		    }
		}
	    }

	  for(int i=0; i<pi_dca_index.size(); i++)
	    {
	      TLorentzVector pi_mom_vec;
	      pi_mom_vec.SetXYZM(rcMomPx[pi_dca_index[i]], rcMomPy[pi_dca_index[i]], rcMomPz[pi_dca_index[i]], gPionMass);
	      for(int j=0; j<k_dca_index.size(); j++)
		{
		  TLorentzVector k_mom_vec;
		  k_mom_vec.SetXYZM(rcMomPx[k_dca_index[j]], rcMomPy[k_dca_index[j]], rcMomPz[k_dca_index[j]], gKaonMass);
		  if(rcTrkqOverP[pi_dca_index[i]]*rcTrkqOverP[k_dca_index[j]]<0)
		    {
		      TLorentzVector parent = pi_mom_vec + k_mom_vec;
		      h3InvMass[1]->Fill(parent.Pt(), parent.Rapidity(), parent.M());
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

  hD0DecayVxVyReco->Write();
  hD0DecayVrVzReco->Write();
  
  hMCD0PtRap->Write();
  hMCD0PtRapReco->Write();
  hD0VtxDist->Write();

  hMcPiPtEta->Write();
  hMcPiPtEtaReco->Write();
  hMcKPtEta->Write();
  hMcKPtEtaReco->Write();
  
  hNRecoVtx->Write();

  for(int ip=0; ip<2; ip++)
    {
      hRcSecPartLocaToMCVtx[ip]->Write();
      hRcSecPartLocbToMCVtx[ip]->Write();
      hRcSecPartLocaToRCVtx[ip]->Write();
      hRcSecPartLocbToRCVtx[ip]->Write();

      hRcPrimPartLocaToMCVtx[ip]->Write();
      hRcPrimPartLocbToMCVtx[ip]->Write();
      hRcPrimPartLocaToRCVtx[ip]->Write();
      hRcPrimPartLocbToRCVtx[ip]->Write();
    }

  for(int i=0; i<2; i++)
    {
      h3InvMass[i]->Write();
    }
  
  
  outfile->Close();

}
