//================================================
// Study the kinematics of produced hadron in 
// certain (x, Q2, nu) bin
//
// To-do:
// 1. Exclude/include feed-down?
//================================================
R__LOAD_LIBRARY(libeicsmear);

#include "bins.h"

using namespace std;


void HF_filter(const char* inFile = "eAu.root", const int select_type = 0)
{

  // Prepare output name
  TString outname = inFile;
  if(select_type==0) outname = outname.ReplaceAll(".root", ".filterD0.root");
  else outname = outname.ReplaceAll(".root", ".filtered.root");
  
  TFile *f = new TFile(inFile);

  //Get EICTree Tree
  TTree *tree = (TTree*)f->Get("EICTree");

  Int_t nEntries = tree->GetEntries();
  cout<<"-------------------------------"<<endl;
  cout<<"Total Number of Events = "<<nEntries<<endl<<endl;

  // Event Class
  erhic::EventBeagle *event_beagle(NULL);

  // Access event Branch
  tree->SetBranchAddress("event",&event_beagle);

  // new tree
  TFile* fout = new TFile(outname.Data(),"recreate");
  fout->SetCompressionLevel(9);
  TTree *newtree = tree->CloneTree(0);
  
  //Loop Over Events
  int nevt = nEntries;
  for(Int_t ievt = 0; ievt < nevt; ievt++)
  {    
    if (ievt%100==0) cout<<"Processing event = "<<ievt<<"/"<<nevt<<endl;

    tree->GetEntry(ievt);

    erhic::EventMC* event = event_beagle;

    //cout << endl << "[i] Event " << ievt << endl;

    bool is_selected = true;
    if(select_type>=0) is_selected = false;

    for(int ipart = 0; ipart < event->GetNTracks(); ipart++)
      {
	erhic::ParticleMC* part = event->GetTrack(ipart);

	if(select_type==0)
	  {
	    // select events with D0 -> pi+k
	    if (abs(part->Id())==421)
	      {
		cout << "[i] Found a D0" << endl;
		int nDaug = part->GetNChildren();
		if(nDaug==2)
		  {
		    // int start = part->GetChild1Index();
		    // int end = part->GetChildNIndex();
		    // cout << part->GetNChildren() << " -> " << start << "  " << end << endl;
		    // for(int i=start; i<=end; i++)
		    //   {
		    //     //const erhic::ParticleMC* daughter = part->GetChild(i);
		    //     erhic::ParticleMC* daughter = event->GetTrack(i-1);
		    //     const erhic::ParticleMC* daughter1 = part->GetChild(i-start);
		    //     cout << "[i] Daughter " << i << ", id = " << daughter->Id() << " =? " << daughter1->Id() << endl;
		    //   }
		    const erhic::ParticleMC* daughter1 = part->GetChild(0);
		    const erhic::ParticleMC* daughter2 = part->GetChild(1);
		    cout << daughter1->Id() << "  " << daughter2->Id() << endl;
		    
		    if( (abs(daughter1->Id())==211 && abs(daughter2->Id())==321) ||
			(abs(daughter1->Id())==321 && abs(daughter2->Id())==211))
		      {
			is_selected = true;
			break;
		      }
		  }
	      }
	  }
      }

    if(is_selected)
      {
	newtree->Fill();
      }
  }
  
  printf("[i] Acceptance fraction %d/%d = %4.2f%%\n",newtree->GetEntries(),nEntries,newtree->GetEntries()*1./nEntries*100);
 
  newtree->AutoSave();
  fout->Close();
}
