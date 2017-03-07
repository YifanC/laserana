////////////////////////////////////////////////////////////////////////
// Class:       HitEvent
// Plugin Type: analyzer (art v2_05_00)
// File:        HitEvent_module.cc
//
// Generated at Wed Mar  1 12:04:00 2017 by yifanch using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#ifndef HitEventTest_Module
#define HitEventTest_Module

//LArSoft includes
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"

//Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"

//Root includes
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TGraph.h"

//C++ includes
#include <memory>
#include <vector>

class HitEvent;


class HitEvent : public art::EDAnalyzer {
public:
  explicit HitEvent(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  HitEvent(HitEvent const &) = delete;
  HitEvent(HitEvent &&) = delete;
  HitEvent & operator = (HitEvent const &) = delete;
  HitEvent & operator = (HitEvent &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void beginRun(art::Run const & r) override;
  void beginSubRun(art::SubRun const & sr) override;
  void endJob() override;
  void endRun(art::Run const & r) override;
  void endSubRun(art::SubRun const & sr) override;
  void reconfigure(fhicl::ParameterSet const & p) override;
  void respondToCloseInputFile(art::FileBlock const & fb) override;
  void respondToCloseOutputFiles(art::FileBlock const & fb) override;
  void respondToOpenInputFile(art::FileBlock const & fb) override;
  void respondToOpenOutputFiles(art::FileBlock const & fb) override;

private:

  // Declare member data here.
  TTree* tree;
  TBranch* PA;
  TBranch* NH;
  short int PeakAmp;
  int NHits;
  TCanvas* can;
  TGraph* gr;
};


HitEvent::HitEvent(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  
  {
  // More initializers here.
  }

void HitEvent::analyze(art::Event const & e)
{
  // Implementation of required member function here.
  art::ValidHandle<std::vector<raw::RawDigit>>DigitVecHandle = e.getValidHandle<std::vector<raw::RawDigit>>("daq");
  int i=0;
  NHits=0;
  for (auto const &RawDigit : *DigitVecHandle){
    raw::ChannelID_t channel = RawDigit.Channel();
    std::vector<short>ADCs = RawDigit.ADCs();
//    std::cout<<"ADCs size"<<RawDigit.ADCs().size()<<std::endl;
    PeakAmp = *std::max_element(ADCs.begin(),ADCs.end());//super naive.....
    std::cout<<"PeakAmp"<<PeakAmp<<std::endl;
    std::cout<<"channel"<<channel<<std::endl;
    if(PeakAmp>1000){
      NHits++;
    }
if(i==1){
    std::cout<<"pass"<<std::endl;
      Int_t ArrADCs[ADCs.size()];
      copy(ADCs.begin(),ADCs.end(),ArrADCs);
      Int_t SampleN[ADCs.size()];
      for(Int_t ii=0;ii<ADCs.size();ii++){
        SampleN[ii]=ii;
      }
      Int_t n = ADCs.size();

      can->cd();
      gr = new TGraph(n,SampleN,ArrADCs);
      gr->Draw("AC*");
      can->Modified();
      can->Update();
}
    tree->Fill();
    //PA->Fill();
    i++;
  }
 // NH->Fill();
  std::cout<<"i:"<<i<<std::endl;
}

void HitEvent::beginJob()
{
  // Implementation of optional member function here.
can = new TCanvas("a","a",600,800);
}

void HitEvent::beginRun(art::Run const & r)
{
  // Implementation of optional member function here.

  art::ServiceHandle<art::TFileService> tfs;

  tree = tfs->make<TTree>("tree","tree1");
  PA = tree->Branch("PeakAmp",&PeakAmp);
  NH = tree->Branch("NHits",&NHits);
  return;

}

void HitEvent::beginSubRun(art::SubRun const & sr)
{
  // Implementation of optional member function here.
}

void HitEvent::endJob()
{
  // Implementation of optional member function here.
  //OutputFile->Close();
}

void HitEvent::endRun(art::Run const & r)
{
  // Implementation of optional member function here.
}

void HitEvent::endSubRun(art::SubRun const & sr)
{
  // Implementation of optional member function here.
}

void HitEvent::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
}

void HitEvent::respondToCloseInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void HitEvent::respondToCloseOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void HitEvent::respondToOpenInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void HitEvent::respondToOpenOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(HitEvent)
#endif
