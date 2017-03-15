

#ifndef GetTracksRootTest_Module
#define GetTracksRootTest_Module

// Always include headers defining everything you use, and only that.
// Starting from LArSoft and ending with standard C++ is a good check
// on LArSoft headers too -- if they can't be loaded by their own, it's a bug!

// LArSoft includes
#include "lardataobj/RecoBase/Track.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"

// uBooNE includes
#include "lardata/Utilities/AssociationUtil.h"

// ROOT includes.
#include "TTree.h"
#include "TFile.h"

#include <fstream>

// Laser Module Classes
#include "LaserObjects/LaserBeam.h"


namespace GetTracksRootTest
{

class GetTracksRootTest : public art::EDProducer
{
public:

    explicit GetTracksRootTest(fhicl::ParameterSet const& parameterSet);

    virtual void beginJob() override;

    virtual void beginRun(art::Run& run);

    virtual void endJob() override;

    virtual void reconfigure(fhicl::ParameterSet const& parameterSet) override;

    // The analysis routine, called once per event. 
    virtual void produce(art::Event& event) override;
   

private:

    bool DEBUG = true;

    // All this goes into the root tree
    TTree* fTrackTree;
    TTree* fLaserTree;

//    std::vector<Double_t> trackx;
//    std::vector<Double_t> tracky;
//    std::vector<Double_t> trackz;
    Double_t trackx;
    Double_t tracky;
    Double_t trackz;
    unsigned int event_id;
    unsigned int track_id;
    unsigned int track_N = 0;

    Double_t laser_entry_x;
    Double_t laser_entry_y;
    Double_t laser_entry_z;

    Double_t laser_exit_x;
    Double_t laser_exit_y;
    Double_t laser_exit_z;

    art::InputTag fTrackLabel;

}; // class GetTracks



GetTracksRootTest::GetTracksRootTest(fhicl::ParameterSet const& pset)
{

    // Read in the parameters from the .fcl file.
    this->reconfigure(pset);
}


//-----------------------------------------------------------------------

void GetTracksRootTest::beginJob()
{

}

void GetTracksRootTest::beginRun(art::Run& run)
{
    
    art::ServiceHandle<art::TFileService> tfs;

    // Initialize time info root file
    fTrackTree = tfs->make<TTree>("Tracks", "Tracks");
    fTrackTree->Branch("event", &event_id);
    fTrackTree->Branch("x", &trackx);
    fTrackTree->Branch("y", &tracky);
    fTrackTree->Branch("z", &trackz);
    fTrackTree->Branch("track_id", &track_id);
    fTrackTree->Branch("track_N", &track_N);

    fLaserTree = tfs->make<TTree>("Laser", "Laser");
    fLaserTree->Branch("event", &event_id);
    fLaserTree->Branch("entry_x", &laser_entry_x);
    fLaserTree->Branch("entry_y", &laser_entry_y);
    fLaserTree->Branch("entry_z", &laser_entry_z);
    fLaserTree->Branch("exit_x", &laser_exit_x);
    fLaserTree->Branch("exit_y", &laser_exit_y);
    fLaserTree->Branch("exit_z", &laser_exit_z);

    return;
}

void GetTracksRootTest::endJob()
{
}

void GetTracksRootTest::reconfigure(fhicl::ParameterSet const& parameterSet)
{
    fTrackLabel = parameterSet.get<art::InputTag>("TrackLabel");
}

//-----------------------------------------------------------------------

void GetTracksRootTest::produce(art::Event& event)
{
    art::Handle<std::vector<recob::Track> > Tracks;
    art::Handle<lasercal::LaserBeam>  Laser;

    event.getByLabel(fTrackLabel, Tracks);

    try {
        event.getByLabel("LaserBeam", Laser);

    	event_id = event.id().event();
        laser_entry_x = Laser->GetEntryPoint().x();
        laser_entry_y = Laser->GetEntryPoint().y();
        laser_entry_z = Laser->GetEntryPoint().z();

        laser_exit_x = Laser->GetExitPoint().x();
        laser_exit_y = Laser->GetExitPoint().y();
        laser_exit_z = Laser->GetExitPoint().z();

        fLaserTree->Fill();
    }
    catch (...){}; // pretty dangerous, but we just ignore writing the laser tree if no laser data is present.


    //auto track = tr.fXYZ;
    for (auto const &Track : *Tracks) {

//        event_id = (unsigned int) event.id().event();
//        track_id = Track.ID();

        size_t track_size = Track.NumberTrajectoryPoints();

        //trackx.resize(track_size); tracky.resize(track_size); trackz.resize(track_size);

        for (size_t i = 0; i < track_size; ++i)  {
            event_id = (unsigned int) event.id().event();
            track_id = Track.ID();
            trackx = Track.LocationAtPoint(i).x();
            tracky = Track.LocationAtPoint(i).y();
            trackz = Track.LocationAtPoint(i).z();
            fTrackTree->Fill();
        }
        track_N++;
        //fTrackTree->Fill();
        //trackx.clear(); tracky.clear(); trackz.clear();
    }
 }
    DEFINE_ART_MODULE(GetTracksRootTest)
} // namespace GetTracks

#endif // GetTracks_Module
