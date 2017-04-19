// LaserDataMerger_module.cc
// A basic "skeleton" to read in art::Event records from a file,
// access their information, and do something with them. 

// See
// <https://cdcvs.fnal.gov/redmine/projects/larsoftsvn/wiki/Using_the_Framework>
// for a description of the ART classes used here.

// Almost everything you see in the code below may have to be changed
// by you to suit your task. The example task is to make histograms
// and n-tuples related to dE/dx of particle tracks in the detector.

// As you try to understand why things are done a certain way in this
// example ("What's all this stuff about 'auto const&'?"), it will help
// to read ADDITIONAL_NOTES.txt in the same directory as this file.

// Also note that, despite our efforts, the documentation and the practices in
// this code may fall out of date. In doubt, ask!
// The last revision of this code was done on July 2015 with LArSoft v04_17_00.

#ifndef LaserDataMerger_Module
#define LaserDataMerger_Module

// Always include headers defining everything you use, and only that.
// Starting from LArSoft and ending with standard C++ is a good check
// on LArSoft headers too -- if they can't be loaded by their own, it's a bug!

// LArSoft includes
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/GeometryCore.h"

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

// Framework includes
#include "canvas/Utilities/Exception.h"
#include "art/Framework/Core/EDProducer.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"

// uBooNE includes
#include "lardata/Utilities/AssociationUtil.h"

// ROOT includes. Note: To look up the properties of the ROOT classes,
// use the ROOT web site; e.g.,
// <http://root.cern.ch/root/html532/ClassIndex.html>
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFile.h"

// C++ Includes
#include <map>
#include <vector>
#include <array>
#include <string>
#include <cmath>
#include <typeinfo>
#include <utility>
#include <memory>
#include <iterator>
#include <iostream>
#include <algorithm>

#include <fstream>
#include <boost/tokenizer.hpp>

// Laser Module Classes
#include "LaserObjects/LaserBeam.h"


namespace LaserDataMerger {

    class LaserDataMerger : public art::EDProducer {
    public:

        explicit LaserDataMerger(fhicl::ParameterSet const &parameterSet);

        virtual void beginJob() override;

        virtual void beginRun(art::Run &run);

        virtual void endJob() override;

        virtual void reconfigure(fhicl::ParameterSet const &parameterSet) override;

        // The analysis routine, called once per event.
        virtual void produce(art::Event &event) override;

        double LinearRawToAngle(double Angles);

        float AttenuatorTickToPercentage(float Tick);

    private:

        bool fDebug = false;

        // All this goes into the root tree
        TTree *fTimeAnalysis;
        int fEvent;
        unsigned int time_s;
        unsigned int time_ms;

        std::map<Long64_t, unsigned int> timemap; ///< Key value: index of event, corresponding index in laser data file
        std::vector<std::vector<double> > laser_values; ///< line by line csv container
        double MaxTimeDelta = 0.1;

        bool fReadTimeMap = false;
        bool fGenerateTimeInfo = false;

        float fTickToAngle;     ///< Conversion constant from linear tick to angle (Heidenhain linear encoder)
        std::array<float, 2> fDirCalLCS1 = {{-999., 999.}};  ///< Position calibration for LCS1 and LCS2:
        ///< The first value stands for the offset of the horizontal
        ///< encoder reading with respect to a straight line along the
        ///< z-axis int the x-z plane (top view).
        ///< The second value stands for the offset of the vertical
        ///< encoder reading with respect to a line along the z-axis in the
        ///< y-z plane (side view).
        TVector2 DirCalLCS1;
        std::array<float, 2> fDirCalLCS2 = {{-999., 999.}};  ///< Same as above but for other laser system
        TVector2 DirCalLCS2;
        std::array<float, 3> fPositionLCS1 = {{-999., -999., -999}};   ///< Mirror Position of LCS1
        TVector3 PositionLCS1;
        std::array<float, 3> fPositionLCS2 = {{-999., -999., -999}};   ///< Mirror Position of LCS2
        TVector3 PositionLCS2;

        std::array<float, 2> fEnergyMinMaxLCS1 = {{-999., 999.}};
        std::array<float, 2> fEnergyMinMaxLCS2 = {{-999., 999.}};

        enum DataStructure {
            LaserSystem, ///< which laser system: 1 or 2
            Status, ///< not defined yet
            RotaryPosition, ///< Position Rotary Heidenhain Encoder
            LinearPosition, ///< Position Linear Heidenhain Encoder
            AttenuatorPosition, ///< Position Attenuator Watt Pilot
            AperturePosition, ///< Position Iris Standa
            TriggerTimeSec, ///< Epoch time (in seconds) of Laser Server at receive of Encoder data
            TriggerTimeUsec, ///< Fraction to add to Epoch time in microseconds
            TriggerCount, ///< Trigger Counter by Heidenhain Encoder
            RunControlStep, ///< Run Counter of step in calibration run
            LaserShotCounter, ///< umber of pulses shot with UV laser (not yet read out)
            MirrorBoxAxis1, ///< Motorized Mirror Zaber T-OMG at box
            MirrorBoxAxis2, ///< Motorized Mirror Zaber T-OMG at box
            MirrorFeedthroughAxis1, ///< Motorized Mirror Zaber T-OMG at flange
            MirrorFeedthroughAxis2 ///< Motorized Mirror Zaber T-OMG at flange
        };

        unsigned int LCS_ID = -1;
        unsigned int RunNumber = 0;

    }; // class LaserDataMerger


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// class implementation

//-----------------------------------------------------------------------
// Constructor

    LaserDataMerger::LaserDataMerger(fhicl::ParameterSet const &pset) {

        // Read in the parameters from the .fcl file.
        this->reconfigure(pset);

        produces<lasercal::LaserBeam>("LaserBeam");
    }


//-----------------------------------------------------------------------

    void LaserDataMerger::beginJob() {

    }

    void LaserDataMerger::beginRun(art::Run &run) {

        art::ServiceHandle<art::TFileService> tfs;
        RunNumber = run.run();
/*
        // read the timemap root file (generated in python)
        //std::string TimemapFile = "TimeMap-" + std::to_string(RunNumber) + ".root";
        //std::cout << "READING TIMEMAP FILE: " << TimemapFile << std::endl;
        //TFile *InputFile = new TFile(TimemapFile.c_str(), "READ");
        //TTree *tree = (TTree *) InputFile->Get("tree");

        unsigned int map_root;
        tree->SetBranchAddress("map", &map_root);
        Long64_t nentries = tree->GetEntries();

        for (Long64_t idx = 0; idx < nentries; idx++) {
            tree->GetEntry(idx);
            timemap.insert(std::pair<Long64_t, unsigned int>(idx, map_root));

            //if (fDebug)
            //{
            //    std::cout << "idx: " << idx << " mapped to: " << timemap.at(idx) << std::endl;
            //}
        }
        delete InputFile;

*/
        // read the laser data file into a vector
        std::string LaserFile = "Run-" + std::to_string(RunNumber) + ".txt";
        std::fstream file(LaserFile, std::ios::in);

        if (file) {
            typedef boost::tokenizer<boost::char_separator<char> > Tokenizer;
            boost::char_separator<char> sep(" ");
            std::string line;

            while (getline(file, line)) {
                Tokenizer info(line, sep); // tokenize the line of data
                std::vector<double> values;

                for (Tokenizer::iterator it = info.begin(); it != info.end(); ++it) {
                    // convert data into double value, and store
                    values.push_back(std::strtod(it->c_str(), 0));
                }

                // store array of values
                laser_values.push_back(values);

                if (fDebug) {
                    for (unsigned int idx = 0; idx < 15; idx++) {
                        //std::cout << std::fixed<<laser_values.back().at(idx) << " ";
                    }
                    //std::cout << std::endl;
                }

            }
            LCS_ID = laser_values.back().at(DataStructure::LaserSystem);
            file.close();
        } else {
            std::cerr << "Error: Unable to open file " << LaserFile << std::endl;
        }
        return;
    }

    void LaserDataMerger::endJob() {
    }

    void LaserDataMerger::reconfigure(fhicl::ParameterSet const &parameterSet) {
        // Read parameters from the .fcl file. The names in the arguments
        // to p.get<TYPE> must match names in the .fcl file.
        fGenerateTimeInfo = parameterSet.get<bool>("GenerateTimeInfo");
        fTickToAngle = parameterSet.get<float>("TickToAngle");
        fDirCalLCS1 = parameterSet.get<std::array<float, 2> >("DirCalLCS1");
        fDirCalLCS2 = parameterSet.get<std::array<float, 2> >("DirCalLCS2");
        fPositionLCS1 = parameterSet.get<std::array<float, 3> >("PositionLCS1");
        fPositionLCS2 = parameterSet.get<std::array<float, 3> >("PositionLCS2");

        fEnergyMinMaxLCS1 = parameterSet.get<std::array<float, 2> >("EnergyMinMaxLCS1");
        fEnergyMinMaxLCS2 = parameterSet.get<std::array<float, 2> >("EnergyMinMaxLCS2");

        // Convert calibration angles to TVector2
        DirCalLCS1.Set(fDirCalLCS1[0],
                       fDirCalLCS1[1]);     ///< LCS1 Calibration values for Rotary and Linear Values (in this order)
        DirCalLCS2.Set(fDirCalLCS2[0],
                       fDirCalLCS2[1]);     ///< LCS2 Calibration values for Rotary and Linear Values (in this order)

        // Convert the Positions to a TVector3
        PositionLCS1.SetXYZ(fPositionLCS1[0], fPositionLCS1[1], fPositionLCS1[2]);
        PositionLCS2.SetXYZ(fPositionLCS2[0], fPositionLCS2[1], fPositionLCS2[2]);

        fDebug = parameterSet.get<bool>("Debug", false);

        //fLaserSystemFile        = parameterSet.get< bool        >("LaserSystemFile");
    }

//-----------------------------------------------------------------------

    void LaserDataMerger::produce(art::Event &event) {
        fEvent = (int) event.id().event();
        int laser_id = -9999;
        std::unique_ptr<lasercal::LaserBeam> LaserAA(new lasercal::LaserBeam());
      
        if(laser_values[0][8]<0){laser_id=0;}//If the laser direction does not change, set event number to -1 (manually).
        else{       
          unsigned int EvtTimeS = event.time().timeHigh();
          unsigned int EvtTimeNS = event.time().timeLow();
          double daq_event_time = EvtTimeS + EvtTimeNS * 1E-9;//ns
          double laser_event_time = -9999;

    
          int search_start = std::max(-fEvent,-20);
          int tempDiff = laser_values.back()[8]-fEvent;
          int search_end = std::min(tempDiff,200);
          double max_laser_time = laser_values.back()[6] + laser_values.back()[7] * 1E-6;
          if(daq_event_time>max_laser_time){
	    std::cout<<"Laser Information is missing. Stoped"<<std::endl;
	    return;
	  }

          for(int j=search_start;j<search_end;j++){
            laser_event_time = laser_values[fEvent+j][6] + laser_values[fEvent+j][7] * 1E-6;//us
            //match the laser time and daq time
            if(TMath::Abs(laser_event_time - daq_event_time) < MaxTimeDelta){
              std::cout<<std::fixed<<"laser: "<<laser_event_time<<"; daq: "<<daq_event_time<<"; Diff: "<<(double)TMath::Abs(laser_event_time - daq_event_time)<<std::endl;
              laser_id = fEvent+j;
              break;
            }
            if(j==search_end-1){
              std::cout<<"Could not find corresponding laser event time for daq event "<<fEvent<<". Something went wrong!"<<std::endl;

              lasercal::LaserBeam Laser(TVector3(-1000, -1000, -1000), -1000, -1000);
              Laser.SetLaserID(9999);

              *LaserAA = Laser;
              event.put(std::move(LaserAA), "LaserBeam");
              return;
            }
          }
        }         
	double Theta;
        double Phi;

        double Theta_raw = laser_values.at(laser_id).at(DataStructure::LinearPosition);
        double Phi_raw = laser_values.at(laser_id).at(DataStructure::RotaryPosition);

        TVector3 Position;
        TVector2 CalibratedAngles;

        if (LCS_ID == 1) { // The downstream laser system (sitting at z = -20)
            Theta = TMath::DegToRad() * (90.0 - LinearRawToAngle(Theta_raw - fDirCalLCS1[1]));
            Phi = TMath::DegToRad() * (Phi_raw - fDirCalLCS1[0]);
            Position = PositionLCS1;
        } else if (LCS_ID == 2) { // The upstream laser system (sitting at z = 1020)
            Theta = TMath::DegToRad() * (LinearRawToAngle(Theta_raw) - 266.60208631 + 60 + 1.161  - 2*0.68) ;
            Phi = -TMath::DegToRad() * (180 + (Phi_raw - fDirCalLCS2[0]));
            Position = PositionLCS2;
        } else {
            std::cerr << "Laser System not recognized " << std::endl;
        }

        if (fDebug) {
            time_s = (unsigned int) event.time().timeHigh();
            time_ms = (unsigned int) event.time().timeLow();

            std::cout << "Positions from entry: " << laser_id << std::endl;
            std::cout << "rot: (raw / calib): " << Theta_raw << " / " << Theta << std::endl;
            std::cout << "lin: (raw / calib): " << Phi_raw << " / " << Phi << std::endl;
            std::cout << "Event Time (low): " << time_s << std::endl;
            std::cout << "Event Time (hig): " << time_ms << std::endl;
        }


        //CalibratedAngles.Set(TMath::DegToRad() * 45, TMath::DegToRad() * 190);
        lasercal::LaserBeam Laser(Position, Phi, Theta);
        Laser.SetLaserID(LCS_ID);
        Laser.SetLaserEventID(laser_values.at(laser_id).at(DataStructure::TriggerCount));
        Laser.SetAssID(laser_id);
        Laser.SetPower(AttenuatorTickToPercentage(laser_values.at(laser_id).at(DataStructure::AttenuatorPosition)));
        Laser.SetTime(laser_values.at(laser_id).at(DataStructure::TriggerTimeSec),
                      laser_values.at(laser_id).at(DataStructure::TriggerTimeUsec));
        if (fDebug) Laser.Print();


        *LaserAA = Laser;

        event.put(std::move(LaserAA), "LaserBeam");

        //for (auto i : fPositionLCS1)
        //    std::cout << i << std::endl;

    }

    float LaserDataMerger::AttenuatorTickToPercentage(float Tick) {
        if (LCS_ID == 1) {
            return (Tick - fEnergyMinMaxLCS1[0]) / (fEnergyMinMaxLCS1[1] - fEnergyMinMaxLCS1[0]);
        } else if (LCS_ID == 2) {
            return (Tick - fEnergyMinMaxLCS2[0]) / (fEnergyMinMaxLCS2[1] - fEnergyMinMaxLCS1[0]);
        } else return -999.;
    }

    double LaserDataMerger::LinearRawToAngle(double RawTicks) {

        float linear2angle = 0.3499;         ///< conversion constant of linear encoder to angle (mm/deg)
        float TickLength = 0.00001;          ///< Tick length in mm
        //float err_linear2angle = 0.0002;     ///< error of conversion factor
        return RawTicks * TickLength / linear2angle;
    }

    DEFINE_ART_MODULE(LaserDataMerger)

} // namespace LaserDataMerger

#endif // LaserDataMerger_Module
