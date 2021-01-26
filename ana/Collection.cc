#include <iostream>
#include <fstream>
#include <stdlib.h>

// check directory
#include <sys/stat.h>
#include <TROOT.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include "TMath.h"
#include "TNtupleD.h"
#include "TLorentzVector.h"

#include "time.h"
#include <dirent.h>
#include <string>
#include <vector>
#include <map>
struct stat sb;

#include "lcio.h"
#include <stdio.h>
#include "IO/LCReader.h"
#include "IMPL/LCTOOLS.h"
#include "EVENT/LCEvent.h"
#include "EVENT/LCRunHeader.h"
#include "EVENT/CalorimeterHit.h"
#include "EVENT/RawCalorimeterHit.h"
#include "EVENT/ReconstructedParticle.h"
#include "IOIMPL/LCFactory.h"
#include "EVENT/MCParticle.h"
#include "EVENT/LCCollection.h"
#include "IMPL/LCEventImpl.h"
#include "UTIL/LCTOOLS.h"
#include "UTIL/Operators.h"
#include "UTIL/LCIterator.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/MCParticleImpl.h"
#include "EVENT/MCParticle.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/Track.h"
#include "EVENT/Cluster.h"
#include "EVENT/CalorimeterHit.h"
#include "EVENT/SimCalorimeterHit.h"
#include "Objects/CartesianVector.h"
#include "Objects/Helix.h"
#include "DetectorGeometrySimple.h"
#include "IDDecoder.h"
#include "Pandora/Pandora.h"

#include "LParticle.h"
#include "CParticle.h"
using namespace std;

const double kPI = TMath::Pi();
const double k2PI = 2 * kPI;
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
using namespace fastjet;

// find all files inside a directory
std::vector<std::string> open(std::string name = "data.in")
{
    vector<std::string> ntup;
    ifstream myfile;
    myfile.open(name.c_str(), ios::in);

    if (!myfile)
    {
        cerr << " -> Can't open input file:  " << name << endl;
        exit(1);
    }
    else
    {
        cout << "-> Read data file=" << name << endl;
    }

    string temp;
    while (myfile >> temp)
    {
        //the following line trims white space from the beginning of the string
        temp.erase(temp.begin(), std::find_if(temp.begin(), temp.end(), not1(ptr_fun<int, int>(isspace))));
        if (temp.find("#") == 0)
            continue;
        ntup.push_back(temp);
    }
    cout << "-> Number of files=" << ntup.size() << endl;
    myfile.close();

    for (unsigned int i = 0; i < ntup.size(); i++)
    {
        cout << ".. file to analyse=" + ntup[i] << endl;
    }
    return ntup;
}
int main(int argc, char **argv)
{
    std::vector<std::string> files = open("data.in");
    int mtype = 2; // single-particles

    string outputfile = "root/output.root";
    cout << "\n -> Output file is =" << outputfile << endl;
    TFile *RootFile = new TFile(outputfile.c_str(), "RECREATE", "Histogram file");
    TH1D *h_debug = new TH1D("debug", "events", 10, 0, 10);

    // read detector geometry for this configuration
    string detector = "./data/rfull009_sifcch7/sifcch7/sifcch7.pandora";
    DetectorGeometrySimple *geom = new DetectorGeometrySimple(detector);
    string caloType = "HAD_BARREL";
    DetectorGeometrySimple::ExtraSubDetectorParameters *xsubdet = geom->getExtraSubDetectorParametersFromType(caloType);
    if (xsubdet == NULL)
    {
        std::cout << "The ExtraSubDetectorParameters for " << caloType << " were not found." << std::endl;
        throw new std::exception;
    }

    string caloType_ecal = "EM_BARREL";
    DetectorGeometrySimple::ExtraSubDetectorParameters *xsubdet_ecal = geom->getExtraSubDetectorParametersFromType(caloType_ecal);
    if (xsubdet_ecal == NULL)
    {
        std::cout << "The ExtraSubDetectorParameters for " << caloType_ecal << " were not found." << std::endl;
        throw new std::exception;
    }
    // Get the decoder.
    IDDecoder *decoder = xsubdet->m_decoder;
    IDDecoder *decoder_ecal = xsubdet_ecal->m_decoder;

    // calculate position
    double layer_size_mm = 27.5 + 5 + 2.5;
    double hcal_inner_R_mm = 2300;

    double layer_size_ecal_mm = 4.0;
    double ecal_inner_R_mm = 2100;

    // calculate weighted average of time
    double XTime = 0;
    double XEnergy = 0;

    int MaxEvents = 100000000;

    int nEvents = 0;
    for (unsigned int mfile = 0; mfile < files.size(); mfile++)
    {
        string Rfile = files[mfile];

        cout << " # File=" << Rfile << endl;

        IO::LCReader *lcReader = IOIMPL::LCFactory::getInstance()->createLCReader();
        lcReader->open(Rfile.c_str());

        EVENT::LCEvent *evt = 0;

        if (nEvents > MaxEvents)
            break;
        while ((evt = lcReader->readNextEvent()) != 0)
        {
            if (nEvents == 0)
                UTIL::LCTOOLS::dumpEvent(evt);

            // UTIL::LCTOOLS::dumpEvent( evt ) ;

            nEvents++;
            if ((nEvents < 100 && nEvents % 10 == 0) || (nEvents > 100 && nEvents % 200 == 0))
                cout << " # Events=" << nEvents << endl;

            if (nEvents > MaxEvents)
                break;
            //std::string mcpName("MCParticle");
            std::string mcpName("HcalBarrelHits");
            // get truth
            IMPL::LCCollectionVec *col = (IMPL::LCCollectionVec *)evt->getCollection("MCParticle");
            //IMPL::LCCollectionVec *col50 = (IMPL::LCCollectionVec *)evt->getCollection("EcalEndcapHits");
            IMPL::LCCollectionVec *col53 = (IMPL::LCCollectionVec *)evt->getCollection("EcalBarrelHits");
            //IMPL::LCCollectionVec *col51 = (IMPL::LCCollectionVec *)evt->getCollection("HAD_BARREL");
            //IMPL::LCCollectionVec *col52 = (IMPL::LCCollectionVec *)evt->getCollection("EM_BARREL");
        }
        lcReader->close();
        delete lcReader;
    }
    RootFile->Write();
    RootFile->Close();

    return 0;
}