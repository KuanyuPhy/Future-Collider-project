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
#include "TNtuple.h"
#include "TTree.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include <fastjet/AreaDefinition.hh>
#include <fastjet/ClusterSequenceArea.hh>
struct stat sb;
const double kPI = TMath::Pi();
const double k2PI = 2 * kPI;
using namespace std;
using namespace fastjet;

float EffectiveRadius(PseudoJet Jet, vector<PseudoJet> constituents, double jetR = 0.5)
{
    float Energy = Jet.Et();
    float numerator = 0;
    int size = constituents.size();
    for (int i = 0; i < size; i++)
    {
        if (Jet.delta_R(constituents[i]) > jetR)
            continue;
        numerator += constituents[i].Et() * Jet.delta_R(constituents[i]);
    }
    //cout << Energy << endl;
    //cout << numerator << endl;
    return numerator / Energy;
}

float eccentricity(PseudoJet Jet, vector<PseudoJet> constituents)
{
    unsigned int num = constituents.size();
    double Dphi[num], Deta[num], E[num];
    double etaSum = 0.;
    double phiSum = 0.;
    double eTot = 0.;

    for (unsigned int j = 0; j < num; j++)
    {
        PseudoJet cp = constituents.at(j);
        E[j] = cp.e();
        Dphi[j] = Jet.phi() - cp.phi();
        // if (Dphi[j]>TMath::Pi()) Dphi[j]=2*TMath::Pi()-Dphi[j];
        if (fabs(Dphi[j] - 2. * TMath::Pi()) < fabs(Dphi[j]))
            Dphi[j] -= 2. * TMath::Pi();
        if (fabs(Dphi[j] + 2. * TMath::Pi()) < fabs(Dphi[j]))
            Dphi[j] += 2. * TMath::Pi();
        Deta[j] = Jet.eta() - cp.eta();
        etaSum = etaSum + Deta[j] * E[j];
        phiSum = phiSum + Dphi[j] * E[j];
        eTot = eTot + E[j];
    }
    etaSum = etaSum / eTot;
    phiSum = phiSum / eTot;
    for (unsigned int j = 0; j < num; j++)
    {
        Deta[j] = Deta[j] - etaSum;
        Dphi[j] = Dphi[j] - phiSum;
    }

    double X1 = 0.;
    double X2 = 0;
    for (unsigned int i = 0; i < num; i++)
    {
        X1 += 2. * E[i] * Deta[i] * Dphi[i];
        X2 += E[i] * (Dphi[i] * Dphi[i] - Deta[i] * Deta[i]);
    }

    // variance calculations
    double Theta = .5 * atan(X1 / X2);
    double sinTheta = TMath::Sin(Theta);
    double cosTheta = TMath::Cos(Theta);
    double Theta2 = Theta + 0.5 * TMath::Pi();
    double sinThetaPrime = TMath::Sin(Theta2);
    double cosThetaPrime = TMath::Cos(Theta2);

    double VarX = 0.;
    double VarY = 0.;
    for (unsigned int i = 0; i < num; i++)
    {
        double X = cosTheta * Deta[i] - sinTheta * Dphi[i];
        double Y = sinTheta * Deta[i] + cosTheta * Dphi[i];
        VarX += E[i] * X * X;
        VarY += E[i] * Y * Y;
    }

    double VarianceMax = VarX;
    double VarianceMin = VarY;
    if (VarianceMax < VarianceMin)
    {
        VarianceMax = VarY;
        VarianceMin = VarX;
    }

    double ECC = 1.0 - (VarianceMin / VarianceMax);

    return ECC;
}
double nsubjettiness(PseudoJet Jet, vector<PseudoJet> constituents, int NSubJets, double jetRad = 0.5)
{
    //vector<CParticle> constit = jet.GetConstituents();
    vector<PseudoJet> jetConstit;

    unsigned int num = constituents.size();

    if (num < (unsigned int)NSubJets)
    {
        return -999;
    }
    TLorentzVector Jet_p;
    Jet_p.SetPxPyPzE(Jet.px(), Jet.py(), Jet.pz(), Jet.E());
    num = 0;
    for (unsigned int i = 0; i < constituents.size(); i++)
    {
        PseudoJet c_i = constituents[i];
        if (c_i.delta_R(Jet_p) > jetRad)
            continue;
        jetConstit.push_back(c_i);
        num++;
    }
    if (num < (unsigned int)NSubJets)
    {
        return -999;
    }
    std::vector<std::vector<fastjet::PseudoJet>> kt_subjets_vec; //a vector of vectors of Pseudojets
    fastjet::JetDefinition *m_jetdef = new fastjet::JetDefinition(fastjet::kt_algorithm, 1.5, fastjet::E_scheme, fastjet::Best);
    fastjet::ClusterSequence kt_seq(jetConstit, *m_jetdef);
    delete m_jetdef;

    for (unsigned int i = 0; i < (unsigned int)NSubJets; i++)
    {
        kt_subjets_vec.push_back(fastjet::sorted_by_pt(kt_seq.exclusive_jets((int)NSubJets)));
    }
    double min_dist = 100000.0;
    double sum_pt = 0.0;
    double sum_dist = 0.0;
    //first find the minimum distance.
    for (unsigned int i = 0; i < jetConstit.size(); i++)
    {
        fastjet::PseudoJet theconstit(jetConstit[i]);
        sum_pt += theconstit.perp();
        float min_dist = 1e10;
        for (unsigned int j = 0; j < (unsigned int)NSubJets; j++)
        {
            const std::vector<fastjet::PseudoJet> *kt_subjets = &(kt_subjets_vec[j]);
            float temp_dist = theconstit.perp() * std::sqrt(kt_subjets->at(j).plain\
_distance(theconstit));
            if (temp_dist < min_dist)
                min_dist = temp_dist;
        } //loop over axis (subjets)
        sum_dist += min_dist;
    } //loop over jet constituents

    sum_dist /= (jetRad * sum_pt);

    double nSubJettiness = sum_dist;
    if (sum_dist > 1.0)
        cout << "uh oh" << sum_dist << endl;

    return nSubJettiness;
}

double splittingscale(PseudoJet Jet)
{
    if (!Jet.has_constituents())
        return -5.;

    vector<PseudoJet> jetConstit = Jet.constituents();

    double dR = Jet.associated_cluster_sequence()->jet_def().R();
    fastjet::JetDefinition *m_jetdef = new fastjet::JetDefinition(fastjet::kt_algorithm, dR, fastjet::E_scheme, fastjet::Best);

    fastjet::ClusterSequence kt_seq(jetConstit, *m_jetdef);
    delete m_jetdef;
    return dR * TMath::Sqrt(kt_seq.exclusive_dmerge(1));
}

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
    //debug
    bool debug = true;
    std::vector<std::string> files = open("data.in");
    int mtype = 1;

    string outputfile = "root/output.root";
    TFile *RootFile = new TFile(outputfile.c_str(), "RECREATE", "Histogram file");
    TH1D *h_debug = new TH1D("debug", "events", 10, 0, 10);
    TH1F *Total_particle_ID_no_cut = new TH1F("Total_particle_ID_no_cut", "Total_particle_ID_no_cut", 20, 0, 20);
    TH1D *Wbosn_nn = new TH1D("Wbosn_nn", "Nr of Wboson", 10, 0, 10);

    TH1F *Timing_Standard = new TH1F("Timing_Standard", "Timing_Standard", 200, 0, 50);
    TH1F *h_Particles_dR_Highest_PT_T_Truth[5];
    TH1F *h_Particles_dR_Highest_PT_PT_Truth[5];
    TH1F *h_Particles_dR_Highest_PT_V_Truth[5];
    TH1F *h_Particles_dR_Highest_PT_Vz_Truth[5];
    for (int j = 0; j < 5; j++)
    {
        h_Particles_dR_Highest_PT_T_Truth[j] = new TH1F(Form("h_Particles_dR_Highest_PT_T_Truth_%i", j), Form("h_Particles_dR_Highest_PT_T_Truth_%i", j), 100, 0, 1);
        h_Particles_dR_Highest_PT_PT_Truth[j] = new TH1F(Form("h_Particles_dR_Highest_PT_PT_Truth_%i", j), Form("h_Particles_dR_Highest_PT_PT_Truth_%i", j), 100, 0, 1);
        h_Particles_dR_Highest_PT_V_Truth[j] = new TH1F(Form("h_Particles_dR_Highest_PT_V_Truth_%i", j), Form("h_Particles_dR_Highest_PT_V_Truth_%i", j), 100, 0, 1);
        h_Particles_dR_Highest_PT_Vz_Truth[j] = new TH1F(Form("h_Particles_dR_Highest_PT_Vz_Truth_%i", j), Form("h_Particles_dR_Highest_PT_Vz_Truth_%i", j), 100, 0, 1);
    }

    TH1D *h_jet_pt_truth_check = new TH1D("h_jet_pt_truth_check", "pT [GeV] ", 270, 0, 2700); // plus Geant4
    h_jet_pt_truth_check->GetXaxis()->SetTitle("P_{T}^{true,jet}");
    //h_jet_pt_truth_check->GetYaxis()->SetTitle("Entries");
    TH1D *h_jet_eta_truth_check = new TH1D("h_jet_eta_truth_check", "#eta", 100, -5, 5); //before eta cut
    TH1D *h_jet_n_truth = new TH1D("h_jet_n_truth", "Nr of truth jets", 10, 0, 10);      //before match
    TH1D *h_jet_nn_truth = new TH1D("h_jet_nn_truth", "Nr of truth jets", 10, 0, 10);    //after match
    TH1D *h_jet_m_truth = new TH1D("jet_m_truth", "Mass [GeV]", 100, 0, 100);
    TH1F *Total_particle_ID_eta_PT_cut = new TH1F("Total_particle_ID_eta_PT_cut", "Total_particle_ID_eta_PT_cut", 20, 0, 20);
    //==============Trailing_particle_ID============//
    TH1F *Timing_detector_dR_Leading_trailing_T = new TH1F("Timing_detector_dR_Leading_trailing_T", "Timing_detector_dR_Leading_trailing_T", 50, 0, 1);
    TH1F *Timing_detector_dR_Leading_next_trailing_T = new TH1F("Timing_detector_dR_Leading_next_trailing_T", "Timing_detector_dR_Leading_next_trailing_T", 50, 0, 1);
    TH1F *Trailing_particle_ID_T = new TH1F("Trailing_particle_ID_T", "Trailing_particle_ID_T", 20, 0, 20);
    TH1F *Trailing_particle_ID_PT = new TH1F("Trailing_particle_ID_PT", "Trailing_particle_ID_PT", 20, 0, 20);
    TH1F *Next_to_trailing_particle_ID_T = new TH1F("Next_to_trailing_particle_ID_T", "Next_to_trailing_particle_ID_T", 20, 0, 20);
    TH1F *Next_to_trailing_particle_ID_PT = new TH1F("Next_to_trailing_particle_ID_PT", "Next_to_trailing_particle_ID_PT", 20, 0, 20);
    TH1F *Timing_detector_Leading = new TH1F("Timing_detector_Leading", "Timing_detector_Leading", 200, 0, 50);
    TH1F *Timing_detector_Trailing = new TH1F("Timing_detector_Trailing", "Timing_detector_Trailing", 200, 0, 50);
    TH1F *Timing_detector_Average = new TH1F("Timing_detector_Average", "Timing_detector_Average", 200, 0, 50);
    TH1F *Timing_detector_next_to_trailing = new TH1F("Timing_detector_next_to_trailing", "Timing_detector_next_to_trailing", 200, 0, 50);
    TH1F *Timing_detector_Trailing_P = new TH1F("Timing_detector_Trailing_P", "Timing_detector_Trailing_P", 200, 0, 100);
    TH1F *Timing_detector_next_to_trailing_P = new TH1F("Timing_detector_next_to_trailing_P", "Timing_detector_next_to_trailing_P", 200, 0, 100);
    TH1F *Timing_detector_Trailing_V = new TH1F("Timing_detector_Trailing_V", "Timing_detector_Trailing_V", 1000, 8, 9);
    TH1F *Timing_detector_next_to_trailing_V = new TH1F("Timing_detector_next_to_trailing_V", "Timing_detector_next_to_trailing_V", 1000, 0.9, 1);
    TH1F *Timing_detector_dR_Leading_trailing_PT = new TH1F("Timing_detector_dR_Leading_trailing_PT", "Timing_detector_dR_Leading_trailing_PT", 50, 0, 1);
    TH1F *Timing_detector_dR_Leading_next_trailing_PT = new TH1F("Timing_detector_dR_Leading_next_trailing_PT", "Timing_detector_dR_Leading_next_trailing_PT", 50, 0, 1);

    Int_t Event;
    Float_t PT_Tr0T_HPt;
    Float_t PT_Tr1T_HPt;
    Float_t PT_Tr2T_HPt;
    Float_t PT_Tr3T_HPt;
    Float_t PT_Tr4T_HPt;
    Float_t PT_Tr0PT_HPt;
    Float_t PT_Tr1PT_HPt;
    Float_t PT_Tr2PT_HPt;
    Float_t PT_Tr3PT_HPt;
    Float_t PT_Tr4PT_HPt;
    Float_t PT_Tr0V_HPt;
    Float_t PT_Tr1V_HPt;
    Float_t PT_Tr2V_HPt;
    Float_t PT_Tr3V_HPt;
    Float_t PT_Tr4V_HPt;
    Float_t PT_Tr0Vz_HPt;
    Float_t PT_Tr1Vz_HPt;
    Float_t PT_Tr2Vz_HPt;
    Float_t PT_Tr3Vz_HPt;
    Float_t PT_Tr4Vz_HPt;
    TTree *T_Truth = new TTree("BDT_variables", "BDT_variables");
    T_Truth->Branch("Event", &Event, "Event/I");
    T_Truth->Branch("PT_Tr0T_HPt", &PT_Tr0T_HPt, "PT_Tr0T_HPt/F");
    T_Truth->Branch("PT_Tr1T_HPt", &PT_Tr1T_HPt, "PT_Tr1T_HPt/F");
    T_Truth->Branch("PT_Tr2T_HPt", &PT_Tr2T_HPt, "PT_Tr2T_HPt/F");
    T_Truth->Branch("PT_Tr3T_HPt", &PT_Tr3T_HPt, "PT_Tr3T_HPt/F");
    T_Truth->Branch("PT_Tr4T_HPt", &PT_Tr4T_HPt, "PT_Tr4T_HPt/F");
    T_Truth->Branch("PT_Tr0PT_HPt", &PT_Tr0PT_HPt, "PT_Tr0PT_HPt/F");
    T_Truth->Branch("PT_Tr1PT_HPt", &PT_Tr1PT_HPt, "PT_Tr1PT_HPt/F");
    T_Truth->Branch("PT_Tr2PT_HPt", &PT_Tr2PT_HPt, "PT_Tr2PT_HPt/F");
    T_Truth->Branch("PT_Tr3PT_HPt", &PT_Tr3PT_HPt, "PT_Tr3PT_HPt/F");
    T_Truth->Branch("PT_Tr4PT_HPt", &PT_Tr4PT_HPt, "PT_Tr4PT_HPt/F");
    T_Truth->Branch("PT_Tr0V_HPt", &PT_Tr0V_HPt, "PT_Tr0V_HPt/F");
    T_Truth->Branch("PT_Tr1V_HPt", &PT_Tr1V_HPt, "PT_Tr1V_HPt/F");
    T_Truth->Branch("PT_Tr2V_HPt", &PT_Tr2V_HPt, "PT_Tr2V_HPt/F");
    T_Truth->Branch("PT_Tr3V_HPt", &PT_Tr3V_HPt, "PT_Tr3V_HPt/F");
    T_Truth->Branch("PT_Tr4V_HPt", &PT_Tr4V_HPt, "PT_Tr4V_HPt/F");
    T_Truth->Branch("PT_Tr0Vz_HPt", &PT_Tr0Vz_HPt, "PT_Tr0Vz_HPt/F");
    T_Truth->Branch("PT_Tr1Vz_HPt", &PT_Tr1Vz_HPt, "PT_Tr1Vz_HPt/F");
    T_Truth->Branch("PT_Tr2Vz_HPt", &PT_Tr2Vz_HPt, "PT_Tr2Vz_HPt/F");
    T_Truth->Branch("PT_Tr3Vz_HPt", &PT_Tr3Vz_HPt, "PT_Tr3Vz_HPt/F");
    T_Truth->Branch("PT_Tr4Vz_HPt", &PT_Tr4Vz_HPt, "PT_Tr4Vz_HPt/F");
    // read detector geometry for this configuration
    //string detector = "./data/rfull009_sifcch7/sifcch7/sifcch7.pandora";
    string detector = "./data/rfull012_sifcch10/sifcch10/sifcch10.pandora";

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

    //calculate the Event pass cut
    Long64_t nPass[20] = {0};

    // calculate position
    double layer_size_mm = 27.5 + 5 + 2.5;
    double hcal_inner_R_mm = 2300;

    double layer_size_ecal_mm = 4.0;
    double ecal_inner_R_mm = 2100;

    // calculate weighted average of time
    double XTime = 0;
    double XEnergy = 0;

    const pandora::SubDetectorType subDetectorType = geom->getPandoraSubDetectorType(caloType);

    double minPtConst = 1.5; // min pT on constituents
    // jets
    double Rparam = 0.4;
    const double ETAmax = 0.6;
    //double ETmin = 0.5;
    if (mtype == 1)
    {
        double ETmin = 200; // increase for boosted
        cout << "mu+mu- mode for boosted jets" << endl;
    }
    cout << "min PT for jets=" << minPtConst << endl;
    cout << "eta max for jets=" << ETAmax << endl;
    cout << "R for jets =" << Rparam << endl;
    // fastjet
    Strategy strategy = fastjet::Best;
    JetDefinition jet_def(fastjet::antikt_algorithm, Rparam, strategy);

    int MaxEvents = 100000000;

    int nEvents = 0;
    int nnEvents = 0;
    // loop over all files
    for (unsigned int mfile = 0; mfile < files.size(); mfile++)
    {
        string Rfile = files[mfile];

        cout << " # File=" << Rfile << endl;

        IO::LCReader *lcReader = IOIMPL::LCFactory::getInstance()->createLCReader();
        lcReader->open(Rfile.c_str());
        bool QQ_event = false;
        bool WW_event = false;
        if ((Rfile).find("_qq_") > 0 and (Rfile).find("_qq_") < 1000)
        {
            cout << "QQ coming!" << endl;
            QQ_event = true;
        }
        if ((Rfile).find("_ww_") > 0 and (Rfile).find("_ww_") < 1000)
        {
            cout << "WW coming!" << endl;
            WW_event = true;
        }
        EVENT::LCEvent *evt = 0;

        if (nEvents > MaxEvents)
            break;
        //----------- the event loop -----------
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

            h_debug->Fill(1.0);
            //create vector
            vector<PseudoJet> avec_truth; // created by generator
            vector<fastjet::PseudoJet> WW_boson;
            // get truth
            IMPL::LCCollectionVec *col = (IMPL::LCCollectionVec *)evt->getCollection("MCParticle");
            int nMCP = col->getNumberOfElements();
            //Particle ID
            int Zprime_pdg = 32;
            int Photon_pdg = 22;
            int Muon_pdg = 13;
            int W_pdg = 24;
            bool Wboson_zero = false;
            //===============================
            //        for Z' -> WW
            //===============================
            for (int i = 0; i < nMCP; ++i)
            {
                EVENT::MCParticle *mcp = (EVENT::MCParticle *)col->getElementAt(i);
                int pdgid = mcp->getPDG();
                double px = mcp->getMomentum()[0];
                double py = mcp->getMomentum()[1];
                double pz = mcp->getMomentum()[2];
                double m = mcp->getMomentum()[3];
                double e = mcp->getEnergy();
                int gs = mcp->getGeneratorStatus();
                bool LWboson = false;
                //WW_boson.clear();
                //================================================
                //                  Check W is hadronic decay
                //================================================
                if (TMath::Abs(pdgid) == W_pdg)
                {
                    if (gs == 2)
                    {
                        for (unsigned int j = 0; j < (mcp->getDaughters().size()); j++)
                        {

                            if ((TMath::Abs(mcp->getDaughters()[j]->getPDG()) < 19) && (TMath::Abs(mcp->getDaughters()[j]->getPDG()) > 10))
                            {
                                LWboson = true;
                            }
                        }
                        if (!LWboson)
                        {
                            WW_boson.push_back(PseudoJet(px, py, pz, e));
                        }
                    }
                } // End W boson ID
            }     // End first MC loop
            Wbosn_nn->Fill(WW_boson.size());
            //===============================
            //Take WW boson 2 (hadronic decay)
            //===============================
            if (WW_event)
            {
                if (WW_boson.size() != 2)
                {
                    continue;
                }
            }
            //===============================
            // only generator level
            //===============================
            for (int i = 0; i < nMCP; ++i)
            {
                EVENT::MCParticle *mcp = (EVENT::MCParticle *)col->getElementAt(i);
                int pdgid = mcp->getPDG();
                double px = mcp->getMomentum()[0];
                double py = mcp->getMomentum()[1];
                double pz = mcp->getMomentum()[2];
                double m = mcp->getMomentum()[3];
                double e = mcp->getEnergy();
                int gs = mcp->getGeneratorStatus();
                //================================================================
                //          exclude Photon decay from muon
                //================================================================
                fastjet::PseudoJet p(px, py, pz, e);
                p.set_user_index(pdgid);
                if (gs == 1)
                {
                    //check photon from mu
                    if (abs(pdgid) == 22 and abs(mcp->getParents()[0]->getPDG()) == 13)
                    {
                        pz = abs(mcp->getParents()[0]->getMomentum()[2]);
                        if (pz > 2500)
                        {
                            continue;
                        }
                    }
                    if (mcp->getCharge() != 0)
                    {
                        if (p.pt() > minPtConst)
                        {
                            avec_truth.push_back(p);
                        }
                    }
                    else
                    {
                        if (p.pt() != 0)
                        {
                            avec_truth.push_back(p);
                        }
                        else
                        {
                            continue;
                        }
                    }
                }
            } // End Second MC loop
            //===================================
            //  Put gs ==1 particle into clusters
            //===================================
            int activeAreaRepeats = 1;
            double ghostArea = 0.01;
            double ghostEtaMax = 7.0;
            fastjet::GhostedAreaSpec fjActiveArea(ghostEtaMax, activeAreaRepeats, ghostArea);
            fastjet::AreaDefinition fjAreaDefinition(fastjet::active_area, fjActiveArea);
            fastjet::ClusterSequenceArea *thisClustering = new fastjet::ClusterSequenceArea(avec_truth, jet_def, fjAreaDefinition);
            vector<fastjet::PseudoJet> sjets_truth = sorted_by_pt(thisClustering->inclusive_jets(25));
            vector<fastjet::PseudoJet> Truthjets;
            vector<vector<double>> jet_time(sjets_truth.size(), vector<double>());
            vector<vector<double>> jet_time_sort(sjets_truth.size(), vector<double>());
            vector<vector<double>> T_sort_number_only(sjets_truth.size(), vector<double>());
            vector<vector<double>> PT_jet(sjets_truth.size(), vector<double>());
            vector<vector<double>> PT_jet_sort(sjets_truth.size(), vector<double>());
            vector<vector<double>> PT_sort_number_only(sjets_truth.size(), vector<double>());
            vector<vector<double>> velocity_jet(sjets_truth.size(), vector<double>());
            vector<vector<double>> velocity_jet_sort(sjets_truth.size(), vector<double>());
            vector<vector<double>> velocity_number_only(sjets_truth.size(), vector<double>());
            vector<vector<double>> velocity_jet_Z(sjets_truth.size(), vector<double>());
            vector<vector<double>> velocity_jet_Z_sort(sjets_truth.size(), vector<double>());
            vector<vector<double>> velocityZ_number_only(sjets_truth.size(), vector<double>());
            vector<vector<TLorentzVector>> FourP_dR_Truth(sjets_truth.size(), vector<TLorentzVector>());
            float mass_average = 0;
            float time_average = 0;
            float event_number = 0;
            float SOL = 3 * TMath::Power(10, 8);
            //==============================================================
            //          Selection: Pt > 25 GeV (include Pt > 1.5)
            //==============================================================
            for (unsigned int k = 0; k < sjets_truth.size(); k++)
            {
                double eta = sjets_truth[k].pseudorapidity();
                double phi = sjets_truth[k].phi();
                if (phi < 0)
                    phi = phi + k2PI;
                double m = sjets_truth[k].m();
                double pt = sjets_truth[k].perp();
                double e = sjets_truth[k].e();
                bool isTrueMatch = false;
                //==============================================================
                // Selection:Jet abs(eta)<1
                //==============================================================
                if (fabs(eta) > 1)
                    continue;
                //===============================================================
                // Matching:  Require Truth jet and W_boson dr < 0.4
                //===============================================================
                TLorentzVector p_using;
                p_using.SetPxPyPzE(sjets_truth[k].px(), sjets_truth[k].py(), sjets_truth[k].pz(), sjets_truth[k].e());
                if (WW_event)
                {
                    for (unsigned int i = 0; i < WW_boson.size(); i++)
                    {
                        TLorentzVector temp_WW(0, 0, 0, 0);
                        temp_WW.SetPxPyPzE(WW_boson[i].px(),
                                           WW_boson[i].py(),
                                           WW_boson[i].pz(),
                                           WW_boson[i].e());
                        TLorentzVector temp_jet(0, 0, 0, 0);
                        temp_jet.SetPxPyPzE(sjets_truth[k].px(),
                                            sjets_truth[k].py(),
                                            sjets_truth[k].pz(),
                                            sjets_truth[k].e());
                        double dr = temp_WW.DeltaR(temp_jet);
                        if (dr < 0.4)
                        {
                            Truthjets.push_back(p_using);
                            isTrueMatch = true;
                            break;
                        }
                    }
                    if (!isTrueMatch)
                    {
                        continue;
                    }
                }
                if (QQ_event)
                {
                    Truthjets.push_back(p_using);
                }
                //================================================================
                //                   Run all Jet constituents
                //================================================================
                vector<PseudoJet> constit = sjets_truth[k].constituents();
                int csize = constit.size();
                for (int i = 0; i < csize; i++)
                {
                    fastjet::PseudoJet constituent(constit[i].px(), constit[i].py(), constit[i].pz(), constit[i].e());
                    TLorentzVector constit_vector;
                    constit_vector.SetPxPyPzE(constit[i].px(), constit[i].py(), constit[i].pz(), constit[i].e());
                    float constit_velocity_vt = TMath::Power((TMath::Power(constit[i].px(), 2) + TMath::Power(constit[i].py(), 2) + TMath::Power(constit[i].pz(), 2)), 0.5) / constit[i].e();
                    float constit_velocity_z = (constit[i].pz() / constit[i].e()); // Magnetic_consideration
                    //Jet Time
                    jet_time[k].push_back(abs(2.3 * TMath::Power(10, 9) / (constit_velocity_z * SOL * (TMath::Tan(constit_vector.Theta())))));
                    jet_time_sort[k].push_back(abs(2.3 * TMath::Power(10, 9) / (constit_velocity_z * SOL * (TMath::Tan(constit_vector.Theta())))));
                    //Velocity Total
                    velocity_jet[k].push_back(constit_velocity_vt);
                    velocity_jet_sort[k].push_back(constit_velocity_vt);
                    //VelocityZ
                    velocity_jet_Z[k].push_back(constit_velocity_z);
                    velocity_jet_Z_sort[k].push_back(constit_velocity_z);
                    //PT
                    PT_jet[k].push_back(constit[i].perp());
                    PT_jet_sort[k].push_back(constit[i].perp());
                    //Four_momentum
                    FourP_dR_Truth[k].push_back(constit_vector);
                }
            } //end of sjets_truth loop
            h_jet_n_truth->Fill(sjets_truth.size());
            h_jet_nn_truth->Fill(Truthjets.size());
            //================================================================
            //                  sort PT & time & velocity
            //================================================================
            for (unsigned int k = 0; k < Truthjets.size(); k++)
            {
                if (PT_jet_sort[k].size() > 0)
                {
                    sort(PT_jet_sort[k].begin(), PT_jet_sort[k].end());
                    sort(jet_time_sort[k].begin(), jet_time_sort[k].end());
                    sort(velocity_jet_sort[k].begin(), velocity_jet_sort[k].end());
                    sort(velocity_jet_Z_sort[k].begin(), velocity_jet_Z_sort[k].end());
                    //================================================================
                    // Finally output PT_sort_number_only == PT_Reco_sort size
                    //================================================================
                    PT_sort_number_only[k].push_back(PT_jet_sort[k][0]);
                    T_sort_number_only[k].push_back(jet_time_sort[k][0]);
                    velocity_number_only[k].push_back(velocity_jet_sort[k][0]);
                    velocityZ_number_only[k].push_back(velocity_jet_Z_sort[k][0]);
                    //===================================================
                    // Only Save different kind of PT, T, Vt
                    //===================================================
                    for (unsigned int j = 0; j < PT_jet_sort[k].size(); j++)
                    {
                        if (PT_jet_sort[k][j] != PT_sort_number_only[k].back())
                        {
                            PT_sort_number_only[k].push_back(PT_jet_sort[k][j]);
                        }
                        if (jet_time_sort[k][j] != T_sort_number_only[k].back())
                        {
                            T_sort_number_only[k].push_back(jet_time_sort[k][j]);
                        }
                        if (velocity_jet_sort[k][j] != velocity_number_only[k].back())
                        {
                            velocity_number_only[k].push_back(velocity_jet_sort[k][j]);
                        }
                        if (velocity_jet_Z_sort[k][j] != velocityZ_number_only[k].back())
                        {
                            velocityZ_number_only[k].push_back(velocity_jet_Z_sort[k][j]);
                        }
                    }
                }
            } //End of Truthjet
            vector<float> dR_Highest_PT_T_Truth;
            vector<float> dR_Highest_PT_PT_Truth;
            vector<float> dR_Highest_PT_V_Truth;
            vector<float> dR_Highest_PT_Vz_Truth;
            //================================================================
            //      Save Highest PT Four momentum
            //================================================================
            vector<vector<TLorentzVector>> Highest_PT_FourP(sjets_truth.size(), vector<TLorentzVector>());
            vector<vector<int>> PT_PDG_Reco(sjets_truth.size(), vector<int>());
            vector<vector<int>> T_PDG_Reco(sjets_truth.size(), vector<int>());
            vector<vector<int>> V_PDG_Reco(sjets_truth.size(), vector<int>());
            for (unsigned int k = 0; k < Truthjets.size(); k++)
            {
                unsigned int Size_T_PT_V_Vz = jet_time[k].size();

                for (unsigned int i = 0; i < Size_T_PT_V_Vz; i++)
                {
                    if (PT_jet_sort[k][Size_T_PT_V_Vz - 1] == PT_jet[k][i])
                    {
                        Highest_PT_FourP[k].push_back(FourP_dR_Truth[k][i]);
                    }
                }
                for (unsigned int j = 0; j < 5; j++)
                {
                    if (j < PT_sort_number_only[k].size())
                    {
                        for (unsigned int i = 0; i < Size_T_PT_V_Vz; i++)
                        {
                            if (PT_sort_number_only[k][j] == PT_jet[k][i])
                            {
                                dR_Highest_PT_PT_Truth.push_back(FourP_dR_Truth[k][i].DeltaR(Highest_PT_FourP[k][0]));
                                h_Particles_dR_Highest_PT_PT_Truth[j]->Fill(FourP_dR_Truth[k][i].DeltaR(Highest_PT_FourP[k][0]));
                            }
                        }
                    }

                    if (j < T_sort_number_only[k].size())
                    {
                        for (unsigned int i = 0; i < Size_T_PT_V_Vz; i++)
                        {
                            if (T_sort_number_only[k][T_sort_number_only[k].size() - 1 - j] == jet_time[k][i])
                            {
                                dR_Highest_PT_T_Truth.push_back(FourP_dR_Truth[k][i].DeltaR(Highest_PT_FourP[k][0]));
                                h_Particles_dR_Highest_PT_T_Truth[j]->Fill(FourP_dR_Truth[k][i].DeltaR(Highest_PT_FourP[k][0]));
                            }
                        }
                    }
                    if (j < velocity_number_only[k].size())
                    {
                        for (unsigned int i = 0; i < Size_T_PT_V_Vz; i++)
                        {
                            if (velocity_number_only[k][j] == velocity_jet[k][i])
                            {
                                dR_Highest_PT_V_Truth.push_back(FourP_dR_Truth[k][i].DeltaR(Highest_PT_FourP[k][0]));
                                h_Particles_dR_Highest_PT_V_Truth[j]->Fill(FourP_dR_Truth[k][i].DeltaR(Highest_PT_FourP[k][0]));
                            }
                        }
                    }
                    if (j < velocityZ_number_only[k].size())
                    {
                        for (unsigned int i = 0; i < Size_T_PT_V_Vz; i++)
                        {
                            if (velocityZ_number_only[k][j] == velocity_jet_Z[k][i])
                            {
                                dR_Highest_PT_Vz_Truth.push_back(FourP_dR_Truth[k][i].DeltaR(Highest_PT_FourP[k][0]));
                                h_Particles_dR_Highest_PT_Vz_Truth[j]->Fill(FourP_dR_Truth[k][i].DeltaR(Highest_PT_FourP[k][0]));
                            }
                        }
                    }
                }
            }

            //================================================================
            //    Save Tree Variables
            //================================================================
            //PT
            if (dR_Highest_PT_PT_Truth.size() >= 1)
            {
                PT_Tr0PT_HPt = dR_Highest_PT_PT_Truth[0];
            }
            if (dR_Highest_PT_PT_Truth.size() >= 2)
            {
                PT_Tr1PT_HPt = dR_Highest_PT_PT_Truth[1];
            }
            if (dR_Highest_PT_PT_Truth.size() >= 3)
            {
                PT_Tr2PT_HPt = dR_Highest_PT_PT_Truth[2];
            }
            if (dR_Highest_PT_PT_Truth.size() >= 4)
            {
                PT_Tr3PT_HPt = dR_Highest_PT_PT_Truth[3];
            }
            if (dR_Highest_PT_PT_Truth.size() >= 5)
            {
                PT_Tr4PT_HPt = dR_Highest_PT_PT_Truth[4];
            }
            //T
            if (dR_Highest_PT_T_Truth.size() >= 1)
            {
                PT_Tr0T_HPt = dR_Highest_PT_T_Truth[0];
            }
            if (dR_Highest_PT_T_Truth.size() >= 2)
            {
                PT_Tr1T_HPt = dR_Highest_PT_T_Truth[1];
            }
            if (dR_Highest_PT_T_Truth.size() >= 3)
            {
                PT_Tr2T_HPt = dR_Highest_PT_T_Truth[2];
            }
            if (dR_Highest_PT_T_Truth.size() >= 4)
            {
                PT_Tr3T_HPt = dR_Highest_PT_T_Truth[3];
            }
            if (dR_Highest_PT_T_Truth.size() >= 5)
            {
                PT_Tr4T_HPt = dR_Highest_PT_T_Truth[4];
            }
            //Vt
            if (dR_Highest_PT_V_Truth.size() >= 1)
            {
                PT_Tr0V_HPt = dR_Highest_PT_V_Truth[0];
            }
            if (dR_Highest_PT_V_Truth.size() >= 2)
            {
                PT_Tr1V_HPt = dR_Highest_PT_V_Truth[1];
            }
            if (dR_Highest_PT_V_Truth.size() >= 3)
            {
                PT_Tr2V_HPt = dR_Highest_PT_V_Truth[2];
            }
            if (dR_Highest_PT_V_Truth.size() >= 4)
            {
                PT_Tr3V_HPt = dR_Highest_PT_V_Truth[3];
            }
            if (dR_Highest_PT_V_Truth.size() >= 5)
            {
                PT_Tr4V_HPt = dR_Highest_PT_V_Truth[4];
            }
            if (dR_Highest_PT_Vz_Truth.size() >= 1)
            {
                PT_Tr0Vz_HPt = dR_Highest_PT_Vz_Truth[0];
            }
            if (dR_Highest_PT_Vz_Truth.size() >= 2)
            {
                PT_Tr1Vz_HPt = dR_Highest_PT_Vz_Truth[1];
            }
            if (dR_Highest_PT_Vz_Truth.size() >= 3)
            {
                PT_Tr2Vz_HPt = dR_Highest_PT_Vz_Truth[2];
            }
            if (dR_Highest_PT_Vz_Truth.size() >= 4)
            {
                PT_Tr3Vz_HPt = dR_Highest_PT_Vz_Truth[3];
            }
            if (dR_Highest_PT_Vz_Truth.size() >= 5)
            {
                PT_Tr4Vz_HPt = dR_Highest_PT_Vz_Truth[4];
            }
            //Timing_detector_Leading->Fill(abs(2.3 * TMath::Power(10, 9) / (Vz_Leading * SOL * TMath::Tan(Theta_Leading))));
            //Timing_detector_Trailing->Fill(abs(2.3 * TMath::Power(10, 9) / (Vz_Trailing * SOL * TMath::Tan(Theta_Trailing))));

            //===============================
            //  Clear the Truth Level vector
            //===============================
            sjets_truth.clear();
            Truthjets.clear();
            WW_boson.clear();
            avec_truth.clear();
            T_Truth->Fill();
        }
        lcReader->close();
        delete lcReader;
    }
    RootFile->Write();
    h_jet_n_truth->Write();
    Timing_Standard->Write();
    h_jet_pt_truth_check->Write();
    h_jet_eta_truth_check->Write();
    h_jet_n_truth->Write();
    h_jet_m_truth->Write();
    h_jet_nn_truth->Write();
    //RootFile->Print();
    RootFile->Close();
    return 0;
}
