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
    TH1F *Timing_Standard = new TH1F("Timing_Standard", "Timing_Standard", 200, 0, 50);
    TH1F *h_Particles_dR_Highest_PT_T[5];
    TH1F *h_Particles_dR_Highest_PT_PT[5];
    for (int j = 0; j < 5; j++)
    {
        h_Particles_dR_Highest_PT_T[j] = new TH1F(Form("h_Particles_dR_Highest_PT_T_%i", j), Form("h_Particles_dR_Highest_PT_T_%i", j), 100, 0, 1);
        h_Particles_dR_Highest_PT_PT[j] = new TH1F(Form("h_Particles_dR_Highest_PT_PT_%i", j), Form("h_Particles_dR_Highest_PT_PT_%i", j), 100, 0, 1);
    }
    TH2F *h_Particles_Rank_T_vs_T[5];
    TH2F *h_Particles_Rank_PT_vs_T[5];
    for (int j = 0; j < 5; j++)
    {
        h_Particles_Rank_T_vs_T[j] = new TH2F(Form("h_Particles_Rank_T_vs_T%i", j), Form("h_Particles_Rank_T_vs_T%i", j), 20, 0, 20, 500, 0, 50);
        h_Particles_Rank_PT_vs_T[j] = new TH2F(Form("h_Particles_Rank_PT_vs_T%i", j), Form("h_Particles_Rank_PT_vs_T%i", j), 20, 0, 20, 500, 0, 50);
    }

    TH2F *h_Particles_Rank_T_vs_PT[5];
    TH2F *h_Particles_Rank_PT_vs_PT[5];
    for (int j = 0; j < 5; j++)
    {
        h_Particles_Rank_T_vs_PT[j] = new TH2F(Form("h_Particles_Rank_T_vs_PT%i", j), Form("h_Particles_Rank_T_vs_PT%i", j), 20, 0, 20, 16, -2, 6);
        h_Particles_Rank_PT_vs_PT[j] = new TH2F(Form("h_Particles_Rank_PT_vs_PT%i", j), Form("h_Particles_Rank_PT_vs_PT%i", j), 20, 0, 20, 16, -2, 6);
    }
    TH1F *h_Particles_Rank_T[5];
    TH1F *h_Particles_Rank_PT[5];
    for (int j = 0; j < 5; j++)
    {
        h_Particles_Rank_T[j] = new TH1F(Form("h_Particles_Rank_T_%i", j), Form("h_Particles_Rank_T_%i", j), 20, 0, 20);
        h_Particles_Rank_PT[j] = new TH1F(Form("h_Particles_Rank_PT_%i", j), Form("h_Particles_Rank_PT_%i", j), 20, 0, 20);
    }
    TH1D *h_jet_pt_truth_check = new TH1D("h_jet_pt_truth_check", "pT [GeV] ", 270, 0, 2700); // plus Geant4
    h_jet_pt_truth_check->GetXaxis()->SetTitle("P_{T}^{true,jet}");
    //h_jet_pt_truth_check->GetYaxis()->SetTitle("Entries");
    TH1D *h_jet_eta_truth_check = new TH1D("h_jet_eta_truth_check", "#eta", 100, -5, 5); //before eta cut
    TH1D *h_jet_n_truth = new TH1D("jet_n_truth", "Nr of truth jets", 10, 0, 10);        //before match
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
    TTree *T = new TTree("BDT_variables", "BDT_variables");
    T->Branch("Event", &Event, "Event/I");
    T->Branch("PT_Tr0T_HPt", &PT_Tr0T_HPt, "PT_Tr0T_HPt/F");
    T->Branch("PT_Tr1T_HPt", &PT_Tr1T_HPt, "PT_Tr1T_HPt/F");
    T->Branch("PT_Tr2T_HPt", &PT_Tr2T_HPt, "PT_Tr2T_HPt/F");
    T->Branch("PT_Tr3T_HPt", &PT_Tr3T_HPt, "PT_Tr3T_HPt/F");
    T->Branch("PT_Tr4T_HPt", &PT_Tr4T_HPt, "PT_Tr4T_HPt/F");
    T->Branch("PT_Tr0PT_HPt", &PT_Tr0PT_HPt, "PT_Tr0PT_HPt/F");
    T->Branch("PT_Tr1PT_HPt", &PT_Tr1PT_HPt, "PT_Tr1PT_HPt/F");
    T->Branch("PT_Tr2PT_HPt", &PT_Tr2PT_HPt, "PT_Tr2PT_HPt/F");
    T->Branch("PT_Tr3PT_HPt", &PT_Tr3PT_HPt, "PT_Tr3PT_HPt/F");
    T->Branch("PT_Tr4PT_HPt", &PT_Tr4PT_HPt, "PT_Tr4PT_HPt/F");

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
    vector<int> Total_particle_kind = {11, 12, 13, 14, 130, 211, 310, 321, 2112, 2212, 3112, 3122, 3312, 3222, 3322, 16, 3334};
    vector<int> Trailing_particle_kind = {11, 13, 130, 211, 321, 2112, 2212, 3122, 3112, 3312, 3222, 3322, 16, 3334, 1000010020, 310};
    vector<int> Trailing_particle_kind_limit = {11, 13, 130, 211, 321, 2112, 2212};
    vector<int> Trailing_particle_kind_limit_pT = {11, 13, 211, 321, 2212};
    // loop over all files
    for (unsigned int mfile = 0; mfile < files.size(); mfile++)
    {
        string Rfile = files[mfile];

        cout << " # File=" << Rfile << endl;

        IO::LCReader *lcReader = IOIMPL::LCFactory::getInstance()->createLCReader();
        lcReader->open(Rfile.c_str());

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
            vector<int> PDG_with_no_charge = {0};
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
            //for Z'->WW
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
                //WW_boson.clear();
                //================================================
                //                  Check W is hadronic decay
                //================================================
                if (mcp->getParents().size() != 0)
                {
                    if (gs == 2)
                    {
                        if (abs(pdgid) == 24)
                        {
                            for (unsigned int j = 0; j < (mcp->getDaughters().size()); j++)
                            {
                                if ((abs(mcp->getDaughters()[j]->getPDG()) < 19) and (abs(mcp->getDaughters()[j]->getPDG()) > 10))
                                {
                                    continue;
                                }
                            }
                            PseudoJet p(px, py, pz, e);
                            //TLorentzVector p;
                            //p.SetPxPyPzE(mcp->getMomentum()[0], mcp->getMomentum()[1], mcp->getMomentum()[2], mcp->getEnergy());
                            WW_boson.push_back(p);
                        }
                    }

                    //================================================================
                    //          exclude Photon decay from muon
                    //================================================================
                    fastjet::PseudoJet p(px, py, pz, e);
                    p.set_user_index(pdgid);
                    if (gs == 1)
                    {
                        if (mcp->getCharge() == 0)
                        {

                            int itt;
                            itt = find(PDG_with_no_charge.begin(), PDG_with_no_charge.end(), abs(pdgid))[0]; //retuen [0] element
                            if (itt != abs(pdgid))
                            {
                                //cout << "abs(pdgid) =" << abs(pdgid) << endl;
                                PDG_with_no_charge.push_back(abs(pdgid));
                            }
                        }
                        //check photon from mu
                        if (abs(pdgid) == 22 and abs(mcp->getParents()[0]->getPDG()) == 13)
                        {
                            pz = abs(mcp->getParents()[0]->getMomentum()[2]);
                            if (pz > 2500)
                            {
                                continue;
                            }
                        }
                        //remove photon
                        if (abs(pdgid) == 22)
                        {
                            continue;
                        }
                        //Pt < 1.5 GeV can not reach ECAL
                        if (p.pt() > minPtConst)
                        {
                            avec_truth.push_back(p);
                        }
                    }
                }
            }
            //cout << WW_boson.size() << endl;
            //===================================
            //  Put gs ==1 particle into clusters
            //===================================
            int activeAreaRepeats = 1;
            double ghostArea = 0.01;
            double ghostEtaMax = 7.0;
            fastjet::GhostedAreaSpec fjActiveArea(ghostEtaMax, activeAreaRepeats, ghostArea);
            fastjet::AreaDefinition fjAreaDefinition(fastjet::active_area, fjActiveArea);
            fastjet::ClusterSequenceArea *thisClustering = new fastjet::ClusterSequenceArea(avec_truth, jet_def, fjAreaDefinition);
            //==============================================================
            //          Selection: Pt > 15 GeV (include Pt > 1.5)
            //==============================================================
            vector<fastjet::PseudoJet> sjets_truth = sorted_by_pt(thisClustering->inclusive_jets(15));
            vector<LParticle> truthjets;
            vector<fastjet::PseudoJet> Truthjets_axis;
            vector<int> Truthjets_axis_index;
            bool isTrueMatch = true;
            //cout << WW_boson.size() << endl;
            for (unsigned int k = 0; k < sjets_truth.size(); k++)
            {
                double eta = sjets_truth[k].pseudorapidity();
                double phi = sjets_truth[k].phi();
                if (phi < 0)
                    phi = phi + k2PI;
                double m = sjets_truth[k].m();
                double pt = sjets_truth[k].perp();
                double e = sjets_truth[k].e();
                h_jet_pt_truth_check->Fill(pt);
                h_jet_eta_truth_check->Fill(eta);
                h_jet_n_truth->Fill(sjets_truth.size());
                h_jet_m_truth->Fill(m);
                //==============================================================
                // Selection: Each Jet constituents & abs(eta)<1
                //==============================================================
                if (TMath::Abs(eta) > 1)
                    continue;
                //===============================================================
                // Matching:  Require Truth jet and W_boson dr < 0.4
                //===============================================================
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

                    //cout << "temp_WW =" << temp_WW[0].size() << endl;
                    //cout << "dr =" << dr << endl;
                    if (dr < 0.4)
                    {
                        TLorentzVector p_using;
                        p_using.SetPxPyPzE(sjets_truth[k].px(), sjets_truth[k].py(), sjets_truth[k].pz(), sjets_truth[k].e());
                        Truthjets_axis.push_back(p_using);
                        Truthjets_axis_index.push_back(k);
                        isTrueMatch = false;
                        //cout << "k =" << k << endl;
                    }
                }
                if (isTrueMatch)
                {
                    continue;
                }
                h_jet_nn_truth->Fill(Truthjets_axis_index.size());
                //================================================================
                //                   Run all Jet constituents
                //================================================================
                vector<float> jet_time;
                vector<float> jet_time_sort;
                vector<float> jet_time_for_rank_sort;
                vector<float> PT_jet;
                vector<float> PT_jet_sort;
                vector<float> momentum_jet;
                vector<float> jet_P_for_rank;
                vector<float> jet_time_for_rank;
                vector<int> Rank_PDGID;
                vector<int> constit_PDG;
                vector<float> velocity_jet;
                vector<float> velocity_jet_sort;
                vector<float> velocity_jet_Z;
                vector<float> velocity_jet_Theta;
                vector<TLorentzVector> FourP;
                float mass_average = 0;
                float time_average = 0;
                float event_number = 0;
                float SOL = 3 * TMath::Power(10, 8);
                vector<PseudoJet> constit = sjets_truth[k].constituents();
                int csize = constit.size();
                for (int i = 0; i < csize; i++)
                {
                    fastjet::PseudoJet constituent(constit[i].px(), constit[i].py(), constit[i].pz(), constit[i].e());
                    TLorentzVector constit_vec;
                    constit_vec.SetPxPyPzE(constit[i].px(), constit[i].py(), constit[i].pz(), constit[i].e());
                    int it;
                    it = find(Total_particle_kind.begin(), Total_particle_kind.end(), abs(constit[i].user_index()))[0];
                    if (it != abs(constit[i].user_index()))
                    {
                        Total_particle_kind.push_back(abs(constit[i].user_index()));
                    }
                    for (unsigned int m = 0; m < Total_particle_kind.size(); m++)
                    {
                        if (abs(constit[i].user_index()) == Total_particle_kind[m])
                            Total_particle_ID_no_cut->Fill(m);
                    }
                    //================================================================
                    // Require each Jet constituent abs(eta) < 1 & cut off pt < 1.5 GeV
                    //================================================================
                    if (abs(constit_vec.Eta()) > 1)
                    {
                        continue;
                    }
                    if (constit[i].perp() < 1.5)
                    {
                        continue;
                    }
                    int itt;
                    itt = find(Total_particle_kind.begin(), Total_particle_kind.end(), abs(constit[i].user_index()))[0];
                    if (itt != abs(constit[i].user_index()))
                    {
                        Total_particle_kind.push_back(abs(constit[i].user_index()));
                    }
                    for (unsigned int m = 0; m < Total_particle_kind.size(); m++)
                    {
                        if (abs(constit[i].user_index()) == Total_particle_kind[m])
                            Total_particle_ID_eta_PT_cut->Fill(m);
                    }

                    mass_average = mass_average + constit_vec.M();
                    float constit_velocity_vt = TMath::Power((TMath::Power(constit[i].px(), 2) + TMath::Power(constit[i].py(), 2) + TMath::Power(constit[i].pz(), 2)), 0.5) / constit[i].e();
                    float constit_velocity_z = (constit[i].pz() / constit[i].e());                                                                 // Magnetic_consideration
                    time_average = time_average + abs(2.3 * TMath::Power(10, 9) / (constit_velocity_z * SOL * (TMath::Tan(constit_vec.Theta())))); //[ns]
                    jet_time.push_back(abs(2.3 * TMath::Power(10, 9) / (constit_velocity_z * SOL * (TMath::Tan(constit_vec.Theta())))));
                    jet_time_sort.push_back(abs(2.3 * TMath::Power(10, 9) / (constit_velocity_z * SOL * (TMath::Tan(constit_vec.Theta())))));
                    jet_time_for_rank.push_back(abs(2.3 * TMath::Power(10, 9) / (constit_velocity_z * SOL * (TMath::Tan(constit_vec.Theta())))));
                    //jet_P_for_rank.push_back(constit_vec.P());
                    Rank_PDGID.push_back(sjets_truth[k].constituents()[i].user_index());
                    velocity_jet.push_back(constit_velocity_vt);
                    velocity_jet_sort.push_back(constit_velocity_vt);
                    velocity_jet_Z.push_back(constit_velocity_z);
                    velocity_jet_Theta.push_back(constit_vec.Theta());
                    Timing_Standard->Fill(abs(2.3 * TMath::Power(10, 9) / (SOL * TMath::Sin(constit_vec.Theta())))); //Suppose all of them are photons.
                    PT_jet.push_back(constit[i].perp());
                    PT_jet_sort.push_back(constit[i].perp());
                    momentum_jet.push_back(constit_vec.P());
                    constit_PDG.push_back(abs(sjets_truth[k].constituents()[i].user_index()));
                    FourP.push_back(constit_vec);
                    event_number = event_number + 1;
                }
                double max_time = *max_element(jet_time.begin(), jet_time.end());
                double min_time = *min_element(jet_time.begin(), jet_time.end());
                double max_perp = *max_element(PT_jet.begin(), PT_jet.end());
                double min_perp = *min_element(PT_jet.begin(), PT_jet.end());
                //================================================================
                //                  sort PT & time & velocity
                //================================================================
                sort(PT_jet_sort.begin(), PT_jet_sort.end());
                sort(jet_time.begin(), jet_time.end());
                sort(jet_time_sort.begin(), jet_time_sort.end());
                sort(velocity_jet.begin(), velocity_jet.end());
                sort(velocity_jet_sort.begin(), velocity_jet_sort.end());
                sort(velocity_jet_Z.begin(), velocity_jet_Z.end());
                sort(jet_time_for_rank.begin(), jet_time_for_rank.end());
                //sort(jet_P_for_rank.begin(), jet_P_for_rank.end());
                int Trailing_ID_T = 0;          //One jet one trailing ID(T)
                int Trailing_ID_PT = 0;         //One jet one trailing ID(PT)
                int Next_to_trailing_ID_T = 0;  //One jet one trailing ID(T)
                int Next_to_trailing_ID_PT = 0; //One jet one trailing ID(PT)
                float Momentum_Trailing = 0;
                float Momentum_Next_to_Trailing = 0;
                //float PT_Trailing = 0;
                float Theta_Leading = 0;
                float Theta_Trailing = 0;
                float Theta_Next_to_Trailing = 0;
                float Vz_Trailing = 0;
                float Vz_Next_to_Trailing = 0;
                float Vz_Leading = 0;
                vector<TLorentzVector> HighestPT_Trailing_and_next_trailing;
                vector<float> Particle_ID_T_PT_Re;
                vector<int> Particle_ID_T;
                vector<int> Particle_ID_T_Re;
                vector<int> Particle_ID_T_T;
                vector<int> Particle_ID_T_T_Re;
                vector<float> Particle_ID_T_PT;
                vector<int> Particle_ID_PT;
                vector<int> Particle_ID_PT_Re;
                vector<float> dR_Highest_PT_T;
                vector<float> dR_Highest_PT_PT;
                vector<float> PT_checking;
                vector<float> Timing_checking;
                vector<float> Particle_ID_PT_T;
                vector<float> Particle_ID_PT_PT;
                vector<float> Particle_ID_PT_PT_Re;
                vector<float> Particle_ID_PT_T_Re;
                vector<float> velocity_checking;
                vector<float> H_Particle_ID_PT;
                vector<int> H_Particle_ID_T;
                int check_timing_number = 0;
                int check_PT_number = 0;
                //================================================================
                //      Save Highest PT Four momentum
                //================================================================
                for (unsigned int i = 0; i < PT_jet.size(); i++)
                {
                    if (max_perp == PT_jet[i])
                    {
                        HighestPT_Trailing_and_next_trailing.push_back(FourP[i]);
                    }
                }
                //================================================================
                //    PT order from small to largest :: Save Variables
                //================================================================
                for (unsigned int j = 0; j < PT_jet.size(); j++)
                {
                    for (unsigned int i = 0; i < PT_jet.size(); i++)
                    {
                        if (PT_jet_sort[j] == PT_jet[i])
                        {
                            //cout << "10" << endl;
                            PT_checking.push_back(PT_jet[i]);
                            check_PT_number = check_PT_number + 1;
                            Particle_ID_PT_PT.push_back(PT_jet[i]);
                            Particle_ID_PT_T.push_back(jet_time[i]);
                            Particle_ID_PT.push_back(constit_PDG[i]);
                            dR_Highest_PT_PT.push_back(HighestPT_Trailing_and_next_trailing[0].DeltaR(FourP[i]));
                            if (j < 5)
                            {
                                h_Particles_dR_Highest_PT_PT[j]->Fill(HighestPT_Trailing_and_next_trailing[0].DeltaR(FourP[i]));
                            }
                        }
                    }
                }
                for (unsigned int j = 0; j < jet_time.size(); j++)
                {
                    for (unsigned int i = 0; i < jet_time.size(); i++)
                    {
                        if (velocity_jet_sort[j] == velocity_jet[i]) //Velocity 由小排到大
                        {
                            check_timing_number = check_timing_number + 1;
                            Timing_checking.push_back(jet_time[i]);
                            velocity_checking.push_back(velocity_jet[i]);
                            Particle_ID_T_PT.push_back(PT_jet[i]);
                            Particle_ID_T_T.push_back(jet_time[i]);
                            Particle_ID_T.push_back(constit_PDG[i]);
                            //cout << "1 =" << HighestPT_Trailing_and_next_trailing[0].DeltaR(FourP[i]);
                            dR_Highest_PT_T.push_back(HighestPT_Trailing_and_next_trailing[0].DeltaR(FourP[i]));
                            if (j < 5)
                            {
                                h_Particles_dR_Highest_PT_T[j]->Fill(HighestPT_Trailing_and_next_trailing[0].DeltaR(FourP[i]));
                            }
                        }
                    }
                }
                int Num_of_Tra_pT = 0;
                for (unsigned int j = 0; j < jet_time.size(); j++)
                {
                    if ((event_number - 1) < j)
                    {
                        continue;
                    }
                    else
                    {
                        for (unsigned int m = 0; m < Trailing_particle_kind_limit_pT.size(); m++)
                        {
                            if (abs(Particle_ID_PT[j]) == Trailing_particle_kind_limit_pT[m] and Num_of_Tra_pT < 5)
                            {
                                Num_of_Tra_pT = Num_of_Tra_pT + 1;
                                //cout << "Input TMath::Log10(Particle_ID_PT_PT[j]):: " << TMath::Log10(Particle_ID_PT_PT[j]) << endl;
                                //cout << "Input (Particle_ID_PT_PT[j]):: " << (Particle_ID_PT_PT[j]) << endl;
                                Particle_ID_PT_PT_Re.push_back((Particle_ID_PT_PT[j]));
                                Particle_ID_PT_T_Re.push_back(Particle_ID_PT_T[j]);
                                Particle_ID_PT_Re.push_back(Particle_ID_PT[j]);
                            }
                        }
                    }
                }
                if (Num_of_Tra_pT != 5)
                {
                    continue;
                }
                int Num_of_Tra_V = 0;
                for (unsigned int j = 0; j < jet_time.size(); j++)
                {
                    if ((event_number - 1) < j)
                    {
                        continue;
                    }
                    if ((Particle_ID_T_PT[j]) == Particle_ID_PT_PT_Re[0] or (Particle_ID_T_PT[j]) == Particle_ID_PT_PT_Re[1] or (Particle_ID_T_PT[j]) == Particle_ID_PT_PT_Re[2] or (Particle_ID_T_PT[j]) == Particle_ID_PT_PT_Re[3] or (Particle_ID_T_PT[j]) == Particle_ID_PT_PT_Re[4])
                    {
                        continue;
                    }
                    else
                    {
                        for (unsigned int m = 0; m < Trailing_particle_kind_limit.size(); m++)
                        {
                            if (abs(Particle_ID_T[j]) == Trailing_particle_kind_limit[m] and Num_of_Tra_V < 5)
                            {
                                //cout << "Used_j: " << j << endl;
                                //cout << "Output (Particle_ID_T_PT[j]): " << (Particle_ID_T_PT[j]) << endl;
                                //cout << "Output TMath::Log10(Particle_ID_T_PT[j]): " << TMath::Log10(Particle_ID_T_PT[j]) << endl;
                                Num_of_Tra_V = Num_of_Tra_V + 1;
                                Particle_ID_T_PT_Re.push_back((Particle_ID_T_PT[j]));
                                Particle_ID_T_T_Re.push_back(Particle_ID_T_T[j]);
                                Particle_ID_T_Re.push_back(Particle_ID_T[j]);
                            }
                        }
                    }
                }
                if (Num_of_Tra_V != 5)
                {
                    continue;
                }
                for (int j = 0; j < 5; j++)
                {
                    if ((event_number - 1) < j)
                    {
                        continue;
                    }
                    else
                    {
                        for (unsigned int m = 0; m < Trailing_particle_kind_limit_pT.size(); m++)
                        {
                            if (abs(Particle_ID_PT_Re[j]) == Trailing_particle_kind_limit_pT[m])
                            {
                                H_Particle_ID_PT.push_back(m);
                                h_Particles_Rank_PT[j]->Fill(m);

                                h_Particles_Rank_PT_vs_T[j]->Fill(m, Particle_ID_PT_T_Re[j]);
                                if (TMath::Log10(Particle_ID_PT_PT[j]) > -1.5)
                                    h_Particles_Rank_PT_vs_PT[j]->Fill(m, TMath::Log10(Particle_ID_PT_PT_Re[j]));
                                if (TMath::Log10(Particle_ID_PT_PT[j]) < -1.5)
                                    h_Particles_Rank_PT_vs_PT[j]->Fill(m, -1.499);
                            }
                        }
                    }
                }
                for (int j = 0; j < 5; j++)
                {
                    if ((event_number - 1) < j)
                    {
                        continue;
                    }
                    else
                    {
                        {
                            for (unsigned int m = 0; m < Trailing_particle_kind_limit.size(); m++)
                            {
                                if (abs(Particle_ID_T_Re[j]) == Trailing_particle_kind_limit[m])
                                {
                                    H_Particle_ID_T.push_back(m);
                                    h_Particles_Rank_T[j]->Fill(m);
                                    h_Particles_Rank_T_vs_T[j]->Fill(m, Particle_ID_T_T_Re[j]);
                                    if (TMath::Log10(Particle_ID_T_PT[j]) > -1.5)
                                        h_Particles_Rank_T_vs_PT[j]->Fill(m, TMath::Log10(Particle_ID_T_PT_Re[j]));
                                    if (TMath::Log10(Particle_ID_T_PT[j]) < -1.5)
                                        h_Particles_Rank_T_vs_PT[j]->Fill(m, -1.499);
                                }
                            }
                        }
                    }
                }
                //================================================================
                //  Save Truth Level Tree
                //================================================================
                if (event_number >= 1)
                {
                    PT_Tr0T_HPt = TMath::Log10(Particle_ID_T_PT_Re[0]);
                    PT_Tr0PT_HPt = TMath::Log10(Particle_ID_PT_PT_Re[0]);
                }
                if (event_number >= 2)
                {
                    PT_Tr1T_HPt = TMath::Log10(Particle_ID_T_PT_Re[1]);
                    PT_Tr1PT_HPt = TMath::Log10(Particle_ID_PT_PT_Re[1]);
                }
                if (event_number >= 3)
                {
                    PT_Tr2T_HPt = TMath::Log10(Particle_ID_T_PT_Re[2]);
                    PT_Tr2PT_HPt = TMath::Log10(Particle_ID_PT_PT_Re[2]);
                }
                if (event_number >= 4)
                {
                    PT_Tr3T_HPt = TMath::Log10(Particle_ID_T_PT_Re[3]);
                    PT_Tr3PT_HPt = TMath::Log10(Particle_ID_PT_PT_Re[3]);
                }
                if (event_number >= 5)
                {
                    PT_Tr4T_HPt = TMath::Log10(Particle_ID_T_PT_Re[4]);
                    PT_Tr4PT_HPt = TMath::Log10(Particle_ID_PT_PT_Re[4]);
                }
                for (unsigned int i = 0; i < jet_time.size(); i++)
                {
                    if (jet_time_sort[jet_time.size() - 1] == jet_time[i])
                    {
                        Trailing_ID_T = constit_PDG[i];
                    }
                    if (jet_time_sort[jet_time.size() - 2] == jet_time[i])
                    {
                        //cout << "debug2" << endl;
                        Next_to_trailing_ID_T = constit_PDG[i];
                    }
                    if (PT_jet_sort[0] == PT_jet[i])
                    {
                        Trailing_ID_PT = constit_PDG[i];
                    }
                    if (PT_jet_sort[1] == PT_jet[i])
                    {
                        Next_to_trailing_ID_PT = constit_PDG[i];
                    }
                }
                int it;
                it = find(Trailing_particle_kind.begin(), Trailing_particle_kind.end(), abs(Trailing_ID_T))[0];
                //cout << "it = " << it << endl;
                //cout << "abs(Trailing_ID_T) = " << abs(Trailing_ID_T) << endl;
                if (it != abs(Trailing_ID_T))
                {
                    Trailing_particle_kind.push_back(abs(Trailing_ID_T));
                }
                else
                {
                    //cout << "" << endl;
                }
                for (unsigned int m = 0; m < Trailing_particle_kind.size(); m++)
                {
                    if (abs(Trailing_ID_T) == Trailing_particle_kind[m])
                        Trailing_particle_ID_T->Fill(m);
                }
                int it2;
                it2 = find(Trailing_particle_kind.begin(), Trailing_particle_kind.end(), abs(Trailing_ID_PT))[0];
                //cout << "it2 = " << it2 << endl;
                //cout << " abs(Trailing_ID_PT) = " << abs(Trailing_ID_PT) << endl;
                if (it2 != abs(Trailing_ID_PT))
                {
                    Trailing_particle_kind.push_back(abs(Trailing_ID_PT));
                }
                else
                {
                    //cout << "" << endl;
                }
                for (unsigned int m = 0; m < Trailing_particle_kind.size(); m++)
                {
                    if (abs(Trailing_ID_PT) == Trailing_particle_kind[m])
                        Trailing_particle_ID_PT->Fill(m);
                }
                int it4;
                it4 = find(Trailing_particle_kind.begin(), Trailing_particle_kind.end(), abs(Next_to_trailing_ID_T))[0];
                //cout << "it4 = " << it4 << endl;
                //cout << " abs(Next_to_trailing_ID_T) = " << abs(Next_to_trailing_ID_T) << endl;
                if (it4 != abs(Next_to_trailing_ID_T))
                {
                    Trailing_particle_kind.push_back(abs(Next_to_trailing_ID_T));
                }
                else
                {
                    //cout << "" << endl;
                }
                for (unsigned int m = 0; m < Trailing_particle_kind.size(); m++)
                {
                    if (abs(Next_to_trailing_ID_T) == Trailing_particle_kind[m])
                        Next_to_trailing_particle_ID_T->Fill(m);
                }
                //cout << "P = " << endl;
                int it3;
                it3 = find(Trailing_particle_kind.begin(), Trailing_particle_kind.end(), abs(Next_to_trailing_ID_PT))[0];
                //cout << "it3 = " << it3 << endl;
                //cout << " abs(Next_to_trailing_ID_PT) = " << abs(Next_to_trailing_ID_PT) << endl;
                if (it3 != abs(Next_to_trailing_ID_PT))
                {
                    Trailing_particle_kind.push_back(abs(Next_to_trailing_ID_PT));
                }
                for (unsigned int m = 0; m < Trailing_particle_kind.size(); m++)
                {
                    //cout << "Trailing_particle_kind[m]: "<< abs(Trailing_particle_kind[m]) << endl;
                    if (abs(Next_to_trailing_ID_PT) == Trailing_particle_kind[m])
                        Next_to_trailing_particle_ID_PT->Fill(m);
                    //cout << "Next_to_trailing_particle_ID_PT->GetBinContent(m): " << Next_to_trailing_particle_ID_PT->GetBinContent(m+1) << endl;
                }
                for (unsigned int i = 0; i < jet_time.size(); i++)
                {
                    if (jet_time_sort[0] == jet_time[i])
                    {
                        Vz_Leading = velocity_jet_Z[i];
                        Theta_Leading = velocity_jet_Theta[i];
                    }
                }
                for (unsigned int i = 0; i < jet_time.size(); i++)
                {
                    if (jet_time_sort[jet_time.size() - 1] == jet_time[i])
                    {
                        Theta_Trailing = velocity_jet_Theta[i];
                        Vz_Trailing = velocity_jet_Z[i];
                        Momentum_Trailing = momentum_jet[i];
                        HighestPT_Trailing_and_next_trailing.push_back(FourP[i]);
                    }
                }
                for (unsigned int i = 0; i < jet_time.size(); i++)
                {
                    if (jet_time_sort[jet_time.size() - 2] == jet_time[i])
                    {
                        Theta_Next_to_Trailing = velocity_jet_Theta[i];
                        Vz_Next_to_Trailing = velocity_jet_Z[i];
                        Momentum_Next_to_Trailing = momentum_jet[i];
                        HighestPT_Trailing_and_next_trailing.push_back(FourP[i]);
                    }
                }
                for (unsigned int i = 0; i < PT_jet.size(); i++)
                {
                    if (PT_jet_sort[0] == PT_jet[i])
                    {
                        HighestPT_Trailing_and_next_trailing.push_back(FourP[i]);
                    }
                }
                for (unsigned int i = 0; i < PT_jet.size(); i++)
                {
                    if (PT_jet_sort[1] == PT_jet[i])
                    {
                        HighestPT_Trailing_and_next_trailing.push_back(FourP[i]);
                    }
                }
                Timing_detector_Average->Fill(abs(time_average));
                Timing_detector_Leading->Fill(abs(2.3 * TMath::Power(10, 9) / (Vz_Leading * SOL * TMath::Tan(Theta_Leading))));
                Timing_detector_Trailing->Fill(abs(2.3 * TMath::Power(10, 9) / (Vz_Trailing * SOL * TMath::Tan(Theta_Trailing))));
                Timing_detector_next_to_trailing->Fill(abs(2.3 * TMath::Power(10, 9) / (Vz_Next_to_Trailing * SOL * TMath::Tan(Theta_Next_to_Trailing))));
                //cout << "P = " << abs(Momentum_Trailing) << endl;
                //Timing_detector_Trailing_P->Fill(abs(Momentum_Trailing));
                Timing_detector_next_to_trailing_P->Fill(abs(Momentum_Next_to_Trailing));
                //cout << "logV =" << 3 * TMath::Power(10, 8) * abs(velocity_jet_sort[0]) << endl;
                //cout << "logV =" << TMath::Log10(3 * TMath::Power(10, 8) * abs(velocity_jet_sort[0])) << endl;
                Timing_detector_Trailing_V->Fill(TMath::Log10(3 * TMath::Power(10, 8) * abs(velocity_jet_sort[0])));
                //cout << "Trailing-V: " << TMath::Log10(3 * TMath::Power(10, 8) * abs(velocity_jet_sort[0])) << endl;
                Timing_detector_next_to_trailing_V->Fill(abs(velocity_jet_sort[1]));
                Timing_detector_dR_Leading_trailing_T->Fill(HighestPT_Trailing_and_next_trailing[0].DeltaR(HighestPT_Trailing_and_next_trailing[1]));
                Timing_detector_dR_Leading_next_trailing_T->Fill(HighestPT_Trailing_and_next_trailing[0].DeltaR(HighestPT_Trailing_and_next_trailing[2]));
                Timing_detector_dR_Leading_trailing_PT->Fill(HighestPT_Trailing_and_next_trailing[0].DeltaR(HighestPT_Trailing_and_next_trailing[3]));
                Timing_detector_dR_Leading_next_trailing_PT->Fill(HighestPT_Trailing_and_next_trailing[0].DeltaR(HighestPT_Trailing_and_next_trailing[4]));
                //velocity_jet.clear();
            } //end of sjets_truth loop
            for (unsigned int i = 0; i < Total_particle_kind.size(); i++)
            {
                Total_particle_ID_no_cut->GetXaxis()->SetBinLabel(i + 1, Form("ID=%i", Total_particle_kind[i]));
                Total_particle_ID_eta_PT_cut->GetXaxis()->SetBinLabel(i + 1, Form("ID=%i", Total_particle_kind[i]));
                Trailing_particle_ID_T->GetXaxis()->SetBinLabel(i + 1, Form("ID=%i", Total_particle_kind[i]));
                Next_to_trailing_particle_ID_T->GetXaxis()->SetBinLabel(i + 1, Form("ID=%i", Total_particle_kind[i]));
                Next_to_trailing_particle_ID_PT->GetXaxis()->SetBinLabel(i + 1, Form("ID=%i", Total_particle_kind[i]));
            }
            for (int j = 0; j < 5; j++)
            {
                for (unsigned int i = 0; i < Trailing_particle_kind_limit.size(); i++)
                {
                    h_Particles_Rank_PT[j]->GetXaxis()->SetBinLabel(i + 1, Form("ID=%i", Trailing_particle_kind_limit[i]));
                    h_Particles_Rank_T[j]->GetXaxis()->SetBinLabel(i + 1, Form("ID=%i", Trailing_particle_kind_limit[i]));
                }
            }
            //===============================
            //  Clear the Truth Level vector
            //===============================
            sjets_truth.clear();
            truthjets.clear();
            Truthjets_axis_index.clear();
            WW_boson.clear();
            avec_truth.clear();
            PDG_with_no_charge.clear();
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
