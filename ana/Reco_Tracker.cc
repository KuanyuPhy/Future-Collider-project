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
#include "TNtuple.h"
#include "TTree.h"

using namespace std;

const double kPI = TMath::Pi();
const double k2PI = 2 * kPI;
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include <fastjet/AreaDefinition.hh>
#include <fastjet/ClusterSequenceArea.hh>
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

    TTree *mytree = new TTree("mytree", "this is a plain tree");
    //Create the histogram
    TH1F *Event_number_check = new TH1F("Event_number_check", "Event_number_check", 12, 0, 12);
    TFile *RootFile = new TFile(outputfile.c_str(), "UPDATE", "Histogram file");
    TH1F *h_Jet_PT = new TH1F("h_Jet_PT", "h_Jet_PT", 100, 0, 3000);
    h_Jet_PT->GetXaxis()->SetTitle("GeV");
    TH1F *Total_particle_ID_no_cut = new TH1F("Total_particle_ID_no_cut", "Total_particle_ID_no_cut", 20, 0, 20);
    TH1F *Eta_plot = new TH1F("Eta_plot", "Eta_plot", 200, -10, 10);
    TH1D *h_pt_clus = new TH1D("clus_pt", "pt RecoClus", 400, 0, 1000);
    TH1D *h_eta_clus = new TH1D("clus_eta", "eta RecoClus", 120, -6, 6);
    TH1D *h_pt_ecal = new TH1D("ecal_pt", "pt ecal", 400, 0, 1000);
    double TTxbins[] = {0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 6.0, 8.0, 10.0, 20.0, 30};
    const int TTBins = sizeof(TTxbins) / sizeof(double);
    TH1D *binsTT = new TH1D("bins_time_ecal", "bins_time_ecal", TTBins - 1, TTxbins);
    binsTT->Sumw2();
    for (Int_t j = 0; j < TTBins - 1; j++)
    {
        float x = TTxbins[j + 1] - TTxbins[j];
        binsTT->Fill(TTxbins[j] + 0.5 * x, x);
    }
    TH1F *Timing_detecto_ECAL_TDif = new TH1F("Timing_detecto_ECAL_TDif", "Timing_detecto_ECAL_TDif", 200, 0, 50);
    //TH1F *Timing_detecto_ECAL_TDif = new TH1F("Timing_detecto_ECAL_TDif", "Timing_detecto_ECAL_TDif", TTBins - 1, TTxbins);

    //TH1D *h_jet_hitstime1D_ECAL = new TH1D("jet_hitstime_Ecal_1D", "1D for time vs hit in HCAL", 100, 0.5, 10);

    TH1F *mass_sum_average_Reco = new TH1F("mass_sum_average_Reco", "mass_sum_average_Reco", 200, 0.5, 4.5);
    TH1F *Timing_detector_Reco_TOF = new TH1F("Timing_detector_Reco_TOF", "Timing_detector_Reco_TOF", 200, 0, 50);
    TH1F *Timing_detector_Reco_TOF_track = new TH1F("Timing_detector_Reco_TOF_track", "Timing_detector_Reco_TOF_track", 200, 0, 50);

    TH1F *h_Particles_Rank_T_Reco[5];
    TH1F *h_Particles_Rank_PT_Reco[5];
    for (int j = 0; j < 5; j++)
    {
        h_Particles_Rank_T_Reco[j] = new TH1F(Form("h_Particles_Rank_T_Reco_%i", j), Form("h_Particles_Rank_T_Reco_%i", j), 30, 0, 30);
        h_Particles_Rank_PT_Reco[j] = new TH1F(Form("h_Particles_Rank_PT_Reco_%i", j), Form("h_Particles_Rank_PT_Reco_%i", j), 30, 0, 30);
    }
    TH1F *h_Particles_dR_Highest_PT_T_Reco[5];
    TH1F *h_Particles_dR_Highest_PT_PT_Reco[5];
    for (int j = 0; j < 5; j++)
    {
        h_Particles_dR_Highest_PT_T_Reco[j] = new TH1F(Form("h_Particles_dR_Highest_PT_T_Reco_%i", j), Form("h_Particles_dR_Highest_PT_T_Reco_%i", j), 200, 0, 1);
        h_Particles_dR_Highest_PT_PT_Reco[j] = new TH1F(Form("h_Particles_dR_Highest_PT_PT_Reco_%i", j), Form("h_Particles_dR_Highest_PT_PT_Reco_%i", j), 200, 0, 1);
    }
    // read detector geometry for this configuration

    Int_t Event;
    Int_t ID_Tr0T;
    Int_t ID_Tr1T;
    Int_t ID_Tr2T;
    Int_t ID_Tr3T;
    Int_t ID_Tr4T;
    Int_t ID_Tr0PT;
    Int_t ID_Tr1PT;
    Int_t ID_Tr2PT;
    Int_t ID_Tr3PT;
    Int_t ID_Tr4PT;
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
    Int_t Event_reco;
    Float_t dR_Tr0T_HPt_Reco;
    Float_t dR_Tr1T_HPt_Reco;
    Float_t dR_Tr2T_HPt_Reco;
    Float_t dR_Tr3T_HPt_Reco;
    Float_t dR_Tr4T_HPt_Reco;
    Float_t dR_Tr0PT_HPt_Reco;
    Float_t dR_Tr1PT_HPt_Reco;
    Float_t dR_Tr2PT_HPt_Reco;
    Float_t dR_Tr3PT_HPt_Reco;
    Float_t dR_Tr4PT_HPt_Reco;
    TTree *T_Reco_T = new TTree("BDT_variables_Reco", "BDT_variables_Reco");
    T_Reco_T->Branch("Event_reco", &Event_reco, "Event_reco/I");
    T_Reco_T->Branch("dR_Tr0T_HPt_Reco", &dR_Tr0T_HPt_Reco, "dR_Tr0T_HPt_Reco/F");
    T_Reco_T->Branch("dR_Tr1T_HPt_Reco", &dR_Tr1T_HPt_Reco, "dR_Tr1T_HPt_Reco/F");
    T_Reco_T->Branch("dR_Tr2T_HPt_Reco", &dR_Tr2T_HPt_Reco, "dR_Tr2T_HPt_Reco/F");
    T_Reco_T->Branch("dR_Tr3T_HPt_Reco", &dR_Tr3T_HPt_Reco, "dR_Tr3T_HPt_Reco/F");
    T_Reco_T->Branch("dR_Tr4T_HPt_Reco", &dR_Tr4T_HPt_Reco, "dR_Tr4T_HPt_Reco/F");
    T_Reco_T->Branch("dR_Tr0PT_HPt_Reco", &dR_Tr0PT_HPt_Reco, "dR_Tr0PT_HPt_Reco/F");
    T_Reco_T->Branch("dR_Tr1PT_HPt_Reco", &dR_Tr1PT_HPt_Reco, "dR_Tr1PT_HPt_Reco/F");
    T_Reco_T->Branch("dR_Tr2PT_HPt_Reco", &dR_Tr2PT_HPt_Reco, "dR_Tr2PT_HPt_Reco/F");
    T_Reco_T->Branch("dR_Tr3PT_HPt_Reco", &dR_Tr3PT_HPt_Reco, "dR_Tr3PT_HPt_Reco/F");
    T_Reco_T->Branch("dR_Tr4PT_HPt_Reco", &dR_Tr4PT_HPt_Reco, "dR_Tr4PT_HPt_Reco/F");

    Int_t Event_reco_track;
    Float_t dR_Tr0T_HPt_Reco_track;
    Float_t dR_Tr1T_HPt_Reco_track;
    Float_t dR_Tr2T_HPt_Reco_track;
    Float_t dR_Tr3T_HPt_Reco_track;
    Float_t dR_Tr4T_HPt_Reco_track;
    Float_t dR_Tr0PT_HPt_Reco_track;
    Float_t dR_Tr1PT_HPt_Reco_track;
    Float_t dR_Tr2PT_HPt_Reco_track;
    Float_t dR_Tr3PT_HPt_Reco_track;
    Float_t dR_Tr4PT_HPt_Reco_track;
    TTree *T_Reco_T_track = new TTree("BDT_variables_Reco_track", "BDT_variables_Reco_track");
    T_Reco_T_track->Branch("Event_reco_track", &Event_reco_track, "Event_reco_track/I");
    T_Reco_T_track->Branch("dR_Tr0T_HPt_Reco_track", &dR_Tr0T_HPt_Reco_track, "dR_Tr0T_HPt_Reco_track/F");
    T_Reco_T_track->Branch("dR_Tr1T_HPt_Reco_track", &dR_Tr1T_HPt_Reco_track, "dR_Tr1T_HPt_Reco_track/F");
    T_Reco_T_track->Branch("dR_Tr2T_HPt_Reco_track", &dR_Tr2T_HPt_Reco_track, "dR_Tr2T_HPt_Reco_track/F");
    T_Reco_T_track->Branch("dR_Tr3T_HPt_Reco_track", &dR_Tr3T_HPt_Reco_track, "dR_Tr3T_HPt_Reco_track/F");
    T_Reco_T_track->Branch("dR_Tr4T_HPt_Reco_track", &dR_Tr4T_HPt_Reco_track, "dR_Tr4T_HPt_Reco_track/F");
    T_Reco_T_track->Branch("dR_Tr0PT_HPt_Reco_track", &dR_Tr0PT_HPt_Reco_track, "dR_Tr0PT_HPt_Reco_track/F");
    T_Reco_T_track->Branch("dR_Tr1PT_HPt_Reco_track", &dR_Tr1PT_HPt_Reco_track, "dR_Tr1PT_HPt_Reco_track/F");
    T_Reco_T_track->Branch("dR_Tr2PT_HPt_Reco_track", &dR_Tr2PT_HPt_Reco_track, "dR_Tr2PT_HPt_Reco_track/F");
    T_Reco_T_track->Branch("dR_Tr3PT_HPt_Reco_track", &dR_Tr3PT_HPt_Reco_track, "dR_Tr3PT_HPt_Reco_track/F");
    T_Reco_T_track->Branch("dR_Tr4PT_HPt_Reco_track", &dR_Tr4PT_HPt_Reco_track, "dR_Tr4PT_HPt_Reco_track/F");

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
    //TNtupleD *tuple = new TNtupleD("MCParticle", "", tupleNames);
    TNtuple *tuple = new TNtuple("tuple", "a n-tuple", "px:py:pz:vx:vy:vz:ex:ey:ez:pdg:np:nd:gs:ss");

    int var[13];

    double minPtConst = 0.2; // min pT on constituents
    // jets
    double Rparam = 0.4;
    const double ETAmax = 0.6;
    //double ETmin = 0.5;
    if (mtype == 1)
    {
        double ETmin = 200; // increase for boosted
        cout << "mu+mu- mode for boosted jets" << endl;
    }
    // fastjet
    Strategy strategy = fastjet::Best;
    JetDefinition jet_def(fastjet::antikt_algorithm, Rparam, strategy);

    int MaxEvents = 100000000;

    int nEvents = 0;
    int nnEvents = 0;
    int Hadronic_decay_total = 0;
    int Leptonic_decay_total = 0;
    vector<int> Total_particle_kind = {11, 12, 13, 14, 22, 130, 211, 310, 321, 2112, 2212, 3112, 3122, 3312, 3222, 3322, 16, 3334};
    vector<int> Trailing_particle_kind = {11, 13, 130, 211, 321, 2112, 2212, 3122, 3112, 3312, 3222, 3322, 16, 3334, 1000010020, 310};
    vector<int> Trailing_particle_kind_limit = {11, 13, 130, 211, 321, 2112, 2212};
    vector<int> Trailing_particle_kind_limit_pT = {11, 13, 211, 321, 2212};
    vector<int> Trailing_particle_kind_T = {11, 12, 13, 14, 22, 130, 211, 310, 321, 2112, 2212, 3112, 3122, 3312, 3222, 3322, 16, 3334, 1000010020};
    vector<int> Trailing_particle_kind_PT = {11, 12, 13, 14, 22, 130, 211, 310, 321, 2112, 2212, 3112, 3122, 3312, 3222, 3322, 16, 3334, 1000010020};
    // loop over all files
    for (unsigned int mfile = 0; mfile < files.size(); mfile++)
    {
        string Rfile = files[mfile];

        cout << " # File=" << Rfile << endl;

        IO::LCReader *lcReader = IOIMPL::LCFactory::getInstance()->createLCReader();
        lcReader->open(Rfile.c_str());
        cout << "File_name" << Rfile.c_str() << endl;

        EVENT::LCEvent *evt = 0;
        if (nEvents > MaxEvents)
            break;
        //----------- the event loop -----------
        nnEvents = 0;
        while ((evt = lcReader->readNextEvent()) != 0)
        {
            if (nEvents == 0)
                UTIL::LCTOOLS::dumpEvent(evt);
            nnEvents++;
            nEvents++;
            int Status1 = 0;
            Event_number_check->Fill(1);
            if (nEvents > MaxEvents)
                break;
            std::string mcpName("MCParticle");
            // get truth
            IMPL::LCCollectionVec *col = (IMPL::LCCollectionVec *)evt->getCollection(mcpName);
            int nMCP = col->getNumberOfElements();

            int neu = 0;
            //list of Particle
            vector<int> PDG_with_no_charge = {0};
            vector<PseudoJet> avec_truth;     // created by generator
            vector<PseudoJet> avec_truth_sim; // also created by geant
            int Forth = 1;
            int Back = -1;
            vector<int> Check_Forth_And_Back = {Forth, Back};
            vector<bool> Check_Forth_And_Back_Bool;

            vector<TLorentzVector> Forth_And_Back_Vector;
            int Zprime_pdg = 32;
            int Photon_pdg = 22;
            int Muon_pdg = 13;
            int W_pdg = 24;

            //Check the leptonic decay and hadronic decay
            for (int CFAB = 0; CFAB < 2; CFAB++)
            {
                int Leptonic_check = 0;

                for (int i = 0; i < nMCP; ++i)
                {
                    EVENT::MCParticle *mcp = (EVENT::MCParticle *)col->getElementAt(i);
                    if (mcp->getParents().size() != 0)
                    {
                        if ((mcp->getParents()[0]->getPDG() == 32) and (mcp->getPDG()) * Check_Forth_And_Back[CFAB] > 0)
                        {
                            if (mcp->getGeneratorStatus() == 3)
                            {
                                TLorentzVector p;
                                p.SetPxPyPzE(mcp->getMomentum()[0], mcp->getMomentum()[1], mcp->getMomentum()[2], mcp->getMomentum()[3]);
                                Forth_And_Back_Vector.push_back(p);
                                for (unsigned int j = 0; j < (mcp->getDaughters().size()); j++)
                                {
                                    if ((abs(mcp->getDaughters()[j]->getPDG()) < 19) and (abs(mcp->getDaughters()[j]->getPDG()) > 10))
                                    {
                                        Leptonic_check = Leptonic_check + 1;
                                    }
                                }
                            }
                        }
                    }
                }
                if (Leptonic_check != 0)
                {
                    Leptonic_decay_total = Leptonic_decay_total + 1;
                    //cout << "Awful leptonic decay:" << endl;
                    Check_Forth_And_Back_Bool.push_back(false);
                }
                else
                {
                    Hadronic_decay_total = Hadronic_decay_total + 1;
                    //cout << "Good Jet! Hadronic decay:" << endl;
                    Check_Forth_And_Back_Bool.push_back(true);
                }
                //cout << "Check_Forth_And_Back_Bool size = " << Check_Forth_And_Back_Bool.size() << endl;
                //cout << "nMCP =" << nMCP << endl;
                if (Check_Forth_And_Back_Bool[0] == false and Check_Forth_And_Back_Bool[1] == false)
                    continue;
                //mytree->Fill();
            }

            for (int i = 0; i < nMCP; ++i)
            {

                EVENT::MCParticle *mcp = (EVENT::MCParticle *)col->getElementAt(i);
                double px = mcp->getMomentum()[0];
                double py = mcp->getMomentum()[1];
                double pz = mcp->getMomentum()[2];
                double m = mcp->getMomentum()[3];
                double e = sqrt(px * px + py * py + pz * pz + m * m);
                int pdgid = mcp->getPDG();
                //cout << "m=" << m << " |p|=" << sqrt(px*px+py*py+pz*pz) << " e=" << e << endl;
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

                fastjet::PseudoJet p(px, py, pz, e);
                //p.SetPxPyPzE(px, py, pz, e);
                p.set_user_index(pdgid);
                double ptt_r = p.pt();
                /*
                if (ptt_r == 0)
                {
                    cout << "ptt_r =" << ptt_r << endl;
                }
                else
                {
                    cout << "ptt_r =" << ptt_r << endl;
                }
                */
                //cout << "test1" << endl;
                // only generator level

                int pdg = mcp->getPDG();
                int np = mcp->getParents().size();
                int nd = mcp->getDaughters().size();
                int gs = mcp->getGeneratorStatus();
                if (abs(pdg) == 22 and gs == 1 and abs(mcp->getParents()[0]->getPDG()) == 13)
                {
                    pz = abs(mcp->getParents()[0]->getMomentum()[2]);
                    if (pz > 2000)
                    {
                        continue;
                    }
                }
                if (gs == 1)
                {
                    Status1 = Status1 + 1;
                    avec_truth.push_back(p);
                    if (avec_truth.size() == 0)
                    {
                        cout << "avec_truth size =" << avec_truth.size() << endl;
                    }
                }
            }
            int activeAreaRepeats = 1;
            double ghostArea = 0.01;
            double ghostEtaMax = 7.0;
            fastjet::GhostedAreaSpec fjActiveArea(ghostEtaMax, activeAreaRepeats, ghostArea);
            //cout << "test7" << endl;
            fastjet::AreaDefinition fjAreaDefinition(fastjet::active_area, fjActiveArea);
            //cout << "test8" << endl;
            fastjet::ClusterSequenceArea *thisClustering = new fastjet::ClusterSequenceArea(avec_truth, jet_def, fjAreaDefinition);
            //cout << "test9" << endl;
            vector<fastjet::PseudoJet> sjets_truth = sorted_by_pt(thisClustering->inclusive_jets(25.0));
            vector<LParticle> truthjets;
            int nn = 0;
            int check = 4;
            int Jet_each_event = 0;
            vector<TLorentzVector> Truthjets_axis;
            for (unsigned int k = 0; k < sjets_truth.size(); k++)
            {
                if (sjets_truth[k].constituents().size() == 1 and abs(sjets_truth[k].constituents()[0].user_index()) == 22)
                {
                    continue;
                }
                double eta = sjets_truth[k].pseudorapidity();
                double phi = sjets_truth[k].phi();
                if (phi < 0)
                    phi = phi + k2PI;
                double m = sjets_truth[k].m();
                double pt = sjets_truth[k].perp();
                double e = sjets_truth[k].e();
                vector<float> velocity_jet;
                vector<float> velocity_jet_sort;
                vector<float> momentum_jet;
                vector<float> PT_jet;
                vector<float> PT_jet_sort;
                vector<int> constit_PDG;
                vector<float> velocity_jet_Z;
                vector<float> velocity_jet_Theta;
                vector<float> jet_time;
                vector<float> jet_time_sort;
                vector<float> P_jet_sort;
                vector<float> jet_time_for_rank_sort;
                vector<float> jet_time_for_rank;
                vector<float> jet_P_for_rank_sort;
                vector<float> jet_P_for_rank;
                vector<int> Rank_PDGID;
                TLorentzVector p_using;

                p_using.SetPxPyPzE(sjets_truth[k].px(), sjets_truth[k].py(), sjets_truth[k].pz(), sjets_truth[k].e());
                //===============
                fastjet::PseudoJet Jet_axis(sjets_truth[k].px(), sjets_truth[k].py(), sjets_truth[k].pz(), sjets_truth[k].e());
                //require jet delta R is less than 0.4
                if (p_using.DeltaR(Forth_And_Back_Vector[0]) < 0.1 * check and Check_Forth_And_Back_Bool[0] == 0)
                {
                    //cout << "backward jet1" << endl;
                    continue;
                }
                if (p_using.DeltaR(Forth_And_Back_Vector[1]) < 0.1 * check and Check_Forth_And_Back_Bool[1] == 0)
                {
                    //cout << "backward jet2" << endl;
                    continue;
                }
                if (p_using.DeltaR(Forth_And_Back_Vector[0]) > 0.1 * check and p_using.DeltaR(Forth_And_Back_Vector[1]) > 0.1 * check)
                {
                    continue;
                }
                float Jet_PT = p_using.Perp();
                Jet_each_event = Jet_each_event + 1;
                h_Jet_PT->Fill(p_using.Perp());
                LParticle p(sjets_truth[k].px(), sjets_truth[k].py(), sjets_truth[k].pz(), sjets_truth[k].e(), 0);
                vector<PseudoJet> constit = sjets_truth[k].constituents();
                int csize = constit.size();
                float event_number = 0;
                float Event_number_out_Eta2P1 = 0;
                float time_average = 0;
                float mass_average = 0;
                float SOL = 3 * TMath::Power(10, 8);
                int Trailing_particle_ID_size = Trailing_particle_kind.size();
                int Check_photon_jet = 0;
                vector<TLorentzVector> FourP;
                vector<PseudoJet> FourP_1;

                for (int i = 0; i < csize; i++)
                {
                    if ((constit[i].user_index() - 22) != 0)
                        Check_photon_jet = Check_photon_jet + 1;
                }
                if (Check_photon_jet == 0)
                {
                    for (int i = 0; i < csize; i++)
                    {
                        continue;
                    }
                }
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
                }
                Eta_plot->Fill(p_using.Eta());
                // fill jet substructure
                p.SetParameter(EffectiveRadius(sjets_truth[k], sjets_truth[k].constituents(), Rparam));  // 0
                p.SetParameter(nsubjettiness(sjets_truth[k], sjets_truth[k].constituents(), 1, Rparam)); // 1
                p.SetParameter(nsubjettiness(sjets_truth[k], sjets_truth[k].constituents(), 2, Rparam)); // 2
                p.SetParameter(nsubjettiness(sjets_truth[k], sjets_truth[k].constituents(), 3, Rparam)); //  3
                p.SetParameter(splittingscale(sjets_truth[k]));                                          // 4
                p.SetParameter(eccentricity(sjets_truth[k], sjets_truth[k].constituents()));             // 5
                truthjets.push_back(p);
                Truthjets_axis.push_back(p_using);
                if (Truthjets_axis.size() <= 0)
                {
                    continue;
                }
                //cout << "Truthjets_axis =" << Truthjets_axis.size() << endl;
            }

            // clusters
            vector<PseudoJet> avec_clus;
            double calo_sum = 0;
            IMPL::LCCollectionVec *col5 = (IMPL::LCCollectionVec *)evt->getCollection("ReconClusters");
            int nCL = col5->getNumberOfElements();
            for (int i = 0; i < nCL; ++i)
            {
                //cout << "nCL=" << nCL << endl;
                EVENT::Cluster *mcp = (EVENT::Cluster *)col5->getElementAt(i);
                const float *pos = mcp->getPosition();
                float x = pos[0];
                float y = pos[1];
                float z = pos[2];
                double e = mcp->getEnergy();
                double _tmp = std::sqrt(x * x + y * y + z * z);
                double px = e * x / _tmp;
                double py = e * y / _tmp;
                double pz = e * z / _tmp;

                calo_sum = calo_sum + e;
                PseudoJet pj(px, py, pz, e);
                double eta_r = pj.pseudorapidity();
                double phi_r = pj.phi();
                double pt_r = pj.pt();
                // fill clusters
                h_pt_clus->Fill(pt_r);
                h_eta_clus->Fill(eta_r);

                //cout << " e=" << e <<  " phi=" << pj.phi() << " eta=" << pj.eta() << endl;
                if (pt_r > minPtConst)
                {
                    avec_clus.push_back(pj);
                }
                //cout << "Calo clus sum=" << calo_sum << endl;
            }
            vector<PseudoJet> avec_hits_raw;
            vector<PseudoJet> avec_hits_raw_sf;
            vector<LParticle> simhits;
            double ecalsum_raw = 0;

            // ECAL simulated raw hits
            IMPL::LCCollectionVec *col53 = (IMPL::LCCollectionVec *)evt->getCollection("EcalBarrelHits");
            nCL = col53->getNumberOfElements();
            for (int i = 0; i < nCL; ++i)
            {
                //cout << "nCL=" << nCL << endl;
                EVENT::SimCalorimeterHit *mcp = (EVENT::SimCalorimeterHit *)col53->getElementAt(i);
                const float *pos = mcp->getPosition();
                float x = pos[0];
                float y = pos[1];
                float z = pos[2];
                double e = mcp->getEnergy();
                double _tmp = std::sqrt(x * x + y * y + z * z);
                double px = e * x / _tmp;
                double py = e * y / _tmp;
                double pz = e * z / _tmp;
                // Get the two 32-bit chunks of the ID.
                int cellId0 = mcp->getCellID0();
                int cellId1 = mcp->getCellID1();
                // Make a 64-bit id for the IDDecoder.  The type MUST be "long long" and not "long".  (from Tony Johnson)
                long long cellId = ((long long)cellId1) << 32 | cellId0;
                int layer = decoder_ecal->getFieldValue("layer", cellId);
                // 1st layer on middle
                double layer_pos_cm = 0.1 * (layer_size_ecal_mm * 0.5 + (layer * layer_size_ecal_mm));
                double Thit = mcp->getTimeCont(0);
                // number of MC contributions to the hit
                // calculate average time  in [ns]
                double avt = 0;
                double ave = 0;
                for (int jj = 0; jj < mcp->getNMCContributions(); jj++)
                {
                    avt = avt + mcp->getEnergyCont(jj) * mcp->getTimeCont(jj);
                    ave = ave + mcp->getEnergyCont(jj);
                }
                avt = avt / ave;
                ecalsum_raw = ecalsum_raw + e;
                PseudoJet pj(px, py, pz, e);
                double eta_r = pj.pseudorapidity();
                double phi_r = pj.phi();
                double pt_r = pj.pt();
                avec_hits_raw.push_back(pj);
                //cout << " e=" << e << " phi=" << pj.phi() << " eta=" << pj.eta() << endl;
                // fill hits
                LParticle p(px, py, pz, e, layer);
                p.SetCharge(layer);
                p.SetType(2); // ECAL
                p.SetStatus(mcp->getNMCContributions());
                p.SetParameter(x * 1000); //*1000???
                p.SetParameter(y * 1000);
                p.SetParameter(z * 1000);
                p.SetParameter(layer_pos_cm);
                // find fastest hit
                float timeCont = mcp->getTimeCont(0);
                EVENT::MCParticle *pmm = mcp->getParticleCont(0);
                int pdg = pmm->getPDG();
                int status = pmm->getSimulatorStatus();
                float rawTime = timeCont;
                float FAM = pmm->getMass();
                int nCont = mcp->getNMCContributions();
                for (int jjj = 0; jjj < nCont; jjj++)
                {
                    if (mcp->getTimeCont(jjj) < rawTime)
                    {
                        rawTime = mcp->getTimeCont(jjj);
                    }
                    EVENT::MCParticle *pmm = mcp->getParticleCont(jjj);
                    pdg = pmm->getPDG();
                    //cout << "pdg =" << pdg << endl;
                    status = pmm->getSimulatorStatus();
                }
                p.SetParameter(rawTime); // fastest hit
                p.SetParameter(avt);     // average hit time
                p.SetParameter((double)pdg);
                p.SetParameter((double)status);
                p.SetParameter((double)FAM);
                // PDG of first constribution
                if (mcp->getNMCContributions() > 0)
                {
                    EVENT::MCParticle *pmm = mcp->getParticleCont(0);
                    int pdg = pmm->getPDG();
                    int status = pmm->getSimulatorStatus();
                    p.SetParameter(pdg);
                    p.SetParameter((double)status);
                    //p.SetParameter( (double)(mcp->getGeneratorStatus()));
                }
                if (layer == 1 or layer == 31)
                {
                    simhits.push_back(p);
                }
                //cout << "simhits.size =" << simhits.size() << endl;
                // correct by SF (1st, 2nd, 3rd layer are 1., 0.0184, 0.0092)
                double ECAL_SF = 0.0184;
                e = e / ECAL_SF;
                px = e * x / _tmp;
                py = e * y / _tmp;
                pz = e * z / _tmp;
                PseudoJet pj_sf(px, py, pz, e);
                avec_hits_raw_sf.push_back(pj_sf); //???
            }
            vector<PseudoJet> avec_hits;
            vector<LParticle> Calhits;
            // ECAL hits after created by slicPandora
            double ecalsum = 0;
            IMPL::LCCollectionVec *col52 = (IMPL::LCCollectionVec *)evt->getCollection("EM_BARREL");
            nCL = col52->getNumberOfElements();
            for (int i = 0; i < nCL; ++i)
            {
                EVENT::CalorimeterHit *mcp = (EVENT::CalorimeterHit *)col52->getElementAt(i);
                const float *pos = mcp->getPosition();
                float x = pos[0];
                float y = pos[1];
                float z = pos[2];
                double e = mcp->getEnergy();
                double _tmp = std::sqrt(x * x + y * y + z * z);
                double px = e * x / _tmp;
                double py = e * y / _tmp;
                double pz = e * z / _tmp;
                int cellId0 = mcp->getCellID0();
                int cellId1 = mcp->getCellID1();
                long long cellId = ((long long)cellId1) << 32 | cellId0;
                int layer = decoder_ecal->getFieldValue("layer", cellId);
                LParticle p(px, py, pz, e, layer);
                ecalsum = ecalsum + e;
                PseudoJet pj(px, py, pz, e);
                double eta_r = pj.pseudorapidity();
                double phi_r = pj.phi();
                double pt_r = pj.pt();
                if (pt_r == 0)
                {
                    cout << "pt_r =" << pt_r << endl;
                }
                p.SetParameter(x * 1000);
                p.SetParameter(y * 1000);
                p.SetParameter(z * 1000);
                p.SetParameter(layer);
                //cout << " e=" << e << " phi=" << pj.phi() << " eta=" << pj.eta() << endl;
                if (e > 0.5)
                {
                    avec_hits.push_back(pj);
                    Calhits.push_back(p);
                    if (avec_hits.size() == 0)
                    {
                        cout << "avec_hits size =" << avec_hits.size() << endl;
                    }
                }
            }
            // ----------------- cluster jets --------------------------
            fastjet::ClusterSequenceArea *thisClustering_reco = new fastjet::ClusterSequenceArea(avec_hits, jet_def, fjAreaDefinition);
            vector<fastjet::PseudoJet> sjets_reco = sorted_by_pt(thisClustering_reco->inclusive_jets(25.0));
            vector<TLorentzVector> Recojets;
            //cout << "sjets_reco.size() =" << sjets_reco.size() << endl;
            //================================================================
            //       Save Recojet
            //================================================================
            for (unsigned int k = 0; k < sjets_reco.size(); k++)
            {
                TLorentzVector p_using_reco;
                p_using_reco.SetPxPyPzE(sjets_reco[k].px(), sjets_reco[k].py(), sjets_reco[k].pz(), sjets_reco[k].e());
                for (unsigned int iii = 0; iii < Truthjets_axis.size(); iii++)
                {
                    if (p_using_reco.DeltaR(Truthjets_axis[iii]) < 0.4)
                    {
                        //cout << "dr =" << p_using_reco.DeltaR(Truthjets_axis[iii]) << endl;
                        Recojets.push_back(p_using_reco);
                    }
                }
            }
            if (Recojets.size() == 0)
                continue;
            vector<vector<TLorentzVector>> FourP_dR_Reco(Recojets.size(), vector<TLorentzVector>());
            vector<vector<int>> PDG_Reco(Recojets.size(), vector<int>());
            vector<vector<double>> PT_Reco_sort(Recojets.size(), vector<double>());
            vector<vector<double>> PT_Reco(Recojets.size(), vector<double>());
            vector<vector<double>> PT_sort_number_only(Recojets.size(), vector<double>());
            vector<vector<double>> T_Reco_sort(Recojets.size(), vector<double>());
            vector<vector<double>> T_Reco(Recojets.size(), vector<double>());
            vector<vector<double>> T_sort_number_only(Recojets.size(), vector<double>());
            vector<vector<double>> T_first_last_average(2, vector<double>());

            vector<double> mass_Reco = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
            vector<int> event_number_Reco = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
            vector<int> Eta_smaller_than_1_event = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
            vector<int> event_number_Reco_for_mass = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
            vector<int> Eta_smaller_than_1_event_for_mass = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
            //cout << "Recojets.size() =" << Recojets.size() << endl;
            for (unsigned int Back_forth = 0; Back_forth < Recojets.size(); Back_forth++)
            {
                for (unsigned j1 = 0; j1 < simhits.size(); j1++)
                {
                    LParticle hit = (LParticle)simhits.at(j1);
                    TLorentzVector LE = hit.GetP();
                    double e_h = LE.E();
                    double phi_h = LE.Phi();
                    double eta_h = LE.PseudoRapidity();
                    vector<double> par = hit.GetParameters();
                    int type = hit.GetType(); // ECAL or HCAL
                    int layer = hit.GetCharge();
                    int x = int(par[0]);
                    int y = int(par[1]);
                    int z = int(par[2]);
                    double layer_pos_cm = par[3];
                    double Thit = par[4];
                    double LogTime = TMath::Log10(Thit);
                    double avt = par[5];
                    int pdg = int(par[6]);
                    int stat = int(par[7]);
                    double Mass = par[8];
                    //cout << "Thit =" << Thit << endl;
                    // Study of the dR
                    //Hadronic decays checking
                    for (unsigned int k1 = 0; k1 < Calhits.size(); k1++) //?????
                    {
                        LParticle p = (LParticle)Calhits.at(k1);
                        //PseudoJet p = (PseudoJet)avec_hits.at(k1);
                        TLorentzVector LE_1 = p.GetP();
                        vector<double> par1 = p.GetParameters();
                        int CellX_Cal = int(par1[0]);
                        int CellY_Cal = int(par1[1]);
                        int CellZ_Cal = int(par1[2]);
                        int LAYER_ECAL_USE = int(par1[3]);
                        //   cout << "Back_forth: " << Back_forth  << endl;
                        if (layer == 1 and (Recojets[Back_forth].DeltaR(LE_1)) < 0.4 and x == CellX_Cal and y == CellY_Cal and z == CellZ_Cal)
                        {
                            T_first_last_average[0].push_back(Thit);
                        }
                        if (layer == 31 and (Recojets[Back_forth].DeltaR(LE_1)) < 0.4 and x == CellX_Cal and y == CellY_Cal and z == CellZ_Cal)
                        {
                            T_first_last_average[1].push_back(Thit);
                            //Mass
                            mass_Reco[Back_forth] = mass_Reco[Back_forth] + Mass;
                            //Timing
                            //cout << "Timing: " << Thit << endl;
                            T_Reco_sort[Back_forth].push_back(Thit);
                            T_Reco[Back_forth].push_back(Thit);
                            //Four_momentum
                            FourP_dR_Reco[Back_forth].push_back(LE);

                            //cout << "PT_sort_number_only2 size = " << FourP_dR_Reco[Back_forth].Pt() << endl;
                            //for (auto x : FourP_dR_Reco[Back_forth].Pt())
                            //{
                            //    cout << x << " ";
                            //}
                            //cout << endl;
                            //PT
                            PT_Reco_sort[Back_forth].push_back(LE.Perp());
                            PT_Reco[Back_forth].push_back(LE.Perp());
                            //PDG
                            PDG_Reco[Back_forth].push_back(abs(pdg));
                            //Eta_selection
                            event_number_Reco[Back_forth] = event_number_Reco[Back_forth] + 1;
                            if (abs(eta_h) < 1)
                            {
                                Eta_smaller_than_1_event[Back_forth] = Eta_smaller_than_1_event[Back_forth] + 1;
                            }
                        }
                    }
                }
            }
            //================================================================
            //       Sort from small to big (T_Reco_sort and PT_Reco_sort)
            //================================================================
            for (unsigned int Back_forth = 0; Back_forth < Recojets.size(); Back_forth++)
            {
                sort(T_Reco_sort[Back_forth].begin(), T_Reco_sort[Back_forth].end());
                sort(PT_Reco_sort[Back_forth].begin(), PT_Reco_sort[Back_forth].end());
            }
            //================================================================
            // Finally output PT_sort_number_only == PT_Reco_sort size
            //================================================================
            vector<vector<TLorentzVector>> Highest_PT_FourP(Recojets.size(), vector<TLorentzVector>());
            for (unsigned int Back_forth = 0; Back_forth < Recojets.size(); Back_forth++)
            {
                //cout << "Back_forth =" << Back_forth << endl;
                if (T_Reco_sort[Back_forth].size() == 0)
                    continue;
                PT_sort_number_only[Back_forth].push_back(PT_Reco_sort[Back_forth][0]);
                T_sort_number_only[Back_forth].push_back(T_Reco_sort[Back_forth][0]);
                //cout << "PT_sort_number_only size = " << PT_sort_number_only[Back_forth].size() << endl;
                for (unsigned int uuu = 0; uuu < T_Reco_sort[Back_forth].size(); uuu++)
                {
                    // << "uuu =" << uuu << endl;
                    if (PT_Reco_sort[Back_forth][uuu] != PT_sort_number_only[Back_forth].back())
                    {
                        //cout << "test = " << PT_Reco_sort[Back_forth][uuu] << endl;
                        PT_sort_number_only[Back_forth].push_back(PT_Reco_sort[Back_forth][uuu]); //只存不同的pt
                    }
                    if (T_Reco_sort[Back_forth][uuu] != T_sort_number_only[Back_forth].back())
                    {
                        T_sort_number_only[Back_forth].push_back(T_Reco_sort[Back_forth][uuu]);
                    }
                }
                //cout << "PT_sort_number_only2 size = " << PT_sort_number_only[Back_forth].size() << endl;
                //for (auto x : PT_sort_number_only[Back_forth])
                //{
                //    cout << x << " ";
                //}
                //cout << endl;
                //cout << "PT_Reco_sort size = " << PT_Reco_sort[Back_forth].size() << endl;
                //for (auto x : PT_Reco_sort[Back_forth])
                //{
                //    cout << x << " ";
                //}
                //cout << endl;
            }
            //
            vector<int> Full_contain;
            unsigned int check_point_eta = 0;
            //=================================
            // make sure that 90% eta are small than 1
            //=================================
            for (unsigned int Back_forth = 0; Back_forth < Recojets.size(); Back_forth++)
            {
                //cout << "Back_forth =" << Back_forth << endl;
                //for (auto x : event_number_Reco)
                //{
                //    cout << x << " ";
                //}
                //cout << endl;
                if (event_number_Reco[Back_forth] > 0)
                {
                    //cout << "Fraction_of_the_event_forth====> " << Eta_smaller_than_1_event[Back_forth]/event_number_Reco[Back_forth] << endl;
                    if ((Eta_smaller_than_1_event[Back_forth] / event_number_Reco[Back_forth]) > 0.9) //???
                    {
                        Full_contain.push_back(1);
                    }
                    else
                    {
                        Full_contain.push_back(0);
                        check_point_eta = check_point_eta + 1;
                    }
                }
                else
                {
                    Full_contain.push_back(0);
                    check_point_eta = check_point_eta + 1;
                }
            }
            if (check_point_eta == Recojets.size())
            {
                continue;
            }
            vector<float> dR_Highest_PT_T_Reco;
            vector<float> dR_Highest_PT_PT_Reco;
            //================================================================
            //  Save the Time_Difference to TH1D
            //================================================================
            vector<vector<int>> PT_PDG_Reco(Recojets.size(), vector<int>());
            vector<vector<int>> T_PDG_Reco(Recojets.size(), vector<int>());

            double T_first = 0;
            double T_last = 0;
            if (T_first_last_average[0].size() != 0 and T_first_last_average[1].size() != 0)
            {
                float Time_Difference = 0;
                for (unsigned int G = 0; G < T_first_last_average[0].size(); G++)
                {
                    T_first = T_first + T_first_last_average[0][G];
                }
                for (unsigned int H = 0; H < T_first_last_average[1].size(); H++)
                {
                    T_last = T_last + T_first_last_average[1][H];
                }
                Time_Difference = (T_last / T_first_last_average[1].size()) - (T_first / T_first_last_average[0].size()); //???
                //cout << "Time_Difference: " << Time_Difference << endl;
                if (Time_Difference > 0)
                {
                    Timing_detecto_ECAL_TDif->Fill(Time_Difference);
                }
            }

            //================================================================
            // Save FourP with Pt is highest to Highest_PT_FourP
            //================================================================

            for (unsigned int Back_forth = 0; Back_forth < Recojets.size(); Back_forth++)
            //for (unsigned int Back_forth = 0; Back_forth <= 0; Back_forth++)
            {
                if (Full_contain[Back_forth] == 1)
                {
                    if (T_Reco_sort[Back_forth].size() > 0)
                    {
                        unsigned int Size_T_PT = T_Reco_sort[Back_forth].size();
                        for (unsigned int jkl = 0; jkl < Size_T_PT; jkl++)
                        {
                            if (PT_Reco_sort[Back_forth][Size_T_PT - 1] == PT_Reco[Back_forth][jkl])
                            {
                                //cout << "PT_Reco_sort size = " << PT_Reco_sort[Back_forth].size() << endl;
                                //cout << "Back_forth =" << Back_forth << " PT =" << PT_Reco_sort[Back_forth][Size_T_PT - 1] << endl;
                                Highest_PT_FourP[Back_forth].push_back(FourP_dR_Reco[Back_forth][jkl]);
                                //cout << "Back_forth =" << Back_forth << " PT =" << FourP_dR_Reco[Back_forth][jkl].Pt() << endl;
                            }
                            // cout << "PT_Reco_sort size = " << Highest_PT_FourP[Back_forth][jkl].Pt() << endl;
                        }

                        //cout << "End of event " << endl;

                        //================================================================
                        // Save T_sort_number_only the order is from largest to small (5 elements)
                        // dR_Highest_PT_T_Reco dR(Highest_PT_FourP, FourP)
                        //================================================================
                        //cout << "T_Reco = " << endl;
                        //for (auto x : T_Reco[0])
                        //{
                        //    cout << x << " ";
                        //}
                        //cout << endl;
                        for (unsigned int jkl = 0; jkl < 5; jkl++)
                        {
                            if (T_sort_number_only[Back_forth].size() > 0 and PT_sort_number_only[Back_forth].size() > 0)
                            {
                                //cout << "T_sort_number_only = " << T_sort_number_only[0].size() << endl;
                                if (jkl < T_sort_number_only[Back_forth].size())
                                {
                                    int identify_1 = 0;
                                    //cout << "Back_forth = " << Back_forth << endl;
                                    //cout << "T_sort_number_only[Back_forth][T_sort_number_only[Back_forth].size() - 1 - jkl] =" << T_sort_number_only[Back_forth][T_sort_number_only[Back_forth].size() - 1 - jkl] << endl;
                                    for (unsigned int ijk = 0; ijk < Size_T_PT; ijk++)
                                    {
                                        //cout << "start check" << endl;
                                        if (T_sort_number_only[Back_forth][T_sort_number_only[Back_forth].size() - 1 - jkl] == T_Reco[Back_forth][ijk])
                                        {
                                            //cout << "T_sort_number_only size = " << T_sort_number_only[Back_forth].size() << endl;
                                            //for (auto x : T_sort_number_only[Back_forth])
                                            //{
                                            //    cout << x << " ";
                                            //}
                                            //cout << endl;
                                            //cout << "start check" << endl;

                                            if (ijk == 0)
                                            {
                                                //cout << "T_Reco[Back_forth][ijk] =" << T_Reco[Back_forth][ijk] << endl;
                                                Timing_detector_Reco_TOF->Fill(T_Reco[Back_forth][ijk]);
                                            }

                                            //cout << "T_Reco[Back_forth][ijk]: " << T_Reco[Back_forth][ijk] << endl;
                                            //cout << "jkl: " << jkl << endl;
                                            //cout << "ijk: " << ijk << endl;
                                            identify_1++;
                                            T_PDG_Reco[Back_forth].push_back(PDG_Reco[Back_forth][ijk]);
                                            //cout << "T_PDG_Reco[Back_forth][ijk]: " << PDG_Reco[Back_forth][ijk] << endl;
                                            dR_Highest_PT_T_Reco.push_back(FourP_dR_Reco[Back_forth][ijk].DeltaR(Highest_PT_FourP[Back_forth][0]));
                                            //
                                            //cout << "FourP_dR_Reco[Back_forth][ijk]: " << FourP_dR_Reco[Back_forth][ijk].Pt() << endl;
                                            //cout << "Highest_PT_FourP[Back_forth]: " << Highest_PT_FourP[Back_forth] << endl;
                                            h_Particles_dR_Highest_PT_T_Reco[jkl]->Fill(FourP_dR_Reco[Back_forth][ijk].DeltaR(Highest_PT_FourP[Back_forth][0]));
                                            //cout << "TTT:FourP_dR_Reco[A  -+
                                            [ijk].DeltaR(Highest_PT_FourP[Back_forth][0]) : " << FourP_dR_Reco[Back_forth][ijk].DeltaR(Highest_PT_FourP[Back_forth][0]) << endl;
                                        }
                                    }
                                    //cout << "FourP_dR_Reco" << endl;
                                    //for (TLorentzVector x : FourP_dR_Reco[0])
                                    //{
                                    //    //cout << x.Print() << " ";
                                    //    x.Print();
                                    //}

                                    if (identify_1 > 1)
                                    {
                                        cout << "Weird! Check!" << endl;
                                    }
                                    //================================================================
                                    // Save PT_sort_number_only the order is from small to larger (5 elements)
                                    // dR_Highest_PT_PT_Reco dR(Highest_PT_FourP, FourP)
                                    //================================================================
                                    if (jkl < PT_sort_number_only[Back_forth].size())
                                    {
                                        //cout << "PT_sort_number_only[Back_forth].size() =" << PT_sort_number_only[Back_forth].size() << endl;
                                        int identify = 0;
                                        for (unsigned int ijk = 0; ijk < Size_T_PT; ijk++)
                                        {
                                            if (PT_sort_number_only[Back_forth][jkl] == PT_Reco[Back_forth][ijk])
                                            {
                                                //cout << "PT_sort_number_only[Back_forth][jkl]: " << PT_sort_number_only[Back_forth][jkl] << "PT_Reco[Back_forth][ijk]: " << PT_Reco[Back_forth][ijk] << endl;
                                                //cout << "ijk: " << ijk << endl;
                                                identify = identify + 1;
                                                PT_PDG_Reco[Back_forth].push_back(PDG_Reco[Back_forth][ijk]);
                                                //cout << "PT_PDG_Reco[Back_forth][ijk]: " << PDG_Reco[Back_forth][ijk] << endl;
                                                dR_Highest_PT_PT_Reco.push_back(FourP_dR_Reco[Back_forth][ijk].DeltaR(Highest_PT_FourP[Back_forth][0]));
                                                h_Particles_dR_Highest_PT_PT_Reco[jkl]->Fill(FourP_dR_Reco[Back_forth][ijk].DeltaR(Highest_PT_FourP[Back_forth][0]));
                                                //cout << "PTPTPT:FourP_dR_Reco[Back_forth][ijk].DeltaR(Highest_PT_FourP[Back_forth][0]): " << FourP_dR_Reco[Back_forth][ijk].DeltaR(Highest_PT_FourP[Back_forth][0]) << endl;
                                            }
                                        }
                                        /*
                                        if (identify > 1)
                                        {
                                            //cout << "Weird! Check_1!" << endl;
                                        }
                                        */
                                    }
                                    //cout << "identify = " << identify_1 << endl;
                                }
                            }
                        }
                    }
                }
            }
            //cout << "FourP_dR_Reco" << endl;
            //for (TLorentzVector x : FourP_dR_Reco[0])
            //{
            //    cout << x.Pt() << " ";
            //}
            //cout << endl;
            //cout << "Highest_PT_FourP" << endl;
            //for (TLorentzVector x : Highest_PT_FourP[0])
            //{
            //    cout << x.Pt() << " ";
            //}
            //cout << endl;
            for (unsigned int Back_forth = 0; Back_forth < Recojets.size(); Back_forth++)
            {
                for (unsigned int j = 0; j < 5; j++)
                {
                    if (T_PDG_Reco[Back_forth].size() > 0)
                    {
                        if ((T_PDG_Reco[Back_forth].size()) > j)
                        {
                            int it;
                            it = find(Trailing_particle_kind_T.begin(), Trailing_particle_kind_T.end(), abs(T_PDG_Reco[Back_forth][j]))[0];
                            //cout << "it =" << it << endl;
                            if (it != abs(T_PDG_Reco[Back_forth][j]))
                            {
                                //       cout << "============================================================================================================================================================" << endl;
                                //cout << "abs(T_PDG_Reco[Back_forth][j]): " << abs(T_PDG_Reco[Back_forth][j]) << endl;
                                Trailing_particle_kind_T.push_back(abs(T_PDG_Reco[Back_forth][j]));
                                for (unsigned m = 0; m < Trailing_particle_kind_T.size(); m++)
                                {
                                    //cout << "Trailing_particle_kind_T[m]: "<< abs(Trailing_particle_kind_T[m]) << endl;
                                    if (abs(T_PDG_Reco[Back_forth][j]) == Trailing_particle_kind_T[m] and abs(T_PDG_Reco[Back_forth][j]) != 22)
                                    {
                                        h_Particles_Rank_T_Reco[j]->Fill(m);
                                    }
                                }
                            }
                            else
                            {
                                for (unsigned int m = 0; m < Trailing_particle_kind_T.size(); m++)
                                {
                                    //cout << "Trailing_particle_kind_T[m]: "<< abs(Trailing_particle_kind_T[m]) << endl;
                                    if (abs(T_PDG_Reco[Back_forth][j]) == Trailing_particle_kind_T[m] and abs(T_PDG_Reco[Back_forth][j]) != 22)
                                    {
                                        h_Particles_Rank_T_Reco[j]->Fill(m);
                                    }
                                    //cout << "Particle_ID_T_Reco[j]->GetBinContent(m): " << h_Particles_Rank_T_Reco[j]->GetBinContent(m + 1) << endl;
                                }
                            }
                        }
                    }
                    if (PT_PDG_Reco[Back_forth].size() > 0)
                    {
                        if ((PT_PDG_Reco[Back_forth].size()) > j)
                        {
                            //cout << "abs(PT_PDG_Reco[Back_forth][j]): " << abs(PT_PDG_Reco[Back_forth][j]) << endl;
                            int it2;
                            it2 = find(Trailing_particle_kind_PT.begin(), Trailing_particle_kind_PT.end(), abs(PT_PDG_Reco[Back_forth][j]))[0];
                            if (it2 != abs(PT_PDG_Reco[Back_forth][j]))
                            {
                                Trailing_particle_kind_PT.push_back(abs(PT_PDG_Reco[Back_forth][j]));
                                for (unsigned int m = 0; m < Trailing_particle_kind_PT.size(); m++)
                                {
                                    if (abs(PT_PDG_Reco[Back_forth][j]) == Trailing_particle_kind_PT[m] and abs(PT_PDG_Reco[Back_forth][j]) != 22)
                                    {
                                        h_Particles_Rank_PT_Reco[j]->Fill(m);
                                    }
                                }
                            }
                            else
                            {
                                for (unsigned int m = 0; m < Trailing_particle_kind_PT.size(); m++)
                                {
                                    //cout << "Trailing_particle_kind_PT[m]: "<< abs(Trailing_particle_kind_PT[m]) << endl;
                                    if (abs(PT_PDG_Reco[Back_forth][j]) == Trailing_particle_kind_PT[m] and abs(PT_PDG_Reco[Back_forth][j]) != 22)
                                    {
                                        h_Particles_Rank_PT_Reco[j]->Fill(m);
                                    }
                                    //cout << "Particle_ID_PT_Reco[j]->GetBinContent(m): " << h_Particles_Rank_PT_Reco[j]->GetBinContent(m+1) << endl;
                                }
                            }
                        }
                    }
                }
            }
            // cout << "99 : ";
            if (dR_Highest_PT_PT_Reco.size() >= 1)
            {
                dR_Tr0PT_HPt_Reco = dR_Highest_PT_PT_Reco[0];
            }
            if (dR_Highest_PT_PT_Reco.size() >= 2)
            {
                dR_Tr1PT_HPt_Reco = dR_Highest_PT_PT_Reco[1];
            }
            if (dR_Highest_PT_PT_Reco.size() >= 3)
            {
                dR_Tr2PT_HPt_Reco = dR_Highest_PT_PT_Reco[2];
            }
            if (dR_Highest_PT_PT_Reco.size() >= 4)
            {
                dR_Tr3PT_HPt_Reco = dR_Highest_PT_PT_Reco[3];
            }
            if (dR_Highest_PT_PT_Reco.size() >= 5)
            {
                dR_Tr4PT_HPt_Reco = dR_Highest_PT_PT_Reco[4];
            }
            if (dR_Highest_PT_T_Reco.size() >= 1)
            {
                dR_Tr0T_HPt_Reco = dR_Highest_PT_T_Reco[0];
            }
            if (dR_Highest_PT_T_Reco.size() >= 2)
            {
                dR_Tr1T_HPt_Reco = dR_Highest_PT_T_Reco[1];
            }
            if (dR_Highest_PT_T_Reco.size() >= 3)
            {
                dR_Tr2T_HPt_Reco = dR_Highest_PT_T_Reco[2];
            }
            if (dR_Highest_PT_T_Reco.size() >= 4)
            {
                dR_Tr3T_HPt_Reco = dR_Highest_PT_T_Reco[3];
            }
            if (dR_Highest_PT_T_Reco.size() >= 5)
            {
                dR_Tr4T_HPt_Reco = dR_Highest_PT_T_Reco[4];
            }
            for (unsigned int Back_forth = 0; Back_forth < Recojets.size(); Back_forth++)
            {
                if (event_number_Reco[Back_forth] > 0)
                {
                    // cout << "event_number_Reco[Back_forth]: " << event_number_Reco[Back_forth] << endl;
                    // cout << "Fraction_of_the_event_forth_for_mass====> " << Eta_smaller_than_1_event[Back_forth]/event_number_Reco[Back_forth] << endl;
                    if ((Eta_smaller_than_1_event[Back_forth] / event_number_Reco[Back_forth]) > 0.9)
                    {
                        mass_sum_average_Reco->Fill(mass_Reco[Back_forth] / event_number_Reco[Back_forth]);
                        //cout << "mass_sum_average_Reco_forth_for_mass====> " << mass_Reco[Back_forth]/event_number_Reco[Back_forth] << endl;
                    }
                }
            }
            //Clear vector
            Full_contain.clear();
            avec_hits.clear();
            avec_clus.clear();
            Calhits.clear();
            avec_hits_raw_sf.clear();
            avec_hits_raw.clear();

            double _bField = 5;
            // Reconstructed tracks
            // look up at: /share/sl6/ilcsoft/slic/release-v05-00-00/slicPandora/HEAD/lcio/v02-04-03/build/include/EVENT
            vector<PseudoJet> avec_tracks;
            vector<LParticle> trackhits;
            IMPL::LCCollectionVec *colT = (IMPL::LCCollectionVec *)evt->getCollection("Tracks");
            unsigned int nTracks = colT->getNumberOfElements();
            for (unsigned int j = 0; j < Forth_And_Back_Vector.size(); j++)
            {
                TLorentzVector Jet_axis_Truth = Forth_And_Back_Vector[j];
                for (unsigned int i = 0; i < nTracks; i++)
                {
                    EVENT::Track *track = (EVENT::Track *)colT->getElementAt(i);
                    float d0 = track->getD0();
                    float z0 = track->getZ0();
                    float omega = track->getOmega();
                    float tanLambda = track->getTanLambda();
                    float phi0 = track->getPhi();
                    float radius = 1.0 / omega;
                    float x0 = radius * cos(phi0 - k2PI);
                    float y0 = radius * sin(phi0 - k2PI);
                    const pandora::Helix helixFit(phi0, d0, z0, omega, tanLambda, _bField);
                    const float recoMomentum(helixFit.GetMomentum().GetMagnitude());
                    double px = helixFit.GetMomentum().GetX();
                    double py = helixFit.GetMomentum().GetY();
                    double pz = helixFit.GetMomentum().GetZ();
                    double m = 0; // mcp->getMomentum()[3];
                    double e = sqrt(px * px + py * py + pz * pz + m * m);
                    PseudoJet pp(px, py, pz, e);
                    double eta_r = pp.pseudorapidity();
                    double phi_r = pp.phi();
                    double pt_r = pp.pt();
                    TLorentzVector TLV;
                    TLV.SetPxPyPzE(px, py, pz, e);
                    // Get the two 32-bit chunks of the ID.
                    //===comment

                    for (int ppp = 0; ppp < track->getTrackerHits().size(); ppp++)
                    {
                        int cellId0 = track->getTrackerHits()[ppp]->getCellID0();
                        int cellId1 = track->getTrackerHits()[ppp]->getCellID1();
                        cout << "cellId0: " << cellId0 << endl;
                        cout << "cellId1: " << cellId1 << endl;
                        // Make a 64-bit id for the IDDecoder.  The type MUST be "long long" and not "long".  (from Tony Johnson)
                        long long cellId = ((long long)cellId1) << 32 | cellId0;
                        int layer = decoder->getFieldValue("layer", cellId);
                        cout << "Tracker_Layer: " << layer << endl;
                        continue;
                    }
                    if (TLV.DeltaR(Jet_axis_Truth) < 0.4 and Check_Forth_And_Back_Bool[j] == 1)
                    {
                        trackhits.push_back(p);
                    }
                } //end matching
            }

            T->Fill();
            T_Reco_T->Fill();
            T_Reco_T_track->Fill();
            //cout << "Final: ";
        }
        cout << "nnEvents = " << nnEvents << endl;
        lcReader->close();
        delete lcReader;
    } //End loop over all files
    RootFile->Write();
    mytree->Write();
    h_pt_clus->Write();
    h_eta_clus->Write();
    h_pt_ecal->Write();
    Timing_detecto_ECAL_TDif->Write();
    Timing_detector_Reco_TOF->Write();
    Timing_detector_Reco_TOF_track->Write();
    RootFile->Close();

    return 0;
}