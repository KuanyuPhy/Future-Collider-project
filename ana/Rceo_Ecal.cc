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
    bool debug = false;
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
    TH1D *h_jet_hitstime1D_ECAL = new TH1D("jet_hitstime_Ecal_1D", "1D for time vs hit in HCAL", 100, 0.5, 10);
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
            }
            if (debug)
            {
                cout << "1" << endl;
            }
            // get Reco from EMhit
            /*
            double ecalsum_raw = 0;
            vector<LParticle> simhits;
            vector<PseudoJet> avec_hits_raw_sf;
            IMPL::LCCollectionVec *col50_2 = (IMPL::LCCollectionVec *)evt->getCollection("EM_BARREL");
            int nCL12 = col50_2->getNumberOfElements();
            for (int i = 0; i < nCL12; ++i)
            {
                if (debug)
                {
                    cout << "1.1" << endl;
                }
                EVENT::SimCalorimeterHit *mcp_1 = (EVENT::SimCalorimeterHit *)col50_2->getElementAt(i);
                const float *pos = mcp_1->getPosition();
                float x = pos[0];
                float y = pos[1];
                float z = pos[2];
                double e = mcp_1->getEnergy();
                double _tmp = std::sqrt(x * x + y * y + z * z);
                double px = e * x / _tmp;
                double py = e * y / _tmp;
                double pz = e * z / _tmp;
                // Get the two 32-bit chunks of the ID.
                int cellId0 = mcp_1->getCellID0();
                int cellId1 = mcp_1->getCellID1();
                if (debug)
                {
                    cout << "1.2" << endl;
                }
                // Make a 64-bit id for the IDDecoder.  The type MUST be "long long" and not "long".  (from Tony Johnson)
                long long cellId = ((long long)cellId1) << 32 | cellId0;
                if (debug)
                {
                    cout << "1.3" << endl;
                }
                int layer = decoder_ecal->getFieldValue("layer", cellId);
                // 1st layer on middle
                if (debug)
                {
                    cout << "1.4" << endl;
                }
                double layer_pos_cm = 0.1 * (layer_size_ecal_mm * 0.5 + (layer * layer_size_ecal_mm));
                if (debug)
                {
                    cout << "1.5" << endl;
                }
                //double Thit = mcp_1->getTimeCont(0);
                if (debug)
                {
                    cout << "1.6" << endl;
                }
                // number of MC contributions to the hit
                // calculate average time  in [ns]
                double avt = 0;
                double ave = 0;
                if (debug)
                {
                    cout << "2" << endl;
                }
                for (int jj = 0; jj < mcp_1->getNMCContributions(); jj++)
                {
                    avt = avt + mcp_1->getEnergyCont(jj) * mcp_1->getTimeCont(jj);
                    ave = ave + mcp_1->getEnergyCont(jj);
                }
                avt = avt / ave;
                ecalsum_raw = ecalsum_raw + e;
                PseudoJet pj(px, py, pz, e);
                double eta_r = pj.pseudorapidity();
                double phi_r = pj.phi();
                double pt_r = pj.pt();
                // fill hits
                LParticle p(px, py, pz, e, layer);
                p.SetCharge(layer);
                p.SetType(2); // ECAL
                p.SetStatus(mcp_1->getNMCContributions());
                p.SetParameter(x * 1000);
                p.SetParameter(y * 1000);
                p.SetParameter(z * 1000);
                p.SetParameter(layer_pos_cm);
                // find fastest hit
                float timeCont = mcp_1->getTimeCont(0);
                EVENT::MCParticle *pmm = mcp_1->getParticleCont(0);
                int pdg = pmm->getPDG();
                int status = pmm->getSimulatorStatus();
                float rawTime = timeCont;
                float FAM = pmm->getMass();
                int nCont = mcp_1->getNMCContributions();
                for (int jjj = 0; jjj < nCont; jjj++)
                {
                    if (mcp_1->getTimeCont(jjj) < rawTime)
                    {
                        rawTime = mcp_1->getTimeCont(jjj);
                        EVENT::MCParticle *pmm = mcp_1->getParticleCont(jjj);
                        pdg = pmm->getPDG();
                        status = pmm->getSimulatorStatus();
                    }
                    p.SetParameter(rawTime); // fastest hit
                    p.SetParameter(avt);     // average hit time
                    p.SetParameter((double)pdg);
                    p.SetParameter((double)status);
                    p.SetParameter((double)FAM);
                    avec_hits_raw_sf.push_back(pj); //need to check
                }
                if (debug)
                {
                    cout << "3" << endl;
                }
                if (layer == 1 or layer == 31)
                {
                    simhits.push_back(p);
                    // correct by SF (1st, 2nd, 3rd layer are 1., 0.0184, 0.0092)
                    double ECAL_SF = 0.0184;
                    e = e / ECAL_SF;
                    px = e * x / _tmp;
                    py = e * y / _tmp;
                    pz = e * z / _tmp;
                    //PseudoJet pj_sf(px, py, pz, e);
                    //avec_hits_raw_sf.push_back(pj_sf);
                }
                if (Truthjets_axis.size() > 0)
                {
                    int activeAreaRepeats_1 = 1;
                    double ghostArea_1 = 0.01;
                    double ghostEtaMax_1 = 7.0;
                    fastjet::GhostedAreaSpec fjActiveArea_1(ghostEtaMax_1, activeAreaRepeats_1, ghostArea_1);
                    fastjet::AreaDefinition fjAreaDefinition_1(fastjet::active_area, fjActiveArea_1);
                    fastjet::ClusterSequenceArea *thisClustering_reco = new fastjet::ClusterSequenceArea(avec_hits_raw_sf, jet_def, fjAreaDefinition_1);
                    vector<fastjet::PseudoJet> sjets_reco = sorted_by_pt(thisClustering_reco->inclusive_jets(25.0));
                    vector<TLorentzVector> Recojets;
                    for (unsigned int k = 0; k < sjets_reco.size(); k++)
                    {
                        TLorentzVector p_using_reco;
                        p_using_reco.SetPxPyPzE(sjets_reco[k].px(), sjets_reco[k].py(), sjets_reco[k].pz(), sjets_reco[k].e());
                        for (unsigned int iii = 0; iii < Truthjets_axis.size(); iii++)
                        {
                            if (p_using_reco.DeltaR(Truthjets_axis[iii]) < 0.4)
                            {
                                Recojets.push_back(p_using_reco);
                            }
                        }
                    }
                }

            }
            */
            /*
            double calo_sum = 0;
            // clusters
            vector<PseudoJet> avec_clus;
            double xsum_cl = 0;
            IMPL::LCCollectionVec *col5 = (IMPL::LCCollectionVec *)evt->getCollection("ReconClusters");
            int nCL = col5->getNumberOfElements();
            for (int i = 0; i < nCL; ++i)
            {
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
                if (pt_r > minPtConst)
                {
                    avec_clus.push_back(pj);
                }
            }*/
            //cout << "Calo clus sum=" << calo_sum << endl;
            /*
            vector<PseudoJet> avec_hits;
            double hcalsum = 0;
            IMPL::LCCollectionVec *col51 = (IMPL::LCCollectionVec *)evt->getCollection("HAD_BARREL");
            nCL = col51->getNumberOfElements();
            for (int i = 0; i < nCL; ++i)
            {
                EVENT::CalorimeterHit *mcp = (EVENT::CalorimeterHit *)col51->getElementAt(i);
                const float *pos = mcp->getPosition();
                float x = pos[0];
                float y = pos[1];
                float z = pos[2];
                double e = mcp->getEnergy();
                double _tmp = std::sqrt(x * x + y * y + z * z);
                double px = e * x / _tmp;
                double py = e * y / _tmp;
                double pz = e * z / _tmp;
                hcalsum = hcalsum + e;
                PseudoJet pj(px, py, pz, e);
                double eta_r = pj.pseudorapidity();
                double phi_r = pj.phi();
                double pt_r = pj.pt();
                //cout << " e=" << e <<  " phi=" << pj.phi() << " eta=" << pj.eta() << endl;
                avec_hits.push_back(pj);

                //cout << "HCal raw sum=" << hcalsum_raw << endl;
                //cout << "HCal sum after SF=" << hcalsum << endl;
            }*/
            vector<LParticle> simhits;
            vector<PseudoJet> avec_hits;
            vector<PseudoJet> Calhits;
            vector<PseudoJet> avec_hits_raw_sf;
            double ecalsum = 0;
            // ECAL hits
            double ecalsum_raw = 0;
            // ECAL hits
            IMPL::LCCollectionVec *col53 = (IMPL::LCCollectionVec *)evt->getCollection("EcalBarrelHits");
            IMPL::LCCollectionVec *col52 = (IMPL::LCCollectionVec *)evt->getCollection("EM_BARREL");
            int nCL = col52->getNumberOfElements();
            for (int i = 0; i < nCL; ++i)
            {
                //EVENT::SimCalorimeterHit *mcp = (EVENT::SimCalorimeterHit *)col52->getElementAt(i);
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
                float time = mcp->getTime();
                //cout << "time = " << time << endl;
                // fill hits

                //Calhits.push_back(p);
                // Get the two 32-bit chunks of the ID.
                int cellId0 = mcp->getCellID0();
                int cellId1 = mcp->getCellID1();
                // Make a 64-bit id for the IDDecoder.  The type MUST be "long long" and not "long".  (from Tony Johnson)
                long long cellId = ((long long)cellId1) << 32 | cellId0;
                int layer = decoder_ecal->getFieldValue("layer", cellId);
                // 1st layer on middle
                double layer_pos_cm = 0.1 * (layer_size_ecal_mm * 0.5 + (layer * layer_size_ecal_mm));
                //double Thit = mcp->getTimeCont();
                double avt = 0;
                double ave = 0;

                LParticle p(px, py, pz, e, layer);
                PseudoJet pj_sf(px, py, pz, e);
                p.SetParameter(x * 1000);
                p.SetParameter(y * 1000);
                p.SetParameter(z * 1000);
                p.SetParameter(layer);
                if (e > 0.5)
                {
                    if (layer == 1 or layer == 31)
                    {
                        avec_hits_raw_sf.push_back(pj_sf);
                    }
                }
            }
            if (Truthjets_axis.size() > 0)
            {
                fastjet::ClusterSequenceArea *thisClustering_reco = new fastjet::ClusterSequenceArea(avec_hits_raw_sf, jet_def, fjAreaDefinition);
                vector<fastjet::PseudoJet> sjets_reco = sorted_by_pt(thisClustering_reco->inclusive_jets(25.0));
                vector<TLorentzVector> Recojets;
                for (unsigned int k = 0; k < sjets_reco.size(); k++)
                {
                    TLorentzVector p_using_reco;
                    p_using_reco.SetPxPyPzE(sjets_reco[k].px(), sjets_reco[k].py(), sjets_reco[k].pz(), sjets_reco[k].e());
                    for (unsigned int iii = 0; iii < Truthjets_axis.size(); iii++)
                    {
                        if (p_using_reco.DeltaR(Truthjets_axis[iii]) < 0.4)
                        {
                            Recojets.push_back(p_using_reco);
                            if (Recojets.size() < 0)
                            {
                                continue;
                            }
                        }
                    }
                }
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
                for (int Back_forth = 0; Back_forth < Recojets.size(); Back_forth++)
                {
                    for (int j1 = 0; j1 < avec_hits_raw_sf.size(); j1++)
                    {
                        cout << "1" << endl;
                    }
                }
            }
            // ----------------- cluster jets --------------------------
            //ClusterSequence clust_seq_clus(avec_clus, jet_def);
        }
        cout << "nnEvents = " << nnEvents << endl;
        lcReader->close();
        delete lcReader;
    } //End loop over all files
    RootFile->Write();
    mytree->Write();
    h_pt_clus->Write();
    h_eta_clus->Write();
    RootFile->Close();
    return 0;
}