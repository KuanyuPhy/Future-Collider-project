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

    TH1D *h_jet_pt_truth_check = new TH1D("h_jet_pt_truth_check", "pT [GeV] ", 100, 0, 2700); // plus Geant4
    h_jet_pt_truth_check->GetXaxis()->SetTitle("P_{T}^{true,jet}");
    //h_jet_pt_truth_check->GetYaxis()->SetTitle("Entries");
    TH1D *h_jet_eta_truth_check = new TH1D("h_jet_eta_truth_check", "#eta", 100, -5, 5); //before eta cut
    TH1D *h_jet_n_truth = new TH1D("h_jet_n_truth", "Nr of truth jets", 10, 0, 10);      //before match
    TH1D *h_jet_nn_truth = new TH1D("h_jet_nn_truth", "Nr of truth jets", 10, 0, 10);    //after match
    TH1D *h_jet_m_truth = new TH1D("jet_m_truth", "Mass [GeV]", 100, 0, 100);

    TH1D *h_jet_pt_reco_check = new TH1D("h_jet_pt_reco_check", "pT [GeV] ", 100, 0, 2700); // plus Geant4
    h_jet_pt_reco_check->GetXaxis()->SetTitle("P_{T}^{Reco,jet}");
    //h_jet_pt_truth_check->GetYaxis()->SetTitle("Entries");
    TH1D *h_jet_eta_reco_check = new TH1D("h_jet_eta_reco_check", "#eta", 100, -5, 5); //before eta cut
    TH1D *h_jet_n_reco = new TH1D("h_jet_n_reco", "Nr of truth jets", 10, 0, 10);      //before match
    TH1D *h_jet_nn_reco = new TH1D("h_jet_nn_reco", "Nr of truth jets", 10, 0, 10);    //after match
    TH1D *h_jet_m_reco = new TH1D("h_jet_m_reco", "Mass [GeV]", 100, 0, 100);

    TH1D *h_pt_ecal = new TH1D("ecal_pt", "pt ecal", 270, 0, 2700);
    TH1D *Timing_detecto_ECAL_TDif = new TH1D("Timing_detecto_ECAL_TDif", "Timing_detecto_ECAL_TDif", 200, 0, 50);
    Timing_detecto_ECAL_TDif->GetXaxis()->SetTitle("T[ns]");

    TH1D *h_effR_j1 = new TH1D("h_effR_j1", "not match true jet", 100, 0., 10.0);
    TH1D *h_effR_j2 = new TH1D("h_effR_j1", "one match true jet", 100, 0., 10.0);
    TH1D *h_effR_j3 = new TH1D("h_effR_j1", "two match true jet", 100, 0., 10.0);

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

    Float_t dr0;
    Float_t dr1;
    TTree *tGEN = new TTree("tGEN", "Tree with vectors");
    tGEN->Branch("dr0", &dr0);
    tGEN->Branch("dr1", &dr1);

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
            //vector<fastjet::PseudoJet> Truthjets_axis;
            vector<TLorentzVector> Truthjets_axis;
            vector<int> Truthjets_axis_index;
            vector<double> deltaRRR;
            double deltaR00, deltaR01, deltaR02;
            bool isTrueMatch = false;
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
                h_jet_m_truth->Fill(m);
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
                        isTrueMatch = true;
                        //cout << "k =" << k << endl;
                    }
                    if (Truthjets_axis.size() == 0)
                    {
                        //deltaR00.push_back(dr);
                        deltaR00 = dr;
                        h_effR_j1->Fill(deltaR00);
                    }
                    if (Truthjets_axis.size() == 1)
                    {
                        //deltaR01.push_back(dr);
                        deltaR01 = dr;
                        h_effR_j2->Fill(deltaR01);
                    }
                    if (Truthjets_axis.size() == 2)
                    {
                        //deltaR02.push_back(dr);
                        deltaR02 = dr;
                        h_effR_j3->Fill(deltaR02);
                    }
                }

                //==============================================================
                // Selection: Each Jet constituents & abs(eta)<1
                //==============================================================
                if (TMath::Abs(eta) > 1)
                    continue;

            } //end of sjets_truth loop
            h_jet_n_truth->Fill(sjets_truth.size());
            h_jet_nn_truth->Fill(Truthjets_axis.size());

            // ECAL hits
            //IMPL::LCCollectionVec *col53 = (IMPL::LCCollectionVec *)evt->getCollection("EcalBarrelHits");
            //int nCL = col53->getNumberOfElements();
            //================================================================
            // Get Collection from EM_BARREL (after sampling function)
            //================================================================
            vector<PseudoJet> avec_hits;
            vector<LParticle> Calhits;
            double ecalsum = 0;
            IMPL::LCCollectionVec *col52 = (IMPL::LCCollectionVec *)evt->getCollection("EM_BARREL");
            int nCL = col52->getNumberOfElements();
            for (int i = 0; i < nCL; i++)
            {
                EVENT::CalorimeterHit *mcp = (EVENT::CalorimeterHit *)col52->getElementAt(i);
                const float *pos = mcp->getPosition();
                float x = pos[0];
                float y = pos[1];
                float z = pos[2];
                double e = mcp->getEnergy();
                //double e = sqrt(px * px + py * py + pz * pz + m * m);
                double _tmp = std::sqrt(x * x + y * y + z * z);
                double px = e * x / _tmp;
                double py = e * y / _tmp;
                double pz = e * z / _tmp;
                //double momentum = TMath::Power((TMath::Power(px, 2) + TMath::Power(py, 2) + TMath::Power(pz, 2)), 0.5);
                int cellId0 = mcp->getCellID0();
                int cellId1 = mcp->getCellID1();
                long long cellId = ((long long)cellId1) << 32 | cellId0;
                int layer = decoder_ecal->getFieldValue("layer", cellId);
                double Thit = mcp->getTime();
                //cout << "CalorimeterHit Thit =" << Thit << endl;
                //cout << "CalorimeterHit layer =" << layer << endl;
                ecalsum = ecalsum + e;
                PseudoJet pj(px, py, pz, e);
                double eta_r = pj.pseudorapidity();
                double phi_r = pj.phi();
                double pt_r = pj.pt();
                if (pt_r == 0)
                {
                    cout << "pt_r =" << pt_r << endl;
                }
                // fill hits
                LParticle p(px, py, pz, e, layer);
                p.SetCharge(layer);
                p.SetType(2);
                p.SetParameter(Thit); //0

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
            vector<fastjet::PseudoJet> Recojetstest;
            vector<int> Recojets_axis_index;
            bool isRecoMatch = false;
            //cout << "sjets_reco.size() =" << sjets_reco.size() << endl;
            //================================================================
            //       Save Recojet
            //================================================================
            for (unsigned int k = 0; k < sjets_reco.size(); k++)
            {
                double eta = sjets_reco[k].pseudorapidity();
                double phi = sjets_reco[k].phi();
                if (phi < 0)
                    phi = phi + k2PI;
                double m = sjets_reco[k].m();
                double pt = sjets_reco[k].perp();
                double e = sjets_reco[k].e();
                h_jet_pt_reco_check->Fill(pt);
                h_jet_eta_reco_check->Fill(eta);

                h_jet_m_reco->Fill(m);
                TLorentzVector p_using_reco;
                p_using_reco.SetPxPyPzE(sjets_reco[k].px(), sjets_reco[k].py(), sjets_reco[k].pz(), sjets_reco[k].e());
                //=====================================================
                // Matching:  Require Truth jet and ReoJet dr < 0.4
                //=====================================================
                for (unsigned int iii = 0; iii < Truthjets_axis.size(); iii++)
                {
                    if (p_using_reco.DeltaR(Truthjets_axis[iii]) < 0.4)
                    {
                        //cout << "dr =" << p_using_reco.DeltaR(Truthjets_axis[iii]) << endl;
                        Recojets.push_back(p_using_reco);
                        Recojets_axis_index.push_back(k); //Pass matching index
                        Recojetstest.push_back(p_using_reco);
                        isRecoMatch = true;
                    }
                    if (isRecoMatch)
                    {
                        break;
                    }
                }
                //if (Recojets.size() == 0)
                //{
                //    cout << "0" << endl;
                //}

            } //end of save Recojets loop
            h_jet_n_reco->Fill(sjets_reco.size());
            h_jet_nn_reco->Fill(Recojets.size());
            for (int i = 0; i < Recojets_axis_index.size(); i++)
            {
                int pass = Recojets_axis_index[i];
                //cout << "pass =" << pass << endl;
                vector<float> ECALhits0, ECALhits31;
                for (unsigned int j1 = 0; j1 < Calhits.size(); j1++)
                {
                    LParticle hit = (LParticle)Calhits.at(j1);
                    TLorentzVector LE = hit.GetP();
                    double phi_h = LE.Phi();
                    double eta_h = LE.PseudoRapidity();
                    vector<double> par = hit.GetParameters();
                    int type = hit.GetType(); // ECAL or HCAL
                    int layer = hit.GetCharge();
                    double Thit = par[0];
                    /*
                    if (layer == 1 && (sjets_reco[pass].DeltaR(LE)) < 0.4)
                    {
                        ECALhits0.push_back(Thit);
                    }
                    */
                    /*
                    if (layer == 31 && (sjets_reco[pass].DeltaR(LE)) < 0.4)
                    {
                        ECALhits31.push_back(Thit);
                        //Timing
                        T_Reco_sort[pass].push_back(Thit);
                        T_Reco[pass].push_back(Thit);
                        //Four_momentum
                        FourP_dR_Reco[pass].push_back(LE);
                        //PT
                        PT_Reco_sort[pass].push_back(LE.Perp());
                        PT_Reco[pass].push_back(LE.Perp());
                        //Velocity
                        vector<PseudoJet> constit = sjets_reco[pass].constituents();
                        int csize = constit.size();
                        for (int j = 0; j < csize; j++)
                        {
                            double px = constit[j].px();
                            double py = constit[j].px();
                            double pz = constit[j].pz();
                            double e = constit[j].e();
                            double constit_velocity_vt = sqrt(px * px + py * py + pz * pz) / e;
                            cout << "Vt = " << constit_velocity_vt << endl;
                        }
                        //Eta_selection
                        event_number_Reco[pass] = event_number_Reco[pass] + 1;
                        if (abs(eta_h) < 1)
                        {
                            Eta_smaller_than_1_event[pass] = Eta_smaller_than_1_event[pass] + 1;
                        }
                    }
                    */
                }
            }

            vector<vector<TLorentzVector>> FourP_dR_Reco(Recojets.size(), vector<TLorentzVector>());
            vector<vector<int>> PDG_Reco(Recojets.size(), vector<int>());
            vector<vector<double>> PT_Reco_sort(Recojets.size(), vector<double>());
            vector<vector<double>> PT_Reco(Recojets.size(), vector<double>());
            vector<vector<double>> PT_sort_number_only(Recojets.size(), vector<double>());
            vector<vector<double>> T_Reco_sort(Recojets.size(), vector<double>());
            vector<vector<double>> T_Reco(Recojets.size(), vector<double>());
            vector<vector<double>> velocity_jet(Recojets.size(), vector<double>());
            vector<vector<double>> T_sort_number_only(Recojets.size(), vector<double>());
            vector<double> mass_Reco = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
            vector<int> event_number_Reco = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
            vector<int> Eta_smaller_than_1_event = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
            vector<int> event_number_Reco_for_mass = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
            vector<int> Eta_smaller_than_1_event_for_mass = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
            for (unsigned int i = 0; i < Recojets.size(); i++)
            {
                double pt = Recojets[i].Pt();
                h_pt_ecal->Fill(pt);
                vector<float> ECALhits0, ECALhits31;
                for (unsigned int j1 = 0; j1 < Calhits.size(); j1++)
                {
                    LParticle hit = (LParticle)Calhits.at(j1);
                    TLorentzVector LE = hit.GetP();
                    double phi_h = LE.Phi();
                    double eta_h = LE.PseudoRapidity();
                    double mass = LE.M();
                    //cout << "mass =" << mass << endl;
                    vector<double> par = hit.GetParameters();
                    int type = hit.GetType(); // ECAL or HCAL
                    int layer = hit.GetCharge();
                    double Thit = par[0];
                    //cout << "Thit =" << Thit << endl;
                    double LogTime = TMath::Log10(Thit);
                    double Momentum = LE.P();
                    double Energy = LE.Energy();
                    //cout << "Momentum = " << Momentum << endl;
                    //cout << "Energy = " << Energy << endl;
                    double Vt = Momentum / Energy;

                    //cout << "Vt" << Vt << endl;
                    if (layer == 1 && (Recojets[i].DeltaR(LE)) < 0.4)
                    {
                        ECALhits0.push_back(Thit);
                    }
                    if (layer == 31 && (Recojets[i].DeltaR(LE)) < 0.4)
                    {
                        ECALhits31.push_back(Thit);
                        //Timing
                        T_Reco_sort[i].push_back(Thit);
                        T_Reco[i].push_back(Thit);
                        //Four_momentum
                        FourP_dR_Reco[i].push_back(LE);
                        //PT
                        PT_Reco_sort[i].push_back(LE.Perp());
                        PT_Reco[i].push_back(LE.Perp());
                        //Velocity

                        //Eta_selection
                        event_number_Reco[i] = event_number_Reco[i] + 1;
                        if (abs(eta_h) < 1)
                        {
                            Eta_smaller_than_1_event[i] = Eta_smaller_than_1_event[i] + 1;
                        }
                    }
                }

                //================================================================
                //  Save the Time_Difference to TH1D
                //================================================================
                double T_first = 0;
                double T_last = 0;
                if (ECALhits31.size() > 0 && ECALhits0.size() > 0)
                {
                    float Time_Difference = 0;
                    for (unsigned int G = 0; G < ECALhits0.size(); G++)
                    {
                        T_first = T_first + ECALhits0[G];
                    }
                    for (unsigned int H = 0; H < ECALhits31.size(); H++)
                    {
                        T_last = T_last + ECALhits31[H];
                    }
                    Time_Difference = (T_last / ECALhits31.size()) - (T_first / ECALhits0.size());
                    if (Time_Difference > 0)
                    {
                        Timing_detecto_ECAL_TDif->Fill(Time_Difference);
                    }
                }
            } //end of loop all Recojets

            //================================================================
            //       Sort from small to big (T_Reco_sort and PT_Reco_sort)
            //================================================================
            for (unsigned int i = 0; i < Recojets.size(); i++)
            {
                sort(T_Reco_sort[i].begin(), T_Reco_sort[i].end());
                sort(PT_Reco_sort[i].begin(), PT_Reco_sort[i].end());
            }
            //================================================================
            // Finally output PT_sort_number_only == PT_Reco_sort size
            //================================================================
            vector<vector<TLorentzVector>> Highest_PT_FourP(Recojets.size(), vector<TLorentzVector>());
            for (unsigned int i = 0; i < Recojets.size(); i++)
            {
                if (T_Reco_sort[i].size() == 0)
                    continue;
                PT_sort_number_only[i].push_back(PT_Reco_sort[i][0]);
                T_sort_number_only[i].push_back(T_Reco_sort[i][0]);
                for (unsigned int uuu = 0; uuu < T_Reco_sort[i].size(); uuu++)
                {
                    if (PT_Reco_sort[i][uuu] != PT_sort_number_only[i].back())
                    {
                        PT_sort_number_only[i].push_back(PT_Reco_sort[i][uuu]); //只存不同的pt
                    }
                    if (T_Reco_sort[i][uuu] != T_sort_number_only[i].back())
                    {
                        T_sort_number_only[i].push_back(T_Reco_sort[i][uuu]);
                    }
                }
            }
            vector<int> Full_contain;
            unsigned int check_point_eta = 0;
            //=================================
            // make sure that 90% eta are small than 1
            //=================================
            for (unsigned int i = 0; i < Recojets.size(); i++)
            {
                if (event_number_Reco[i] > 0)
                {
                    if ((Eta_smaller_than_1_event[i] / event_number_Reco[i]) > 0.9) //???
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

            //===============================
            //  Clear the Truth Level vector
            //===============================
            sjets_truth.clear();
            truthjets.clear();
            Truthjets_axis_index.clear();
            WW_boson.clear();
            avec_truth.clear();
            PDG_with_no_charge.clear();
            tGEN->Fill();
        }
        lcReader->close();
        delete lcReader;
    }
    RootFile->Write();
    h_jet_n_truth->Write();
    h_jet_pt_truth_check->Write();
    h_jet_eta_truth_check->Write();
    h_jet_n_truth->Write();
    h_jet_m_truth->Write();
    h_jet_nn_truth->Write();
    h_jet_n_truth->Write();
    h_jet_pt_reco_check->Write();
    h_jet_eta_reco_check->Write();
    h_jet_n_reco->Write();
    h_jet_m_reco->Write();
    h_jet_nn_reco->Write();
    h_pt_ecal->Write();
    Timing_detecto_ECAL_TDif->Write();
    h_effR_j1->Write();
    h_effR_j2->Write();
    h_effR_j3->Write();
    //RootFile->Print();
    RootFile->Close();
    return 0;
}
