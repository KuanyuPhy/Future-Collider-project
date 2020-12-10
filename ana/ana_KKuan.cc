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
  std::vector<std::string> files = open("data.in");
  int mtype = 1;

  //for (auto i = files.begin(); i != files.end(); ++i)
  //  cout << "files = " << i << endl;

  string outputfile = "root/output.root";
  //cout << "\n -> Output file is =" << outputfile << endl;
  int gg;
  TTree *mytree = new TTree("mytree", "this is a plain tree");
  mytree->Branch("gg", &gg, "gg/I");
  TH1F *Event_number_check = new TH1F("Event_number_check", "Event_number_check", 12, 0, 12);
  TH1F *Event_number_check_jet = new TH1F("Event_number_check_jet", "Event_number_check_jet", 12, 0, 12);
  TH1F *Eta_plot = new TH1F("Eta_plot", "Eta_plot", 200, -10, 10);
  TH1F *Eta_plot_after_cut = new TH1F("Eta_plot_after_cut", "Eta_plot_after_cut", 200, -10, 10);
  TH1F *Eta_plot_check = new TH1F("Eta_plot_check", "Eta_plot_check", 200, -10, 10);
  TH1F *Timing_Standard = new TH1F("Timing_Standard", "Timing_Standard", 200, 0, 50);
  TFile *RootFile = new TFile(outputfile.c_str(), "UPDATE", "Histogram file");
  TH1D *h_debug = new TH1D("debug", "events", 10, 0, 10);
  TH1D *h_e_jet = new TH1D("jet_energy", "energy jets", 1000, 0, 100);

  double TTxbins[] = {0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 6.0, 8.0, 10.0, 20.0, 30};
  const int TTBins = sizeof(TTxbins) / sizeof(double);

  TH1D *firstLayerTime = new TH1D("firstLayerTime", "First Layer Time", 100, 7.0, 8.0);
  firstLayerTime->GetXaxis()->SetTitle("T[ns]");
  firstLayerTime->GetYaxis()->SetTitle("Entries");

  TH1D *lastLayerTime = new TH1D("lastLayerTime", "last Layer Time", 100, 7.0, 8.0);
  lastLayerTime->GetXaxis()->SetTitle("T[ns]");
  lastLayerTime->GetYaxis()->SetTitle("Entries");

  TH1F *mass_sum_average = new TH1F("mass_sum_average", "mass_sum_average", 100, 0.5, 2.5);

  TH1D *h_jet_pt_truth_mass = new TH1D("Mass", "Jet mass", 100, 0, 10000);
  TH1F *Total_particle_ID_no_cut = new TH1F("Total_particle_ID_no_cut", "Total_particle_ID_no_cut", 20, 0, 20);
  TH1F *h_Jet_PT = new TH1F("h_Jet_PT", "h_Jet_PT", 100, 0, 3000);
  h_Jet_PT->GetXaxis()->SetTitle("GeV");
  TH1D *OwnLayerTime = new TH1D("OwnLayerTime", "Time difference between 31 nd 0 layers for 1 GeV", 1000, 0.0, 1.0);

  TH1F *Total_particle_ID_eta_cut = new TH1F("Total_particle_ID_eta_cut", "Total_particle_ID_eta_cut", 20, 0, 20);
  TH1D *TimeEcalDiff1 = new TH1D("time_ecal_diff1gev", "Time difference between 31 nd 0 layers for 1 GeV", TTBins - 1, TTxbins);
  TH1D *TimeEcalDiff10 = new TH1D("time_ecal_diff10gev", "Time difference between 31 nd 0 layers for 10 GeV", TTBins - 1, TTxbins);
  TProfile *TimeEcalDiffPt = new TProfile("time_ecal_diff vs pt", "Time difference between 31 nd 0 layers vs pT", 20, 0, 20, 's');

  TH1F *Total_particle_ID_eta_PT_cut = new TH1F("Total_particle_ID_eta_PT_cut", "Total_particle_ID_eta_PT_cut", 20, 0, 20);
  TH1F *check_jet_particle_number = new TH1F("check_jet_particle_number", "check_jet_particle_number", 100, 0, 100);
  TH1D *binsTT = new TH1D("bins_time_ecal", "bins_time_ecal", TTBins - 1, TTxbins);
  binsTT->Sumw2();
  for (Int_t j = 0; j < TTBins - 1; j++)
  {
    float x = TTxbins[j + 1] - TTxbins[j];
    binsTT->Fill(TTxbins[j] + 0.5 * x, x);
  }
  TH1F *Timing_detector_Leading = new TH1F("Timing_detector_Leading", "Timing_detector_Leading", 200, 0, 50);
  TH1F *Timing_detector_Trailing = new TH1F("Timing_detector_Trailing", "Timing_detector_Trailing", 200, 0, 50);
  TH1F *Timing_detector_Average = new TH1F("Timing_detector_Average", "Timing_detector_Average", 200, 0, 50);
  TH1F *Timing_detector_next_to_trailing = new TH1F("Timing_detector_next_to_trailing", "Timing_detector_next_to_trailing", 200, 0, 50);
  TH1F *Timing_detector_Trailing_P = new TH1F("Timing_detector_Trailing_P", "Timing_detector_Trailing_P", 50, 0, 200);
  TH1F *Timing_detector_next_to_trailing_P = new TH1F("Timing_detector_next_to_trailing_P", "Timing_detector_next_to_trailing_P", 50, 0, 100);
  TH1F *Timing_detector_Trailing_V = new TH1F("Timing_detector_Trailing_V", "Timing_detector_Trailing_V", 1000, 8, 9);
  TH1F *Timing_detector_next_to_trailing_V = new TH1F("Timing_detector_next_to_trailing_V", "Timing_detector_next_to_trailing_V", 1000, 0.9, 1);
  TH1F *Timing_detector_dR_Leading_trailing_PT = new TH1F("Timing_detector_dR_Leading_trailing_PT", "Timing_detector_dR_Leading_trailing_PT", 50, 0, 1);
  TH1F *Timing_detector_dR_Leading_next_trailing_PT = new TH1F("Timing_detector_dR_Leading_next_trailing_PT", "Timing_detector_dR_Leading_next_trailing_PT", 50, 0, 1);
  TH1F *Timing_detector_Average_tem = new TH1F("Timing_detector_Average_tem", "Timing_detector_Average_tem", 1000, 7, 12);
  TH1F *Timing_detector_dR_Leading_Proton_PT = new TH1F("Timing_detector_dR_Leading_Proton_PT", "Timing_detector_dR_Leading_Proton_PT", 50, 0, 1);

  TH2F *Timing_momentum_correlation = new TH2F("Timing_momentum_correlation", "Timing_momentum_correlation", 20, 0, 1, 20, 0, 1);
  TH2F *Timing_P_rank_difference_momentum_correlation = new TH2F("Timing_P_rank_difference_momentum_correlation", "Timing_P_rank_difference_momentum_correlation", 6, 0, 3, 40, -2, 2);
  //==============Trailing_particle_ID============//
  TH1F *Timing_detector_dR_Leading_trailing_T = new TH1F("Timing_detector_dR_Leading_trailing_T", "Timing_detector_dR_Leading_trailing_T", 50, 0, 1);
  TH1F *Timing_detector_dR_Leading_next_trailing_T = new TH1F("Timing_detector_dR_Leading_next_trailing_T", "Timing_detector_dR_Leading_next_trailing_T", 50, 0, 1);
  TH1F *Trailing_particle_ID_T = new TH1F("Trailing_particle_ID_T", "Trailing_particle_ID_T", 20, 0, 20);
  TH1F *Trailing_particle_ID_PT = new TH1F("Trailing_particle_ID_PT", "Trailing_particle_ID_PT", 20, 0, 20);
  TH1F *Next_to_trailing_particle_ID_T = new TH1F("Next_to_trailing_particle_ID_T", "Next_to_trailing_particle_ID_T", 20, 0, 20);
  TH1F *Next_to_trailing_particle_ID_PT = new TH1F("Next_to_trailing_particle_ID_PT", "Next_to_trailing_particle_ID_PT", 20, 0, 20);
  //=============TH2F_Trailing_ID_PT==============//

  TH1F *Timing_detector_Reco_TOF = new TH1F("Timing_detector_Reco_TOF", "Timing_detector_Reco_TOF", 200, 0, 50);
  TH1F *Timing_detector_Reco_TOF_track = new TH1F("Timing_detector_Reco_TOF_track", "Timing_detector_Reco_TOF_track", 200, 0, 50);

  TH1F *h_Particles_Rank_T[5];
  TH1F *h_Particles_Rank_PT[5];
  for (int j = 0; j < 5; j++)
  {
    h_Particles_Rank_T[j] = new TH1F(Form("h_Particles_Rank_T_%i", j), Form("h_Particles_Rank_T_%i", j), 20, 0, 20);
    h_Particles_Rank_PT[j] = new TH1F(Form("h_Particles_Rank_PT_%i", j), Form("h_Particles_Rank_PT_%i", j), 20, 0, 20);
  }

  TH2F *h_Particles_Rank_T_vs_PT[5];
  TH2F *h_Particles_Rank_PT_vs_PT[5];
  for (int j = 0; j < 5; j++)
  {
    h_Particles_Rank_T_vs_PT[j] = new TH2F(Form("h_Particles_Rank_T_vs_PT%i", j), Form("h_Particles_Rank_T_vs_PT%i", j), 20, 0, 20, 16, -2, 6);
    h_Particles_Rank_PT_vs_PT[j] = new TH2F(Form("h_Particles_Rank_PT_vs_PT%i", j), Form("h_Particles_Rank_PT_vs_PT%i", j), 20, 0, 20, 16, -2, 6);
  }
  TH2F *h_Particles_Rank_T_vs_T[5];
  TH2F *h_Particles_Rank_PT_vs_T[5];
  for (int j = 0; j < 5; j++)
  {
    h_Particles_Rank_T_vs_T[j] = new TH2F(Form("h_Particles_Rank_T_vs_T%i", j), Form("h_Particles_Rank_T_vs_T%i", j), 20, 0, 20, 500, 0, 50);
    h_Particles_Rank_PT_vs_T[j] = new TH2F(Form("h_Particles_Rank_PT_vs_T%i", j), Form("h_Particles_Rank_PT_vs_T%i", j), 20, 0, 20, 500, 0, 50);
  }

  TH1D *h_jet_pt_truth_sim = new TH1D("jet_pt_truth_sim", "pT [GeV] plus Geant4", 100, 0, 3000);
  TH1D *h_eta_jet = new TH1D("jet_eta", "eta jets", 100, -4, 4);
  TH1D *pz_gamma = new TH1D("pz_gamma", "#gamma pz [GeV]", 10, 0, 3000);
  pz_gamma->GetXaxis()->SetTitle("[GeV]");
  pz_gamma->GetYaxis()->SetTitle("Entries");

  TH1F *h_Jet_PT_HI_PT_C_P = new TH1F("h_Jet_PT_HI_PT_C_P", "h_Jet_PT_HI_PT_C_P", 100, 0, 1);
  TH1F *h_Jet_SM_PT_PT = new TH1F("h_Jet_SM_PT_PT", "h_Jet_SM_PT_PT", 80, -2, 6);
  TH1F *h_Jet_SM_PT_T = new TH1F("h_Jet_SM_PT_T", "h_Jet_SM_PT_T", 5000, 0, 50);
  TH1F *h_Jet_T_PT = new TH1F("h_Jet_T_PT", "h_Jet_T_PT", 80, -2, 6);
  TH1F *h_Jet_T_T = new TH1F("h_Jet_T_T", "h_Jet_T_T", 5000, 0, 50);
  TH1D *checkgs = new TH1D("checkgs", "checkgs", 40, 0, 4);

  TH1F *h_Particles_dR_Highest_PT_T[5];
  TH1F *h_Particles_dR_Highest_PT_PT[5];
  for (int j = 0; j < 5; j++)
  {
    h_Particles_dR_Highest_PT_T[j] = new TH1F(Form("h_Particles_dR_Highest_PT_T_%i", j), Form("h_Particles_dR_Highest_PT_T_%i", j), 100, 0, 1);
    h_Particles_dR_Highest_PT_PT[j] = new TH1F(Form("h_Particles_dR_Highest_PT_PT_%i", j), Form("h_Particles_dR_Highest_PT_PT_%i", j), 100, 0, 1);
  }
  TH1F *check_Pion_DZ = new TH1F("check_Pion_DZ", "check_Pion_DZ", 100, -5, 5);
  TH1F *check_Proton_DZ = new TH1F("check_Proton_DZ", "check_Proton_DZ", 100, -5, 5);
  TH1F *check_Pion_VZ = new TH1F("check_Pion_VZ", "check_Pion_VZ", 100, 0.0, 1);
  TH1F *check_Proton_V = new TH1F("check_Proton_VZ", "check_Proton_Vz", 100, 0.0, 1);
  TH1F *check_Pion_T = new TH1F("check_Pion_T", "check_Pion_T", 200, 0, 50);
  TH1F *check_Proton_T = new TH1F("check_Proton_T", "check_Proton_T", 200, 0, 50);

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

  // *tupleNames = "px:py:pz:vx:vy:vz:ex:ey:ez:pdg:np:nd:gs:ss";
  //TNtupleD *tuple = new TNtupleD("MCParticle", "", tupleNames);
  TNtuple *tuple = new TNtuple("tuple", "a n-tuple", "px:py:pz:vx:vy:vz:ex:ey:ez:pdg:np:nd:gs:ss");

  int var[13];

  double minPtConst = 0.2; // min pT on constituents
  // jets
  double Rparam = 0.4;
  const double ETAmax = 0.6;
  double ETmin = 0.5;
  if (mtype == 1)
    ETmin = 200; // increase for boosted

  if (mtype == 1)
    cout << "mu+mu- mode for boosted jets" << endl;

  cout << "min PT for jets=" << ETmin << endl;
  cout << "eta max for jets=" << ETAmax << endl;
  cout << "R for jets =" << Rparam << endl;

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

      // UTIL::LCTOOLS::dumpEvent( evt ) ;
      nnEvents++;
      nEvents++;
      int Status1 = 0;
      //if ((nEvents < 100 && nEvents % 10 == 0) || (nEvents > 100 && nEvents % 200 == 0))
      //cout << "nEvents: " << nEvents << endl;
      //cout << " # Events=" << nEvents << endl;
      Event_number_check->Fill(1);

      if (nEvents > MaxEvents)
        break;

      h_debug->Fill(1.0);

      std::string mcpName("MCParticle");
      // get truth
      IMPL::LCCollectionVec *col = (IMPL::LCCollectionVec *)evt->getCollection(mcpName);
      int nMCP = col->getNumberOfElements();

      int neu = 0;
      //list of Particle
      vector<int> PDG_with_no_charge = {0};
      vector<PseudoJet> avec_truth;     // created by generator
      vector<PseudoJet> avec_truth_sim; // also created by geant
                                        //==============================================================//
      int Forth = 1;
      int Back = -1;
      vector<int> Check_Forth_And_Back = {Forth, Back};
      vector<bool> Check_Forth_And_Back_Bool;
      //cout << "Check_Forth_And_Back_Bool size = " << Check_Forth_And_Back_Bool.size() << endl;
      vector<TLorentzVector> Forth_And_Back_Vector;
      int Zprime_pdg = 32;
      int Photon_pdg = 22;
      int Muon_pdg = 13;
      int W_pdg = 24;
      //cout << "nMCP = " << nMCP << endl;
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
              int gg = mcp->getGeneratorStatus();
              checkgs->Fill(gg);
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
      /*
      for (int i = 0; i < 10; i++)
      {
        cout << "i = " << i << "Check_Forth_And_Back_Bool  = " << Check_Forth_And_Back_Bool[i] << endl;
      }*/

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
          pz_gamma->Fill(pz);
          if (pz > 2000)
          {
            continue;
          }
        }
        /*if (abs(pdg) == 22 and gs == 1 and abs(mcp->getParents()[0]->getPDG()) == 13 and abs(mcp->getParents()[0]->getMomentum()[2]) > 2000)
        {
          continue;
        }*/
        if (gs == 1)
        {
          Status1 = Status1 + 1;
          avec_truth.push_back(p);
        }
      }
      //cout << "Status1 : " << Status1 << endl;
      // assume remnant for mu+mu-
      //cout << "test6" << endl;
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
      //cout << sjets_truth.size() << endl;
      for (unsigned int k = 0; k < sjets_truth.size(); k++)
      {
        //cout << "sjets_truth.size() = " << sjets_truth.size() << endl;
        //cout << "test8" << endl;
        if (sjets_truth[k].constituents().size() == 1 and abs(sjets_truth[k].constituents()[0].user_index()) == 22)
        { //cout << "Fucking ISR Photon" << endl;
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
        //
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
        //======================================
        //Since to remove constit is photon ??
        for (int i = 0; i < csize; i++)
        {
          if ((constit[i].user_index() - 22) != 0)
            Check_photon_jet = Check_photon_jet + 1;
        }
        if (Check_photon_jet == 0)
        {
          for (int i = 0; i < csize; i++)
          {
            //           cout << "constit_Photon jet:" << constit[i].user_index() << endl;
            continue;
          }
        }
        //======================================
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
        if (abs(p_using.Eta()) > 1.7)
        {
          //cout << "Particles of jet outside Eta==1" << endl;
          continue;
        }
        Eta_plot_after_cut->Fill(p_using.Eta());
        //======================================
        // fill jet substructure
        p.SetParameter(EffectiveRadius(sjets_truth[k], sjets_truth[k].constituents(), Rparam));  // 0
        p.SetParameter(nsubjettiness(sjets_truth[k], sjets_truth[k].constituents(), 1, Rparam)); // 1
        p.SetParameter(nsubjettiness(sjets_truth[k], sjets_truth[k].constituents(), 2, Rparam)); // 2
        p.SetParameter(nsubjettiness(sjets_truth[k], sjets_truth[k].constituents(), 3, Rparam)); //  3
        p.SetParameter(splittingscale(sjets_truth[k]));                                          // 4
        p.SetParameter(eccentricity(sjets_truth[k], sjets_truth[k].constituents()));             // 5
        truthjets.push_back(p);
        Truthjets_axis.push_back(p_using);

        for (int i = 0; i < csize; i++)
        {
          fastjet::PseudoJet constituent(constit[i].px(), constit[i].py(), constit[i].pz(), constit[i].e());
          TLorentzVector constit_vec;
          constit_vec.SetPxPyPzE(constit[i].px(), constit[i].py(), constit[i].pz(), constit[i].e());
          /*
          if (abs(constit_vec.Eta()) > 1)
          {
            Event_number_out_Eta2P1 = Event_number_out_Eta2P1 + 1;
            continue;
          }*/
          int it;
          it = find(Total_particle_kind.begin(), Total_particle_kind.end(), abs(constit[i].user_index()))[0];
          if (it != abs(constit[i].user_index()))
          {
            Total_particle_kind.push_back(abs(constit[i].user_index()));
          }
          for (unsigned int m = 0; m < Total_particle_kind.size(); m++)
          {
            if (abs(constit[i].user_index()) == Total_particle_kind[m])
              Total_particle_ID_eta_cut->Fill(m);
          }
        }
        //=========================================Cut and Find the information we want=================================//
        for (int i = 0; i < csize; i++)
        {
          fastjet::PseudoJet constituent(constit[i].px(), constit[i].py(), constit[i].pz(), constit[i].e());
          TLorentzVector constit_vec;
          constit_vec.SetPxPyPzE(constit[i].px(), constit[i].py(), constit[i].pz(), constit[i].e());
          int ID = 0;
          ID = find(PDG_with_no_charge.begin(), PDG_with_no_charge.end(), abs(constit[i].user_index()))[0];
          if (constit[i].perp() < 1.5 and ID != abs(constit[i].user_index()))
          {
            continue;
          }
          /*if (abs(constit_vec.Eta()) > 1)
          {
            Event_number_out_Eta2P1 = Event_number_out_Eta2P1 + 1;
            continue;
          }*/
          int it;
          it = find(Total_particle_kind.begin(), Total_particle_kind.end(), abs(constit[i].user_index()))[0];
          if (it != abs(constit[i].user_index()))
          {
            Total_particle_kind.push_back(abs(constit[i].user_index()));
          }
          for (unsigned int m = 0; m < Total_particle_kind.size(); m++)
          {
            if (abs(constit[i].user_index()) == Total_particle_kind[m])
              Total_particle_ID_eta_PT_cut->Fill(m);
          }
          mass_average = mass_average + constit_vec.M();

          int ID1 = 0;
          ID1 = find(PDG_with_no_charge.begin(), PDG_with_no_charge.end(), abs(constit[i].user_index()))[0];
          if (constit[i].perp() < 1.5 and ID != abs(constit[i].user_index()))
          {

            continue;
          }
          //cut off pt < 1.5 GeV
          if (constit[i].perp() < 1.5)
          {
            continue;
          }
          float constit_velocity = TMath::Power((TMath::Power(constit[i].px(), 2) + TMath::Power(constit[i].py(), 2) + TMath::Power(constit[i].pz(), 2)), 0.5) / constit[i].e(); //Beta
          float constit_velocity_vt = TMath::Power((TMath::Power(constit[i].px(), 2) + TMath::Power(constit[i].py(), 2) + TMath::Power(constit[i].pz(), 2)), 0.5) / constit[i].e();
          float constit_velocity_z = (constit[i].pz() / constit[i].e()); // Magnetic_consideration
          if (abs(sjets_truth[k].constituents()[i].user_index()) == 211)
          {
            //float dz = 2.3 / (TMath::Tan(constit_vec.Theta()));
            //cout << "pion dz = " << dz << endl;
            check_Pion_DZ->Fill(2.3 / (TMath::Tan(constit_vec.Theta())));
            //float vz = abs(constit_velocity_z);
            //cout << "pion vz = " << vz << endl;
            check_Pion_VZ->Fill(abs(constit_velocity_z));
            //float t = abs(2.3 * TMath::Power(10, 9) / (constit_velocity_z * SOL * (TMath::Tan(constit_vec.Theta()))));
            //cout << "pion t = " << t << endl;
            check_Pion_T->Fill(abs(2.3 * TMath::Power(10, 9) / (constit_velocity_z * SOL * (TMath::Tan(constit_vec.Theta())))));
          }
          if (abs(sjets_truth[k].constituents()[i].user_index()) == 2212)
          {
            check_Proton_DZ->Fill(2.3 / (TMath::Tan(constit_vec.Theta())));
            //float v = abs(constit_velocity_z);
            //cout << "Proton v = " << v << endl;
            check_Proton_V->Fill(abs(constit_velocity_z));
            check_Proton_T->Fill(abs(2.3 * TMath::Power(10, 9) / (constit_velocity_z * SOL * (TMath::Tan(constit_vec.Theta())))));
          }
          if (abs(sjets_truth[k].constituents()[i].user_index()) != 22)
          {
            time_average = time_average + abs(2.3 * TMath::Power(10, 9) / (constit_velocity_z * SOL * (TMath::Tan(constit_vec.Theta())))); //[ns]
            jet_time.push_back(abs(2.3 * TMath::Power(10, 9) / (constit_velocity_z * SOL * (TMath::Tan(constit_vec.Theta())))));
            jet_time_sort.push_back(abs(2.3 * TMath::Power(10, 9) / (constit_velocity_z * SOL * (TMath::Tan(constit_vec.Theta())))));

            jet_time_for_rank_sort.push_back(abs(2.3 * TMath::Power(10, 9) / (constit_velocity_z * SOL * (TMath::Tan(constit_vec.Theta())))));
            jet_time_for_rank.push_back(abs(2.3 * TMath::Power(10, 9) / (constit_velocity_z * SOL * (TMath::Tan(constit_vec.Theta())))));
            jet_P_for_rank.push_back(constit_vec.P());
            jet_P_for_rank_sort.push_back(constit_vec.P());
            Rank_PDGID.push_back(sjets_truth[k].constituents()[i].user_index());

            velocity_jet_sort.push_back(constit_velocity_vt);
            velocity_jet.push_back(constit_velocity_vt);
            velocity_jet_Z.push_back(constit[i].pz() / constit[i].e());
            velocity_jet_Theta.push_back(constit_vec.Theta());
            /*
            if (constit_vec.Eta() == 0)
            {
              //cout << "Timing_eta_0: " << abs(2.3 * TMath./A  ::Power(10, 9) / (SOL * TMath::Sin(constit_vec.Theta()))) << endl;
            }*/
            Timing_Standard->Fill(abs(2.3 * TMath::Power(10, 9) / (SOL * TMath::Sin(constit_vec.Theta())))); //Suppose all of them are photons.
            PT_jet.push_back(constit[i].perp());
            PT_jet_sort.push_back(constit[i].perp());
            momentum_jet.push_back(constit_vec.P());
            P_jet_sort.push_back(constit_vec.P());
            constit_PDG.push_back(abs(sjets_truth[k].constituents()[i].user_index()));
            FourP.push_back(constit_vec);
            FourP_1.push_back(constituent);
            event_number = event_number + 1;
            Eta_plot_check->Fill(constit_vec.Eta());
            h_e_jet->Fill(constit[i].perp());
          }
        }
        //cout << "Jet_each_event: " << Jet_each_event << endl;
        Event_number_check_jet->Fill(Jet_each_event);

        if (event_number == 0)
          continue;
        check_jet_particle_number->Fill(jet_time_for_rank_sort.size());
        //cout << "event_number: " << event_number << endl;
        //cout << "Event_number_out_Eta2P1: " << Event_number_out_Eta2P1 << endl;
        /*
        if (event_number > 0 and Event_number_out_Eta2P1 > 0) ///Event_number_out_Eta2P1 almost = 0
        {
          cout << "Event_number_out_Eta2P1: " << Event_number_out_Eta2P1 << endl;
          cout << "event_number: " << (event_number + Event_number_out_Eta2P1) << endl;
        }*/
        mass_average = mass_average / event_number;
        time_average = time_average / event_number;
        mass_sum_average->Fill(mass_average);
        double max_time = *max_element(jet_time.begin(), jet_time.end());
        double min_time = *min_element(jet_time.begin(), jet_time.end());
        double max_perp = *max_element(PT_jet.begin(), PT_jet.end());
        double min_perp = *min_element(PT_jet.begin(), PT_jet.end());
        //////////////////////////////////////////////
        sort(PT_jet_sort.begin(), PT_jet_sort.end());
        sort(jet_time_sort.begin(), jet_time_sort.end());
        sort(velocity_jet_sort.begin(), velocity_jet_sort.end());
        sort(P_jet_sort.begin(), P_jet_sort.end());
        sort(jet_time_for_rank_sort.begin(), jet_time_for_rank_sort.end());
        sort(jet_P_for_rank_sort.begin(), jet_P_for_rank_sort.end());
        int Trailing_ID_T = 0;          //One jet one trailing ID(T)
        int Trailing_ID_PT = 0;         //One jet one trailing ID(PT)
        int Next_to_trailing_ID_T = 0;  //One jet one trailing ID(T)
        int Next_to_trailing_ID_PT = 0; //One jet one trailing ID(PT)
        float Momentum_Trailing = 0;
        float Momentum_Next_to_Trailing = 0;
        float PT_Trailing = 0;
        float Theta_Leading = 0;
        float Theta_Trailing = 0;
        float Theta_Next_to_Trailing = 0;
        float Vz_Trailing = 0;
        float Vz_Next_to_Trailing = 0;
        float Vz_Leading = 0;
        vector<int> Particle_ID_T;
        vector<int> Particle_ID_T_Re;
        vector<int> Particle_ID_PT;
        vector<int> Particle_ID_PT_Re;
        vector<int> H_Particle_ID_T;
        vector<int> H_Particle_ID_PT;

        vector<float> Particle_ID_T_PT;
        vector<float> Particle_ID_T_PT_Re;
        vector<float> Particle_ID_T_T;
        vector<float> Particle_ID_T_T_Re;
        vector<float> Particle_ID_PT_T;
        vector<float> Particle_ID_PT_T_Re;
        vector<float> Particle_ID_PT_PT;
        vector<float> Particle_ID_PT_PT_Re;

        vector<float> dR_Highest_PT_T;
        vector<float> dR_Highest_PT_PT;
        vector<float> Timing_checking;
        vector<float> PT_checking;
        int check_timing_number = 0;
        int check_PT_number = 0;
        vector<TLorentzVector> HighestPT_Trailing_and_next_trailing;
        //===============================Find_minimum_velocity_in_particle=========================//
        for (unsigned int i = 0; i < PT_jet.size(); i++)
        {
          if (max_perp == PT_jet[i])
          {
            HighestPT_Trailing_and_next_trailing.push_back(FourP[i]);
          }
        }
        for (unsigned int j = 0; j < PT_jet.size(); j++)
        {
          for (unsigned int i = 0; i < PT_jet.size(); i++)
          {
            if (PT_jet_sort[j] == PT_jet[i])
            {
              //cout << "PT_jet_sort[j] = " << PT_jet_sort[j] << endl;
              //cout << "PT_jet[i] = " << PT_jet[i] << endl;
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
            if (velocity_jet_sort[j] == velocity_jet[i])
            {
              check_timing_number = check_timing_number + 1;
              Timing_checking.push_back(jet_time[i]);
              Particle_ID_T_PT.push_back(PT_jet[i]);
              Particle_ID_T_T.push_back(jet_time[i]);
              Particle_ID_T.push_back(constit_PDG[i]);
              dR_Highest_PT_T.push_back(HighestPT_Trailing_and_next_trailing[0].DeltaR(FourP[i]));

              if (j < 5)
              {
                h_Particles_dR_Highest_PT_T[j]->Fill(HighestPT_Trailing_and_next_trailing[0].DeltaR(FourP[i]));
              }
            }
          }
        }
        if (PT_jet_sort.size() > 0)
        {
          //cout << "Log10(Particle_ID_PT_PT[0]) " << TMath::Log10(Particle_ID_PT_PT[0]) << endl;
          //cout << "Log10(Particle_ID_PT_PT[1]) " << TMath::Log10(Particle_ID_PT_PT[1]) << endl;
          h_Jet_SM_PT_PT->Fill(TMath::Log10(Particle_ID_PT_PT[0]));
          h_Jet_SM_PT_T->Fill(Particle_ID_PT_T[0]);
          h_Jet_T_PT->Fill(TMath::Log10(Particle_ID_T_PT[0]));
          h_Jet_T_T->Fill(Particle_ID_T_T[0]);
          //??
          for (unsigned int joke = 0; joke < PT_jet_sort.size(); joke++)
          {
            int Check = 0;

            for (unsigned int joke_1 = 0; joke_1 < PT_jet_sort.size(); joke_1++)
            {
              //cout << "PT_jet_sort.size() - 1 = " << PT_jet_sort.size() << endl;
              if (PT_jet_sort[PT_jet_sort.size() - 1 - joke] == PT_jet[joke_1])
              {
                int ID3 = 0;
                ID3 = find(PDG_with_no_charge.begin(), PDG_with_no_charge.end(), abs(constit_PDG[joke_1]))[0];
                if (ID3 != abs(constit_PDG[joke_1]))
                {
                  //cout << "abs(constit_PDG[joke_1]) = " << abs(constit_PDG[joke_1]) << endl;
                  //cout << "PT_jet_sort[PT_jet_sort.size()-1-joke]/Jet_PT: " << PT_jet_sort[PT_jet_sort.size()-1-joke]/Jet_PT << endl;
                  h_Jet_PT_HI_PT_C_P->Fill(PT_jet_sort[PT_jet_sort.size() - 1 - joke] / Jet_PT);
                  //cout << "2: " << endl;
                  Check = Check + 1;
                  break;
                }
                else
                {
                  {
                    continue;
                  }
                }
              }
              if (Check == 1)
              {
                break;
              }
            }
          } //???
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
              if (abs(Particle_ID_PT[j]) == Trailing_particle_kind_limit_pT[m] and abs(Particle_ID_PT[j]) != 22 and Num_of_Tra_pT < 5)
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
          else
          {
            if ((Particle_ID_T_PT[j]) == Particle_ID_PT_PT_Re[0] or (Particle_ID_T_PT[j]) == Particle_ID_PT_PT_Re[1] or (Particle_ID_T_PT[j]) == Particle_ID_PT_PT_Re[2] or (Particle_ID_T_PT[j]) == Particle_ID_PT_PT_Re[3] or (Particle_ID_T_PT[j]) == Particle_ID_PT_PT_Re[4])
            {
              //cout << "??? " << endl;
              //cout << "Overlap Output TMath::Log10(Particle_ID_T_PT[j]): " << TMath::Log10(Particle_ID_T_PT[j]) << endl;
              continue;
            }
            else
            {
              {
                for (unsigned int m = 0; m < Trailing_particle_kind_limit.size(); m++)
                {
                  if (abs(Particle_ID_T[j]) == Trailing_particle_kind_limit[m] and abs(Particle_ID_T[j]) != 22 and Num_of_Tra_V < 5)
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
          }
        }
        if (Num_of_Tra_V != 5)
        {
          continue;
        }
        for (int j = 0; j < 5; j++)
        {
          //cout << "Particle_ID_PT_PT_Re: " << Particle_ID_PT_PT_Re[j] << endl;
          //cout << "TMath::Log10(Particle_ID_PT_PT_Re: " << TMath::Log10(Particle_ID_PT_PT_Re[j]) << endl;
          //cout << "Particle_ID_T_PT_Re: " << Particle_ID_T_PT_Re[j] << endl;
          //cout << "TMath::Log10(Particle_ID_T_PT_Re: " << TMath::Log10(Particle_ID_T_PT_Re[j]) << endl;
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
        //==============================================================================
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
          if (abs(Trailing_ID_T) == Trailing_particle_kind[m] and abs(Trailing_ID_T) != 22)
            Trailing_particle_ID_T->Fill(m);
        }
        int it2;
        it2 = find(Trailing_particle_kind.begin(), Trailing_particle_kind.end(), abs(Trailing_ID_PT))[0];
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
          if (abs(Trailing_ID_PT) == Trailing_particle_kind[m] and abs(Trailing_ID_PT) != 22)
            Trailing_particle_ID_PT->Fill(m);
        }
        int it4;
        it4 = find(Trailing_particle_kind.begin(), Trailing_particle_kind.end(), abs(Next_to_trailing_ID_T))[0];
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
          if (abs(Next_to_trailing_ID_T) == Trailing_particle_kind[m] and abs(Next_to_trailing_ID_T) != 22)
            Next_to_trailing_particle_ID_T->Fill(m);
        }
        int it3;
        it3 = find(Trailing_particle_kind.begin(), Trailing_particle_kind.end(), abs(Next_to_trailing_ID_PT))[0];
        if (it3 != abs(Next_to_trailing_ID_PT))
        {
          Trailing_particle_kind.push_back(abs(Next_to_trailing_ID_PT));
        }
        else
          //cout << "" << endl; //

          for (unsigned int m = 0; m < Trailing_particle_kind.size(); m++)
          {
            //cout << "Trailing_particle_kind[m]: "<< abs(Trailing_particle_kind[m]) << endl;
            if (abs(Next_to_trailing_ID_PT) == Trailing_particle_kind[m] and abs(Next_to_trailing_ID_PT) != 22)
              Next_to_trailing_particle_ID_PT->Fill(m);
            //cout << "Next_to_trailing_particle_ID_PT->GetBinContent(m): " << Next_to_trailing_particle_ID_PT->GetBinContent(m+1) << endl;
          }
        int event_checkk = 0;
        for (unsigned int j = 0; j < jet_time_for_rank_sort.size(); j++)
        {
          for (unsigned int i = 0; i < jet_time_for_rank.size(); i++)
          {
            if (jet_time_for_rank_sort[(jet_time_for_rank_sort.size() - 1) - j] == jet_time_for_rank[i])
            {
              for (unsigned int kk = 0; kk < jet_P_for_rank.size(); kk++)
              {
                if (jet_P_for_rank[i] == jet_P_for_rank_sort[kk])
                {
                  Timing_momentum_correlation->Fill(float((j + 1)) / float((jet_time_for_rank_sort.size())), float((kk + 1)) / float((jet_time_for_rank_sort.size())));

                  Timing_P_rank_difference_momentum_correlation->Fill(TMath::Log10(jet_P_for_rank[kk]), float(j - kk) / float((jet_time_for_rank_sort.size())));
                  event_checkk = event_checkk + 1;
                }
                else
                {
                  continue;
                }
              }
            }
            else
            {
              continue;
            }
          }
        }
        for (unsigned int i = 0; i < constit_PDG.size(); i++)
        {
          if (constit_PDG[i] == 2212)
          {
            Timing_detector_dR_Leading_Proton_PT->Fill(HighestPT_Trailing_and_next_trailing[0].DeltaR(FourP[i]));
          }
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
            PT_Trailing = PT_jet[i];
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
        Timing_detector_Trailing_P->Fill(abs(Momentum_Trailing));
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
        velocity_jet.clear();
        Timing_detector_Average_tem->Clear();
        if (pt >= 10000)
        {
          for (int i = 0; i < csize; i++)
          {
            // std::cout << "Truth Jets const.Eti20(): " << (constit[i].Et()) << std::endl;
            //            h_truth_jets_const_Et10->Fill(constit[i].Et());
          }
        }
      }

      //cout << "Final: " << endl;

    } //end loop
    cout << "nnEvents = " << nnEvents << endl;
    lcReader->close();
    delete lcReader;
  }

  RootFile->Write();
  mytree->Write();
  Event_number_check->Write();
  Eta_plot_after_cut->Write();
  Eta_plot_check->Write();
  Eta_plot->Write();
  check_Pion_DZ->Write();
  check_Pion_VZ->Write();
  check_Pion_T->Write();
  check_Proton_DZ->Write();
  check_Proton_V->Write();
  check_Proton_T->Write();
  Timing_Standard->Write();
  //h_Particles_dR_Highest_PT_PT->Write();

  //RootFile->Print();
  RootFile->Close();
  cout << "end" << endl;

  return 0;
}