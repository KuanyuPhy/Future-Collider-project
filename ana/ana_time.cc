/***************************************************************************
 *  How to use SLCIO files from HepSim, and how to build anti-KT jets 
 *  S.Chekanov (ANL) chekanov@anl.gov
 *  A library for HEP events storage and processing based on Google's PB   
 *  The project web site: http://atlaswww.hep.anl.gov/hepsim/
****************************************************************************/

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
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include <fastjet/AreaDefinition.hh>
#include <fastjet/ClusterSequenceArea.hh>
using namespace std;

const double kPI = TMath::Pi();
const double k2PI = 2 * kPI;

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

// main example
int main(int argc, char **argv)
{

  std::vector<std::string> files = open("data.in");
  int mtype = 2; // single-particles

  string outputfile = "root/output.root";
  cout << "\n -> Output file is =" << outputfile << endl;
  TFile *RootFile = new TFile(outputfile.c_str(), "RECREATE", "Histogram file");
  TH1D *h_debug = new TH1D("debug", "events", 10, 0, 10);
  TH1D *h_jet_pt_truth_check = new TH1D("h_jet_pt_truth_check", "pT [GeV] ", 270, 0, 2700); // plus Geant4
  h_jet_pt_truth_check->GetXaxis()->SetTitle("P_{T}^{true,jet}");

  TH1D *h_jet_eta_truth_check = new TH1D("h_jet_eta_truth_check", "#eta", 100, -5, 5); //before eta cut
  TH1D *h_jet_n_truth = new TH1D("h_jet_n_truth", "Nr of truth jets", 10, 0, 10);      //before match
  TH1D *h_jet_nn_truth = new TH1D("h_jet_nn_truth", "Nr of truth jets", 10, 0, 10);    //after match
  TH1D *h_jet_m_truth = new TH1D("jet_m_truth", "Mass [GeV]", 100, 0, 100);

  // read detector geometry for this configuration
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

  // calculate position
  double layer_size_mm = 27.5 + 5 + 2.5;
  double hcal_inner_R_mm = 2300;

  double layer_size_ecal_mm = 4.0;
  double ecal_inner_R_mm = 2100;

  // calculate weighted average of time
  double XTime = 0;
  double XEnergy = 0;

  //pandora::Pandora* m_pandora = new pandora::Pandora();
  //m_pandora->ReadSettings(detector);

  const pandora::SubDetectorType subDetectorType = geom->getPandoraSubDetectorType(caloType);
  //const pandora::SubDetector &subdet(pandora.GetGeometry()->GetSubDetector(subDetectorType));

  // pandora::Pandora* pandora = new pandora::Pandora();
  //const pandora::SubDetector &subdet(pandora->GetGeometry()->GetSubDetector(subDetectorType));
  //  PandoraApi::Geometry::LayerParametersList layer_par= xsubdet->m_layerParametersList;

  const char *tupleNames = "px:py:pz:vx:vy:vz:ex:ey:ez:pdg:np:nd:gs:ss";
  TNtupleD *tuple = new TNtupleD("MCParticle", "", tupleNames);

  double minPtConst = 0.2; // min pT on constituents
  // jets
  double Rparam = 0.5;
  const double ETAmax = 1.0;
  double ETmin = 0.5;
  if (mtype == 1)
    ETmin = 200; // increase for boosted
  if (mtype == 31)
    ETmin = 50;
  if (mtype == 32)
    ETmin = 2000;
  if (mtype == 33)
    ETmin = 8000;

  if (mtype == 1)
    cout << "mu+mu- mode for boosted jets" << endl;
  if (mtype == 2)
    cout << "single-particle mode for pgun" << endl;
  if (mtype == 3)
    cout << "pp mode for jet resolutions" << endl;

  cout << "min PT for jets=" << ETmin << endl;
  cout << "eta max for jets=" << ETAmax << endl;
  cout << "R for jets =" << Rparam << endl;

  Float_t Timediff1tree;
  Float_t Timediff31tree;
  TTree *tGEN = new TTree("tGEN", "Tree with vectors");
  tGEN->Branch("Timediff1tree", &Timediff1tree);
  tGEN->Branch("Timediff31tree", &Timediff31tree);
  // fastjet
  Strategy strategy = fastjet::Best;
  JetDefinition jet_def(fastjet::antikt_algorithm, Rparam, strategy);

  //double m,px,py,pz,e,vx,vy,vz,ex,ey,ez ;
  //int pdg,np,nd,gs,ss ;

  int MaxEvents = 100000000;

  int nEvents = 0;

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
      // get truth
      vector<fastjet::PseudoJet> WW_boson;
      vector<PseudoJet> avec_truth;     // created by generator
      vector<PseudoJet> avec_truth_sim; // also created by geant
      vector<int> PDG_with_no_charge = {0};
      IMPL::LCCollectionVec *col = (IMPL::LCCollectionVec *)evt->getCollection("MCParticle");
      int nMCP = col->getNumberOfElements();
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
      // assume remnant for mu+mu-
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
      bool isTrueMatch = false;
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
            isTrueMatch = true;
            //cout << "k =" << k << endl;
          }
          if (isTrueMatch)
          {
            break;
          }
        }
      }
      //==============================================================
      //          Get Collection from EcalBarrelHits
      //==============================================================
      double ecalsum_raw = 0;
      vector<PseudoJet> avec_hits_raw;
      vector<LParticle> simhits;
      IMPL::LCCollectionVec *col53 = (IMPL::LCCollectionVec *)evt->getCollection("EcalBarrelHits");
      int nCL = col53->getNumberOfElements();
      for (int i = 0; i < nCL; ++i)
      {
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
        // fill hits
        LParticle p(px, py, pz, e, layer);
        p.SetCharge(layer);
        p.SetType(2);         // ECAL
        p.SetParameter(Thit); //0 fastest hit
        if (layer == 1 or layer == 31)
        {
          simhits.push_back(p);
        }
      } //end of EcalBarrelHits loop
      // ECAL hits after created by slicPandora
      double ecalsum = 0;
      vector<PseudoJet> avec_hits;
      vector<LParticle> Calhits;
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
        double Thit = mcp->getTime();
        LParticle p(px, py, pz, e, layer);
        ecalsum = ecalsum + e;
        PseudoJet pj(px, py, pz, e);
        p.SetCharge(layer);
        p.SetType(2);         // ECAL
        p.SetParameter(Thit); //0
        if (layer == 1 or layer == 31)
        {
          Calhits.push_back(p);
        }
        if (e > 0.5)
        {
          avec_hits.push_back(pj);
        }
      } //end of EM_BARREL loop
      //================================================================
      //  Check Hit time
      //================================================================
      for (unsigned i = 0; i < simhits.size(); i++)
      {
        LParticle hit = (LParticle)simhits.at(i);
        vector<double> par = hit.GetParameters();
        int Simlayer = hit.GetCharge();
        double SimThit = par[0];
        for (unsigned j = 0; j < Calhits.size(); j++)
        {
          LParticle p = (LParticle)Calhits.at(j);
          vector<double> par1 = p.GetParameters();
          int Callayer = p.GetCharge();
          double CalThit = par1[0];
          double Timediff = SimThit - CalThit;
          if (Simlayer == 1 && Callayer == 1)
          {
            double Timediff1 = SimThit - CalThit;
            Timediff1tree = Timediff1;
            //cout << "Timediff =" << Timediff << endl;
          }
          if (Simlayer == 31 && Callayer == 31)
          {
            double Timediff31 = SimThit - CalThit;
            Timediff31tree = Timediff31;
            //cout << "Timediff =" << Timediff << endl;
          }
        } //end of Calhits loop
      }   //end of simhits loop
      tGEN->Fill();
    }
    lcReader->close();
    delete lcReader;
    // close the file
  }
  //m_ntuple->Fill();
  RootFile->Write();
  h_jet_pt_truth_check->Write();
  h_jet_eta_truth_check->Write();
  h_jet_n_truth->Write();
  h_jet_m_truth->Write();
  //RootFile->Print();
  RootFile->Close();
  return 0;
}
