// This code is modefied to calculate asymmetry for individual beams and averaged for total asymmetry

#include <iostream>
#include <fstream>
#include <cmath>
//#include "std_lib_facilities.h"
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH1.h"
#include "TMath.h"
#include "TGraphErrors.h"
using namespace std;
ofstream Output;
// void asymmetryVspT_averaged(const char *ifile = "Ntuple_RawTree_18063008_iff_data.root")
void asymmetryVspT_final(const char *ifile = "../Ntuple_RawTree_polCut_TPCorTOF_P123.root")
// void asymmetryVspT_averaged(const char *ifile = "NTuple_P123_V1_July18.root")
{
    // Histograms to get mean polarization for yellow and blue beam
    TH1D *hpolB = new TH1D("hpolB", "", 100, 0, 1);
    TH1D *hpolY = new TH1D("hpolY", "", 100, 0, 1);
    // Chi-sqare distribution
    TH1D *h_chi_BG = new TH1D("h_chi_BG", "h_chi_BG", 5, 0, 5);
    TH1D *h_chi_BL = new TH1D("h_chi_BL", "h_chi_BL", 5, 0, 5);
    TH1D *h_chi_YG = new TH1D("h_chi_YG", "h_chi_YG", 5, 0, 5);
    TH1D *h_chi_YL = new TH1D("h_chi_YL", "h_chi_YL", 5, 0, 5);
    // pT distribution for each pT bin for average pT
    TH1D *hMinv[5];
    TH1D *hMinvGtB[5];
    TH1D *hMinvLtB[5];
    TH1D *hMinvGtY[5];
    TH1D *hMinvLtY[5];

    for (Int_t n = 0; n < 5; n++)
    {
        hMinv[n] = new TH1D(Form("Minv_mbin_%i", n), "", 100, 0, 5);
        hMinvGtB[n] = new TH1D(Form("MinvGtB_mbin_%i", n), "", 100, 0, 5);
        hMinvLtB[n] = new TH1D(Form("MinvLtB_mbin_%i", n), "", 100, 0, 5);
        hMinvGtY[n] = new TH1D(Form("MinvGtY_mbin_%i", n), "", 100, 0, 5);
        hMinvLtY[n] = new TH1D(Form("MinvLtY_mbin_%i", n), "", 100, 0, 5);
    }

    TFile *f = new TFile(ifile);
    TTree *ntuple1 = (TTree *)f->Get("ntuple1_TPCorTOF"); // get trees
    TTree *ntuple2 = (TTree *)f->Get("ntuple2_TPCorTOF");
    TTree *ntuple4 = (TTree *)f->Get("ntuple4_TPCorTOF");
    TTree *ntuple5 = (TTree *)f->Get("ntuple5_TPCorTOF");
    // TTree *ntuple6 = (TTree *)f->Get("ntuple6_TPCorTOF");
    // define variables to hold trees variables
    Output.open("AUT_Vs_pT_9Bin_trigBias_polCut_TPCorTOF.txt");

    float eta_pair;
    float fspinconfig;
    float cone;
    float pT_pair;
    float Minv;
    float PhiRSB;
    float PhiRSY;
    float fitPts_min_pair;
    float polB_corr, polY_corr;
    double pi = 3.14159265359;

    double Npairs[9] = {0};
    double pTpairs[9] = {0};
    double Mpairs[9] = {0};
    double pT1[9] = {0}; // for pT  bin average
    double pT2[9] = {0};
    double pT3[9] = {0};
    double pT4[9] = {0};
    double pT5[9] = {0};
    double avg_Minv[5] = {0};
    double NpairsBGt[5][9] = {0};
    double NpairsBLt[5][9] = {0};
    double pTpairsBGt[5][9] = {0};
    double pTpairsBLt[5][9] = {0};
    double MpairsBGt[5][9] = {0};
    double MpairsBLt[5][9] = {0};
    double EtapairsBGt[5][9] = {0};
    double EtapairsBLt[5][9] = {0};

    double NpairsYGt[5][9] = {0};
    double NpairsYLt[5][9] = {0};
    double pTpairsYGt[5][9] = {0};
    double pTpairsYLt[5][9] = {0};
    double MpairsYGt[5][9] = {0};
    double MpairsYLt[5][9] = {0};
    double EtapairsYGt[5][9] = {0};
    double EtapairsYLt[5][9] = {0};

    for (int k = 0; k < 9; k++)
    {
        Npairs[k] = 0.;
        pTpairs[k] = 0.;
        Mpairs[k] = 0.;
        pT1[k] = 0;
        pT2[k] = 0;
        pT3[k] = 0;
        pT4[k] = 0;
        pT5[k] = 0;
    }

    // get the values of variables from tree

    ntuple1->SetBranchAddress("fspinconfig", &fspinconfig);
    ntuple1->SetBranchAddress("cone", &cone);
    ntuple2->SetBranchAddress("Minv", &Minv);
    ntuple2->SetBranchAddress("pT_pair", &pT_pair);
    ntuple2->SetBranchAddress("eta_pair", &eta_pair);
    ntuple4->SetBranchAddress("PhiRSB", &PhiRSB);
    ntuple4->SetBranchAddress("PhiRSY", &PhiRSY);
    ntuple5->SetBranchAddress("fitPts_min_pair", &fitPts_min_pair);
    // ntuple6->SetBranchAddress("polB_corr", &polB_corr);

    // To store average polarization values from histograms
    double avgPolB, avgPolY, rmsB, rmsY, avgPolT, rmsT;

    Int_t nentries = (Int_t)ntuple1->GetEntries(); // entries of trees
    // add friend to the tree. no root file to add friend means the tree is on the same root file
    ntuple1->AddFriend("ntuple2_TPCorTOF");
    ntuple1->AddFriend("ntuple4_TPCorTOF");
    ntuple1->AddFriend("ntuple5_TPCorTOF");
    // double pT[10]={2.80, 3.47,  3.79, 4.12, 4.49, 4.938,  5.505, 6.300, 7.66, 15 };
    // double eta_range[10]={-1.200,   -0.668,  -0.469, -0.281, -0.096, 0.089, 0.275, 0.47,   0.675,  1.2 };
    Double_t avgMinv[5] = {0};
    // double pT[6]={2.80, 3.727, 4.343, 5.157, 6.529, 15.00};
    // double M[10]={0.250, 0.403, 0.516, 2, 0.711, 0.803, 0.921, 1.070, 1.286, 4.000};
    double M[6] = {0.20, 0.4848, 0.6565, 0.8431, 1.1710, 4.000}; // Navagyan's bining for 90% of Run17 data
    // double M[6]={0.20, 0.4795, 0.6714, 0.8452, 1.1100, 4.00}; // Babu binning
    // pT binning for each Minv-Bin
    double p1[10] = {2.60, 3.982, 4.453, 4.914, 5.407, 5.982, 6.7080, 7.740, 9.540, 25};
    double p2[10] = {2.600, 4.1065, 4.620, 5.1220, 5.6850, 6.350, 7.200, 8.400, 10.500, 25};
    double p3[10] = {2.600, 4.183, 4.737, 5.299, 5.920, 6.657, 7.598, 8.939, 11.260, 25};
    double p4[10] = {2.600, 4.160, 4.715, 5.285, 5.933, 6.710, 7.710, 9.120, 11.584, 25};
    double p5[10] = {2.600, 4.615, 5.286, 5.957, 6.705, 7.590, 8.708, 10.280, 12.990, 25};

    /*double p1[10]={2.80, 3.432, 3.709, 3.9884, 4.3098, 4.6781, 5.147, 5.809, 6.966, 15};
      double p2[10]={2.80, 3.432, 3.7192, 4.0147, 4.3474, 4.745, 5.255, 5.973, 7.217, 15};
      double p3[10]={2.80, 3.4294, 3.734, 4.0515, 4.4108, 4.8396, 5.386, 6.153, 7.468, 15};
      double p4[10]={2.80, 3.3631, 3.6523, 3.9558, 4.3046, 4.729, 5.2784, 6.056, 7.392, 15};
      double p5[10]={2.80, 3.7464, 4.144, 4.538, 4.9705, 5.4766, 6.111, 6.9860, 8.459, 15};
      */

    // int ne = (int)ntuple6->GetEntries();
    /*for (int pol=0; pol<ne;pol++)
    //for (int pol=0; pol<10000;pol++)
    {
    ntuple6->GetEntry(pol);
    hpolB->Fill(polB_corr);
    hpolY->Fill(polY_corr);
    }
    */
    // TCanvas *fitCanvas[5];
    TCanvas *FitCanv_Blue_Gt = new TCanvas("FitCanv_Blue_Gt", "FitCanv_Blue_Gt", 900, 900);
    TCanvas *FitCanv_Blue_Lt = new TCanvas("FitCanv_Blue_Lt", "FitCanv_Blue_Lt", 900, 900);
    TCanvas *FitCanv_Yellow_Gt = new TCanvas("FitCanv_Yellow_Gt", "FitCanv_YLllow_Gt", 900, 900);
    TCanvas *FitCanv_Yellow_Lt = new TCanvas("FitCanv_Yellow_Lt", "FitCanv_Yellow_Lt", 900, 900);
    FitCanv_Blue_Gt->Divide(3, 3);
    FitCanv_Blue_Lt->Divide(3, 3);
    FitCanv_Yellow_Gt->Divide(3, 3);
    FitCanv_Yellow_Lt->Divide(3, 3);

    double A_pT1bg[9] = {0};
    double deltaA_pT1bg[9] = {0};
    double A_pT1bl[9] = {0};
    double deltaA_pT1bl[9] = {0};
    double A_pT1yg[9] = {0};
    double deltaA_pT1yg[9] = {0};
    double A_pT1yl[9] = {0};
    double deltaA_pT1yl[9] = {0};

    double A_pT2bg[9] = {0};
    double deltaA_pT2bg[9] = {0};
    double A_pT2bl[9] = {0};
    double deltaA_pT2bl[9] = {0};
    double A_pT2yg[9] = {0};
    double deltaA_pT2yg[9] = {0};
    double A_pT2yl[9] = {0};
    double deltaA_pT2yl[9] = {0};

    double A_pT3bg[9] = {0};
    double deltaA_pT3bg[9] = {0};
    double A_pT3bl[9] = {0};
    double deltaA_pT3bl[9] = {0};
    double A_pT3yg[9] = {0};
    double deltaA_pT3yg[9] = {0};
    double A_pT3yl[9] = {0};
    double deltaA_pT3yl[9] = {0};

    double A_pT4bg[9] = {0};
    double deltaA_pT4bg[9] = {0};
    double A_pT4bl[9] = {0};
    double deltaA_pT4bl[9] = {0};
    double A_pT4yg[9] = {0};
    double deltaA_pT4yg[9] = {0};
    double A_pT4yl[9] = {0};
    double deltaA_pT4yl[9] = {0};

    double A_pT5bg[9] = {0};
    double deltaA_pT5bg[9] = {0};
    double A_pT5bl[9] = {0};
    double deltaA_pT5bl[9] = {0};
    double A_pT5yg[9] = {0};
    double deltaA_pT5yg[9] = {0};
    double A_pT5yl[9] = {0};
    double deltaA_pT5yl[9] = {0};

    for (Int_t mbin = 0; mbin < 5; mbin++)
    {
        int NpTGtUpB[165] = {0};
        int NpTLtUpB[165] = {0};
        int NpTGtDnB[165] = {0};
        int NpTLtDnB[165] = {0};
        int NMinvGtUpB[165] = {0};
        int NMinvLtUpB[165] = {0};
        int NMinvGtDnB[165] = {0};
        int NMinvLtDnB[165] = {0};
        int NetaUpB[165] = {0};
        int NetaDnB[165] = {0};
        int NpTGtUpY[165] = {0};
        int NpTLtUpY[165] = {0};
        int NpTGtDnY[165] = {0};
        int NpTLtDnY[165] = {0};
        int NMinvGtUpY[165] = {0};
        int NMinvLtUpY[165] = {0};
        int NMinvGtDnY[165] = {0};
        int NMinvLtDnY[165] = {0};
        int NetaUpY[165] = {0};
        int NetaDnY[165] = {0};
        int NpTGtUpT[165] = {0};
        int NpTLtUpT[165] = {0};
        int NpTGtDnT[165] = {0};
        int NpTLtDnT[165] = {0};
        int NMinvGtUpT[165] = {0};
        int NMinvLtUpT[165] = {0};
        int NMinvGtDnT[165] = {0};
        int NMinvLtDnT[165] = {0};
        int NetaUpT[165] = {0};
        int NetaDnT[165] = {0};

        double pT[10] = {0};
        // choose correspondig pT bin for 5 different  Minv-Bin
        if (mbin == 0)
            for (int k = 0; k < 10; k++)
            {
                pT[k] = p1[k];
            }
        if (mbin == 1)
            for (int k = 0; k < 10; k++)
            {
                pT[k] = p2[k];
            }
        if (mbin == 2)
            for (int k = 0; k < 10; k++)
            {
                pT[k] = p3[k];
            }
        if (mbin == 3)
            for (int k = 0; k < 10; k++)
            {
                pT[k] = p4[k];
            }
        if (mbin == 4)
            for (int k = 0; k < 10; k++)
            {
                pT[k] = p5[k];
            }

        cout << "mbin: " << mbin << ", pT[10]={" << pT[0] << ", " << pT[1] << ", " << pT[2] << ", " << pT[3] << ", " << pT[4] << ", " << pT[5] << ", " << pT[6] << ", " << pT[7] << ", " << pT[8] << ", " << pT[9] << "}" << endl;
        ;
        cout << "Entries: " << nentries << endl;
        for (Int_t i = 0; i < nentries; i++)
        // for (Int_t i = 0; i < 10000; i++)
        {

            ntuple1->GetEntry(i);
            if (Minv < M[mbin] || Minv >= M[mbin + 1])
                continue;
            // if(pT_pair<3.75) continue;
            if (cone >= .7)
                continue;
            if (Minv > 4.)
                continue;
            //	if(!(Minv>0.5076 || Minv<0.4876))  continue;//dodge K0 mass range, doesn't cause asymmetry. Implemention was not working before.This is fine now.
            if (fitPts_min_pair < 15)
                continue;

            hMinv[mbin]->Fill(Minv);
            // cout << "polarizztion values   "<< polB_corr << "  "<<polY_corr <<endl;
            // cout << "selection cuts are good ...working on phi loop.... " << endl;

            //-------------BLUE-------------------------------------
            // Phi
            for (int phi = 0; phi < 16; phi++)
            {
                if (PhiRSB >= (phi - 8.) / 8. * pi && PhiRSB <= (phi - 7.) / 8. * pi)
                {
                    // avg pT pair loop
                    for (int m = 0; m < 9; m++)
                    {
                        if (pT_pair >= pT[m] && pT_pair < pT[m + 1])
                        {
                            Npairs[m] = Npairs[m] + 1;
                            pTpairs[m] = pTpairs[m] + pT_pair;
                            Mpairs[m] = Mpairs[m] + Minv;
                            if (eta_pair > 0)
                            {
                                hMinvGtB[mbin]->Fill(Minv);
                                NpairsBGt[mbin][m] = NpairsBGt[mbin][m] + 1;
                                pTpairsBGt[mbin][m] = pTpairsBGt[mbin][m] + pT_pair;
                                MpairsBGt[mbin][m] = MpairsBGt[mbin][m] + Minv;
                                EtapairsBGt[mbin][m] = EtapairsBGt[mbin][m] + eta_pair;
                            }
                            if (eta_pair < 0)
                            {
                                hMinvLtB[mbin]->Fill(Minv);
                                NpairsBLt[mbin][m] = NpairsBLt[mbin][m] + 1;
                                pTpairsBLt[mbin][m] = pTpairsBLt[mbin][m] + pT_pair;
                                MpairsBLt[mbin][m] = MpairsBLt[mbin][m] + Minv;
                                EtapairsBLt[mbin][m] = EtapairsBLt[mbin][m] + eta_pair;
                            }

                            if (fspinconfig == 51 || fspinconfig == 53)
                            {
                                if (eta_pair > 0)
                                {
                                    NpTGtUpB[m * 16 + phi]++;
                                }
                                if (eta_pair < 0)
                                {
                                    NpTLtUpB[m * 16 + phi]++;
                                }
                            }
                            if (fspinconfig == 83 || fspinconfig == 85)
                            {
                                if (eta_pair > 0)
                                {
                                    NpTGtDnB[m * 16 + phi]++;
                                }
                                if (eta_pair < 0)
                                {
                                    NpTLtDnB[m * 16 + phi]++;
                                }
                            }
                        }
                    } // pT loop

                }                                                           // phi-range loop
            }                                                               // phi loop
            /**************** BLUE BEAM ENDS *****************************/ ///////

            ////***************** YELLOW BEAM *****************************////////////
            // Phi loop
            for (int phi = 0; phi < 16; phi++)
            {
                if (PhiRSY >= (phi - 8.) / 8. * pi && PhiRSY <= (phi - 7.) / 8. * pi)
                {
                    // pT loop
                    for (int m = 0; m < 9; m++)
                    {
                        if (pT_pair >= pT[m] && pT_pair < pT[m + 1])
                        {

                            if (eta_pair > 0)
                            {
                                hMinvGtY[mbin]->Fill(Minv);
                                NpairsYGt[mbin][m] = NpairsYGt[mbin][m] + 1;
                                pTpairsYGt[mbin][m] = pTpairsYGt[mbin][m] + pT_pair;
                                MpairsYGt[mbin][m] = MpairsYGt[mbin][m] + Minv;
                                EtapairsYGt[mbin][m] = EtapairsYGt[mbin][m] + eta_pair;
                            }
                            if (eta_pair < 0)
                            {
                                hMinvLtY[mbin]->Fill(Minv);
                                NpairsYLt[mbin][m] = NpairsYLt[mbin][m] + 1;
                                pTpairsYLt[mbin][m] = pTpairsYLt[mbin][m] + pT_pair;
                                MpairsYLt[mbin][m] = MpairsYLt[mbin][m] + Minv;
                                EtapairsYLt[mbin][m] = EtapairsYLt[mbin][m] + eta_pair;
                            }

                            if (fspinconfig == 51 || fspinconfig == 83)
                            {
                                // if(eta_pair>0)
                                if (eta_pair < 0)
                                {
                                    NpTGtUpY[m * 16 + phi]++;
                                }
                                // if(eta_pair<0)
                                if (eta_pair > 0)
                                {
                                    NpTLtUpY[m * 16 + phi]++;
                                }
                            }
                            if (fspinconfig == 53 || fspinconfig == 85)
                            {
                                // if(eta_pair>0)
                                if (eta_pair < 0)
                                {
                                    NpTGtDnY[m * 16 + phi]++;
                                }
                                // if(eta_pair<0)
                                if (eta_pair > 0)
                                {
                                    NpTLtDnY[m * 16 + phi]++;
                                }
                            }
                        }
                    }                                                    // pT loop
                }                                                        // PhiRS range loop
            }                                                            // phi loop
        }                                                                // entries  for loop
        /******************** YELLOW BEAM ENDS ************************/ /////

        // Get average polarizations and rms from histograms
        // avgPolB = hpolB->GetMean();
        // avgPolY = hpolY->GetMean();
        // rmsB    = hpolB->GetRMS();
        // rmsY    = hpolY->GetRMS();
        avgPolB = 0.59;
        avgPolY = 0.5948;
        rmsB = 0.03148;
        rmsY = 0.03697;

        // Store pT bin averages for each pT bin in an array for plotting purpose

        // Here you are taking average pT for Blue and Yellow and Average all the same and plot them....might need to reconsider this one......
        if (mbin == 0)
        {
            for (int k0 = 0; k0 < 9; k0++)
            {
                pT1[k0] = pTpairs[k0] / (double)Npairs[k0];
            }
        }
        if (mbin == 1)
        {
            for (int k1 = 0; k1 < 9; k1++)
            {
                pT2[k1] = pTpairs[k1] / (double)Npairs[k1];
            }
        }
        if (mbin == 2)
        {
            for (int k2 = 0; k2 < 9; k2++)
            {
                pT3[k2] = pTpairs[k2] / (double)Npairs[k2];
            }
        }
        if (mbin == 3)
        {
            for (int k3 = 0; k3 < 9; k3++)
            {
                pT4[k3] = pTpairs[k3] / (double)Npairs[k3];
            }
        }
        if (mbin == 4)
        {
            for (int k4 = 0; k4 < 9; k4++)
            {
                pT5[k4] = pTpairs[k4] / (double)Npairs[k4];
            }
        }

        double dAdB, dAdC, dAdD, dAdE, dAdP; // variables for error calculation
        // Use rms for polarization error from distribution!
        double dP_B = rmsB; // Polarization errors blue
        double dP_Y = rmsY; // Polarization errors yellow
        double dA[165];     // Asymmetry amplitude from sine fit
        double B, C, D, E;
        double a, b;
        double Asym[165];
        // double pi = 3.14159265359;

        // fitCanvas[mbin] = new TCanvas(Form("fitCanvas_MBin%i", mbin), "", 650, 550);
        // fitCanvas[mbin]->Print(Form("./FitPlots_Mbin_%i_AsymVspT.pdf(", mbin));

        double Abg[9] = {0}; // Here Abg refers to Asymmetry for Blue in eta >0
        double deltaAbg[9] = {0};
        double Ayg[9] = {0};
        double deltaAyg[9] = {0};
        // Asymmetry for BLUE eta >0

        for (int m = 0; m < 9; m++)
        {
            for (int ang = 0; ang < 16; ang++)
            {
                if (ang < 8)
                {
                    a = sqrt(double(NpTGtUpB[m * 16 + ang]) * double(NpTGtDnB[m * 16 + ang + 8]));
                    b = sqrt(double(NpTGtDnB[m * 16 + ang]) * double(NpTGtUpB[m * 16 + ang + 8]));
                    B = double(NpTGtUpB[m * 16 + ang]);
                    C = double(NpTGtUpB[m * 16 + ang + 8]);
                    D = double(NpTGtDnB[m * 16 + ang]);
                    E = double(NpTGtDnB[m * 16 + ang + 8]);
                }
                if (ang > 7)
                {
                    a = sqrt(double(NpTGtUpB[m * 16 + ang]) * double(NpTGtDnB[m * 16 + ang - 8]));
                    b = sqrt(double(NpTGtDnB[m * 16 + ang]) * double(NpTGtUpB[m * 16 + ang - 8]));
                    B = double(NpTGtUpB[m * 16 + ang]);
                    C = double(NpTGtUpB[m * 16 + ang - 8]);
                    D = double(NpTGtDnB[m * 16 + ang]);
                    E = double(NpTGtDnB[m * 16 + ang - 8]);
                }
                Asym[ang] = (1. / avgPolB) * ((a - b) / (a + b));
                dAdB = (1. / avgPolB) * E * sqrt(D * C) / (sqrt(B * E) * (pow((sqrt(B * E) + sqrt(D * C)), 2)));
                dAdE = (1. / avgPolB) * B * sqrt(D * C) / (sqrt(B * E) * (pow((sqrt(B * E) + sqrt(D * C)), 2)));
                dAdD = (-1. / avgPolB) * C * sqrt(B * E) / (sqrt(D * C) * (pow((sqrt(B * E) + sqrt(D * C)), 2)));
                dAdC = (-1. / avgPolB) * D * sqrt(B * E) / (sqrt(D * C) * (pow((sqrt(B * E) + sqrt(D * C)), 2)));
                dAdP = (-1. / (avgPolB * avgPolB)) * (sqrt(B * E) - sqrt(D * C)) / (sqrt(B * E) + sqrt(D * C));
                dA[ang] = sqrt(pow((fabs(dAdB) * sqrt(B)), 2) + pow((fabs(dAdC) * sqrt(C)), 2) + pow((fabs(dAdD) * sqrt(D)), 2) + pow((fabs(dAdE) * sqrt(E)), 2) + pow((fabs(dAdP) * dP_B), 2));
                // cout <<"polB:"<<avgPolB<<", polY: "<<avgPolY<< ", a="<<a<<", b="<<b<<", B="<<B<<", C="<<C<<", D="<<D<<", E="<<E<<", Asym["<<ang<<"]="<<Asym[ang]<<", dA["<<ang<<"]="<<dA[ang]<<endl;
            } // angle loop
            FitCanv_Blue_Gt->cd(m + 1);
            char name[600];
            char title[600];
            double chi2Ndf[9];
            double angle[16] = {-15. / 16. * pi, -13. / 16. * pi, -11. / 16. * pi, -9. / 16. * pi, -7. / 16. * pi, -5. / 16. * pi, -3. / 16. * pi, -1. / 16. * pi, 1. / 16. * pi, 3. / 16. * pi, 5. / 16. * pi, 7. / 16. * pi, 9. / 16. * pi, 11. / 16. * pi, 13. / 16. * pi, 15. / 16. * pi};
            double ex[16] = {0};
            gStyle->SetOptDate(0);
            gStyle->SetOptFit(1);
            auto grt = new TGraphErrors(16, angle, Asym, ex, dA);
            grt->SetMarkerStyle(20);
            grt->Draw("AP");
            TF1 *fit = new TF1("fit", "[0]*sin(x)", -3.14159265359, 0);
            // fit->SetParameter(0,0.0001);
            grt->Fit(fit, "R");
            grt->SetMarkerColor(4);
            grt->GetXaxis()->SetTitle("#Phi_{RS}");
            grt->GetXaxis()->SetTitleOffset(1);
            sprintf(title, "Minv bin %i, pT bin %i, BLUE, #eta > 0", mbin, m);
            grt->SetTitle(title);
            chi2Ndf[m] = fit->GetChisquare() / fit->GetNDF();
            cout << "ChiSquare = " << chi2Ndf[m] << endl;
            grt->GetYaxis()->SetTitle("A_{UT}(#Phi_{RS})");
            grt->GetYaxis()->SetTitleOffset(1);
            grt->GetYaxis()->CenterTitle();
            grt->GetYaxis()->SetRangeUser(-0.08, 0.08);
            Abg[m] = fit->GetParameter(0); // This is where it is extracting A_UT from the sin fit for each pT bin
            deltaAbg[m] = fit->GetParError(0);
            h_chi_BG->Fill(chi2Ndf[m]);
            // sprintf(name,"Minv%i_pTbin%i_#etaGt_BLUE_9bin.png",mbin, m);
            // c1->SaveAs(name);
            // fitCanvas[mbin]->Print(Form("./TriggerBias/FitPlots_Mbin_%i_AsymVspT.pdf", mbin));
            FitCanv_Blue_Gt->Update();

        } // 9 pT bin loop
        if (mbin == 0)
        {
            FitCanv_Blue_Gt->Print("./TriggerBias/FitPlot_Mbin_Blue_Gt.pdf(", "pdf");
        }
        else if (mbin == 4)
        {
            FitCanv_Blue_Gt->Print("./TriggerBias/FitPlot_Mbin_Blue_Gt.pdf)", "pdf");
        }
        else
        {
            FitCanv_Blue_Gt->Print("./TriggerBias/FitPlot_Mbin_Blue_Gt.pdf", "pdf");
        }
        // FitCanv_Blue->Print("./TriggerBias/FitPlot_Blue_Gt.pdf");
        // YELLOW eta > 0
        for (int m = 0; m < 9; m++)
        {
            for (int ang = 0; ang < 16; ang++)
            {
                if (ang < 8)
                {
                    a = sqrt(double(NpTGtUpY[m * 16 + ang]) * double(NpTGtDnY[m * 16 + ang + 8]));
                    b = sqrt(double(NpTGtDnY[m * 16 + ang]) * double(NpTGtUpY[m * 16 + ang + 8]));
                    B = double(NpTGtUpY[m * 16 + ang]);
                    C = double(NpTGtUpY[m * 16 + ang + 8]);
                    D = double(NpTGtDnY[m * 16 + ang]);
                    E = double(NpTGtDnY[m * 16 + ang + 8]);
                }
                if (ang > 7)
                {
                    a = sqrt(double(NpTGtUpY[m * 16 + ang]) * double(NpTGtDnY[m * 16 + ang - 8]));
                    b = sqrt(double(NpTGtDnY[m * 16 + ang]) * double(NpTGtUpY[m * 16 + ang - 8]));
                    B = double(NpTGtUpY[m * 16 + ang]);
                    C = double(NpTGtUpY[m * 16 + ang - 8]);
                    D = double(NpTGtDnY[m * 16 + ang]);
                    E = double(NpTGtDnY[m * 16 + ang - 8]);
                }
                Asym[ang] = (1. / avgPolY) * ((a - b) / (a + b));

                dAdB = (1. / avgPolY) * E * sqrt(D * C) / (sqrt(B * E) * (pow((sqrt(B * E) + sqrt(D * C)), 2)));
                dAdE = (1. / avgPolY) * B * sqrt(D * C) / (sqrt(B * E) * (pow((sqrt(B * E) + sqrt(D * C)), 2)));
                dAdD = (-1. / avgPolY) * C * sqrt(B * E) / (sqrt(D * C) * (pow((sqrt(B * E) + sqrt(D * C)), 2)));
                dAdC = (-1. / avgPolY) * D * sqrt(B * E) / (sqrt(D * C) * (pow((sqrt(B * E) + sqrt(D * C)), 2)));
                dAdP = (-1. / (avgPolY * avgPolY)) * (sqrt(B * E) - sqrt(D * C)) / (sqrt(B * E) + sqrt(D * C));
                dA[ang] = sqrt(pow((fabs(dAdB) * sqrt(B)), 2) + pow((fabs(dAdC) * sqrt(C)), 2) + pow((fabs(dAdD) * sqrt(D)), 2) + pow((fabs(dAdE) * sqrt(E)), 2) + pow((fabs(dAdP) * dP_Y), 2));
            }
            FitCanv_Yellow_Gt->cd(m + 1);
            char name[600];
            char title[600];
            double chi2Ndf[9];
            double angle[16] = {-15. / 16. * pi, -13. / 16. * pi, -11. / 16. * pi, -9. / 16. * pi, -7. / 16. * pi, -5. / 16. * pi, -3. / 16. * pi, -1. / 16. * pi, 1. / 16. * pi, 3. / 16. * pi, 5. / 16. * pi, 7. / 16. * pi, 9. / 16. * pi, 11. / 16. * pi, 13. / 16. * pi, 15. / 16. * pi};
            double ex[16] = {0};
            gStyle->SetOptDate(0);
            gStyle->SetOptFit(1);
            auto gr = new TGraphErrors(16, angle, Asym, ex, dA);
            gr->SetMarkerStyle(20);
            gr->Draw("AP");
            TF1 *fit = new TF1("fit", "[0]*sin(x)", -3.14159265359, 0);
            // fit->SetParameter(0,0.0001);
            gr->Fit(fit, "R");
            gr->SetMarkerColor(4);
            gr->GetXaxis()->SetTitle("#Phi_{RS}");
            gr->GetXaxis()->SetTitleOffset(1);
            sprintf(title, "Minv bin %i, pT bin %i, YELLOW, #eta > 0", mbin, m);
            gr->SetTitle(title);
            chi2Ndf[m] = fit->GetChisquare() / fit->GetNDF();
            cout << "ChiSquare= " << chi2Ndf[m] << endl;
            gr->GetYaxis()->SetTitle("A_{UT}(#Phi_{RS})");
            gr->GetYaxis()->SetTitleOffset(1);
            gr->GetYaxis()->CenterTitle(kTRUE);
            gr->GetYaxis()->SetRangeUser(-0.08, 0.08);
            Ayg[m] = fit->GetParameter(0); // This is where it extract it A_UT for Yellow beam eta>0 for each pT pair bin
            deltaAyg[m] = fit->GetParError(0);
            h_chi_YG->Fill(chi2Ndf[m]);
            FitCanv_Yellow_Gt->Update();
            // sprintf(name,"Minv%i_pTbin%i_EtaGt_YELLOW_9bin.png",mbin, m);
            // c1->SaveAs(name);
            // fitCanvas[mbin]->Print(Form("./FitPlots_Mbin_%i_AsymVspT.pdf", mbin));
        } // pT loop
        if (mbin == 0)
        {
            FitCanv_Yellow_Gt->Print("./TriggerBias/FitPlot_Mbin_Yellow_Gt.pdf(", "pdf");
        }
        else if (mbin == 4)
        {
            FitCanv_Yellow_Gt->Print("./TriggerBias/FitPlot_Mbin_Yellow_Gt.pdf)", "pdf");
        }
        else
        {
            FitCanv_Yellow_Gt->Print("./TriggerBias/FitPlot_Mbin_Yellow_Gt.pdf", "pdf");
        }
        // BLUE beam in eta < 0

        double Abl[9] = {0};
        double deltaAbl[9] = {0}; /// bl refers to blue beam eta less than 0
        double Ayl[9] = {0};
        double deltaAyl[9] = {0};
        for (int m = 0; m < 9; m++)
        {
            for (int ang = 0; ang < 16; ang++)
            {
                if (ang < 8)
                {
                    a = sqrt(double(NpTLtUpB[m * 16 + ang]) * double(NpTLtDnB[m * 16 + ang + 8]));
                    b = sqrt(double(NpTLtDnB[m * 16 + ang]) * double(NpTLtUpB[m * 16 + ang + 8]));
                    B = double(NpTLtUpB[m * 16 + ang]);
                    C = double(NpTLtUpB[m * 16 + ang + 8]);
                    D = double(NpTLtDnB[m * 16 + ang]);
                    E = double(NpTLtDnB[m * 16 + ang + 8]);
                }
                if (ang > 7)
                {
                    a = sqrt(double(NpTLtUpB[m * 16 + ang]) * double(NpTLtDnB[m * 16 + ang - 8]));
                    b = sqrt(double(NpTLtDnB[m * 16 + ang]) * double(NpTLtUpB[m * 16 + ang - 8]));
                    B = double(NpTLtUpB[m * 16 + ang]);
                    C = double(NpTLtUpB[m * 16 + ang - 8]);
                    D = double(NpTLtDnB[m * 16 + ang]);
                    E = double(NpTLtDnB[m * 16 + ang - 8]);
                }
                Asym[ang] = (1. / avgPolB) * ((a - b) / (a + b));
                dAdB = (1. / avgPolB) * E * sqrt(D * C) / (sqrt(B * E) * (pow((sqrt(B * E) + sqrt(D * C)), 2)));
                dAdE = (1. / avgPolB) * B * sqrt(D * C) / (sqrt(B * E) * (pow((sqrt(B * E) + sqrt(D * C)), 2)));
                dAdD = (-1. / avgPolB) * C * sqrt(B * E) / (sqrt(D * C) * (pow((sqrt(B * E) + sqrt(D * C)), 2)));
                dAdC = (-1. / avgPolB) * D * sqrt(B * E) / (sqrt(D * C) * (pow((sqrt(B * E) + sqrt(D * C)), 2)));
                dAdP = (-1. / (avgPolB * avgPolB)) * (sqrt(B * E) - sqrt(D * C)) / (sqrt(B * E) + sqrt(D * C));
                dA[ang] = sqrt(pow((fabs(dAdB) * sqrt(B)), 2) + pow((fabs(dAdC) * sqrt(C)), 2) + pow((fabs(dAdD) * sqrt(D)), 2) + pow((fabs(dAdE) * sqrt(E)), 2) + pow((fabs(dAdP) * dP_B), 2));
            } // angle loop
            FitCanv_Blue_Lt->cd(m + 1);
            char name[600];
            char title[600];
            double chi2Ndf[9];
            double angle[16] = {-15. / 16. * pi, -13. / 16. * pi, -11. / 16. * pi, -9. / 16. * pi, -7. / 16. * pi, -5. / 16. * pi, -3. / 16. * pi, -1. / 16. * pi, 1. / 16. * pi, 3. / 16. * pi, 5. / 16. * pi, 7. / 16. * pi, 9. / 16. * pi, 11. / 16. * pi, 13. / 16. * pi, 15. / 16. * pi};
            double ex[16] = {0};
            gStyle->SetOptDate(0);
            gStyle->SetOptFit(1);
            auto grt = new TGraphErrors(16, angle, Asym, ex, dA);
            grt->SetMarkerStyle(20);
            grt->Draw("AP");
            TF1 *fit = new TF1("fit", "[0]*sin(x)", -3.14159265359, 0);
            // fit->SetParameter(0,0.0001);
            grt->Fit(fit, "R");
            grt->SetMarkerColor(4);
            grt->GetXaxis()->SetTitle("#Phi_{RS}");
            grt->GetXaxis()->SetTitleOffset(1);
            sprintf(title, "Minv bin %i, pT bin %i, BLUE, #eta < 0", mbin, m);
            grt->SetTitle(title);
            chi2Ndf[m] = fit->GetChisquare() / fit->GetNDF();
            cout << "ChiSquare = " << chi2Ndf[m] << endl;
            grt->GetYaxis()->SetTitle("A_{UT}(#Phi_{RS})");
            grt->GetYaxis()->SetTitleOffset(1);
            grt->GetYaxis()->CenterTitle(kTRUE);
            grt->GetYaxis()->SetTitleOffset(1);
            grt->GetYaxis()->SetRangeUser(-0.08, 0.08);
            Abl[m] = fit->GetParameter(0); //    This is where it extract the A_UT for blue beam for eta<0 for each pT bin
            deltaAbl[m] = fit->GetParError(0);
            h_chi_BL->Fill(chi2Ndf[m]);

            FitCanv_Blue_Lt->Update();
            // sprintf(name,"Minv%i_pTbin%i_#etaLt_BLUE_9bin.png",mbin, m);
            // c1->Print(name);
            // fitCanvas[mbin]->Print(Form("./FitPlots_Mbin_%i_AsymVspT.pdf", mbin));
        } // invariant mass loop
        if (mbin == 0)
        {
            FitCanv_Blue_Lt->Print("./TriggerBias/FitPlot_Mbin_Blue_Lt.pdf(", "pdf");
        }
        else if (mbin == 4)
        {
            FitCanv_Blue_Lt->Print("./TriggerBias/FitPlot_Mbin_Blue_Lt.pdf)", "pdf");
        }
        else
        {
            FitCanv_Blue_Lt->Print("./TriggerBias/FitPlot_Mbin_Blue_Lt.pdf", "pdf");
        }
        // YELLOW Beam Eta<0
        for (int m = 0; m < 9; m++)
        {
            for (int ang = 0; ang < 16; ang++)
            {
                if (ang < 8)
                {
                    a = sqrt(double(NpTLtUpY[m * 16 + ang]) * double(NpTLtDnY[m * 16 + ang + 8]));
                    b = sqrt(double(NpTLtDnY[m * 16 + ang]) * double(NpTLtUpY[m * 16 + ang + 8]));
                    B = double(NpTLtUpY[m * 16 + ang]);
                    C = double(NpTLtUpY[m * 16 + ang + 8]);
                    D = double(NpTLtDnY[m * 16 + ang]);
                    E = double(NpTLtDnY[m * 16 + ang + 8]);
                }
                if (ang > 7)
                {
                    a = sqrt(double(NpTLtUpY[m * 16 + ang]) * double(NpTLtDnY[m * 16 + ang - 8]));
                    b = sqrt(double(NpTLtDnY[m * 16 + ang]) * double(NpTLtUpY[m * 16 + ang - 8]));
                    B = double(NpTLtUpY[m * 16 + ang]);
                    C = double(NpTLtUpY[m * 16 + ang - 8]);
                    D = double(NpTLtDnY[m * 16 + ang]);
                    E = double(NpTLtDnY[m * 16 + ang - 8]);
                }
                Asym[ang] = (1. / avgPolY) * ((a - b) / (a + b));

                dAdB = (1. / avgPolY) * E * sqrt(D * C) / (sqrt(B * E) * (pow((sqrt(B * E) + sqrt(D * C)), 2)));
                dAdE = (1. / avgPolY) * B * sqrt(D * C) / (sqrt(B * E) * (pow((sqrt(B * E) + sqrt(D * C)), 2)));
                dAdD = (-1. / avgPolY) * C * sqrt(B * E) / (sqrt(D * C) * (pow((sqrt(B * E) + sqrt(D * C)), 2)));
                dAdC = (-1. / avgPolY) * D * sqrt(B * E) / (sqrt(D * C) * (pow((sqrt(B * E) + sqrt(D * C)), 2)));
                dAdP = (-1. / (avgPolY * avgPolY)) * (sqrt(B * E) - sqrt(D * C)) / (sqrt(B * E) + sqrt(D * C));
                dA[ang] = sqrt(pow((fabs(dAdB) * sqrt(B)), 2) + pow((fabs(dAdC) * sqrt(C)), 2) + pow((fabs(dAdD) * sqrt(D)), 2) + pow((fabs(dAdE) * sqrt(E)), 2) + pow((fabs(dAdP) * dP_Y), 2));
            }

            FitCanv_Yellow_Lt->cd(m + 1);
            char name[600];
            char title[600];
            double chi2Ndf[9];
            // double angle[16]={-15./16.*pi,-13./16.*pi,-11./16.*pi,-9./16.*pi,-7./16.*pi,-5./16.*pi,-3./16.*pi,-1./16.*pi};
            double angle[16] = {-15. / 16. * pi, -13. / 16. * pi, -11. / 16. * pi, -9. / 16. * pi, -7. / 16. * pi, -5. / 16. * pi, -3. / 16. * pi, -1. / 16. * pi, 1. / 16. * pi, 3. / 16. * pi, 5. / 16. * pi, 7. / 16. * pi, 9. / 16. * pi, 11. / 16. * pi, 13. / 16. * pi, 15. / 16. * pi};
            double ex[16] = {0};
            gStyle->SetOptDate(0);
            gStyle->SetOptFit(1);
            auto gr = new TGraphErrors(16, angle, Asym, ex, dA);
            gr->SetMarkerStyle(20);
            gr->Draw("AP");
            TF1 *fit = new TF1("fit", "[0]*sin(x)", -3.14159265359, 0);
            fit->SetParameter(0, 0.0001);
            gr->Fit(fit, "R");
            gr->SetMarkerColor(4);
            gr->GetXaxis()->SetTitle("#Phi_{RS}");
            gr->GetXaxis()->SetTitleOffset(1);
            sprintf(title, "Minv bin %i, pT bin %i, YELLOW, #eta < 0", mbin, m);
            chi2Ndf[m] = fit->GetChisquare() / fit->GetNDF();
            cout << "ChiSquare= " << chi2Ndf[m] << endl;
            gr->SetTitle(title);
            gr->GetYaxis()->SetTitle("A_{UT}(#Phi_{RS})");
            gr->GetYaxis()->CenterTitle(kTRUE);
            gr->GetYaxis()->SetTitleOffset(1);
            gr->GetYaxis()->SetRangeUser(-0.08, 0.08);
            Ayl[m] = fit->GetParameter(0); // This is where it extrat it AUT for yellow beam for Eta<0 for each pT bin
            deltaAyl[m] = fit->GetParError(0);
            h_chi_YL->Fill(chi2Ndf[m]);
            FitCanv_Yellow_Lt->Update();
            // sprintf(name,"Minv%i_pTbin%i_etaLt0_YELLOW_9bin.png",mbin, m);
            // c1->Print(name);
            // fitCanvas[mbin]->Print(Form("./FitPlots_Mbin_%i_AsymVspT.pdf", mbin));
        } // 9 bin loop
        if (mbin == 0)
        {
            FitCanv_Yellow_Lt->Print("./TriggerBias/FitPlot_Mbin_Yellow_Lt.pdf(", "pdf");
        }
        else if (mbin == 4)
        {
            FitCanv_Yellow_Lt->Print("./TriggerBias/FitPlot_Mbin_Yellow_Lt.pdf)", "pdf");
        }
        else
        {
            FitCanv_Yellow_Lt->Print("./TriggerBias/FitPlot_Mbin_Yellow_Lt.pdf", "pdf");
        }
        // Yellow beam asymmetry ends
        // fitCanvas[mbin]->Print(Form("./FitPlots_Mbin_%i_AsymVspT.pdf)", mbin));

        // mbin is invarient mass bin where all data set is divided equally in 5 bins, the whole code is in this loop
        if (mbin == 0)
        { // for Blue beam eta > 0

            for (int kkk = 0; kkk < 9; kkk++)

            {
                A_pT1bg[kkk] = Abg[kkk];
                deltaA_pT1bg[kkk] = deltaAbg[kkk];
                A_pT1bl[kkk] = Abl[kkk];
                deltaA_pT1bl[kkk] = deltaAbl[kkk];
                A_pT1yg[kkk] = Ayg[kkk];
                deltaA_pT1yg[kkk] = deltaAyg[kkk];
                A_pT1yl[kkk] = Ayl[kkk];
                deltaA_pT1yl[kkk] = deltaAyl[kkk];
            }

            // double A_pT1bg[9] = {Abg[0], Abg[1], Abg[2], Abg[3], Abg[4], Abg[5], Abg[6], Abg[7], Abg[8]};
            // double deltaA_pT1bg[9] = {deltaAbg[0], deltaAbg[1], deltaAbg[2], deltaAbg[3], deltaAbg[4], deltaAbg[5], deltaAbg[6], deltaAbg[7], deltaAbg[8]};
            //// for Blue beam eta < 0
            // double A_pT1bl[9] = {Abl[0], Abl[1], Abl[2], Abl[3], Abl[4], Abl[5], Abl[6], Abl[7], Abl[8]};
            // double deltaA_pT1bl[9] = {deltaAbl[0], deltaAbl[1], deltaAbl[2], deltaAbl[3], deltaAbl[4], deltaAbl[5], deltaAbl[6], deltaAbl[7], deltaAbl[8]};
            //// for Yellow beam eta >0
            // double A_pT1yg[9] = {Ayg[0], Ayg[1], Ayg[2], Ayg[3], Ayg[4], Ayg[5], Ayg[6], Ayg[7], Ayg[8]};
            // double deltaA_pT1yg[9] = {deltaAyg[0], deltaAyg[1], deltaAyg[2], deltaAyg[3], deltaAyg[4], deltaAyg[5], deltaAyg[6], deltaAyg[7], deltaAyg[8]};
            //
            // double A_pT1yl[9] = {Ayl[0], Ayl[1], Ayl[2], Ayl[3], Ayl[4], Ayl[5], Ayl[6], Ayl[7], Ayl[8]};
            // double deltaA_pT1yl[9] = {deltaAyl[0], deltaAyl[1], deltaAyl[2], deltaAyl[3], deltaAyl[4], deltaAyl[5], deltaAyl[6], deltaAyl[7], deltaAyl[8]};
        }

        if (mbin == 1)
        { // for Blue beam eta > 0
            for (int kkk = 0; kkk < 9; kkk++)

            {
                A_pT2bg[kkk] = Abg[kkk];
                deltaA_pT2bg[kkk] = deltaAbg[kkk];
                A_pT2bl[kkk] = Abl[kkk];
                deltaA_pT2bl[kkk] = deltaAbl[kkk];
                A_pT2yg[kkk] = Ayg[kkk];
                deltaA_pT2yg[kkk] = deltaAyg[kkk];
                A_pT2yl[kkk] = Ayl[kkk];
                deltaA_pT2yl[kkk] = deltaAyl[kkk];
            }

            // double A_pT2bg[9] = {Abg[0], Abg[1], Abg[2], Abg[3], Abg[4], Abg[5], Abg[6], Abg[7], Abg[8]};
            // double deltaA_pT2bg[9] = {deltaAbg[0], deltaAbg[1], deltaAbg[2], deltaAbg[3], deltaAbg[4], deltaAbg[5], deltaAbg[6], deltaAbg[7], deltaAbg[8]};

            // double A_pT2bl[9] = {Abl[0], Abl[1], Abl[2], Abl[3], Abl[4], Abl[5], Abl[6], Abl[7], Abl[8]};
            // double deltaA_pT2bl[9] = {deltaAbl[0], deltaAbl[1], deltaAbl[2], deltaAbl[3], deltaAbl[4], deltaAbl[5], deltaAbl[6], deltaAbl[7], deltaAbl[8]};

            // double A_pT2yg[9] = {Ayg[0], Ayg[1], Ayg[2], Ayg[3], Ayg[4], Ayg[5], Ayg[6], Ayg[7], Ayg[8]};
            // double deltaA_pT2yg[9] = {deltaAyg[0], deltaAyg[1], deltaAyg[2], deltaAyg[3], deltaAyg[4], deltaAyg[5], deltaAyg[6], deltaAyg[7], deltaAyg[8]};

            // double A_pT2yl[9] = {Ayl[0], Ayl[1], Ayl[2], Ayl[3], Ayl[4], Ayl[5], Ayl[6], Ayl[7], Ayl[8]};
            // double deltaA_pT2yl[9] = {deltaAyl[0], deltaAyl[1], deltaAyl[2], deltaAyl[3], deltaAyl[4], deltaAyl[5], deltaAyl[6], deltaAyl[7], deltaAyl[8]};
        }
        if (mbin == 2)
        {

            for (int kkk = 0; kkk < 9; kkk++)

            {
                A_pT3bg[kkk] = Abg[kkk];
                deltaA_pT3bg[kkk] = deltaAbg[kkk];
                A_pT3bl[kkk] = Abl[kkk];
                deltaA_pT3bl[kkk] = deltaAbl[kkk];
                A_pT3yg[kkk] = Ayg[kkk];
                deltaA_pT3yg[kkk] = deltaAyg[kkk];
                A_pT3yl[kkk] = Ayl[kkk];
                deltaA_pT3yl[kkk] = deltaAyl[kkk];
            }

            // double A_pT3bg[9] = {Abg[0], Abg[1], Abg[2], Abg[3], Abg[4], Abg[5], Abg[6], Abg[7], Abg[8]};
            // double deltaA_pT3bg[9] = {deltaAbg[0], deltaAbg[1], deltaAbg[2], deltaAbg[3], deltaAbg[4], deltaAbg[5], deltaAbg[6], deltaAbg[7], deltaAbg[8]};
            //
            // double A_pT3bl[9] = {Abl[0], Abl[1], Abl[2], Abl[3], Abl[4], Abl[5], Abl[6], Abl[7], Abl[8]};
            // double deltaA_pT3bl[9] = {deltaAbl[0], deltaAbl[1], deltaAbl[2], deltaAbl[3], deltaAbl[4], deltaAbl[5], deltaAbl[6], deltaAbl[7], deltaAbl[8]};
            //
            // double A_pT3yg[9] = {Ayg[0], Ayg[1], Ayg[2], Ayg[3], Ayg[4], Ayg[5], Ayg[6], Ayg[7], Ayg[8]};
            // double deltaA_pT3yg[9] = {deltaAyg[0], deltaAyg[1], deltaAyg[2], deltaAyg[3], deltaAyg[4], deltaAyg[5], deltaAyg[6], deltaAyg[7], deltaAyg[8]};
            //
            // double A_pT3yl[9] = {Ayl[0], Ayl[1], Ayl[2], Ayl[3], Ayl[4], Ayl[5], Ayl[6], Ayl[7], Ayl[8]};
            // double deltaA_pT3yl[9] = {deltaAyl[0], deltaAyl[1], deltaAyl[2], deltaAyl[3], deltaAyl[4], deltaAyl[5], deltaAyl[6], deltaAyl[7], deltaAyl[8]};
        }
        if (mbin == 3)
        {
            for (int kkk = 0; kkk < 9; kkk++)

            {
                A_pT4bg[kkk] = Abg[kkk];
                deltaA_pT4bg[kkk] = deltaAbg[kkk];
                A_pT4bl[kkk] = Abl[kkk];
                deltaA_pT4bl[kkk] = deltaAbl[kkk];
                A_pT4yg[kkk] = Ayg[kkk];
                deltaA_pT4yg[kkk] = deltaAyg[kkk];
                A_pT4yl[kkk] = Ayl[kkk];
                deltaA_pT4yl[kkk] = deltaAyl[kkk];
            }

            // double A_pT4bg[9] = {Abg[0], Abg[1], Abg[2], Abg[3], Abg[4], Abg[5], Abg[6], Abg[7], Abg[8]};
            // double deltaA_pT4bg[9] = {deltaAbg[0], deltaAbg[1], deltaAbg[2], deltaAbg[3], deltaAbg[4], deltaAbg[5], deltaAbg[6], deltaAbg[7], deltaAbg[8]};
            //
            // double A_pT4bl[9] = {Abl[0], Abl[1], Abl[2], Abl[3], Abl[4], Abl[5], Abl[6], Abl[7], Abl[8]};
            // double deltaA_pT4bl[9] = {deltaAbl[0], deltaAbl[1], deltaAbl[2], deltaAbl[3], deltaAbl[4], deltaAbl[5], deltaAbl[6], deltaAbl[7], deltaAbl[8]};
            //
            // double A_pT4yg[9] = {Ayg[0], Ayg[1], Ayg[2], Ayg[3], Ayg[4], Ayg[5], Ayg[6], Ayg[7], Ayg[8]};
            // double deltaA_pT4yg[9] = {deltaAyg[0], deltaAyg[1], deltaAyg[2], deltaAyg[3], deltaAyg[4], deltaAyg[5], deltaAyg[6], deltaAyg[7], deltaAyg[8]};
            //
            // double A_pT4yl[9] = {Ayl[0], Ayl[1], Ayl[2], Ayl[3], Ayl[4], Ayl[5], Ayl[6], Ayl[7], Ayl[8]};
            // double deltaA_pT4yl[9] = {deltaAyl[0], deltaAyl[1], deltaAyl[2], deltaAyl[3], deltaAyl[4], deltaAyl[5], deltaAyl[6], deltaAyl[7], deltaAyl[8]};
        }
        if (mbin == 4)
        {

            for (int kkk = 0; kkk < 9; kkk++)

            {
                A_pT5bg[kkk] = Abg[kkk];
                deltaA_pT5bg[kkk] = deltaAbg[kkk];
                A_pT5bl[kkk] = Abl[kkk];
                deltaA_pT5bl[kkk] = deltaAbl[kkk];
                A_pT5yg[kkk] = Ayg[kkk];
                deltaA_pT5yg[kkk] = deltaAyg[kkk];
                A_pT5yl[kkk] = Ayl[kkk];
                deltaA_pT5yl[kkk] = deltaAyl[kkk];
            }

            // double A_pT5bg[9] = {Abg[0], Abg[1], Abg[2], Abg[3], Abg[4], Abg[5], Abg[6], Abg[7], Abg[8]};
            // double deltaA_pT5bg[9] = {deltaAbg[0], deltaAbg[1], deltaAbg[2], deltaAbg[3], deltaAbg[4], deltaAbg[5], deltaAbg[6], deltaAbg[7], deltaAbg[8]};
            //
            // double A_pT5bl[9] = {Abl[0], Abl[1], Abl[2], Abl[3], Abl[4], Abl[5], Abl[6], Abl[7], Abl[8]};
            // double deltaA_pT5bl[9] = {deltaAbl[0], deltaAbl[1], deltaAbl[2], deltaAbl[3], deltaAbl[4], deltaAbl[5], deltaAbl[6], deltaAbl[7], deltaAbl[8]};
            //
            // double A_pT5yg[9] = {Ayg[0], Ayg[1], Ayg[2], Ayg[3], Ayg[4], Ayg[5], Ayg[6], Ayg[7], Ayg[8]};
            // double deltaA_pT5yg[9] = {deltaAyg[0], deltaAyg[1], deltaAyg[2], deltaAyg[3], deltaAyg[4], deltaAyg[5], deltaAyg[6], deltaAyg[7], deltaAyg[8]};
            //
            // double A_pT5yl[9] = {Ayl[0], Ayl[1], Ayl[2], Ayl[3], Ayl[4], Ayl[5], Ayl[6], Ayl[7], Ayl[8]};
            // double deltaA_pT5yl[9] = {deltaAyl[0], deltaAyl[1], deltaAyl[2], deltaAyl[3], deltaAyl[4], deltaAyl[5], deltaAyl[6], deltaAyg[7], deltaAyg[8]};
        }

    } // mbin bin loop

    double MinvBGt[5][9] = {0};
    double MinvBLt[5][9] = {0};
    double MinvYGt[5][9] = {0};
    double MinvYLt[5][9] = {0};
    double MinvAvgGt[5][9] = {0};
    double MinvAvgLt[5][9] = {0};

    double Eta_pairBGt[5][9] = {0};
    double Eta_pairBLt[5][9] = {0};
    double Eta_pairYGt[5][9] = {0};
    double Eta_pairYLt[5][9] = {0};
    double Eta_pairAvgGt[5][9] = {0};
    double Eta_pairAvgLt[5][9] = {0};

    double pTBGt[5][9] = {0};
    double pTBLt[5][9] = {0};
    double pTYGt[5][9] = {0};
    double pTYLt[5][9] = {0};
    double pTAvgGt[5][9] = {0};
    double pTAvgLt[5][9] = {0};

    double avgMinvGt[5] = {0};
    double avgMinvLt[5] = {0};
    double avg_MinvGtB[5] = {0};
    double avg_MinvLtB[5] = {0};
    double avg_MinvGtY[5] = {0};
    double avg_MinvLtY[5] = {0};

    for (Int_t n = 0; n < 5; n++)
    {
        avg_Minv[n] = hMinv[n]->GetMean();
        avg_MinvGtB[n] = hMinvGtB[n]->GetMean();
        avg_MinvLtB[n] = hMinvLtB[n]->GetMean();
        avg_MinvGtY[n] = hMinvGtY[n]->GetMean();
        avg_MinvLtY[n] = hMinvLtY[n]->GetMean();
    }
    for (int pbin = 0; pbin < 5; pbin++)
    {
        avgMinvGt[pbin] = 0.5 * (avg_MinvGtB[pbin] + avg_MinvGtY[pbin]);
        avgMinvLt[pbin] = 0.5 * (avg_MinvLtB[pbin] + avg_MinvLtY[pbin]);
        for (int mbin = 0; mbin < 9; mbin++)
        {
            MinvBGt[pbin][mbin] = MpairsBGt[pbin][mbin] / (double)NpairsBGt[pbin][mbin];
            pTBGt[pbin][mbin] = pTpairsBGt[pbin][mbin] / (double)NpairsBGt[pbin][mbin];
            Eta_pairBGt[pbin][mbin] = EtapairsBGt[pbin][mbin] / (double)NpairsBGt[pbin][mbin];

            MinvYGt[pbin][mbin] = MpairsYGt[pbin][mbin] / (double)NpairsYGt[pbin][mbin];
            pTYGt[pbin][mbin] = pTpairsYGt[pbin][mbin] / (double)NpairsYGt[pbin][mbin];
            Eta_pairYGt[pbin][mbin] = EtapairsYGt[pbin][mbin] / (double)NpairsYGt[pbin][mbin];

            MinvAvgGt[pbin][mbin] = 0.5 * (MinvBGt[pbin][mbin] + MinvYGt[pbin][mbin]);
            pTAvgGt[pbin][mbin] = 0.5 * (pTBGt[pbin][mbin] + pTYGt[pbin][mbin]);
            Eta_pairAvgGt[pbin][mbin] = 0.5 * (Eta_pairBGt[pbin][mbin] + Eta_pairYGt[pbin][mbin]);

            MinvBLt[pbin][mbin] = MpairsBLt[pbin][mbin] / (double)NpairsBLt[pbin][mbin];
            pTBLt[pbin][mbin] = pTpairsBLt[pbin][mbin] / (double)NpairsBLt[pbin][mbin];
            Eta_pairBLt[pbin][mbin] = EtapairsBLt[pbin][mbin] / (double)NpairsBLt[pbin][mbin];

            MinvYLt[pbin][mbin] = MpairsYLt[pbin][mbin] / (double)NpairsYLt[pbin][mbin];
            pTYLt[pbin][mbin] = pTpairsYLt[pbin][mbin] / (double)NpairsYLt[pbin][mbin];
            Eta_pairYLt[pbin][mbin] = EtapairsYLt[pbin][mbin] / (double)NpairsBLt[pbin][mbin];

            MinvAvgLt[pbin][mbin] = 0.5 * (MinvBLt[pbin][mbin] + MinvYLt[pbin][mbin]);
            pTAvgLt[pbin][mbin] = 0.5 * (pTBLt[pbin][mbin] + pTYLt[pbin][mbin]);
            Eta_pairAvgLt[pbin][mbin] = 0.5 * (Eta_pairBLt[pbin][mbin] + Eta_pairYLt[pbin][mbin]);
        }
    }

    Output << "avgMinvGtB[5]={" << avg_MinvGtB[0] << "," << avg_MinvGtB[1] << "," << avg_MinvGtB[2] << "," << avg_MinvGtB[3] << "," << avg_MinvGtB[4] << "}" << endl;
    Output << "avgMinvGtY[5]={" << avg_MinvGtY[0] << "," << avg_MinvGtY[1] << "," << avg_MinvGtY[2] << "," << avg_MinvGtY[3] << "," << avg_MinvGtY[4] << "}" << endl;
    Output << "avgMinvLtB[5]={" << avg_MinvLtB[0] << "," << avg_MinvLtB[1] << "," << avg_MinvLtB[2] << "," << avg_MinvLtB[3] << "," << avg_MinvLtB[4] << "}" << endl;
    Output << "avgMinvLtY[5]={" << avg_MinvLtY[0] << "," << avg_MinvLtY[1] << "," << avg_MinvLtY[2] << "," << avg_MinvLtY[3] << "," << avg_MinvLtY[4] << "}" << endl;
    Output << "avgMinvGt[5]={" << avgMinvGt[0] << "," << avgMinvGt[1] << "," << avgMinvGt[2] << "," << avgMinvGt[3] << "," << avgMinvGt[4] << "}" << endl;
    Output << "avgMinvLt[5]={" << avgMinvLt[0] << "," << avgMinvLt[1] << "," << avgMinvLt[2] << "," << avgMinvLt[3] << "," << avgMinvLt[4] << "}" << endl;

    Output << "MinvAvgGt[0][9]={" << MinvAvgGt[0][0] << "," << MinvAvgGt[0][1] << "," << MinvAvgGt[0][2] << "," << MinvAvgGt[0][3] << "," << MinvAvgGt[0][4] << "," << MinvAvgGt[0][5] << "," << MinvAvgGt[0][6] << "," << MinvAvgGt[0][7] << "," << MinvAvgGt[0][8] << "}" << endl;
    Output << "MinvAvgGt[1][9]={" << MinvAvgGt[1][0] << "," << MinvAvgGt[1][1] << "," << MinvAvgGt[1][2] << "," << MinvAvgGt[1][3] << "," << MinvAvgGt[1][4] << "," << MinvAvgGt[1][5] << "," << MinvAvgGt[1][6] << "," << MinvAvgGt[1][7] << "," << MinvAvgGt[1][8] << "}" << endl;
    Output << "MinvAvgGt[2][9]={" << MinvAvgGt[2][0] << "," << MinvAvgGt[2][1] << "," << MinvAvgGt[2][2] << "," << MinvAvgGt[2][3] << "," << MinvAvgGt[2][4] << "," << MinvAvgGt[2][5] << "," << MinvAvgGt[2][6] << "," << MinvAvgGt[2][7] << "," << MinvAvgGt[2][8] << "}" << endl;
    Output << "MinvAvgGt[3][9]={" << MinvAvgGt[3][0] << "," << MinvAvgGt[3][1] << "," << MinvAvgGt[3][2] << "," << MinvAvgGt[3][3] << "," << MinvAvgGt[3][4] << "," << MinvAvgGt[3][5] << "," << MinvAvgGt[3][6] << "," << MinvAvgGt[3][7] << "," << MinvAvgGt[3][8] << "}" << endl;
    Output << "MinvAvgGt[4][9]={" << MinvAvgGt[4][0] << "," << MinvAvgGt[4][1] << "," << MinvAvgGt[4][2] << "," << MinvAvgGt[4][3] << "," << MinvAvgGt[4][4] << "," << MinvAvgGt[4][5] << "," << MinvAvgGt[4][6] << "," << MinvAvgGt[4][7] << "," << MinvAvgGt[4][8] << "}" << endl;

    Output << "pTAvgGt[0][9]={" << pTAvgGt[0][0] << "," << pTAvgGt[0][1] << "," << pTAvgGt[0][2] << "," << pTAvgGt[0][3] << "," << pTAvgGt[0][4] << "," << pTAvgGt[0][5] << "," << pTAvgGt[0][6] << "," << pTAvgGt[0][7] << "," << pTAvgGt[0][8] << "}" << endl;
    Output << "pTAvgGt[1][9]={" << pTAvgGt[1][0] << "," << pTAvgGt[1][1] << "," << pTAvgGt[1][2] << "," << pTAvgGt[1][3] << "," << pTAvgGt[1][4] << "," << pTAvgGt[1][5] << "," << pTAvgGt[1][6] << "," << pTAvgGt[1][7] << "," << pTAvgGt[1][8] << "}" << endl;
    Output << "pTAvgGt[2][9]={" << pTAvgGt[2][0] << "," << pTAvgGt[2][1] << "," << pTAvgGt[2][2] << "," << pTAvgGt[2][3] << "," << pTAvgGt[2][4] << "," << pTAvgGt[2][5] << "," << pTAvgGt[2][6] << "," << pTAvgGt[2][7] << "," << pTAvgGt[2][8] << "}" << endl;
    Output << "pTAvgGt[3][9]={" << pTAvgGt[3][0] << "," << pTAvgGt[3][1] << "," << pTAvgGt[3][2] << "," << pTAvgGt[3][3] << "," << pTAvgGt[3][4] << "," << pTAvgGt[3][5] << "," << pTAvgGt[3][6] << "," << pTAvgGt[3][7] << "," << pTAvgGt[3][8] << "}" << endl;
    Output << "pTAvgGt[4][9]={" << pTAvgGt[4][0] << "," << pTAvgGt[4][1] << "," << pTAvgGt[4][2] << "," << pTAvgGt[4][3] << "," << pTAvgGt[4][4] << "," << pTAvgGt[4][5] << "," << pTAvgGt[4][6] << "," << pTAvgGt[4][7] << "," << pTAvgGt[4][8] << "}" << endl;

    //===================================
    Output << "Eta_pairAvgGt[0][9]={" << Eta_pairAvgGt[0][0] << "," << Eta_pairAvgGt[0][1] << "," << Eta_pairAvgGt[0][2] << "," << Eta_pairAvgGt[0][3] << "," << Eta_pairAvgGt[0][4] << "," << Eta_pairAvgGt[0][5] << "," << Eta_pairAvgGt[0][6] << "," << Eta_pairAvgGt[0][7] << "," << Eta_pairAvgGt[0][8] << "}" << endl;
    Output << "Eta_pairAvgGt[1][9]={" << Eta_pairAvgGt[1][0] << "," << Eta_pairAvgGt[1][1] << "," << Eta_pairAvgGt[1][2] << "," << Eta_pairAvgGt[1][3] << "," << Eta_pairAvgGt[1][4] << "," << Eta_pairAvgGt[1][5] << "," << Eta_pairAvgGt[1][6] << "," << Eta_pairAvgGt[1][7] << "," << Eta_pairAvgGt[1][8] << "}" << endl;
    Output << "Eta_pairAvgGt[2][9]={" << Eta_pairAvgGt[2][0] << "," << Eta_pairAvgGt[2][1] << "," << Eta_pairAvgGt[2][2] << "," << Eta_pairAvgGt[2][3] << "," << Eta_pairAvgGt[2][4] << "," << Eta_pairAvgGt[2][5] << "," << Eta_pairAvgGt[2][6] << "," << Eta_pairAvgGt[2][7] << "," << Eta_pairAvgGt[2][8] << "}" << endl;
    Output << "Eta_pairAvgGt[3][9]={" << Eta_pairAvgGt[3][0] << "," << Eta_pairAvgGt[3][1] << "," << Eta_pairAvgGt[3][2] << "," << Eta_pairAvgGt[3][3] << "," << Eta_pairAvgGt[3][4] << "," << Eta_pairAvgGt[3][5] << "," << Eta_pairAvgGt[3][6] << "," << Eta_pairAvgGt[3][7] << "," << Eta_pairAvgGt[3][8] << "}" << endl;
    Output << "Eta_pairAvgGt[4][9]={" << Eta_pairAvgGt[4][0] << "," << Eta_pairAvgGt[4][1] << "," << Eta_pairAvgGt[4][2] << "," << Eta_pairAvgGt[4][3] << "," << Eta_pairAvgGt[4][4] << "," << Eta_pairAvgGt[4][5] << "," << Eta_pairAvgGt[4][6] << "," << Eta_pairAvgGt[4][7] << "," << Eta_pairAvgGt[4][8] << "}" << endl;
    //===================================
    Output << "MinvAvgLt[0][9]={" << MinvAvgLt[0][0] << "," << MinvAvgLt[0][1] << "," << MinvAvgLt[0][2] << "," << MinvAvgLt[0][3] << "," << MinvAvgLt[0][4] << "," << MinvAvgLt[0][5] << "," << MinvAvgLt[0][6] << "," << MinvAvgLt[0][7] << "," << MinvAvgLt[0][8] << "}" << endl;
    Output << "MinvAvgLt[1][9]={" << MinvAvgLt[1][0] << "," << MinvAvgLt[1][1] << "," << MinvAvgLt[1][2] << "," << MinvAvgLt[1][3] << "," << MinvAvgLt[1][4] << "," << MinvAvgLt[1][5] << "," << MinvAvgLt[1][6] << "," << MinvAvgLt[1][7] << "," << MinvAvgLt[1][8] << "}" << endl;
    Output << "MinvAvgLt[2][9]={" << MinvAvgLt[2][0] << "," << MinvAvgLt[2][1] << "," << MinvAvgLt[2][2] << "," << MinvAvgLt[2][3] << "," << MinvAvgLt[2][4] << "," << MinvAvgLt[2][5] << "," << MinvAvgLt[2][6] << "," << MinvAvgLt[2][7] << "," << MinvAvgLt[2][8] << "}" << endl;
    Output << "MinvAvgLt[3][9]={" << MinvAvgLt[3][0] << "," << MinvAvgLt[3][1] << "," << MinvAvgLt[3][2] << "," << MinvAvgLt[3][3] << "," << MinvAvgLt[3][4] << "," << MinvAvgLt[3][5] << "," << MinvAvgLt[3][6] << "," << MinvAvgLt[3][7] << "," << MinvAvgLt[3][8] << "}" << endl;
    Output << "MinvAvgLt[4][9]={" << MinvAvgLt[4][0] << "," << MinvAvgLt[4][1] << "," << MinvAvgLt[4][2] << "," << MinvAvgLt[4][3] << "," << MinvAvgLt[4][4] << "," << MinvAvgLt[4][5] << "," << MinvAvgLt[4][6] << "," << MinvAvgLt[4][7] << "," << MinvAvgLt[4][8] << "}" << endl;

    Output << "pTAvgLt[0][9]={" << pTAvgLt[0][0] << "," << pTAvgLt[0][1] << "," << pTAvgLt[0][2] << "," << pTAvgLt[0][3] << "," << pTAvgLt[0][4] << "," << pTAvgLt[0][5] << "," << pTAvgLt[0][6] << "," << pTAvgLt[0][7] << "," << pTAvgLt[0][8] << "}" << endl;
    Output << "pTAvgLt[1][9]={" << pTAvgLt[1][0] << "," << pTAvgLt[1][1] << "," << pTAvgLt[1][2] << "," << pTAvgLt[1][3] << "," << pTAvgLt[1][4] << "," << pTAvgLt[1][5] << "," << pTAvgLt[1][6] << "," << pTAvgLt[1][7] << "," << pTAvgLt[1][8] << "}" << endl;
    Output << "pTAvgLt[2][9]={" << pTAvgLt[2][0] << "," << pTAvgLt[2][1] << "," << pTAvgLt[2][2] << "," << pTAvgLt[2][3] << "," << pTAvgLt[2][4] << "," << pTAvgLt[2][5] << "," << pTAvgLt[2][6] << "," << pTAvgLt[2][7] << "," << pTAvgLt[2][8] << "}" << endl;
    Output << "pTAvgLt[3][9]={" << pTAvgLt[3][0] << "," << pTAvgLt[3][1] << "," << pTAvgLt[3][2] << "," << pTAvgLt[3][3] << "," << pTAvgLt[3][4] << "," << pTAvgLt[3][5] << "," << pTAvgLt[3][6] << "," << pTAvgLt[3][7] << "," << pTAvgLt[3][8] << "}" << endl;
    Output << "pTAvgLt[4][9]={" << pTAvgLt[4][0] << "," << pTAvgLt[4][1] << "," << pTAvgLt[4][2] << "," << pTAvgLt[4][3] << "," << pTAvgLt[4][4] << "," << pTAvgLt[4][5] << "," << pTAvgLt[4][6] << "," << pTAvgLt[4][7] << "," << pTAvgLt[4][8] << "}" << endl;

    //===================================
    Output << "Eta_pairAvgLt[0][9]={" << Eta_pairAvgLt[0][0] << "," << Eta_pairAvgLt[0][1] << "," << Eta_pairAvgLt[0][2] << "," << Eta_pairAvgLt[0][3] << "," << Eta_pairAvgLt[0][4] << "," << Eta_pairAvgLt[0][5] << "," << Eta_pairAvgLt[0][6] << "," << Eta_pairAvgLt[0][7] << "," << Eta_pairAvgLt[0][8] << "}" << endl;
    Output << "Eta_pairAvgLt[1][9]={" << Eta_pairAvgLt[1][0] << "," << Eta_pairAvgLt[1][1] << "," << Eta_pairAvgLt[1][2] << "," << Eta_pairAvgLt[1][3] << "," << Eta_pairAvgLt[1][4] << "," << Eta_pairAvgLt[1][5] << "," << Eta_pairAvgLt[1][6] << "," << Eta_pairAvgLt[1][7] << "," << Eta_pairAvgLt[1][8] << "}" << endl;
    Output << "Eta_pairAvgLt[2][9]={" << Eta_pairAvgLt[2][0] << "," << Eta_pairAvgLt[2][1] << "," << Eta_pairAvgLt[2][2] << "," << Eta_pairAvgLt[2][3] << "," << Eta_pairAvgLt[2][4] << "," << Eta_pairAvgLt[2][5] << "," << Eta_pairAvgLt[2][6] << "," << Eta_pairAvgLt[2][7] << "," << Eta_pairAvgLt[2][8] << "}" << endl;
    Output << "Eta_pairAvgLt[3][9]={" << Eta_pairAvgLt[3][0] << "," << Eta_pairAvgLt[3][1] << "," << Eta_pairAvgLt[3][2] << "," << Eta_pairAvgLt[3][3] << "," << Eta_pairAvgLt[3][4] << "," << Eta_pairAvgLt[3][5] << "," << Eta_pairAvgLt[3][6] << "," << Eta_pairAvgLt[3][7] << "," << Eta_pairAvgLt[3][8] << "}" << endl;
    Output << "Eta_pairAvgLt[4][9]={" << Eta_pairAvgLt[4][0] << "," << Eta_pairAvgLt[4][1] << "," << Eta_pairAvgLt[4][2] << "," << Eta_pairAvgLt[4][3] << "," << Eta_pairAvgLt[4][4] << "," << Eta_pairAvgLt[4][5] << "," << Eta_pairAvgLt[4][6] << "," << Eta_pairAvgLt[4][7] << "," << Eta_pairAvgLt[4][8] << "}" << endl;

    //===================================

    // calculate average asymmetry
    double avgA_pT1g[9] = {0};
    double avgA_pT1l[9] = {0};
    double avgA_pT2g[9] = {0};
    double avgA_pT2l[9] = {0};
    double avgA_pT3g[9] = {0};
    double avgA_pT3l[9] = {0};
    double avgA_pT4g[9] = {0};
    double avgA_pT4l[9] = {0};
    double avgA_pT5g[9] = {0};
    double avgA_pT5l[9] = {0};
    double errA_pT1g[9] = {0};
    double errA_pT1l[9] = {0};
    double errA_pT2g[9] = {0};
    double errA_pT2l[9] = {0};
    double errA_pT3g[9] = {0};
    double errA_pT3l[9] = {0};
    double errA_pT4g[9] = {0};
    double errA_pT4l[9] = {0};
    double errA_pT5g[9] = {0};
    double errA_pT5l[9] = {0};

    double WavgA_pT1g[9] = {0};
    double WerrA_pT1g[9] = {0};
    double WavgA_pT1l[9] = {0};
    double WerrA_pT1l[9] = {0};
    double WavgA_pT2g[9] = {0};
    double WerrA_pT2g[9] = {0};
    double WavgA_pT2l[9] = {0};
    double WerrA_pT2l[9] = {0};
    double WavgA_pT3g[9] = {0};
    double WerrA_pT3g[9] = {0};
    double WavgA_pT3l[9] = {0};
    double WerrA_pT3l[9] = {0};
    double WavgA_pT4g[9] = {0};
    double WerrA_pT4g[9] = {0};
    double WavgA_pT4l[9] = {0};
    double WerrA_pT4l[9] = {0};
    double WavgA_pT5g[9] = {0};
    double WerrA_pT5g[9] = {0};
    double WavgA_pT5l[9] = {0};
    double WerrA_pT5l[9] = {0};

    for (Int_t ii = 0; ii < 9; ii++)
    {
        // Simple Average
        avgA_pT1g[ii] = (A_pT1bg[ii] + A_pT1yg[ii]) / 2.;
        avgA_pT1l[ii] = (A_pT1bl[ii] + A_pT1yl[ii]) / 2.;
        avgA_pT2g[ii] = (A_pT2bg[ii] + A_pT2yg[ii]) / 2.;
        avgA_pT2l[ii] = (A_pT2bl[ii] + A_pT2yl[ii]) / 2.;
        avgA_pT3g[ii] = (A_pT3bg[ii] + A_pT3yg[ii]) / 2.;
        avgA_pT3l[ii] = (A_pT3bl[ii] + A_pT3yl[ii]) / 2.;
        avgA_pT4g[ii] = (A_pT4bg[ii] + A_pT4yg[ii]) / 2.;
        avgA_pT4l[ii] = (A_pT4bl[ii] + A_pT4yl[ii]) / 2.;
        avgA_pT5g[ii] = (A_pT5bg[ii] + A_pT5yg[ii]) / 2.;
        avgA_pT5l[ii] = (A_pT5bl[ii] + A_pT5yl[ii]) / 2.;

        errA_pT1g[ii] = .5 * sqrt(pow(deltaA_pT1bg[ii], 2) + pow(deltaA_pT1yg[ii], 2)); // total asym err, pT bin 1, eta > 0
        errA_pT1l[ii] = .5 * sqrt(pow(deltaA_pT1bl[ii], 2) + pow(deltaA_pT1yl[ii], 2)); // total asym err, pT bin 1, eta < 0
        errA_pT2g[ii] = .5 * sqrt(pow(deltaA_pT2bg[ii], 2) + pow(deltaA_pT2yg[ii], 2));
        errA_pT2l[ii] = .5 * sqrt(pow(deltaA_pT2bl[ii], 2) + pow(deltaA_pT2yl[ii], 2));
        errA_pT3g[ii] = .5 * sqrt(pow(deltaA_pT3bg[ii], 2) + pow(deltaA_pT3yg[ii], 2));
        errA_pT3l[ii] = .5 * sqrt(pow(deltaA_pT3bl[ii], 2) + pow(deltaA_pT3yl[ii], 2));
        errA_pT4g[ii] = .5 * sqrt(pow(deltaA_pT4bg[ii], 2) + pow(deltaA_pT4yg[ii], 2));
        errA_pT4l[ii] = .5 * sqrt(pow(deltaA_pT4bl[ii], 2) + pow(deltaA_pT4yl[ii], 2));
        errA_pT5g[ii] = .5 * sqrt(pow(deltaA_pT5bg[ii], 2) + pow(deltaA_pT5yg[ii], 2));
        errA_pT5l[ii] = .5 * sqrt(pow(deltaA_pT5bl[ii], 2) + pow(deltaA_pT5yl[ii], 2));

        // Weighted Average
        WavgA_pT1g[ii] = (A_pT1bg[ii] * (1 / pow(deltaA_pT1bg[ii], 2)) + A_pT1yg[ii] * (1 / pow(deltaA_pT1yg[ii], 2))) / ((1 / pow(deltaA_pT1bg[ii], 2)) + (1 / pow(deltaA_pT1yg[ii], 2)));
        WavgA_pT2g[ii] = (A_pT2bg[ii] * (1 / pow(deltaA_pT2bg[ii], 2)) + A_pT2yg[ii] * (1 / pow(deltaA_pT2yg[ii], 2))) / ((1 / pow(deltaA_pT2bg[ii], 2)) + (1 / pow(deltaA_pT2yg[ii], 2)));
        WavgA_pT3g[ii] = (A_pT3bg[ii] * (1 / pow(deltaA_pT3bg[ii], 2)) + A_pT3yg[ii] * (1 / pow(deltaA_pT3yg[ii], 2))) / ((1 / pow(deltaA_pT3bg[ii], 2)) + (1 / pow(deltaA_pT3yg[ii], 2)));
        WavgA_pT4g[ii] = (A_pT4bg[ii] * (1 / pow(deltaA_pT4bg[ii], 2)) + A_pT4yg[ii] * (1 / pow(deltaA_pT4yg[ii], 2))) / ((1 / pow(deltaA_pT4bg[ii], 2)) + (1 / pow(deltaA_pT4yg[ii], 2)));
        WavgA_pT5g[ii] = (A_pT5bg[ii] * (1 / pow(deltaA_pT5bg[ii], 2)) + A_pT5yg[ii] * (1 / pow(deltaA_pT5yg[ii], 2))) / ((1 / pow(deltaA_pT5bg[ii], 2)) + (1 / pow(deltaA_pT5yg[ii], 2)));

        WerrA_pT1g[ii] = 1 / sqrt((1 / pow(deltaA_pT1bg[ii], 2)) + (1 / pow(deltaA_pT1yg[ii], 2))); // total asym err, pT bin 1, eta > 0
        WerrA_pT2g[ii] = 1 / sqrt((1 / pow(deltaA_pT2bg[ii], 2)) + (1 / pow(deltaA_pT2yg[ii], 2)));
        WerrA_pT3g[ii] = 1 / sqrt((1 / pow(deltaA_pT3bg[ii], 2)) + (1 / pow(deltaA_pT3yg[ii], 2)));
        WerrA_pT4g[ii] = 1 / sqrt((1 / pow(deltaA_pT4bg[ii], 2)) + (1 / pow(deltaA_pT4yg[ii], 2)));
        WerrA_pT5g[ii] = 1 / sqrt((1 / pow(deltaA_pT5bg[ii], 2)) + (1 / pow(deltaA_pT5yg[ii], 2)));

        WavgA_pT1l[ii] = (A_pT1bl[ii] * (1 / pow(deltaA_pT1bl[ii], 2)) + A_pT1yl[ii] * (1 / pow(deltaA_pT1yl[ii], 2))) / ((1 / pow(deltaA_pT1bl[ii], 2)) + (1 / pow(deltaA_pT1yl[ii], 2)));
        WavgA_pT2l[ii] = (A_pT2bl[ii] * (1 / pow(deltaA_pT2bl[ii], 2)) + A_pT2yl[ii] * (1 / pow(deltaA_pT2yl[ii], 2))) / ((1 / pow(deltaA_pT2bl[ii], 2)) + (1 / pow(deltaA_pT2yl[ii], 2)));
        WavgA_pT3l[ii] = (A_pT3bl[ii] * (1 / pow(deltaA_pT3bl[ii], 2)) + A_pT3yl[ii] * (1 / pow(deltaA_pT3yl[ii], 2))) / ((1 / pow(deltaA_pT3bl[ii], 2)) + (1 / pow(deltaA_pT3yl[ii], 2)));
        WavgA_pT4l[ii] = (A_pT4bl[ii] * (1 / pow(deltaA_pT4bl[ii], 2)) + A_pT4yl[ii] * (1 / pow(deltaA_pT4yl[ii], 2))) / ((1 / pow(deltaA_pT4bl[ii], 2)) + (1 / pow(deltaA_pT4yl[ii], 2)));
        WavgA_pT5l[ii] = (A_pT5bl[ii] * (1 / pow(deltaA_pT5bl[ii], 2)) + A_pT5yl[ii] * (1 / pow(deltaA_pT5yl[ii], 2))) / ((1 / pow(deltaA_pT5bl[ii], 2)) + (1 / pow(deltaA_pT5yl[ii], 2)));

        WerrA_pT1l[ii] = 1 / sqrt((1 / pow(deltaA_pT1bl[ii], 2)) + (1 / pow(deltaA_pT1yl[ii], 2))); // total asym err, pT bin 1, eta < 0
        WerrA_pT2l[ii] = 1 / sqrt((1 / pow(deltaA_pT2bl[ii], 2)) + (1 / pow(deltaA_pT2yl[ii], 2)));
        WerrA_pT3l[ii] = 1 / sqrt((1 / pow(deltaA_pT3bl[ii], 2)) + (1 / pow(deltaA_pT3yl[ii], 2)));
        WerrA_pT4l[ii] = 1 / sqrt((1 / pow(deltaA_pT4bl[ii], 2)) + (1 / pow(deltaA_pT4yl[ii], 2)));
        WerrA_pT5l[ii] = 1 / sqrt((1 / pow(deltaA_pT5bl[ii], 2)) + (1 / pow(deltaA_pT5yl[ii], 2)));

        cout << "Blue and Yellow beam asymmetry and its error" << endl;
        cout << A_pT1bg[ii] << "\t" << deltaA_pT1bg[ii] << "<==Blue" << A_pT1yg[ii] << "\t" << deltaA_pT1yg[ii] << "<==Yellow" << endl;
        cout << "Simple Ave and its error" << endl;
        cout << avgA_pT1g[ii] << "\t" << errA_pT1g[ii] << endl;
        cout << "Weighted Ave and its error" << endl;
        cout << WavgA_pT1g[ii] << "\t" << WerrA_pT1g[ii] << endl;
        cout << "   " << endl;
        cout << "   " << endl;
        cout << "   " << endl;
    }
    // total asymmetry error calculation
    // since two asymmetry(BLUE+YELLOW) are averaged, error propagation formaula is used for total asymmetry error
    // if function, f = (a+b)/2, Then Error, df =1/2 √{(∂f/∂a)^2*(∆a)^2+(∂f/∂b)^2*(∆b)^2}

    Output << "                                                                                                  " << endl;
    Output << "                                                                                                  " << endl;
    Output << "                                                                                                  " << endl;
    Output << "                                      BLUE BEAM                                                   " << endl;
    Output << "********************   FINAL ASYMMETRY , CONE > 0.7, Avg pT Binning for BLUE BEAM Only****************" << endl;

    Output << "<<<<<<<<<<<<<<<   pT Bin Averages for asymemtry plot  >>>>>>>>>>>>>>" << endl;
    Output << "Minv-1, pTBGt[0][9]= {" << pTBGt[0][0] << ", " << pTBGt[0][1] << ", " << pTBGt[0][2] << ", " << pTBGt[0][3] << ", " << pTBGt[0][4] << ", " << pTBGt[0][5] << ", " << pTBGt[0][6] << ", " << pTBGt[0][7] << ", " << pTBGt[0][8] << "}" << endl;
    Output << "Minv-2, pTBGt[1][9]= {" << pTBGt[1][0] << ", " << pTBGt[1][1] << ", " << pTBGt[1][2] << ", " << pTBGt[1][3] << ", " << pTBGt[1][4] << ", " << pTBGt[1][5] << ", " << pTBGt[1][6] << ", " << pTBGt[1][7] << ", " << pTBGt[1][8] << "}" << endl;
    Output << "Minv-3, pTBGt[2][9]= {" << pTBGt[2][0] << ", " << pTBGt[2][1] << ", " << pTBGt[2][2] << ", " << pTBGt[2][3] << ", " << pTBGt[2][4] << ", " << pTBGt[2][5] << ", " << pTBGt[2][6] << ", " << pTBGt[2][7] << ", " << pTBGt[2][8] << "}" << endl;
    Output << "Minv-4, pTBGt[3][9]= {" << pTBGt[3][0] << ", " << pTBGt[3][1] << ", " << pTBGt[3][2] << ", " << pTBGt[3][3] << ", " << pTBGt[3][4] << ", " << pTBGt[3][5] << ", " << pTBGt[3][6] << ", " << pTBGt[3][7] << ", " << pTBGt[3][8] << "}" << endl;
    Output << "Minv-5, pTBGt[4][9]= {" << pTBGt[4][0] << ", " << pTBGt[4][1] << ", " << pTBGt[4][2] << ", " << pTBGt[4][3] << ", " << pTBGt[4][4] << ", " << pTBGt[4][5] << ", " << pTBGt[4][6] << ", " << pTBGt[4][7] << ", " << pTBGt[4][8] << "}" << endl;

    Output << "Minv-1, pTBLt[0][9]= {" << pTBLt[0][0] << ", " << pTBLt[0][1] << ", " << pTBLt[0][2] << ", " << pTBLt[0][3] << ", " << pTBLt[0][4] << ", " << pTBLt[0][5] << ", " << pTBLt[0][6] << ", " << pTBLt[0][7] << ", " << pTBLt[0][8] << "}" << endl;
    Output << "Minv-2, pTBLt[1][9]= {" << pTBLt[1][0] << ", " << pTBLt[1][1] << ", " << pTBLt[1][2] << ", " << pTBLt[1][3] << ", " << pTBLt[1][4] << ", " << pTBLt[1][5] << ", " << pTBLt[1][6] << ", " << pTBLt[1][7] << ", " << pTBLt[1][8] << "}" << endl;
    Output << "Minv-3, pTBLt[2][]9= {" << pTBLt[2][0] << ", " << pTBLt[2][1] << ", " << pTBLt[2][2] << ", " << pTBLt[2][3] << ", " << pTBLt[2][4] << ", " << pTBLt[2][5] << ", " << pTBLt[2][6] << ", " << pTBLt[2][7] << ", " << pTBLt[2][8] << "}" << endl;
    Output << "Minv-4, pTBLt[3][9]= {" << pTBLt[3][0] << ", " << pTBLt[3][1] << ", " << pTBLt[3][2] << ", " << pTBLt[3][3] << ", " << pTBLt[3][4] << ", " << pTBLt[3][5] << ", " << pTBLt[3][6] << ", " << pTBLt[3][7] << ", " << pTBLt[3][8] << "}" << endl;
    Output << "Minv-5, pTBLt[4][9]= {" << pTBLt[4][0] << ", " << pTBLt[4][1] << ", " << pTBLt[4][2] << ", " << pTBLt[4][3] << ", " << pTBLt[4][4] << ", " << pTBLt[4][5] << ", " << pTBLt[4][6] << ", " << pTBLt[4][7] << ", " << pTBLt[4][8] << "}" << endl;

    Output << "<<<<<<<< Average Polarization Values >>>>>>>>" << endl;
    Output << "BLUE <P>: " << avgPolB << ", Err_P: " << rmsB << endl;
    Output << "YELLOW <P>: " << avgPolY << ", Err_P: " << rmsY << endl;

    Output << "                                                                                                  " << endl;
    Output << "                                                                                                  " << endl;
    Output << "<<<<<<<<<<<       Eta < 0, BLUE BEAM A_{UT}   >>>>>>>>>>>>>>>>>" << endl;
    Output << "Minv-1,<A_{UT}> ={" << A_pT1bl[0] << ", " << A_pT1bl[1] << ", " << A_pT1bl[2] << ", " << A_pT1bl[3] << ", " << A_pT1bl[4] << ", " << A_pT1bl[5] << ", " << A_pT1bl[6] << ", " << A_pT1bl[7] << ", " << A_pT1bl[8] << "}" << endl;
    Output << "Minv-1,<Err>={" << deltaA_pT1bl[0] << ", " << deltaA_pT1bl[1] << ", " << deltaA_pT1bl[2] << ", " << deltaA_pT1bl[3] << ", " << deltaA_pT1bl[4] << ", " << deltaA_pT1bl[5] << ", " << deltaA_pT1bl[6] << ", " << deltaA_pT1bl[7] << ", " << deltaA_pT1bl[8] << "}" << endl;
    Output << "Minv-2,<A_{UT}>={" << A_pT2bl[0] << ", " << A_pT2bl[1] << ", " << A_pT2bl[2] << ", " << A_pT2bl[3] << ", " << A_pT2bl[4] << ", " << A_pT2bl[5] << ", " << A_pT2bl[6] << ", " << A_pT2bl[7] << ", " << A_pT2bl[8] << "}" << endl;
    Output << "Minv-2,<Err> ={" << deltaA_pT2bl[0] << ", " << deltaA_pT2bl[1] << ", " << deltaA_pT2bl[2] << ", " << deltaA_pT2bl[3] << ", " << deltaA_pT2bl[4] << ", " << deltaA_pT2bl[5] << ", " << deltaA_pT2bl[6] << ", " << deltaA_pT2bl[7] << ", " << deltaA_pT2bl[8] << "}" << endl;
    Output << "Minv-3,<A_{UT}> ={" << A_pT3bl[0] << ", " << A_pT3bl[1] << ", " << A_pT3bl[2] << ", " << A_pT3bl[3] << ", " << A_pT3bl[4] << ", " << A_pT3bl[5] << ", " << A_pT3bl[6] << ", " << A_pT3bl[7] << ", " << A_pT3bl[8] << "}" << endl;
    Output << "Minv-3,<Err> ={" << deltaA_pT3bl[0] << ", " << deltaA_pT3bl[1] << ", " << deltaA_pT3bl[2] << ", " << deltaA_pT3bl[3] << ", " << deltaA_pT3bl[4] << ", " << deltaA_pT3bl[5] << ", " << deltaA_pT3bl[6] << ", " << deltaA_pT3bl[7] << ", " << deltaA_pT3bl[8] << "}" << endl;
    Output << "Minv-4,<A_{UT}> ={" << A_pT4bl[0] << ", " << A_pT4bl[1] << ", " << A_pT4bl[2] << ", " << A_pT4bl[3] << ", " << A_pT4bl[4] << ", " << A_pT4bl[5] << ", " << A_pT4bl[6] << ", " << A_pT4bl[7] << ", " << A_pT4bl[8] << "}" << endl;
    Output << "Minv-4,<Err> ={" << deltaA_pT4bl[0] << ", " << deltaA_pT4bl[1] << ", " << deltaA_pT4bl[2] << ", " << deltaA_pT4bl[3] << ", " << deltaA_pT4bl[4] << ", " << deltaA_pT4bl[5] << ", " << deltaA_pT4bl[6] << ", " << deltaA_pT4bl[7] << ", " << deltaA_pT4bl[8] << "}" << endl;
    Output << "Minv-5,<A_{UT}> ={" << A_pT5bl[0] << ", " << A_pT5bl[1] << ", " << A_pT5bl[2] << ", " << A_pT5bl[3] << ", " << A_pT5bl[4] << ", " << A_pT5bl[5] << ", " << A_pT5bl[6] << ", " << A_pT5bl[7] << ", " << A_pT5bl[8] << "}" << endl;
    Output << "Minv-5,<Err> ={" << deltaA_pT5bl[0] << ", " << deltaA_pT5bl[1] << ", " << deltaA_pT5bl[2] << ", " << deltaA_pT5bl[3] << ", " << deltaA_pT5bl[4] << ", " << deltaA_pT5bl[5] << ", " << deltaA_pT5bl[6] << ", " << deltaA_pT5bl[7] << ", " << deltaA_pT5bl[8] << "}" << endl;

    Output << "                                                                                                  " << endl;
    Output << "                                                                                                  " << endl;
    Output << "<<<<<<<<<<<       Eta > 0, BLUE BEAM A_{UT}   >>>>>>>>>>>>>>>>>" << endl;
    Output << "Minv-1,<A_{UT}> ={" << A_pT1bg[0] << ", " << A_pT1bg[1] << ", " << A_pT1bg[2] << ", " << A_pT1bg[3] << ", " << A_pT1bg[4] << ", " << A_pT1bg[5] << ", " << A_pT1bg[6] << ", " << A_pT1bg[7] << ", " << A_pT1bg[8] << "}" << endl;
    Output << "Minv-1,<Err> ={" << deltaA_pT1bg[0] << ", " << deltaA_pT1bg[1] << ", " << deltaA_pT1bg[2] << ", " << deltaA_pT1bg[3] << ", " << deltaA_pT1bg[4] << ", " << deltaA_pT1bg[5] << ", " << deltaA_pT1bg[6] << ", " << deltaA_pT1bg[7] << ", " << deltaA_pT1bg[8] << "}" << endl;
    Output << "Minv-2,<A_{UT}> ={" << A_pT2bg[0] << ", " << A_pT2bg[1] << ", " << A_pT2bg[2] << ", " << A_pT2bg[3] << ", " << A_pT2bg[4] << ", " << A_pT2bg[5] << ", " << A_pT2bg[6] << ", " << A_pT2bg[7] << ", " << A_pT2bg[8] << "}" << endl;
    Output << "Minv-2,<Err> ={" << deltaA_pT2bg[0] << ", " << deltaA_pT2bg[1] << ", " << deltaA_pT2bg[2] << ", " << deltaA_pT2bg[3] << ", " << deltaA_pT2bg[4] << ", " << deltaA_pT2bg[5] << ", " << deltaA_pT2bg[6] << ", " << deltaA_pT2bg[7] << ", " << deltaA_pT2bg[8] << "}" << endl;
    Output << "Minv-3,<A_{UT}> ={" << A_pT3bg[0] << ", " << A_pT3bg[1] << ", " << A_pT3bg[2] << ", " << A_pT3bg[3] << ", " << A_pT3bg[4] << ", " << A_pT3bg[5] << ", " << A_pT3bg[6] << ", " << A_pT3bg[7] << ", " << A_pT3bg[8] << "}" << endl;
    Output << "Minv-3,<Err> ={" << deltaA_pT3bg[0] << ", " << deltaA_pT3bg[1] << ", " << deltaA_pT3bg[2] << ", " << deltaA_pT3bg[3] << ", " << deltaA_pT3bg[4] << ", " << deltaA_pT3bg[5] << ", " << deltaA_pT3bg[6] << ", " << deltaA_pT3bg[7] << ", " << deltaA_pT3bg[8] << "}" << endl;
    Output << "Minv-4,<A_{UT}> ={" << A_pT4bg[0] << ", " << A_pT4bg[1] << ", " << A_pT4bg[2] << ", " << A_pT4bg[3] << ", " << A_pT4bg[4] << ", " << A_pT4bg[5] << ", " << A_pT4bg[6] << ", " << A_pT4bg[7] << ", " << A_pT4bg[8] << "}" << endl;
    Output << "Minv-4,<Err> ={" << deltaA_pT4bg[0] << ", " << deltaA_pT4bg[1] << ", " << deltaA_pT4bg[2] << ", " << deltaA_pT4bg[3] << ", " << deltaA_pT4bg[4] << ", " << deltaA_pT4bg[5] << ", " << deltaA_pT4bg[6] << ", " << deltaA_pT4bg[7] << ", " << deltaA_pT4bg[8] << "}" << endl;
    Output << "Minv-5,<A_{UT}> ={" << A_pT5bg[0] << ", " << A_pT5bg[1] << ", " << A_pT5bg[2] << ", " << A_pT5bg[3] << ", " << A_pT5bg[4] << ", " << A_pT5bg[5] << ", " << A_pT5bg[6] << ", " << A_pT5bg[7] << ", " << A_pT5bg[8] << "}" << endl;
    Output << "Minv-5,<Err> ={" << deltaA_pT5bg[0] << ", " << deltaA_pT5bg[1] << ", " << deltaA_pT5bg[2] << ", " << deltaA_pT5bg[3] << ", " << deltaA_pT5bg[4] << ", " << deltaA_pT5bg[5] << ", " << deltaA_pT5bg[6] << ", " << deltaA_pT5bg[7] << ", " << deltaA_pT5bg[8] << "}" << endl;

    Output << "                                                                                                  " << endl;
    Output << "                                                                                                  " << endl;
    Output << "                                                                                                  " << endl;
    Output << "                                      YELLOW BEAM                                                 " << endl;
    Output << "********************   FINAL ASYMMETRY , CONE > 0.7, pT Binning for YELLOW BEAM Only****************" << endl;

    Output << "<<<<<<<<<<<<<<<   pT Bin Averages for asymemtry plot  >>>>>>>>>>>>>>" << endl;
    Output << "Minv-1, pTYGt[0][9]= {" << pTYGt[0][0] << ", " << pTYGt[0][1] << ", " << pTYGt[0][2] << ", " << pTYGt[0][3] << ", " << pTYGt[0][4] << ", " << pTYGt[0][5] << ", " << pTYGt[0][6] << ", " << pTYGt[0][7] << ", " << pTYGt[0][8] << "}" << endl;
    Output << "Minv-2, pTYGt[1][9]= {" << pTYGt[1][0] << ", " << pTYGt[1][1] << ", " << pTYGt[1][2] << ", " << pTYGt[1][3] << ", " << pTYGt[1][4] << ", " << pTYGt[1][5] << ", " << pTYGt[1][6] << ", " << pTYGt[1][7] << ", " << pTYGt[1][8] << "}" << endl;
    Output << "Minv-3, pTYGt[2][]9= {" << pTYGt[2][0] << ", " << pTYGt[2][1] << ", " << pTYGt[2][2] << ", " << pTYGt[2][3] << ", " << pTYGt[2][4] << ", " << pTYGt[2][5] << ", " << pTYGt[2][6] << ", " << pTYGt[2][7] << ", " << pTYGt[2][8] << "}" << endl;
    Output << "Minv-4, pTYGt[3][9]= {" << pTYGt[3][0] << ", " << pTYGt[3][1] << ", " << pTYGt[3][2] << ", " << pTYGt[3][3] << ", " << pTYGt[3][4] << ", " << pTYGt[3][5] << ", " << pTYGt[3][6] << ", " << pTYGt[3][7] << ", " << pTYGt[3][8] << "}" << endl;
    Output << "Minv-5, pTYGt[4][9]= {" << pTYGt[4][0] << ", " << pTYGt[4][1] << ", " << pTYGt[4][2] << ", " << pTYGt[4][3] << ", " << pTYGt[4][4] << ", " << pTYGt[4][5] << ", " << pTYGt[4][6] << ", " << pTYGt[4][7] << ", " << pTYGt[4][8] << "}" << endl;

    Output << "Minv-1, pTYLt[0][9]= {" << pTYLt[0][0] << ", " << pTYLt[0][1] << ", " << pTYLt[0][2] << ", " << pTYLt[0][3] << ", " << pTYLt[0][4] << ", " << pTYLt[0][5] << ", " << pTYLt[0][6] << ", " << pTYLt[0][7] << ", " << pTYLt[0][8] << "}" << endl;
    Output << "Minv-2, pTYLt[1][9]= {" << pTYLt[1][0] << ", " << pTYLt[1][1] << ", " << pTYLt[1][2] << ", " << pTYLt[1][3] << ", " << pTYLt[1][4] << ", " << pTYLt[1][5] << ", " << pTYLt[1][6] << ", " << pTYLt[1][7] << ", " << pTYLt[1][8] << "}" << endl;
    Output << "Minv-3, pTYLt[2][]9= {" << pTYLt[2][0] << ", " << pTYLt[2][1] << ", " << pTYLt[2][2] << ", " << pTYLt[2][3] << ", " << pTYLt[2][4] << ", " << pTYLt[2][5] << ", " << pTYLt[2][6] << ", " << pTYLt[2][7] << ", " << pTYLt[2][8] << "}" << endl;
    Output << "Minv-4, pTYLt[3][9]= {" << pTYLt[3][0] << ", " << pTYLt[3][1] << ", " << pTYLt[3][2] << ", " << pTYLt[3][3] << ", " << pTYLt[3][4] << ", " << pTYLt[3][5] << ", " << pTYLt[3][6] << ", " << pTYLt[3][7] << ", " << pTYLt[3][8] << "}" << endl;
    Output << "Minv-5, pTYLt[4][9]= {" << pTYLt[4][0] << ", " << pTYLt[4][1] << ", " << pTYLt[4][2] << ", " << pTYLt[4][3] << ", " << pTYLt[4][4] << ", " << pTYLt[4][5] << ", " << pTYLt[4][6] << ", " << pTYLt[4][7] << ", " << pTYLt[4][8] << "}" << endl;

    Output << "                                                                                                  " << endl;
    Output << "                                                                                                  " << endl;
    Output << "<<<<<<<<<<<       Eta < 0, YELLOW BEAM A_{UT}   >>>>>>>>>>>>>>>>>" << endl;
    Output << "Minv-1,<A_{UT}> ={" << A_pT1yl[0] << ", " << A_pT1yl[1] << ", " << A_pT1yl[2] << ", " << A_pT1yl[3] << ", " << A_pT1yl[4] << ", " << A_pT1yl[5] << ", " << A_pT1yl[6] << ", " << A_pT1yl[7] << ", " << A_pT1yl[8] << "}" << endl;
    Output << "Minv-1,<Err>={" << deltaA_pT1yl[0] << ", " << deltaA_pT1yl[1] << ", " << deltaA_pT1yl[2] << ", " << deltaA_pT1yl[3] << ", " << deltaA_pT1yl[4] << ", " << deltaA_pT1yl[5] << ", " << deltaA_pT1yl[6] << ", " << deltaA_pT1yl[7] << ", " << deltaA_pT1yl[8] << "}" << endl;
    Output << "Minv-2,<A_{UT}>={" << A_pT2yl[0] << ", " << A_pT2yl[1] << ", " << A_pT2yl[2] << ", " << A_pT2yl[3] << ", " << A_pT2yl[4] << ", " << A_pT2yl[5] << ", " << A_pT2yl[6] << ", " << A_pT2yl[7] << ", " << A_pT2yl[8] << "}" << endl;
    Output << "Minv-2,<Err> ={" << deltaA_pT2yl[0] << ", " << deltaA_pT2yl[1] << ", " << deltaA_pT2yl[2] << ", " << deltaA_pT2yl[3] << ", " << deltaA_pT2yl[4] << ", " << deltaA_pT2yl[5] << ", " << deltaA_pT2yl[6] << ", " << deltaA_pT2yl[7] << ", " << deltaA_pT2yl[8] << "}" << endl;
    Output << "Minv-3,<A_{UT}> ={" << A_pT3yl[0] << ", " << A_pT3yl[1] << ", " << A_pT3yl[2] << ", " << A_pT3yl[3] << ", " << A_pT3yl[4] << ", " << A_pT3yl[5] << ", " << A_pT3yl[6] << ", " << A_pT3yl[7] << ", " << A_pT3yl[8] << "}" << endl;
    Output << "Minv-3,<Err> ={" << deltaA_pT3yl[0] << ", " << deltaA_pT3yl[1] << ", " << deltaA_pT3yl[2] << ", " << deltaA_pT3yl[3] << ", " << deltaA_pT3yl[4] << ", " << deltaA_pT3yl[5] << ", " << deltaA_pT3yl[6] << ", " << deltaA_pT3yl[7] << ", " << deltaA_pT3yl[8] << "}" << endl;
    Output << "Minv-4,<A_{UT}> ={" << A_pT4yl[0] << ", " << A_pT4yl[1] << ", " << A_pT4yl[2] << ", " << A_pT4yl[3] << ", " << A_pT4yl[4] << ", " << A_pT4yl[5] << ", " << A_pT4yl[6] << ", " << A_pT4yl[7] << ", " << A_pT4yl[8] << "}" << endl;
    Output << "Minv-4,<Err> ={" << deltaA_pT4yl[0] << ", " << deltaA_pT4yl[1] << ", " << deltaA_pT4yl[2] << ", " << deltaA_pT4yl[3] << ", " << deltaA_pT4yl[4] << ", " << deltaA_pT4yl[5] << ", " << deltaA_pT4yl[6] << ", " << deltaA_pT4yl[7] << ", " << deltaA_pT4yl[8] << "}" << endl;
    Output << "Minv-5,<A_{UT}> ={" << A_pT5yl[0] << ", " << A_pT5yl[1] << ", " << A_pT5yl[2] << ", " << A_pT5yl[3] << ", " << A_pT5yl[4] << ", " << A_pT5yl[5] << ", " << A_pT5yl[6] << ", " << A_pT5yl[7] << ", " << A_pT5yl[8] << "}" << endl;
    Output << "Minv-5,<Err> ={" << deltaA_pT5yl[0] << ", " << deltaA_pT5yl[1] << ", " << deltaA_pT5yl[2] << ", " << deltaA_pT5yl[3] << ", " << deltaA_pT5yl[4] << ", " << deltaA_pT5yl[5] << ", " << deltaA_pT5yl[6] << ", " << deltaA_pT5yl[7] << ", " << deltaA_pT5yl[8] << "}" << endl;

    Output << "                                                                                                  " << endl;
    Output << "                                                                                                  " << endl;
    Output << "<<<<<<<<<<<       Eta > 0, YELLOw BEAM A_{UT}   >>>>>>>>>>>>>>>>>" << endl;
    Output << "Minv-1,<A_{UT}> ={" << A_pT1yg[0] << ", " << A_pT1yg[1] << ", " << A_pT1yg[2] << ", " << A_pT1yg[3] << ", " << A_pT1yg[4] << ", " << A_pT1yg[5] << ", " << A_pT1yg[6] << ", " << A_pT1yg[7] << ", " << A_pT1yg[8] << "}" << endl;
    Output << "Minv-1,<Err> ={" << deltaA_pT1yg[0] << ", " << deltaA_pT1yg[1] << ", " << deltaA_pT1yg[2] << ", " << deltaA_pT1yg[3] << ", " << deltaA_pT1yg[4] << ", " << deltaA_pT1yg[5] << ", " << deltaA_pT1yg[6] << ", " << deltaA_pT1yg[7] << ", " << deltaA_pT1yg[8] << "}" << endl;
    Output << "Minv-2,<A_{UT}> ={" << A_pT2yg[0] << ", " << A_pT2yg[1] << ", " << A_pT2yg[2] << ", " << A_pT2yg[3] << ", " << A_pT2yg[4] << ", " << A_pT2yg[5] << ", " << A_pT2yg[6] << ", " << A_pT2yg[7] << ", " << A_pT2yg[8] << "}" << endl;
    Output << "Minv-2,<Err> ={" << deltaA_pT2yg[0] << ", " << deltaA_pT2yg[1] << ", " << deltaA_pT2yg[2] << ", " << deltaA_pT2yg[3] << ", " << deltaA_pT2yg[4] << ", " << deltaA_pT2yg[5] << ", " << deltaA_pT2yg[6] << ", " << deltaA_pT2yg[7] << ", " << deltaA_pT2yg[8] << "}" << endl;
    Output << "Minv-3,<A_{UT}> ={" << A_pT3yg[0] << ", " << A_pT3yg[1] << ", " << A_pT3yg[2] << ", " << A_pT3yg[3] << ", " << A_pT3yg[4] << ", " << A_pT3yg[5] << ", " << A_pT3yg[6] << ", " << A_pT3yg[7] << ", " << A_pT3yg[8] << "}" << endl;
    Output << "Minv-3,<Err> ={" << deltaA_pT3yg[0] << ", " << deltaA_pT3yg[1] << ", " << deltaA_pT3yg[2] << ", " << deltaA_pT3yg[3] << ", " << deltaA_pT3yg[4] << ", " << deltaA_pT3yg[5] << ", " << deltaA_pT3yg[6] << ", " << deltaA_pT3yg[7] << ", " << deltaA_pT3yg[8] << "}" << endl;
    Output << "Minv-4,<A_{UT}> ={" << A_pT4yg[0] << ", " << A_pT4yg[1] << ", " << A_pT4yg[2] << ", " << A_pT4yg[3] << ", " << A_pT4yg[4] << ", " << A_pT4yg[5] << ", " << A_pT4yg[6] << ", " << A_pT4yg[7] << ", " << A_pT4yg[8] << "}" << endl;
    Output << "Minv-4,<Err> ={" << deltaA_pT4yg[0] << ", " << deltaA_pT4yg[1] << ", " << deltaA_pT4yg[2] << ", " << deltaA_pT4yg[3] << ", " << deltaA_pT4yg[4] << ", " << deltaA_pT4yg[5] << ", " << deltaA_pT4yg[6] << ", " << deltaA_pT4yg[7] << ", " << deltaA_pT4yg[8] << "}" << endl;
    Output << "Minv-5,<A_{UT}> ={" << A_pT5yg[0] << ", " << A_pT5yg[1] << ", " << A_pT5yg[2] << ", " << A_pT5yg[3] << ", " << A_pT5yg[4] << ", " << A_pT5yg[5] << ", " << A_pT5yg[6] << ", " << A_pT5yg[7] << ", " << A_pT5yg[8] << "}" << endl;
    Output << "Minv-5,<Err> ={" << deltaA_pT5yg[0] << ", " << deltaA_pT5yg[1] << ", " << deltaA_pT5yg[2] << ", " << deltaA_pT5yg[3] << ", " << deltaA_pT5yg[4] << ", " << deltaA_pT5yg[5] << ", " << deltaA_pT5yg[6] << ", " << deltaA_pT5yg[7] << ", " << deltaA_pT5yg[8] << "}" << endl;

    Output << "                                                                                                  " << endl;
    Output << "                                                                                                  " << endl;
    Output << "                                                                                                  " << endl;
    Output << "                                      AVERAGE A_UT                                                " << endl;
    Output << "********************   FINAL ASYMMETRY , CONE > 0.7, pT Binning for Average BLUE and YELLOW BEAM****************" << endl;
    Output << "<<<<<<<<<<<<<<<   pT Bin Averages for asymemtry plot  >>>>>>>>>>>>>>" << endl;

    Output << "                                                                                                  " << endl;
    Output << "                                                                                                  " << endl;

    Output << "Minv-1 pTAvgGt[0][9]={" << pTAvgGt[0][0] << "," << pTAvgGt[0][1] << "," << pTAvgGt[0][2] << "," << pTAvgGt[0][3] << "," << pTAvgGt[0][4] << "," << pTAvgGt[0][5] << "," << pTAvgGt[0][6] << "," << pTAvgGt[0][7] << "," << pTAvgGt[0][8] << "}" << endl;
    Output << "Minv-2 pTAvgGt[1][9]={" << pTAvgGt[1][0] << "," << pTAvgGt[1][1] << "," << pTAvgGt[1][2] << "," << pTAvgGt[1][3] << "," << pTAvgGt[1][4] << "," << pTAvgGt[1][5] << "," << pTAvgGt[1][6] << "," << pTAvgGt[1][7] << "," << pTAvgGt[1][8] << "}" << endl;
    Output << "Minv-3  pTAvgGt[2][9]={" << pTAvgGt[2][0] << "," << pTAvgGt[2][1] << "," << pTAvgGt[2][2] << "," << pTAvgGt[2][3] << "," << pTAvgGt[2][4] << "," << pTAvgGt[2][5] << "," << pTAvgGt[2][6] << "," << pTAvgGt[2][7] << "," << pTAvgGt[2][8] << "}" << endl;
    Output << "Minv-4  pTAvgGt[3][9]={" << pTAvgGt[3][0] << "," << pTAvgGt[3][1] << "," << pTAvgGt[3][2] << "," << pTAvgGt[3][3] << "," << pTAvgGt[3][4] << "," << pTAvgGt[3][5] << "," << pTAvgGt[3][6] << "," << pTAvgGt[3][7] << "," << pTAvgGt[3][8] << "}" << endl;
    Output << "Minv-5 pTAvgGt[4][9]={" << pTAvgGt[4][0] << "," << pTAvgGt[4][1] << "," << pTAvgGt[4][2] << "," << pTAvgGt[4][3] << "," << pTAvgGt[4][4] << "," << pTAvgGt[4][5] << "," << pTAvgGt[4][6] << "," << pTAvgGt[4][7] << "," << pTAvgGt[4][8] << "}" << endl;

    Output << "                                                                                                  " << endl;
    Output << "                                                                                                  " << endl;

    Output << "Minv-1 pTAvgLt[0][9]={" << pTAvgLt[0][0] << "," << pTAvgLt[0][1] << "," << pTAvgLt[0][2] << "," << pTAvgLt[0][3] << "," << pTAvgLt[0][4] << "," << pTAvgLt[0][5] << "," << pTAvgLt[0][6] << "," << pTAvgLt[0][7] << "," << pTAvgLt[0][8] << "}" << endl;
    Output << "Minv-2 pTAvgLt[1][9]={" << pTAvgLt[1][0] << "," << pTAvgLt[1][1] << "," << pTAvgLt[1][2] << "," << pTAvgLt[1][3] << "," << pTAvgLt[1][4] << "," << pTAvgLt[1][5] << "," << pTAvgLt[1][6] << "," << pTAvgLt[1][7] << "," << pTAvgLt[1][8] << "}" << endl;
    Output << "Minv-3 pTAvgLt[2][9]={" << pTAvgLt[2][0] << "," << pTAvgLt[2][1] << "," << pTAvgLt[2][2] << "," << pTAvgLt[2][3] << "," << pTAvgLt[2][4] << "," << pTAvgLt[2][5] << "," << pTAvgLt[2][6] << "," << pTAvgLt[2][7] << "," << pTAvgLt[2][8] << "}" << endl;
    Output << "Minv-4 pTAvgLt[3][9]={" << pTAvgLt[3][0] << "," << pTAvgLt[3][1] << "," << pTAvgLt[3][2] << "," << pTAvgLt[3][3] << "," << pTAvgLt[3][4] << "," << pTAvgLt[3][5] << "," << pTAvgLt[3][6] << "," << pTAvgLt[3][7] << "," << pTAvgLt[3][8] << "}" << endl;
    Output << "Minv-5 pTAvgLt[4][9]={" << pTAvgLt[4][0] << "," << pTAvgLt[4][1] << "," << pTAvgLt[4][2] << "," << pTAvgLt[4][3] << "," << pTAvgLt[4][4] << "," << pTAvgLt[4][5] << "," << pTAvgLt[4][6] << "," << pTAvgLt[4][7] << "," << pTAvgLt[4][8] << "}" << endl;

    Output << "                                                                                                  " << endl;
    Output << "                                                                                                  " << endl;

    Output << "<<<<<<<<<<<       Eta < 0, Average A_{UT}   >>>>>>>>>>>>>>>>>" << endl;
    Output << "Minv-1,<A_{UT}> ={" << avgA_pT1l[0] << ", " << avgA_pT1l[1] << ", " << avgA_pT1l[2] << ", " << avgA_pT1l[3] << ", " << avgA_pT1l[4] << ", " << avgA_pT1l[5] << ", " << avgA_pT1l[6] << ", " << avgA_pT1l[7] << ", " << avgA_pT1l[8] << "}" << endl;
    Output << "Minv-1,<Err>={" << errA_pT1l[0] << ", " << errA_pT1l[1] << ", " << errA_pT1l[2] << ", " << errA_pT1l[3] << ", " << errA_pT1l[4] << ", " << errA_pT1l[5] << ", " << errA_pT1l[6] << ", " << errA_pT1l[7] << ", " << errA_pT1l[8] << "}" << endl;
    Output << "Minv-2,<A_{UT}>={" << avgA_pT2l[0] << ", " << avgA_pT2l[1] << ", " << avgA_pT2l[2] << ", " << avgA_pT2l[3] << ", " << avgA_pT2l[4] << ", " << avgA_pT2l[5] << ", " << avgA_pT2l[6] << ", " << avgA_pT2l[7] << ", " << avgA_pT2l[8] << "}" << endl;
    Output << "Minv-2,<Err> ={" << errA_pT2l[0] << ", " << errA_pT2l[1] << ", " << errA_pT2l[2] << ", " << errA_pT2l[3] << ", " << errA_pT2l[4] << ", " << errA_pT2l[5] << ", " << errA_pT2l[6] << ", " << errA_pT2l[7] << ", " << errA_pT2l[8] << "}" << endl;
    Output << "Minv-3,<A_{UT}> ={" << avgA_pT3l[0] << ", " << avgA_pT3l[1] << ", " << avgA_pT3l[2] << ", " << avgA_pT3l[3] << ", " << avgA_pT3l[4] << ", " << avgA_pT3l[5] << ", " << avgA_pT3l[6] << ", " << avgA_pT3l[7] << ", " << avgA_pT3l[8] << "}" << endl;
    Output << "Minv-3,<Err> ={" << errA_pT3l[0] << ", " << errA_pT3l[1] << ", " << errA_pT3l[2] << ", " << errA_pT3l[3] << ", " << errA_pT3l[4] << ", " << errA_pT3l[5] << ", " << errA_pT3l[6] << ", " << errA_pT3l[7] << ", " << errA_pT3l[8] << "}" << endl;
    Output << "Minv-4,<A_{UT}> ={" << avgA_pT4l[0] << ", " << avgA_pT4l[1] << ", " << avgA_pT4l[2] << ", " << avgA_pT4l[3] << ", " << avgA_pT4l[4] << ", " << avgA_pT4l[5] << ", " << avgA_pT4l[6] << ", " << avgA_pT4l[7] << ", " << avgA_pT4l[8] << "}" << endl;
    Output << "Minv-4,<Err> ={" << errA_pT4l[0] << ", " << errA_pT4l[1] << ", " << errA_pT4l[2] << ", " << errA_pT4l[3] << ", " << errA_pT4l[4] << ", " << errA_pT4l[5] << ", " << errA_pT4l[6] << ", " << errA_pT4l[7] << ", " << errA_pT4l[8] << "}" << endl;
    Output << "Minv-5,<A_{UT}> ={" << avgA_pT5l[0] << ", " << avgA_pT5l[1] << ", " << avgA_pT5l[2] << ", " << avgA_pT5l[3] << ", " << avgA_pT5l[4] << ", " << avgA_pT5l[5] << ", " << avgA_pT5l[6] << ", " << avgA_pT5l[7] << ", " << avgA_pT5l[8] << "}" << endl;
    Output << "Minv-5,<Err> ={" << errA_pT5l[0] << ", " << errA_pT5l[1] << ", " << errA_pT5l[2] << ", " << errA_pT5l[3] << ", " << errA_pT5l[4] << ", " << errA_pT5l[5] << ", " << errA_pT5l[6] << ", " << errA_pT5l[7] << ", " << errA_pT5l[8] << "}" << endl;

    Output << "                                                                                                  " << endl;
    Output << "                                                                                                  " << endl;
    Output << "<<<<<<<<<<<       Eta > 0, Average A_{UT}   >>>>>>>>>>>>>>>>>" << endl;
    Output << "Minv-1,<A_{UT}> ={" << avgA_pT1g[0] << ", " << avgA_pT1g[1] << ", " << avgA_pT1g[2] << ", " << avgA_pT1g[3] << ", " << avgA_pT1g[4] << ", " << avgA_pT1g[5] << ", " << avgA_pT1g[6] << ", " << avgA_pT1g[7] << ", " << avgA_pT1g[8] << "}" << endl;
    Output << "Minv-1,<Err> ={" << errA_pT1g[0] << ", " << errA_pT1g[1] << ", " << errA_pT1g[2] << ", " << errA_pT1g[3] << ", " << errA_pT1g[4] << ", " << errA_pT1g[5] << ", " << errA_pT1g[6] << ", " << errA_pT1g[7] << ", " << errA_pT1g[8] << "}" << endl;
    Output << "Minv-2,<A_{UT}> ={" << avgA_pT2g[0] << ", " << avgA_pT2g[1] << ", " << avgA_pT2g[2] << ", " << avgA_pT2g[3] << ", " << avgA_pT2g[4] << ", " << avgA_pT2g[5] << ", " << avgA_pT2g[6] << ", " << avgA_pT2g[7] << ", " << avgA_pT2g[8] << "}" << endl;
    Output << "Minv-2,<Err> ={" << errA_pT2g[0] << ", " << errA_pT2g[1] << ", " << errA_pT2g[2] << ", " << errA_pT2g[3] << ", " << errA_pT2g[4] << ", " << errA_pT2g[5] << ", " << errA_pT2g[6] << ", " << errA_pT2g[7] << ", " << errA_pT2g[8] << "}" << endl;
    Output << "Minv-3,<A_{UT}> ={" << avgA_pT3g[0] << ", " << avgA_pT3g[1] << ", " << avgA_pT3g[2] << ", " << avgA_pT3g[3] << ", " << avgA_pT3g[4] << ", " << avgA_pT3g[5] << ", " << avgA_pT3g[6] << ", " << avgA_pT3g[7] << ", " << avgA_pT3g[8] << "}" << endl;
    Output << "Minv-3,<Err> ={" << errA_pT3g[0] << ", " << errA_pT3g[1] << ", " << errA_pT3g[2] << ", " << errA_pT3g[3] << ", " << errA_pT3g[4] << ", " << errA_pT3g[5] << ", " << errA_pT3g[6] << ", " << errA_pT3g[7] << ", " << errA_pT3g[8] << "}" << endl;
    Output << "Minv-4,<A_{UT}> ={" << avgA_pT4g[0] << ", " << avgA_pT4g[1] << ", " << avgA_pT4g[2] << ", " << avgA_pT4g[3] << ", " << avgA_pT4g[4] << ", " << avgA_pT4g[5] << ", " << avgA_pT4g[6] << ", " << avgA_pT4g[7] << ", " << avgA_pT4g[8] << "}" << endl;
    Output << "Minv-4,<Err> ={" << errA_pT4g[0] << ", " << errA_pT4g[1] << ", " << errA_pT4g[2] << ", " << errA_pT4g[3] << ", " << errA_pT4g[4] << ", " << errA_pT4g[5] << ", " << errA_pT4g[6] << ", " << errA_pT4g[7] << ", " << errA_pT4g[8] << "}" << endl;
    Output << "Minv-5,<A_{UT}> ={" << avgA_pT5g[0] << ", " << avgA_pT5g[1] << ", " << avgA_pT5g[2] << ", " << avgA_pT5g[3] << ", " << avgA_pT5g[4] << ", " << avgA_pT5g[5] << ", " << avgA_pT5g[6] << ", " << avgA_pT5g[7] << ", " << avgA_pT5g[8] << "}" << endl;
    Output << "Minv-5,<Err> ={" << errA_pT5g[0] << ", " << errA_pT5g[1] << ", " << errA_pT5g[2] << ", " << errA_pT5g[3] << ", " << errA_pT5g[4] << ", " << errA_pT5g[5] << ", " << errA_pT5g[6] << ", " << errA_pT5g[7] << ", " << errA_pT5g[8] << "}" << endl;

    // Average Asymmetry
    double AvgA_G[9] = {0};
    double dA_G[9] = {0};
    double AvgA_L[9] = {0};
    double dA_L[9] = {0};

    double WAvgA_G[9] = {0};
    double WdA_G[9] = {0};
    double WAvgA_L[9] = {0};
    double WdA_L[9] = {0};

    double A_BG[9] = {0};
    double dA_BG[9] = {0};
    double A_BL[9] = {0};
    double dA_BL[9] = {0};
    double A_YG[9] = {0};
    double dA_YG[9] = {0};
    double A_YL[9] = {0};
    double dA_YL[9] = {0};

    double AvgpT_G[9] = {0};
    double AvgpT_L[9] = {0};
    double pT_BG[9] = {0};
    double pT_BL[9] = {0};
    double pT_YG[9] = {0};
    double pT_YL[9] = {0};

    vector<double> A_UT[3][2];
    vector<double> dA_UT[3][2];
    vector<double> PT[3][2];
    vector<double> PIDSYS[3][2];

    vector<double> WA_UT[3][2];
    vector<double> WdA_UT[3][2];
    vector<double> WPIDSYS[3][2];

    vector<double> WPIDSYS_TPC[3][2];
    vector<double> WPIDSYS_TPCorTOF[3][2];
    vector<double> WPIDSYS_TOF[3][2];
    vector<double> WTOTSYS[3][2];
    // double A_UT[3][2][9]={0};
    // double dA_UT[3][2][9]={0};

    // TGraphErrors *gr[3][2][5];
    TGraphErrors *GRAPH[3][2][5];
    TGraphErrors *gr_sys[3][2][5];

    TGraphErrors *WGRAPH[3][2][5];
    TGraphErrors *Wgr_sys[3][2][5];

    TGraphErrors *Wgr_sys_TPC[3][2][5];
    TGraphErrors *Wgr_sys_TPCorTOF[3][2][5];
    TGraphErrors *Wgr_sys_TOF[3][2][5];
    TGraphErrors *Wgr_TOTSYS[3][2][5];

    TGraphErrors *gr_BG[5];
    TGraphErrors *gr_YG[5];
    TGraphErrors *gr_BL[5];
    TGraphErrors *gr_YL[5];
    TGraphErrors *gr_AvgG[5];
    TGraphErrors *gr_AvgL[5];
    TGraphErrors *GRAPH_R11[5];

    // To Compare with Run11 Result;;
    double pT_R11_G_P1[9] = {4.212, 5.085, 6.248, 7.994, 10.13, 12.16, 16.87};
    double A_R11_G_P1[9] = {0.003871, 0.004516, 0.006774, -0.003871, 0.01129, -0.01516, 0.01226};
    double dA_R11_G_P1[9] = {0.006452, 0.005806, 0.006129, 0.006774, 0.0129, 0.01484, 0.02161};
    double pT_R11_G_P2[9] = {4.15, 5.06, 6.257, 7.982, 10.09, 12.25, 17.08};
    double A_R11_G_P2[9] = {0.005484, 0.0003226, 0.01613, -0.005806, 0.01645, 0.0371, 0.02968};
    double dA_R11_G_P2[9] = {0.007097, 0.006452, 0.006129, 0.006452, 0.01161, 0.01194, 0.01677};
    double pT_R11_G_P3[9] = {4.174, 5.104, 6.278, 8.04, 10.14, 12.25, 17.48};
    double A_R11_G_P3[9] = {0.005161, 0.01287, 0.0003226, 0.01065, 0.03097, 0.003871, 0.02871};
    double dA_R11_G_P3[9] = {0.007097, 0.006129, 0.006452, 0.006452, 0.01129, 0.01194, 0.01516};
    double pT_R11_G_P4[9] = {4.164, 5.085, 6.248, 7.994, 10.08, 12.21, 17.69};
    double A_R11_G_P4[9] = {0.003981, -0.003148, 0.01338, 0.005278, 0.01662, 0.02472, 0.03639};
    double dA_R11_G_P4[9] = {0.005833, 0.005509, 0.005509, 0.005509, 0.009722, 0.009722, 0.01199};
    double pT_R11_G_P5[9] = {4.16, 5.127, 6.287, 8.076, 10.15, 12.28, 18.08};
    double A_R11_G_P5[9] = {0.008843, 0.004306, -0.001852, 0.01241, 0.01565, 0.01759, 0.03347};
    double dA_R11_G_P5[9] = {0.002106, 0.003565, 0.003241, 0.003241, 0.005509, 0.005185, 0.005833};

    double pT_R11_G[9] = {0};
    double A_R11_G[9] = {0};
    double dA_R11_G[9] = {0};

    for (int pad = 0; pad < 5; pad++)
    {
        if (pad == 0)
        {
            for (int mbin = 0; mbin < 9; mbin++)
            {
                pT_R11_G[mbin] = pT_R11_G_P1[mbin];
                A_R11_G[mbin] = A_R11_G_P1[mbin];
                dA_R11_G[mbin] = dA_R11_G_P1[mbin];
            }
        }
        if (pad == 1)
        {
            for (int mbin = 0; mbin < 9; mbin++)
            {
                // 2nd pad <Minv> of Run17 and comparable to 3rd pad of Run11.
                pT_R11_G[mbin] = pT_R11_G_P3[mbin];
                A_R11_G[mbin] = A_R11_G_P3[mbin];
                dA_R11_G[mbin] = dA_R11_G_P3[mbin];
            }
        }
        if (pad == 2)
        {
            for (int mbin = 0; mbin < 9; mbin++)
            {
                // 3rd pad <Minv> of Run17 and comparable to 4th pad of Run11.
                pT_R11_G[mbin] = pT_R11_G_P4[mbin];
                A_R11_G[mbin] = A_R11_G_P4[mbin];
                dA_R11_G[mbin] = dA_R11_G_P4[mbin];
            }
        }
        if (pad == 3)
        {
            for (int mbin = 0; mbin < 9; mbin++)
            {
                // 4th pad <Minv> of Run17 and comparable to 5th pad of Run11.
                pT_R11_G[mbin] = pT_R11_G_P5[mbin];
                A_R11_G[mbin] = A_R11_G_P5[mbin];
                dA_R11_G[mbin] = dA_R11_G_P5[mbin];
            }
        }
        if (pad == 4)
        {
            for (int mbin = 0; mbin < 9; mbin++)
            {
                // pT_R11_G[mbin] = pT_R11_G_P5[mbin];
                // A_R11_G[mbin] = A_R11_G_P5[mbin];
                // dA_R11_G[mbin] = dA_R11_G_P5[mbin];
            }
        }
        GRAPH_R11[pad] = new TGraphErrors(9, pT_R11_G, A_R11_G, 0, dA_R11_G);
    }

    //=============================Compare with JAMDiFF theory ============================

    TGraphErrors *GRAPH_JAMDiFF_G[5];
    TGraph *GRAPH_JAMDiFF_CentralVal_G[5];

    TGraphErrors *GRAPH_JAMDiFF_L[5];
    TGraph *GRAPH_JAMDiFF_CentralVal_L[5];

    double pT_JAMDiFF_G[9];
    double A_JAMDiFF_G[9];
    double E_JAMDiFF_G[9];

    double pT_JAMDiFF_L[9];
    double A_JAMDiFF_L[9];
    double E_JAMDiFF_L[9];

    double JAM_pT510_P1_G[9] = {3.6, 4.2, 4.7, 5.2, 5.7, 6.3, 7.2, 8.5, 12.2};
    double JAM_pT510_P2_G[9] = {3.7, 4.4, 4.9, 5.4, 6, 6.8, 7.8, 9.3, 13.6};
    double JAM_pT510_P3_G[9] = {3.8, 4.5, 5, 5.6, 6.3, 7.1, 8.2, 10, 14.5};
    double JAM_pT510_P4_G[9] = {3.8, 4.4, 5, 5.6, 6.3, 7.2, 8.4, 10.2, 14.9};
    double JAM_pT510_P5_G[9] = {4.2, 5, 5.6, 6.3, 7.1, 8.1, 9.4, 11.5, 16.4};

    double JAM_A510_P1_G[9] = {0.001312132, 0.001238796, 0.001316385, 0.001451633, 0.001614641, 0.001823855, 0.002168886, 0.00276044, 0.00462898};
    double JAM_A510_P2_G[9] = {0.004050971, 0.00395322, 0.004224003, 0.004685851, 0.00527832, 0.006234868, 0.007523822, 0.009603016, 0.016019271};
    double JAM_A510_P3_G[9] = {0.005413927, 0.005325898, 0.005741395, 0.00647483, 0.007510655, 0.008787163, 0.010674915, 0.014245735, 0.023277272};
    double JAM_A510_P4_G[9] = {0.008246897, 0.008095234, 0.008923283, 0.010043087, 0.011646116, 0.014163971, 0.017672042, 0.0240481, 0.042727704};
    double JAM_A510_P5_G[9] = {0.006514535, 0.007760767, 0.008928565, 0.008817759, 0.010285121, 0.009460434, 0.011717572, 0.011031579, 0.019346454};

    double JAM_E510_P1_G[9] = {0.000791104, 0.000717357, 0.000737903, 0.000787978, 0.000847544, 0.000918129, 0.001023543, 0.001185906, 0.001605443};
    double JAM_E510_P2_G[9] = {0.001831177, 0.001698468, 0.001746152, 0.001863752, 0.002003393, 0.002215442, 0.002456975, 0.002770208, 0.003439282};
    double JAM_E510_P3_G[9] = {0.002098447, 0.001951959, 0.002021477, 0.002179107, 0.002394783, 0.002623155, 0.002895131, 0.003283833, 0.003857956};
    double JAM_E510_P4_G[9] = {0.003172949, 0.00292236, 0.002988946, 0.003194102, 0.003544317, 0.004083616, 0.004683846, 0.005471877, 0.006732057};
    double JAM_E510_P5_G[9] = {0.002629226, 0.003442705, 0.004472552, 0.004884617, 0.005630013, 0.00513854, 0.0062335, 0.005678914, 0.008856724};
    // ======= JAM for eta<0 ======
    double JAM_pT510_P1_L[9] = {3.6, 4.2, 4.7, 5.2, 5.7, 6.3, 7.2, 8.5, 12.2};
    double JAM_pT510_P2_L[9] = {3.7, 4.4, 4.9, 5.4, 6, 6.8, 7.8, 9.3, 13.6};
    double JAM_pT510_P3_L[9] = {3.8, 4.5, 5, 5.6, 6.3, 7.1, 8.2, 10, 14.5};
    double JAM_pT510_P4_L[9] = {3.8, 4.4, 5, 5.6, 6.3, 7.2, 8.4, 10.2, 14.9};
    double JAM_pT510_P5_L[9] = {4.2, 5, 5.6, 6.3, 7.1, 8.1, 9.4, 11.5, 16.4};

    double JAM_A510_P1_L[9] = {0.000153878, 0.000217861, 0.000250227, 0.000273722, 0.00029959, 0.000321068, 0.000361832, 0.000416118, 0.000581353};
    double JAM_A510_P2_L[9] = {0.000559319, 0.000768354, 0.000855695, 0.000923982, 0.001008683, 0.001103997, 0.001226636, 0.001416006, 0.001969838};
    double JAM_A510_P3_L[9] = {0.000846702, 0.001153174, 0.001263544, 0.00137806, 0.001491551, 0.001624896, 0.001813505, 0.002087751, 0.002866723};
    double JAM_A510_P4_L[9] = {0.001557014, 0.002078998, 0.002331342, 0.002563129, 0.002795331, 0.003048915, 0.003457702, 0.004023966, 0.00567386};
    double JAM_A510_P5_L[9] = {0.001736956, 0.002292128, 0.002527055, 0.002370012, 0.002622776, 0.002216065, 0.002547217, 0.002115148, 0.002908497};

    double JAM_E510_P1_L[9] = {0.000132126, 0.000158689, 0.000161433, 0.000159638, 0.00016041, 0.000159728, 0.000167354, 0.000182709, 0.000239407};
    double JAM_E510_P2_L[9] = {0.000427339, 0.00046253, 0.0004529, 0.000438344, 0.000426474, 0.000414084, 0.000413681, 0.000435835, 0.000536128};
    double JAM_E510_P3_L[9] = {0.000611319, 0.000650335, 0.000622759, 0.000586807, 0.000548351, 0.000520117, 0.000508975, 0.000530919, 0.000684843};
    double JAM_E510_P4_L[9] = {0.001259802, 0.001394831, 0.001355718, 0.001296315, 0.001212347, 0.001104517, 0.001013501, 0.000973637, 0.001259322};
    double JAM_E510_P5_L[9] = {0.001265564, 0.001505364, 0.001659398, 0.001528869, 0.001543, 0.001193386, 0.001269962, 0.001006418, 0.001392541};

    for (int pad = 0; pad < 5; pad++)
    {
        if (pad == 0)
        {
            for (int mbin = 0; mbin < sizeof(JAM_pT510_P1_G) / sizeof(JAM_pT510_P1_G[0]); mbin++)
            {

                pT_JAMDiFF_G[mbin] = JAM_pT510_P1_G[mbin];
                A_JAMDiFF_G[mbin] = JAM_A510_P1_G[mbin];
                E_JAMDiFF_G[mbin] = JAM_E510_P1_G[mbin];

                pT_JAMDiFF_L[mbin] = JAM_pT510_P1_L[mbin];
                A_JAMDiFF_L[mbin] = JAM_A510_P1_L[mbin];
                E_JAMDiFF_L[mbin] = JAM_E510_P1_L[mbin];
            }
        }
        if (pad == 1)
        {
            for (int mbin = 0; mbin < sizeof(JAM_pT510_P2_G) / sizeof(JAM_pT510_P2_G[0]); mbin++)
            {

                pT_JAMDiFF_G[mbin] = JAM_pT510_P2_G[mbin];
                A_JAMDiFF_G[mbin] = JAM_A510_P2_G[mbin];
                E_JAMDiFF_G[mbin] = JAM_E510_P2_G[mbin];

                pT_JAMDiFF_L[mbin] = JAM_pT510_P2_L[mbin];
                A_JAMDiFF_L[mbin] = JAM_A510_P2_L[mbin];
                E_JAMDiFF_L[mbin] = JAM_E510_P2_L[mbin];
            }
        }
        if (pad == 2)
        {
            for (int mbin = 0; mbin < sizeof(JAM_pT510_P3_G) / sizeof(JAM_pT510_P3_G[0]); mbin++)
            {

                pT_JAMDiFF_G[mbin] = JAM_pT510_P3_G[mbin];
                A_JAMDiFF_G[mbin] = JAM_A510_P3_G[mbin];
                E_JAMDiFF_G[mbin] = JAM_E510_P3_G[mbin];

                pT_JAMDiFF_L[mbin] = JAM_pT510_P3_L[mbin];
                A_JAMDiFF_L[mbin] = JAM_A510_P3_L[mbin];
                E_JAMDiFF_L[mbin] = JAM_E510_P3_L[mbin];
            }
        }
        if (pad == 3)
        {
            for (int mbin = 0; mbin < sizeof(JAM_pT510_P4_G) / sizeof(JAM_pT510_P4_G[0]); mbin++)
            {

                pT_JAMDiFF_G[mbin] = JAM_pT510_P4_G[mbin];
                A_JAMDiFF_G[mbin] = JAM_A510_P4_G[mbin];
                E_JAMDiFF_G[mbin] = JAM_E510_P4_G[mbin];

                pT_JAMDiFF_L[mbin] = JAM_pT510_P4_L[mbin];
                A_JAMDiFF_L[mbin] = JAM_A510_P4_L[mbin];
                E_JAMDiFF_L[mbin] = JAM_E510_P4_L[mbin];
            }
        }
        if (pad == 4)
        {
            for (int mbin = 0; mbin < sizeof(JAM_pT510_P5_G) / sizeof(JAM_pT510_P5_G[0]); mbin++)
            {

                pT_JAMDiFF_G[mbin] = JAM_pT510_P5_G[mbin];
                A_JAMDiFF_G[mbin] = JAM_A510_P5_G[mbin];
                E_JAMDiFF_G[mbin] = JAM_E510_P5_G[mbin];

                pT_JAMDiFF_L[mbin] = JAM_pT510_P5_L[mbin];
                A_JAMDiFF_L[mbin] = JAM_A510_P5_L[mbin];
                E_JAMDiFF_L[mbin] = JAM_E510_P5_L[mbin];
            }
        }

        GRAPH_JAMDiFF_G[pad] = new TGraphErrors(sizeof(pT_JAMDiFF_G) / sizeof(pT_JAMDiFF_G[0]), pT_JAMDiFF_G, A_JAMDiFF_G, 0, E_JAMDiFF_G);
        GRAPH_JAMDiFF_CentralVal_G[pad] = new TGraph(sizeof(pT_JAMDiFF_G) / sizeof(pT_JAMDiFF_G[0]), pT_JAMDiFF_G, A_JAMDiFF_G);

        GRAPH_JAMDiFF_L[pad] = new TGraphErrors(sizeof(pT_JAMDiFF_L) / sizeof(pT_JAMDiFF_L[0]), pT_JAMDiFF_L, A_JAMDiFF_L, 0, E_JAMDiFF_L);
        GRAPH_JAMDiFF_CentralVal_L[pad] = new TGraph(sizeof(pT_JAMDiFF_L) / sizeof(pT_JAMDiFF_L[0]), pT_JAMDiFF_L, A_JAMDiFF_L);
    }
    //============================= Compare with JAMDiFF theory End ============================

    // Proton and Pion separate
    // Pion Pair purity for sysyematic study;
    // double pionPairPtGt[5] = {0.827833, 0.832974, 0.810826, 0.781036, 0.711847};
    // double pionPairPtLt[5] = {0.86172, 0.835669, 0.809801, 0.775498, 0.688176};

    // Kaon and Proton combine
    double pionPairMGt[5] = {0.838334, 0.842922, 0.845917, 0.813829, 0.780844};
    double pionPairMLt[5] = {0.839661, 0.845568, 0.848364, 0.816115, 0.783557};
    double pionPairPtGt[5] = {0.833136, 0.854063, 0.86806, 0.869287, 0.83688};
    double pionPairPtLt[5] = {0.833513, 0.85424, 0.867603, 0.8696, 0.840349};

    // double pionPairPtGt_TPC[5] = {0.890242, 0.891946, 0.899971, 0.900113, 0.878957};
    // double pionPairPtLt_TPC[5] = {0.882257, 0.88228, 0.893677, 0.89351, 0.87064};
    // Final Result

    // the values in pionPairPtGt_TPC array are for A_UT vs Minv plots 5 different pT Bin. eg pionPairPtGt_TPC[0] is for 2.6<pT_pair<4.628 and pionPairGt_TPC[1] is for 4.628<pT_pair<5.643....and so on

    double pionPairPtGt_TPC[5] = {0.819141, 0.834208, 0.851712, 0.85491, 0.824758};
    double pionPairPtLt_TPC[5] = {0.826245, 0.83548, 0.852501, 0.85595, 0.827991};
    double pionPairMGt_TPC[5] = {0.829238, 0.83258, 0.834903, 0.800073, 0.765804};
    double pionPairMLt_TPC[5] = {0.83052, 0.834775, 0.83705, 0.802637, 0.768361};
    double pionPairEta_TPC[9] = {0.862276, 0.826708, 0.797671, 0.77381, 0.771638, 0.790094, 0.808581, 0.831038, 0.856093};

    // double pionPairPtGt_TPCorTOF[5] = {0.9416, 0.924323, 0.922676, 0.914645, 0.888728};
    // double pionPairPtLt_TPCorTOF[5] = {0.934431, 0.917851, 0.917512, 0.910225, 0.881531};
    // Final Result

    // the values in pionPairPtGt_TPCorTOF array are for A_UT vs Minv plots 5 different pT Bin. eg pionPairPtGt_TPCorTOF[0] is for 2.6<pT_pair<4.628 and pionPairPtGt_TPCorTOF[1] is for 4.628<pT_pair<5.643....and so on
    // the values in pionPairMGt_TPCorTOF array are for A_UT vs pT plots 5 different Minv Bin. eg pionPairMGt_TPCorTOF[0] is for 0.2<pT_pair<0.48 and pionPairMGt_TPCorTOF[1] is for 0.48<pT_pair<0.65....and so on

    //================================== Assigining PID sys from Relative difference of A_UT ===========================//

    double pionPairPtGt_TPCorToF[5] = {0};
    double pionPairPtLt_TPCorToF[5] = {0};
    double pionPairMGt_TPCorToF[5] = {0};
    double pionPairMLt_TPCorToF[5] = {0};
    double pionPairEta_TPCorToF[9] = {0};
    // double pionPairPtGt_TPCorToF[5] = {0.885761, 0.879613, 0.884539, 0.879013, 0.843039};
    // double pionPairPtLt_TPCorToF[5] = {0.890978, 0.880087, 0.88465, 0.879312, 0.845038};
    // double pionPairMGt_TPCorToF[5] = {0.860666, 0.862725, 0.865509, 0.839492, 0.807601};
    // double pionPairMLt_TPCorToF[5] = {0.86062, 0.863598, 0.866285, 0.840665, 0.808823};
    // double pionPairEta_TPCorToF[9] = {0.888436, 0.861023, 0.834911, 0.809617, 0.807478, 0.8277, 0.846784, 0.866677, 0.884602};
    // For Relative diff from 2D Bining;;
    // double TPCorTOF_rel_diff_PtGt_p0[5] = {-0.009152, -0.141025, -0.050961, -0.033589, 0.004620};
    // double TPCorTOF_rel_diff_PtGt_p0_err[5] = {0.301869, 0.182559, 0.125431, 0.078658, 0.044062};
    // double TOF_rel_diff_PtGt_p0[5] = {0.068537, -0.277958, -0.0199939, 0.065080, -0.077160};
    // double TOF_rel_diff_PtGt_p0_err[5] = {0.425255, 0.265452, 0.173306, 0.104850, 0.060484};
    // double TPCorTOF_rel_diff_PtLt_p0[5] = {0.354138, 0.253378, 0.212848, -0.232608, 0.024217};
    // double TPCorTOF_rel_diff_PtLt_p0_err[5] = {0.404620, 0.464401, 0.47733, 0.492759, 0.290594};
    // double TOF_rel_diff_PtLt_p0[5] = {0.282023, 0.310293, 0.389626, -0.296303, 0.443480};
    // double TOF_rel_diff_PtLt_p0_err[5] = {0.697597, 0.718014, 0.708123, 0.715590, 0.376371};

    double TPC_TPCandTOF_rel_diff_intPt_p0 = -0.078693;
    double TPC_TPCandTOF_rel_diff_intPt_p0_err = 0.055447;

    double TPC_TPCandTOF_rel_diff_intPt_p0_Lt = 0.143213;
    double TPC_TPCandTOF_rel_diff_intPt_p0_err_Lt = 0.350029;

    double TPC_TPCandTOF_rel_diff_intPt_p0_final = TMath::Max(fabs(TPC_TPCandTOF_rel_diff_intPt_p0), fabs(TPC_TPCandTOF_rel_diff_intPt_p0_err));
    double TPC_TPCandTOF_rel_diff_intPt_p0_final_Lt = TMath::Max(fabs(TPC_TPCandTOF_rel_diff_intPt_p0_Lt), fabs(TPC_TPCandTOF_rel_diff_intPt_p0_err_Lt));

    double TPC_TPCandTOF_rel_diff_intMinv_p0 = -0.057273;
    double TPC_TPCandTOF_rel_diff_intMinv_p0_err = 0.035243;

    double TPC_TPCandTOF_rel_diff_intMinv_p0_Lt = -0.089911;
    double TPC_TPCandTOF_rel_diff_intMinv_p0_err_Lt = 0.241005;

    double TPC_TPCandTOF_rel_diff_intPt_final = TMath::Max(fabs(TPC_TPCandTOF_rel_diff_intPt_p0), fabs(TPC_TPCandTOF_rel_diff_intPt_p0_err));
    double TPC_TPCandTOF_rel_diff_intPt_final_Lt = TMath::Max(fabs(TPC_TPCandTOF_rel_diff_intPt_p0_Lt), fabs(TPC_TPCandTOF_rel_diff_intPt_p0_err_Lt));

    double TPC_TPCandTOF_rel_diff_intMinv_final = TMath::Max(fabs(TPC_TPCandTOF_rel_diff_intMinv_p0), fabs(TPC_TPCandTOF_rel_diff_intMinv_p0_err));
    double TPC_TPCandTOF_rel_diff_intMinv_final_Lt = TMath::Max(fabs(TPC_TPCandTOF_rel_diff_intMinv_p0_Lt), fabs(TPC_TPCandTOF_rel_diff_intMinv_p0_err_Lt));

    double TPC_TPCandTOF_rel_Eta_diff_p0 = -0.109716;
    double TPC_TPCandTOF_rel_Eta_diff_p0_err = 0.062488;

    double TPC_TPCandTOF_rel_Eta_diff_final = TMath::Max(fabs(TPC_TPCandTOF_rel_Eta_diff_p0), fabs(TPC_TPCandTOF_rel_Eta_diff_p0_err));

    for (int jj = 0; jj < 5; jj++)
    {

        pionPairPtGt_TPCorToF[jj] = 1 + fabs(TPC_TPCandTOF_rel_diff_intPt_final);
        pionPairPtLt_TPCorToF[jj] = 1 + fabs(TPC_TPCandTOF_rel_diff_intPt_final_Lt);
        pionPairMGt_TPCorToF[jj] = 1 + fabs(TPC_TPCandTOF_rel_diff_intMinv_final);
        pionPairMLt_TPCorToF[jj] = 1 + fabs(TPC_TPCandTOF_rel_diff_intMinv_final_Lt);
        // cout << pionPairPtGt_TPCorToF[jj] << "\t PID sys from Relative difference of AUT between TOFandTPC and TPC sample for #eta>0\t" << endl;
        // cout << pionPairPtLt_TPCorToF[jj] << "\t PID sys from Relative difference of AUT between TOFandTPC and TPC sample for #eta<0\t" << endl;
        // cout << pionPairMGt_TPCorToF[jj] << "\t PID sys from Relative difference of AUT between TOFandTPC and TPC sample for #eta>0\t" << endl;
        // cout << pionPairMLt_TPCorToF[jj] << "\t PID sys from Relative difference of AUT between TOFandTPC and TPC sample for #eta<0\t" << endl;
    }

    double TPC_TPCandTOF_diff_intpT_Gt_p0_final = 0.000730;
    double TPC_TPCandTOF_diff_intpT_Lt_p0_final = 0.000066;
    double TPC_TPCandTOF_diff_intMinv_Gt_p0_final = 0.000437;
    double TPC_TPCandTOF_diff_intMinv_Lt_p0_final = 0.000195;
    double TPC_TPCandTOF_diff_Eta_final_p0_final = 0.000193;

    for (int kk = 0; kk < 9; kk++)
    {
        pionPairEta_TPCorToF[kk] = 1 + fabs(TPC_TPCandTOF_rel_Eta_diff_final);
    }

    // double pionPairPtGt_TOF[5] = {0.979031, 0.949769, 0.941655, 0.928794, 0.897969};
    // double pionPairPtLt_TOF[5] = {0.974835, 0.946528, 0.936604, 0.924064, 0.891258};
    // Final result
    double pionPairPtGt_ToF[5] = {0.944208, 0.917944, 0.912098, 0.899861, 0.859474};
    double pionPairPtLt_ToF[5] = {0.948695, 0.919545, 0.912388, 0.900562, 0.862505};
    double pionPairMGt_ToF[5] = {0.884668, 0.88769, 0.88976, 0.871818, 0.84386};
    double pionPairMLt_ToF[5] = {0.884153, 0.888802, 0.890719, 0.873281, 0.846215};
    double pionPairEta_ToF[9] = {0.902961, 0.886886, 0.865055, 0.843727, 0.840466, 0.85691, 0.874384, 0.891725, 0.901808};

    // the values in TrigBias_pTGt array are for A_UT vs pT plots for 5 different Minv Bins, eg TrigBias_pTGt[0] is for 0.2<Minv<0.48, TrigBias_pTGt[1] for 0.48<Minv<0.65 TrigBias_Gt[1] is for 0.65<Minv<0.84 and so on.
    // becare full, the name scheme is different than that of PID sysyematic//
    //===============================Assiging Trigger Bias Max between p0 or p0_err of Trigger Bias============================//
    double TrigBias_pTGt_p0[5] = {0.8857, 0.9317, 0.8950, 0.9661, 0.9451};
    double TrigBias_pTGt_p0_err[5] = {0.0951, 0.0841, 0.076, 0.0778, 0.0911};

    double TrigBias_pTLt_p0[5] = {1.05, 1.057, 0.9392, 0.9812, 0.9166};
    double TrigBias_pTLt_p0_err[5] = {0.09, 0.077, 0.0684, 0.0753, 0.0682};

    double TrigBias_MinvGt_p0[5] = {0.8906, 0.8057, 0.9887, 1.004, 0.9465};
    double TrigBias_MinvGt_p0_err[5] = {0.1826, 0.1086, 0.0978, 0.091, 0.0620};

    double TrigBias_MinvLt_p0[5] = {0.8803, 0.9203, 0.9651, 0.9992, 0.9665};
    double TrigBias_MinvLt_p0_err[5] = {0.1360, 0.1075, 0.0962, 0.0763, 0.0588};

    double TrigBias_pTGt[5] = {0};
    double TrigBias_pTLt[5] = {0};

    double TrigBias_MinvGt[5] = {0};
    double TrigBias_MinvLt[5] = {0};

    for (int ii = 0; ii < 5; ii++)
    {
        double pTGt_p0 = fabs(1 - TrigBias_pTGt_p0[ii]);
        double pTGt_p0_err = fabs(TrigBias_pTGt_p0_err[ii]);

        double pTLt_p0 = fabs(1 - TrigBias_pTLt_p0[ii]);
        double pTLt_p0_err = fabs(TrigBias_pTLt_p0_err[ii]);

        double MinvGt_p0 = fabs(1 - TrigBias_MinvGt_p0[ii]);
        double MinvGt_p0_err = fabs(TrigBias_MinvGt_p0_err[ii]);

        double MinvLt_p0 = fabs(1 - TrigBias_MinvLt_p0[ii]);
        double MinvLt_p0_err = fabs(TrigBias_MinvLt_p0_err[ii]);

        TrigBias_pTGt[ii] = (pTGt_p0 > pTGt_p0_err) ? TrigBias_pTGt_p0[ii] : 1 + pTGt_p0_err;
        TrigBias_pTLt[ii] = (pTLt_p0 > pTLt_p0_err) ? TrigBias_pTLt_p0[ii] : 1 + pTLt_p0_err;

        TrigBias_MinvGt[ii] = (MinvGt_p0 > MinvGt_p0_err) ? TrigBias_MinvGt_p0[ii] : 1 + MinvGt_p0_err;
        TrigBias_MinvLt[ii] = (MinvLt_p0 > MinvLt_p0_err) ? TrigBias_MinvLt_p0_err[ii] : 1 + MinvLt_p0_err;

        // cout << TrigBias_MinvGt[ii] << "\t Trig Bias for AUT vs Minv #eta>0 for pTbin" << ii << endl;
        // cout << TrigBias_MinvLt[ii] << "\t Trig Bias for AUT vs Minv #eta<0 for pTbin" << ii << endl;
        // cout << TrigBias_pTGt[ii] << "\t Trig Bias for AUT vs pT #eta>0 for Minvbin" << ii << endl;
        // cout << TrigBias_pTLt[ii] << "\t Trig Bias for AUT vs pT #eta<0 for Minvbin" << ii << endl;
    }

    double trigBias_eta_p0 = fabs(1 - 1.036);
    double trigBias_eta_p0_err = 0.076;
    double TrigBias_Eta = (trigBias_eta_p0 > trigBias_eta_p0_err) ? 1 + trigBias_eta_p0 : 1 + trigBias_eta_p0_err;
    // cout << TrigBias_Eta << "\t Trigger Bias Eta\t" << endl;

    //===============================Assiging Trigger Bias Max between p0 or p0_err of Trigger Bias Ends Here============================//

    double errXsys[9] = {0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15};
    // double errXsys[9] = {0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17};

    double pidsys_Gt[9] = {0};
    double pidsys_Lt[9] = {0};
    double Wpidsys_Gt[9] = {0};
    double Wpidsys_Lt[9] = {0};

    double WpidsysTPC_Gt[9] = {0};
    double WpidsysTPC_Lt[9] = {0};
    double WpidsysTPCorTOF_Gt[9] = {0};
    double WpidsysTPCorTOF_Lt[9] = {0};
    double WpidsysTOF_Gt[9] = {0};
    double WpidsysTOF_Lt[9] = {0};

    // Trigger Bias Sys
    double WTrigsys_Gt[9] = {0};
    double WTrigsys_Lt[9] = {0};
    // Total sys
    double WTotsys_Gt[9] = {0};
    double WTotsys_Lt[9] = {0};

    // ofstream AUT_Stat_PIDSys;
    // const char *outtxtfile = "AUT_Vs_pT_Stats_PIDSys_trigBias.txt";
    //
    // if (remove(outtxtfile) != 0)
    //{
    //    cout << "Unable to delete file" << endl;
    //}
    // else
    //{
    //    cout << outtxtfile << "\t file deleted\t" << endl;
    //}
    // AUT_Stat_PIDSys.open(outtxtfile, ios::app);

    ofstream AUT_Stat_LatexTable_Gt;
    ofstream AUT_Stat_LatexTable_Lt;
    const char *outtxtfile_Table_Gt = "AUT_Vs_pT_LatexTable_Gt.tex";
    if (remove(outtxtfile_Table_Gt) != 0)
    {
        cout << "Unable to delete file" << endl;
    }
    else
    {
        cout << outtxtfile_Table_Gt << "\t file deleted\t" << endl;
    }

    const char *outtxtfile_Table_Lt = "AUT_Vs_pT_LatexTable_Lt.tex";
    if (remove(outtxtfile_Table_Lt) != 0)
    {
        cout << "Unable to delete file" << endl;
    }
    else
    {
        cout << outtxtfile_Table_Lt << "\t file deleted\t" << endl;
    }
    // AUT_Stat_PIDSys.open(outtxtfile, ios::app);
    AUT_Stat_LatexTable_Gt.open(outtxtfile_Table_Gt, ios::app);
    AUT_Stat_LatexTable_Lt.open(outtxtfile_Table_Lt, ios::app);

    // Do make an arry for different graph and loop over this to make code short//
    // double pT_avg[9] = {0};
    for (int pad = 0; pad < 5; pad++)
    {
        if (pad == 0)
        {
            for (int mbin = 0; mbin < 9; mbin++)
            {
                AvgpT_G[mbin] = pTAvgGt[pad][mbin];
                AvgpT_L[mbin] = pTAvgLt[pad][mbin];
                pT_BG[mbin] = pTBGt[pad][mbin];
                pT_BL[mbin] = pTBLt[pad][mbin];
                pT_YG[mbin] = pTYGt[pad][mbin];
                pT_YL[mbin] = pTYLt[pad][mbin];
                // pT_avg[mbin] = pT1[mbin];
                AvgA_G[mbin] = avgA_pT1g[mbin];
                dA_G[mbin] = errA_pT1g[mbin];
                AvgA_L[mbin] = avgA_pT1l[mbin];
                dA_L[mbin] = errA_pT1l[mbin];

                WAvgA_G[mbin] = WavgA_pT1g[mbin];
                WdA_G[mbin] = WerrA_pT1g[mbin];
                WAvgA_L[mbin] = WavgA_pT1l[mbin];
                WdA_L[mbin] = WerrA_pT1l[mbin];

                A_BG[mbin] = A_pT1bg[mbin];
                dA_BG[mbin] = deltaA_pT1bg[mbin];
                A_BL[mbin] = A_pT1bl[mbin];
                dA_BL[mbin] = deltaA_pT1bl[mbin];
                A_YG[mbin] = A_pT1yg[mbin];
                dA_YG[mbin] = deltaA_pT1yg[mbin];
                A_YL[mbin] = A_pT1yl[mbin];
                dA_YL[mbin] = deltaA_pT1yl[mbin];

                pidsys_Gt[mbin] = fabs((1 - pionPairMGt[0]) * TMath::Max(fabs(AvgA_G[mbin]), fabs(dA_G[mbin])));
                pidsys_Lt[mbin] = fabs((1 - pionPairMLt[0]) * TMath::Max(fabs(AvgA_L[mbin]), fabs(dA_L[mbin])));
                Wpidsys_Gt[mbin] = fabs((1 - pionPairMGt[0]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin])));
                Wpidsys_Lt[mbin] = fabs((1 - pionPairMLt[0]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin])));
                // TPC
                WpidsysTPC_Gt[mbin] = fabs((1 - pionPairMGt_TPC[0]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin])));
                WpidsysTPC_Lt[mbin] = fabs((1 - pionPairMLt_TPC[0]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin])));
                // TPCorTOF
                WpidsysTPCorTOF_Gt[mbin] = fabs((1 - pionPairMGt_TPCorToF[0]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin])));
                WpidsysTPCorTOF_Lt[mbin] = fabs((1 - pionPairMLt_TPCorToF[0]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin])));
                // TOF
                WpidsysTOF_Gt[mbin] = fabs((1 - pionPairMGt_ToF[0]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin])));
                WpidsysTOF_Lt[mbin] = fabs((1 - pionPairMLt_ToF[0]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin])));

                // Trigger Bias
                WTrigsys_Gt[mbin] = fabs((1 - TrigBias_pTGt[0]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin])));
                WTrigsys_Lt[mbin] = fabs((1 - TrigBias_pTLt[0]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin])));

                // Total sys
                // WTotsys_Gt[mbin] = sqrt(pow(WpidsysTPCorTOF_Gt[mbin], 2) + pow(WTrigsys_Gt[mbin], 2));
                // WTotsys_Lt[mbin] = sqrt(pow(WpidsysTPCorTOF_Lt[mbin], 2) + pow(WTrigsys_Lt[mbin], 2));
                WTotsys_Gt[mbin] = sqrt(pow(TPC_TPCandTOF_diff_intMinv_Gt_p0_final, 2) + pow(WTrigsys_Gt[mbin], 2));
                WTotsys_Lt[mbin] = sqrt(pow(TPC_TPCandTOF_diff_intMinv_Lt_p0_final, 2) + pow(WTrigsys_Lt[mbin], 2));
            }
        }
        if (pad == 1)
        {
            for (int mbin = 0; mbin < 9; mbin++)
            {
                AvgpT_G[mbin] = pTAvgGt[pad][mbin];
                AvgpT_L[mbin] = pTAvgLt[pad][mbin];
                pT_BG[mbin] = pTBGt[pad][mbin];
                pT_BL[mbin] = pTBLt[pad][mbin];
                pT_YG[mbin] = pTYGt[pad][mbin];
                pT_YL[mbin] = pTYLt[pad][mbin];
                // pT_avg[mbin] = pT2[mbin];
                AvgA_G[mbin] = avgA_pT2g[mbin];
                dA_G[mbin] = errA_pT2g[mbin];
                AvgA_L[mbin] = avgA_pT2l[mbin];
                dA_L[mbin] = errA_pT2l[mbin];

                WAvgA_G[mbin] = WavgA_pT2g[mbin];
                WdA_G[mbin] = WerrA_pT2g[mbin];
                WAvgA_L[mbin] = WavgA_pT2l[mbin];
                WdA_L[mbin] = WerrA_pT2l[mbin];

                A_BG[mbin] = A_pT2bg[mbin];
                dA_BG[mbin] = deltaA_pT2bg[mbin];
                A_BL[mbin] = A_pT2bl[mbin];
                dA_BL[mbin] = deltaA_pT2bl[mbin];
                A_YG[mbin] = A_pT2yg[mbin];
                dA_YG[mbin] = deltaA_pT2yg[mbin];
                A_YL[mbin] = A_pT2yl[mbin];
                dA_YL[mbin] = deltaA_pT2yl[mbin];

                pidsys_Gt[mbin] = fabs((1 - pionPairMGt[1]) * TMath::Max(fabs(AvgA_G[mbin]), fabs(dA_G[mbin])));
                pidsys_Lt[mbin] = fabs((1 - pionPairMLt[1]) * TMath::Max(fabs(AvgA_L[mbin]), fabs(dA_L[mbin])));
                Wpidsys_Gt[mbin] = fabs((1 - pionPairMGt[1]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin])));
                Wpidsys_Lt[mbin] = fabs((1 - pionPairMLt[1]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin])));
                // TPC
                WpidsysTPC_Gt[mbin] = fabs((1 - pionPairMGt_TPC[1]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin])));
                WpidsysTPC_Lt[mbin] = fabs((1 - pionPairMLt_TPC[1]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin])));
                // TPCorTOF
                WpidsysTPCorTOF_Gt[mbin] = fabs((1 - pionPairMGt_TPCorToF[1]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin])));
                WpidsysTPCorTOF_Lt[mbin] = fabs((1 - pionPairMLt_TPCorToF[1]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin])));
                // TOF
                WpidsysTOF_Gt[mbin] = fabs((1 - pionPairMGt_ToF[1]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin])));
                WpidsysTOF_Lt[mbin] = fabs((1 - pionPairMLt_ToF[1]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin])));

                // Trigger Bias
                WTrigsys_Gt[mbin] = fabs((1 - TrigBias_pTGt[1]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin])));
                WTrigsys_Lt[mbin] = fabs((1 - TrigBias_pTLt[1]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin])));

                // Total sys
                // WTotsys_Gt[mbin] = sqrt(pow(WpidsysTPCorTOF_Gt[mbin], 2) + pow(WTrigsys_Gt[mbin], 2));
                // WTotsys_Lt[mbin] = sqrt(pow(WpidsysTPCorTOF_Lt[mbin], 2) + pow(WTrigsys_Lt[mbin], 2));
                WTotsys_Gt[mbin] = sqrt(pow(TPC_TPCandTOF_diff_intMinv_Gt_p0_final, 2) + pow(WTrigsys_Gt[mbin], 2));
                WTotsys_Lt[mbin] = sqrt(pow(TPC_TPCandTOF_diff_intMinv_Lt_p0_final, 2) + pow(WTrigsys_Lt[mbin], 2));
            }
        }
        if (pad == 2)
        {
            for (int mbin = 0; mbin < 9; mbin++)
            {
                AvgpT_G[mbin] = pTAvgGt[pad][mbin];
                AvgpT_L[mbin] = pTAvgLt[pad][mbin];
                pT_BG[mbin] = pTBGt[pad][mbin];
                pT_BL[mbin] = pTBLt[pad][mbin];
                pT_YG[mbin] = pTYGt[pad][mbin];
                pT_YL[mbin] = pTYLt[pad][mbin];

                // pT_avg[mbin] = pT3[mbin];
                AvgA_G[mbin] = avgA_pT3g[mbin];
                dA_G[mbin] = errA_pT3g[mbin];
                AvgA_L[mbin] = avgA_pT3l[mbin];
                dA_L[mbin] = errA_pT3l[mbin];

                WAvgA_G[mbin] = WavgA_pT3g[mbin];
                WdA_G[mbin] = WerrA_pT3g[mbin];
                WAvgA_L[mbin] = WavgA_pT3l[mbin];
                WdA_L[mbin] = WerrA_pT3l[mbin];

                A_BG[mbin] = A_pT3bg[mbin];
                dA_BG[mbin] = deltaA_pT3bg[mbin];
                A_BL[mbin] = A_pT3bl[mbin];
                dA_BL[mbin] = deltaA_pT3bl[mbin];
                A_YG[mbin] = A_pT3yg[mbin];
                dA_YG[mbin] = deltaA_pT3yg[mbin];
                A_YL[mbin] = A_pT3yl[mbin];
                dA_YL[mbin] = deltaA_pT3yl[mbin];

                pidsys_Gt[mbin] = fabs((1 - pionPairMGt[2]) * TMath::Max(fabs(AvgA_G[mbin]), fabs(dA_G[mbin])));
                pidsys_Lt[mbin] = fabs((1 - pionPairMLt[2]) * TMath::Max(fabs(AvgA_L[mbin]), fabs(dA_L[mbin])));
                Wpidsys_Gt[mbin] = fabs((1 - pionPairMGt[2]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin])));
                Wpidsys_Lt[mbin] = fabs((1 - pionPairMLt[2]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin])));
                // TPC
                WpidsysTPC_Gt[mbin] = fabs((1 - pionPairMGt_TPC[2]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin])));
                WpidsysTPC_Lt[mbin] = fabs((1 - pionPairMLt_TPC[2]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin])));
                // TPCorTOF
                WpidsysTPCorTOF_Gt[mbin] = fabs((1 - pionPairMGt_TPCorToF[2]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin])));
                WpidsysTPCorTOF_Lt[mbin] = fabs((1 - pionPairMLt_TPCorToF[2]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin])));
                // TOF
                WpidsysTOF_Gt[mbin] = fabs((1 - pionPairMGt_ToF[2]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin])));
                WpidsysTOF_Lt[mbin] = fabs((1 - pionPairMLt_ToF[2]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin])));

                // Trigger Bias
                WTrigsys_Gt[mbin] = fabs((1 - TrigBias_pTGt[2]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin])));
                WTrigsys_Lt[mbin] = fabs((1 - TrigBias_pTLt[2]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin])));

                // Total sys
                // WTotsys_Gt[mbin] = sqrt(pow(WpidsysTPCorTOF_Gt[mbin], 2) + pow(WTrigsys_Gt[mbin], 2));
                // WTotsys_Lt[mbin] = sqrt(pow(WpidsysTPCorTOF_Lt[mbin], 2) + pow(WTrigsys_Lt[mbin], 2));
                WTotsys_Gt[mbin] = sqrt(pow(TPC_TPCandTOF_diff_intMinv_Gt_p0_final, 2) + pow(WTrigsys_Gt[mbin], 2));
                WTotsys_Lt[mbin] = sqrt(pow(TPC_TPCandTOF_diff_intMinv_Lt_p0_final, 2) + pow(WTrigsys_Lt[mbin], 2));
            }
        }
        if (pad == 3)
        {
            for (int mbin = 0; mbin < 9; mbin++)
            {
                AvgpT_G[mbin] = pTAvgGt[pad][mbin];
                AvgpT_L[mbin] = pTAvgLt[pad][mbin];
                pT_BG[mbin] = pTBGt[pad][mbin];
                pT_BL[mbin] = pTBLt[pad][mbin];
                pT_YG[mbin] = pTYGt[pad][mbin];
                pT_YL[mbin] = pTYLt[pad][mbin];
                // pT_avg[mbin] = pT4[mbin];
                AvgA_G[mbin] = avgA_pT4g[mbin];
                dA_G[mbin] = errA_pT4g[mbin];
                AvgA_L[mbin] = avgA_pT4l[mbin];
                dA_L[mbin] = errA_pT4l[mbin];

                WAvgA_G[mbin] = WavgA_pT4g[mbin];
                WdA_G[mbin] = WerrA_pT4g[mbin];
                WAvgA_L[mbin] = WavgA_pT4l[mbin];
                WdA_L[mbin] = WerrA_pT4l[mbin];

                A_BG[mbin] = A_pT4bg[mbin];
                dA_BG[mbin] = deltaA_pT4bg[mbin];
                A_BL[mbin] = A_pT4bl[mbin];
                dA_BL[mbin] = deltaA_pT4bl[mbin];
                A_YG[mbin] = A_pT4yg[mbin];
                dA_YG[mbin] = deltaA_pT4yg[mbin];
                A_YL[mbin] = A_pT4yl[mbin];
                dA_YL[mbin] = deltaA_pT4yl[mbin];

                pidsys_Gt[mbin] = fabs((1 - pionPairMGt[3]) * TMath::Max(fabs(AvgA_G[mbin]), fabs(dA_G[mbin])));
                pidsys_Lt[mbin] = fabs((1 - pionPairMLt[3]) * TMath::Max(fabs(AvgA_L[mbin]), fabs(dA_L[mbin])));
                Wpidsys_Gt[mbin] = fabs((1 - pionPairMGt[3]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin])));
                Wpidsys_Lt[mbin] = fabs((1 - pionPairMLt[3]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin])));
                // TPC
                WpidsysTPC_Gt[mbin] = fabs((1 - pionPairMGt_TPC[3]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin])));
                WpidsysTPC_Lt[mbin] = fabs((1 - pionPairMLt_TPC[3]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin])));
                // TPCorTOF
                WpidsysTPCorTOF_Gt[mbin] = fabs((1 - pionPairMGt_TPCorToF[3]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin])));
                WpidsysTPCorTOF_Lt[mbin] = fabs((1 - pionPairMLt_TPCorToF[3]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin])));
                // TOF
                WpidsysTOF_Gt[mbin] = fabs((1 - pionPairMGt_ToF[3]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin])));
                WpidsysTOF_Lt[mbin] = fabs((1 - pionPairMLt_ToF[3]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin])));

                // Trigger Bias
                WTrigsys_Gt[mbin] = fabs((1 - TrigBias_pTGt[3]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin])));
                WTrigsys_Lt[mbin] = fabs((1 - TrigBias_pTLt[3]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin])));
                // Total sys
                // WTotsys_Gt[mbin] = sqrt(pow(WpidsysTPCorTOF_Gt[mbin], 2) + pow(WTrigsys_Gt[mbin], 2));
                // WTotsys_Lt[mbin] = sqrt(pow(WpidsysTPCorTOF_Lt[mbin], 2) + pow(WTrigsys_Lt[mbin], 2));
                WTotsys_Gt[mbin] = sqrt(pow(TPC_TPCandTOF_diff_intMinv_Gt_p0_final, 2) + pow(WTrigsys_Gt[mbin], 2));
                WTotsys_Lt[mbin] = sqrt(pow(TPC_TPCandTOF_diff_intMinv_Lt_p0_final, 2) + pow(WTrigsys_Lt[mbin], 2));
            }
        }
        if (pad == 4)
        {
            for (int mbin = 0; mbin < 9; mbin++)
            {
                AvgpT_G[mbin] = pTAvgGt[pad][mbin];
                AvgpT_L[mbin] = pTAvgLt[pad][mbin];
                pT_BG[mbin] = pTBGt[pad][mbin];
                pT_BL[mbin] = pTBLt[pad][mbin];
                pT_YG[mbin] = pTYGt[pad][mbin];
                pT_YL[mbin] = pTYLt[pad][mbin];
                // pT_avg[mbin] = pT5[mbin];
                AvgA_G[mbin] = avgA_pT5g[mbin];
                dA_G[mbin] = errA_pT5g[mbin];
                AvgA_L[mbin] = avgA_pT5l[mbin];
                dA_L[mbin] = errA_pT5l[mbin];

                WAvgA_G[mbin] = WavgA_pT5g[mbin];
                WdA_G[mbin] = WerrA_pT5g[mbin];
                WAvgA_L[mbin] = WavgA_pT5l[mbin];
                WdA_L[mbin] = WerrA_pT5l[mbin];

                A_BG[mbin] = A_pT5bg[mbin];
                dA_BG[mbin] = deltaA_pT5bg[mbin];
                A_BL[mbin] = A_pT5bl[mbin];
                dA_BL[mbin] = deltaA_pT5bl[mbin];
                A_YG[mbin] = A_pT5yg[mbin];
                dA_YG[mbin] = deltaA_pT5yg[mbin];
                A_YL[mbin] = A_pT5yl[mbin];
                dA_YL[mbin] = deltaA_pT5yl[mbin];

                pidsys_Gt[mbin] = fabs((1 - pionPairMGt[4]) * TMath::Max(fabs(AvgA_G[mbin]), fabs(dA_G[mbin])));
                pidsys_Lt[mbin] = fabs((1 - pionPairMLt[4]) * TMath::Max(fabs(AvgA_L[mbin]), fabs(dA_L[mbin])));
                Wpidsys_Gt[mbin] = fabs((1 - pionPairMGt[4]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin])));
                Wpidsys_Lt[mbin] = fabs((1 - pionPairMLt[4]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin])));
                // TPC
                WpidsysTPC_Gt[mbin] = fabs((1 - pionPairMGt_TPC[4]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin])));
                WpidsysTPC_Lt[mbin] = fabs((1 - pionPairMLt_TPC[4]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin])));
                // TPCorTOF
                WpidsysTPCorTOF_Gt[mbin] = fabs((1 - pionPairMGt_TPCorToF[4]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin])));
                WpidsysTPCorTOF_Lt[mbin] = fabs((1 - pionPairMLt_TPCorToF[4]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin])));
                // TOF
                WpidsysTOF_Gt[mbin] = fabs((1 - pionPairMGt_ToF[4]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin])));
                WpidsysTOF_Lt[mbin] = fabs((1 - pionPairMLt_ToF[4]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin])));

                // Trigger Bias
                WTrigsys_Gt[mbin] = fabs((1 - TrigBias_pTGt[4]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin])));
                WTrigsys_Lt[mbin] = fabs((1 - TrigBias_pTLt[4]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin])));

                // Total sys
                // WTotsys_Gt[mbin] = sqrt(pow(WpidsysTPCorTOF_Gt[mbin], 2) + pow(WTrigsys_Gt[mbin], 2));
                // WTotsys_Lt[mbin] = sqrt(pow(WpidsysTPCorTOF_Lt[mbin], 2) + pow(WTrigsys_Lt[mbin], 2));
                WTotsys_Gt[mbin] = sqrt(pow(TPC_TPCandTOF_diff_intMinv_Gt_p0_final, 2) + pow(WTrigsys_Gt[mbin], 2));
                WTotsys_Lt[mbin] = sqrt(pow(TPC_TPCandTOF_diff_intMinv_Lt_p0_final, 2) + pow(WTrigsys_Lt[mbin], 2));
            }
        }

        for (int BYA = 0; BYA < 3; BYA++)
        {
            for (int GL = 0; GL < 2; GL++)
            {
                for (int mbin = 0; mbin < 9; mbin++)
                {
                    if (BYA == 0 && GL == 0)
                    {
                        PT[BYA][GL].push_back(pT_BG[mbin]);
                        A_UT[BYA][GL].push_back(A_BG[mbin]);
                        dA_UT[BYA][GL].push_back(dA_BG[mbin]);
                    }
                    if (BYA == 0 && GL == 1)
                    {
                        PT[BYA][GL].push_back(pT_BL[mbin]);
                        A_UT[BYA][GL].push_back(A_BL[mbin]);
                        dA_UT[BYA][GL].push_back(dA_BL[mbin]);
                    }
                    if (BYA == 1 && GL == 0)
                    {
                        PT[BYA][GL].push_back(pT_YG[mbin]);
                        A_UT[BYA][GL].push_back(A_YG[mbin]);
                        dA_UT[BYA][GL].push_back(dA_YG[mbin]);
                    }
                    if (BYA == 1 && GL == 1)
                    {
                        PT[BYA][GL].push_back(pT_YL[mbin]);
                        A_UT[BYA][GL].push_back(A_YL[mbin]);
                        dA_UT[BYA][GL].push_back(dA_YL[mbin]);
                    }
                    if (BYA == 2 && GL == 0)
                    {
                        PT[BYA][GL].push_back(AvgpT_G[mbin]);
                        A_UT[BYA][GL].push_back(AvgA_G[mbin]);
                        dA_UT[BYA][GL].push_back(dA_G[mbin]);
                        PIDSYS[BYA][GL].push_back(pidsys_Gt[mbin]);
                        WA_UT[BYA][GL].push_back(WAvgA_G[mbin]);
                        WdA_UT[BYA][GL].push_back(WdA_G[mbin]);
                        WPIDSYS[BYA][GL].push_back(Wpidsys_Gt[mbin]);

                        WPIDSYS_TPC[BYA][GL].push_back(WpidsysTPC_Gt[mbin]);
                        WPIDSYS_TPCorTOF[BYA][GL].push_back(WpidsysTPCorTOF_Gt[mbin]);
                        WPIDSYS_TOF[BYA][GL].push_back(WpidsysTOF_Gt[mbin]);
                        // Total Systematic
                        WTOTSYS[BYA][GL].push_back(WTotsys_Gt[mbin]);
                    }
                    if (BYA == 2 && GL == 1)
                    {
                        PT[BYA][GL].push_back(AvgpT_L[mbin]);
                        A_UT[BYA][GL].push_back(AvgA_L[mbin]);
                        dA_UT[BYA][GL].push_back(dA_L[mbin]);
                        PIDSYS[BYA][GL].push_back(pidsys_Lt[mbin]);
                        WA_UT[BYA][GL].push_back(WAvgA_L[mbin]);
                        WdA_UT[BYA][GL].push_back(WdA_L[mbin]);
                        WPIDSYS[BYA][GL].push_back(Wpidsys_Lt[mbin]);

                        WPIDSYS_TPC[BYA][GL].push_back(WpidsysTPC_Lt[mbin]);
                        WPIDSYS_TPCorTOF[BYA][GL].push_back(WpidsysTPCorTOF_Lt[mbin]);
                        WPIDSYS_TOF[BYA][GL].push_back(WpidsysTOF_Lt[mbin]);
                        // Total Systematic
                        WTOTSYS[BYA][GL].push_back(WTotsys_Lt[mbin]);
                    }

                    // if(BYA==0&&GL==0){PT[BYA][GL].push_back(pT[mbin]); A_UT[BYA][GL].push_back(A_BG[mbin]); dA_UT[BYA][GL].push_back(dA_BG[mbin]);}
                    // if(BYA==0&&GL==1){PT[BYA][GL].push_back(pT[mbin]); A_UT[BYA][GL].push_back(A_BL[mbin]); dA_UT[BYA][GL].push_back(dA_BL[mbin]);}
                    // if(BYA==1&&GL==0){PT[BYA][GL].push_back(pT[mbin]); A_UT[BYA][GL].push_back(A_YG[mbin]); dA_UT[BYA][GL].push_back(dA_YG[mbin]);}
                    // if(BYA==1&&GL==1){PT[BYA][GL].push_back(pT[mbin]); A_UT[BYA][GL].push_back(A_YL[mbin]); dA_UT[BYA][GL].push_back(dA_YL[mbin]);}
                    // if(BYA==2&&GL==0){PT[BYA][GL].push_back(pT[mbin]); A_UT[BYA][GL].push_back(AvgA_G[mbin]); dA_UT[BYA][GL].push_back(dA_G[mbin]);}
                    // if(BYA==2&&GL==1){PT[BYA][GL].push_back(pT[mbin]); A_UT[BYA][GL].push_back(AvgA_L[mbin]); dA_UT[BYA][GL].push_back(dA_L[mbin]);}
                } // mbin loop ends here
                if (BYA == 2 && GL == 0)
                {
                    if (pad == 0)
                    {
                        AUT_Stat_LatexTable_Gt << "\\newgeometry{top=1cm,bottom=0.5cm}" << endl;
                        AUT_Stat_LatexTable_Gt << "\\section*{Table for $A_{UT}$ vs $p^{\\pi ^ +\\pi ^ -}_{T}$ for $\\eta^{\\pi ^ +\\pi ^ -}>0$}" << endl;
                        AUT_Stat_LatexTable_Gt << " \\begin{tabular}{| c | c | c | c | c | c |}" << endl;
                        AUT_Stat_LatexTable_Gt << "\\hline" << endl;
                        AUT_Stat_LatexTable_Gt << "\\textbf{$p^{\\pi ^ +\\pi ^ -}_{T}}$} &\\textbf{$A_{UT}$} &\\textbf{$\\sigma_{stat}$}&\\textbf{$\\sigma_{PID}$}  &\\textbf{$\\sigma_{Trig}$} &\\textbf{$\\sigma_{Tot}$}"
                                               << "\\"
                                               << "\\" << endl;
                        AUT_Stat_LatexTable_Gt << "\\hline" << endl;
                    }
                    if (pad == 1 || pad == 2 || pad == 3 || pad == 4)
                    {
                        AUT_Stat_LatexTable_Gt
                            << "\\hdashline" << endl;
                    }
                    for (int nBin = 0; nBin < 9; nBin++)
                    {
                        AUT_Stat_LatexTable_Gt << std::setprecision(2) << AvgpT_G[nBin] << "&" << std::setprecision(4) << WAvgA_G[nBin] << "&" << WdA_G[nBin] << "&" << fabs(TPC_TPCandTOF_diff_intMinv_Gt_p0_final) << "&" << fabs(WTrigsys_Gt[nBin]) << "&" << WTotsys_Gt[nBin] << "\\"
                                               << "\\" << endl;
                    }
                    if (pad == 4)
                    {
                        AUT_Stat_LatexTable_Gt
                            << "\\hline" << endl;
                        AUT_Stat_LatexTable_Gt << " \\end{tabular}" << endl;
                    }

                    // AUT_Stat_PIDSys << "=================================== Weighted Average Forward Dir =============================" << endl;
                    // AUT_Stat_PIDSys << "             " << endl;
                    // AUT_Stat_PIDSys << "pT-" << pad << "<WA_UT> = "
                    //                 << "{" << WAvgA_G[0] << "," << WAvgA_G[1] << "," << WAvgA_G[2] << "," << WAvgA_G[3] << "," << WAvgA_G[4] << "," << WAvgA_G[5] << "," << WAvgA_G[6] << "," << WAvgA_G[7] << "," << WAvgA_G[8] << "}" << endl;
                    // AUT_Stat_PIDSys << "pT-" << pad << "<WdAErr> = "
                    //                 << "{" << WdA_G[0] << "," << WdA_G[1] << "," << WdA_G[2] << "," << WdA_G[3] << "," << WdA_G[4] << "," << WdA_G[5] << "," << WdA_G[6] << "," << WdA_G[7] << "," << WdA_G[8] << "}" << endl;
                    // AUT_Stat_PIDSys << "pT-" << pad << "<WPidSys_TPCorTOF> = "
                    //                 << "{" << WpidsysTPCorTOF_Gt[0] << "," << WpidsysTPCorTOF_Gt[1] << "," << WpidsysTPCorTOF_Gt[2] << "," << WpidsysTPCorTOF_Gt[3] << "," << WpidsysTPCorTOF_Gt[4] << "," << WpidsysTPCorTOF_Gt[5] << "," << WpidsysTPCorTOF_Gt[6] << "," << WpidsysTPCorTOF_Gt[7] << "," << WpidsysTPCorTOF_Gt[8] << "}" << endl;
                    // AUT_Stat_PIDSys << "pT-" << pad << "<WTOTSys_TPCorTOF> = "
                    //                 << "{" << WTotsys_Gt[0] << "," << WTotsys_Gt[1] << "," << WTotsys_Gt[2] << "," << WTotsys_Gt[3] << "," << WTotsys_Gt[4] << "," << WTotsys_Gt[5] << "," << WTotsys_Gt[6] << "," << WTotsys_Gt[7] << "," << WTotsys_Gt[8] << "}" << endl;
                    // AUT_Stat_PIDSys << "             " << endl;
                    // AUT_Stat_PIDSys << "             " << endl;

                } // Average and Forward direction
                if (BYA == 2 && GL == 1)
                {

                    if (pad == 0)
                    {
                        AUT_Stat_LatexTable_Lt << "\\newgeometry{top=1cm,bottom=0.5cm}" << endl;
                        AUT_Stat_LatexTable_Lt << "\\section*{Table for $A_{UT}$ vs $p^{\\pi ^ +\\pi ^ -}_{T}$ for $\\eta^{\\pi ^ +\\pi ^ -}<0$}" << endl;
                        AUT_Stat_LatexTable_Lt << " \\begin{tabular}{| c | c | c | c | c | c |}" << endl;
                        AUT_Stat_LatexTable_Lt << "\\hline" << endl;
                        AUT_Stat_LatexTable_Lt << "\\textbf{$p^{\\pi ^ +\\pi ^ -}_{T}}$} &\\textbf{$A_{UT}$} &\\textbf{$\\sigma_{stat}$}&\\textbf{$\\sigma_{PID}$}  &\\textbf{$\\sigma_{Trig}$}  & \\textbf{$\\sigma_{Tot}$}"
                                               << "\\"
                                               << "\\" << endl;
                        AUT_Stat_LatexTable_Lt << "\\hline" << endl;
                    }
                    if (pad == 1 || pad == 2 || pad == 3 || pad == 4)
                    {
                        AUT_Stat_LatexTable_Lt
                            << "\\hdashline" << endl;
                    }
                    for (int nBin = 0; nBin < 9; nBin++)
                    {
                        AUT_Stat_LatexTable_Lt << std::setprecision(2) << AvgpT_L[nBin] << "&" << std::setprecision(4) << WAvgA_L[nBin] << "&" << WdA_L[nBin] << "&" << fabs(TPC_TPCandTOF_diff_intMinv_Lt_p0_final) << "&" << WTrigsys_Lt[nBin] << "&" << WTotsys_Lt[nBin] << "\\"
                                               << "\\" << endl;
                    }
                    if (pad == 4)
                    {
                        AUT_Stat_LatexTable_Lt
                            << "\\hline" << endl;
                        AUT_Stat_LatexTable_Lt << " \\end{tabular}" << endl;
                    }

                    // AUT_Stat_PIDSys << "=================================== Weighted Average Backward Dir =============================" << endl;
                    // AUT_Stat_PIDSys << "             " << endl;
                    // AUT_Stat_PIDSys << "pT-" << pad << "<WA_UT> = "
                    //                 << "{" << WAvgA_L[0] << "," << WAvgA_L[1] << "," << WAvgA_L[2] << "," << WAvgA_L[3] << "," << WAvgA_L[4] << "," << WAvgA_L[5] << "," << WAvgA_L[6] << "," << WAvgA_L[7] << "," << WAvgA_L[8] << "}" << endl;
                    // AUT_Stat_PIDSys << "pT-" << pad << "<WdAErr> = "
                    //                 << "{" << WdA_L[0] << "," << WdA_L[1] << "," << WdA_L[2] << "," << WdA_L[3] << "," << WdA_L[4] << "," << WdA_L[5] << "," << WdA_L[6] << "," << WdA_L[7] << "," << WdA_L[8] << "}" << endl;
                    // AUT_Stat_PIDSys << "pT-" << pad << "<WPidSys_TPCorTOF> = "
                    //                 << "{" << WpidsysTPCorTOF_Lt[0] << "," << WpidsysTPCorTOF_Lt[1] << "," << WpidsysTPCorTOF_Lt[2] << "," << WpidsysTPCorTOF_Lt[3] << "," << WpidsysTPCorTOF_Lt[4] << "," << WpidsysTPCorTOF_Lt[5]
                    //                 << "," << WpidsysTPCorTOF_Lt[6] << "," << WpidsysTPCorTOF_Lt[7] << "," << WpidsysTPCorTOF_Lt[8] << "}" << endl;
                    // AUT_Stat_PIDSys << "pT-" << pad << "<WTOTSys_TPCorTOF> = "
                    //                 << "{" << WTotsys_Lt[0] << "," << WTotsys_Lt[1] << "," << WTotsys_Lt[2] << "," << WTotsys_Lt[3] << "," << WTotsys_Lt[4] << "," << WTotsys_Lt[5] << "," << WTotsys_Lt[6] << "," << WTotsys_Lt[7] << "," << WTotsys_Lt[8] << "}" << endl;
                    // AUT_Stat_PIDSys << "             " << endl;
                    // AUT_Stat_PIDSys << "             " << endl;
                } // Average and Backward dir

            } // GL loops ends here
        }     // BYA loops ends here

        for (int BYA = 0; BYA < 3; BYA++)
        {
            for (int GL = 0; GL < 2; GL++)
            {
                //   cout << PT[BYA][GL].size() << " =pT size for "<< BYA << "\t BYA\t and  "<< GL << "\t GL\t "<< endl;
                //   cout << A_UT[BYA][GL].size() << " =A_UT size for "<< BYA << "\t BYA\t and  "<< GL << "\t GL\t "<< endl;
                //   cout << dA_UT[BYA][GL].size()<< " =dA_UT size for "<< BYA << "\t BYA\t and  "<< GL << "\t GL\t "<< endl;

                GRAPH[BYA][GL][pad] = new TGraphErrors(PT[BYA][GL].size(), &PT[BYA][GL][0], &A_UT[BYA][GL][0], 0, &dA_UT[BYA][GL][0]);
                if (BYA == 2)
                {
                    gr_sys[BYA][GL][pad] = new TGraphErrors(PT[BYA][GL].size(), &PT[BYA][GL][0], &A_UT[BYA][GL][0], errXsys, &PIDSYS[BYA][GL][0]);

                    Wgr_sys[BYA][GL][pad] = new TGraphErrors(PT[BYA][GL].size(), &PT[BYA][GL][0], &WA_UT[BYA][GL][0], errXsys, &WPIDSYS[BYA][GL][0]);
                    Wgr_sys_TPC[BYA][GL][pad] = new TGraphErrors(PT[BYA][GL].size(), &PT[BYA][GL][0], &WA_UT[BYA][GL][0], errXsys, &WPIDSYS_TPC[BYA][GL][0]);
                    Wgr_sys_TPCorTOF[BYA][GL][pad] = new TGraphErrors(PT[BYA][GL].size(), &PT[BYA][GL][0], &WA_UT[BYA][GL][0], errXsys, &WPIDSYS_TPCorTOF[BYA][GL][0]);
                    Wgr_sys_TOF[BYA][GL][pad] = new TGraphErrors(PT[BYA][GL].size(), &PT[BYA][GL][0], &WA_UT[BYA][GL][0], errXsys, &WPIDSYS_TOF[BYA][GL][0]);

                    // Wgr_sys[BYA][GL][pad] = new TGraphErrors(PT[BYA][GL].size(), &PT[BYA][GL][0], &WA_UT[BYA][GL][0], 0, 0);
                    Wgr_TOTSYS[BYA][GL][pad] = new TGraphErrors(PT[BYA][GL].size(), &PT[BYA][GL][0], &WA_UT[BYA][GL][0], errXsys, &WTOTSYS[BYA][GL][0]);
                    WGRAPH[BYA][GL][pad] = new TGraphErrors(PT[BYA][GL].size(), &PT[BYA][GL][0], &WA_UT[BYA][GL][0], 0, &WdA_UT[BYA][GL][0]);
                }
                // GRAPH[BYA][GL][pad] = new TGraphErrors(PT[BYA][GL].size(), &PT[BYA][GL][0], &A_UT[BYA][GL][0], 0, &dA_UT[BYA][GL][0]);
            }
        }

        for (int BYA = 0; BYA < 3; BYA++)
        {
            for (int GL = 0; GL < 2; GL++)
            {
                PT[BYA][GL].clear();
                A_UT[BYA][GL].clear();
                dA_UT[BYA][GL].clear();
                PIDSYS[BYA][GL].clear();
                WA_UT[BYA][GL].clear();
                WdA_UT[BYA][GL].clear();
                WPIDSYS[BYA][GL].clear();

                WPIDSYS_TPC[BYA][GL].clear();
                WPIDSYS_TPCorTOF[BYA][GL].clear();
                WPIDSYS_TOF[BYA][GL].clear();
                WTOTSYS[BYA][GL].clear();
            }
        }
        // cout << A_UT[BYA][GL].size() << "size of A_UT " << endl;
        // cout << dA_UT[BYA][GL].size() << "size of dA_UT" << endl;

        for (int BYA = 0; BYA < 3; BYA++)
        {
            for (int GL = 0; GL < 2; GL++)
            {
                GRAPH[BYA][GL][pad]->GetYaxis()->SetTitle("#font[22]{A_{UT}^{Sin(#Phi_{S}-#Phi_{R})}}");
                GRAPH[BYA][GL][pad]->GetYaxis()->CenterTitle(kTRUE);
                GRAPH[BYA][GL][pad]->SetTitle("");
                GRAPH[BYA][GL][pad]->GetYaxis()->SetTitleOffset(1.4);
                GRAPH[BYA][GL][pad]->GetYaxis()->SetLabelOffset(0.017);
                GRAPH[BYA][GL][pad]->GetYaxis()->SetTitleSize(0.06);
                GRAPH[BYA][GL][pad]->GetYaxis()->SetLabelSize(0.06);
                GRAPH[BYA][GL][pad]->GetYaxis()->SetLabelFont(22);
                GRAPH[BYA][GL][pad]->GetXaxis()->SetLabelFont(22);
                GRAPH[BYA][GL][pad]->GetXaxis()->SetTitle("#font[22]{p_{T}^{#pi^{+}#pi^{-}}(GeV/c)}");
                GRAPH[BYA][GL][pad]->GetXaxis()->SetTitleSize(0.06);
                GRAPH[BYA][GL][pad]->GetXaxis()->SetLabelSize(0.06);
                GRAPH[BYA][GL][pad]->GetYaxis()->SetRangeUser(-0.0252, 0.062);
                GRAPH[BYA][GL][pad]->GetXaxis()->SetLimits(2.45, 19.45);
                GRAPH[BYA][GL][pad]->GetXaxis()->SetNdivisions(505);
                GRAPH[BYA][GL][pad]->GetYaxis()->SetNdivisions(505);

                if (BYA == 2)
                {
                    // Wgr_sys[BYA][GL][pad]->GetYaxis()->SetTitle("#font[22]{A_{UT}^{Sin(#Phi_{S}-#Phi_{R})}}");
                    // Wgr_sys[BYA][GL][pad]->SetTitle("");
                    // Wgr_sys[BYA][GL][pad]->GetYaxis()->SetTitleOffset(1.4);
                    // Wgr_sys[BYA][GL][pad]->GetYaxis()->SetLabelOffset(0.017);
                    // Wgr_sys[BYA][GL][pad]->GetYaxis()->SetTitleSize(0.06);
                    // Wgr_sys[BYA][GL][pad]->GetYaxis()->SetLabelSize(0.06);
                    // Wgr_sys[BYA][GL][pad]->GetYaxis()->SetLabelFont(22);
                    // Wgr_sys[BYA][GL][pad]->GetXaxis()->SetLabelFont(22);
                    // Wgr_sys[BYA][GL][pad]->GetXaxis()->SetTitle("#font[22]{p_{T}^{#pi^{+}#pi^{-}}(GeV/c)}");
                    // Wgr_sys[BYA][GL][pad]->GetXaxis()->SetTitleSize(0.06);
                    // Wgr_sys[BYA][GL][pad]->GetXaxis()->SetLabelSize(0.06);
                    // Wgr_sys[BYA][GL][pad]->GetYaxis()->SetRangeUser(-0.0252, 0.062);
                    // Wgr_sys[BYA][GL][pad]->GetXaxis()->SetLimits(2.45, 19.45);
                    // Wgr_sys[BYA][GL][pad]->GetXaxis()->SetNdivisions(505);
                    // Wgr_sys[BYA][GL][pad]->GetYaxis()->SetNdivisions(505);
                    // Wgr_sys[BYA][GL][pad]->SetLineColor(1);
                    // Wgr_sys[BYA][GL][pad]->SetLineWidth(1);
                    // Wgr_sys[BYA][GL][pad]->SetFillStyle(0);
                    // PID systematics
                    // Wgr_sys_TPCorTOF[BYA][GL][pad]->GetYaxis()->SetTitle("#font[22]{A_{UT}^{Sin(#Phi_{S}-#Phi_{R})}}");
                    // Wgr_sys_TPCorTOF[BYA][GL][pad]->SetTitle("");
                    // Wgr_sys_TPCorTOF[BYA][GL][pad]->GetYaxis()->SetTitleOffset(1.4);
                    // Wgr_sys_TPCorTOF[BYA][GL][pad]->GetYaxis()->SetLabelOffset(0.017);
                    // Wgr_sys_TPCorTOF[BYA][GL][pad]->GetYaxis()->SetTitleSize(0.06);
                    // Wgr_sys_TPCorTOF[BYA][GL][pad]->GetYaxis()->SetLabelSize(0.06);
                    // Wgr_sys_TPCorTOF[BYA][GL][pad]->GetYaxis()->SetLabelFont(22);
                    // Wgr_sys_TPCorTOF[BYA][GL][pad]->GetXaxis()->SetLabelFont(22);
                    // Wgr_sys_TPCorTOF[BYA][GL][pad]->GetXaxis()->SetTitle("#font[22]{p_{T}^{#pi^{+}#pi^{-}}(GeV/c)}");
                    // Wgr_sys_TPCorTOF[BYA][GL][pad]->GetXaxis()->SetTitleSize(0.06);
                    // Wgr_sys_TPCorTOF[BYA][GL][pad]->GetXaxis()->SetLabelSize(0.06);
                    // Wgr_sys_TPCorTOF[BYA][GL][pad]->GetYaxis()->SetRangeUser(-0.0252, 0.062);
                    // Wgr_sys_TPCorTOF[BYA][GL][pad]->GetXaxis()->SetLimits(2.45, 19.45);
                    // Wgr_sys_TPCorTOF[BYA][GL][pad]->GetXaxis()->SetNdivisions(505);
                    // Wgr_sys_TPCorTOF[BYA][GL][pad]->GetYaxis()->SetNdivisions(505);
                    // Wgr_sys_TPCorTOF[BYA][GL][pad]->SetLineColor(1);
                    // Wgr_sys_TPCorTOF[BYA][GL][pad]->SetLineWidth(1);
                    // Wgr_sys_TPCorTOF[BYA][GL][pad]->SetFillStyle(0);
                    // Total Sytematic
                    Wgr_TOTSYS[BYA][GL][pad]->GetYaxis()->SetTitle("#font[22]{A_{UT}^{Sin(#Phi_{S}-#Phi_{R})}}");
                    Wgr_TOTSYS[BYA][GL][pad]->GetYaxis()->CenterTitle(kTRUE);
                    Wgr_TOTSYS[BYA][GL][pad]->SetTitle("");
                    Wgr_TOTSYS[BYA][GL][pad]->GetYaxis()->SetTitleOffset(1.4);
                    Wgr_TOTSYS[BYA][GL][pad]->GetYaxis()->SetLabelOffset(0.017);
                    Wgr_TOTSYS[BYA][GL][pad]->GetYaxis()->SetTitleSize(0.06);
                    Wgr_TOTSYS[BYA][GL][pad]->GetYaxis()->SetLabelSize(0.06);
                    Wgr_TOTSYS[BYA][GL][pad]->GetYaxis()->SetLabelFont(22);
                    Wgr_TOTSYS[BYA][GL][pad]->GetXaxis()->SetLabelFont(22);
                    Wgr_TOTSYS[BYA][GL][pad]->GetXaxis()->SetTitle("#font[22]{p_{T}^{#pi^{+}#pi^{-}}(GeV/c)}");
                    Wgr_TOTSYS[BYA][GL][pad]->GetXaxis()->SetTitleSize(0.06);
                    Wgr_TOTSYS[BYA][GL][pad]->GetXaxis()->SetLabelSize(0.06);
                    // Wgr_TOTSYS[BYA][GL][pad]->GetYaxis()->SetRangeUser(-0.0252, 0.062);
                    Wgr_TOTSYS[BYA][GL][pad]->GetYaxis()->SetRangeUser(-0.0122, 0.0452);
                    Wgr_TOTSYS[BYA][GL][pad]->GetXaxis()->SetLimits(2.45, 19.45);
                    Wgr_TOTSYS[BYA][GL][pad]->GetXaxis()->SetNdivisions(505);
                    Wgr_TOTSYS[BYA][GL][pad]->GetYaxis()->SetNdivisions(505);
                    Wgr_TOTSYS[BYA][GL][pad]->SetLineColor(1);
                    Wgr_TOTSYS[BYA][GL][pad]->SetLineWidth(1);
                    Wgr_TOTSYS[BYA][GL][pad]->SetFillStyle(0);

                    if (GL == 0)
                    {
                        // GL ==0 /eta>0, GL==1 /eta<0
                        GRAPH[BYA][GL][pad]->SetMarkerStyle(27);
                        GRAPH[BYA][GL][pad]->SetMarkerColor(1);
                        GRAPH[BYA][GL][pad]->SetLineColor(1);

                        WGRAPH[BYA][GL][pad]->SetMarkerStyle(27);
                        WGRAPH[BYA][GL][pad]->SetMarkerColor(1);
                        WGRAPH[BYA][GL][pad]->SetLineColor(1);
                        // GRAPH[BYA][GL][pad]->SetLineWidth(1);
                    }
                    if (GL == 1)
                    {
                        GRAPH[BYA][GL][pad]->SetMarkerStyle(20);
                        GRAPH[BYA][GL][pad]->SetMarkerColor(2);
                        GRAPH[BYA][GL][pad]->SetLineColor(2);

                        WGRAPH[BYA][GL][pad]->SetMarkerStyle(20);
                        WGRAPH[BYA][GL][pad]->SetMarkerColor(2);
                        WGRAPH[BYA][GL][pad]->SetLineColor(2);
                        // GRAPH[BYA][GL][pad]->SetLineWidth(1);
                    }
                } // If BYA==2
            }     // GL for loop
        }         // BYA for loop

    }                               // pad loop ends here
    AUT_Stat_LatexTable_Gt.close(); // Txt file with AUT_pT_PIDSys_TrigSys
    AUT_Stat_LatexTable_Lt.close(); // Txt file with AUT_pT_PIDSys_TrigSys

    // pT_boundary line;
    double b_line[5][10];
    for (int i = 0; i < 5; i++)
    {
        if (i == 0)
        {
            for (int j = 0; j < 10; j++)
            {
                b_line[i][j] = p1[j];
            }
        }
        if (i == 1)
        {
            for (int j = 0; j < 10; j++)
            {
                b_line[i][j] = p2[j];
            }
        }
        if (i == 2)
        {
            for (int j = 0; j < 10; j++)
            {
                b_line[i][j] = p3[j];
            }
        }
        if (i == 3)
        {
            for (int j = 0; j < 10; j++)
            {
                b_line[i][j] = p4[j];
            }
        }
        if (i == 4)
        {
            for (int j = 0; j < 10; j++)
            {
                b_line[i][j] = p5[j];
            }
        }
    }

    TLine *line[5]; // line passing through 0
    TLegend *leg;
    TLatex tex[5];

    TCanvas *myCanvA[3];
    const char *FBA_name[] = {"Forward", "BackWard", "Avg_BL"};

    // FBA= Forward(Blue and Yellow), BackWard(Blue and Yellow) and Averaged(Average of Blue and Yellow for forward and Backward)
    TLine *bound_line = new TLine();
    for (int FBA = 0; FBA < 3; FBA++)
    {
        myCanvA[FBA] = new TCanvas(Form("myCanvA_%s", FBA_name[FBA]), Form("myCanvA_%s", FBA_name[FBA]), 150, 10, 1100, 700);
        // gStyle->SetOptStat(0);
        myCanvA[FBA]->SetGrid(0, 0);
        myCanvA[FBA]->Divide(3, 2);
        // The 6th PAD is to write stuffs
        for (int PAD = 0; PAD < 6; PAD++)
        {
            // cout << FBA_name[FBA] << endl;
            // cout << PAD << "PAD Number"<< endl;
            if (PAD == 2)
            {
                myCanvA[FBA]->cd(PAD + 1);
                GRAPH[FBA][0][PAD]->GetYaxis()->SetTitle("");
                GRAPH[FBA][0][PAD]->GetYaxis()->SetLabelSize(0);
                GRAPH[FBA][1][PAD]->GetYaxis()->SetTitle("");
                GRAPH[FBA][1][PAD]->GetYaxis()->SetLabelSize(0);
                gPad->SetTopMargin(0.137);
                gPad->SetLeftMargin(0.02);
                gPad->SetBottomMargin(0.15);
                gPad->SetRightMargin(0.02);
                gPad->SetFillStyle(4000);
                gPad->SetPad(0.6634, 0.4197, 0.97, 0.9897);
            }

            if (PAD == 5)
            {
                myCanvA[FBA]->cd(PAD + 1);
                gPad->SetTopMargin(0.4);
                gPad->SetLeftMargin(0.01);
                gPad->SetBottomMargin(0.15);
                gPad->SetRightMargin(0.02);
                gPad->SetPad(0.6634, 0.01, 0.97, 0.49);
                gPad->SetFillStyle(4000);
                gPad->SetFrameFillColor(0);
                gPad->SetFrameFillStyle(0);
                gPad->SetFrameBorderMode(0);
                TLatex text;
                text.SetTextSize(0.08);
                text.DrawLatex(0.1, 0.7, "#font[22]{#color[2]{STAR Preliminary 2017}}");
                text.SetTextSize(0.045);
                text.DrawLatex(0.1, 0.6, "#font[22]{p^{#uparrow}+p #rightarrow #pi^{+}#pi^{-} + X at #sqrt{s} = 510 GeV, L_{int} = 350 pb^{-1}}");
                text.DrawLatex(0.1, 0.50, "#font[22]{#pm 1.4% scale uncertanity from beam}");
                text.DrawLatex(0.1, 0.45, "#font[22]{polarization (not shown)}");
            }

            if (PAD == 1)
            {
                myCanvA[FBA]->cd(PAD + 1);
                GRAPH[FBA][0][PAD]->GetYaxis()->SetTitle("");
                GRAPH[FBA][0][PAD]->GetYaxis()->SetLabelSize(0);
                GRAPH[FBA][1][PAD]->GetYaxis()->SetTitle("");
                GRAPH[FBA][1][PAD]->GetYaxis()->SetLabelSize(0);
                gPad->SetTopMargin(0.16);
                gPad->SetLeftMargin(0.0);
                gPad->SetBottomMargin(0.01);
                gPad->SetRightMargin(0.0);
                gPad->SetPad(0.37, 0.5, 0.67, 0.99);
            }

            if (PAD == 0)
            {
                myCanvA[FBA]->cd(PAD + 1);
                gPad->SetTopMargin(0.16);
                gPad->SetLeftMargin(0.2);
                gPad->SetBottomMargin(0.01);
                gPad->SetRightMargin(0.0);
                gPad->SetPad(0, 0.5, 0.37, 0.99);
            }

            if (PAD == 3)
            {
                myCanvA[FBA]->cd(PAD + 1);
                gPad->SetTopMargin(0.0);
                gPad->SetLeftMargin(0.2);
                gPad->SetBottomMargin(0.15);
                gPad->SetRightMargin(0.0);
                gPad->SetPad(0.0, 0.02052, 0.37, 0.5052);
            }

            if (PAD == 4)
            {
                myCanvA[FBA]->cd(PAD + 1);
                GRAPH[FBA][0][PAD]->GetYaxis()->SetTitle("");
                GRAPH[FBA][0][PAD]->GetYaxis()->SetLabelSize(0);
                GRAPH[FBA][1][PAD]->GetYaxis()->SetTitle("");
                GRAPH[FBA][1][PAD]->GetYaxis()->SetLabelSize(0);
                gPad->SetTopMargin(0.0);
                gPad->SetLeftMargin(0.0);
                gPad->SetBottomMargin(0.15);
                gPad->SetRightMargin(0.0035);
                gPad->SetPad(0.37, 0.02052, 0.67, 0.5052);
            }

            if (PAD < 5)
            {

                if (FBA == 0)
                {
                    // Plot Blue and Yello Eta_pair>0
                    // GRAPH[BYA][GL][pad]
                    GRAPH[0][0][PAD]->SetMarkerStyle(24);      // Blue Eta_pair>0;
                    GRAPH[0][0][PAD]->SetMarkerColor(4);       // Blue Eta_pair>0;
                    GRAPH[0][0][PAD]->SetLineColor(4);         // Blue Eta_pair>0;
                    GRAPH[1][0][PAD]->SetMarkerStyle(20);      // Yellow Eta_pair>0;
                    GRAPH[1][0][PAD]->SetMarkerColor(kYellow); // Yellow Eta_pair>0;
                    GRAPH[1][0][PAD]->SetLineColor(41);        // Yellow Eta_pair>0;
                    GRAPH[0][0][PAD]->Draw("AP");
                    GRAPH[1][0][PAD]->Draw("P same");

                    for (int kk = 0; kk < 10; kk++)
                    {
                        // bound_line->DrawLine(b_line[PAD][kk], 0.062 - 0.005, b_line[PAD][kk], 0.062);
                    }
                    if (PAD == 1)
                    {
                        leg = new TLegend(0.60, 0.70, 1.0, 0.84);
                        leg->AddEntry(GRAPH[0][0][PAD], "#font[22]{#eta^{#pi^{+}#pi^{-}}>0,BLUE}", "lp");
                        leg->AddEntry(GRAPH[1][0][PAD], "#font[22]{#eta^{#pi^{+}#pi^{-}}>0,YELLOW}", "lp");
                        leg->Draw();
                    }
                }
                if (FBA == 1)
                {
                    // GRAPH[BYA][GL][pad]
                    GRAPH[0][1][PAD]->SetMarkerStyle(24);      // Blue Eta_pair<0;
                    GRAPH[0][1][PAD]->SetMarkerColor(4);       // Blue Eta_pair<0;
                    GRAPH[0][1][PAD]->SetLineColor(4);         // Blue Eta_pair<0;
                    GRAPH[1][1][PAD]->SetMarkerStyle(20);      // Yellow Eta_pair<0;
                    GRAPH[1][1][PAD]->SetMarkerColor(kYellow); // Yellow Eta_pair<0;
                    GRAPH[1][1][PAD]->SetLineColor(41);        // Yellow Eta_pair<0;
                    GRAPH[0][1][PAD]->Draw("AP");
                    GRAPH[1][1][PAD]->Draw("P same");
                    for (int kk = 0; kk < 10; kk++)
                    {
                        // bound_line->DrawLine(b_line[PAD][kk], 0.062 - 0.005, b_line[PAD][kk], 0.062);
                    }
                    if (PAD == 1)
                    {

                        leg = new TLegend(0.60, 0.70, 1.0, 0.84);
                        leg->AddEntry(GRAPH[0][1][PAD], "#font[22]{#eta^{#pi^{+}#pi^{-}}<0,BLUE}", "lp");
                        leg->AddEntry(GRAPH[1][1][PAD], "#font[22]{#eta^{#pi^{+}#pi^{-}}<0,YELLOW}", "lp");
                        leg->Draw();
                    }
                }
                if (FBA == 2)
                {
                    // GRAPH[BYA][GL][pad]
                    // gr_sys[2][1][PAD]->Draw("AE2");   // BackWard
                    // GRAPH[2][1][PAD]->Draw("P same"); // Backward

                    // gr_sys[2][0][PAD]->Draw("AE2 same"); // Forward
                    // GRAPH[2][0][PAD]->Draw("P same");    // Forward

                    // Wgr_sys[2][1][PAD]->Draw("AE2 same"); // BackWard
                    // WGRAPH[2][1][PAD]->Draw("P same");    // Backward
                    // Wgr_sys[2][0][PAD]->SetFillStyle(0);
                    // Wgr_sys[2][0][PAD]->SetLineColor(3); // Green Color
                    // Wgr_sys[2][0][PAD]->SetLineWidth(2);
                    // Wgr_sys[2][0][PAD]->Draw("A2"); // Forward
                    // gPad->Update();

                    // Wgr_sys_TPC[2][0][PAD]->SetFillStyle(0);
                    // Wgr_sys_TPC[2][0][PAD]->SetLineColor(6); // pink color
                    // Wgr_sys_TPC[2][0][PAD]->SetLineWidth(2);
                    // Wgr_sys_TPC[2][0][PAD]->Draw("A2"); // Forward
                    //  Wgr_sys_TPC[2][0][PAD]->Draw("2 same"); // Forward
                    gPad->Update();
                    // Wgr_sys_TPCorTOF[2][0][PAD]->SetFillStyle(0);
                    // Wgr_sys_TPCorTOF[2][0][PAD]->SetLineColor(1); // Black color
                    // // Wgr_sys_TPCorTOF[2][0][PAD]->SetLineWidth(2);
                    // //  Wgr_sys_TPCorTOF[2][0][PAD]->Draw("2 same");
                    // Wgr_sys_TPCorTOF[2][0][PAD]->Draw("A2"); // Forward
                    //
                    // Wgr_sys_TPCorTOF[2][1][PAD]->SetFillStyle(0);
                    // Wgr_sys_TPCorTOF[2][1][PAD]->SetLineColor(2); // Black color
                    // // Wgr_sys_TPCorTOF[2][1][PAD]->SetLineWidth(2);
                    // //  Wgr_sys_TPCorTOF[2][0][PAD]->Draw("2 same");
                    // Wgr_sys_TPCorTOF[2][1][PAD]->Draw("2 same"); // Forward

                    // Apply Total Systematic
                    Wgr_TOTSYS[2][0][PAD]->SetFillStyle(0);
                    Wgr_TOTSYS[2][0][PAD]->SetLineColor(1); // Black color

                    Wgr_TOTSYS[2][0][PAD]->Draw("A2"); // Forward

                    Wgr_TOTSYS[2][1][PAD]->SetFillStyle(0);
                    Wgr_TOTSYS[2][1][PAD]->SetLineColor(1); // Red color

                    Wgr_TOTSYS[2][1][PAD]->Draw("2 same"); // Forward

                    // Wgr_sys_TOF[2][0][PAD]->SetFillStyle(0);
                    // Wgr_sys_TOF[2][0][PAD]->SetLineColor(kBlue); // Blue Color
                    // Wgr_sys_TOF[2][0][PAD]->SetLineWidth(2);
                    // Wgr_sys_TOF[2][0][PAD]->Draw("2 same");

                    WGRAPH[2][1][PAD]->Draw("P same"); // Backward
                    WGRAPH[2][0][PAD]->Draw("P same"); // Forward
                    // if (PAD < 4)
                    //{
                    //     GRAPH_R11[PAD]->SetMarkerStyle(3);
                    //     GRAPH_R11[PAD]->SetMarkerColor(4);
                    //     GRAPH_R11[PAD]->SetLineColor(4);
                    //     GRAPH_R11[PAD]->Draw("P same"); // Run 11
                    // }
                    GRAPH_JAMDiFF_G[PAD]->SetTitle("");
                    GRAPH_JAMDiFF_G[PAD]->SetFillStyle(1001);
                    GRAPH_JAMDiFF_G[PAD]->SetFillColorAlpha(8, 0.6);
                    GRAPH_JAMDiFF_G[PAD]->Draw("E3");
                    GRAPH_JAMDiFF_CentralVal_G[PAD]->Draw("LP");

                    GRAPH_JAMDiFF_L[PAD]->SetTitle("");
                    GRAPH_JAMDiFF_L[PAD]->SetFillStyle(1001);
                    GRAPH_JAMDiFF_L[PAD]->SetFillColorAlpha(9, 0.6);
                    GRAPH_JAMDiFF_L[PAD]->Draw("E3");
                    GRAPH_JAMDiFF_CentralVal_L[PAD]->Draw("LP");

                    for (int kk = 0; kk < 10; kk++)
                    {
                        // bound_line->DrawLine(b_line[PAD][kk], 0.0452 - 0.002, b_line[PAD][kk], 0.0452);
                    }

                    // GRAPH[2][0][PAD]->Draw("AP2");
                    // GRAPH[2][1][PAD]->Draw("same P");
                    if (PAD == 1)
                    {
                        leg = new TLegend(0.45, 0.75, 1.0, 0.82);
                        leg->AddEntry(WGRAPH[2][0][PAD], "#font[22]{#eta^{#pi^{+}#pi^{-}} > 0}", "lep");
                        // leg->AddEntry(GRAPH_R11[PAD], "Run11 #font[22]{#eta^{#pi^{+}#pi^{-}}>0}", "lep");
                        leg->AddEntry(WGRAPH[2][1][PAD], "#font[22]{#eta^{#pi^{+}#pi^{-}} < 0}", "lep");
                        leg->AddEntry(Wgr_TOTSYS[2][0][PAD], "#font[22]{Tot Sys.}", "f");
                        leg->SetNColumns(3);
                        // leg->AddEntry(GRAPH_R11[PAD], "Run11 #font[22]{#eta^{#pi^{+}#pi^{-}}>0}", "lp ");
                        // leg->AddEntry(Wgr_sys_TPC[2][0][PAD], "TPC");
                        // leg->AddEntry(Wgr_sys_TPCorTOF[2][0][PAD], "TPCorTOF");
                        // leg->AddEntry(Wgr_sys_TOF[2][0][PAD], "TOF");

                        TLegend *leg_JAMDiFF = new TLegend(0.45, 0.65, 1.0, 0.75);
                        leg_JAMDiFF->SetHeader("#font[22]{C. Cocuzza et. al. (JAMDiFF)}");
                        // leg->AddEntry((TObject *)0, "C. Cocuzza et. al. (JAMDiFF)", ""); // Placeholder entry for the text
                        leg_JAMDiFF->AddEntry(GRAPH_JAMDiFF_G[PAD], "#font[22]{#eta^{#pi^{+}#pi^{-}} > 0}", "f");
                        leg_JAMDiFF->AddEntry(GRAPH_JAMDiFF_L[PAD], "#font[22]{#eta^{#pi^{+}#pi^{-}} < 0}", "f");
                        leg_JAMDiFF->SetNColumns(2);

                        leg->SetBorderSize(0);
                        leg->Draw();
                        leg_JAMDiFF->SetBorderSize(0);
                        leg_JAMDiFF->Draw();
                    }
                }
                gPad->Update();
                line[PAD] = new TLine(myCanvA[FBA]->cd(PAD + 1)->GetUxmin(), 0, myCanvA[FBA]->cd(PAD + 1)->GetUxmax(), 0);
                line[PAD]->SetLineStyle(2);
                line[PAD]->SetLineWidth(1);
                line[PAD]->Draw();

                tex[PAD].SetTextSize(0.06);
                // tex[PAD].SetTextAlign(13);
                double Avg_Minv = std::round((avg_Minv[PAD]) * pow(10, 2)) / pow(10, 2);
                tex[PAD].DrawLatex(3, 0.028, Form("#color[1]{#font[22]{<M^{#pi^{+}#pi^{-}}_{inv}> = %g GeV/c^{2}}}", Avg_Minv));
                if (PAD == 1)
                {
                    // tex[PAD].DrawLatex(7, 0.048, "p_{T} bin boundires");
                }
                gPad->Update();
            } // PAD<5 control statement

        } // PAD Loop
        myCanvA[FBA]->Print(Form("./TriggerBias/AUT_Vs_pT_%s_Cone0.7_polCut_TPCorTOF.pdf", FBA_name[FBA]));
    } // FBA Loop ends here

    TCanvas *c_chi = new TCanvas("c_chi", "canv_chi", 900, 700);
    c_chi->Divide(2, 2);
    c_chi->cd(1);
    h_chi_BG->Draw();
    c_chi->cd(2);
    h_chi_BL->Draw();
    c_chi->cd(3);
    h_chi_YG->Draw();
    c_chi->cd(4);
    h_chi_YL->Draw();
    c_chi->Print("./TriggerBias/Chi_NDF_pT.pdf");
} // main
