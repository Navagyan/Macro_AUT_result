//////////////////////////////////////////////////////////////
// Edited : July 27 , 2021 Navagyan 				    //
//  This code produce histograms for TSSA as a function of Eta  //
//////////////////////////////////////////////////////////////
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

ofstream outfile;
ofstream output_LatexTable;
// void asymmetryVsEta_9bin(const char *ifile)
// void asymmetryVsEta_9bin(const char *ifile = "NTuple_P123_V1_July18.root")
// void asymmetryVsEta_9bin_trigBias_Working_LatexTable_TPCorTOF(const char *ifile = "../Ntuple_RawTree_polCut_TPCorTOF_P123.root")//
void asymmetryVsEta_9bin_TPCorTOF_final(const char *mode = "integrated_pT", const char *ifile = "../Ntuple_RawTree_polCut_TPCorTOF_P123.root")
{
    TH1D *hpolB = new TH1D("hpolB", "", 100, 0, 1);
    TH1D *hpolY = new TH1D("hpolY", "", 100, 0, 1);

    TH1D *h_chi_B = new TH1D("h_chi_B", "h_chi_B", 5, 0, 5);
    TH1D *h_chi_Y = new TH1D("h_chi_Y", "h_chi_Y", 5, 0, 5);

    TFile *f = new TFile(ifile);
    if (std::strcmp(mode, "integrated_pT") == 0)
    {
        outfile.open("AUT_Vs_Eta_trigBias_LatexTable_integrated_pT_TPCorTOF.txt");
        output_LatexTable.open("AUT_Vs_Eta_trigBias_LatexTable_integrated_pT_TPCorTOF.tex");
    }
    if (std::strcmp(mode, "highest_pT") == 0)
    {
        outfile.open("AUT_Vs_Eta_trigBias_LatexTable_highestpT.txt");
        output_LatexTable.open("AUT_Vs_Eta_trigBias_LatexTable_highestpT.tex");
    }

    //      get the trees
    TTree *ntuple1 = (TTree *)f->Get("ntuple1_TPCorTOF");
    TTree *ntuple2 = (TTree *)f->Get("ntuple2_TPCorTOF");
    TTree *ntuple4 = (TTree *)f->Get("ntuple4_TPCorTOF");
    TTree *ntuple5 = (TTree *)f->Get("ntuple5_TPCorTOF");
    // TTree *ntuple6 = (TTree *)f->Get("ntuple6_TPCorTOF");
    //  Define variables to hold leaves content
    float eta_pair;
    float PhiRS;
    float fspinconfig;
    float cone;
    float pT_pair;
    float Minv;
    float trigger;
    float PhiRSB;
    float PhiRSY;
    float fitPts_min_pair;
    float polB_corr, polY_corr;
    // To store average polarization values from histograms
    double avgPolB, avgPolY, rmsB, rmsY, avgPolT, rmsT;
    double pi = 3.14159265359;

    // variables for eta bin average
    double Npairs[9] = {0};
    double pTpairs[9] = {0};
    double etapairs[9] = {0};
    double Mpairs[9] = {0};
    double avg_eta[9] = {0};
    double avg_Minv[9] = {0};
    double avg_pT[9] = {0};
    double Npairs_Y[9] = {0};
    double pTpairs_Y[9] = {0};
    double etapairs_Y[9] = {0};
    double Mpairs_Y[9] = {0};
    double avg_eta_B[9] = {0};
    double avg_Minv_B[9] = {0};
    double avg_pT_B[9] = {0};
    double avg_eta_Y[9] = {0};
    double avg_Minv_Y[9] = {0};
    double avg_pT_Y[9] = {0};

    // Get the leaves content and store on variables
    // structure:: tree -> SetBanchAddress("variable to store leaf content", &Leaf)
    ntuple1->SetBranchAddress("fspinconfig", &fspinconfig);
    ntuple1->SetBranchAddress("cone", &cone);
    ntuple2->SetBranchAddress("Minv", &Minv);
    ntuple2->SetBranchAddress("pT_pair", &pT_pair);
    ntuple2->SetBranchAddress("eta_pair", &eta_pair);
    ntuple4->SetBranchAddress("PhiRSB", &PhiRSB);
    ntuple4->SetBranchAddress("PhiRSY", &PhiRSY);
    ntuple5->SetBranchAddress("fitPts_min_pair", &fitPts_min_pair);
    // ntuple6->SetBranchAddress("polB_corr", &polB_corr);

    // pT bin boundaries (boundaries are set to ensure roughlyequal stat in all bins)
    double pT[6] = {2.60, 4.628, 5.643, 6.979, 9.265, 25}; // Navagyan's Binnig for 90% of dataset
    double M[6] = {0.20, 0.4848, 0.6565, 0.8431, 1.1710, 4.000};
    // eta bin boundaries
    // double eta_range[10]={-1.200, -0.668, -0.469, -0.281, -0.096, 0.089, 0.275, 0.470, 0.675, 1.200};
    // double eta_range[10]={-1.200, -0.7310, -0.5575, -0.3490, -0.06350, 0.23725, 0.4583, 0.62350, 0.7694, 1.200};
    // double eta_range[10]={-1.200, -0.7202, -0.5444, -0.3190, 0.0115, 0.3033, 0.4951, 0.6402, 0.7722, 1.200};
    // double eta_range[10]={-1.200, -0.7243, -0.5438, -0.3183, 0.0145, 0.3055, 0.4961, 0.6402, 0.7717, 1.200};
    double eta_range[10] = {-1.200, -0.7243, -0.5539, -0.3410, -0.0225, 0.2833, 0.4852, 0.6346, 0.7695, 1.200};
    // double eta_range[10]={-1.200, -0.7210, -0.5442, -0.3185, 0.0110, 0.3042, 0.4981, 0.6446, 0.7775, 1.200};

    // Get number of entries in the tree
    int nentries = (int)ntuple1->GetEntries();
    cout << "nentries = " << nentries << endl;
    // add friend to the tree
    // structure :: tree_where_to_add -> AddFriend("tree_name_to_add", "root_file_where_the_tree_t0_add_is"
    // if no file name is given it is assumed that the friend tree is on the same root file
    ntuple1->AddFriend("ntuple2_TPCorTOF");
    ntuple1->AddFriend("ntuple4_TPCorTOF");
    ntuple1->AddFriend("ntuple5_TPCorTOF");
    int NpTGtUpB[175] = {0};
    int NpTLtUpB[175] = {0};
    int NpTGtDnB[175] = {0};
    int NpTLtDnB[175] = {0};
    int NMinvGtUpB[175] = {0};
    int NMinvLtUpB[175] = {0};
    int NMinvGtDnB[175] = {0};
    int NMinvLtDnB[175] = {0};
    int NetaUpB[175] = {0};
    int NetaDnB[175] = {0};
    int NpTGtUpY[175] = {0};
    int NpTLtUpY[175] = {0};
    int NpTGtDnY[175] = {0};
    int NpTLtDnY[175] = {0};
    int NMinvGtUpY[175] = {0};
    int NMinvLtUpY[175] = {0};
    int NMinvGtDnY[175] = {0};
    int NMinvLtDnY[175] = {0};
    int NetaUpY[175] = {0};
    int NetaDnY[175] = {0};
    int NpTGtUpT[175] = {0};
    int NpTLtUpT[175] = {0};
    int NpTGtDnT[175] = {0};
    int NpTLtDnT[175] = {0};
    int NMinvGtUpT[175] = {0};
    int NMinvLtUpT[175] = {0};
    int NMinvGtDnT[175] = {0};
    int NMinvLtDnT[175] = {0};
    int NetaUpT[175] = {0};
    int NetaDnT[175] = {0};

    for (int j = 0; j < nentries; j++)
    // for (int j = 0; j < 1000; j++)
    {
        ntuple1->GetEntry(j);
        if (cone >= .7)
            continue;
        // if (pT_pair < 3.75)
        //     continue;
        //	if(cone>=.3)continue;
        if (Minv > 4.)
            continue;
        if (fitPts_min_pair < 15)
            continue;
        //------------------for highest pT pair----------------
        if (std::strcmp(mode, "highest_pT") == 0)
        {
            if (pT_pair > 25 || pT_pair < 9.265)
                continue;
            // cout << pT_pair << endl;
        }
        else
        {
            // cout << "integrated pT mode  " << pT_pair << endl;
        }
        //---------for highest pT pair end--------

        //   if (pT_pair > 25 || pT_pair < 6.98)
        //      continue;
        //   if (pT_pair > 25 || pT_pair < 5.64)
        //      continue;
        //   if (pT_pair > 25 || pT_pair < 4.63)
        //      continue;
        //   if (pT_pair > 25 || pT_pair < 2.60)
        //      continue;
        //   if (Minv > 4.00 || Minv < 1.17) // Highest Minv Bin
        //       continue;
        //   if (Minv > 1.1710 || Minv < 0.8431) // 2nd Highest Minv Bin
        //      continue;
        //   if (Minv > 0.84 || Minv < 0.65) // 3nd Highest Minv Bin
        //      continue;
        //   if (Minv > 0.65 || Minv < 0.48) // 4th  Highest Minv Bin
        //      continue;
        //   if (Minv > 0.48 || Minv < 0.20) // 5th  Highest Minv Bin
        //   continue;
        //  if (Minv > 1.1710 || Minv < 0.65) // 3rd highest to 2nd highest Minv Bin
        //     continue;
        //  if (Minv > 1.1710 || Minv < 0.48) // 4th highest to 2nd highest Minv Bin
        //      continue;
        //    if (Minv > 1.1710 || Minv < 0.20) // 5th highest to 2nd highest Minv Bin
        //       continue;
        //    if (Minv > 4.00 || Minv < 0.20) // 5th highest to 1nd highest Minv Bin
        //       continue;
        //	if(Minv>0.4876 && Minv<0.5076) continue;//dodge K0 mass range, doesn't cause asymmetry.

        //-------------BLUE-------------------------------------
        // Phi
        for (int phi = 0; phi < 16; phi++)
        {
            if (PhiRSB >= (phi - 8.) / 8. * pi && PhiRSB <= (phi - 7.) / 8. * pi)
            {
                // Eta
                for (int eta = 0; eta < 9; eta++)
                {
                    if (eta_pair >= eta_range[eta] && eta_pair < eta_range[eta + 1])
                    {
                        Npairs[eta] = Npairs[eta] + 1;
                        pTpairs[eta] = pTpairs[eta] + pT_pair;
                        etapairs[eta] = etapairs[eta] + eta_pair;
                        Mpairs[eta] = Mpairs[eta] + Minv;
                        if (fspinconfig == 51 || fspinconfig == 53)
                        {
                            NetaUpB[eta * 16 + phi]++;
                        }
                        if (fspinconfig == 83 || fspinconfig == 85)
                        {
                            NetaDnB[eta * 16 + phi]++;
                        }
                    }
                } // Eta loop
            }     // Phi control
        }         // Phi loop

        //--------------end BLUE---------------------------------------

        //-----------------YELLOW---------------------------------------
        // Phi
        for (int phi = 0; phi < 16; phi++)
        {
            if (PhiRSY >= (phi - 8.) / 8. * pi && PhiRSY <= (phi - 7.) / 8. * pi)
            {
                // Eta
                for (int eta = 0; eta < 9; eta++)
                {
                    // Since beam circulating counter clock wise is color coded as the yellow beam which is in opposite direction to blue direction. So its eta direction is opposite to that of blue.
                    /// I made a change here
                    // double etapair = (-1) * eta_pair;
                    if (eta_pair >= eta_range[eta] && eta_pair < eta_range[eta + 1])
                    // if (etapair >= eta_range[eta] && etapair < eta_range[eta + 1])
                    {
                        Npairs_Y[eta] = Npairs_Y[eta] + 1;
                        pTpairs_Y[eta] = pTpairs_Y[eta] + pT_pair;
                        etapairs_Y[eta] = etapairs_Y[eta] + eta_pair;
                        // etapairs_Y[eta] = etapairs_Y[eta] + etapair;
                        Mpairs_Y[eta] = Mpairs_Y[eta] + Minv;
                        if (fspinconfig == 51 || fspinconfig == 83)
                        {
                            // NetaUpY[eta * 16 + phi]++; //if you choose to fill array with this method, you need to flip the etapair direction for Yellow, results are slightly different.
                            if (eta == 0)
                                NetaUpY[(eta + 8) * 16 + phi]++;
                            if (eta == 1)
                                NetaUpY[(eta + 6) * 16 + phi]++;
                            if (eta == 2)
                                NetaUpY[(eta + 4) * 16 + phi]++;
                            if (eta == 3)
                                NetaUpY[(eta + 2) * 16 + phi]++;
                            if (eta == 4)
                                NetaUpY[(eta)*16 + phi]++;
                            if (eta == 5)
                                NetaUpY[(eta - 2) * 16 + phi]++;
                            if (eta == 6)
                                NetaUpY[(eta - 4) * 16 + phi]++;
                            if (eta == 7)
                                NetaUpY[(eta - 6) * 16 + phi]++;
                            if (eta == 8)
                                NetaUpY[(eta - 8) * 16 + phi]++;
                        }
                        if (fspinconfig == 53 || fspinconfig == 85)
                        {
                            // NetaDnY[eta * 16 + phi]++;
                            if (eta == 0)
                                NetaDnY[(eta + 8) * 16 + phi]++;
                            if (eta == 1)
                                NetaDnY[(eta + 6) * 16 + phi]++;
                            if (eta == 2)
                                NetaDnY[(eta + 4) * 16 + phi]++;
                            if (eta == 3)
                                NetaDnY[(eta + 2) * 16 + phi]++;
                            if (eta == 4)
                                NetaDnY[(eta)*16 + phi]++;
                            if (eta == 5)
                                NetaDnY[(eta - 2) * 16 + phi]++;
                            if (eta == 6)
                                NetaDnY[(eta - 4) * 16 + phi]++;
                            if (eta == 7)
                                NetaDnY[(eta - 6) * 16 + phi]++;
                            if (eta == 8)
                                NetaDnY[(eta - 8) * 16 + phi]++;
                        }
                    }
                } // Eta loop
            }
        } // phi loop
        //---------------end YELLOW--------------------------------------
    }

    for (int k = 0; k < 9; k++)
    {
        avg_eta_B[k] = (double)etapairs[k] / (double)Npairs[k];
        avg_Minv_B[k] = (double)Mpairs[k] / (double)Npairs[k];
        avg_pT_B[k] = (double)pTpairs[k] / (double)Npairs[k];
        avg_eta_Y[k] = (double)etapairs_Y[k] / (double)Npairs_Y[k];
        avg_Minv_Y[k] = (double)Mpairs_Y[k] / (double)Npairs_Y[k];
        avg_pT_Y[k] = (double)pTpairs_Y[k] / (double)Npairs_Y[k];

        avg_eta[k] = 0.5 * (avg_eta_B[k] + avg_eta_Y[k]);
        avg_Minv[k] = 0.5 * (avg_Minv_B[k] + avg_Minv_Y[k]);
        avg_pT[k] = 0.5 * (avg_pT_B[k] + avg_pT_Y[k]);
        // cout<<"no. of pairs["<<k<<"]="<<Npairs[k]<<endl;
        // cout<<"pT of pairs["<<k<<"]="<<pTpairs[k]<<endl;
        // cout<<"no. of pairs["<<k<<"]="<<Npairs[k]<<endl;
        // cout<<"pT of pairs["<<k<<"]="<<pTpairs[k]<<endl;
        // cout<<"eta of pairs["<<k<<"]="<<etapairs[k]<<endl;
        // cout<<"<pT>["<<k<<"]="<<pTpairs[k]/Npairs[k]<<endl;
        // cout<<"<eta>["<<k<<"]="<<etapairs[k]/Npairs[k]<<endl;
        // cout<<"<Minv>["<<k<<"]="<<Mpairs[k]/Npairs[k]<<endl;
        if (outfile.is_open())
        {
            outfile << "Blue:=>"
                    << "avg_eta_B[" << k << "] = " << avg_eta_B[k] << ",  avg_Minv_B[" << k << "] = " << avg_Minv_B[k] << ", avg_pT_B[" << k << "] = " << avg_pT_B[k] << endl;
            outfile << "Yellow:=>"
                    << "avg_eta_Y[" << k << "] = " << avg_eta_Y[k] << ",  avg_Minv_Y[" << k << "] = " << avg_Minv_Y[k] << ", avg_pT_Y[" << k << "] = " << avg_pT_Y[k] << endl;
            outfile << "BY avg:=>"
                    << "avg_eta[" << k << "] = " << avg_eta[k] << ",  avg_Minv[" << k << "] = " << avg_Minv[k] << ", avg_pT[" << k << "] = " << avg_pT[k] << endl;
        }
    }

    //-----------------end total counts-----------------------------

    // Get polarization mean and rms from histograms
    // avgPolB = hpolB->GetMean();
    // avgPolY = hpolY->GetMean();
    // rmsB    = hpolB->GetRMS();
    // rmsY    = hpolY->GetRMS();
    avgPolB = 0.59;
    avgPolY = 0.5948;
    rmsB = 0.03148;
    rmsY = 0.03697;

    double a, b;
    double dAdB, dAdC, dAdD, dAdE, dAdP;
    double dP_B = rmsB; // use rms from total distribution for now. Should do independent calculatioon for blue and yellow then averaged.
    double dP_Y = rmsY;
    double dA[175], Asym[175]; // 16 array would be enough
    double B, C, D, E;
    // double pi = 3.14159265359;
    double BE, DC;
    // asymmetry for BLUE eta > 0
    double AB[9] = {0};
    double deltaAB[9] = {0};
    double AY[9] = {0};
    double deltaAY[9] = {0};

    TCanvas *canv_Blue = new TCanvas("Blue_Canv", "Blue_Canv", 900, 900);
    canv_Blue->Divide(3, 3);
    for (int eta = 0; eta < 9; eta++)
    {
        for (int ang = 0; ang < 16; ang++)
        {
            if (ang < 8)
            {
                B = NetaUpB[eta * 16 + ang];
                C = NetaUpB[eta * 16 + ang + 8];
                D = NetaDnB[eta * 16 + ang];
                E = NetaDnB[eta * 16 + ang + 8];
                BE = (double)(B * E);
                DC = (double)(D * C);
                a = sqrt(BE);
                b = sqrt(DC);
            }
            if (ang > 7)
            {
                B = NetaUpB[eta * 16 + ang];
                C = NetaUpB[eta * 16 + ang - 8];
                D = NetaDnB[eta * 16 + ang];
                E = NetaDnB[eta * 16 + ang - 8];
                BE = (double)(B * E);
                DC = (double)(D * C);
                a = sqrt(BE);
                b = sqrt(DC);
            }
            Asym[ang] = (1. / avgPolB) * ((a - b) / (a + b));
            dAdB = (1. / avgPolB) * (E * sqrt(DC)) / (sqrt(BE) * (pow((sqrt(BE) + sqrt(DC)), 2)));
            dAdE = (1. / avgPolB) * (B * sqrt(DC)) / (sqrt(BE) * (pow((sqrt(BE) + sqrt(DC)), 2)));
            dAdD = (-1. / avgPolB) * (C * sqrt(BE)) / (sqrt(DC) * (pow((sqrt(BE) + sqrt(DC)), 2)));
            dAdC = (-1. / avgPolB) * (D * sqrt(BE)) / (sqrt(DC) * (pow((sqrt(BE) + sqrt(DC)), 2)));
            dAdP = (-1. / (avgPolB * avgPolB)) * (sqrt(BE) - sqrt(DC)) / (sqrt(BE) + sqrt(DC));
            dA[ang] = sqrt(pow((fabs(dAdB) * sqrt(B)), 2) + pow((fabs(dAdC) * sqrt(C)), 2) + pow((fabs(dAdD) * sqrt(D)), 2) + pow((fabs(dAdE) * sqrt(E)), 2) + pow((fabs(dAdP) * dP_B), 2));

            // cout << "***************** eta bin: " << eta << ", phi bin: "<< ang << "  **********************"<< endl;
            // cout <<"B: "<< B << ", C: "<< C << ", D: "<< D << ", E: "<< E << ", a: "<< a << ", b: "<< b <<", BE: "<<BE << ", DC:  "<<DC<<endl;
            // cout << "dAdB: "<< dAdB <<  ", dAdE: "<< dAdE << ", dAdD: "<< dAdD << ", dAdC: "<< dAdC << ", dAdP:  "<< dAdP << ", Asym["<<ang << "]: "<<Asym[ang]<< ", dA["<<ang<<"]: "<< dA[ang]<<", avgPolB: "<< avgPolB << endl;
        } // angle loop

        char name[600];
        char title[600];
        double chi2Ndf[9];
        double chi2[9];
        double angle[16] = {-15. / 16. * pi, -13. / 16. * pi, -11. / 16. * pi, -9. / 16. * pi, -7. / 16. * pi, -5. / 16. * pi, -3. / 16. * pi, -1. / 16. * pi, pi / 16, 3 * pi / 16, 5 * pi / 16, 7 * pi / 16, 9 * pi / 16, 11 * pi / 16, 13 * pi / 16, 15 * pi / 16};
        canv_Blue->cd(eta + 1);
        double ex[16] = {0};
        gStyle->SetOptDate(0);
        gStyle->SetOptFit(1);
        auto grt = new TGraphErrors(16, angle, Asym, ex, dA);
        grt->SetMarkerStyle(20);
        sprintf(title, "eta bin %i, BLUE", eta);
        grt->SetTitle(title);
        grt->Draw("AP");
        TF1 *fit = new TF1("fit", "[0]*sin(x)", -3.14159265359, 0);

        fit->SetParameter(0, 0.0001);
        grt->Fit(fit, "RQ");
        grt->SetMarkerColor(4);
        grt->GetXaxis()->SetTitle("#Phi_{RS}");
        grt->GetYaxis()->SetRangeUser(-0.06, 0.06);
        grt->GetXaxis()->SetTitleOffset(1);

        chi2Ndf[eta] = fit->GetChisquare() / fit->GetNDF();
        chi2[eta] = fit->GetChisquare();
        cout << "chisquare: " << chi2Ndf[eta] << endl;
        grt->GetYaxis()->SetTitle("A_{UT}(#Phi_{RS})");
        grt->GetYaxis()->CenterTitle(kTRUE);
        grt->GetYaxis()->SetTitleOffset(1);
        AB[eta] = fit->GetParameter(0);
        deltaAB[eta] = fit->GetParError(0);
        // cout << AB[eta] << "\t Blue A_UT Value for \t" << eta << "\t Bin" << endl;
        // cout << deltaAB[eta] << "\t Blue delta A_UT Value for \t" << eta << "\t Bin" << endl;
        h_chi_B->Fill(chi2Ndf[eta]);
        canv_Blue->Update();
        // canv_Blue->Print("./TriggerBias/FitPlot_EtaBin_Blue_integratedpT_TPCorTOF.pdf");
        if (std::strcmp(mode, "highest_pT") == 0)
        {
            canv_Blue->Print("./TriggerBias/FitPlot_EtaBin_Blue_highestpT_TPCorTOF.pdf");
        }
        if (std::strcmp(mode, "integrated_pT") == 0)
        {
            canv_Blue->Print("./TriggerBias/FitPlot_EtaBin_Blue_integratedpT_TPCorTOF.pdf");
        }
        // sprintf(name, "eta_bin%i_Blue_9bin.png", eta);

    } // eta loop

    // asymmetry for yellow beam
    TCanvas *canv_Yellow = new TCanvas("Yellow_Beam", "Yellow_Beam", 900, 900);
    canv_Yellow->Divide(3, 3);
    for (int eta = 0; eta < 9; eta++)
    {
        for (int ang = 0; ang < 16; ang++)
        {
            if (ang < 8)
            {
                B = NetaUpY[eta * 16 + ang];
                C = NetaUpY[eta * 16 + ang + 8];
                D = NetaDnY[eta * 16 + ang];
                E = NetaDnY[eta * 16 + ang + 8];
                BE = (double)(B * E);
                DC = (double)(D * C);
                a = sqrt(BE);
                b = sqrt(DC);
            }
            if (ang > 7)
            {
                B = NetaUpY[eta * 16 + ang];
                C = NetaUpY[eta * 16 + ang - 8];
                D = NetaDnY[eta * 16 + ang];
                E = NetaDnY[eta * 16 + ang - 8];
                BE = (double)(B * E);
                DC = (double)(D * C);
                a = sqrt(BE);
                b = sqrt(DC);
            }
            Asym[ang] = (1. / avgPolB) * ((a - b) / (a + b));
            dAdB = (1. / avgPolB) * (E * sqrt(DC)) / (sqrt(BE) * (pow((sqrt(BE) + sqrt(DC)), 2)));
            dAdE = (1. / avgPolB) * (B * sqrt(DC)) / (sqrt(BE) * (pow((sqrt(BE) + sqrt(DC)), 2)));
            dAdD = (-1. / avgPolB) * (C * sqrt(BE)) / (sqrt(DC) * (pow((sqrt(BE) + sqrt(DC)), 2)));
            dAdC = (-1. / avgPolB) * (D * sqrt(BE)) / (sqrt(DC) * (pow((sqrt(BE) + sqrt(DC)), 2)));
            dAdP = (-1. / (avgPolB * avgPolB)) * (sqrt(BE) - sqrt(DC)) / (sqrt(BE) + sqrt(DC));
            dA[ang] = sqrt(pow((fabs(dAdB) * sqrt(B)), 2) + pow((fabs(dAdC) * sqrt(C)), 2) + pow((fabs(dAdD) * sqrt(D)), 2) + pow((fabs(dAdE) * sqrt(E)), 2) + pow((fabs(dAdP) * dP_Y), 2));

            // cout << "***************** eta bin: " << eta << ", phi bin: "<< ang << " Yellow Beam **********************"<< endl;
            // cout <<"B: "<< B << ", C: "<< C << ", D: "<< D << ", E: "<< E << ", a: "<< a << ", b: "<< b <<", BE: "<<BE << ", DC:  "<<DC<<endl;
            // cout << "dAdB: "<< dAdB <<  ", dAdE: "<< dAdE << ", dAdD: "<< dAdD << ", dAdC: "<< dAdC << ", dAdP:  "<< dAdP << ", Asym["<<ang << "]: "<<Asym[ang]<< ", dA["<<ang<<"]: "<< dA[ang]<<", avgPolB: "<< avgPolB << endl;

        } // angle loop
        canv_Yellow->cd(eta + 1);
        char name[600];
        char title[600];
        double chi2Ndf[9];
        double chi2[9];
        double angle[16] = {-15. / 16. * pi, -13. / 16. * pi, -11. / 16. * pi, -9. / 16. * pi, -7. / 16. * pi, -5. / 16. * pi, -3. / 16. * pi, -1. / 16. * pi, pi / 16, 3 * pi / 16, 5 * pi / 16, 7 * pi / 16, 9 * pi / 16, 11 * pi / 16, 13 * pi / 16, 15 * pi / 16};

        double ex[16] = {0};
        gStyle->SetOptDate(0);
        auto grt = new TGraphErrors(16, angle, Asym, ex, dA);
        sprintf(title, "eta bin %i, YELLOW", eta);
        grt->SetTitle(title);
        grt->SetMarkerStyle(20);
        grt->Draw("AP");
        TF1 *fit = new TF1("fit", "[0]*sin(x)", -3.14159265359, 0);

        fit->SetParameter(0, 0.0001);
        grt->Fit(fit, "R Q");
        grt->SetMarkerColor(4);
        grt->GetXaxis()->SetTitle("#Phi_{RS}");
        grt->GetXaxis()->SetTitleOffset(1);

        grt->GetYaxis()->SetRangeUser(-0.06, 0.06);
        chi2Ndf[eta] = fit->GetChisquare() / fit->GetNDF();
        chi2[eta] = fit->GetChisquare();
        // cout<<"chisquare: "<<chi2Ndf[eta]<<endl;
        grt->GetYaxis()->SetTitle("A_{UT}(#Phi_{RS})");
        grt->GetYaxis()->CenterTitle(kTRUE);
        grt->GetYaxis()->SetTitleOffset(1);
        AY[eta] = fit->GetParameter(0);
        deltaAY[eta] = fit->GetParError(0);
        cout << AY[eta] << "\t Yellow A_UT Value for \t" << eta << "\t Bin" << endl;
        cout << deltaAY[eta] << "\t Yellow delta A_UT Value for \t" << eta << "\t Bin" << endl;
        sprintf(name, "eta_bin%i_YELLOW_9bin.png", eta);
        // canv_yellow->Print(name);
        canv_Yellow->Update();
        h_chi_Y->Fill(chi2Ndf[eta]);
        if (std::strcmp(mode, "highest_pT") == 0)
        {
            canv_Yellow->Print("./TriggerBias/FitPlot_EtaBin_Yellow_highestpT_TPCorTOF.pdf");
        }
        if (std::strcmp(mode, "integrated_pT") == 0)
        {
            canv_Yellow->Print("./TriggerBias/FitPlot_EtaBin_Yellow_integratedpT_TPCorTOF.pdf");
        }

    } // eta loop

    double A_B[9] = {0};
    double deltaA_B[9] = {0};
    double A_Y[9] = {0};
    double deltaA_Y[9] = {0};

    for (int ii = 0; ii < 9; ii++)
    {
        A_B[ii] = AB[ii];
        deltaA_B[ii] = deltaAB[ii];
        A_Y[ii] = AY[ii];
        deltaA_Y[ii] = deltaAY[ii];
    }

    // double A_B[9] = {AB[0], AB[1], AB[2], AB[3], AB[4], AB[5], AB[6], AB[7], AB[8]};
    // double deltaA_B[9] = {deltaAB[0], deltaAB[1], deltaAB[2], deltaAB[3], deltaAB[4], deltaAB[5], deltaAB[6], deltaAB[7], deltaAB[8]};
    // double A_Y[9] = {AY[0], AY[1], AY[2], AY[3], AY[4], AY[5], AY[6], AY[7], AY[8]};
    // double deltaA_Y[9] = {deltaAY[0], deltaAY[1], deltaAY[2], deltaAY[3], deltaAY[4], deltaAY[5], deltaAY[6], deltaAY[7], deltaAY[8]};

    double avgA[9] = {0};
    double errA[9] = {0};
    // Weighted Average
    double WavgA[9] = {0};
    double WerrA[9] = {0};
    // calculate avaerage asymmetry and total error
    double SumpT = 0;
    double SumMinv = 0;
    for (Int_t l = 0; l < 9; l++)
    {
        avgA[l] = (A_B[l] + A_Y[l]) / 2.;
        errA[l] = .5 * sqrt(pow(deltaA_B[l], 2) + pow(deltaA_Y[l], 2));

        WavgA[l] = (A_B[l] * (1 / pow(deltaA_B[l], 2)) + A_Y[l] * (1 / pow(deltaA_Y[l], 2))) / ((1 / pow(deltaA_B[l], 2)) + (1 / pow(deltaA_Y[l], 2)));
        WerrA[l] = 1 / sqrt((1 / pow(deltaA_B[l], 2)) + (1 / pow(deltaA_Y[l], 2)));

        SumpT = SumpT + avg_pT[l];
        SumMinv = SumMinv + avg_Minv[l];
        // outfile << "avgA["<<l<<"]="<< avgA[l]<< endl;
        // outfile << "errA["<<l<<"]="<< errA[l]<< endl;
    }

    double Avg_pT = std::round(SumpT / 9);
    // set precission upto 2 decimal number
    double Avg_Minv = std::round((SumMinv / 9) * pow(10, 2)) / pow(10, 2);
    outfile << "                                                       " << endl;
    outfile << "                                                       " << endl;
    outfile << "                                                       " << endl;
    outfile << "/////////////////////Average Asymmetry (eta) Values/////////////////////" << endl;
    outfile << "<pT>=" << SumpT / 9 << endl;
    outfile << "<Minv>=" << SumMinv / 9 << endl;
    outfile << " avg <pT>=" << Avg_pT << endl;
    outfile << "avg <Minv>=" << Avg_Minv << endl;
    outfile << "avgEta[9]={" << avg_eta[0] << ", " << avg_eta[1] << "," << avg_eta[2] << "," << avg_eta[3] << "," << avg_eta[4] << avg_eta[5] << ", " << avg_eta[6] << "," << avg_eta[7] << ", " << avg_eta[8] << "};" << endl;
    outfile << "avgAsymEta[9]={" << avgA[0] << ", " << avgA[1] << "," << avgA[2] << "," << avgA[3] << "," << avgA[4] << ", " << avgA[5] << ", " << avgA[6] << ", " << avgA[7] << ", " << avgA[8] << "};" << endl;
    outfile << "avgErrEta[9]={" << errA[0] << ", " << errA[1] << "," << errA[2] << "," << errA[3] << "," << errA[4] << ", " << errA[5] << ", " << errA[6] << ", " << errA[7] << ", " << errA[8] << "};" << endl;
    outfile << "                                                       " << endl;
    outfile << "                                                       " << endl;
    outfile << "*******Weighted Average Asymmetry  Values **************" << endl;
    outfile << "WavgAsymEta[9]={" << WavgA[0] << ", " << WavgA[1] << "," << WavgA[2] << "," << WavgA[3] << "," << WavgA[4] << ", " << WavgA[5] << ", " << WavgA[6] << ", " << WavgA[7] << ", " << WavgA[8] << "};" << endl;
    outfile << "WavgErrEta[9]={" << WerrA[0] << ", " << WerrA[1] << "," << WerrA[2] << "," << WerrA[3] << "," << WerrA[4] << ", " << WerrA[5] << ", " << WerrA[6] << ", " << WerrA[7] << ", " << WerrA[8] << "};" << endl;
    outfile << "                                                       " << endl;
    outfile << "                                                       " << endl;

    // Run2006 asymmetry values for comparison
    // double asym06[4] = {-0.00020, .0075, .013, .031};
    // double err06[4] = {.0060, .0061, .0061, .0060};
    // double eta06[4] = {-0.75, -0.26, 0.26, 0.75};

    // Systematics
    // Kaon and Proton Seperate
    // double pionPairEta[9] = {0.712991, 0.845311, 0.784701, 0.763526, 0.762675, 0.777183, 0.790308, 0.85319, 0.766311};
    // Kaon and Proton Combine
    // P20ic.SL20c data TPC PID
    double pionPairEta[9] = {0.851266, 0.833897, 0.821624, 0.817317, 0.81957, 0.820086, 0.825542, 0.836172, 0.848973};
    double errXsys[9] = {0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03};

    // P20ic.SL22b data TPC PID
    // double pionPairEta_TPC[9] = {0.903331, 0.86843, 0.84591, 0.832521, 0.837576, 0.85384, 0.869885, 0.886617, 0.909103};
    // Final PID result
    double pionPairEta_TPC[9] = {0.862276, 0.826708, 0.797671, 0.77381, 0.771638, 0.790094, 0.808581, 0.831038, 0.856093};
    // double errXsys[9] = {0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03};
    //  P20ic.SL22b data TPCorTOF PID
    // double pionPairEta_TPCorTOF[9] = {0.921772, 0.895675, 0.875577, 0.861274, 0.863155, 0.879875, 0.896546, 0.910883, 0.92669};
    // Final PID result
    // double pionPairEta_TPCorTOF[9] = {0.888436, 0.861023, 0.834911, 0.809617, 0.807478, 0.8277, 0.846784, 0.866677, 0.884602};

    //================================== Assigining PID sys from Relative difference of A_UT ===========================//

    double pionPairPtGt_TPCorToF[5] = {0};
    double pionPairPtLt_TPCorToF[5] = {0};
    double pionPairMGt_TPCorToF[5] = {0};
    double pionPairMLt_TPCorToF[5] = {0};
    double pionPairEta_TPCorTOF[9] = {0};
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

    double TPC_TPCandTOF_rel_diff_intMinv_p0_Lt = -0.89911;
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
    }
    for (int kk = 0; kk < 9; kk++)
    {
        pionPairEta_TPCorTOF[kk] = 1 + fabs(TPC_TPCandTOF_rel_Eta_diff_final);
        cout << pionPairEta_TPCorTOF[kk] << "\t PID sys pionPairEta purity\t" << endl;
    }

    // double errXsys[9] = {0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03};
    // P20ic.SL22b TOF only PID
    // double pionPairEta_TOF[9] = {0.927569, 0.914154, 0.8984, 0.884785, 0.887742, 0.901717, 0.916273, 0.926812, 0.933777};
    // Final PID result
    double pionPairEta_TOF[9] = {0.902961, 0.886886, 0.865055, 0.843727, 0.840466, 0.85691, 0.874384, 0.891725, 0.901808};
    // double errXsys[9] = {0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03};

    double TPC_TPCandTOF_diff_intpT_Gt_p0_final = 0.000730;
    double TPC_TPCandTOF_diff_intpT_Lt_p0_final = 0.000066;
    double TPC_TPCandTOF_diff_intMinv_Gt_p0_final = 0.000437;
    double TPC_TPCandTOF_diff_intMinv_Lt_p0_final = 0.000195;
    double TPC_TPCandTOF_diff_Eta_final_p0_final = 0.000193;

    double pidsys[9] = {0};
    double Wpidsys[9] = {0};
    double pidsys_TPC[9] = {0};
    double Wpidsys_TPC[9] = {0};
    double pidsys_TPCorTOF[9] = {0};
    double Wpidsys_TPCorTOF[9] = {0};
    double pidsys_TOF[9] = {0};
    double Wpidsys_TOF[9] = {0};

    double trigBias_eta_p0 = fabs(1 - 1.036);
    double trigBias_eta_p0_err = 0.076;
    double TrigBias_Eta = (trigBias_eta_p0 > trigBias_eta_p0_err) ? 1 + trigBias_eta_p0 : 1 + trigBias_eta_p0_err;
    cout << TrigBias_Eta << "\t Trigger Bias Eta\t" << endl;

    // Trig Bias Sys
    double WTrigsys[9] = {0};
    // Total Sys
    double WTot_sys[9] = {0};

    for (int etabin = 0; etabin < 9; etabin++)
    {
        // pidsys[etabin] = (1 - pionPairEta[etabin]) * TMath::Max(avgA[etabin], errA[etabin]);
        Wpidsys[etabin] = (1 - pionPairEta[etabin]) * TMath::Max(fabs(WavgA[etabin]), fabs(WerrA[etabin]));
        pidsys_TPC[etabin] = (1 - pionPairEta_TPC[etabin]) * TMath::Max(fabs(avgA[etabin]), fabs(errA[etabin]));
        Wpidsys_TPC[etabin] = (1 - pionPairEta_TPC[etabin]) * TMath::Max(fabs(WavgA[etabin]), fabs(WerrA[etabin]));

        pidsys_TPCorTOF[etabin] = (1 - pionPairEta_TPCorTOF[etabin]) * TMath::Max(fabs(avgA[etabin]), fabs(errA[etabin]));
        Wpidsys_TPCorTOF[etabin] = (1 - pionPairEta_TPCorTOF[etabin]) * TMath::Max(fabs(WavgA[etabin]), fabs(WerrA[etabin]));

        pidsys_TOF[etabin] = (1 - pionPairEta_TOF[etabin]) * TMath::Max(fabs(avgA[etabin]), fabs(errA[etabin]));
        Wpidsys_TOF[etabin] = (1 - pionPairEta_TOF[etabin]) * TMath::Max(fabs(WavgA[etabin]), fabs(WerrA[etabin]));
        // WTrigsys[etabin] = 0.2 * TMath::Max(WavgA[etabin], WerrA[etabin]);
        WTrigsys[etabin] = fabs((1 - TrigBias_Eta) * TMath::Max(fabs(WavgA[etabin]), fabs(WerrA[etabin])));
        // WTot_sys[etabin] = sqrt(pow(Wpidsys_TPCorTOF[etabin], 2) + pow(WTrigsys[etabin], 2));
        WTot_sys[etabin] = sqrt(pow(TPC_TPCandTOF_diff_Eta_final_p0_final, 2) + pow(WTrigsys[etabin], 2));
    }

    // outfile << "WPIDSysErrEtaTOF[9]={" << Wpidsys_TPCorTOF[0] << ", " << Wpidsys_TPCorTOF[1] << "," << Wpidsys_TPCorTOF[2] << "," << Wpidsys_TPCorTOF[3] << "," << Wpidsys_TPCorTOF[4] << ", " << Wpidsys_TPCorTOF[5] << ", " << Wpidsys_TPCorTOF[6] << ", " << Wpidsys_TPCorTOF[7] << ", " << Wpidsys_TPCorTOF[8] << "};" << endl;

    outfile << "WPIDSysErrEtaTOF[9]={" << TPC_TPCandTOF_diff_Eta_final_p0_final << ", " << TPC_TPCandTOF_diff_Eta_final_p0_final << "," << TPC_TPCandTOF_diff_Eta_final_p0_final << "," << TPC_TPCandTOF_diff_Eta_final_p0_final << "," << TPC_TPCandTOF_diff_Eta_final_p0_final << ", " << TPC_TPCandTOF_diff_Eta_final_p0_final << ", " << TPC_TPCandTOF_diff_Eta_final_p0_final << ", " << TPC_TPCandTOF_diff_Eta_final_p0_final << ", " << TPC_TPCandTOF_diff_Eta_final_p0_final << "};" << endl;

    outfile << "                                                       " << endl;

    outfile << "WTrigBiasTOF[9]={" << WTrigsys[0] << ", " << WTrigsys[1] << "," << WTrigsys[2] << "," << WTrigsys[3] << "," << WTrigsys[4] << ", " << WTrigsys[5] << ", " << WTrigsys[6] << ", " << WTrigsys[7] << ", " << WTrigsys[8] << "};" << endl;

    outfile << "                                                       " << endl;
    outfile << " WTot_sys[9]={" << WTot_sys[0] << ", " << WTot_sys[1] << "," << WTot_sys[2] << "," << WTot_sys[3] << "," << WTot_sys[4] << ", " << WTot_sys[5] << ", " << WTot_sys[6] << ", " << WTot_sys[7] << ", " << WTot_sys[8] << "};" << endl;
    outfile << "                                                       " << endl;
    outfile << "                                                       " << endl;

    output_LatexTable << "\\newgeometry{top=1cm,bottom=0.5cm}" << endl;

    if (std::strcmp(mode, "highest_pT") == 0)
    {
        output_LatexTable << "\\section*{Table for $A_{UT}$ vs $\\eta^{\\pi ^ +\\pi ^ -}$ for $9<p^{\\pi^+\\pi^-}_{T}<25$ }" << endl;
    }
    if (std::strcmp(mode, "integrated_pT") == 0)
    {
        output_LatexTable << "\\section*{Table for $A_{UT}$ vs $\\eta^{\\pi ^ +\\pi ^ -}$ for integrated $p^{\\pi^+\\pi^-}_{T}$ }" << endl;
    }

    // output_LatexTable << "\\section*{Table for $A_{UT}$ vs $\\eta^{\\pi ^ +\\pi ^ -}$ for $9<p^{\\pi^+\\pi^-}_{T}<25$ }" << endl;
    //  output_LatexTable << "\\section*{Table for $A_{UT}$ vs $\\eta^{\\pi ^ +\\pi ^ -}$ for integrated $p^{\\pi^+\\pi^-}_{T}$ }" << endl;
    output_LatexTable << " \\begin{tabular}{| c | c | c | c | c|c|}" << endl;
    output_LatexTable << "\\hline" << endl;
    output_LatexTable << "\\textbf{$\\eta^{\\pi ^ +\\pi ^ -}$}} &\\textbf{$A_{UT}$} &\\textbf{$\\sigma_{stat}$}&\\textbf{$\\sigma_{PID}$}  &\\textbf{$\\sigma_{Trig}$} & \\textbf{$\\sigma_{Tot}$}"
                      << "\\"
                      << "\\" << endl;
    output_LatexTable << "\\hline" << endl;

    for (int nBin = 0; nBin < 9; nBin++)
    {
        output_LatexTable << std::setprecision(2) << avg_eta[nBin] << "&" << std::setprecision(4) << WavgA[nBin] << "&" << WerrA[nBin] << "&" << fabs(TPC_TPCandTOF_diff_Eta_final_p0_final) << "&" << fabs(WTrigsys[nBin]) << "&" << WTot_sys[nBin] << "\\"
                          << "\\" << endl;
    }
    output_LatexTable
        << "\\hline" << endl;
    output_LatexTable << " \\end{tabular}" << endl;

    output_LatexTable.close();

    gStyle->SetOptDate(0);
    gROOT->ForceStyle(0);

    TLine *bound_line = new TLine();

    // Blue and yelllow asymmetry
    TCanvas *myCan = new TCanvas("myCan", "myCan", 500, 450);
    gStyle->SetTitleX(.35);
    gStyle->SetLegendBorderSize(0);
    myCan->cd();
    myCan->SetGrid(0, 0);
    gPad->SetLeftMargin(.12);

    auto gr1 = new TGraphErrors(9, avg_eta, A_B, 0, deltaA_B);
    gr1->SetMarkerStyle(21);
    gr1->SetMarkerSize(1);
    gr1->SetTitle("A_{UT} vs #eta");
    // gr1->GetXaxis()->SetRangeUser(-1.22, 1.22);
    // gr1->GetXaxis()->SetLimits(-1.22, 1.22);
    gr1->GetXaxis()->SetLabelSize(0.05);
    gr1->GetYaxis()->SetRangeUser(-.0252, .050);
    // gr1->GetYaxis()->SetRangeUser(-.0252, .070);
    gr1->SetMarkerColor(kBlue);
    gr1->GetYaxis()->SetTitle("A_{UT}");
    gr1->GetXaxis()->SetTitle("#eta^{#pi^{+}#pi^{-}}");
    gr1->GetXaxis()->CenterTitle();
    gr1->GetYaxis()->CenterTitle();
    gr1->GetYaxis()->SetTitleOffset(1.43);
    gr1->GetXaxis()->SetTitleOffset(1.);
    gr1->GetYaxis()->SetLabelSize(0.035);
    gr1->GetXaxis()->SetLabelSize(0.035);
    gr1->GetYaxis()->SetTitleSize(0.04);
    gr1->GetXaxis()->SetTitleSize(0.04);
    gr1->GetXaxis()->SetNdivisions(505);
    gr1->GetYaxis()->SetNdivisions(505);
    gr1->SetMarkerStyle(20);
    gr1->Draw("AP");
    myCan->Update();

    auto gr2 = new TGraphErrors(9, avg_eta, A_Y, 0, deltaA_Y);
    gr2->SetMarkerStyle(20);
    gr2->SetMarkerColor(kRed);
    gr2->Draw("same P");

    for (int i = 0; i < 10; i++)
    {
        bound_line->DrawLine(eta_range[i], 0.05 - .002, eta_range[i], 0.05);
    }

    myCan->Update();

    TLine *line1 = new TLine(gPad->GetUxmin(), 0., gPad->GetUxmax(), 0.);
    line1->SetLineStyle(2);
    line1->Draw();
    TLegend *leg1 = new TLegend(0.2, .68, 0.5, 0.88);
    leg1->AddEntry("", " Cone < 0.7 ", "");
    leg1->AddEntry(gr1, " Blue Beam  ", "lp");
    leg1->AddEntry(gr2, " Yellow Beam  ", "lp");
    leg1->SetTextSize(0.04);
    leg1->Draw();
    myCan->Update();
    // myCan->SaveAs("AsymVsEta_ConeLt7_9bin.pdf");
    if (std::strcmp(mode, "highest_pT") == 0)
    {
        myCan->SaveAs("./TriggerBias/AsymVsEta_ConeLt7_9bin_highestpT_TPCorTOF.pdf");
    }
    if (std::strcmp(mode, "integrated_pT") == 0)
    {
        myCan->SaveAs("./TriggerBias/AsymVsEta_ConeLt7_9bin_integratedpT_TPCorTOF.pdf");
    }

    // average of yellow and blue asymmetry
    TCanvas *avg = new TCanvas("avg", "avg", 500, 450);
    gStyle->SetTitle("");
    gStyle->SetTitleX(.4);
    gStyle->SetLegendBorderSize(0);
    avg->cd();
    avg->SetGrid(0, 0);
    gPad->SetLeftMargin(.12);

    auto gr_sys = new TGraphErrors(9, avg_eta, avgA, errXsys, pidsys);
    auto grA = new TGraphErrors(9, avg_eta, avgA, 0, errA);

    auto Wgr_sys = new TGraphErrors(9, avg_eta, WavgA, errXsys, Wpidsys);

    auto Wgr_sys_TPC = new TGraphErrors(9, avg_eta, WavgA, errXsys, Wpidsys_TPC);
    auto Wgr_sys_TPCorTOF = new TGraphErrors(9, avg_eta, WavgA, errXsys, Wpidsys_TPCorTOF);
    auto Wgr_sys_TOF = new TGraphErrors(9, avg_eta, WavgA, errXsys, Wpidsys_TOF);
    //  Total systematics;
    auto Wgr_TOTSYS = new TGraphErrors(9, avg_eta, WavgA, errXsys, WTot_sys);

    auto WgrA = new TGraphErrors(9, avg_eta, WavgA, 0, WerrA);

    double avg_eta_Run11[7] = {-0.85, -0.57, -0.29, 0.00, 0.29, 0.57, 0.85};
    double avgA_Run11[7] = {0.007, 0.003, 0.012, 0.020, 0.013, 0.026, 0.028};
    double errA_Run11[7] = {0.005, 0.004, 0.004, 0.005, 0.006, 0.006, 0.007};
    auto gr_Run11 = new TGraphErrors(7, avg_eta_Run11, avgA_Run11, 0, errA_Run11);

    // Wgr_sys->SetTitle(" < A_{UT} > vs #eta");
    //  Wgr_sys->GetXaxis()->SetRangeUser(-1.22, 1.22);
    //  Wgr_sys->GetXaxis()->SetLimits(-1.22, 1.22);
    // Wgr_sys_TPCorTOF->GetXaxis()->SetLabelSize(0.05);
    // Wgr_sys_TPCorTOF->GetYaxis()->SetTitle("#font[22]{A_{UT}^{Sin(#Phi_{S}-#Phi_{R})}}");
    // Wgr_sys_TPCorTOF->GetXaxis()->SetTitle("#eta^{#pi^{+}#pi^{-}}");
    // Wgr_sys_TPCorTOF->SetTitle("");
    // Wgr_sys_TPCorTOF->GetXaxis()->CenterTitle();
    // Wgr_sys_TPCorTOF->GetYaxis()->CenterTitle();
    // Wgr_sys_TPCorTOF->GetYaxis()->SetTitleOffset(1.43);
    // Wgr_sys_TPCorTOF->GetXaxis()->SetTitleOffset(1.);
    // Wgr_sys_TPCorTOF->GetYaxis()->SetLabelSize(0.035);
    // Wgr_sys_TPCorTOF->GetXaxis()->SetLabelSize(0.035);
    // Wgr_sys_TPCorTOF->GetYaxis()->SetTitleSize(0.04);
    // Wgr_sys_TPCorTOF->GetXaxis()->SetTitleSize(0.04);
    // Wgr_sys_TPCorTOF->GetYaxis()->SetRangeUser(-.0252, .050);
    // Wgr_sys_TPCorTOF->SetFillStyle(0);
    // Wgr_sys_TPCorTOF->SetLineColor(1);
    // Wgr_sys_TPCorTOF->SetLineWidth(1);

    Wgr_TOTSYS->GetXaxis()->SetLabelSize(0.05);
    Wgr_TOTSYS->GetYaxis()->SetTitle("#font[22]{A_{UT}^{Sin(#Phi_{S}-#Phi_{R})}}");
    Wgr_TOTSYS->GetYaxis()->CenterTitle(kTRUE);
    Wgr_TOTSYS->GetXaxis()->SetTitle("#eta^{#pi^{+}#pi^{-}}");
    Wgr_TOTSYS->SetTitle("");
    Wgr_TOTSYS->GetXaxis()->CenterTitle();
    Wgr_TOTSYS->GetYaxis()->SetTitleOffset(1.43);
    Wgr_TOTSYS->GetXaxis()->SetTitleOffset(1.);
    Wgr_TOTSYS->GetYaxis()->SetLabelSize(0.035);
    Wgr_TOTSYS->GetXaxis()->SetLabelSize(0.035);
    Wgr_TOTSYS->GetYaxis()->SetTitleSize(0.04);
    Wgr_TOTSYS->GetXaxis()->SetTitleSize(0.04);
    // For JAM and Run11 Calculation;
    Wgr_TOTSYS->GetYaxis()->SetRangeUser(-0.0048, .070);
    // This is preliminanty
    // Wgr_TOTSYS->GetYaxis()->SetRangeUser(-0.0048, .030);
    Wgr_TOTSYS->GetXaxis()->SetNdivisions(505);
    Wgr_TOTSYS->GetYaxis()->SetNdivisions(505);
    Wgr_TOTSYS->SetFillStyle(0);
    Wgr_TOTSYS->SetLineColor(1);
    Wgr_TOTSYS->SetLineWidth(1);

    // Wgr_sys_TPCorTOF->SetFillStyle(0);
    // Wgr_sys_TPCorTOF->SetLineColor(6);
    // Wgr_sys_TPCorTOF->SetLineWidth(1);

    // Wgr_sys_TPCorTOF->Draw("AE2");
    Wgr_TOTSYS->Draw("AE2");

    //------------------------ JAM theory Calc------------------------
    double JAM_Eta510_int[9] = {-0.83, -0.64, -0.45, -0.19, 0.14, 0.39, 0.56, 0.7, 0.86};
    double JAM_A510_int[9] = {0.001176242, 0.001684061, 0.00315167, 0.005069976, 0.008783057, 0.012789528, 0.015910953, 0.01505275, 0.01898724};
    double JAM_E510_int[9] = {0.00048114, 0.000570162, 0.000954386, 0.001363856, 0.002408741, 0.003543845, 0.00439211, 0.004030387, 0.004830157};

    double JAM_Eta510_high[9] = {-0.83, -0.64, -0.45, -0.19, 0.13, 0.39, 0.56, 0.7, 0.86};
    double JAM_A510_high[9] = {0.002342537, 0.003563637, 0.005341299, 0.008929921, 0.016462086, 0.026441394, 0.036377589, 0.046075776, 0.058210743};
    double JAM_E510_high[9] = {0.000497216, 0.000773359, 0.001212715, 0.002234751, 0.003980716, 0.005795919, 0.006592458, 0.007745962, 0.008092689};

    TGraphErrors *JAM_th510;
    TGraph *JAM_th510_CentralVal;
    if (std::strcmp(mode, "highest_pT") == 0)
    {
        JAM_th510 = new TGraphErrors(9, JAM_Eta510_high, JAM_A510_high, 0, JAM_E510_high);
        JAM_th510_CentralVal = new TGraph(9, JAM_Eta510_high, JAM_A510_high);
    }
    if (std::strcmp(mode, "integrated_pT") == 0)
    {
        JAM_th510 = new TGraphErrors(9, JAM_Eta510_int, JAM_A510_int, 0, JAM_E510_int);
        JAM_th510_CentralVal = new TGraph(9, JAM_Eta510_int, JAM_A510_int);
    }
    JAM_th510->SetTitle("");
    JAM_th510->SetFillStyle(1001);
    JAM_th510->SetFillColorAlpha(8, 0.6);
    JAM_th510->Draw("E3");
    JAM_th510_CentralVal->Draw("LP");
    //------------------------ JAM theory Calc End------------------------

    // Wgr_sys_TPC->GetXaxis()->SetLabelSize(0.05);
    // Wgr_sys_TPC->GetYaxis()->SetTitle("#font[22]{A_{UT}^{Sin(#Phi_{S}-#Phi_{R})}}");
    // Wgr_sys_TPC->GetXaxis()->SetTitle("#eta^{#pi^{+}#pi^{-}}");
    // Wgr_sys_TPC->GetXaxis()->CenterTitle();
    // Wgr_sys_TPC->GetYaxis()->CenterTitle();
    // Wgr_sys_TPC->GetYaxis()->SetTitleOffset(1.43);
    // Wgr_sys_TPC->GetXaxis()->SetTitleOffset(1.);
    // Wgr_sys_TPC->GetYaxis()->SetLabelSize(0.035);
    // Wgr_sys_TPC->GetXaxis()->SetLabelSize(0.035);
    // Wgr_sys_TPC->GetYaxis()->SetTitleSize(0.04);
    // Wgr_sys_TPC->GetXaxis()->SetTitleSize(0.04);
    // Wgr_sys_TPC->GetYaxis()->SetRangeUser(-.0252, .050);
    //  Wgr_sys_TPC->SetFillStyle(0);
    //  Wgr_sys_TPC->SetLineColor(6);
    //  Wgr_sys_TPC->SetLineWidth(1);

    // Wgr_sys_TPC->SetFillStyle(0);
    // Wgr_sys_TPC->SetLineColor(6);
    //// Wgr_sys_TPC->SetLineWidth(1);
    // Wgr_sys_TPC->Draw("A2");
    // gPad->Update();

    // Wgr_sys_TPCorTOF->SetFillStyle(3002);
    // Wgr_sys_TPCorTOF->SetFillColor(kRed);
    // Wgr_sys_TPCorTOF->Draw("2 same");
    //// gPad->Update();
    // Wgr_sys_TOF->SetFillStyle(3001);
    // Wgr_sys_TOF->SetFillColor(8);
    // Wgr_sys_TOF->Draw("2 same");

    // Wgr_sys->SetFillStyle(0);
    // Wgr_sys->SetLineColor(3); // Green Color
    // Wgr_sys->SetLineWidth(2);
    // Wgr_sys->Draw("A2"); // Forward
    // gPad->Update();
    // Wgr_sys_TPC->SetFillStyle(0);
    // Wgr_sys_TPC->SetLineColor(6); // pink color
    // Wgr_sys_TPC->SetLineWidth(2);
    // Wgr_sys_TPC->Draw("A2"); // Forward
    // Wgr_sys_TPC->Draw("2 same"); // Forward
    //  gPad->Update();

    // Wgr_sys_TPCorTOF->SetFillStyle(0);
    // Wgr_sys_TPCorTOF->SetLineColor(1); // Cyan colo
    //// Wgr_sys_TPCorTOF->SetLineWidth(2);
    ////  Wgr_sys_TPCorTOF->Draw("2 same");
    // Wgr_sys_TPCorTOF->Draw("AE2");

    // Wgr_sys_TOF->SetFillStyle(0);
    // Wgr_sys_TOF->SetLineColor(kBlue); // Blue Color
    // Wgr_sys_TOF->SetLineWidth(2);
    // Wgr_sys_TOF->Draw("2 same");

    WgrA->SetMarkerStyle(20);
    WgrA->SetMarkerColor(2);
    WgrA->SetTitle("");
    WgrA->Draw("P same");
    // WgrA->SetTitle("Run17");
    gr_Run11->SetMarkerStyle(27);
    // gr_Run11->SetMarkerColor(1);
    // gr_Run11->SetLineColor(4);
    gr_Run11->Draw("P same");
    // gr_Run11->SetTitle("Run11");
    //  Wgr_sys->SetTitle("");
    //  Wgr_sys_TPC->SetTitle("");
    avg->Update();

    for (int i = 0; i < 10; i++)
    {
        if (i == 0)
            continue;
        if (i == 9)
            continue;
        // bound_line->DrawLine(eta_range[i], 0.03 - .001, eta_range[i], 0.03);
    }

    TLatex *tex = new TLatex();
    tex->SetTextSize(0.06);
    // tex->DrawLatex(-0.25, .0302, "#scale[0.5]{#eta^{#pi^{+}#pi^{-}}bin boundary}");
    //  tex->DrawLatex(-0.50, 0.062, "#scale[0.5]{#eta^{#pi^{+}#pi^{-}}bin boundary}");
    //   tex->DrawLatex();

    TLine *avgL = new TLine(gPad->GetUxmin(), 0., gPad->GetUxmax(), 0.);
    avgL->SetLineStyle(2);
    avgL->Draw();

    TLegend *leg2 = new TLegend(0.18, 0.60, 0.48, 0.88);
    leg2->SetTextSize(0.01);
    leg2->AddEntry("", "#font[22]{#color[2]{STAR Preliminary 2017}}", "");
    leg2->SetTextSize(0.08);
    leg2->AddEntry(JAM_th510, "#font[22]{C. Cocuzza et. al. (JAMDiFF), #sqrt{s} = 510 GeV}", "f");
    leg2->AddEntry("", "#font[22]{p^{#uparrow}+p #rightarrow #pi^{+}#pi^{-}+ X at #sqrt{s} = 510 GeV, L_{int} = 350 pb^{-1}}", "");
    leg2->SetTextSize(0.06);
    // leg2->AddEntry("", "#font[22]{Cone < 0.7}", "");
    if (std::strcmp(mode, "highest_pT") == 0)
    {
        leg2->AddEntry("", "#font[22]{9 < p^{#pi^{+}#pi^{-}}_{T}(GeV/c)< 25, 0.2 < M^{#pi^{+}#pi^{-}}_{inv}(GeV/c^{2})< 4}", "");
    }
    if (std::strcmp(mode, "integrated_pT") == 0)
    {
        leg2->AddEntry("", "#font[22]{2.6 < p^{#pi^{+}#pi^{-}}_{T}(GeV/c) < 25, 0.2 < M^{#pi^{+}#pi^{-}}_{inv}(GeV/c^{2}) < 4}", "");
    }

    leg2->AddEntry("", Form("#font[22]{<p^{#pi^{+}#pi^{-}}_{T}> = %g GeV/c, <M^{#pi^{+}#pi^{-}}_{inv}> = %g GeV/c^{2}}", Avg_pT, Avg_Minv), "");
    leg2->SetTextSize(0.03);
    leg2->AddEntry(Wgr_TOTSYS, "#font[22]{Tot Sys.}", "f");

    //  leg2->AddEntry(WgrA, "Run17");
    //  leg2->AddEntry(Wgr_sys_TPC, "TPC");
    //  leg2->AddEntry(Wgr_sys_TPCorTOF, "TPCorTOF");
    //  leg2->AddEntry(Wgr_sys_TOF, "TOF");
    leg2->AddEntry(gr_Run11, "Run11");

    leg2->Draw();
    avg->Update();

    TLatex text;
    text.SetTextSize(0.025);
    text.SetTextFont(22);
    text.DrawLatex(-0.50, -0.0025, "#pm 1.4% scale uncertanity from beam polarization (not shown)");

    if (std::strcmp(mode, "highest_pT") == 0)
    {
        avg->SaveAs("./TriggerBias/AsymVsEta_ConeLt7_9bin_Average_highestpT_TOF.pdf");
    }
    if (std::strcmp(mode, "integrated_pT") == 0)
    {
        avg->SaveAs("./TriggerBias/AsymVsEta_ConeLt7_9bin_integratedpT_TPCorTOF.pdf");
    }

    TCanvas *c_chi = new TCanvas("c_chi", "canv_chi", 900, 700);
    c_chi->Divide(1, 2);
    c_chi->cd(1);
    h_chi_B->Draw();
    c_chi->cd(2);
    h_chi_Y->Draw();

    if (std::strcmp(mode, "highest_pT") == 0)
    {
        c_chi->Print("./TriggerBias/Chi_NDF_Eta_highestpT_TOF.pdf");
    }
    if (std::strcmp(mode, "integrated_pT") == 0)
    {
        c_chi->Print("./TriggerBias/Chi_NDF_Eta_integratedpT_TPCorTOF.pdf");
    }
}
