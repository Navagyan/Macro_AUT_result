// This code is modefied to calculate asymmetry for individual beams and averaged for total asymmetry
#include <iostream>
#include <fstream>
#include <cstdio>
#include <iomanip>
//#include "std_lib_facilities.h"
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH1.h"
#include "TMath.h"
#include "TGraphErrors.h"
using namespace std;
ofstream Output;

// void asymmetryVsinvMass_averaged(const char *ifile = "Ntuple_V1_P123.root")
void asymmetryVsinvMass_final(const char *ifile = "../Ntuple_RawTree_polCut_TPCorTOF_P123.root")
{
	// Histograms to get mean polarization for yellow and blue beam
	TH1D *hpolB = new TH1D("hpolB", "", 100, 0, 1);
	TH1D *hpolY = new TH1D("hpolY", "", 100, 0, 1);

	TH1D *h_chi_BG = new TH1D("h_chi_BG", "h_chi_BG", 5, 0, 5);
	TH1D *h_chi_BL = new TH1D("h_chi_BL", "h_chi_BL", 5, 0, 5);
	TH1D *h_chi_YG = new TH1D("h_chi_YG", "h_chi_YG", 5, 0, 5);
	TH1D *h_chi_YL = new TH1D("h_chi_YL", "h_chi_YL", 5, 0, 5);
	// pT distribution for each pT bin for average pT
	TH1D *hpT_pair[5];
	TH1D *hpT_pairGtB[5];
	TH1D *hpT_pairLtB[5];
	TH1D *hpT_pairGtY[5];
	TH1D *hpT_pairLtY[5];

	for (Int_t n = 0; n < 5; n++)
	{
		hpT_pair[n] = new TH1D(Form("pT_pair_pbin_%i", n), "", 100, 0, 25);
		hpT_pairGtB[n] = new TH1D(Form("pT_pairGtB_pbin_%i", n), "", 100, 0, 25);
		hpT_pairLtB[n] = new TH1D(Form("pT_pairLtB_pbin_%i", n), "", 100, 0, 25);
		hpT_pairGtY[n] = new TH1D(Form("pT_pairGtY_pbin_%i", n), "", 100, 0, 25);
		hpT_pairLtY[n] = new TH1D(Form("pT_pairLtY_pbin_%i", n), "", 100, 0, 25);
	}

	TFile *f = new TFile(ifile);
	TTree *ntuple1 = (TTree *)f->Get("ntuple1_TPCorTOF"); // get trees
	TTree *ntuple2 = (TTree *)f->Get("ntuple2_TPCorTOF");
	TTree *ntuple4 = (TTree *)f->Get("ntuple4_TPCorTOF");
	TTree *ntuple5 = (TTree *)f->Get("ntuple5_TPCorTOF");
	// TTree *ntuple6 = (TTree *)f->Get("ntuple6_TPCorTOF");

	Output.open("AUT_Vs_Minv_9Bin_trigBias_polCut_TPCorTOF.txt");

	// define variables to hold trees variables

	float eta_pair;
	float PhiRS;
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
	double Minv1[9] = {0}; // for invariant mass bin average
	double Minv2[9] = {0};
	double Minv3[9] = {0};
	double Minv4[9] = {0};
	double Minv5[9] = {0};

	double avg_pT_pair[5] = {0};
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
		Minv1[k] = 0;
		Minv2[k] = 0;
		Minv3[k] = 0;
		Minv4[k] = 0;
		Minv5[k] = 0;
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
	//	ntuple6->SetBranchAddress("polB_corr", &polB_corr);

	// To store average polarization values from histograms
	double avgPolB, avgPolY, rmsB, rmsY, avgPolT, rmsT;

	Int_t nentries = (Int_t)ntuple1->GetEntries(); // entries of trees
	// add friend to the tree. no root file to add friend means the tree is on the same root file
	ntuple1->AddFriend("ntuple2_TPCorTOF");
	ntuple1->AddFriend("ntuple4_TPCorTOF");
	ntuple1->AddFriend("ntuple5_TPCorTOF");
	// double pT[10]={2.80, 3.47,  3.79, 4.12, 4.49, 4.938,  5.505, 6.300, 7.66, 15 };
	// double eta_range[10]={-1.200,   -0.668,  -0.469, -0.281, -0.096, 0.089, 0.275, 0.47,   0.675,  1.2 };
	// double pT[6]={2.80, 3.727, 4.343, 5.157, 6.529, 15.00};
	// double M[10]={0.250, 0.403, 0.516, 0.612, 0.711, 0.803, 0.921, 1.070, 1.286, 4.000};
	// double pT[6]={2.80,   3.71,  4.303 ,   5.084,   6.404,  15.00};//Babu's binning
	double pT[6] = {2.60, 4.628, 5.643, 6.979, 9.265, 25}; // Navagyan's Binnig for 90% of dataset
	// Invariant mass binning for each pT-Bin
	double M1[10] = {0.20, 0.3740, 0.4554, 0.5327, 0.6207, 0.7190, 0.8235, 0.9685, 1.175, 4.000};
	double M2[10] = {0.20, 0.3916, 0.4799, 0.5640, 0.6597, 0.7591, 0.8793, 1.0550, 1.327, 4.000};
	double M3[10] = {0.20, 0.4049, 0.4967, 0.5874, 0.6877, 0.7857, 0.9171, 1.1138, 1.4390, 4.000};
	double M4[10] = {0.200, 0.4249, 0.5217, 0.6222, 0.7244, 0.8255, 0.9678, 1.1913, 1.578, 4.000};
	double M5[10] = {0.2000, 0.4688, 0.5805, 0.6900, 0.7870, 0.9041, 1.067, 1.313, 1.768, 4.00};

	/* Babu Bining for his data set
	   double M1[10]={0.20, 0.3891, 0.4792, 0.5852, 0.6742, 0.7623, 0.8566, 0.9689,  1.098, 4.000};
	   double M2[10]={0.20, 0.3900, 0.4786, 0.5858, 0.6792, 0.7719, 0.8799, 1.019, 1.207, 4.000};
	   double M3[10]= {0.20, 0.3973, 0.4900, 0.6004, 0.6985, 0.7921, 0.912, 1.075, 1.304, 4.000};
	   double M4[10]= {0.200, 0.4096, 0.524, 0.6246, 0.7246, 0.8218, 0.9521, 1.1418, 1.418, 4.000};
	   double M5[10]= {0.20, 0.4397, 0.5648, 0.6778, 0.7744, 0.8823, 1.031, 1.2582, 1.6275, 4.00};
	   */

	// int ne = (int)ntuple6->GetEntries();
	//  for (int pol=0; pol<ne;pol++)
	//  for (int pol=0; pol<10000;pol++)
	//{
	//	ntuple6->GetEntry(pol);
	//	hpolB->Fill(polB_corr);
	//	hpolY->Fill(polY_corr);
	//  }

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

	for (Int_t pbin = 0; pbin < 5; pbin++)
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

		double M[10] = {0};
		// choose correspondig invariant mass bin for different  pT-Bin
		if (pbin == 0)
			for (int k = 0; k < 10; k++)
			{
				M[k] = M1[k];
			}
		if (pbin == 1)
			for (int k = 0; k < 10; k++)
			{
				M[k] = M2[k];
			}
		if (pbin == 2)
			for (int k = 0; k < 10; k++)
			{
				M[k] = M3[k];
			}
		if (pbin == 3)
			for (int k = 0; k < 10; k++)
			{
				M[k] = M4[k];
			}
		if (pbin == 4)
			for (int k = 0; k < 10; k++)
			{
				M[k] = M5[k];
			}

		cout << "Pbin: " << pbin << ", M[10]={" << M[0] << ", " << M[1] << ", " << M[2] << ", " << M[3] << ", " << M[4] << ", " << M[5] << ", " << M[6] << ", " << M[7] << ", " << M[8] << ", " << M[9] << "}" << endl;
		;
		cout << "Entries: " << nentries << endl;
		for (Int_t i = 0; i < nentries; i++)
		// for (Int_t i = 0; i < 100000; i++)
		{

			ntuple1->GetEntry(i);
			if (pT_pair < pT[pbin] || pT_pair >= pT[pbin + 1])
				continue;
			// if (pT_pair < 3.75)
			//	continue;
			if (cone >= .7)
				continue;
			if (Minv > 4.)
				continue;
			// if(!(Minv>0.5076 || Minv<0.4876))  continue;//dodge K0 mass range, doesn't cause asymmetry. Implemention was not working before.This is fine now.
			if (fitPts_min_pair < 15)
				continue;

			hpT_pair[pbin]->Fill(pT_pair);

			//-------------BLUE-------------------------------------
			// Phi
			for (int phi = 0; phi < 16; phi++)
			{
				if (PhiRSB >= (phi - 8.) / 8. * pi && PhiRSB <= (phi - 7.) / 8. * pi)
				{
					// Invariant mass loop
					for (int m = 0; m < 9; m++)
					{
						if (Minv >= M[m] && Minv < M[m + 1])
						{
							Npairs[m] = Npairs[m] + 1;
							pTpairs[m] = pTpairs[m] + pT_pair;
							Mpairs[m] = Mpairs[m] + Minv;

							if (eta_pair > 0)
							{
								hpT_pairGtB[pbin]->Fill(Minv);
								NpairsBGt[pbin][m] = NpairsBGt[pbin][m] + 1;
								pTpairsBGt[pbin][m] = pTpairsBGt[pbin][m] + pT_pair;
								MpairsBGt[pbin][m] = MpairsBGt[pbin][m] + Minv;
								EtapairsBGt[pbin][m] = EtapairsBGt[pbin][m] + eta_pair;
							}
							if (eta_pair < 0)
							{
								hpT_pairLtB[pbin]->Fill(Minv);
								NpairsBLt[pbin][m] = NpairsBLt[pbin][m] + 1;
								pTpairsBLt[pbin][m] = pTpairsBLt[pbin][m] + pT_pair;
								MpairsBLt[pbin][m] = MpairsBLt[pbin][m] + Minv;
								EtapairsBLt[pbin][m] = EtapairsBLt[pbin][m] + eta_pair;
							}

							if (fspinconfig == 51 || fspinconfig == 53)
							{
								if (eta_pair > 0)
								{
									NMinvGtUpB[m * 16 + phi]++;
								}
								if (eta_pair < 0)
								{
									NMinvLtUpB[m * 16 + phi]++;
								}
							}
							if (fspinconfig == 83 || fspinconfig == 85)
							{
								if (eta_pair > 0)
								{
									NMinvGtDnB[m * 16 + phi]++;
								}
								if (eta_pair < 0)
								{
									NMinvLtDnB[m * 16 + phi]++;
								}
							}
						} // Minv loop bin control
					}	  // Minv loop

				} // phi-range loop

			}																// phi loop
			/**************** BLUE BEAM ENDS *****************************/ ///////

			////***************** YELLOW BEAM *****************************////////////
			// Phi loop
			for (int phi = 0; phi < 16; phi++)
			{
				if (PhiRSY >= (phi - 8.) / 8. * pi && PhiRSY <= (phi - 7.) / 8. * pi)
				{
					// cout << " Yellow Phi loop" << phi << endl;
					//  Minv loop
					for (int m = 0; m < 9; m++)
					{
						if (Minv >= M[m] && Minv < M[m + 1])
						{
							// cout << " Yellow Minv loop" << m << endl;
							// cout << " Yellow Minv " << Minv << endl;
							if (eta_pair > 0)
							{
								hpT_pairGtY[pbin]->Fill(Minv);
								NpairsYGt[pbin][m] = NpairsYGt[pbin][m] + 1;
								pTpairsYGt[pbin][m] = pTpairsYGt[pbin][m] + pT_pair;
								MpairsYGt[pbin][m] = MpairsYGt[pbin][m] + Minv;
								EtapairsYGt[pbin][m] = EtapairsYGt[pbin][m] + eta_pair;
							}
							if (eta_pair < 0)
							{
								hpT_pairLtY[pbin]->Fill(Minv);
								NpairsYLt[pbin][m] = NpairsYLt[pbin][m] + 1;
								pTpairsYLt[pbin][m] = pTpairsYLt[pbin][m] + pT_pair;
								MpairsYLt[pbin][m] = MpairsYLt[pbin][m] + Minv;
								EtapairsYLt[pbin][m] = EtapairsYLt[pbin][m] + eta_pair;
							}
							if (fspinconfig == 51 || fspinconfig == 83)
							{
								// if(eta_pair>0)
								if (eta_pair < 0)
								{
									NMinvGtUpY[m * 16 + phi]++;
								}
								// if(eta_pair<0)
								if (eta_pair > 0)
								{
									NMinvLtUpY[m * 16 + phi]++;
								}
							}
							if (fspinconfig == 53 || fspinconfig == 85)
							{
								// if(eta_pair>0)
								if (eta_pair < 0)
								{
									NMinvGtDnY[m * 16 + phi]++;
								}
								// if(eta_pair<0)
								if (eta_pair > 0)
								{
									NMinvLtDnY[m * 16 + phi]++;
								}
							}
						}
					} // Minv loop
				}	  // PhiRS range loop
			}		  // phi loop
		}			  // entries  for loop

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

		// Store Minv bin averages for each pT bin in an array for plotting purpose
		// Here you are taking average Minv for Blue and Yellow and Average all the same and plot them......might need to reconsider this one...

		if (pbin == 0)
		{
			for (int k0 = 0; k0 < 9; k0++)
			{
				Minv1[k0] = Mpairs[k0] / (double)Npairs[k0];
			}
		}
		if (pbin == 1)
		{
			for (int k1 = 0; k1 < 9; k1++)
			{
				Minv2[k1] = Mpairs[k1] / (double)Npairs[k1];
			}
		}
		if (pbin == 2)
		{
			for (int k2 = 0; k2 < 9; k2++)
			{
				Minv3[k2] = Mpairs[k2] / (double)Npairs[k2];
			}
		}
		if (pbin == 3)
		{
			for (int k3 = 0; k3 < 9; k3++)
			{
				Minv4[k3] = Mpairs[k3] / (double)Npairs[k3];
			}
		}
		if (pbin == 4)
		{
			for (int k4 = 0; k4 < 9; k4++)
			{
				Minv5[k4] = Mpairs[k4] / (double)Npairs[k4];
			}
		}
		gStyle->SetOptDate(0);
		gStyle->SetOptFit(1);
		// Asymmetry calculation for BLUE eta >0
		double dAdB, dAdC, dAdD, dAdE, dAdP; // variables for error calculation
		// Use rms for polarization error from distribution!
		double dP_B = rmsB; // Polarization errors blue
		double dP_Y = rmsY; // Polarization errors yellow
		double dA[165];		// Asymmetry amplitude from sine fit
		double B, C, D, E;
		double a, b;
		double Asym[165];
		// double pi = 3.14159265359;

		double Abg[9] = {0}; // Here Abg refers to Asymmetry for Blue in eta >0
		double deltaAbg[9] = {0};
		double Ayg[9] = {0};
		double deltaAyg[9] = {0};

		for (int m = 0; m < 9; m++)
		{
			for (int ang = 0; ang < 16; ang++)
			{
				if (ang < 8)
				{
					a = sqrt(double(NMinvGtUpB[m * 16 + ang]) * double(NMinvGtDnB[m * 16 + ang + 8]));
					b = sqrt(double(NMinvGtDnB[m * 16 + ang]) * double(NMinvGtUpB[m * 16 + ang + 8]));
					B = double(NMinvGtUpB[m * 16 + ang]);
					C = double(NMinvGtUpB[m * 16 + ang + 8]);
					D = double(NMinvGtDnB[m * 16 + ang]);
					E = double(NMinvGtDnB[m * 16 + ang + 8]);
				}
				if (ang > 7)
				{
					a = sqrt(double(NMinvGtUpB[m * 16 + ang]) * double(NMinvGtDnB[m * 16 + ang - 8]));
					b = sqrt(double(NMinvGtDnB[m * 16 + ang]) * double(NMinvGtUpB[m * 16 + ang - 8]));
					B = double(NMinvGtUpB[m * 16 + ang]);
					C = double(NMinvGtUpB[m * 16 + ang - 8]);
					D = double(NMinvGtDnB[m * 16 + ang]);
					E = double(NMinvGtDnB[m * 16 + ang - 8]);
				}
				Asym[ang] = (1. / avgPolB) * ((a - b) / (a + b)); // A_UT for each phiRS bin for each Minv Bin
				// cout << "asymmetry for blue beam for bin "<< ang <<" --> "<< Asym[ang] << endl;
				dAdB = (1. / avgPolB) * E * sqrt(D * C) / (sqrt(B * E) * (pow((sqrt(B * E) + sqrt(D * C)), 2)));
				dAdE = (1. / avgPolB) * B * sqrt(D * C) / (sqrt(B * E) * (pow((sqrt(B * E) + sqrt(D * C)), 2)));
				dAdD = (-1. / avgPolB) * C * sqrt(B * E) / (sqrt(D * C) * (pow((sqrt(B * E) + sqrt(D * C)), 2)));
				dAdC = (-1. / avgPolB) * D * sqrt(B * E) / (sqrt(D * C) * (pow((sqrt(B * E) + sqrt(D * C)), 2)));
				dAdP = (-1. / (avgPolB * avgPolB)) * (sqrt(B * E) - sqrt(D * C)) / (sqrt(B * E) + sqrt(D * C));
				dA[ang] = sqrt(pow((fabs(dAdB) * sqrt(B)), 2) + pow((fabs(dAdC) * sqrt(C)), 2) + pow((fabs(dAdD) * sqrt(D)), 2) + pow((fabs(dAdE) * sqrt(E)), 2) + pow((fabs(dAdP) * dP_B), 2));
				// cout <<"polB:"<<avgPolB<<", polY: "<<avgPolY<< ", a="<<a<<", b="<<b<<", B="<<B<<", C="<<C<<", D="<<D<<", E="<<E<<", Asym["<<ang<<"]="<<Asym[ang]<<", dA["<<ang<<"]="<<dA[ang]<<endl;
			} // angle loop
			gStyle->SetOptFit(1);
			FitCanv_Blue_Gt->cd(m + 1);
			char name[600];
			char title[600];
			double chi2Ndf[9];
			double angle[16] = {-15. / 16. * pi, -13. / 16. * pi, -11. / 16. * pi, -9. / 16. * pi, -7. / 16. * pi, -5. / 16. * pi, -3. / 16. * pi, -1. / 16. * pi, 1. / 16. * pi, 3. / 16. * pi, 5. / 16. * pi, 7. / 16. * pi, 9. / 16. * pi, 11. / 16. * pi, 13. / 16. * pi, 15. / 16. * pi};
			double ex[16] = {0};

			auto grt = new TGraphErrors(16, angle, Asym, ex, dA); // Vital Step this grt Pointer will be used to plot sin Fit and then amplitude of that sin fit will be used to determine the asymetry for that particular Minv Bin and store in Abg[m]
			grt->SetMarkerStyle(20);
			grt->Draw("AP");
			TF1 *fit = new TF1("fit", "[0]*sin(x)", -3.14159265359, 0);
			// fit->SetParameter(0,0.0001);
			grt->Fit(fit, "R");
			grt->SetMarkerColor(4);
			grt->GetXaxis()->SetTitle("#Phi_{RS}");
			grt->GetXaxis()->SetTitleOffset(1);
			// sprintf(title, "pT bin %i, Minv bin %i, BLUE, #eta > 0", pbin, m);
			snprintf(title, sizeof(title), "pT bin %i, Minv bin %i, BLUE, #eta > 0", pbin, m);
			grt->SetTitle(title);
			chi2Ndf[m] = fit->GetChisquare() / fit->GetNDF();
			cout << "ChiSquare = " << chi2Ndf[m] << endl;
			grt->GetYaxis()->SetTitle("A_{UT}(#Phi_{RS})");
			grt->GetYaxis()->CenterTitle();
			grt->GetYaxis()->SetRangeUser(-0.08, 0.08);
			Abg[m] = fit->GetParameter(0); // Extractring the asymmetry
			deltaAbg[m] = fit->GetParError(0);
			h_chi_BG->Fill(chi2Ndf[m]);
			FitCanv_Blue_Gt->Update();
			// sprintf(name,"pT%i_Minvbin%i_#etaGt_BLUE_9bin.png",pbin, m);
			// c1->SaveAs(name);

		} // invariant mass loop
		if (pbin == 0)
		{
			FitCanv_Blue_Gt->Print("./TriggerBias/FitPlot_pBin_Blue_Gt.pdf(", "pdf");
		}
		else if (pbin == 4)
		{
			FitCanv_Blue_Gt->Print("./TriggerBias/FitPlot_pBin_Blue_Gt.pdf)", "pdf");
		}
		else
		{
			FitCanv_Blue_Gt->Print("./TriggerBias/FitPlot_pBin_Blue_Gt.pdf", "pdf");
		}
		// YELLOW eta > 0

		for (int m = 0; m < 9; m++)
		{
			for (int ang = 0; ang < 16; ang++)
			{
				if (ang < 8)
				{
					a = sqrt(double(NMinvGtUpY[m * 16 + ang]) * double(NMinvGtDnY[m * 16 + ang + 8]));
					b = sqrt(double(NMinvGtDnY[m * 16 + ang]) * double(NMinvGtUpY[m * 16 + ang + 8]));
					B = double(NMinvGtUpY[m * 16 + ang]);
					C = double(NMinvGtUpY[m * 16 + ang + 8]);
					D = double(NMinvGtDnY[m * 16 + ang]);
					E = double(NMinvGtDnY[m * 16 + ang + 8]);
				}
				if (ang > 7)
				{
					a = sqrt(double(NMinvGtUpY[m * 16 + ang]) * double(NMinvGtDnY[m * 16 + ang - 8]));
					b = sqrt(double(NMinvGtDnY[m * 16 + ang]) * double(NMinvGtUpY[m * 16 + ang - 8]));
					B = double(NMinvGtUpY[m * 16 + ang]);
					C = double(NMinvGtUpY[m * 16 + ang - 8]);
					D = double(NMinvGtDnY[m * 16 + ang]);
					E = double(NMinvGtDnY[m * 16 + ang - 8]);
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
			auto gr = new TGraphErrors(16, angle, Asym, ex, dA);
			gr->SetMarkerStyle(20);
			gr->Draw("AP");
			TF1 *fit = new TF1("fit", "[0]*sin(x)", -3.14159265359, 0);
			// fit->SetParameter(0,0.0001);
			gr->Fit(fit, "R");
			gr->SetMarkerColor(4);
			gr->GetXaxis()->SetTitle("#Phi_{RS}");
			gr->GetXaxis()->SetTitleOffset(1);
			// sprintf(title, "pT bin %i, Minv bin %i, YELLOW, #eta > 0", pbin, m);
			snprintf(title, sizeof(title), "pT bin %i, Minv bin %i, YELLOW, #eta > 0", pbin, m);
			gr->SetTitle(title);
			chi2Ndf[m] = fit->GetChisquare() / fit->GetNDF();
			cout << "ChiSquare= " << chi2Ndf[m] << endl;
			gr->GetYaxis()->SetTitle("A_{UT}(#Phi_{RS})");
			gr->GetYaxis()->CenterTitle(kTRUE);
			gr->GetYaxis()->SetRangeUser(-0.08, 0.08);
			Ayg[m] = fit->GetParameter(0);
			deltaAyg[m] = fit->GetParError(0);
			h_chi_YG->Fill(chi2Ndf[m]);
			FitCanv_Yellow_Gt->Update();

			//	sprintf(name,"pT%i_Minvbin%i_#etaGt_YELLOW_9bin.png",pbin, m);
			//	c1->SaveAs(name);
		} // Minv loop
		if (pbin == 0)
		{
			FitCanv_Yellow_Gt->Print("./TriggerBias/FitPlot_pBin_Yellow_Gt.pdf(", "pdf");
		}
		else if (pbin == 4)
		{
			FitCanv_Yellow_Gt->Print("./TriggerBias/FitPlot_pBin_Yellow_Gt.pdf)", "pdf");
		}
		else
		{
			FitCanv_Yellow_Gt->Print("./TriggerBias/FitPlot_pBin_Yellow_Gt.pdf", "pdf");
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
					a = sqrt(double(NMinvLtUpB[m * 16 + ang]) * double(NMinvLtDnB[m * 16 + ang + 8]));
					b = sqrt(double(NMinvLtDnB[m * 16 + ang]) * double(NMinvLtUpB[m * 16 + ang + 8]));
					B = double(NMinvLtUpB[m * 16 + ang]);
					C = double(NMinvLtUpB[m * 16 + ang + 8]);
					D = double(NMinvLtDnB[m * 16 + ang]);
					E = double(NMinvLtDnB[m * 16 + ang + 8]);
				}
				if (ang > 7)
				{
					a = sqrt(double(NMinvLtUpB[m * 16 + ang]) * double(NMinvLtDnB[m * 16 + ang - 8]));
					b = sqrt(double(NMinvLtDnB[m * 16 + ang]) * double(NMinvLtUpB[m * 16 + ang - 8]));
					B = double(NMinvLtUpB[m * 16 + ang]);
					C = double(NMinvLtUpB[m * 16 + ang - 8]);
					D = double(NMinvLtDnB[m * 16 + ang]);
					E = double(NMinvLtDnB[m * 16 + ang - 8]);
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
			auto grt = new TGraphErrors(16, angle, Asym, ex, dA);
			grt->SetMarkerStyle(20);
			grt->Draw("AP");
			TF1 *fit = new TF1("fit", "[0]*sin(x)", -3.14159265359, 0);
			// fit->SetParameter(0,0.0001);
			grt->Fit(fit, "R");
			grt->SetMarkerColor(4);
			grt->GetXaxis()->SetTitle("#Phi_{RS}");
			grt->GetXaxis()->SetTitleOffset(1);
			// sprintf(title, "pT bin %i, Minv bin %i, BLUE, #eta < 0", pbin, m);
			snprintf(title, sizeof(title), "pT bin %i, Minv bin %i, BLUE, #eta < 0", pbin, m);
			grt->SetTitle(title);
			chi2Ndf[m] = fit->GetChisquare() / fit->GetNDF();
			cout << "ChiSquare = " << chi2Ndf[m] << endl;
			grt->GetYaxis()->SetTitle("A_{UT}(#Phi_{RS})");
			grt->GetYaxis()->SetTitleOffset(1);
			grt->GetYaxis()->CenterTitle(kTRUE);
			grt->GetYaxis()->SetRangeUser(-0.08, 0.08);
			Abl[m] = fit->GetParameter(0);
			deltaAbl[m] = fit->GetParError(0);
			h_chi_BL->Fill(chi2Ndf[m]);
			FitCanv_Blue_Lt->Update();
			//	sprintf(name,"pT%i_Minvbin%i_#etaLt_BLUE_9bin.png",pbin, m);
			//	c1->Print(name);

		} // invariant mass loop
		if (pbin == 0)
		{
			FitCanv_Blue_Lt->Print("./TriggerBias/FitPlot_pBin_Blue_Lt.pdf(", "pdf");
		}
		else if (pbin == 4)
		{
			FitCanv_Blue_Lt->Print("./TriggerBias/FitPlot_pBin_Blue_Lt.pdf)", "pdf");
		}
		else
		{
			FitCanv_Blue_Lt->Print("./TriggerBias/FitPlot_pBin_Blue_Lt.pdf", "pdf");
		}
		// YELLOW Beam Eta<0

		for (int m = 0; m < 9; m++)
		{
			for (int ang = 0; ang < 16; ang++)
			{
				if (ang < 8)
				{
					a = sqrt(double(NMinvLtUpY[m * 16 + ang]) * double(NMinvLtDnY[m * 16 + ang + 8]));
					b = sqrt(double(NMinvLtDnY[m * 16 + ang]) * double(NMinvLtUpY[m * 16 + ang + 8]));
					B = double(NMinvLtUpY[m * 16 + ang]);
					C = double(NMinvLtUpY[m * 16 + ang + 8]);
					D = double(NMinvLtDnY[m * 16 + ang]);
					E = double(NMinvLtDnY[m * 16 + ang + 8]);
				}
				if (ang > 7)
				{
					a = sqrt(double(NMinvLtUpY[m * 16 + ang]) * double(NMinvLtDnY[m * 16 + ang - 8]));
					b = sqrt(double(NMinvLtDnY[m * 16 + ang]) * double(NMinvLtUpY[m * 16 + ang - 8]));
					B = double(NMinvLtUpY[m * 16 + ang]);
					C = double(NMinvLtUpY[m * 16 + ang - 8]);
					D = double(NMinvLtDnY[m * 16 + ang]);
					E = double(NMinvLtDnY[m * 16 + ang - 8]);
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
			auto gr = new TGraphErrors(16, angle, Asym, ex, dA);
			gr->SetMarkerStyle(20);
			gr->Draw("AP");
			TF1 *fit = new TF1("fit", "[0]*sin(x)", -3.14159265359, 0);
			fit->SetParameter(0, 0.0001);
			gr->Fit(fit, "R");
			gr->SetMarkerColor(4);
			gr->GetXaxis()->SetTitle("#Phi_{RS}");
			gr->GetXaxis()->SetTitleOffset(1);
			// sprintf(title, "pT bin %i, Minv bin %i, YELLOW, #eta < 0", pbin, m);
			snprintf(title, sizeof(title), "pT bin %i, Minv bin %i, YELLOW, #eta < 0", pbin, m);
			chi2Ndf[m] = fit->GetChisquare() / fit->GetNDF();
			cout << "ChiSquare= " << chi2Ndf[m] << endl;
			gr->SetTitle(title);
			gr->GetYaxis()->SetTitle("A_{UT}(#Phi_{RS})");
			gr->GetYaxis()->CenterTitle(kTRUE);
			gr->GetYaxis()->SetRangeUser(-0.08, 0.08);
			Ayl[m] = fit->GetParameter(0);
			deltaAyl[m] = fit->GetParError(0);
			h_chi_YL->Fill(chi2Ndf[m]);
			FitCanv_Yellow_Lt->Update();
			//	sprintf(name,"pTbin%i_Minvbin%i_etaLt0_YELLOW_9bin.png",pbin, m);
			//	c1->Print(name);

		} // Minv loop
		if (pbin == 0)
		{
			FitCanv_Yellow_Lt->Print("./TriggerBias/FitPlot_pBin_Yellow_Lt.pdf(", "pdf");
		}
		else if (pbin == 4)
		{
			FitCanv_Yellow_Lt->Print("./TriggerBias/FitPlot_pBin_Yellow_Lt.pdf)", "pdf");
		}
		else
		{
			FitCanv_Yellow_Lt->Print("./TriggerBias/FitPlot_pBin_Yellow_Lt.pdf", "pdf");
		}
		// Yellow beam asymmetry ends

		if (pbin == 0)
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

			// double A_pT1yl[9] = {Ayl[0], Ayl[1], Ayl[2], Ayl[3], Ayl[4], Ayl[5], Ayl[6], Ayl[7], Ayl[8]};
			// double deltaA_pT1yl[9] = {deltaAyl[0], deltaAyl[1], deltaAyl[2], deltaAyl[3], deltaAyl[4], deltaAyl[5], deltaAyl[6], deltaAyl[7], deltaAyl[8]};
		}

		if (pbin == 1)
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
		if (pbin == 2)
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
		if (pbin == 3)
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
		if (pbin == 4)
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

	} // pbin bin loop

	double MinvBGt[5][9] = {0};
	double MinvBLt[5][9] = {0};
	double MinvYGt[5][9] = {0};
	double MinvYLt[5][9] = {0};
	double MinvAvgGt[5][9] = {0};
	double MinvAvgLt[5][9] = {0};

	double pTBGt[5][9] = {0};
	double pTBLt[5][9] = {0};
	double pTYGt[5][9] = {0};
	double pTYLt[5][9] = {0};
	double pTAvgGt[5][9] = {0};
	double pTAvgLt[5][9] = {0};

	double Eta_pairBGt[5][9] = {0};
	double Eta_pairBLt[5][9] = {0};
	double Eta_pairYGt[5][9] = {0};
	double Eta_pairYLt[5][9] = {0};
	double Eta_pairAvgGt[5][9] = {0};
	double Eta_pairAvgLt[5][9] = {0};

	double avgpT_pairGt[5] = {0};
	double avgpT_pairLt[5] = {0};
	double avg_pT_pairGtB[5] = {0};
	double avg_pT_pairLtB[5] = {0};
	double avg_pT_pairGtY[5] = {0};
	double avg_pT_pairLtY[5] = {0};

	for (Int_t n = 0; n < 5; n++)
	{
		avg_pT_pair[n] = hpT_pair[n]->GetMean();
		// cout << "avg_pT_pair[" << n << "]=" << avg_pT_pair[n] << endl;
		avg_pT_pairGtB[n] = hpT_pairGtB[n]->GetMean();
		avg_pT_pairLtB[n] = hpT_pairLtB[n]->GetMean();
		avg_pT_pairGtY[n] = hpT_pairGtY[n]->GetMean();
		avg_pT_pairLtY[n] = hpT_pairLtY[n]->GetMean();
	}

	for (int pbin = 0; pbin < 5; pbin++)
	{
		avgpT_pairGt[pbin] = 0.5 * (avg_pT_pairGtB[pbin] + avg_pT_pairGtY[pbin]);
		avgpT_pairLt[pbin] = 0.5 * (avg_pT_pairLtB[pbin] + avg_pT_pairLtY[pbin]);
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

	Output << "avgpT_pairGtB[5]={" << avg_pT_pairGtB[0] << "," << avg_pT_pairGtB[1] << "," << avg_pT_pairGtB[2] << "," << avg_pT_pairGtB[3] << "," << avg_pT_pairGtB[4] << "}" << endl;
	Output << "avgpT_pairGtY[5]={" << avg_pT_pairGtY[0] << "," << avg_pT_pairGtY[1] << "," << avg_pT_pairGtY[2] << "," << avg_pT_pairGtY[3] << "," << avg_pT_pairGtY[4] << "}" << endl;
	Output << "avgpT_pairLtB[5]={" << avg_pT_pairLtB[0] << "," << avg_pT_pairLtB[1] << "," << avg_pT_pairLtB[2] << "," << avg_pT_pairLtB[3] << "," << avg_pT_pairLtB[4] << "}" << endl;
	Output << "avgpT_pairLtY[5]={" << avg_pT_pairLtY[0] << "," << avg_pT_pairLtY[1] << "," << avg_pT_pairLtY[2] << "," << avg_pT_pairLtY[3] << "," << avg_pT_pairLtY[4] << "}" << endl;
	Output << "avgpT_pairGt[5]={" << avgpT_pairGt[0] << "," << avgpT_pairGt[1] << "," << avgpT_pairGt[2] << "," << avgpT_pairGt[3] << "," << avgpT_pairGt[4] << "}" << endl;
	Output << "avgpT_pairLt[5]={" << avgpT_pairLt[0] << "," << avgpT_pairLt[1] << "," << avgpT_pairLt[2] << "," << avgpT_pairLt[3] << "," << avgpT_pairLt[4] << "}" << endl;

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
		// total asymmetry error calculation
		// since two asymmetry(BLUE+YELLOW) are averaged, error propagation formaula is used for total asymmetry error
		// if function, f = (a+b)/2, Then Error, df =1/2 √{(∂f/∂a)^2*(∆a)^2+(∂f/∂b)^2*(∆b)^2}

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
	}

	Output << "                                                                                                  " << endl;
	Output << "                                                                                                  " << endl;
	Output << "                                                                                                  " << endl;
	Output << "                                      BLUE BEAM                                                   " << endl;
	Output << "********************   FINAL ASYMMETRY , CONE > 0.7, Avg Minv  Binning for BLUE BEAM Only****************" << endl;

	Output << "<<<<<<<<<<<<<<<   Minv Bin Averages for asymemtry plot  >>>>>>>>>>>>>>" << endl;
	Output << "pT-1, MinvBGt[0][9]= {" << MinvBGt[0][0] << ", " << MinvBGt[0][1] << ", " << MinvBGt[0][2] << ", " << MinvBGt[0][3] << ", " << MinvBGt[0][4] << ", " << MinvBGt[0][5] << ", " << MinvBGt[0][6] << ", " << MinvBGt[0][7] << ", " << MinvBGt[0][8] << "}" << endl;
	Output << "pT-2, MinvBGt[1][9]= {" << MinvBGt[1][0] << ", " << MinvBGt[1][1] << ", " << MinvBGt[1][2] << ", " << MinvBGt[1][3] << ", " << MinvBGt[1][4] << ", " << MinvBGt[1][5] << ", " << MinvBGt[1][6] << ", " << MinvBGt[1][7] << ", " << MinvBGt[1][8] << "}" << endl;
	Output << "pT-3, MinvBGt[2][9]= {" << MinvBGt[2][0] << ", " << MinvBGt[2][1] << ", " << MinvBGt[2][2] << ", " << MinvBGt[2][3] << ", " << MinvBGt[2][4] << ", " << MinvBGt[2][5] << ", " << MinvBGt[2][6] << ", " << MinvBGt[2][7] << ", " << MinvBGt[2][8] << "}" << endl;
	Output << "pT-4, MinvBGt[3][9]= {" << MinvBGt[3][0] << ", " << MinvBGt[3][1] << ", " << MinvBGt[3][2] << ", " << MinvBGt[3][3] << ", " << MinvBGt[3][4] << ", " << MinvBGt[3][5] << ", " << MinvBGt[3][6] << ", " << MinvBGt[3][7] << ", " << MinvBGt[3][8] << "}" << endl;
	Output << "pT-5, MinvBGt[4][9]= {" << MinvBGt[4][0] << ", " << MinvBGt[4][1] << ", " << MinvBGt[4][2] << ", " << MinvBGt[4][3] << ", " << MinvBGt[4][4] << ", " << MinvBGt[4][5] << ", " << MinvBGt[4][6] << ", " << MinvBGt[4][7] << ", " << MinvBGt[4][8] << "}" << endl;

	Output << "<<<<<<<<<<<<<<<   Minv Bin Averages for asymemtry plot  >>>>>>>>>>>>>>" << endl;
	Output << "pT-1, MinvBLt[0][9]= {" << MinvBLt[0][0] << ", " << MinvBLt[0][1] << ", " << MinvBLt[0][2] << ", " << MinvBLt[0][3] << ", " << MinvBLt[0][4] << ", " << MinvBLt[0][5] << ", " << MinvBLt[0][6] << ", " << MinvBLt[0][7] << ", " << MinvBLt[0][8] << "}" << endl;
	Output << "pT-2, MinvBLt[1][9]= {" << MinvBLt[1][0] << ", " << MinvBLt[1][1] << ", " << MinvBLt[1][2] << ", " << MinvBLt[1][3] << ", " << MinvBLt[1][4] << ", " << MinvBLt[1][5] << ", " << MinvBLt[1][6] << ", " << MinvBLt[1][7] << ", " << MinvBLt[1][8] << "}" << endl;
	Output << "pT-3, MinvBLt[2][9]= {" << MinvBLt[2][0] << ", " << MinvBLt[2][1] << ", " << MinvBLt[2][2] << ", " << MinvBLt[2][3] << ", " << MinvBLt[2][4] << ", " << MinvBLt[2][5] << ", " << MinvBLt[2][6] << ", " << MinvBLt[2][7] << ", " << MinvBLt[2][8] << "}" << endl;
	Output << "pT-4, MinvBLt[3][9]= {" << MinvBLt[3][0] << ", " << MinvBLt[3][1] << ", " << MinvBLt[3][2] << ", " << MinvBLt[3][3] << ", " << MinvBLt[3][4] << ", " << MinvBLt[3][5] << ", " << MinvBLt[3][6] << ", " << MinvBLt[3][7] << ", " << MinvBLt[3][8] << "}" << endl;
	Output << "pT-5, MinvBLt[4][9]= {" << MinvBLt[4][0] << ", " << MinvBLt[4][1] << ", " << MinvBLt[4][2] << ", " << MinvBLt[4][3] << ", " << MinvBLt[4][4] << ", " << MinvBLt[4][5] << ", " << MinvBLt[4][6] << ", " << MinvBLt[4][7] << ", " << MinvBLt[4][8] << "}" << endl;

	Output << "<<<<<<<< Average Polarization Values >>>>>>>>" << endl;
	Output << "BLUE <P>: " << avgPolB << ", Err_P: " << rmsB << endl;
	Output << "YELLOW <P>: " << avgPolY << ", Err_P: " << rmsY << endl;

	Output << "                                                                                                  " << endl;
	Output << "                                                                                                  " << endl;
	Output << "<<<<<<<<<<<       Eta < 0, BLUE BEAM A_{UT}   >>>>>>>>>>>>>>>>>" << endl;
	Output << "pT-1,<A_{UT}> ={" << A_pT1bl[0] << ", " << A_pT1bl[1] << ", " << A_pT1bl[2] << ", " << A_pT1bl[3] << ", " << A_pT1bl[4] << ", " << A_pT1bl[5] << ", " << A_pT1bl[6] << ", " << A_pT1bl[7] << ", " << A_pT1bl[8] << "}" << endl;
	Output << "pT-1,<Err>={" << deltaA_pT1bl[0] << ", " << deltaA_pT1bl[1] << ", " << deltaA_pT1bl[2] << ", " << deltaA_pT1bl[3] << ", " << deltaA_pT1bl[4] << ", " << deltaA_pT1bl[5] << ", " << deltaA_pT1bl[6] << ", " << deltaA_pT1bl[7] << ", " << deltaA_pT1bl[8] << "}" << endl;
	Output << "pT-2,<A_{UT}>={" << A_pT2bl[0] << ", " << A_pT2bl[1] << ", " << A_pT2bl[2] << ", " << A_pT2bl[3] << ", " << A_pT2bl[4] << ", " << A_pT2bl[5] << ", " << A_pT2bl[6] << ", " << A_pT2bl[7] << ", " << A_pT2bl[8] << "}" << endl;
	Output << "pT-2,<Err> ={" << deltaA_pT2bl[0] << ", " << deltaA_pT2bl[1] << ", " << deltaA_pT2bl[2] << ", " << deltaA_pT2bl[3] << ", " << deltaA_pT2bl[4] << ", " << deltaA_pT2bl[5] << ", " << deltaA_pT2bl[6] << ", " << deltaA_pT2bl[7] << ", " << deltaA_pT2bl[8] << "}" << endl;
	Output << "pT-3,<A_{UT}> ={" << A_pT3bl[0] << ", " << A_pT3bl[1] << ", " << A_pT3bl[2] << ", " << A_pT3bl[3] << ", " << A_pT3bl[4] << ", " << A_pT3bl[5] << ", " << A_pT3bl[6] << ", " << A_pT3bl[7] << ", " << A_pT3bl[8] << "}" << endl;
	Output << "pT-3,<Err> ={" << deltaA_pT3bl[0] << ", " << deltaA_pT3bl[1] << ", " << deltaA_pT3bl[2] << ", " << deltaA_pT3bl[3] << ", " << deltaA_pT3bl[4] << ", " << deltaA_pT3bl[5] << ", " << deltaA_pT3bl[6] << ", " << deltaA_pT3bl[7] << ", " << deltaA_pT3bl[8] << "}" << endl;
	Output << "pT-4,<A_{UT}> ={" << A_pT4bl[0] << ", " << A_pT4bl[1] << ", " << A_pT4bl[2] << ", " << A_pT4bl[3] << ", " << A_pT4bl[4] << ", " << A_pT4bl[5] << ", " << A_pT4bl[6] << ", " << A_pT4bl[7] << ", " << A_pT4bl[8] << "}" << endl;
	Output << "pT-4,<Err> ={" << deltaA_pT4bl[0] << ", " << deltaA_pT4bl[1] << ", " << deltaA_pT4bl[2] << ", " << deltaA_pT4bl[3] << ", " << deltaA_pT4bl[4] << ", " << deltaA_pT4bl[5] << ", " << deltaA_pT4bl[6] << ", " << deltaA_pT4bl[7] << ", " << deltaA_pT4bl[8] << "}" << endl;
	Output << "pT-5,<A_{UT}> ={" << A_pT5bl[0] << ", " << A_pT5bl[1] << ", " << A_pT5bl[2] << ", " << A_pT5bl[3] << ", " << A_pT5bl[4] << ", " << A_pT5bl[5] << ", " << A_pT5bl[6] << ", " << A_pT5bl[7] << ", " << A_pT5bl[8] << "}" << endl;
	Output << "pT-5,<Err> ={" << deltaA_pT5bl[0] << ", " << deltaA_pT5bl[1] << ", " << deltaA_pT5bl[2] << ", " << deltaA_pT5bl[3] << ", " << deltaA_pT5bl[4] << ", " << deltaA_pT5bl[5] << ", " << deltaA_pT5bl[6] << ", " << deltaA_pT5bl[7] << ", " << deltaA_pT5bl[8] << "}" << endl;

	Output << "                                                                                                  " << endl;
	Output << "                                                                                                  " << endl;
	Output << "<<<<<<<<<<<       Eta > 0, BLUE BEAM A_{UT}   >>>>>>>>>>>>>>>>>" << endl;
	Output << "pT-1,<A_{UT}> ={" << A_pT1bg[0] << ", " << A_pT1bg[1] << ", " << A_pT1bg[2] << ", " << A_pT1bg[3] << ", " << A_pT1bg[4] << ", " << A_pT1bg[5] << ", " << A_pT1bg[6] << ", " << A_pT1bg[7] << ", " << A_pT1bg[8] << "}" << endl;
	Output << "pT-1,<Err>={" << deltaA_pT1bg[0] << ", " << deltaA_pT1bg[1] << ", " << deltaA_pT1bg[2] << ", " << deltaA_pT1bg[3] << ", " << deltaA_pT1bg[4] << ", " << deltaA_pT1bg[5] << ", " << deltaA_pT1bg[6] << ", " << deltaA_pT1bg[7] << ", " << deltaA_pT1bg[8] << "}" << endl;
	Output << "pT-2,<A_{UT}>={" << A_pT2bg[0] << ", " << A_pT2bg[1] << ", " << A_pT2bg[2] << ", " << A_pT2bg[3] << ", " << A_pT2bg[4] << ", " << A_pT2bg[5] << ", " << A_pT2bg[6] << ", " << A_pT2bg[7] << ", " << A_pT2bg[8] << "}" << endl;
	Output << "pT-2,<Err> ={" << deltaA_pT2bg[0] << ", " << deltaA_pT2bg[1] << ", " << deltaA_pT2bg[2] << ", " << deltaA_pT2bg[3] << ", " << deltaA_pT2bg[4] << ", " << deltaA_pT2bg[5] << ", " << deltaA_pT2bg[6] << ", " << deltaA_pT2bg[7] << ", " << deltaA_pT2bg[8] << "}" << endl;
	Output << "pT-3,<A_{UT}> ={" << A_pT3bg[0] << ", " << A_pT3bg[1] << ", " << A_pT3bg[2] << ", " << A_pT3bg[3] << ", " << A_pT3bg[4] << ", " << A_pT3bg[5] << ", " << A_pT3bg[6] << ", " << A_pT3bg[7] << ", " << A_pT3bg[8] << "}" << endl;
	Output << "pT-3,<Err> ={" << deltaA_pT3bg[0] << ", " << deltaA_pT3bg[1] << ", " << deltaA_pT3bg[2] << ", " << deltaA_pT3bg[3] << ", " << deltaA_pT3bg[4] << ", " << deltaA_pT3bg[5] << ", " << deltaA_pT3bg[6] << ", " << deltaA_pT3bg[7] << ", " << deltaA_pT3bg[8] << "}" << endl;
	Output << "pT-4,<A_{UT}> ={" << A_pT4bg[0] << ", " << A_pT4bg[1] << ", " << A_pT4bg[2] << ", " << A_pT4bg[3] << ", " << A_pT4bg[4] << ", " << A_pT4bg[5] << ", " << A_pT4bg[6] << ", " << A_pT4bg[7] << ", " << A_pT4bg[8] << "}" << endl;
	Output << "pT-4,<Err> ={" << deltaA_pT4bg[0] << ", " << deltaA_pT4bg[1] << ", " << deltaA_pT4bg[2] << ", " << deltaA_pT4bg[3] << ", " << deltaA_pT4bg[4] << ", " << deltaA_pT4bg[5] << ", " << deltaA_pT4bg[6] << ", " << deltaA_pT4bg[7] << ", " << deltaA_pT4bg[8] << "}" << endl;
	Output << "pT-5,<A_{UT}> ={" << A_pT5bg[0] << ", " << A_pT5bg[1] << ", " << A_pT5bg[2] << ", " << A_pT5bg[3] << ", " << A_pT5bg[4] << ", " << A_pT5bg[5] << ", " << A_pT5bg[6] << ", " << A_pT5bg[7] << ", " << A_pT5bg[8] << "}" << endl;
	Output << "pT-5,<Err> ={" << deltaA_pT5bg[0] << ", " << deltaA_pT5bg[1] << ", " << deltaA_pT5bg[2] << ", " << deltaA_pT5bg[3] << ", " << deltaA_pT5bg[4] << ", " << deltaA_pT5bg[5] << ", " << deltaA_pT5bg[6] << ", " << deltaA_pT5bg[7] << ", " << deltaA_pT5bg[8] << "}" << endl;

	Output << "                                                                                                  " << endl;
	Output << "                                                                                                  " << endl;
	Output << "                                                                                                  " << endl;
	Output << "                                      YELLOW BEAM                                                   " << endl;
	Output << "********************   FINAL ASYMMETRY , CONE > 0.7, Avg Minv  Binning for Yellow BEAM Only****************" << endl;

	Output << "<<<<<<<<<<<<<<<   Minv Bin Averages for asymemtry plot  >>>>>>>>>>>>>>" << endl;
	Output << "pT-1, MinvYGt[0][9]= {" << MinvYGt[0][0] << ", " << MinvYGt[0][1] << ", " << MinvYGt[0][2] << ", " << MinvYGt[0][3] << ", " << MinvYGt[0][4] << ", " << MinvYGt[0][5] << ", " << MinvYGt[0][6] << ", " << MinvYGt[0][7] << ", " << MinvYGt[0][8] << "}" << endl;
	Output << "pT-2, MinvYGt[1][9]= {" << MinvYGt[1][0] << ", " << MinvYGt[1][1] << ", " << MinvYGt[1][2] << ", " << MinvYGt[1][3] << ", " << MinvYGt[1][4] << ", " << MinvYGt[1][5] << ", " << MinvYGt[1][6] << ", " << MinvYGt[1][7] << ", " << MinvYGt[1][8] << "}" << endl;
	Output << "pT-3, MinvYGt[2][9]= {" << MinvYGt[2][0] << ", " << MinvYGt[2][1] << ", " << MinvYGt[2][2] << ", " << MinvYGt[2][3] << ", " << MinvYGt[2][4] << ", " << MinvYGt[2][5] << ", " << MinvYGt[2][6] << ", " << MinvYGt[2][7] << ", " << MinvYGt[2][8] << "}" << endl;
	Output << "pT-4, MinvYGt[3][9]= {" << MinvYGt[3][0] << ", " << MinvYGt[3][1] << ", " << MinvYGt[3][2] << ", " << MinvYGt[3][3] << ", " << MinvYGt[3][4] << ", " << MinvYGt[3][5] << ", " << MinvYGt[3][6] << ", " << MinvYGt[3][7] << ", " << MinvYGt[3][8] << "}" << endl;
	Output << "pT-5, MinvYGt[4][9]= {" << MinvYGt[4][0] << ", " << MinvYGt[4][1] << ", " << MinvYGt[4][2] << ", " << MinvYGt[4][3] << ", " << MinvYGt[4][4] << ", " << MinvYGt[4][5] << ", " << MinvYGt[4][6] << ", " << MinvYGt[4][7] << ", " << MinvYGt[4][8] << "}" << endl;

	Output << "<<<<<<<<<<<<<<<   Minv Bin Averages for asymemtry plot  >>>>>>>>>>>>>>" << endl;
	Output << "pT-1, MinvYLt[0][9]= {" << MinvYLt[0][0] << ", " << MinvYLt[0][1] << ", " << MinvYLt[0][2] << ", " << MinvYLt[0][3] << ", " << MinvYLt[0][4] << ", " << MinvYLt[0][5] << ", " << MinvYLt[0][6] << ", " << MinvYLt[0][7] << ", " << MinvYLt[0][8] << "}" << endl;
	Output << "pT-2, MinvYLt[1][9]= {" << MinvYLt[1][0] << ", " << MinvYLt[1][1] << ", " << MinvYLt[1][2] << ", " << MinvYLt[1][3] << ", " << MinvYLt[1][4] << ", " << MinvYLt[1][5] << ", " << MinvYLt[1][6] << ", " << MinvYLt[1][7] << ", " << MinvYLt[1][8] << "}" << endl;
	Output << "pT-3, MinvYLt[2][9]= {" << MinvYLt[2][0] << ", " << MinvYLt[2][1] << ", " << MinvYLt[2][2] << ", " << MinvYLt[2][3] << ", " << MinvYLt[2][4] << ", " << MinvYLt[2][5] << ", " << MinvYLt[2][6] << ", " << MinvYLt[2][7] << ", " << MinvYLt[2][8] << "}" << endl;
	Output << "pT-4, MinvYLt[3][9]= {" << MinvYLt[3][0] << ", " << MinvYLt[3][1] << ", " << MinvYLt[3][2] << ", " << MinvYLt[3][3] << ", " << MinvYLt[3][4] << ", " << MinvYLt[3][5] << ", " << MinvYLt[3][6] << ", " << MinvYLt[3][7] << ", " << MinvYLt[3][8] << "}" << endl;
	Output << "pT-5, MinvYLt[4][9]= {" << MinvYLt[4][0] << ", " << MinvYLt[4][1] << ", " << MinvYLt[4][2] << ", " << MinvYLt[4][3] << ", " << MinvYLt[4][4] << ", " << MinvYLt[4][5] << ", " << MinvYLt[4][6] << ", " << MinvYLt[4][7] << ", " << MinvYLt[4][8] << "}" << endl;

	Output << "                                                                                                  " << endl;
	Output << "                                                                                                  " << endl;
	Output << "<<<<<<<<<<<       Eta < 0, YELLOW BEAM A_{UT}   >>>>>>>>>>>>>>>>>" << endl;
	Output << "pT-1,<A_{UT}> ={" << A_pT1yl[0] << ", " << A_pT1yl[1] << ", " << A_pT1yl[2] << ", " << A_pT1yl[3] << ", " << A_pT1yl[4] << ", " << A_pT1yl[5] << ", " << A_pT1yl[6] << ", " << A_pT1yl[7] << ", " << A_pT1yl[8] << "}" << endl;
	Output << "pT-1,<Err>={" << deltaA_pT1yl[0] << ", " << deltaA_pT1yl[1] << ", " << deltaA_pT1yl[2] << ", " << deltaA_pT1yl[3] << ", " << deltaA_pT1yl[4] << ", " << deltaA_pT1yl[5] << ", " << deltaA_pT1yl[6] << ", " << deltaA_pT1yl[7] << ", " << deltaA_pT1yl[8] << "}" << endl;
	Output << "pT-2,<A_{UT}>={" << A_pT2yl[0] << ", " << A_pT2yl[1] << ", " << A_pT2yl[2] << ", " << A_pT2yl[3] << ", " << A_pT2yl[4] << ", " << A_pT2yl[5] << ", " << A_pT2yl[6] << ", " << A_pT2yl[7] << ", " << A_pT2yl[8] << "}" << endl;
	Output << "pT-2,<Err> ={" << deltaA_pT2yl[0] << ", " << deltaA_pT2yl[1] << ", " << deltaA_pT2yl[2] << ", " << deltaA_pT2yl[3] << ", " << deltaA_pT2yl[4] << ", " << deltaA_pT2yl[5] << ", " << deltaA_pT2yl[6] << ", " << deltaA_pT2yl[7] << ", " << deltaA_pT2yl[8] << "}" << endl;
	Output << "pT-3,<A_{UT}> ={" << A_pT3yl[0] << ", " << A_pT3yl[1] << ", " << A_pT3yl[2] << ", " << A_pT3yl[3] << ", " << A_pT3yl[4] << ", " << A_pT3yl[5] << ", " << A_pT3yl[6] << ", " << A_pT3yl[7] << ", " << A_pT3yl[8] << "}" << endl;
	Output << "pT-3,<Err> ={" << deltaA_pT3yl[0] << ", " << deltaA_pT3yl[1] << ", " << deltaA_pT3yl[2] << ", " << deltaA_pT3yl[3] << ", " << deltaA_pT3yl[4] << ", " << deltaA_pT3yl[5] << ", " << deltaA_pT3yl[6] << ", " << deltaA_pT3yl[7] << ", " << deltaA_pT3yl[8] << "}" << endl;
	Output << "pT-4,<A_{UT}> ={" << A_pT4yl[0] << ", " << A_pT4yl[1] << ", " << A_pT4yl[2] << ", " << A_pT4yl[3] << ", " << A_pT4yl[4] << ", " << A_pT4yl[5] << ", " << A_pT4yl[6] << ", " << A_pT4yl[7] << ", " << A_pT4yl[8] << "}" << endl;
	Output << "pT-4,<Err> ={" << deltaA_pT4yl[0] << ", " << deltaA_pT4yl[1] << ", " << deltaA_pT4yl[2] << ", " << deltaA_pT4yl[3] << ", " << deltaA_pT4yl[4] << ", " << deltaA_pT4yl[5] << ", " << deltaA_pT4yl[6] << ", " << deltaA_pT4yl[7] << ", " << deltaA_pT4yl[8] << "}" << endl;
	Output << "pT-5,<A_{UT}> ={" << A_pT5yl[0] << ", " << A_pT5yl[1] << ", " << A_pT5yl[2] << ", " << A_pT5yl[3] << ", " << A_pT5yl[4] << ", " << A_pT5yl[5] << ", " << A_pT5yl[6] << ", " << A_pT5yl[7] << ", " << A_pT5yl[8] << "}" << endl;
	Output << "pT-5,<Err> ={" << deltaA_pT5yl[0] << ", " << deltaA_pT5yl[1] << ", " << deltaA_pT5yl[2] << ", " << deltaA_pT5yl[3] << ", " << deltaA_pT5yl[4] << ", " << deltaA_pT5yl[5] << ", " << deltaA_pT5yl[6] << ", " << deltaA_pT5yl[7] << ", " << deltaA_pT5yl[8] << "}" << endl;

	Output << "                                                                                                  " << endl;
	Output << "                                                                                                  " << endl;
	Output << "<<<<<<<<<<<       Eta > 0, BLUE BEAM A_{UT}   >>>>>>>>>>>>>>>>>" << endl;
	Output << "pT-1,<A_{UT}> ={" << A_pT1yg[0] << ", " << A_pT1yg[1] << ", " << A_pT1yg[2] << ", " << A_pT1yg[3] << ", " << A_pT1yg[4] << ", " << A_pT1yg[5] << ", " << A_pT1yg[6] << ", " << A_pT1yg[7] << ", " << A_pT1yg[8] << "}" << endl;
	Output << "pT-1,<Err>={" << deltaA_pT1yg[0] << ", " << deltaA_pT1yg[1] << ", " << deltaA_pT1yg[2] << ", " << deltaA_pT1yg[3] << ", " << deltaA_pT1yg[4] << ", " << deltaA_pT1yg[5] << ", " << deltaA_pT1yg[6] << ", " << deltaA_pT1yg[7] << ", " << deltaA_pT1yg[8] << "}" << endl;
	Output << "pT-2,<A_{UT}>={" << A_pT2yg[0] << ", " << A_pT2yg[1] << ", " << A_pT2yg[2] << ", " << A_pT2yg[3] << ", " << A_pT2yg[4] << ", " << A_pT2yg[5] << ", " << A_pT2yg[6] << ", " << A_pT2yg[7] << ", " << A_pT2yg[8] << "}" << endl;
	Output << "pT-2,<Err> ={" << deltaA_pT2yg[0] << ", " << deltaA_pT2yg[1] << ", " << deltaA_pT2yg[2] << ", " << deltaA_pT2yg[3] << ", " << deltaA_pT2yg[4] << ", " << deltaA_pT2yg[5] << ", " << deltaA_pT2yg[6] << ", " << deltaA_pT2yg[7] << ", " << deltaA_pT2yg[8] << "}" << endl;
	Output << "pT-3,<A_{UT}> ={" << A_pT3yg[0] << ", " << A_pT3yg[1] << ", " << A_pT3yg[2] << ", " << A_pT3yg[3] << ", " << A_pT3yg[4] << ", " << A_pT3yg[5] << ", " << A_pT3yg[6] << ", " << A_pT3yg[7] << ", " << A_pT3yg[8] << "}" << endl;
	Output << "pT-3,<Err> ={" << deltaA_pT3yg[0] << ", " << deltaA_pT3yg[1] << ", " << deltaA_pT3yg[2] << ", " << deltaA_pT3yg[3] << ", " << deltaA_pT3yg[4] << ", " << deltaA_pT3yg[5] << ", " << deltaA_pT3yg[6] << ", " << deltaA_pT3yg[7] << ", " << deltaA_pT3yg[8] << "}" << endl;
	Output << "pT-4,<A_{UT}> ={" << A_pT4yg[0] << ", " << A_pT4yg[1] << ", " << A_pT4yg[2] << ", " << A_pT4yg[3] << ", " << A_pT4yg[4] << ", " << A_pT4yg[5] << ", " << A_pT4yg[6] << ", " << A_pT4yg[7] << ", " << A_pT4yg[8] << "}" << endl;
	Output << "pT-4,<Err> ={" << deltaA_pT4yg[0] << ", " << deltaA_pT4yg[1] << ", " << deltaA_pT4yg[2] << ", " << deltaA_pT4yg[3] << ", " << deltaA_pT4yg[4] << ", " << deltaA_pT4yg[5] << ", " << deltaA_pT4yg[6] << ", " << deltaA_pT4yg[7] << ", " << deltaA_pT4yg[8] << "}" << endl;
	Output << "pT-5,<A_{UT}> ={" << A_pT5yg[0] << ", " << A_pT5yg[1] << ", " << A_pT5yg[2] << ", " << A_pT5yg[3] << ", " << A_pT5yg[4] << ", " << A_pT5yg[5] << ", " << A_pT5yg[6] << ", " << A_pT5yg[7] << ", " << A_pT5yg[8] << "}" << endl;
	Output << "pT-5,<Err> ={" << deltaA_pT5yg[0] << ", " << deltaA_pT5yg[1] << ", " << deltaA_pT5yg[2] << ", " << deltaA_pT5yg[3] << ", " << deltaA_pT5yg[4] << ", " << deltaA_pT5yg[5] << ", " << deltaA_pT5yg[6] << ", " << deltaA_pT5yg[7] << ", " << deltaA_pT5yg[8] << "}" << endl;

	Output << "                                                                                                  " << endl;
	Output << "                                                                                                  " << endl;
	Output << "                                                                                                  " << endl;
	Output << "                                      AVERAGE A_UT                                                    " << endl;
	Output << "********************   FINAL ASYMMETRY , CONE > 0.7, Avg Minv  Binning for BLUE and YELLOW BEAM ****************" << endl;

	Output << "<<<<<<<<<<<<<<<   Minv Bin Averages for asymemtry plot  >>>>>>>>>>>>>>" << endl;
	Output << "pT-1 MinvAvgGt[0][9]={" << MinvAvgGt[0][0] << "," << MinvAvgGt[0][1] << "," << MinvAvgGt[0][2] << "," << MinvAvgGt[0][3] << "," << MinvAvgGt[0][4] << "," << MinvAvgGt[0][5] << "," << MinvAvgGt[0][6] << "," << MinvAvgGt[0][7] << "," << MinvAvgGt[0][8] << "}" << endl;
	Output << "pT-2 MinvAvgGt[1][9]={" << MinvAvgGt[1][0] << "," << MinvAvgGt[1][1] << "," << MinvAvgGt[1][2] << "," << MinvAvgGt[1][3] << "," << MinvAvgGt[1][4] << "," << MinvAvgGt[1][5] << "," << MinvAvgGt[1][6] << "," << MinvAvgGt[1][7] << "," << MinvAvgGt[1][8] << "}" << endl;
	Output << "pT-3 MinvAvgGt[2][9]={" << MinvAvgGt[2][0] << "," << MinvAvgGt[2][1] << "," << MinvAvgGt[2][2] << "," << MinvAvgGt[2][3] << "," << MinvAvgGt[2][4] << "," << MinvAvgGt[2][5] << "," << MinvAvgGt[2][6] << "," << MinvAvgGt[2][7] << "," << MinvAvgGt[2][8] << "}" << endl;
	Output << "pT-4 MinvAvgGt[3][9]={" << MinvAvgGt[3][0] << "," << MinvAvgGt[3][1] << "," << MinvAvgGt[3][2] << "," << MinvAvgGt[3][3] << "," << MinvAvgGt[3][4] << "," << MinvAvgGt[3][5] << "," << MinvAvgGt[3][6] << "," << MinvAvgGt[3][7] << "," << MinvAvgGt[3][8] << "}" << endl;
	Output << "pT-5 MinvAvgGt[4][9]={" << MinvAvgGt[4][0] << "," << MinvAvgGt[4][1] << "," << MinvAvgGt[4][2] << "," << MinvAvgGt[4][3] << "," << MinvAvgGt[4][4] << "," << MinvAvgGt[4][5] << "," << MinvAvgGt[4][6] << "," << MinvAvgGt[4][7] << "," << MinvAvgGt[4][8] << "}" << endl;

	Output << "pT-1 MinvAvgLt[0][9]={" << MinvAvgLt[0][0] << "," << MinvAvgLt[0][1] << "," << MinvAvgLt[0][2] << "," << MinvAvgLt[0][3] << "," << MinvAvgLt[0][4] << "," << MinvAvgLt[0][5] << "," << MinvAvgLt[0][6] << "," << MinvAvgLt[0][7] << "," << MinvAvgLt[0][8] << "}" << endl;
	Output << "pT-2 MinvAvgLt[1][9]={" << MinvAvgLt[1][0] << "," << MinvAvgLt[1][1] << "," << MinvAvgLt[1][2] << "," << MinvAvgLt[1][3] << "," << MinvAvgLt[1][4] << "," << MinvAvgLt[1][5] << "," << MinvAvgLt[1][6] << "," << MinvAvgLt[1][7] << "," << MinvAvgLt[1][8] << "}" << endl;
	Output << "pT-3 MinvAvgLt[2][9]={" << MinvAvgLt[2][0] << "," << MinvAvgLt[2][1] << "," << MinvAvgLt[2][2] << "," << MinvAvgLt[2][3] << "," << MinvAvgLt[2][4] << "," << MinvAvgLt[2][5] << "," << MinvAvgLt[2][6] << "," << MinvAvgLt[2][7] << "," << MinvAvgLt[2][8] << "}" << endl;
	Output << "pT-4 MinvAvgLt[3][9]={" << MinvAvgLt[3][0] << "," << MinvAvgLt[3][1] << "," << MinvAvgLt[3][2] << "," << MinvAvgLt[3][3] << "," << MinvAvgLt[3][4] << "," << MinvAvgLt[3][5] << "," << MinvAvgLt[3][6] << "," << MinvAvgLt[3][7] << "," << MinvAvgLt[3][8] << "}" << endl;
	Output << "pT-5 MinvAvgLt[4][9]={" << MinvAvgLt[4][0] << "," << MinvAvgLt[4][1] << "," << MinvAvgLt[4][2] << "," << MinvAvgLt[4][3] << "," << MinvAvgLt[4][4] << "," << MinvAvgLt[4][5] << "," << MinvAvgLt[4][6] << "," << MinvAvgLt[4][7] << "," << MinvAvgLt[4][8] << "}" << endl;

	Output << "                                                                                                  " << endl;
	Output << "                                                                                                  " << endl;
	Output << "<<<<<<<<<<<       Eta < 0, Average A_{UT}   >>>>>>>>>>>>>>>>>" << endl;
	Output << "pT-1,<A_{UT}> ={" << avgA_pT1l[0] << ", " << avgA_pT1l[1] << ", " << avgA_pT1l[2] << ", " << avgA_pT1l[3] << ", " << avgA_pT1l[4] << ", " << avgA_pT1l[5] << ", " << avgA_pT1l[6] << ", " << avgA_pT1l[7] << ", " << avgA_pT1l[8] << "}" << endl;
	Output << "pT-1,<Err>={" << errA_pT1l[0] << ", " << errA_pT1l[1] << ", " << errA_pT1l[2] << ", " << errA_pT1l[3] << ", " << errA_pT1l[4] << ", " << errA_pT1l[5] << ", " << errA_pT1l[6] << ", " << errA_pT1l[7] << ", " << errA_pT1l[8] << "}" << endl;
	Output << "pT-2,<A_{UT}>={" << avgA_pT2l[0] << ", " << avgA_pT2l[1] << ", " << avgA_pT2l[2] << ", " << avgA_pT2l[3] << ", " << avgA_pT2l[4] << ", " << avgA_pT2l[5] << ", " << avgA_pT2l[6] << ", " << avgA_pT2l[7] << ", " << avgA_pT2l[8] << "}" << endl;
	Output << "pT-2,<Err> ={" << errA_pT2l[0] << ", " << errA_pT2l[1] << ", " << errA_pT2l[2] << ", " << errA_pT2l[3] << ", " << errA_pT2l[4] << ", " << errA_pT2l[5] << ", " << errA_pT2l[6] << ", " << errA_pT2l[7] << ", " << errA_pT2l[8] << "}" << endl;
	Output << "pT-3,<A_{UT}> ={" << avgA_pT3l[0] << ", " << avgA_pT3l[1] << ", " << avgA_pT3l[2] << ", " << avgA_pT3l[3] << ", " << avgA_pT3l[4] << ", " << avgA_pT3l[5] << ", " << avgA_pT3l[6] << ", " << avgA_pT3l[7] << ", " << avgA_pT3l[8] << "}" << endl;
	Output << "pT-3,<Err> ={" << errA_pT3l[0] << ", " << errA_pT3l[1] << ", " << errA_pT3l[2] << ", " << errA_pT3l[3] << ", " << errA_pT3l[4] << ", " << errA_pT3l[5] << ", " << errA_pT3l[6] << ", " << errA_pT3l[7] << ", " << errA_pT3l[8] << "}" << endl;
	Output << "pT-4,<A_{UT}> ={" << avgA_pT4l[0] << ", " << avgA_pT4l[1] << ", " << avgA_pT4l[2] << ", " << avgA_pT4l[3] << ", " << avgA_pT4l[4] << ", " << avgA_pT4l[5] << ", " << avgA_pT4l[6] << ", " << avgA_pT4l[7] << ", " << avgA_pT4l[8] << "}" << endl;
	Output << "pT-4,<Err> ={" << errA_pT4l[0] << ", " << errA_pT4l[1] << ", " << errA_pT4l[2] << ", " << errA_pT4l[3] << ", " << errA_pT4l[4] << ", " << errA_pT4l[5] << ", " << errA_pT4l[6] << ", " << errA_pT4l[7] << ", " << errA_pT4l[8] << "}" << endl;
	Output << "pT-5,<A_{UT}> ={" << avgA_pT5l[0] << ", " << avgA_pT5l[1] << ", " << avgA_pT5l[2] << ", " << avgA_pT5l[3] << ", " << avgA_pT5l[4] << ", " << avgA_pT5l[5] << ", " << avgA_pT5l[6] << ", " << avgA_pT5l[7] << ", " << avgA_pT5l[8] << "}" << endl;
	Output << "pT-5,<Err> ={" << errA_pT5l[0] << ", " << errA_pT5l[1] << ", " << errA_pT5l[2] << ", " << errA_pT5l[3] << ", " << errA_pT5l[4] << ", " << errA_pT5l[5] << ", " << errA_pT5l[6] << ", " << errA_pT5l[7] << ", " << errA_pT5l[8] << "}" << endl;

	Output << "                                                                                                  " << endl;
	Output << "                                                                                                  " << endl;
	Output << "<<<<<<<<<<<       Eta > 0, Average A_{UT}   >>>>>>>>>>>>>>>>>" << endl;
	Output << "pT-1,<A_{UT}> ={" << avgA_pT1g[0] << ", " << avgA_pT1g[1] << ", " << avgA_pT1g[2] << ", " << avgA_pT1g[3] << ", " << avgA_pT1g[4] << ", " << avgA_pT1g[5] << ", " << avgA_pT1g[6] << ", " << avgA_pT1g[7] << ", " << avgA_pT1g[8] << "}" << endl;
	Output << "pT-1,<Err> ={" << errA_pT1g[0] << ", " << errA_pT1g[1] << ", " << errA_pT1g[2] << ", " << errA_pT1g[3] << ", " << errA_pT1g[4] << ", " << errA_pT1g[5] << ", " << errA_pT1g[6] << ", " << errA_pT1g[7] << ", " << errA_pT1g[8] << "}" << endl;
	Output << "pT-2,<A_{UT}> ={" << avgA_pT2g[0] << ", " << avgA_pT2g[1] << ", " << avgA_pT2g[2] << ", " << avgA_pT2g[3] << ", " << avgA_pT2g[4] << ", " << avgA_pT2g[5] << ", " << avgA_pT2g[6] << ", " << avgA_pT2g[7] << ", " << avgA_pT2g[8] << "}" << endl;
	Output << "pT-2,<Err> ={" << errA_pT2g[0] << ", " << errA_pT2g[1] << ", " << errA_pT2g[2] << ", " << errA_pT2g[3] << ", " << errA_pT2g[4] << ", " << errA_pT2g[5] << ", " << errA_pT2g[6] << ", " << errA_pT2g[7] << ", " << errA_pT2g[8] << "}" << endl;
	Output << "pT-3,<A_{UT}> ={" << avgA_pT3g[0] << ", " << avgA_pT3g[1] << ", " << avgA_pT3g[2] << ", " << avgA_pT3g[3] << ", " << avgA_pT3g[4] << ", " << avgA_pT3g[5] << ", " << avgA_pT3g[6] << ", " << avgA_pT3g[7] << ", " << avgA_pT3g[8] << "}" << endl;
	Output << "pT-3,<Err> ={" << errA_pT3g[0] << ", " << errA_pT3g[1] << ", " << errA_pT3g[2] << ", " << errA_pT3g[3] << ", " << errA_pT3g[4] << ", " << errA_pT3g[5] << ", " << errA_pT3g[6] << ", " << errA_pT3g[7] << ", " << errA_pT3g[8] << "}" << endl;
	Output << "pT-4,<A_{UT}> ={" << avgA_pT4g[0] << ", " << avgA_pT4g[1] << ", " << avgA_pT4g[2] << ", " << avgA_pT4g[3] << ", " << avgA_pT4g[4] << ", " << avgA_pT4g[5] << ", " << avgA_pT4g[6] << ", " << avgA_pT4g[7] << ", " << avgA_pT4g[8] << "}" << endl;
	Output << "pT-4,<Err> ={" << errA_pT4g[0] << ", " << errA_pT4g[1] << ", " << errA_pT4g[2] << ", " << errA_pT4g[3] << ", " << errA_pT4g[4] << ", " << errA_pT4g[5] << ", " << errA_pT4g[6] << ", " << errA_pT4g[7] << ", " << errA_pT4g[8] << "}" << endl;
	Output << "pT-5,<A_{UT}> ={" << avgA_pT5g[0] << ", " << avgA_pT5g[1] << ", " << avgA_pT5g[2] << ", " << avgA_pT5g[3] << ", " << avgA_pT5g[4] << ", " << avgA_pT5g[5] << ", " << avgA_pT5g[6] << ", " << avgA_pT5g[7] << ", " << avgA_pT5g[8] << "}" << endl;
	Output << "pT-5,<Err> ={" << errA_pT5g[0] << ", " << errA_pT5g[1] << ", " << errA_pT5g[2] << ", " << errA_pT5g[3] << ", " << errA_pT5g[4] << ", " << errA_pT5g[5] << ", " << errA_pT5g[6] << ", " << errA_pT5g[7] << ", " << errA_pT5g[8] << "}" << endl;

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

	double AvgM_G[9] = {0};
	double AvgM_L[9] = {0};
	double M_BG[9] = {0};
	double M_BL[9] = {0};
	double M_YG[9] = {0};
	double M_YL[9] = {0};

	vector<double> MINV[3][2];

	vector<double> A_UT[3][2];
	vector<double> dA_UT[3][2];
	vector<double> PIDSYS[3][2];

	vector<double> WA_UT[3][2];
	vector<double> WdA_UT[3][2];
	vector<double> WPIDSYS[3][2];

	vector<double> WPIDSYS_TPC[3][2];
	vector<double> WPIDSYS_TPCorTOF[3][2];
	vector<double> WPIDSYS_TOF[3][2];
	vector<double> WTOTSYS[3][2];

	TGraphErrors *GRAPH[3][2][5];
	TGraphErrors *WGRAPH[3][2][5];
	TGraphErrors *gr_sys[3][2][5];
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

	// To Compare with Run11
	double M_R11_G_P1[9] = {0.4, 0.5, 0.6, 0.7, 1.0, 1.3};
	double A_R11_G_P1[9] = {0.0041, 0.0057, 0.0055, 0.0045, 0.0073, 0.0126};
	double dA_R11_G_P1[9] = {0.0065, 0.0069, 0.0071, 0.0065, 0.0055, 0.0068};
	double M_R11_G_P2[9] = {0.4, 0.5, 0.6, 0.7, 0.9, 1.3, 1.7};
	double A_R11_G_P2[9] = {0.0047, 0.0009, 0.0140, -0.0028, 0.0043, 0.0066, 0.0000};
	double dA_R11_G_P2[9] = {0.0060, 0.0061, 0.0064, 0.0057, 0.0049, 0.0061, 0.0113};
	double M_R11_G_P3[9] = {0.4, 0.5, 0.6, 0.7, 0.9, 1.3, 1.9};
	double A_R11_G_P3[9] = {0.0070, 0.0163, 0.0005, 0.0135, 0.0059, -0.0086, -0.0103};
	double dA_R11_G_P3[9] = {0.0063, 0.0062, 0.0064, 0.0057, 0.0049, 0.0062, 0.0081};
	double M_R11_G_P4[9] = {0.4, 0.5, 0.6, 0.7, 0.9, 1.3, 2.0};
	double A_R11_G_P4[9] = {-0.0038, -0.0056, 0.0107, 0.0055, 0.0204, 0.0068, 0.0039};
	double dA_R11_G_P4[9] = {0.0070, 0.0065, 0.0066, 0.0057, 0.0049, 0.0062, 0.0069};
	double M_R11_G_P5[9] = {0.4, 0.5, 0.6, 0.7, 0.9, 1.3, 2.2};
	double A_R11_G_P5[9] = {0.0023, 0.0269, 0.0207, 0.0249, 0.0298, 0.0196, 0.0104};
	double dA_R11_G_P5[9] = {0.0090, 0.0075, 0.0072, 0.0060, 0.0050, 0.0060, 0.0061};
	double M_R11_G[9] = {0};
	double A_R11_G[9] = {0};
	double dA_R11_G[9] = {0};

	for (int pad = 0; pad < 5; pad++)
	{
		if (pad == 0)
		{
			for (int mbin = 0; mbin < 9; mbin++)
			{

				M_R11_G[mbin] = M_R11_G_P1[mbin];
				A_R11_G[mbin] = A_R11_G_P1[mbin];
				dA_R11_G[mbin] = dA_R11_G_P1[mbin];
			}
		}
		if (pad == 1)
		{
			for (int mbin = 0; mbin < 9; mbin++)
			{

				M_R11_G[mbin] = M_R11_G_P2[mbin];
				A_R11_G[mbin] = A_R11_G_P2[mbin];
				dA_R11_G[mbin] = dA_R11_G_P2[mbin];
			}
		}
		if (pad == 2)
		{
			for (int mbin = 0; mbin < 9; mbin++)
			{

				M_R11_G[mbin] = M_R11_G_P3[mbin];
				A_R11_G[mbin] = A_R11_G_P3[mbin];
				dA_R11_G[mbin] = dA_R11_G_P3[mbin];
			}
		}
		if (pad == 3)
		{
			for (int mbin = 0; mbin < 9; mbin++)
			{

				M_R11_G[mbin] = M_R11_G_P4[mbin];
				A_R11_G[mbin] = A_R11_G_P4[mbin];
				dA_R11_G[mbin] = dA_R11_G_P4[mbin];
			}
		}
		if (pad == 4)
		{
			for (int mbin = 0; mbin < 9; mbin++)
			{

				M_R11_G[mbin] = M_R11_G_P5[mbin];
				A_R11_G[mbin] = A_R11_G_P5[mbin];
				dA_R11_G[mbin] = dA_R11_G_P5[mbin];
			}
		}

		GRAPH_R11[pad] = new TGraphErrors(9, M_R11_G, A_R11_G, 0, dA_R11_G);
	}

	//=============================Compare with JAMDiFF theory ============================

	TGraphErrors *GRAPH_JAMDiFF_G[5];
	TGraph *GRAPH_JAMDiFF_CentralVal_G[5];

	TGraphErrors *GRAPH_JAMDiFF_L[5];
	TGraph *GRAPH_JAMDiFF_CentralVal_L[5];

	double M_JAMDiFF_G[9];
	double A_JAMDiFF_G[9];
	double E_JAMDiFF_G[9];

	double M_JAMDiFF_L[9];
	double A_JAMDiFF_L[9];
	double E_JAMDiFF_L[9];

	double JAM_M510_P1_G[9] = {0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.3};
	double JAM_M510_P2_G[9] = {0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1, 1.2, 1.6};
	double JAM_M510_P3_G[9] = {0.4, 0.5, 0.5, 0.6, 0.7, 0.8, 1, 1.3, 1.8};
	double JAM_M510_P4_G[8] = {0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.4};
	double JAM_M510_P5_G[8] = {0.4, 0.5, 0.6, 0.7, 0.8, 1, 1.2, 1.5};

	double JAM_A510_P1_G[9] = {0.00124005, 0.00127839, 0.002220933, 0.003932849, 0.005467239, 0.006989092, 0.008221294, 0.007738282, 0.006514535};
	double JAM_A510_P2_G[9] = {0.001321984, 0.001415474, 0.002474664, 0.004365926, 0.005875476, 0.00749873, 0.009124923, 0.008192817, 0.006917135};
	double JAM_A510_P3_G[9] = {0.001823855, 0.003205427, 0.003152915, 0.00562829, 0.007510655, 0.009553294, 0.011646116, 0.009558897, 0.004426999};
	double JAM_A510_P4_G[8] = {0.002481064, 0.004377054, 0.007792102, 0.010321801, 0.013352078, 0.016472733, 0.015696895, 0.014135812};
	double JAM_A510_P5_G[8] = {0.004472185, 0.007788241, 0.01420061, 0.019224171, 0.024984835, 0.033189751, 0.029596744, 0.029117996};

	double JAM_E510_P1_G[9] = {0.001063449, 0.000750701, 0.001039825, 0.001741493, 0.002067357, 0.002396818, 0.002945778, 0.003064155, 0.002629226};
	double JAM_E510_P2_G[9] = {0.001071201, 0.000773406, 0.001071701, 0.001777762, 0.002052738, 0.002371788, 0.003024682, 0.002806216, 0.003963517};
	double JAM_E510_P3_G[9] = {0.000918129, 0.001265935, 0.001247165, 0.002084682, 0.002394783, 0.002800845, 0.003544317, 0.003368451, 0.002500325};
	double JAM_E510_P4_G[8] = {0.001112399, 0.001504264, 0.002502025, 0.002849374, 0.00344219, 0.004328849, 0.004577453, 0.00582897};
	double JAM_E510_P5_G[8] = {0.001572482, 0.00196801, 0.003276316, 0.003631747, 0.00438972, 0.006166589, 0.006599213, 0.012416455};
	// ======JAM for eta<0======
	double JAM_M510_P1_L[9] = {0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.3};
	double JAM_M510_P2_L[9] = {0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1, 1.2, 1.6};
	double JAM_M510_P3_L[9] = {0.4, 0.5, 0.5, 0.6, 0.7, 0.8, 1, 1.3, 1.8};
	double JAM_M510_P4_L[8] = {0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.4};
	double JAM_M510_P5_L[8] = {0.4, 0.5, 0.6, 0.7, 0.8, 1, 1.2, 1.5};

	double JAM_A510_P1_L[9] = {0.00015529, 0.0001962, 0.000352972, 0.000656859, 0.00094409, 0.001323591, 0.001723591, 0.001773787, 0.001736956};
	double JAM_A510_P2_L[9] = {0.000235613, 0.000283578, 0.000498525, 0.000905394, 0.001278706, 0.001736395, 0.002293429, 0.002149785, 0.002157643};
	double JAM_A510_P3_L[9] = {0.000321068, 0.000563789, 0.000573739, 0.00106234, 0.001491551, 0.002050178, 0.002795331, 0.002409782, 0.001273614};
	double JAM_A510_P4_L[8] = {0.000390638, 0.000687907, 0.001251554, 0.001778813, 0.002396393, 0.003182554, 0.003371645, 0.003170512};
	double JAM_A510_P5_L[8] = {0.00056776, 0.000968409, 0.001816028, 0.002539988, 0.003446136, 0.004948846, 0.004645352, 0.004822397};

	double JAM_E510_P1_L[9] = {0.000149564, 0.00015151, 0.000243897, 0.000441512, 0.000620406, 0.000885953, 0.001214953, 0.00135102, 0.001265564};
	double JAM_E510_P2_L[9] = {0.000195313, 0.000166357, 0.000250285, 0.00045442, 0.00061529, 0.000852319, 0.001320272, 0.001269788, 0.001604628};
	double JAM_E510_P3_L[9] = {0.000159728, 0.000230462, 0.000233942, 0.000426397, 0.000548351, 0.00073869, 0.001212347, 0.001123241, 0.000843578};
	double JAM_E510_P4_L[8] = {0.000174782, 0.000239233, 0.000415314, 0.00050878, 0.000638265, 0.000883973, 0.001126112, 0.001370875};
	double JAM_E510_P5_L[8] = {0.000234959, 0.000296188, 0.000509964, 0.000621307, 0.00080402, 0.001112517, 0.001165571, 0.002126754};

	for (int pad = 0; pad < 5; pad++)
	{
		if (pad == 0)
		{
			for (int mbin = 0; mbin < sizeof(JAM_M510_P1_G) / sizeof(JAM_M510_P1_G[0]); mbin++)
			{

				M_JAMDiFF_G[mbin] = JAM_M510_P1_G[mbin];
				A_JAMDiFF_G[mbin] = JAM_A510_P1_G[mbin];
				E_JAMDiFF_G[mbin] = JAM_E510_P1_G[mbin];

				M_JAMDiFF_L[mbin] = JAM_M510_P1_L[mbin];
				A_JAMDiFF_L[mbin] = JAM_A510_P1_L[mbin];
				E_JAMDiFF_L[mbin] = JAM_E510_P1_L[mbin];
			}
		}
		if (pad == 1)
		{
			for (int mbin = 0; mbin < sizeof(JAM_M510_P2_G) / sizeof(JAM_M510_P2_G[0]); mbin++)
			{

				M_JAMDiFF_G[mbin] = JAM_M510_P2_G[mbin];
				A_JAMDiFF_G[mbin] = JAM_A510_P2_G[mbin];
				E_JAMDiFF_G[mbin] = JAM_E510_P2_G[mbin];

				M_JAMDiFF_L[mbin] = JAM_M510_P2_L[mbin];
				A_JAMDiFF_L[mbin] = JAM_A510_P2_L[mbin];
				E_JAMDiFF_L[mbin] = JAM_E510_P2_L[mbin];
			}
		}
		if (pad == 2)
		{
			for (int mbin = 0; mbin < sizeof(JAM_M510_P3_G) / sizeof(JAM_M510_P3_G[0]); mbin++)
			{

				M_JAMDiFF_G[mbin] = JAM_M510_P3_G[mbin];
				A_JAMDiFF_G[mbin] = JAM_A510_P3_G[mbin];
				E_JAMDiFF_G[mbin] = JAM_E510_P3_G[mbin];

				M_JAMDiFF_L[mbin] = JAM_M510_P3_L[mbin];
				A_JAMDiFF_L[mbin] = JAM_A510_P3_L[mbin];
				E_JAMDiFF_L[mbin] = JAM_E510_P3_L[mbin];
			}
		}
		if (pad == 3)
		{
			for (int mbin = 0; mbin < sizeof(JAM_M510_P4_G) / sizeof(JAM_M510_P4_G[0]); mbin++)
			{

				M_JAMDiFF_G[mbin] = JAM_M510_P4_G[mbin];
				A_JAMDiFF_G[mbin] = JAM_A510_P4_G[mbin];
				E_JAMDiFF_G[mbin] = JAM_E510_P4_G[mbin];

				M_JAMDiFF_L[mbin] = JAM_M510_P4_L[mbin];
				A_JAMDiFF_L[mbin] = JAM_A510_P4_L[mbin];
				E_JAMDiFF_L[mbin] = JAM_E510_P4_L[mbin];
			}
		}
		if (pad == 4)
		{
			for (int mbin = 0; mbin < sizeof(JAM_M510_P5_G) / sizeof(JAM_M510_P5_G[0]); mbin++)
			{

				M_JAMDiFF_G[mbin] = JAM_M510_P5_G[mbin];
				A_JAMDiFF_G[mbin] = JAM_A510_P5_G[mbin];
				E_JAMDiFF_G[mbin] = JAM_E510_P5_G[mbin];

				M_JAMDiFF_L[mbin] = JAM_M510_P5_L[mbin];
				A_JAMDiFF_L[mbin] = JAM_A510_P5_L[mbin];
				E_JAMDiFF_L[mbin] = JAM_E510_P5_L[mbin];
			}
		}

		GRAPH_JAMDiFF_G[pad] = new TGraphErrors(sizeof(M_JAMDiFF_G) / sizeof(M_JAMDiFF_G[0]), M_JAMDiFF_G, A_JAMDiFF_G, 0, E_JAMDiFF_G);
		GRAPH_JAMDiFF_CentralVal_G[pad] = new TGraph(sizeof(M_JAMDiFF_G) / sizeof(M_JAMDiFF_G[0]), M_JAMDiFF_G, A_JAMDiFF_G);

		GRAPH_JAMDiFF_L[pad] = new TGraphErrors(sizeof(M_JAMDiFF_L) / sizeof(M_JAMDiFF_L[0]), M_JAMDiFF_L, A_JAMDiFF_L, 0, E_JAMDiFF_L);
		GRAPH_JAMDiFF_CentralVal_L[pad] = new TGraph(sizeof(M_JAMDiFF_L) / sizeof(M_JAMDiFF_L[0]), M_JAMDiFF_L, A_JAMDiFF_L);

		// JAMDiFF have only 8 points for pad 3 and 4, so I manually reduce the number of points on the graph by substrating by 1 which is not a good practice should use dynamic allocator such as vector here than fix memory array.
		if (pad == 4 || pad == 3)
		{
			GRAPH_JAMDiFF_G[pad] = new TGraphErrors(sizeof(M_JAMDiFF_G) / sizeof(M_JAMDiFF_G[0]) - 1, M_JAMDiFF_G, A_JAMDiFF_G, 0, E_JAMDiFF_G);
			GRAPH_JAMDiFF_CentralVal_G[pad] = new TGraph(sizeof(M_JAMDiFF_G) / sizeof(M_JAMDiFF_G[0]) - 1, M_JAMDiFF_G, A_JAMDiFF_G);

			GRAPH_JAMDiFF_L[pad] = new TGraphErrors(sizeof(M_JAMDiFF_L) / sizeof(M_JAMDiFF_L[0]) - 1, M_JAMDiFF_L, A_JAMDiFF_L, 0, E_JAMDiFF_L);
			GRAPH_JAMDiFF_CentralVal_L[pad] = new TGraph(sizeof(M_JAMDiFF_L) / sizeof(M_JAMDiFF_L[0]) - 1, M_JAMDiFF_L, A_JAMDiFF_L);
		}
	}
	//============================= Compare with JAMDiFF theory End ============================

	// Kaon and Proton seperate
	// double pionPairMGt[5] = {0.817349, 0.822488, 0.825511, 0.791272, 0.728912};
	// double pionPairMLt[5] = {0.8287, 0.819703, 0.821156, 0.791168, 0.730723};
	// Kaon and Proton Combine
	double pionPairMGt[5] = {0.838334, 0.842922, 0.845917, 0.813829, 0.780844};
	double pionPairMLt[5] = {0.839661, 0.845568, 0.848364, 0.816115, 0.783557};
	double pionPairPtGt[5] = {0.833136, 0.854063, 0.86806, 0.869287, 0.83688};
	double pionPairPtLt[5] = {0.833513, 0.85424, 0.867603, 0.8696, 0.840349};

	// Kaon and Proton Combine P20ic.SL22b TPC, TPCorTOF and TOF pid
	// double pionPairMGt_TPC[5] = {0.877089, 0.887985, 0.890326, 0.864459, 0.837882};
	// double pionPairMLt_TPC[5] = {0.869475, 0.877599, 0.881896, 0.85364, 0.828147};
	// Final TPC result
	double pionPairPtGt_TPC[5] = {0.819141, 0.834208, 0.851712, 0.85491, 0.824758};
	double pionPairPtLt_TPC[5] = {0.826245, 0.83548, 0.852501, 0.85595, 0.827991};
	double pionPairMGt_TPC[5] = {0.829238, 0.83258, 0.834903, 0.800073, 0.765804};
	double pionPairMLt_TPC[5] = {0.83052, 0.834775, 0.83705, 0.802637, 0.768361};
	double pionPairEta_TPC[9] = {0.862276, 0.826708, 0.797671, 0.77381, 0.771638, 0.790094, 0.808581, 0.831038, 0.856093};

	// Final TPCorTOF result
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

	double TPC_TPCandTOF_diff_intpT_Gt_p0_final = 0.000730;
	double TPC_TPCandTOF_diff_intpT_Lt_p0_final = 0.000066;
	double TPC_TPCandTOF_diff_intMinv_Gt_p0_final = 0.000437;
	double TPC_TPCandTOF_diff_intMinv_Lt_p0_final = 0.000195;
	double TPC_TPCandTOF_diff_Eta_final_p0_final = 0.000193;

	for (int kk = 0; kk < 9; kk++)
	{
		pionPairEta_TPCorToF[kk] = 1 + fabs(TPC_TPCandTOF_rel_Eta_diff_final);
	}
	// Final TOF result
	double pionPairPtGt_ToF[5] = {0.944208, 0.917944, 0.912098, 0.899861, 0.859474};
	double pionPairPtLt_ToF[5] = {0.948695, 0.919545, 0.912388, 0.900562, 0.862505};
	double pionPairMGt_ToF[5] = {0.884668, 0.88769, 0.88976, 0.871818, 0.84386};
	double pionPairMLt_ToF[5] = {0.884153, 0.888802, 0.890719, 0.873281, 0.846215};
	double pionPairEta_ToF[9] = {0.902961, 0.886886, 0.865055, 0.843727, 0.840466, 0.85691, 0.874384, 0.891725, 0.901808};

	double errXsys[9] = {0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03};
	// double errXsys[9] = {0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17};

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

	double MinV[9] = {0};

	// Trigger Bias Sys
	double WTrigsys_Gt[9] = {0};
	double WTrigsys_Lt[9] = {0};
	// Total sys
	double WTotsys_Gt[9] = {0};
	double WTotsys_Lt[9] = {0};

	// ofstream AUT_Stat_PIDSys;
	ofstream AUT_Stat_LatexTable;
	ofstream AUT_Stat_LatexTable_Lt;

	// const char *outtxtfile = "AUT_Vs_Minv_Stats_PIDSys_trigBias.txt";
	//
	// if (remove(outtxtfile) != 0)
	//{
	//	cout << "Unable to delete file" << endl;
	//}
	// else
	//{
	//	cout << outtxtfile << "\t file deleted\t" << endl;
	//}
	// AUT_Stat_PIDSys.open(outtxtfile, ios::app);

	const char *outtxtfile_Table = "AUT_Vs_Minv_LatexTable_Gt_polCut_TPCorTOF.tex";
	if (remove(outtxtfile_Table) != 0)
	{
		cout << "Unable to delete file" << endl;
	}
	else
	{
		cout << outtxtfile_Table << "\t file deleted\t" << endl;
	}

	const char *outtxtfile_Table_Lt = "AUT_Vs_Minv_LatexTable_Lt_polCut_TPCorTOF.tex";
	if (remove(outtxtfile_Table_Lt) != 0)
	{
		cout << "Unable to delete file" << endl;
	}
	else
	{
		cout << outtxtfile_Table_Lt << "\t file deleted\t" << endl;
	}
	// AUT_Stat_PIDSys.open(outtxtfile, ios::app);
	AUT_Stat_LatexTable.open(outtxtfile_Table, ios::app);
	AUT_Stat_LatexTable_Lt.open(outtxtfile_Table_Lt, ios::app);
	for (int pad = 0; pad < 5; pad++)
	{
		if (pad == 0)
		{
			for (int mbin = 0; mbin < 9; mbin++)
			{
				AvgM_G[mbin] = MinvAvgGt[pad][mbin];
				AvgM_L[mbin] = MinvAvgLt[pad][mbin];
				M_BG[mbin] = MinvBGt[pad][mbin];
				M_BL[mbin] = MinvBLt[pad][mbin];
				M_YG[mbin] = MinvYGt[pad][mbin];
				M_YL[mbin] = MinvYLt[pad][mbin];
				MinV[mbin] = Minv1[mbin];

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

				pidsys_Gt[mbin] = (1 - pionPairPtGt[0]) * TMath::Max(fabs(AvgA_G[mbin]), fabs(dA_G[mbin]));
				pidsys_Lt[mbin] = (1 - pionPairPtLt[0]) * TMath::Max(fabs(AvgA_L[mbin]), fabs(dA_L[mbin]));
				Wpidsys_Gt[mbin] = (1 - pionPairPtGt[0]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin]));
				Wpidsys_Lt[mbin] = (1 - pionPairPtLt[0]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin]));
				// TPC
				WpidsysTPC_Gt[mbin] = (1 - pionPairPtGt_TPC[0]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin]));
				WpidsysTPC_Lt[mbin] = (1 - pionPairPtLt_TPC[0]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin]));
				// TPCorTOF
				WpidsysTPCorTOF_Gt[mbin] = (1 - pionPairPtGt_TPCorToF[0]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin]));
				WpidsysTPCorTOF_Lt[mbin] = (1 - pionPairPtLt_TPCorToF[0]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin]));
				// TOF
				WpidsysTOF_Gt[mbin] = (1 - pionPairPtGt_ToF[0]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin]));
				WpidsysTOF_Lt[mbin] = (1 - pionPairPtLt_ToF[0]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin]));

				// Trigger Bias
				WTrigsys_Gt[mbin] = (1 - TrigBias_MinvGt[0]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin]));
				WTrigsys_Lt[mbin] = (1 - TrigBias_MinvLt[0]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin]));

				// Total sys
				// WTotsys_Gt[mbin] = sqrt(pow(WpidsysTPCorTOF_Gt[mbin], 2) + pow(WTrigsys_Gt[mbin], 2));
				// WTotsys_Lt[mbin] = sqrt(pow(WpidsysTPCorTOF_Lt[mbin], 2) + pow(WTrigsys_Lt[mbin], 2));
				WTotsys_Gt[mbin] = sqrt(pow(TPC_TPCandTOF_diff_intpT_Gt_p0_final, 2) + pow(WTrigsys_Gt[mbin], 2));
				WTotsys_Lt[mbin] = sqrt(pow(TPC_TPCandTOF_diff_intpT_Lt_p0_final, 2) + pow(WTrigsys_Lt[mbin], 2));
			}
		}

		if (pad == 1)
		{
			for (int mbin = 0; mbin < 9; mbin++)
			{
				AvgM_G[mbin] = MinvAvgGt[pad][mbin];
				AvgM_L[mbin] = MinvAvgLt[pad][mbin];
				M_BG[mbin] = MinvBGt[pad][mbin];
				M_BL[mbin] = MinvBLt[pad][mbin];
				M_YG[mbin] = MinvYGt[pad][mbin];
				M_YL[mbin] = MinvYLt[pad][mbin];
				MinV[mbin] = Minv2[mbin];

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

				pidsys_Gt[mbin] = (1 - pionPairPtGt[1]) * TMath::Max(fabs(AvgA_G[mbin]), fabs(dA_G[mbin]));
				pidsys_Lt[mbin] = (1 - pionPairPtLt[1]) * TMath::Max(fabs(AvgA_L[mbin]), fabs(dA_L[mbin]));
				Wpidsys_Gt[mbin] = (1 - pionPairPtGt[1]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin]));
				Wpidsys_Lt[mbin] = (1 - pionPairPtLt[1]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin]));

				// TPC
				WpidsysTPC_Gt[mbin] = (1 - pionPairPtGt_TPC[1]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin]));
				WpidsysTPC_Lt[mbin] = (1 - pionPairPtLt_TPC[1]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin]));
				// TPCorTOF
				WpidsysTPCorTOF_Gt[mbin] = (1 - pionPairPtGt_TPCorToF[1]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin]));
				WpidsysTPCorTOF_Lt[mbin] = (1 - pionPairPtLt_TPCorToF[1]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin]));
				// TOF
				WpidsysTOF_Gt[mbin] = (1 - pionPairPtGt_ToF[1]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin]));
				WpidsysTOF_Lt[mbin] = (1 - pionPairPtLt_ToF[1]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin]));
				// Trigger Bias
				WTrigsys_Gt[mbin] = (1 - TrigBias_MinvGt[1]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin]));
				WTrigsys_Lt[mbin] = (1 - TrigBias_MinvLt[1]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin]));

				// Total sys
				// WTotsys_Lt[mbin] = sqrt(pow(WpidsysTPCorTOF_Lt[mbin], 2) + pow(WTrigsys_Lt[mbin], 2));
				// WTotsys_Gt[mbin] = sqrt(pow(WpidsysTPCorTOF_Gt[mbin], 2) + pow(WTrigsys_Gt[mbin], 2));
				WTotsys_Gt[mbin] = sqrt(pow(TPC_TPCandTOF_diff_intpT_Gt_p0_final, 2) + pow(WTrigsys_Gt[mbin], 2));
				WTotsys_Lt[mbin] = sqrt(pow(TPC_TPCandTOF_diff_intpT_Lt_p0_final, 2) + pow(WTrigsys_Lt[mbin], 2));
			}
		}
		if (pad == 2)
		{
			for (int mbin = 0; mbin < 9; mbin++)
			{
				AvgM_G[mbin] = MinvAvgGt[pad][mbin];
				AvgM_L[mbin] = MinvAvgLt[pad][mbin];
				M_BG[mbin] = MinvBGt[pad][mbin];
				M_BL[mbin] = MinvBLt[pad][mbin];
				M_YG[mbin] = MinvYGt[pad][mbin];
				M_YL[mbin] = MinvYLt[pad][mbin];
				MinV[mbin] = Minv3[mbin];

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

				pidsys_Gt[mbin] = (1 - pionPairPtGt[2]) * TMath::Max(fabs(AvgA_G[mbin]), fabs(dA_G[mbin]));
				pidsys_Lt[mbin] = (1 - pionPairPtLt[2]) * TMath::Max(fabs(AvgA_L[mbin]), fabs(dA_L[mbin]));
				Wpidsys_Gt[mbin] = (1 - pionPairPtGt[2]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin]));
				Wpidsys_Lt[mbin] = (1 - pionPairPtLt[2]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin]));

				// TPC
				WpidsysTPC_Gt[mbin] = (1 - pionPairPtGt_TPC[2]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin]));
				WpidsysTPC_Lt[mbin] = (1 - pionPairPtLt_TPC[2]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin]));
				// TPCorTOF
				WpidsysTPCorTOF_Gt[mbin] = (1 - pionPairPtGt_TPCorToF[2]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin]));
				WpidsysTPCorTOF_Lt[mbin] = (1 - pionPairPtLt_TPCorToF[2]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin]));
				// TOF
				WpidsysTOF_Gt[mbin] = (1 - pionPairPtGt_ToF[2]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin]));
				WpidsysTOF_Lt[mbin] = (1 - pionPairPtLt_ToF[2]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin]));
				// Trigger Bias
				WTrigsys_Gt[mbin] = (1 - TrigBias_MinvGt[2]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin]));
				WTrigsys_Lt[mbin] = (1 - TrigBias_MinvLt[2]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin]));

				// Total sys
				// WTotsys_Gt[mbin] = sqrt(pow(WpidsysTPCorTOF_Gt[mbin], 2) + pow(WTrigsys_Gt[mbin], 2));
				// WTotsys_Lt[mbin] = sqrt(pow(WpidsysTPCorTOF_Lt[mbin], 2) + pow(WTrigsys_Lt[mbin], 2));
				WTotsys_Gt[mbin] = sqrt(pow(TPC_TPCandTOF_diff_intpT_Gt_p0_final, 2) + pow(WTrigsys_Gt[mbin], 2));
				WTotsys_Lt[mbin] = sqrt(pow(TPC_TPCandTOF_diff_intpT_Lt_p0_final, 2) + pow(WTrigsys_Lt[mbin], 2));
			}
		}
		if (pad == 3)
		{
			for (int mbin = 0; mbin < 9; mbin++)
			{
				AvgM_G[mbin] = MinvAvgGt[pad][mbin];
				AvgM_L[mbin] = MinvAvgLt[pad][mbin];
				M_BG[mbin] = MinvBGt[pad][mbin];
				M_BL[mbin] = MinvBLt[pad][mbin];
				M_YG[mbin] = MinvYGt[pad][mbin];
				M_YL[mbin] = MinvYLt[pad][mbin];
				MinV[mbin] = Minv4[mbin];

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

				pidsys_Gt[mbin] = (1 - pionPairPtGt[3]) * TMath::Max(fabs(AvgA_G[mbin]), fabs(dA_G[mbin]));
				pidsys_Lt[mbin] = (1 - pionPairPtLt[3]) * TMath::Max(fabs(AvgA_L[mbin]), fabs(dA_L[mbin]));
				Wpidsys_Gt[mbin] = (1 - pionPairPtGt[3]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin]));
				Wpidsys_Lt[mbin] = (1 - pionPairPtLt[3]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin]));

				// TPC
				WpidsysTPC_Gt[mbin] = (1 - pionPairPtGt_TPC[3]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin]));
				WpidsysTPC_Lt[mbin] = (1 - pionPairPtLt_TPC[3]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin]));
				// TPCorTOF
				WpidsysTPCorTOF_Gt[mbin] = (1 - pionPairPtGt_TPCorToF[3]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin]));
				WpidsysTPCorTOF_Lt[mbin] = (1 - pionPairPtLt_TPCorToF[3]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin]));
				// TOF
				WpidsysTOF_Gt[mbin] = (1 - pionPairPtGt_ToF[3]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin]));
				WpidsysTOF_Lt[mbin] = (1 - pionPairPtLt_ToF[3]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin]));
				// Trigger Bias
				// WTrigsys_Gt[mbin] = 0.2 * TMath::Max(WAvgA_G[mbin], WdA_G[mbin]);
				// WTrigsys_Lt[mbin] = 0.2 * TMath::Max(WAvgA_L[mbin], WdA_L[mbin]);
				WTrigsys_Gt[mbin] = (1 - TrigBias_MinvGt[3]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin]));
				WTrigsys_Lt[mbin] = (1 - TrigBias_MinvLt[3]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin]));

				// Total sys
				// WTotsys_Gt[mbin] = sqrt(pow(WpidsysTPCorTOF_Gt[mbin], 2) + pow(WTrigsys_Gt[mbin], 2));
				// WTotsys_Lt[mbin] = sqrt(pow(WpidsysTPCorTOF_Lt[mbin], 2) + pow(WTrigsys_Lt[mbin], 2));
				WTotsys_Gt[mbin] = sqrt(pow(TPC_TPCandTOF_diff_intpT_Gt_p0_final, 2) + pow(WTrigsys_Gt[mbin], 2));
				WTotsys_Lt[mbin] = sqrt(pow(TPC_TPCandTOF_diff_intpT_Lt_p0_final, 2) + pow(WTrigsys_Lt[mbin], 2));
			}
		}
		if (pad == 4)
		{
			for (int mbin = 0; mbin < 9; mbin++)
			{
				AvgM_G[mbin] = MinvAvgGt[pad][mbin];
				AvgM_L[mbin] = MinvAvgLt[pad][mbin];
				M_BG[mbin] = MinvBGt[pad][mbin];
				M_BL[mbin] = MinvBLt[pad][mbin];
				M_YG[mbin] = MinvYGt[pad][mbin];
				M_YL[mbin] = MinvYLt[pad][mbin];
				MinV[mbin] = Minv5[mbin];

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

				pidsys_Gt[mbin] = (1 - pionPairPtGt[4]) * TMath::Max(fabs(AvgA_G[mbin]), fabs(dA_G[mbin]));
				pidsys_Lt[mbin] = (1 - pionPairPtLt[4]) * TMath::Max(fabs(AvgA_L[mbin]), fabs(dA_L[mbin]));
				Wpidsys_Gt[mbin] = (1 - pionPairPtGt[4]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin]));
				Wpidsys_Lt[mbin] = (1 - pionPairPtLt[4]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin]));

				// TPC
				WpidsysTPC_Gt[mbin] = (1 - pionPairPtGt_TPC[4]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin]));
				WpidsysTPC_Lt[mbin] = (1 - pionPairPtLt_TPC[4]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin]));
				// TPCorTOF
				WpidsysTPCorTOF_Gt[mbin] = (1 - pionPairPtGt_TPCorToF[4]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin]));
				WpidsysTPCorTOF_Lt[mbin] = (1 - pionPairPtLt_TPCorToF[4]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin]));
				// TOF
				WpidsysTOF_Gt[mbin] = (1 - pionPairPtGt_ToF[4]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin]));
				WpidsysTOF_Lt[mbin] = (1 - pionPairPtLt_ToF[4]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin]));

				// Trigger Bias
				// WTrigsys_Gt[mbin] = 0.2 * TMath::Max(WAvgA_G[mbin], WdA_G[mbin]);
				// WTrigsys_Lt[mbin] = 0.2 * TMath::Max(WAvgA_L[mbin], WdA_L[mbin]);
				WTrigsys_Gt[mbin] = (1 - TrigBias_MinvGt[4]) * TMath::Max(fabs(WAvgA_G[mbin]), fabs(WdA_G[mbin]));
				WTrigsys_Lt[mbin] = (1 - TrigBias_MinvLt[4]) * TMath::Max(fabs(WAvgA_L[mbin]), fabs(WdA_L[mbin]));

				// Total sys
				// WTotsys_Gt[mbin] = sqrt(pow(WpidsysTPCorTOF_Gt[mbin], 2) + pow(WTrigsys_Gt[mbin], 2));
				// WTotsys_Lt[mbin] = sqrt(pow(WpidsysTPCorTOF_Lt[mbin], 2) + pow(WTrigsys_Lt[mbin], 2));

				WTotsys_Gt[mbin] = sqrt(pow(TPC_TPCandTOF_diff_intpT_Gt_p0_final, 2) + pow(WTrigsys_Gt[mbin], 2));
				WTotsys_Lt[mbin] = sqrt(pow(TPC_TPCandTOF_diff_intpT_Lt_p0_final, 2) + pow(WTrigsys_Lt[mbin], 2));
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
						MINV[BYA][GL].push_back(M_BG[mbin]);
						A_UT[BYA][GL].push_back(A_BG[mbin]);
						dA_UT[BYA][GL].push_back(dA_BG[mbin]);
					}
					if (BYA == 0 && GL == 1)
					{
						MINV[BYA][GL].push_back(M_BL[mbin]);
						A_UT[BYA][GL].push_back(A_BL[mbin]);
						dA_UT[BYA][GL].push_back(dA_BL[mbin]);
					}
					if (BYA == 1 && GL == 0)
					{
						MINV[BYA][GL].push_back(M_YG[mbin]);
						A_UT[BYA][GL].push_back(A_YG[mbin]);
						dA_UT[BYA][GL].push_back(dA_YG[mbin]);
					}
					if (BYA == 1 && GL == 1)
					{
						MINV[BYA][GL].push_back(M_YL[mbin]);
						A_UT[BYA][GL].push_back(A_YL[mbin]);
						dA_UT[BYA][GL].push_back(dA_YL[mbin]);
					}
					if (BYA == 2 && GL == 0)
					{
						MINV[BYA][GL].push_back(AvgM_G[mbin]);
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
						MINV[BYA][GL].push_back(AvgM_L[mbin]);
						A_UT[BYA][GL].push_back(AvgA_L[mbin]);
						dA_UT[BYA][GL].push_back(dA_L[mbin]);
						PIDSYS[BYA][GL].push_back(pidsys_Lt[mbin]);
						WA_UT[BYA][GL].push_back(WAvgA_L[mbin]);
						WdA_UT[BYA][GL].push_back(WdA_L[mbin]);
						WPIDSYS[BYA][GL].push_back(Wpidsys_Lt[mbin]);
						// TPC
						WPIDSYS_TPC[BYA][GL].push_back(WpidsysTPC_Lt[mbin]);
						WPIDSYS_TPCorTOF[BYA][GL].push_back(WpidsysTPCorTOF_Lt[mbin]);
						WPIDSYS_TOF[BYA][GL].push_back(WpidsysTOF_Lt[mbin]);
						// Total Systematic
						WTOTSYS[BYA][GL].push_back(WTotsys_Lt[mbin]);
					}

					// if(BYA==0&&GL==0){MINV[BYA][GL].push_back(MinV[mbin]); A_UT[BYA][GL].push_back(A_BG[mbin]); dA_UT[BYA][GL].push_back(dA_BG[mbin]);}
					// if(BYA==0&&GL==1){MINV[BYA][GL].push_back(MinV[mbin]); A_UT[BYA][GL].push_back(A_BL[mbin]); dA_UT[BYA][GL].push_back(dA_BL[mbin]);}
					// if(BYA==1&&GL==0){MINV[BYA][GL].push_back(MinV[mbin]); A_UT[BYA][GL].push_back(A_YG[mbin]); dA_UT[BYA][GL].push_back(dA_YG[mbin]);}
					// if(BYA==1&&GL==1){MINV[BYA][GL].push_back(MinV[mbin]); A_UT[BYA][GL].push_back(A_YL[mbin]); dA_UT[BYA][GL].push_back(dA_YL[mbin]);}
					// if(BYA==2&&GL==0){MINV[BYA][GL].push_back(MinV[mbin]); A_UT[BYA][GL].push_back(AvgA_G[mbin]); dA_UT[BYA][GL].push_back(dA_G[mbin]);}
					// if(BYA==2&&GL==1){MINV[BYA][GL].push_back(MinV[mbin]); A_UT[BYA][GL].push_back(AvgA_L[mbin]); dA_UT[BYA][GL].push_back(dA_L[mbin]);}
				} // mbin loop ends her//Here we have a problem instead of pT1 I need to store M1 and so on....

				if (BYA == 2 && GL == 0)
				{

					if (pad == 0)
					{
						AUT_Stat_LatexTable << "\\newgeometry{top=1cm,bottom=0.5cm}" << endl;
						AUT_Stat_LatexTable << "\\section*{Table for $A_{UT}$ vs $M^{\\pi ^ +\\pi ^ -}_{inv}$ for $\\eta^{\\pi ^ +\\pi ^ -}>0$}" << endl;
						AUT_Stat_LatexTable << " \\begin{tabular}{| c | c | c | c | c | c |}" << endl;
						AUT_Stat_LatexTable << "\\hline" << endl;
						AUT_Stat_LatexTable << "\\textbf{$M^{\\pi ^ +\\pi ^ -}_{inv}}$} &\\textbf{$A_{UT}$} &\\textbf{$\\sigma_{stat}$}&\\textbf{$\\sigma_{PID}$}  &\\textbf{$\\sigma_{Trig}$}&\\textbf{$\\sigma_{Tot}$}"
											<< "\\"
											<< "\\" << endl;
						AUT_Stat_LatexTable << "\\hline" << endl;
					}
					if (pad == 1 || pad == 2 || pad == 3 || pad == 4)
					{
						AUT_Stat_LatexTable
							<< "\\hdashline" << endl;
					}
					for (int nBin = 0; nBin < 9; nBin++)
					{
						// AUT_Stat_LatexTable << std::setprecision(2) << AvgM_G[nBin] << "&" << std::setprecision(4) << WAvgA_G[nBin] << "&" << WdA_G[nBin] << "&" << fabs(WpidsysTPCorTOF_Gt[nBin]) << "&" << fabs(WTrigsys_Gt[nBin]) << "\\"
						//<< "\\" << endl;

						AUT_Stat_LatexTable << std::setprecision(2) << AvgM_G[nBin] << "&" << std::setprecision(4) << WAvgA_G[nBin] << "&" << WdA_G[nBin] << "&" << fabs(TPC_TPCandTOF_diff_intpT_Gt_p0_final) << "&" << fabs(WTrigsys_Gt[nBin]) << "&" << WTotsys_Gt[nBin] << "\\"
											<< "\\" << endl;
					}
					if (pad == 4)
					{
						AUT_Stat_LatexTable
							<< "\\hline" << endl;
						AUT_Stat_LatexTable << " \\end{tabular}" << endl;
					}

					// AUT_Stat_PIDSys << "=================================== Weighted Average Forward Dir =============================" << endl;
					// AUT_Stat_PIDSys << "             " << endl;
					// AUT_Stat_PIDSys << "pT-" << pad << "<WA_UT> = "
					//				<< "{" << WAvgA_G[0] << "," << WAvgA_G[1] << "," << WAvgA_G[2] << "," << WAvgA_G[3] << "," << WAvgA_G[4] << "," << WAvgA_G[5] << "," << WAvgA_G[6] << "," << WAvgA_G[7] << "," << WAvgA_G[8] << "}" << endl;
					// AUT_Stat_PIDSys << "pT-" << pad << "<WdAErr> = "
					//				<< "{" << WdA_G[0] << "," << WdA_G[1] << "," << WdA_G[2] << "," << WdA_G[3] << "," << WdA_G[4] << "," << WdA_G[5] << "," << WdA_G[6] << "," << WdA_G[7] << "," << WdA_G[8] << "}" << endl;
					// AUT_Stat_PIDSys << "pT-" << pad << "<WPidSys_TPCorTOF> = "
					//				<< "{" << WpidsysTPCorTOF_Gt[0] << "," << WpidsysTPCorTOF_Gt[1] << "," << WpidsysTPCorTOF_Gt[2] << "," << WpidsysTPCorTOF_Gt[3] << "," << WpidsysTPCorTOF_Gt[4] << "," << WpidsysTPCorTOF_Gt[5] << "," << WpidsysTPCorTOF_Gt[6] << "," << WpidsysTPCorTOF_Gt[7] << "," << WpidsysTPCorTOF_Gt[8] << "}" << endl;
					// AUT_Stat_PIDSys << "pT-" << pad << "<WTOTSys_TPCorTOF> = "
					//				<< "{" << WTotsys_Gt[0] << "," << WTotsys_Gt[1] << "," << WTotsys_Gt[2] << "," << WTotsys_Gt[3] << "," << WTotsys_Gt[4] << "," << WTotsys_Gt[5] << "," << WTotsys_Gt[6] << "," << WTotsys_Gt[7] << "," << WTotsys_Gt[8] << "}" << endl;
					// AUT_Stat_PIDSys << "             " << endl;
					// AUT_Stat_PIDSys << "             " << endl;

				} // Average and Forward direction
				if (BYA == 2 && GL == 1)
				{

					if (pad == 0)
					{
						AUT_Stat_LatexTable_Lt << "\\newpage" << endl;
						AUT_Stat_LatexTable_Lt << "\\newgeometry{top=1cm,bottom=0.5cm}" << endl;
						AUT_Stat_LatexTable_Lt << "\\section*{Table for $A_{UT}$ vs $M^{\\pi ^ +\\pi ^ -}_{inv}$ for $\\eta^{\\pi ^ +\\pi ^ -}<0$}" << endl;
						AUT_Stat_LatexTable_Lt << " \\begin{tabular}{ |c | c | c | c | c|c|}" << endl;
						AUT_Stat_LatexTable_Lt << "\\hline" << endl;
						AUT_Stat_LatexTable_Lt << "\\textbf{$M^{\\pi ^ +\\pi ^ -}_{inv}$}&\\textbf{$A_{UT}$} &\\textbf{$\\sigma_{stat}$}&\\textbf{$\\sigma_{PID}$}  &\\textbf{$\\sigma_{Trig}$}&\\textbf{$\\sigma_{Tot}$}"
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
						AUT_Stat_LatexTable_Lt << std::setprecision(2) << AvgM_L[nBin] << "&" << std::setprecision(4) << WAvgA_L[nBin] << "&" << WdA_L[nBin] << "&" << fabs(TPC_TPCandTOF_diff_intpT_Lt_p0_final) << "&" << fabs(WTrigsys_Lt[nBin]) << "&" << WTotsys_Lt[nBin]
											   << "\\"
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
					//				<< "{" << WAvgA_L[0] << "," << WAvgA_L[1] << "," << WAvgA_L[2] << "," << WAvgA_L[3] << "," << WAvgA_L[4] << "," << WAvgA_L[5] << "," << WAvgA_L[6] << "," << WAvgA_L[7] << "," << WAvgA_L[8] << "}" << endl;
					// AUT_Stat_PIDSys << "pT-" << pad << "<WdAErr> = "
					//				<< "{" << WdA_L[0] << "," << WdA_L[1] << "," << WdA_L[2] << "," << WdA_L[3] << "," << WdA_L[4] << "," << WdA_L[5] << "," << WdA_L[6] << "," << WdA_L[7] << "," << WdA_L[8] << "}" << endl;
					// AUT_Stat_PIDSys << "pT-" << pad << "<WPidSys_TPCorTOF> = "
					//				<< "{" << WpidsysTPCorTOF_Lt[0] << "," << WpidsysTPCorTOF_Lt[1] << "," << WpidsysTPCorTOF_Lt[2] << "," << WpidsysTPCorTOF_Lt[3] << "," << WpidsysTPCorTOF_Lt[4] << "," << WpidsysTPCorTOF_Lt[5]
					//				<< "," << WpidsysTPCorTOF_Lt[6] << "," << WpidsysTPCorTOF_Lt[7] << "," << WpidsysTPCorTOF_Lt[8] << "}" << endl;
					// AUT_Stat_PIDSys << "pT-" << pad << "<WTOTSys_TPCorTOF> = "
					//				<< "{" << Wpidsys_Lt[0] << "," << Wpidsys_Lt[1] << "," << Wpidsys_Lt[2] << "," << Wpidsys_Lt[3] << "," << Wpidsys_Lt[4] << "," << Wpidsys_Lt[5]
					//				<< "," << Wpidsys_Lt[6] << "," << Wpidsys_Lt[7] << "," << Wpidsys_Lt[8] << "}" << endl;
					//
					// AUT_Stat_PIDSys << "             " << endl;
					// AUT_Stat_PIDSys << "             " << endl;
				} // Average and Backward dir
			}	  // GL loop ends here

		} // BYA loops ends here

		for (int BYA = 0; BYA < 3; BYA++)
		{
			for (int GL = 0; GL < 2; GL++)
			{
				// GL==0 forward GL==1 backward
				//    cout << PT[BYA][GL].size() << " =pT size for "<< BYA << "\t BYA\t and  "<< GL << "\t GL\t "<< endl;
				//    cout << A_UT[BYA][GL].size() << " =A_UT size for "<< BYA << "\t BYA\t and  "<< GL << "\t GL\t "<< endl;
				//    cout << dA_UT[BYA][GL].size()<< " =dA_UT size for "<< BYA << "\t BYA\t and  "<< GL << "\t GL\t "<< endl;
				GRAPH[BYA][GL][pad] = new TGraphErrors(MINV[BYA][GL].size(), &MINV[BYA][GL][0], &A_UT[BYA][GL][0], 0, &dA_UT[BYA][GL][0]);
				if (BYA == 2)
				{
					gr_sys[BYA][GL][pad] = new TGraphErrors(MINV[BYA][GL].size(), &MINV[BYA][GL][0], &A_UT[BYA][GL][0], errXsys, &PIDSYS[BYA][GL][0]);

					Wgr_sys[BYA][GL][pad] = new TGraphErrors(MINV[BYA][GL].size(), &MINV[BYA][GL][0], &WA_UT[BYA][GL][0], errXsys, &WPIDSYS[BYA][GL][0]);

					Wgr_sys_TPC[BYA][GL][pad] = new TGraphErrors(MINV[BYA][GL].size(), &MINV[BYA][GL][0], &WA_UT[BYA][GL][0], errXsys, &WPIDSYS_TPC[BYA][GL][0]);
					Wgr_sys_TPCorTOF[BYA][GL][pad] = new TGraphErrors(MINV[BYA][GL].size(), &MINV[BYA][GL][0], &WA_UT[BYA][GL][0], errXsys, &WPIDSYS_TPCorTOF[BYA][GL][0]);
					Wgr_sys_TOF[BYA][GL][pad] = new TGraphErrors(MINV[BYA][GL].size(), &MINV[BYA][GL][0], &WA_UT[BYA][GL][0], errXsys, &WPIDSYS_TOF[BYA][GL][0]);
					// Only not to draw systematic
					Wgr_TOTSYS[BYA][GL][pad] = new TGraphErrors(MINV[BYA][GL].size(), &MINV[BYA][GL][0], &WA_UT[BYA][GL][0], errXsys, &WTOTSYS[BYA][GL][0]);
					// Wgr_sys[BYA][GL][pad] = new TGraphErrors(MINV[BYA][GL].size(), &MINV[BYA][GL][0], &WA_UT[BYA][GL][0], 0, 0);
					WGRAPH[BYA][GL][pad] = new TGraphErrors(MINV[BYA][GL].size(), &MINV[BYA][GL][0], &WA_UT[BYA][GL][0], 0, &WdA_UT[BYA][GL][0]);
				}
			}
		}

		for (int BYA = 0; BYA < 3; BYA++)
		{
			for (int GL = 0; GL < 2; GL++)
			{
				MINV[BYA][GL].clear();
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
				GRAPH[BYA][GL][pad]->SetTitle("");
				GRAPH[BYA][GL][pad]->GetYaxis()->SetTitleOffset(1.4);
				GRAPH[BYA][GL][pad]->GetYaxis()->SetLabelOffset(0.017);
				GRAPH[BYA][GL][pad]->GetYaxis()->SetTitleSize(0.06);
				GRAPH[BYA][GL][pad]->GetYaxis()->SetLabelSize(0.06);
				GRAPH[BYA][GL][pad]->GetYaxis()->SetLabelFont(22);
				GRAPH[BYA][GL][pad]->GetXaxis()->SetLabelFont(22);
				GRAPH[BYA][GL][pad]->GetXaxis()->SetTitle("#font[22]{M_{inv}(GeV/c^{2})}");
				GRAPH[BYA][GL][pad]->GetXaxis()->SetTitleSize(0.06);
				GRAPH[BYA][GL][pad]->GetXaxis()->SetLabelSize(0.06);
				GRAPH[BYA][GL][pad]->GetYaxis()->SetRangeUser(-0.0252, 0.058);
				GRAPH[BYA][GL][pad]->GetXaxis()->SetLimits(0.02, 2.48);
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
					// Wgr_sys[BYA][GL][pad]->GetXaxis()->SetTitle("#font[22]{M_{inv}^{#pi^{+}#pi^{-}}(GeV/c^{2})}");
					// Wgr_sys[BYA][GL][pad]->GetXaxis()->SetTitleSize(0.06);
					// Wgr_sys[BYA][GL][pad]->GetXaxis()->SetLabelSize(0.06);
					// Wgr_sys[BYA][GL][pad]->GetYaxis()->SetRangeUser(-0.0252, 0.058);
					// Wgr_sys[BYA][GL][pad]->GetXaxis()->SetLimits(0.02, 2.48);
					// Wgr_sys[BYA][GL][pad]->GetXaxis()->SetNdivisions(505);
					// Wgr_sys[BYA][GL][pad]->GetYaxis()->SetNdivisions(505);
					// Wgr_sys[BYA][GL][pad]->SetLineColor(1);
					// Wgr_sys[BYA][GL][pad]->SetLineWidth(1);
					// Wgr_sys[BYA][GL][pad]->SetFillStyle(0);

					// PID sys
					// Wgr_sys_TPCorTOF[BYA][GL][pad]->GetYaxis()->SetTitle("#font[22]{A_{UT}^{Sin(#Phi_{S}-#Phi_{R})}}");
					// Wgr_sys_TPCorTOF[BYA][GL][pad]->SetTitle("");
					// Wgr_sys_TPCorTOF[BYA][GL][pad]->GetYaxis()->SetTitleOffset(1.4);
					// Wgr_sys_TPCorTOF[BYA][GL][pad]->GetYaxis()->SetLabelOffset(0.017);
					// Wgr_sys_TPCorTOF[BYA][GL][pad]->GetYaxis()->SetTitleSize(0.06);
					// Wgr_sys_TPCorTOF[BYA][GL][pad]->GetYaxis()->SetLabelSize(0.06);
					// Wgr_sys_TPCorTOF[BYA][GL][pad]->GetYaxis()->SetLabelFont(22);
					// Wgr_sys_TPCorTOF[BYA][GL][pad]->GetXaxis()->SetLabelFont(22);
					// Wgr_sys_TPCorTOF[BYA][GL][pad]->GetXaxis()->SetTitle("#font[22]{M_{inv}^{#pi^{+}#pi^{-}}(GeV/c^{2})}");
					// Wgr_sys_TPCorTOF[BYA][GL][pad]->GetXaxis()->SetTitleSize(0.06);
					// Wgr_sys_TPCorTOF[BYA][GL][pad]->GetXaxis()->SetLabelSize(0.06);
					// Wgr_sys_TPCorTOF[BYA][GL][pad]->GetYaxis()->SetRangeUser(-0.0252, 0.058);
					// Wgr_sys_TPCorTOF[BYA][GL][pad]->GetXaxis()->SetLimits(0.02, 2.48);
					// Wgr_sys_TPCorTOF[BYA][GL][pad]->GetXaxis()->SetNdivisions(505);
					// Wgr_sys_TPCorTOF[BYA][GL][pad]->GetYaxis()->SetNdivisions(505);
					// Wgr_sys_TPCorTOF[BYA][GL][pad]->SetLineColor(1);
					// Wgr_sys_TPCorTOF[BYA][GL][pad]->SetLineWidth(1);
					// Wgr_sys_TPCorTOF[BYA][GL][pad]->SetFillStyle(0);

					// Total Systematic
					Wgr_TOTSYS[BYA][GL][pad]->GetYaxis()->SetTitle("#font[22]{A_{UT}^{Sin(#Phi_{S}-#Phi_{R})}}");
					Wgr_TOTSYS[BYA][GL][pad]->GetYaxis()->CenterTitle(kTRUE);
					Wgr_TOTSYS[BYA][GL][pad]->SetTitle("");
					Wgr_TOTSYS[BYA][GL][pad]->GetYaxis()->SetTitleOffset(1.4);
					Wgr_TOTSYS[BYA][GL][pad]->GetYaxis()->SetLabelOffset(0.017);
					Wgr_TOTSYS[BYA][GL][pad]->GetYaxis()->SetTitleSize(0.06);
					Wgr_TOTSYS[BYA][GL][pad]->GetYaxis()->SetLabelSize(0.06);
					Wgr_TOTSYS[BYA][GL][pad]->GetYaxis()->SetLabelFont(22);
					Wgr_TOTSYS[BYA][GL][pad]->GetXaxis()->SetLabelFont(22);
					Wgr_TOTSYS[BYA][GL][pad]->GetXaxis()->SetTitle("#font[22]{M_{inv}^{#pi^{+}#pi^{-}}(GeV/c^{2})}");
					Wgr_TOTSYS[BYA][GL][pad]->GetXaxis()->SetTitleSize(0.06);
					Wgr_TOTSYS[BYA][GL][pad]->GetXaxis()->SetLabelSize(0.06);
					// Wgr_TOTSYS[BYA][GL][pad]->GetYaxis()->SetRangeUser(-0.0252, 0.058);
					Wgr_TOTSYS[BYA][GL][pad]->GetYaxis()->SetRangeUser(-0.0122, 0.0452);
					Wgr_TOTSYS[BYA][GL][pad]->GetXaxis()->SetLimits(0.02, 2.48);
					Wgr_TOTSYS[BYA][GL][pad]->GetXaxis()->SetNdivisions(505);
					Wgr_TOTSYS[BYA][GL][pad]->GetYaxis()->SetNdivisions(505);
					Wgr_TOTSYS[BYA][GL][pad]->SetLineColor(1);
					Wgr_TOTSYS[BYA][GL][pad]->SetLineWidth(1);
					Wgr_TOTSYS[BYA][GL][pad]->SetFillStyle(0);

					if (GL == 0)
					{
						GRAPH[BYA][GL][pad]->SetMarkerStyle(27);
						GRAPH[BYA][GL][pad]->SetMarkerColor(1);
						GRAPH[BYA][GL][pad]->SetLineColor(1);

						WGRAPH[BYA][GL][pad]->SetMarkerStyle(27);
						WGRAPH[BYA][GL][pad]->SetMarkerColor(1);
						WGRAPH[BYA][GL][pad]->SetLineColor(1);
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

			} // GL loop ends here
		}	  // BYA loop ends here

	} // pad loop ends here
	// AUT_Stat_PIDSys.close();		// Txt file with AUT_StatError_PIDSysError close here.
	AUT_Stat_LatexTable.close();	// Txt file with AUT_Minv_PIDSys_TrigSys
	AUT_Stat_LatexTable_Lt.close(); // Txt file with AUT_Minv_PIDSys_TrigSys
	// pT_boundary line;
	double b_line[5][10];
	for (int i = 0; i < 5; i++)
	{
		if (i == 0)
		{
			for (int j = 0; j < 10; j++)
			{
				b_line[i][j] = M1[j];
			}
		}
		if (i == 1)
		{
			for (int j = 0; j < 10; j++)
			{
				b_line[i][j] = M2[j];
			}
		}
		if (i == 2)
		{
			for (int j = 0; j < 10; j++)
			{
				b_line[i][j] = M3[j];
			}
		}
		if (i == 3)
		{
			for (int j = 0; j < 10; j++)
			{
				b_line[i][j] = M4[j];
			}
		}
		if (i == 4)
		{
			for (int j = 0; j < 10; j++)
			{
				b_line[i][j] = M5[j];
			}
		}
	}
	TLine *line[5];
	TLegend *leg;
	TLatex tex[5];

	TCanvas *myCanvA[3];
	const char *FBA_name[] = {"Forward", "BackWard", "Avg_BL"};
	TLine *bound_line = new TLine();

	// FBA= Forward(Blue and Yellow), BackWard(Blue and Yellow) and Averaged(Average of Blue and Yellow for forward and Backward)
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
					GRAPH[0][0][PAD]->SetMarkerStyle(24);	   // Blue Eta_pair>0;
					GRAPH[0][0][PAD]->SetMarkerColor(4);	   // Blue Eta_pair>0;
					GRAPH[0][0][PAD]->SetLineColor(4);		   // Blue Eta_pair>0;
					GRAPH[1][0][PAD]->SetMarkerStyle(20);	   // Yellow Eta_pair>0;
					GRAPH[1][0][PAD]->SetMarkerColor(kYellow); // Yellow Eta_pair>0;
					GRAPH[1][0][PAD]->SetLineColor(41);		   // Yellow Eta_pair>0;
					GRAPH[0][0][PAD]->Draw("AP");
					GRAPH[1][0][PAD]->Draw("P same");
					for (int kk = 0; kk < 10; kk++)
					{
						bound_line->DrawLine(b_line[PAD][kk], 0.058 - 0.005, b_line[PAD][kk], 0.058);
					}
					if (PAD == 1)
					{
						leg = new TLegend(0.55, 0.66, 0.97, 0.84);
						leg->AddEntry(GRAPH[0][0][PAD], "#font[22]{#eta^{#pi^{+}#pi^{-}}>0,BLUE}", "lp");
						leg->AddEntry(GRAPH[1][0][PAD], "#font[22]{#eta^{#pi^{+}#pi^{-}}>0,YELLOW}", "lp");
						leg->Draw();
					}
				}
				if (FBA == 1)
				{
					// GRAPH[BYA][GL][pad]
					GRAPH[0][1][PAD]->SetMarkerStyle(24);	   // Blue Eta_pair<0;
					GRAPH[0][1][PAD]->SetMarkerColor(4);	   // Blue Eta_pair<0;
					GRAPH[0][1][PAD]->SetLineColor(4);		   // Blue Eta_pair<0;
					GRAPH[1][1][PAD]->SetMarkerStyle(20);	   // Yellow Eta_pair<0;
					GRAPH[1][1][PAD]->SetMarkerColor(kYellow); // Yellow Eta_pair<0;
					GRAPH[1][1][PAD]->SetLineColor(41);		   // Yellow Eta_pair<0;
					GRAPH[0][1][PAD]->Draw("AP");
					GRAPH[1][1][PAD]->Draw("P same");
					for (int kk = 0; kk < 10; kk++)
					{
						bound_line->DrawLine(b_line[PAD][kk], 0.058 - 0.005, b_line[PAD][kk], 0.058);
					}
					if (PAD == 1)
					{
						leg = new TLegend(0.55, 0.66, 0.97, 0.84);
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
					// WGRAPH[2][1][PAD]->Draw("P same");	  // Backward

					// Wgr_sys[2][0][PAD]->SetFillStyle(0);
					// Wgr_sys[2][0][PAD]->SetLineColor(3); // Green Color
					// Wgr_sys[2][0][PAD]->SetLineWidth(2);
					// Wgr_sys[2][0][PAD]->Draw("A2"); // Forward
					// gPad->Update();

					// Wgr_sys_TPC[2][0][PAD]->SetFillStyle(0);
					// Wgr_sys_TPC[2][0][PAD]->SetLineColor(6); // pink color
					// Wgr_sys_TPC[2][0][PAD]->SetLineWidth(2);
					//  Wgr_sys_TPC[2][0][PAD]->Draw("2 same"); // Forward
					// Wgr_sys_TPC[2][0][PAD]->Draw("A2"); // Forward
					// gPad->Update();

					// TPCorTOF PID systematic;
					// Wgr_sys_TPCorTOF[2][0][PAD]->SetFillStyle(0);
					// Wgr_sys_TPCorTOF[2][0][PAD]->SetLineColor(1); // Black color
					// Wgr_sys_TPCorTOF[2][0][PAD]->Draw("A2");	  // Forward
					//
					// Wgr_sys_TPCorTOF[2][1][PAD]->SetFillStyle(0);
					// Wgr_sys_TPCorTOF[2][1][PAD]->SetLineColor(2); // Red color
					//
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
													   // Wgr_sys[2][0][PAD]->Draw("AE2");   // Forward
													   // WGRAPH[2][0][PAD]->Draw("P same"); // Forward

					// GRAPH_R11[PAD]->SetMarkerStyle(3);
					// GRAPH_R11[PAD]->SetMarkerColor(4);
					// GRAPH_R11[PAD]->SetLineColor(4);
					// GRAPH_R11[PAD]->Draw("P same"); // Draw Run 11
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
						// To draw Bin Boundaries
						// bound_line->DrawLine(b_line[PAD][kk], 0.0452 - 0.002, b_line[PAD][kk], 0.0452);
					}
					if (PAD == 1)
					{
						// leg = new TLegend(0.55, 0.56, 1.0, 0.72);
						leg = new TLegend(0.45, 0.75, 1.0, 0.82);
						leg->AddEntry(WGRAPH[2][0][PAD], "#font[22]{#eta^{#pi^{+}#pi^{-}} > 0}", "lep");

						// leg->AddEntry(GRAPH_R11[PAD], "Run11 #font[22]{#eta^{#pi^{+}#pi^{-}}>0}", "lp");
						leg->AddEntry(WGRAPH[2][1][PAD], "#font[22]{#eta^{#pi^{+}#pi^{-}} < 0}", "lep");
						leg->AddEntry(Wgr_TOTSYS[2][0][PAD], "#font[22]{Tot Sys.}", "f");
						leg->SetNColumns(3);

						TLegend *leg_JAMDiFF = new TLegend(0.45, 0.65, 1.0, 0.75);
						leg_JAMDiFF->SetHeader("#font[22]{C. Cocuzza et. al. (JAMDiFF)}");
						// leg->AddEntry((TObject *)0, "C. Cocuzza et. al. (JAMDiFF)", ""); // Placeholder entry for the text
						leg_JAMDiFF->AddEntry(GRAPH_JAMDiFF_G[PAD], "#font[22]{#eta^{#pi^{+}#pi^{-}} > 0}", "f");
						leg_JAMDiFF->AddEntry(GRAPH_JAMDiFF_L[PAD], "#font[22]{#eta^{#pi^{+}#pi^{-}} < 0}", "f");
						leg_JAMDiFF->SetNColumns(2);
						// leg->AddEntry(GRAPH_JAMDiFF_G[PAD], "C. Cocuzza et. al. (JAMDiFF)", "f");
						// leg->AddEntry(GRAPH_JAMDiFF_L[PAD], "C. Cocuzza et. al. (JAMDiFF)", "f");

						// leg->AddEntry("", "#font[22]{cone > 0.7}", "");

						// leg->AddEntry(Wgr_sys_TPC[2][0][PAD], "TPC");
						// leg->AddEntry(Wgr_sys_TPCorTOF[2][0][PAD], "PID sys.");
						// leg->AddEntry(Wgr_sys_TOF[2][0][PAD], "TOF");
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
				tex[PAD].DrawLatex(1.25, 0.022, Form("#color[1]{#font[22]{<p^{#pi^{+}#pi^{-}}_{T}> = %g GeV/c}}", std::round(avg_pT_pair[PAD])));
				gPad->Update();
				if (PAD == 1)
				{
					// tex[PAD].DrawLatex(0.7, 0.048, "M_{inv} bin boundaries");
				}
			} // PAD<5 Control statement
		}	  // PAD Loop
		myCanvA[FBA]->Print(Form("./TriggerBias/AUT_Vs_Minv_%s_Cone0.7.pdf", FBA_name[FBA]));
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
	c_chi->Print("./TriggerBias/Chi_NDF_Minv.pdf");

} // Main function ends here
