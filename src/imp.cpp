#define imp_cxx
#include "imp.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include<sstream>
#include <TSpectrum.h>
#include<vector>
#include<TGraph.h>
#include<TF1.h>
void imp::MakeTH(TFile *File, TTree *tree)
{
	// 查看LaBr数量，Ge数量，2LaBr填入到三维矩阵，2Ge填入到二维矩阵
	tree -> Branch("hit", &nXIA, "hit/I"); 
	tree -> Branch("cfd", Cfd, "cfd[hit]/D"); 
	tree -> Branch("rise", Rise, "rise[hit]/D"); 
	tree -> Branch("OrginalCfd", data_XIA_CFD, "OrginalCfd[hit]/D"); 
	tree -> Branch("evte", Evte, "evte[hit]/D"); 
	tree -> Branch("ch", Ch, "ch[hit]/s"); 
	tree -> Branch("cid", data_XIA_crate_id, "cid[hit]/s"); 
	tree -> Branch("ts", data_XIA_Event_ts, "ts[hit]/L"); 
	tree -> Branch("e", e, "e[hit]/D"); 

	for(Short_t i = 0; i < kMaxdata_XIA; i++)
	{
		hLa[i] = new TH1I(Form("hLa%02d", i), Form("hLa%02d", i), 4000, -100, 100); 
	}
	hHitCh = new TH1I("hHitCh", "hHitCh", LaNum, 0, LaNum); 
	hTimeCh = new TH2I("hTimeCh%02d", "hTimeCh%02d", 2e4, -1000, 1000, LaNum - 1, 1, LaNum); 
	f -> cd(); 
	fChain -> GetEntry(0); 
	Double_t Ts_Begin = data_XIA_Event_ts[0]; 
	Long64_t Nentry = fChain -> GetEntries(); 
	b_data_XIA_Event_ts -> GetEntry(Nentry - 1); 
	Double_t Ts_End = data_XIA_Event_ts[0]; 

	hTimeMeasure = new TH2I("hTimeMeasure01", "hTimeMeasure01", 2000, -50, 50, 2000, Ts_Begin, Ts_End); 

	std::cout << Nentry << std::endl; 
	Short_t k = 0; 
	for(Long64_t Ientry = 0; Ientry < Nentry; Ientry++)
	{
		b_nXIA -> GetEntry(Ientry); 
		if(Ientry * 100. / Nentry >= k)
		{
			std::cout << Form("%.2f", Ientry * 100. / Double_t(Nentry)) << "%" << std::endl; 
			k = k + 10; 
		}
		// b_data_XIA_ -> GetEntry(Ientry); 
		if(nXIA < 2) continue; 
		//b_data_XIA_channel -> GetEntry(Ientry); 
		fChain -> GetEntry(Ientry); 
		Short_t GeNum = 0; 
		for(Short_t i = 0; i < nXIA; i++)
		{
			if(data_XIA_channel[i] > 199)
			{
				GeNum++; 
				Ch[i] = data_XIA_channel[i] / 100 * 16 - 16 + data_XIA_channel[i] % 100; 
				hHitCh -> Fill(Ch[i]); 
				Evte[i] = data_XIA_Energy[i] * kbHPGe[0][Ch[i] - 32] + kbHPGe[1][Ch[i] - 32]; 
				Cfd[i] = 0; 
				Rise[i] = 0; 
				e[i] = 0; 
				continue; 
			}
			Ch[i] = data_XIA_channel[i] / 100 * 16 - 16 + data_XIA_channel[i] % 100; 
			hHitCh -> Fill(Ch[i]); 
			Evte[i] = data_XIA_Energy[i] * data_XIA_Energy[i] * kbLaBr[0][Ch[i]] + data_XIA_Energy[i] * kbLaBr[1][Ch[i]] + kbLaBr[2][Ch[i]]; 
			Filter(i); 
			Double_t WaveMax = FindMax(); 
			Double_t BaseLine = FindBaseLine(); 
			// Rise[i] = CalculateCFD(WaveMax, .9, BaseLine) / 1000 - CalculateCFD(WaveMax, .1, BaseLine) / 1000; 
			Cfd[i] = CalculateCFD(WaveMax, Rate, BaseLine) / 1000; 
			}
		// if(GeNum > 1)	FillGe(Ientry); 
		tree -> Fill(); 
				for(Short_t i = 0; i < nXIA; i++) if(Ch[i] == 0)	//if(Evte[i] < CoPeak[1] + Sigma && Evte[i] > CoPeak[1] - Sigma)
				for(Short_t j = 0; j < nXIA; j++)
				{
					if(Ch[j] != 0)	//if(Evte[j] < CoPeak[0] + Sigma && Evte[j] > CoPeak[0] - Sigma)
					{
						hLa[Ch[j]] -> Fill(Cfd[i] - Cfd[j] + data_XIA_Event_ts[i] - data_XIA_Event_ts[j]);
						hTimeCh -> Fill(data_XIA_Event_ts[i] - data_XIA_Event_ts[j] + Cfd[i] - Cfd[j], Ch[i]); 
						if(Ch[j] == 1)	hTimeMeasure -> Fill(Cfd[i] - Cfd[j] + data_XIA_Event_ts[i] - data_XIA_Event_ts[j], data_XIA_Event_ts[j]); 
					}
				}
		// if(nXIA - GeNum > 1)	FillLa(Ientry); 
	}
	TGraph* g = new TGraph(); 
	TF1 *f = new TF1("f", "gaus"); 
	for(Short_t i = 0; i < kMaxdata_XIA - 1; i++)
	{
		if(hLa[i] -> GetEntries())
		{
			File -> WriteObject(hLa[i], Form("hLa%02d", i)); 
			hLa[i] -> Fit(f, "SQ+", "", hLa[i] -> GetBinCenter(hLa[i] -> GetMaximumBin()) - 0.5, hLa[i] -> GetBinCenter(hLa[i] -> GetMaximumBin()) + 0.5); 
			g -> SetPoint(i, i + 1, f -> GetParameter(1)); 
		}
		else
		{
			g -> SetPoint(i, i + 1, 5000); 
		}
	}
	File -> WriteObject(g, "g"); 
	File -> WriteObject(hTimeMeasure, "hTimeMeasure"); 
	File -> WriteObject(hTimeCh, "hTimeCh"); 
	File -> WriteObject(hHitCh, "hHitCh"); 
	File -> WriteObject(hGe, "HPGe"); 
	File -> WriteObject(tree, "tree"); 
	//File -> WriteObject(hLa, "LaBr"); 
}
void imp::FillGe(Long64_t Ientry)
{
	b_data_XIA_Energy -> GetEntry(Ientry); 
	// std::cout << Ientry << std::endl; 
	for(Short_t i = 0; i < nXIA; i++)	if(data_XIA_channel[i] > 299)	for(Short_t j = 0; j < nXIA; j++)	if(i != j && data_XIA_channel[j] > 299)
	{
		Short_t Indexi = data_XIA_channel[i] % 100; 
		Short_t Indexj = data_XIA_channel[j] % 100; 
		hGe -> Fill(kbHPGe[0][Indexi] * data_XIA_Energy[i] + kbHPGe[1][Indexi], kbHPGe[0][Indexj] * data_XIA_Energy[j] + kbHPGe[1][Indexj]); 
	}
}
void imp::OriginalCFD()
{
	Long64_t N = 0; 
	TFile* File = new TFile("impsys.root", "update"); 
	TH1I* hTime = new TH1I(Form("Time%05d", FileNum), "Time", 100, 450, 480); 

	Long64_t Nentry = fChain->GetEntries(); 
	for(Long64_t Ientry = 0; Ientry < Nentry; Ientry++)
	{
		fChain -> GetEntry(Ientry); 
		if(nXIA < 2)	continue; 
		N++; 
		for(Short_t i = 0; i < nXIA; i++)
		{
			Double_t Energy = data_XIA_Energy[i] * kb[0][(data_XIA_channel[i] - 100) / 2] + kb[1][(data_XIA_channel[i] - 100) / 2]; 
			if(Energy < CoPeak[1] + Sigma && Energy > CoPeak[1] - Sigma)
			{
				for (Short_t j = 0; j < nXIA; j++)
				{
					// if(i == j)	continue; 
					Double_t Energy2 = data_XIA_Energy[j] * kb[0][(data_XIA_channel[j] - 100) / 2] + kb[1][(data_XIA_channel[j] - 100) / 2]; 
					if(Energy2 < CoPeak[0] + Sigma && Energy2 > CoPeak[0] - Sigma)
					{
						hTime -> Fill(data_XIA_Event_ts[i] + data_XIA_CFD[i] - data_XIA_CFD[j] - data_XIA_Event_ts[j]); 
					}
				}
			}
		}
	}
	std::cout << N << std::endl; 
	File -> WriteObject(hTime, "hTime");
	File -> Close();	
}

Double_t imp::TSINC(Short_t m)
{
	if(m == 0)	return 1; 
	return NSINC * TMath::Sin(m * TMath::Pi() / NSINC) * TMath::Gaus(m, 0, 2 * NSINC)/ m /TMath::Pi(); 
}
void imp::Filter(Short_t ihit)
{
	//for(Short_t i = 0; i < MaxTrace - 2; i++)	Wave[i] = data_XIA_trace[ihit][i] + data_XIA_trace[ihit][i + 1] + data_XIA_trace[ihit][i + 2]; 
	for(Short_t i = 0; i < MaxTrace; i++)	Wave[i] = data_XIA_trace[ihit][i]; 
}
Double_t imp::FindMax()
{
	IndexMax = 0; 
	Double_t CFDed[2 * NSINC + 1]; 
	CFDed[2 * NSINC] = Wave[0]; 

	for(Short_t i = LBase; i < MaxTrace - 2; i++)
	{
		if(Wave[i] > Wave[IndexMax])	IndexMax = i; 
	}
	Double_t M = Wave[IndexMax]; 
	for(Short_t i = 0; i < 2; i++)
	{
		for(Short_t j = 0; j < NSINC; j++)
		{
			CFDed[NSINC * i + j] = 0; 
			for(Short_t k = 0; k < LSINC; k++)
			{
				CFDed[NSINC * i + j] += Wave[(i + IndexMax - 1) - k] * TSINC(k * NSINC + j) + Wave[i + IndexMax + k] * TSINC((k + 1) * NSINC - j); 
			}
			if(CFDed[NSINC * i + j] > M)
			{
				M = CFDed[NSINC * i + j]; 
			}
		}
	}
	return M; 
}
Double_t imp::FindBaseLine()
{
	Double_t BaseLine; 
	BaseLine = 0; 
	for(Short_t i = 0; i < LBase; i++)	BaseLine += Wave[i]; 
	BaseLine /= LBase; 
	return BaseLine; 
}
Double_t imp::CalculateCFD(Double_t M, Double_t RiseRate, Double_t BaseLine){
	Double_t	CFD = 0; 
	Double_t CFDed[NSINC]; 

	M = (M - BaseLine) * RiseRate; 
	for(Short_t i = LBase; i < IndexMax + 1; i++)	if(Wave[i] - BaseLine > M)
	{
		IndexMin = i - 1; 
		break; 
	}
	Short_t IndexSINC = NSINC - 1; 
	for(Short_t j = 1; j < NSINC; j++)
	{
		CFDed[j] = 0; 
		for(Short_t k = 0; k < LSINC; k++)
		{
			CFDed[j] += Wave[IndexMin - k] * TSINC(k * NSINC + j) + Wave[IndexMin + k + 1] * TSINC((k + 1) * NSINC - j); 
		}
		if(CFDed[j] - BaseLine > M)
		{
			IndexSINC = j - 1; 
			break; 
		}
	}
	if(IndexSINC == 0)	CFD = IndexMin * N + (Wave[IndexMin] - BaseLine - M) / (Wave[IndexMin] - CFDed[1]) * (N / NSINC); 
	else	if(IndexSINC == NSINC - 1)	CFD = IndexMin * N + IndexSINC * (N / NSINC) + (CFDed[IndexSINC] - BaseLine - M) / (CFDed[IndexSINC] - Wave[IndexMin + 1]) * (N / NSINC); 
	else	CFD = IndexMin * N + IndexSINC * (N / NSINC) + (CFDed[IndexSINC] - BaseLine - M) / (CFDed[IndexSINC] - CFDed[IndexSINC + 1]) * (N / NSINC); 
	return CFD; 
}

void imp::ReadOut(std::ifstream &LaFile, std::ifstream &GeFile, std::ifstream &TimeFile)
{
	for(Short_t i = 0; i < LaNum; i++)
	{
		kbLaBr[0][i] = 1; 
		kbLaBr[1][i] = 0; 
		kbLaBr[2][i] = 0; 
	}
	for(Short_t i = 0; i < GeNum; i++)
	{
		kbHPGe[0][i] = 1; 
		kbHPGe[1][i] = 0; 
	}
	std::string line; 
	Int_t j = 0; 
		std::getline(LaFile, line); 
	while (std::getline(LaFile, line)) { // 循环读取每一行
		Double_t Number; // 创建整数变量
		std::stringstream ss(line); // 将字符串转换为数据流
		Int_t i = 0; 
				ss >> Number; 
		while (ss >> Number) { // 循环读取每一个数字
			kbLaBr[i][j] = Number; // 将数字添加到一维数组中
			i += 1; 
		}
		j += 1; 
	}

	std::string Line; 
	j = 0;	
	while (std::getline(GeFile, Line)) { // 循环读取每一行
		Double_t Number; // 创建整数变量
		std::stringstream ss(line); // 将字符串转换为数据流
		Int_t i = 0; 
		while (ss >> Number) { // 循环读取每一个数字
			kbHPGe[i][j] = Number; // 将数字添加到一维数组中
			i += 1; 
		}
		j += 1; 
	}

	std::string LIne; 
	j = 0; 
		std::getline(TimeFile, LIne); 
	while (std::getline(TimeFile, LIne))
	{
		Double_t Number; 
		std::stringstream ss(LIne); 
				ss >> Number; 
		ss >> Number; 
		TimeLa[j] = Number; 
		j++; 
	}
	
}
