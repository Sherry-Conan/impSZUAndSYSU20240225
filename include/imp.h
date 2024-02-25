#ifndef imp_h
#define imp_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1I.h>
#include <iostream>
#include <fstream>
#include <THnSparse.h>
#include <TH2D.h>

#include "TClonesArray.h"
#include "TObject.h"

class imp {
public :
	TTree *fChain = nullptr; 
	Int_t		fCurrent; 
	TFile *f = nullptr; 

	static constexpr Int_t 		kMaxdata_XIA = 4; // 探测器数量
	static const		 Short_t	MaxTrace = 124; //波形长度
	static const		 Short_t	GeNum = 16; // HPGe+Colover*4
	static const		 Short_t	LaNum = 28; // 数字电路：LaBr

// Declaration of leaf types
	Int_t           nXIA;
	Int_t           data_XIA_;
	UInt_t          data_XIA_fUniqueID[kMaxdata_XIA];   //[data_XIA_]
	UInt_t          data_XIA_fBits[kMaxdata_XIA];   //[data_XIA_]
	Long64_t        data_XIA_Event_ts[kMaxdata_XIA];   //[data_XIA_]
	Long64_t        data_XIA_Ext_ts[kMaxdata_XIA];   //[data_XIA_]
	UShort_t        data_XIA_finish_code[kMaxdata_XIA];   //[data_XIA_]
	UShort_t        data_XIA_crate_id[kMaxdata_XIA];   //[data_XIA_]
	UShort_t        data_XIA_channel[kMaxdata_XIA];   //[data_XIA_]
	UShort_t        data_XIA_Energy[kMaxdata_XIA];   //[data_XIA_]
	Double_t        data_XIA_CFD[kMaxdata_XIA];   //[data_XIA_]
	UShort_t        data_XIA_CFD_Bit[kMaxdata_XIA];   //[data_XIA_]
	UInt_t          data_XIA_nESum[kMaxdata_XIA];   //[data_XIA_]
	UInt_t         *data_XIA_ESum[kMaxdata_XIA];   //[data_XIA_nESum]
	UInt_t          data_XIA_nQDC[kMaxdata_XIA];   //[data_XIA_]
	UInt_t         *data_XIA_QDC[kMaxdata_XIA];   //[data_XIA_nQDC]
	UInt_t          data_XIA_trace_length[kMaxdata_XIA];   //[data_XIA_]
	UShort_t       *data_XIA_trace[kMaxdata_XIA];   //[data_XIA_trace_length]


	// List of branches
	TBranch		*b_nXIA;	//!
	TBranch		*b_data_XIA_;	//!
	TBranch		*b_data_XIA_fUniqueID;	//!
	TBranch		*b_data_XIA_fBits;	//!
	TBranch		*b_data_XIA_Event_ts;	//!
	TBranch		*b_data_XIA_Ext_ts;	//!
	TBranch		*b_data_XIA_finish_code;	//!
	TBranch		*b_data_XIA_crate_id;	//!
	TBranch		*b_data_XIA_channel;	//!
	TBranch		*b_data_XIA_Energy;	//!
	TBranch		*b_data_XIA_CFD;	//!
	TBranch		*b_data_XIA_CFD_Bit;	//!
	TBranch		*b_data_XIA_nESum;	//!
	TBranch		*b_data_XIA_ESum;	//!
	TBranch		*b_data_XIA_nQDC;	//!
	TBranch		*b_data_XIA_QDC;	//!
	TBranch		*b_data_XIA_trace_length;	//!
	TBranch		*b_data_XIA_trace;	//!

	TH1I *hEnergy[kMaxdata_XIA]; 
	Short_t FileNum; 
	Double_t kb[2][kMaxdata_XIA]; 
	Double_t CoPeak[2] = {1173, 1332}; 
	Short_t Sigma = 30; 

	Short_t 	N = 4000; // 两个点间隔多少皮秒，用以处理单位
	Short_t		LSINC = 12; //sinc插值函数，积分核宽度。
	const Short_t		NSINC = 8; // sinc插入点数，尽量别改
	Short_t 	IndexMin = 800, IndexMax = 800; // 初始化，别动
	Short_t		LBase = 30; // 基线长度
	Double_t	Rate = .25; // 上升比例
	UInt_t Wave[MaxTrace]; 

	Double_t Cfd[kMaxdata_XIA]; 
	Double_t Rise[kMaxdata_XIA]; 
	Double_t Evte[kMaxdata_XIA]; 
	UShort_t Ch[kMaxdata_XIA]; 

	Double_t	kbLaBr[2][LaNum]; 
	Double_t	kbHPGe[2][GeNum]; 
	Double_t	TimeLa[27]; 
	Double_t	kbAnal[2][2]; 

	TH1I			*hLa[kMaxdata_XIA - 1]; 
	TH2I			*hTimeMeasure; // 添加时间差-测量时间，只留0和1的符合。
	TH2I			*hTimeCh; // 添加时间差-探测器序号
	TH1I			*hHitCh; // 添加hit-探测器序号
	TH2D 			*hGe = new TH2D("Ge", "GeEngryEngry", 1600, 0, 1600, 1600, 0, 1600); 
	//Int_t			bins[3] = {1500, 1500, 2000}; 
	//Double_t	xmin[3] = {0, 0, -50}; 
	//Double_t	xmax[3] = {15000, 15000, 50}; 
	//THnSparseD* hLa = new THnSparseD("EEdT", "LaEnergy_Energy_Time", 3, bins, xmin, xmax); 

	imp(Short_t FileNum); 
	virtual ~imp();
	virtual void	Init(); 
	Double_t			CalculateCFD(Double_t M, Double_t RiseRate); 
	void					OriginalCFD(); //没什么用了，XIA原来处理出的cfd
	Double_t 			TSINC(Short_t m); 
	void 					ReadOut(std::ifstream &LaFile, std::ifstream &GeFile, std::ifstream &TimeFile); 
	void					MakeTH(TFile *File, TTree *tree); 
	void					FillGe(Long64_t Ientry); 
	Double_t			FindMax(); 
	void					Filter(Short_t ihit); 
	
	//void					FillLa(Long64_t Ientry); 
};

#endif

#ifdef imp_cxx
imp::imp(Short_t tFileNum) : FileNum(tFileNum)
{
	f = new TFile(Form("run%05d_final.root", FileNum)); 

	f->GetObject("tr", fChain); 
	Init(); 
}

imp::~imp()
{
}
void imp::Init()
{
	fChain->SetBranchAddress("nXIA", &nXIA, &b_nXIA);
	fChain->SetBranchAddress("data_XIA", &data_XIA_, &b_data_XIA_);
	fChain->SetBranchAddress("data_XIA.fUniqueID", data_XIA_fUniqueID, &b_data_XIA_fUniqueID);
	fChain->SetBranchAddress("data_XIA.fBits", data_XIA_fBits, &b_data_XIA_fBits);
	fChain->SetBranchAddress("data_XIA.Event_ts", data_XIA_Event_ts, &b_data_XIA_Event_ts);
	fChain->SetBranchAddress("data_XIA.Ext_ts", data_XIA_Ext_ts, &b_data_XIA_Ext_ts);
	fChain->SetBranchAddress("data_XIA.finish_code", data_XIA_finish_code, &b_data_XIA_finish_code);
	fChain->SetBranchAddress("data_XIA.crate_id", data_XIA_crate_id, &b_data_XIA_crate_id);
	fChain->SetBranchAddress("data_XIA.channel", data_XIA_channel, &b_data_XIA_channel);
	fChain->SetBranchAddress("data_XIA.Energy", data_XIA_Energy, &b_data_XIA_Energy);
	fChain->SetBranchAddress("data_XIA.CFD", data_XIA_CFD, &b_data_XIA_CFD);
	fChain->SetBranchAddress("data_XIA.CFD_Bit", data_XIA_CFD_Bit, &b_data_XIA_CFD_Bit);
	fChain->SetBranchAddress("data_XIA.nESum", data_XIA_nESum, &b_data_XIA_nESum);
	fChain->SetBranchAddress("data_XIA.ESum", data_XIA_ESum, &b_data_XIA_ESum);
	fChain->SetBranchAddress("data_XIA.nQDC", data_XIA_nQDC, &b_data_XIA_nQDC);
	fChain->SetBranchAddress("data_XIA.QDC", data_XIA_QDC, &b_data_XIA_QDC);
	fChain->SetBranchAddress("data_XIA.trace_length", data_XIA_trace_length, &b_data_XIA_trace_length);
	fChain->SetBranchAddress("data_XIA.trace", data_XIA_trace, &b_data_XIA_trace);
}
#endif // #ifdef imp_cxx