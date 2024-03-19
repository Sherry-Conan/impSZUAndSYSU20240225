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

	static constexpr Int_t 		kMaxdata_XIA = 14; // 探测器数量
	static const		 Short_t	MaxTrace = 150; //波形长度
	static const		 Short_t	GeNum = 16; // HPGe+Colover*4
	static const		 Short_t	LaNum = 28; // 数字电路：LaBr

// Declaration of leaf types
	Int_t           nXIA;
	Int_t           data_XIA_;
	Long64_t        data_XIA_Event_ts[kMaxdata_XIA];   //[data_XIA_]
	UShort_t        data_XIA_crate_id[kMaxdata_XIA];   //[data_XIA_]
	UShort_t        data_XIA_channel[kMaxdata_XIA];   //[data_XIA_]
	UShort_t        data_XIA_Energy[kMaxdata_XIA];   //[data_XIA_]
	Double_t        data_XIA_CFD[kMaxdata_XIA];   //[data_XIA_]
	UInt_t          data_XIA_trace_length[kMaxdata_XIA];   //[data_XIA_]
	UShort_t       *data_XIA_trace[kMaxdata_XIA];   //[data_XIA_trace_length]


	// List of branches
	TBranch		*b_nXIA;	//!
	TBranch		*b_data_XIA_;	//!
	TBranch		*b_data_XIA_Event_ts;	//!
	TBranch		*b_data_XIA_crate_id;	//!
	TBranch		*b_data_XIA_channel;	//!
	TBranch		*b_data_XIA_Energy;	//!
	TBranch		*b_data_XIA_CFD;	//!
	TBranch		*b_data_XIA_CFD_Bit;	//!
	TBranch		*b_data_XIA_trace_length;	//!
	TBranch		*b_data_XIA_trace;	//!

	TH1I *hEnergy[kMaxdata_XIA]; 
	Short_t FileNum; 
	Double_t kb[2][kMaxdata_XIA]; 
	Double_t CoPeak[2] = {1173, 1332}; 
	Short_t Sigma = 30; 

	Short_t 	N = 2000; // 两个点间隔多少皮秒，用以处理单位
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
	Double_t e[kMaxdata_XIA]; 

	Double_t	kbLaBr[3][LaNum]; 
	Double_t	kbHPGe[2][GeNum]; 
	Double_t	TimeLa[LaNum - 1]; 
	Double_t	kbAnal[2][2]; 

	TH1I			*hLa[kMaxdata_XIA]; 
	TH2I			*hTimeMeasure; // 添加时间差-测量时间，只留0和1的符合。
	TH2I			*hTimeCh; // 添加时间差-探测器序号
	TH1I			*hHitCh; // 添加hit-探测器序号
	TH2D 			*hGe = new TH2D("Ge", "GeEngryEngry", 1600, 0, 1600, 1600, 0, 1600); 
	//Int_t			bins[3] = {1500, 1500, 2000}; 
	//Double_t	xmin[3] = {0, 0, -50}; 
	//Double_t	xmax[3] = {15000, 15000, 50}; 
	//THnSparseD* hLa = new THnSparseD("EEdT", "LaEnergy_Energy_Time", 3, bins, xmin, xmax); 

	imp(Short_t tFileNum, char *FileName); 
	virtual ~imp();
	virtual void	Init(); 
	Double_t			CalculateCFD(Double_t M, Double_t RiseRate, Double_t BaseLine); 
	void					OriginalCFD(); //没什么用了，XIA原来处理出的cfd
	Double_t 			TSINC(Short_t m); 
	void 					ReadOut(std::ifstream &LaFile, std::ifstream &GeFile, std::ifstream &TimeFile); 
	void					MakeTH(TFile *File, TTree *tree); 
	void					FillGe(Long64_t Ientry); 
	Double_t			FindMax(); 
	void					Filter(Short_t ihit); 
	Double_t			FindBaseLine(); 
	
	//void					FillLa(Long64_t Ientry); 
};

#endif

#ifdef imp_cxx
imp::imp(Short_t tFileNum, char *FileName) : FileNum(tFileNum)
{
	f = new TFile(Form("../data/%s%05d_final.root", FileName, FileNum)); 

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
	fChain->SetBranchAddress("data_XIA.Event_ts", data_XIA_Event_ts, &b_data_XIA_Event_ts);
	fChain->SetBranchAddress("data_XIA.crate_id", data_XIA_crate_id, &b_data_XIA_crate_id);
	fChain->SetBranchAddress("data_XIA.channel", data_XIA_channel, &b_data_XIA_channel);
	fChain->SetBranchAddress("data_XIA.Energy", data_XIA_Energy, &b_data_XIA_Energy);
	fChain->SetBranchAddress("data_XIA.CFD", data_XIA_CFD, &b_data_XIA_CFD);
	fChain->SetBranchAddress("data_XIA.trace_length", data_XIA_trace_length, &b_data_XIA_trace_length);
	fChain->SetBranchAddress("data_XIA.trace", data_XIA_trace, &b_data_XIA_trace);
}
#endif // #ifdef imp_cxx
