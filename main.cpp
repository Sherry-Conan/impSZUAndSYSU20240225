#include"imp.h"

int main(int argc, char const *argv[])
{
	Int_t	FileNum = std::stoi(argv[1]); 
	char* FileName = (char*)argv[2]; 

	imp *impsysu = new imp(FileNum, FileName); 

	std::ifstream LaFile(Form("LaFile.txt")); 
	std::ifstream GeFile(Form("GeFile.txt")); 
	std::ifstream TimeFile(Form("TimeFile.txt")); 

	impsysu -> ReadOut(LaFile, GeFile, TimeFile); 

	TFile *File = new TFile(Form("../data/Ana_%s%04d.root", argv[2], FileNum), "recreate"); 
	TTree *tree = new TTree("tree", "tree"); 

	impsysu -> MakeTH(File, tree); 
	File -> Close(); 

	return 0;
}
