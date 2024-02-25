#include"imp.h"

Int_t main(int argc, char const *argv[])
{
	Int_t	FileNum = std::stoi(argv[1]); 

	imp *impsysu = new imp(FileNum); 

	std::ifstream LaFile(Form("LaFile.txt")); 
	std::ifstream GeFile(Form("GeFile.txt")); 
	std::ifstream TimeFile(Form("TimeFile.txt")); 

	impsysu -> ReadOut(LaFile, GeFile, TimeFile); 

	TFile *File = new TFile(Form("ana%04d.root", FileNum), "recreate"); 
	TTree *tree = new TTree("tree", "tree"); 

	impsysu -> MakeTH(File, tree); 
	File -> Close(); 

	return 0;
}