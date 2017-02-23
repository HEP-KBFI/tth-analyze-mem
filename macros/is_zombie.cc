// call it with
//   root -b -l -q "is_zombie.cc(\"$PATH_TO_ROOT_FILE\")"
// where $PATH_TO_ROOT_FILE is the path to the ROOT file
// you want to test the macro against (do not use BASH
// variables, though)
int
is_zombie(const std::string & filename)
{
  if(! gSystem -> AccessPathName(filename.c_str()))
  {
    TFile file(filename.c_str());
    if(file.IsZombie())
    {
      std::cerr << "File " << filename << " is corrupted\n";
      return EXIT_FAILURE;
    }
    if(file.TestBit(TFile::kRecovered))
    {
      std::cout << "File " << filename << " was corrputed, but now recovered\n";
    }
    return EXIT_SUCCESS;
  }
  else
  {
    std::cerr << "No such file: " << filename << '\n';
    return EXIT_FAILURE;
  }
}
