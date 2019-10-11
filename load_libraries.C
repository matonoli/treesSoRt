void load_libraries() {
  gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/ -I$ALICE_ROOT/include -I$ALICE_PHYSICS/PWGLF/SPECTRA/PiKaPr/TPCTOFfits/ -I$PWD/Scripts/");
  
  TString libraries = "OADB;PWGPP;PWGLFspectra";
  TObjArray *libTokens = libraries.Tokenize(";");
  
  for (Int_t i=0; i<libTokens->GetEntries(); i++)
    gSystem->Load(Form("lib%s", libTokens->At(i)->GetName()));
}
