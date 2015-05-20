//------------------------------------------------------------------------------
//
// index = -2 only loads the macro
// index = -1 runs on all events
// index >  0 runs on _index.root file
//
//------------------------------------------------------------------------------
void runFiducialXS(Int_t index = -2)
{
  gInterpreter->LoadMacro("FiducialXS.C+");

  if (index > -2)
    {
      FiducialXS fxs;

      fxs.Loop(index);
    }
}
