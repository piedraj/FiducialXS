//------------------------------------------------------------------------------
//
// index < -1 only loads the macro
// index = -1 runs on all events
// index >  0 runs on _index.root file
//
//------------------------------------------------------------------------------
void runFiducialXS(int index = -999)
{
  gInterpreter->LoadMacro("FiducialXS.C+");

  if (index < -1) return;

  FiducialXS fxs(index);

  fxs.Loop(index);
}
