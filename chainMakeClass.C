void chainMakeClass()
{
  TChain chain("demo/Tree");

  chain.Add("/gpfs/csic_projects/tier3data/TreesPHYS14/PU20bx25/Tree_TTJets_MadSpin_NoSkim_0.root");
  chain.Add("/gpfs/csic_projects/tier3data/TreesPHYS14/PU20bx25/Tree_TTJets_MadSpin_NoSkim_1.root");
  chain.Add("/gpfs/csic_projects/tier3data/TreesPHYS14/PU20bx25/Tree_TTJets_MadSpin_NoSkim_2.root");
  chain.Add("/gpfs/csic_projects/tier3data/TreesPHYS14/PU20bx25/Tree_TTJets_MadSpin_NoSkim_3.root");
  chain.Add("/gpfs/csic_projects/tier3data/TreesPHYS14/PU20bx25/Tree_TTJets_MadSpin_NoSkim_4.root");
  chain.Add("/gpfs/csic_projects/tier3data/TreesPHYS14/PU20bx25/Tree_TTJets_MadSpin_NoSkim_5.root");

  chain.MakeClass("FiducialXS");
}
