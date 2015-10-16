void extractFiducialXS()
{
  TFile* file = new TFile("rootfiles/fiducial.root", "read");

  TH1F* h_events                     = (TH1F*)file->Get("h_events");
  TH1F* h_ttbar                      = (TH1F*)file->Get("h_ttbar");
  TH1F* h_ttbar_fiducial             = (TH1F*)file->Get("h_ttbar_fiducial");
  TH1F* h_ttbar_fiducial_selected    = (TH1F*)file->Get("h_ttbar_fiducial_selected");
  TH1F* h_ttbar_nonfiducial_selected = (TH1F*)file->Get("h_ttbar_nonfiducial_selected");


  // Summary
  //----------------------------------------------------------------------------
  float n_events                     = h_events                    ->Integral();
  float n_ttbar                      = h_ttbar                     ->Integral();
  float n_ttbar_fiducial             = h_ttbar_fiducial            ->Integral();
  float n_ttbar_fiducial_selected    = h_ttbar_fiducial_selected   ->Integral();
  float n_ttbar_nonfiducial_selected = h_ttbar_nonfiducial_selected->Integral();

  float eff      = (n_ttbar > 0) ? (n_ttbar_fiducial_selected + n_ttbar_nonfiducial_selected) / n_ttbar : -999;
  float eff_fid  = (n_ttbar_fiducial > 0) ? n_ttbar_fiducial_selected / n_ttbar_fiducial : -999;
  float f        = (n_ttbar_fiducial_selected > 0) ? n_ttbar_nonfiducial_selected / n_ttbar_fiducial_selected : -999;
  float xs_ratio = (n_ttbar > 0) ? n_ttbar_fiducial / n_ttbar : -999;

  printf("\n");
  printf("--------------------------------------------------\n");
  printf(" N(events)                        = %.0f\n", n_events);
  printf("--------------------------------------------------\n");
  printf(" N(ttbar)                         = %.0f\n", n_ttbar);
  printf(" N(ttbar selected)                = %.0f\n", n_ttbar_fiducial_selected+n_ttbar_nonfiducial_selected);
  printf(" total efficiency eff             = %5.2f%s\n", 1e2 * eff, "%");
  printf("--------------------------------------------------\n");
  printf(" N(fiducial)                      = %.0f\n", n_ttbar_fiducial);
  printf(" N(fiducial selected)             = %.0f\n", n_ttbar_fiducial_selected);
  printf(" N(non-fiducial selected)         = %.0f\n", n_ttbar_nonfiducial_selected);
  printf(" fiducial efficiency eff_fid      = %5.2f%s\n", 1e2 * eff_fid, "%");
  printf(" contamination fraction f         = %5.2f%s\n", 1e2 * f, "%");
  printf("--------------------------------------------------\n");
  printf(" xs_fid/xs = N(fiducial)/N(ttvar) = %.4f\n", xs_ratio);
  printf("--------------------------------------------------\n");
  printf("\n");
}
