Get ROOT
====

    ssh -Y gridui.ifca.es -o ServerAliveInterval=240
    source /cvmfs/cms.cern.ch/cmsset_default.sh
    pushd /gpfs/csic_users/piedra/CMSSW_7_3_0/src
    cmsenv
    popd


Get the material
====

    git clone https://github.com/piedraj/FiducialXS


Submit the jobs
====

    root -l -b -q runFiducialXS.C
    rm -rf rootfiles/fiducial_*.root
    qsub submitFiducialXS.sge


Extract the final values
====

    hadd -f rootfiles/fiducial.root rootfiles/fiducial_*.root
    root -l -b -q extractFiducialXS.C


Results
====

    --------------------------------------------------
     N(events)                        = 4983143
    --------------------------------------------------
     N(ttbar)                         = 4983143
     N(ttbar selected)                = 33904
     total efficiency eff             =  0.68%
    --------------------------------------------------
     N(fiducial)                      = 82861
     N(fiducial selected)             = 33138
     N(non-fiducial selected)         = 766
     fiducial efficiency eff_fid      = 39.99%
     contamination fraction f         =  2.31%
    --------------------------------------------------
     xs_fid/xs = N(ttbar)/N(fiducial) = 0.0166
    --------------------------------------------------
