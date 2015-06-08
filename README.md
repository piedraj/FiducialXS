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
    qstat -u piedra


Extract the final values
====

    hadd -f rootfiles/fiducial.root rootfiles/fiducial_*.root
    root -l -b -q extractFiducialXS.C


Results
====

    --------------------------------------------------
     N(events)                        = 4995994
    --------------------------------------------------
     N(ttbar)                         = 4995994
     N(ttbar selected)                = 33211
     total efficiency eff             =  0.66%
    --------------------------------------------------
     N(fiducial)                      = 83107
     N(fiducial selected)             = 32497
     N(non-fiducial selected)         = 714
     fiducial efficiency eff_fid      = 39.10%
     contamination fraction f         =  2.20%
    --------------------------------------------------
     xs_fid/xs = N(ttbar)/N(fiducial) = 0.0166
    --------------------------------------------------
