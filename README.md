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

    root -l -b -q 'runFiducialXS.C(-999)'
    qsub submitFiducialXS.sge


Extract the final values
====

    hadd -f rootfiles/fiducial.root rootfiles/fiducial_*.root
    root -l -b -q extractFiducialXS.C


Results
====

    --------------------------------------------------
     N(POWHEG events)                      = 1283520
    --------------------------------------------------
     N(ttbar inclusive)                    = 1283520
     N(ttbar inclusive selected)           = 8241
     N(ttbar -> emu selected)              = 8154 (not used)
     total efficiency eff                  =  0.64%
    --------------------------------------------------
     N(ttbar -> emu)                       = 41473
     N(ttbar -> emu fiducial)              = 21444
     N(ttbar -> emu fiducial selected)     = 8052
     N(ttbar -> emu non-fiducial selected) = 102
     fiducial efficiency eff_fid           = 37.55%
     contamination fraction f              =  1.27%
    --------------------------------------------------
     xs_fid/xs = eff/eff_fid/(1+f)         = 0.0169
    --------------------------------------------------

