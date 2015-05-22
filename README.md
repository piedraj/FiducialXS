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

    --------------------------------------------------
     N(events)                             = 25446992
    --------------------------------------------------
     N(ttbar inclusive)                    = 25446992
     N(ttbar inclusive selected)           = 189234
     N(ttbar -> emu selected)              = 187063 (not used)
     total efficiency eff                  =  0.74%
    --------------------------------------------------
     N(ttbar -> emu)                       = 870888
     N(ttbar -> emu fiducial)              = 455738
     N(ttbar -> emu fiducial selected)     = 185285
     N(ttbar -> emu non-fiducial selected) = 1778
     fiducial efficiency eff_fid           = 40.66%
     contamination fraction f              =  0.96%
    --------------------------------------------------
     xs_fid/xs = eff/eff_fid/(1+f)         = 0.0181
    --------------------------------------------------

