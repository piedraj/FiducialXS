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

