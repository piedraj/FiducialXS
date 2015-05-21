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


Load FiducialXS
====

    root -l -b -q 'runFiducialXS.C(-999)'


Submit the jobs
====

    qsub submitFiducialXS.sge
    hadd -f rootfiles/fiducial.root rootfiles/fiducial_*.root

