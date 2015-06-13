1. Get ROOT
====

    ssh -Y gridui.ifca.es -o ServerAliveInterval=240
    cd /gpfs/csic_projects/cms/$USER
    mkdir fiducial
    cd fiducial

    source /cvmfs/cms.cern.ch/cmsset_default.sh
    export SCRAM_ARCH=slc6_amd64_gcc491
    cmsrel CMSSW_7_3_0
    cd CMSSW_7_3_0/src
    cmsenv


2. Get the material
====

    git clone https://github.com/piedraj/FiducialXS


3. Submit the jobs
====

    cd FiducialXS
    root -l -b -q runFiducialXS.C
    rm -rf rootfiles/fiducial_*.root
    qsub submitFiducialXS.sge
    qstat -u $USER

Alternatively one can login to a node and run interactively

    qlogin -P l.gaes
    source /cvmfs/cms.cern.ch/cmsset_default.sh
    cd /gpfs/csic_projects/cms/$USER/fiducial/CMSSW_7_3_0/src
    cmsenv

    cd FiducialXS
    rm -rf rootfiles/fiducial_*.root
    root -l -b -q 'runFiducialXS.C(0)'
    root -l -b -q 'runFiducialXS.C(1)'
    root -l -b -q 'runFiducialXS.C(2)'
    root -l -b -q 'runFiducialXS.C(3)'
    root -l -b -q 'runFiducialXS.C(4)'


4. Extract the final values
====

    hadd -f rootfiles/fiducial.root rootfiles/fiducial_*.root
    root -l -b -q extractFiducialXS.C


5. Results
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

