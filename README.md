# L1RateComputation
Scripts to be run in cmssw to compute the L1 rate for VBF+Mu

These scripts are used to compute the L1 rates from official Zero Bias L1 Ntuples as function of
- jet pt 1
- jet pt 2
- mjj
- mu pt

# Instructions

    cmsrel CMSSW_13_0_0_pre2
    cd CMSSW_13_0_0_pre2/src
    cmsenv
    git cms-init
    git remote add cms-l1t-offline git@github.com:cms-l1t-offline/cmssw.git
    git fetch cms-l1t-offline l1t-integration-CMSSW_13_0_0_pre2
    git cms-merge-topic -u cms-l1t-offline:l1t-integration-v141-CMSSW_13_0_0_pre2
    git cms-addpkg L1Trigger/L1TCommon
    git cms-addpkg L1Trigger/L1TMuon
    git clone https://github.com/cms-l1t-offline/L1Trigger-L1TMuon.git L1Trigger/L1TMuon/data
    git cms-addpkg L1Trigger/L1TCalorimeter
    git cms-addpkg L1Trigger/L1TCaloLayer1
    git cms-addpkg L1Trigger/L1TNtuples
    git clone https://github.com/cms-l1t-offline/L1Trigger-L1TCalorimeter.git L1Trigger/L1TCalorimeter/data

    scram b -j 8

Clone the repository (then move the script so that they are in the folder 'CMSSW_13_0_0_pre2/src/L1Trigger/L1TNtuples').

    root -l -b 
    .L Rate_ZeroBias_Run362616.C+
    Rate()

    root -l -b 
    .L Plot_Rates_2D.C+
    Plot_2D()