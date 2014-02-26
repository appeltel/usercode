import FWCore.ParameterSet.Config as cms



def getPy6ProcParameters(genTypePtHatRange):
    para = {
        # NSD (beyon pt>5~10 GeV/c, it's safe to say NSD = MB, but below, it's not)
        'NSD_0_to_5':['MSEL=1','CKIN(3)=0','CKIN(4)=5'], 
        'NSD_5_to_10':['MSEL=1','CKIN(3)=5','CKIN(4)=10'],
        'NSD_10_to_20':['MSEL=1','CKIN(3)=10','CKIN(4)=20'],
        'NSD_0_to_10':['MSEL=1','CKIN(3)=0','CKIN(4)=10'],
        'NSD_0_to_15':['MSEL=1','CKIN(3)=0','CKIN(4)=15'],
        'NSD_0_to_20':['MSEL=1','CKIN(3)=0','CKIN(4)=20'],
        'NSD_15_to_20':['MSEL=1','CKIN(3)=15','CKIN(4)=20'],
        'NSD_20_to_30':['MSEL=1','CKIN(3)=20','CKIN(4)=30'],
        'NSD_30_to_50':['MSEL=1','CKIN(3)=30','CKIN(4)=50'],
        'NSD_30_to_inf':['MSEL=1','CKIN(3)=30'],
        'NSD_50_to_80':['MSEL=1','CKIN(3)=50','CKIN(4)=80'],
        'NSD_80_to_inf':['MSEL=1','CKIN(3)=80'],
        'NSD_80_to_120':['MSEL=1','CKIN(3)=80','CKIN(4)=120'],
        'NSD_120_to_170':['MSEL=1','CKIN(3)=120','CKIN(4)=170'],
        'NSD_170_to_230':['MSEL=1','CKIN(3)=170','CKIN(4)=230'],
        'NSD_230_to_300':['MSEL=1','CKIN(3)=230','CKIN(4)=300'],
        'NSD_300_to_380':['MSEL=1','CKIN(3)=300','CKIN(4)=380'],
        'NSD_380_to_470':['MSEL=1','CKIN(3)=380','CKIN(4)=470'],
        'NSD_380_to_inf':['MSEL=1','CKIN(3)=380'],
        'NSD_470_to_inf':['MSEL=1','CKIN(3)=470'],
        'NSD_470_to_570':['MSEL=1','CKIN(3)=470','CKIN(4)=570'],
        'NSD_570_to_670':['MSEL=1','CKIN(3)=570','CKIN(4)=670'],
        'NSD_670_to_800':['MSEL=1','CKIN(3)=670','CKIN(4)=800'],
        'NSD_800_to_930':['MSEL=1','CKIN(3)=800','CKIN(4)=930'],
        'NSD_930_to_inf':['MSEL=1','CKIN(3)=930'],
        'NSD_930_to_1000':['MSEL=1','CKIN(3)=930','CKIN(4)=1000'],
        'NSD_1000_to_1100':['MSEL=1','CKIN(3)=1000','CKIN(4)=1100'],
        'NSD_1100_to_1200':['MSEL=1','CKIN(3)=1100','CKIN(4)=1200'],
        'NSD_1200_to_1400':['MSEL=1','CKIN(3)=1200','CKIN(4)=1400'],
        'NSD_1400_to_1600':['MSEL=1','CKIN(3)=1400','CKIN(4)=1600'],
        'NSD_1600_to_1800':['MSEL=1','CKIN(3)=1600','CKIN(4)=1800'],
        'NSD_1800_to_2000':['MSEL=1','CKIN(3)=1800','CKIN(4)=2000'],
        'NSD_2000_to_2200':['MSEL=1','CKIN(3)=2000','CKIN(4)=2200'],
        'NSD_2000_to_inf':['MSEL=1','CKIN(3)=2000'],
        'NSD_2200_to_2400':['MSEL=1','CKIN(3)=2200','CKIN(4)=2400'],
        'NSD_2400_to_2600':['MSEL=1','CKIN(3)=2400','CKIN(4)=2600'],
        'NSD_2600_to_2800':['MSEL=1','CKIN(3)=2600','CKIN(4)=2800'],
        'NSD_2800_to_3000':['MSEL=1','CKIN(3)=2800','CKIN(4)=3000'],
        'NSD_2800_to_2900':['MSEL=1','CKIN(3)=2800','CKIN(4)=2900'],
        'NSD_3000_to_3200':['MSEL=1','CKIN(3)=3000','CKIN(4)=3200'],
        'NSD_2900_to_inf':['MSEL=1','CKIN(3)=2900'],
        'NSD_3000_to_inf':['MSEL=1','CKIN(3)=3000'],
        'NSD_3200_to_inf':['MSEL=1','CKIN(3)=3200'],
        'NSD_3200_to_3600':['MSEL=1','CKIN(3)=3200','CKIN(4)=3600'],
        'NSD_1200_to_inf':['MSEL=1','CKIN(3)=1200'],
        'NSD_0_to_inf':['MSEL=1','CKIN(3)=0'],
        'TOP':['MSEL=6','PMAS(6,1) = 173.5'],
        # explicit NSD
        'ENSD_0_to_5':['MSEL=0','MSUB(11)=1','MSUB(12)=1','MSUB(13)=1','MSUB(28)=1','MSUB(53)=1','MSUB(68)=1','MSUB(94)=1','MSUB(95)=1','CKIN(3)=0','CKIN(4)=5'],
        'ENSD_5_to_10':['MSEL=0','MSUB(11)=1','MSUB(12)=1','MSUB(13)=1','MSUB(28)=1','MSUB(53)=1','MSUB(68)=1','MSUB(94)=1','MSUB(95)=1','CKIN(3)=5','CKIN(4)=10'],
        'ENSD_10_to_20':['MSEL=0','MSUB(11)=1','MSUB(12)=1','MSUB(13)=1','MSUB(28)=1','MSUB(53)=1','MSUB(68)=1','MSUB(94)=1','MSUB(95)=1','CKIN(3)=10','CKIN(4)=20'],
        'ENSD_0_to_10':['MSEL=0','MSUB(11)=1','MSUB(12)=1','MSUB(13)=1','MSUB(28)=1','MSUB(53)=1','MSUB(68)=1','MSUB(94)=1','MSUB(95)=1','CKIN(3)=0','CKIN(4)=10'],
        'ENSD_0_to_15':['MSEL=0','MSUB(11)=1','MSUB(12)=1','MSUB(13)=1','MSUB(28)=1','MSUB(53)=1','MSUB(68)=1','MSUB(94)=1','MSUB(95)=1','CKIN(3)=0','CKIN(4)=15'],
        'ENSD_0_to_20':['MSEL=0','MSUB(11)=1','MSUB(12)=1','MSUB(13)=1','MSUB(28)=1','MSUB(53)=1','MSUB(68)=1','MSUB(94)=1','MSUB(95)=1','CKIN(3)=0','CKIN(4)=20'],
        'ENSD_15_to_20':['MSEL=0','MSUB(11)=1','MSUB(12)=1','MSUB(13)=1','MSUB(28)=1','MSUB(53)=1','MSUB(68)=1','MSUB(94)=1','MSUB(95)=1','CKIN(3)=15','CKIN(4)=20'],
        'ENSD_20_to_30':['MSEL=0','MSUB(11)=1','MSUB(12)=1','MSUB(13)=1','MSUB(28)=1','MSUB(53)=1','MSUB(68)=1','MSUB(94)=1','MSUB(95)=1','CKIN(3)=20','CKIN(4)=30'],
        'ENSD_30_to_50':['MSEL=0','MSUB(11)=1','MSUB(12)=1','MSUB(13)=1','MSUB(28)=1','MSUB(53)=1','MSUB(68)=1','MSUB(94)=1','MSUB(95)=1','CKIN(3)=30','CKIN(4)=50'],
        'ENSD_50_to_80':['MSEL=0','MSUB(11)=1','MSUB(12)=1','MSUB(13)=1','MSUB(28)=1','MSUB(53)=1','MSUB(68)=1','MSUB(94)=1','MSUB(95)=1','CKIN(3)=50','CKIN(4)=80'],
        'ENSD_80_to_inf':['MSEL=0','MSUB(11)=1','MSUB(12)=1','MSUB(13)=1','MSUB(28)=1','MSUB(53)=1','MSUB(68)=1','MSUB(94)=1','MSUB(95)=1','CKIN(3)=80'],
        'ENSD_80_to_120':['MSEL=0','MSUB(11)=1','MSUB(12)=1','MSUB(13)=1','MSUB(28)=1','MSUB(53)=1','MSUB(68)=1','MSUB(94)=1','MSUB(95)=1','CKIN(3)=80','CKIN(4)=120'],
        'ENSD_120_to_170':['MSEL=0','MSUB(11)=1','MSUB(12)=1','MSUB(13)=1','MSUB(28)=1','MSUB(53)=1','MSUB(68)=1','MSUB(94)=1','MSUB(95)=1','CKIN(3)=120','CKIN(4)=170'],
        'ENSD_170_to_230':['MSEL=0','MSUB(11)=1','MSUB(12)=1','MSUB(13)=1','MSUB(28)=1','MSUB(53)=1','MSUB(68)=1','MSUB(94)=1','MSUB(95)=1','CKIN(3)=170','CKIN(4)=230'],
        'ENSD_230_to_300':['MSEL=0','MSUB(11)=1','MSUB(12)=1','MSUB(13)=1','MSUB(28)=1','MSUB(53)=1','MSUB(68)=1','MSUB(94)=1','MSUB(95)=1','CKIN(3)=230','CKIN(4)=300'],
        'ENSD_300_to_380':['MSEL=0','MSUB(11)=1','MSUB(12)=1','MSUB(13)=1','MSUB(28)=1','MSUB(53)=1','MSUB(68)=1','MSUB(94)=1','MSUB(95)=1','CKIN(3)=300','CKIN(4)=380'],
        'ENSD_380_to_470':['MSEL=0','MSUB(11)=1','MSUB(12)=1','MSUB(13)=1','MSUB(28)=1','MSUB(53)=1','MSUB(68)=1','MSUB(94)=1','MSUB(95)=1','CKIN(3)=380','CKIN(4)=470'],
        'ENSD_470_to_inf':['MSEL=0','MSUB(11)=1','MSUB(12)=1','MSUB(13)=1','MSUB(28)=1','MSUB(53)=1','MSUB(68)=1','MSUB(94)=1','MSUB(95)=1','CKIN(3)=470'],
        # MB
        'MB_0_to_5':['MSEL=0','MSUB(11)=1','MSUB(12)=1','MSUB(13)=1','MSUB(28)=1','MSUB(53)=1','MSUB(68)=1','MSUB(92)=1','MSUB(93)=1','MSUB(94)=1','MSUB(95)=1','CKIN(3)=0','CKIN(4)=5'],
        'MB_0_to_10':['MSEL=0','MSUB(11)=1','MSUB(12)=1','MSUB(13)=1','MSUB(28)=1','MSUB(53)=1','MSUB(68)=1','MSUB(92)=1','MSUB(93)=1','MSUB(94)=1','MSUB(95)=1','CKIN(3)=0','CKIN(4)=10'],
        'MB_10_to_inf':['MSEL=0','MSUB(11)=1','MSUB(12)=1','MSUB(13)=1','MSUB(28)=1','MSUB(53)=1','MSUB(68)=1','MSUB(92)=1','MSUB(93)=1','MSUB(94)=1','MSUB(95)=1','CKIN(3)=10'],
        'MB_10_to_20':['MSEL=0','MSUB(11)=1','MSUB(12)=1','MSUB(13)=1','MSUB(28)=1','MSUB(53)=1','MSUB(68)=1','MSUB(92)=1','MSUB(93)=1','MSUB(94)=1','MSUB(95)=1','CKIN(3)=10','CKIN(4)=20'],
        'MB_0_to_15':['MSEL=0','MSUB(11)=1','MSUB(12)=1','MSUB(13)=1','MSUB(28)=1','MSUB(53)=1','MSUB(68)=1','MSUB(92)=1','MSUB(93)=1','MSUB(94)=1','MSUB(95)=1','CKIN(3)=0','CKIN(4)=15'],
        'MB_0_to_20':['MSEL=0','MSUB(11)=1','MSUB(12)=1','MSUB(13)=1','MSUB(28)=1','MSUB(53)=1','MSUB(68)=1','MSUB(92)=1','MSUB(93)=1','MSUB(94)=1','MSUB(95)=1','CKIN(3)=0','CKIN(4)=20'],
        'MB_0_to_25':['MSEL=0','MSUB(11)=1','MSUB(12)=1','MSUB(13)=1','MSUB(28)=1','MSUB(53)=1','MSUB(68)=1','MSUB(92)=1','MSUB(93)=1','MSUB(94)=1','MSUB(95)=1','CKIN(3)=0','CKIN(4)=25']
        }
    print "PYTHIA process parameters: ",para[genTypePtHatRange]
    return para[genTypePtHatRange]

def getPy8ProcParameters(genTypePtHatRange):
    para = {
        # 0 to X (MB)
        'MB_0_to_10':['Main:timesAllowErrors = 10000','ParticleDecays:limitTau0 = on','ParticleDecays:tauMax = 10','SoftQCD:minBias = on','SoftQCD:singleDiffractive = on','SoftQCD:doubleDiffractive = on','PhaseSpace:pTHatMin = 0','PhaseSpace:pTHatMax = 10','Tune:pp 2','Tune:ee 3'],
        'MB_0_to_15':['Main:timesAllowErrors = 10000','ParticleDecays:limitTau0 = on','ParticleDecays:tauMax = 10','SoftQCD:minBias = on','SoftQCD:singleDiffractive = on','SoftQCD:doubleDiffractive = on','PhaseSpace:pTHatMin = 0','PhaseSpace:pTHatMax = 15','Tune:pp 2','Tune:ee 3'],
        # 0 to X (Explicit)
        'ENSD_0_to_5':['Main:timesAllowErrors = 10000','ParticleDecays:limitTau0 = on','ParticleDecays:tauMax = 10','SoftQCD:minBias = on','SoftQCD:singleDiffractive = off','SoftQCD:doubleDiffractive = on','PhaseSpace:pTHatMin = 0','PhaseSpace:pTHatMax = 5','Tune:pp 2','Tune:ee 3'],
        'ENSD_0_to_10':['Main:timesAllowErrors = 10000','ParticleDecays:limitTau0 = on','ParticleDecays:tauMax = 10','SoftQCD:minBias = on','SoftQCD:singleDiffractive = off','SoftQCD:doubleDiffractive = on','PhaseSpace:pTHatMin = 0','PhaseSpace:pTHatMax = 10','Tune:pp 2','Tune:ee 3'],
        'ENSD_0_to_15':['Main:timesAllowErrors = 10000','ParticleDecays:limitTau0 = on','ParticleDecays:tauMax = 10','SoftQCD:minBias = on','SoftQCD:singleDiffractive = off','SoftQCD:doubleDiffractive = on','PhaseSpace:pTHatMin = 0','PhaseSpace:pTHatMax = 15','Tune:pp 2','Tune:ee 3'],
        'ENSD_0_to_20':['Main:timesAllowErrors = 10000','ParticleDecays:limitTau0 = on','ParticleDecays:tauMax = 10','SoftQCD:minBias = on','SoftQCD:singleDiffractive = off','SoftQCD:doubleDiffractive = on','PhaseSpace:pTHatMin = 0','PhaseSpace:pTHatMax = 20','Tune:pp 2','Tune:ee 3'],
        'ENSD_0_to_30':['Main:timesAllowErrors = 10000','ParticleDecays:limitTau0 = on','ParticleDecays:tauMax = 10','SoftQCD:minBias = on','SoftQCD:singleDiffractive = off','SoftQCD:doubleDiffractive = on','PhaseSpace:pTHatMin = 0','PhaseSpace:pTHatMax = 30','Tune:pp 2','Tune:ee 3'],
        'ENSD_0_to_40':['Main:timesAllowErrors = 10000','ParticleDecays:limitTau0 = on','ParticleDecays:tauMax = 10','SoftQCD:minBias = on','SoftQCD:singleDiffractive = off','SoftQCD:doubleDiffractive = on','PhaseSpace:pTHatMin = 0','PhaseSpace:pTHatMax = 40','Tune:pp 2','Tune:ee 3'],
        # X to Y
        'NSD_5_to_10':['Main:timesAllowErrors = 10000','ParticleDecays:limitTau0 = on','ParticleDecays:tauMax = 10','HardQCD:all = on','PhaseSpace:pTHatMin = 5','PhaseSpace:pTHatMax = 10','Tune:pp 2','Tune:ee 3'],
        'NSD_15_to_20':['Main:timesAllowErrors = 10000','ParticleDecays:limitTau0 = on','ParticleDecays:tauMax = 10','HardQCD:all = on','PhaseSpace:pTHatMin = 15','PhaseSpace:pTHatMax = 20','Tune:pp 2','Tune:ee 3'],
        'NSD_10_to_20':['Main:timesAllowErrors = 10000','ParticleDecays:limitTau0 = on','ParticleDecays:tauMax = 10','HardQCD:all = on','PhaseSpace:pTHatMin = 10','PhaseSpace:pTHatMax = 20','Tune:pp 2','Tune:ee 3'],
        'NSD_20_to_30':['Main:timesAllowErrors = 10000','ParticleDecays:limitTau0 = on','ParticleDecays:tauMax = 10','HardQCD:all = on','PhaseSpace:pTHatMin = 20','PhaseSpace:pTHatMax = 30','Tune:pp 2','Tune:ee 3'],
        'NSD_30_to_50':['Main:timesAllowErrors = 10000','ParticleDecays:limitTau0 = on','ParticleDecays:tauMax = 10','HardQCD:all = on','PhaseSpace:pTHatMin = 30','PhaseSpace:pTHatMax = 50','Tune:pp 2','Tune:ee 3'],
        'NSD_50_to_80':['Main:timesAllowErrors = 10000','ParticleDecays:limitTau0 = on','ParticleDecays:tauMax = 10','HardQCD:all = on','PhaseSpace:pTHatMin = 50','PhaseSpace:pTHatMax = 80','Tune:pp 2','Tune:ee 3'],
        'NSD_80_to_inf':['Main:timesAllowErrors = 10000','ParticleDecays:limitTau0 = on','ParticleDecays:tauMax = 10','HardQCD:all = on','PhaseSpace:pTHatMin = 80','Tune:pp 2','Tune:ee 3'],
        'NSD_80_to_120':['Main:timesAllowErrors = 10000','ParticleDecays:limitTau0 = on','ParticleDecays:tauMax = 10','HardQCD:all = on','PhaseSpace:pTHatMin = 80','PhaseSpace:pTHatMax = 120','Tune:pp 2','Tune:ee 3'],
        'NSD_120_to_170':['Main:timesAllowErrors = 10000','ParticleDecays:limitTau0 = on','ParticleDecays:tauMax = 10','HardQCD:all = on','PhaseSpace:pTHatMin = 120','PhaseSpace:pTHatMax = 170','Tune:pp 2','Tune:ee 3'],
        'NSD_170_to_230':['Main:timesAllowErrors = 10000','ParticleDecays:limitTau0 = on','ParticleDecays:tauMax = 10','HardQCD:all = on','PhaseSpace:pTHatMin = 170','PhaseSpace:pTHatMax = 230','Tune:pp 2','Tune:ee 3'],
        'NSD_230_to_300':['Main:timesAllowErrors = 10000','ParticleDecays:limitTau0 = on','ParticleDecays:tauMax = 10','HardQCD:all = on','PhaseSpace:pTHatMin = 230','PhaseSpace:pTHatMax = 300','Tune:pp 2','Tune:ee 3'],
        'NSD_300_to_380':['Main:timesAllowErrors = 10000','ParticleDecays:limitTau0 = on','ParticleDecays:tauMax = 10','HardQCD:all = on','PhaseSpace:pTHatMin = 300','PhaseSpace:pTHatMax = 380','Tune:pp 2','Tune:ee 3'],
        'NSD_380_to_470':['Main:timesAllowErrors = 10000','ParticleDecays:limitTau0 = on','ParticleDecays:tauMax = 10','HardQCD:all = on','PhaseSpace:pTHatMin = 380','PhaseSpace:pTHatMax = 470','Tune:pp 2','Tune:ee 3'],
        'NSD_380_to_inf':['Main:timesAllowErrors = 10000','ParticleDecays:limitTau0 = on','ParticleDecays:tauMax = 10','HardQCD:all = on','PhaseSpace:pTHatMin = 380','Tune:pp 2','Tune:ee 3'],
        'NSD_470_to_inf':['Main:timesAllowErrors = 10000','ParticleDecays:limitTau0 = on','ParticleDecays:tauMax = 10','HardQCD:all = on','PhaseSpace:pTHatMin = 470','Tune:pp 2','Tune:ee 3']
        }
    print "PYTHIA process parameters: ",para[genTypePtHatRange]
    return para[genTypePtHatRange]
        

def updatePy6ProcParameters(gen,genTypePtHatRange,sqrtS):
    print "Center-of-mass energy: ",sqrtS
    gen.PythiaParameters.processParameters = getPy6ProcParameters(genTypePtHatRange)
    gen.comEnergy = sqrtS

def updatePy8ProcParameters(gen,genTypePtHatRange,sqrtS):
    print "Center-of-mass energy: ",sqrtS
    gen.PythiaParameters.processParameters = getPy8ProcParameters(genTypePtHatRange)
    gen.comEnergy = sqrtS
    


