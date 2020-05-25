from .cantera_tools import *

def twoStage_ideal(phi_global,phi_main,tau_global,tau_sec,airSplit=1,phiSec=None,T_fuel=300,T_ox=650,P=25, mech="gri30.xml",trace=False):
    P *= ct.one_atm
    tau_global *= milliseconds 
    tau_sec *= milliseconds 
    tau_main = tau_global - max(0,tau_sec)
    mfm, mam, mfs, mas = solve_mass_airsplit(phi_global, phi_main, airSplit=airSplit)
    if ((mfm < 0) | (mfs < 0) | (mas < 0) | (mam < 0)):
        return np.hstack([0, 0, 100000, 100000, 10000, 10000, 0, 0, 0, np.full(53, 0, dtype=np.float64), np.full(53, 0, dtype=np.float64), 100000, 100000, 100000])
    m = [mfm + mam, mfs + mas] 
    mainMassPerc = 100*(mfm + mam)/(mfm + mam + mfs + mas) 
    secMassPerc = 100*(mfs + mas)/(mfm + mam + mfs + mas)
    condColumns = ['phi_global', 'phi_main', 'tau_global', 'tau_main', 'tau_sec', 'tau_mix', 'dt', 'T_fuel', 'T_ox', 'main_mass_perc', 'sec_mass_perc', 'mfm', 'mam', 'mfs', 'mas', 'airSplit']

    if (secMassPerc <= 0.001):
        tau_main = tau_global;
        tau_sec = 0;
    if (tau_sec <= 0):
        phi_main = phi_global;

    # RUN MAIN BURNER
    tic = time.time() 

    # pdb.set_trace()
    vitReactor, mainBurnerDF = run_main_burner(phi_main, tau_main, T_fuel, T_ox, P=P, mech=mech)
    massFracs = mainBurnerDF.columns[mainBurnerDF.columns.str.contains('Y_')] 
    moleFracs = mainBurnerDF.columns[mainBurnerDF.columns.str.contains('X_')]       
    flameTime = mainBurnerDF.index.values

    toc = time.time() 
    # print("Time taken to get flame:", (toc-tic), "secs")
    if (tau_sec <= 0):
        T = mainBurnerDF['T'].iloc[-1]
        Y_values = mainBurnerDF[massFracs].iloc[-1]

        NOCOppmvd = correct_nox(np.array(mainBurnerDF[['X_NO', 'X_CO', 'X_NO2', 'X_N2O']].iloc[-1], dtype=np.float64), mainBurnerDF['X_H2O'].iloc[-1], mainBurnerDF['X_O2'].iloc[-1])
        if (trace==True): 
            mainDF = mainBurnerDF[['T', 'NOppmvd', 'COppmvd']];
            mainDF.index.name = 'time';
            return mainDF;
        return np.hstack([flameTime[-1], mainBurnerDF['T'].iloc[-1], NOCOppmvd[0], NOCOppmvd[1], NOCOppmvd[2], NOCOppmvd[3], np.array(mainBurnerDF.iloc[-1],dtype=np.float64)])
    # CREATE FUEL PARTICLE
    fuelStream = ct.Solution(mech) 
    oxStream = ct.Solution(mech) 
    fuelStream.TPX = [T_fuel, P, {'CH4':1.0}]
    oxStream.TPX = [T_ox, P, {'O2':0.21, 'N2':0.79}]     
    secondaryReactor = ct.ConstPressureReactor(mix([fuelStream, oxStream, vitReactor.thermo], [mfs, mas, mfm + mam], mech=mech)) 
    tpArray = np.array([vitReactor.thermo, secondaryReactor.thermo])
    rArray = np.array([vitReactor, secondaryReactor])

    rn = ct.ReactorNet([secondaryReactor])
    if (trace == False):
        # RUN SECONDARY REACTOR
        rn.advance(tau_sec) 

        # GET ENDPOINT DATA
        NOCOppmvd = correct_nox(secondaryReactor.thermo['NO', 'CO', 'NO2', 'N2O'].X, secondaryReactor.thermo['H2O'].X, secondaryReactor.thermo['O2'].X)
        out = np.hstack([rn.time, secondaryReactor.thermo.T, NOCOppmvd[0], NOCOppmvd[1], NOCOppmvd[2], NOCOppmvd[3], secondaryReactor.thermo.T, 0, secondaryReactor.thermo.mean_molecular_weight, secondaryReactor.thermo.Y, secondaryReactor.thermo.X, secondaryReactor.thermo.P, NOCOppmvd[0], NOCOppmvd[1]])   
        # out = np.hstack([phi_main, tau_sec, NOCOppmvd[0], NOCOppmvd[1], secondaryReactor.thermo.T])
        return out
    else:
        dt = 0.00001 * 1e-3
        sec_tList = np.arange(0, tau_sec, dt)
        columnNames = ['time', 'T', 'NOppmvd', 'COppmvd']
        secArray = np.array([None] * len(sec_tList) * len(columnNames)).reshape(len(sec_tList), len(columnNames))
        # performance note: we don't need to do this. we can simply just advance to the desired tau_main
        mainDF = mainBurnerDF[['T', 'NOppmvd', 'COppmvd']];
        mainDF.index.name = 'time';
        for i in range(0, len(sec_tList)):
            MWi = secondaryReactor.thermo.mean_molecular_weight
            NOCOppmvd = correct_nox(secondaryReactor.thermo['NO', 'CO', 'NO2', 'N2O'].X, secondaryReactor.thermo['H2O'].X, secondaryReactor.thermo['O2'].X)
            # secArray[i, :] = np.hstack([vel_final * dt, vel_final, secondaryReactor.thermo.T, 0, MWi, secondaryReactor.thermo.Y, secondaryReactor.thermo.X, secondaryReactor.thermo.P])
            rn.advance(sec_tList[i])
            secArray[i,:] = np.hstack([sec_tList[i] + tau_main, secondaryReactor.thermo.T, NOCOppmvd[0], NOCOppmvd[1]])
        sec_tList += tau_main             
        secDF = pd.DataFrame(data=secArray, columns=columnNames)
        secDF = secDF.set_index('time');
        secondaryReactor.syncState() 
        stagedDF = pd.concat([mainDF, secDF]) 
        return stagedDF;

# if __name__ == "__main__":
#     parser = ArgumentParser()
#     # parser.add_argument("enttype", type=str)  Simplifying, only assuming time for now
#     parser.add_argument("ent_main", type=float)
#     parser.add_argument("ent_sec", type=float)
#     parser.add_argument("out_dir", type=str)
#     # parser.add_argument("tau_sec", type=float)
#     parser.add_argument("phi_jet_norm", nargs='?', type=float, default=1.0)
#     parser.add_argument("phi_global", nargs='?', type=float, default=0.635)
#     parser.add_argument("phi_main", nargs='?', type=float, default=0.3719)
#     args = parser.parse_args()
#     # enttype = args.enttype
#     ent_main=args.ent_main
#     ent_sec=args.ent_sec
#     out_dir=args.out_dir
#     # tau_sec=args.tau_sec
#     phi_jet_norm=args.phi_jet_norm
#     phi_global=args.phi_global
#     phi_main=args.phi_main