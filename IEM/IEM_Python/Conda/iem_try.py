# for the glory of God
from collections import Counter
from cantera import Solution
from cantera._cantera import ConstPressureReactor, ReactorNet
from matplotlib.pyplot import plot, show
from numpy.ma import array, arange


def main():
# Create Reactor Networks
    m, tpArray, rArray, rn = createReactorNetwork()

# Create t vector
    milliseconds = 0.001;
    dt = 0.0005*milliseconds;
    tau_mix = 0.05*milliseconds;
    end_time = 1*milliseconds;
    t = arange(0, end_time, dt)

# Iterative Solution
    T_1 = [None]*len(t)
    for i in range(0,len(t)-1):
        T_1[i] = tpArray[1].T;
        rn.advance(t[i]);
        iem(m, tpArray, rArray, rn, dt, 1/tau_mix);

# Return values
    return t, T_1;

def iem(m, tpArray, rArray, rn, dt, omega):
    # Constant k:
    C_phi = 1;
    k = -C_phi * omega * 0.5 * dt;
    # Calculate average
    m_total_r = 1/sum(m)
    H = 0
    mX_total = Counter(tpArray[0].mass_fraction_dict().keys());
    for i in range(0,len(tpArray)-1):

        X = tpArray[i].mass_fraction_dict();
        for keys in X:
            X[keys] *= m[i];

        H += tpArray[i].enthalpy_mass * m[i]
        mX_total += Counter(X)

    for keys in mX_total:
        mX_total[keys] *= m_total_r;

    H *= m_total_r;

    for i in range(0, len(tpArray)-1):
        h = tpArray[i].enthalpy_mass;
        h = h + (k * (h - H));
        Y = tpArray[i].mass_fraction_dict();

        for keys in Y:
            Y[keys] = Y[keys] + (k * (Y[keys] - mX_total[keys]))

        tpArray[i].HPY = [h, tpArray[i].P, Y]
        rArray[i].syncState();

    rn.reinitialize();

def createReactorNetwork():
    m = array([6, 2, 3]);
    s1 = Solution('h2o2.cti');
    s2 = Solution('h2o2.cti');
    s3 = Solution('h2o2.cti');

    s1.TPY = 300, 101325, 'H2:0.75, O2:0.25';
    s2.TPY = [300, 101325, 'O2:0.75, H2:0.25'];
    s3.TPY = [500, 101325, 'H2:0.5, O2:0.5'];
    s3.equilibrate();
    tpArray = array([s1, s2, s3])
    r1 = ConstPressureReactor(s1);
    r2 = ConstPressureReactor(s2);
    r3 = ConstPressureReactor(s3);
    rArray = array([r1, r2, r3])
    
    rn = ReactorNet([r1, r2, r3]);
    r1.syncState(); 
    r2.syncState(); 
    r3.syncState();
    rn.reinitialize()
    return m, tpArray, rArray, rn;

t, T = main()
print("Done")
plot(t, T)
show()