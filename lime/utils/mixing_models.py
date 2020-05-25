def iem(m, tpArray, rArray, rn, dt, omega):
    # Constant k:
    C_phi = 2
    k = -C_phi * omega * 0.5 * dt
    # Calculate average:
    m_total_r = 1 / sum(m)
    M_species_total = sum([m[i] * tpArray[i].Y for i in range(0, len(tpArray))])
    H_total = sum([m[i] * tpArray[i].enthalpy_mass for i in range(0, len(tpArray))])
    Y_avg = M_species_total * m_total_r  # Y_species_avg = (M_total_species)/(M_total_system)
    h_avg = H_total * m_total_r  # H_avg is the specific mass-weighted average across all reactors of the total enthalpy.
    # Adjust reactor state:
    for i in range(0, len(tpArray)):
        Y_current = tpArray[i].Y
        Y_new = Y_current + k * (Y_current - Y_avg)
        h = tpArray[i].enthalpy_mass
        h_new = h + k * (h - h_avg)
        tpArray[i].HPY = [h_new, tpArray[i].P, Y_new]
        rArray[i].syncState()
    # Reinitialize reactor network solver:
    rn.reinitialize()
    return Y_avg, h_avg