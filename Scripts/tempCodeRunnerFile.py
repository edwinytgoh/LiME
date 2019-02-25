
    for i in range(len(t)):
        tnow = t[i]
        reactorNetwork.advance(tnow)
        timeHistory.loc[tnow] = reactorNetwork.get_state()
    #print(timeHistory['OH'])
