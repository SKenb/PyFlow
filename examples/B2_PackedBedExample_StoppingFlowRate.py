import numpy as np

from PyFlow.PackedBedFlowReactor import PackedBedFlowReactor

def reactionNetworkCallbackWithPackedBed(Cs, index, theta, packedBedQuality):
    # C1, C2 + Packed Bed Quality
    # C1 + Packed Bed -> C2
    assert Cs.shape[0] == 2

    theta_K = theta + 273.15
    rate = Cs[0, :] * packedBedQuality * 10 * np.exp(-10e3 / theta_K / 8.314)
    
    if index == 0: return -.2*rate # Consumption of C1
    if index == 1: return rate # Formation of C2

    if index == -1: return -.2*rate  # Packed bed quality decrease

    raise ValueError("Invalid species index")


def main():
    reactor = PackedBedFlowReactor(length=10, diameter=.8*1e-3)
    reactor.setReactionNetworkCallback(reactionNetworkCallbackWithPackedBed)
    print(reactor)


    timeVec = np.linspace(0, 30, 100)
    Cin = np.zeros((2, len(timeVec)))
    flowRate = np.ones_like(timeVec)
    temperature = 150

    Cin[0, 10:30] = 1.0  # Species 1 pulse from t=1s to t=3s
    flowRate[timeVec > 10] = 0 

    reactor.simulate(Cin=Cin, timeVec=timeVec, flowRate=flowRate, temperature=temperature, isVolumetricFlowRate=False)
    reactor.plot(Cin=Cin, CinTimeVec=timeVec, showPlot=True)


# python -m examples.A0_SimpleExample 
if __name__ == "__main__":
    main()