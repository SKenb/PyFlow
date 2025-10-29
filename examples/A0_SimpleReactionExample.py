import numpy as np

from PyFlow.FlowReactor import FlowReactor

def reactionNetworkCallback(Cs, index, theta):
    # C1, C2, C3
    # C1 + C2 --> C3
    assert Cs.shape[0] == 3

    theta_K = theta + 273.15
    rate = Cs[0, :] * Cs[1, :] * 10 * np.exp(-10e3 / theta_K / 8.314)
    return rate if index == 2 else -1 * rate

def main():
    reactor = FlowReactor(length=10, diameter=.8*1e-3)
    reactor.setReactionNetworkCallback(reactionNetworkCallback)
    print(reactor)


    timeVec = np.linspace(0, 20, 100)
    Cin = np.zeros((3, len(timeVec)))
    flowRate = 1
    temperature = 150

    Cin[0, 10:30] = 1.0  # Species 1 pulse from t=1s to t=3s
    Cin[1, :10] = 0.5  # Species 2 pulse from t=5s to t=7s

    reactor.simulate(Cin=Cin, timeVec=timeVec, flowRate=flowRate, temperature=temperature, isVolumetricFlowRate=False)
    reactor.plot(Cin=Cin, CinTimeVec=timeVec, showPlot=True)


# python -m examples.A0_SimpleExample 
if __name__ == "__main__":
    main()