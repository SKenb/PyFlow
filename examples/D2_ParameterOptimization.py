import os
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import minimize
from PyFlow.Measurements import loadFromPickle, plotData, reduceDataToTimeRange
from PyFlow.FlowReactor import FlowReactor

def buchwaldHartwigReaction(C, i, temperature):
    temperature_K = temperature + 273.15
    reactionRate = lambda Ca, Cb, A, E: A * np.exp(-E/(8.314*temperature_K)) * Ca * Cb
    
    reactionMap = {
        0: { "consumed": {1: 1}, "produced": {} },
        1: { "consumed": {2: 1}, "produced": {} },
        2: { "consumed": {3: 1}, "produced": {} },
        3: { "consumed": {1: .1}, "produced": {4: 1} },
        4: { "consumed": {2: 1}, "produced": {1: 1} },
        5: { "consumed": {3: 1}, "produced": {2: 1} },
        6: { "consumed": {4: 1}, "produced": {3: 1} },
        7: { "consumed": {}, "produced": {4: 1} },
    }
    
    # Reaction 1 - Catalyst (PdL) + 2-Brom --> PdInt1
    R1 = reactionRate(C[3, :], C[0, :], 1.8*1e3, 34_000)
    R1[np.where(np.isnan(R1))] = 0

    # Reaction 2 - PdInt1 + Thiophene --> PdInt2
    R2 = reactionRate(C[4, :], C[1, :], 4.95*1e3, 40_382)
    R2[np.where(np.isnan(R2))] = 0

    # Reaction 3 - PdInt2 + DBU --> PdInt3 + HBr
    R3 = reactionRate(C[5, :], C[2, :], 265.83, 60_000)
    R3[np.where(np.isnan(R3))] = 0

    # Reaction 4 - PdInt3 --> Catalyst (PdL) + Product
    R4 = reactionRate(C[6, :], 1, 160.7, 60_000)
    R4[np.where(np.isnan(R4))] = 0

    RMap = { 1: R1, 2: R2, 3: 0*R3, 4: 0*R4}

    # Total rates of change
    reactionInfos = reactionMap[i]

    R = np.zeros_like(C[0, :])
    for idx, coeff in reactionInfos["produced"].items(): R += coeff * RMap[idx]
    for idx, coeff in reactionInfos["consumed"].items(): R -= coeff * RMap[idx]
    if np.isnan(R).any(): raise ValueError(f"Reaction rate for species index {i} contains NaN values.")
    
    return R

def buchwaldHartwigReactionParams(C, i, temperature, A1, Ea1, A2, Ea2, A3, Ea3, A4, Ea4, k11, k12, k21, k22, k31, k32, k41):
    temperature_K = temperature + 273.15
    reactionRate = lambda Ca, Cb, A, E: A * np.exp(-E/(8.314*temperature_K)) * Ca * Cb
    
    reactionMap = {
        0: { "consumed": {1: k11}, "produced": {} },
        1: { "consumed": {2: k21}, "produced": {} },
        2: { "consumed": {3: k31}, "produced": {} },
        3: { "consumed": {1: k12}, "produced": {4: 1} },
        4: { "consumed": {2: k22}, "produced": {1: 1} },
        5: { "consumed": {3: k32}, "produced": {2: 1} },
        6: { "consumed": {4: k41}, "produced": {3: 1} },
        7: { "consumed": {}, "produced": {4: 1} },
    }
    
    # Reaction 1 - Catalyst (PdL) + 2-Brom --> PdInt1
    R1 = reactionRate(C[3, :], C[0, :], A1, Ea1)
    R1[np.where(np.isnan(R1))] = 0

    # Reaction 2 - PdInt1 + Thiophene --> PdInt2
    R2 = reactionRate(C[4, :], C[1, :], A2, Ea2)
    R2[np.where(np.isnan(R2))] = 0

    # Reaction 3 - PdInt2 + DBU --> PdInt3 + HBr
    R3 = reactionRate(C[5, :], C[2, :], A3, Ea3)
    R3[np.where(np.isnan(R3))] = 0

    # Reaction 4 - PdInt3 --> Catalyst (PdL) + Product
    R4 = reactionRate(C[6, :], 1, A4, Ea4)
    R4[np.where(np.isnan(R4))] = 0

    RMap = { 1: R1, 2: R2, 3: 0*R3, 4: 0*R4}

    # Total rates of change
    reactionInfos = reactionMap[i]

    R = np.zeros_like(C[0, :])
    for idx, coeff in reactionInfos["produced"].items(): R += coeff * RMap[idx]
    for idx, coeff in reactionInfos["consumed"].items(): R -= coeff * RMap[idx]
    if np.isnan(R).any(): raise ValueError(f"Reaction rate for species index {i} contains NaN values.")
    
    return R


def main():
    pickleFilePath = "../Tmp/"
    pickleFiles = os.listdir(pickleFilePath)
    dataSets = {}

    flowReactor = FlowReactor(
        length=10, diameter=0.8*1e-3, 
        reactionNetworkCallback=buchwaldHartwigReaction,
        setSpaceSamples=20
    )
    flowReactor.setVolumeIn_mL(4.67)

    print(flowReactor)

    for pf in pickleFiles:
        if pf.endswith(".pkl"):

            fullPath = os.path.join(pickleFilePath, pf)
            data = loadFromPickle(fullPath)
            dataSets[pf] = data
            print(f"Loaded data from {pf} with keys: {list(data.keys())} and {len(next(iter(data.values())))} entries.")

            data = reduceDataToTimeRange(data, "RelTime", 0, 5000)

            timeVec = np.array([t - data["RelTime"][0] for t in data["RelTime"]])
            temperature = np.array(data["Temp"])
            totalFlowRate = np.array(data["TotalFlowRate"])

            Cin = np.zeros((8, len(data["RelTime"])))

            Cin[0, :] = np.array(data["Conc 2-bromonitrobenze"])
            Cin[1, :] = np.array(data["Conc Thiophene"])
            Cin[2, :] = np.array(data["Conc DBU"])
            Cin[3, :] = np.array(data["Conc Xantphos"])

            ## Test on reactor
            #time, Cout, _ = flowReactor.simulate(
            #    timeVec, 
            #    Cin, 
            #    totalFlowRate, 
            #    temperature,
            #)

            ## Do optimization
            def objectiveFunction(x, flowReactor, timeVec, Cin, totalFlowRate, temperature, measuredData, timeMinMax):

                flowReactor.reactionNetworkCallback = lambda C, i, T: buchwaldHartwigReactionParams(
                    C, i, T, 
                    A1=x[0], Ea1=x[1], A2=x[2], Ea2=x[3], A3=x[4], Ea3=x[5], A4=x[6], Ea4=x[7], 
                    k11=x[8], k12=x[9], k21=x[10], k22=x[11], k31=x[12], k32=x[13], k41=x[14]
                )

                timeSim, CoutSim, _ = flowReactor.simulate(
                    timeVec, 
                    Cin, 
                    totalFlowRate, 
                    temperature,
                )

                error0 = CoutSim[0, :] - np.array(measuredData["Meas Bromonitrobenzene (IR)"])
                error0 = error0[np.where((timeSim >= timeMinMax[0]) & (timeSim <= timeMinMax[1]))]
                error0 = 1e3*np.sum(error0**2)

                error1 = CoutSim[1, :] - np.array(measuredData["Meas Thiophene (IR)"])
                error1 = error1[np.where((timeSim >= timeMinMax[0]) & (timeSim <= timeMinMax[1]))]
                error1 = 1e3*np.sum(error1**2)

                error2 = CoutSim[4, :] - np.array(measuredData["Meas Product (IR)"])
                error2 = error2[np.where((timeSim >= timeMinMax[0]) & (timeSim <= timeMinMax[1]))]
                error2 = 1e3*np.sum(error2**2)

                error = error0 + error1 + error2
                print(f"Simulated with params: {x} --> Objective Error: {error}")
                return error


            ABound = (1, 1e3)
            EaBound = (1e3, 1e5)
            kBound = (0, 300)

            bounds = [ABound, EaBound, ABound, EaBound, ABound, EaBound, ABound, EaBound, kBound, kBound, kBound, kBound, kBound, kBound, kBound, kBound]
            x0 = [1, 3e4, 1, 3e4, 1, 3e4, 1, 3e4, 1, 1, 1, 1, 1, 1, 1, 1]

            objective = lambda x: objectiveFunction(
                x, flowReactor, timeVec, Cin, totalFlowRate, temperature, data,
                timeMinMax=(2000, 2700)
            )
            result = minimize(objective, x0, method='SLSQP', bounds=bounds, constraints=[])

            print("Optimal x:", result.x)
            print("Objective value:", result.fun)

            x = result.x

            
            ## Plot data
            #plotData(data,
            #    showPlot=False,
            #    plotDataDictArray=[
            #        {
            #            "xKey": "RelTime",
            #            "xKeyFormatter": lambda tInSeconds: tInSeconds/3600,
            #            "yKeys": [key for key in data.keys() if key.startswith("Meas")],
            #            "title": f"Measured Concentrations over Time from {pf}",
            #            "xLabel": "Time in hours",
            #            "yLabel": "Concentration in mol/L",
            #        }
            #    ]
            #)

            # Plot outputs

            print(f"Simulation for {pf} completed.")
            exit()


if __name__ == "__main__":
    main()

# Optimal x: [1.56508575e+01 2.99364938e+04 1.97952575e+01 3.00788008e+04
# 1.00000000e+00 3.00000000e+04 1.00000000e+00 3.00000000e+04
# 2.32389337e+01 4.36997442e+00 2.69346357e+00 1.65610402e+00
# 1.00000000e+00 1.00000000e+00 1.00000000e+00 1.00000000e+00]
#Objective value: 1569.4289177386713
# Simulation for PreparedMeasurements_10_6.pkl completed.