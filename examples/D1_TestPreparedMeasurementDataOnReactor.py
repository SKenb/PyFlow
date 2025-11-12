import os
import numpy as np
import matplotlib.pyplot as plt

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


def main():
    pickleFilePath = "C:\\Users\\sebkno\\SynologyDrive\\PhD\\11_Hydrogenation\\B0_BuchwaldHartwigReaction\\PreparedData\\"
    pickleFiles = os.listdir(pickleFilePath)
    dataSets = {}

    flowReactor = FlowReactor(
        length=10, diameter=0.8*1e-3, 
        reactionNetworkCallback=buchwaldHartwigReaction,
        setSpaceSamples=40
    )
    flowReactor.setVolumeIn_mL(4.67)

    print(flowReactor)

    for pf in pickleFiles:
        if pf.endswith(".pkl"):
            fullPath = os.path.join(pickleFilePath, pf)
            data = loadFromPickle(fullPath)
            dataSets[pf] = data
            print(f"Loaded data from {pf} with keys: {list(data.keys())} and {len(next(iter(data.values())))} entries.")

            data = reduceDataToTimeRange(data, "RelTime", 0, 10000)

            timeVec = np.array([t - data["RelTime"][0] for t in data["RelTime"]])
            temperature = np.array(data["Temp"])
            totalFlowRate = np.array(data["TotalFlowRate"])

            Cin = np.zeros((8, len(data["RelTime"])))

            Cin[0, :] = np.array(data["Conc 2-bromonitrobenze"])
            Cin[1, :] = np.array(data["Conc Thiophene"])
            Cin[2, :] = np.array(data["Conc DBU"])
            Cin[3, :] = np.array(data["Conc Xantphos"])

            # Test on reactor
            time, Cout, Cspatial = flowReactor.simulate(
                timeVec, 
                Cin, 
                totalFlowRate, 
                temperature,
            )


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
            figure = plt.figure(figsize=(10, 6))

            toPlotData = [
                {"CIdx": 0, "label": "Cout 2-bromonitrobenze", "UHPLC_Key": "Meas Bromonitrobenzene (UHPLC)" , "IR_Key": "Meas Bromonitrobenzene (IR)"},
                {"CIdx": 1, "label": "Cout Thiophene", "UHPLC_Key": "Meas Thiophene (UHPLC)" , "IR_Key": "Meas Thiophene (IR)"},
                {"CIdx": 4, "label": "Cout Product", "UHPLC_Key": "Meas Product (UHPLC)" , "IR_Key": "Meas Product (IR)"},
            ]

            for idx, plotInfo in enumerate(toPlotData):

                ax = figure.add_subplot(len(toPlotData), 1, idx + 1)

                ax.plot(timeVec, Cin[plotInfo["CIdx"], :], '--', label=f"Input: {plotInfo['label']}")
                ax.plot(time, Cout[plotInfo["CIdx"], :], label=plotInfo["label"])
                ax.plot(timeVec, data[plotInfo["UHPLC_Key"]], label=f"Measured {plotInfo['label']} (UHPLC)")
                ax.plot(timeVec, data[plotInfo["IR_Key"]], label=f"Measured {plotInfo['label']} (IR)")
                ax.set_title(f"Outlet Concentration vs Measured Data from {pf}")
                ax.set_xlabel("Time in seconds")
                ax.set_ylabel("Concentration in mol/L")
                ax.legend()
            
            flowReactor.plot(showPlot=False)

            plt.show()

            print(f"Simulation for {pf} completed.")


if __name__ == "__main__":
    main()