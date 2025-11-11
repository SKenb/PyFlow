import os
import numpy as np
import matplotlib.pyplot as plt

from PyFlow.Measurements import loadFromPickle, plotData
from PyFlow.FlowReactor import FlowReactor

def main():
    pickleFilePath = "C:\\Users\\sebkno\\SynologyDrive\\PhD\\11_Hydrogenation\\B0_BuchwaldHartwigReaction\\PreparedData\\"
    pickleFiles = os.listdir(pickleFilePath)
    dataSets = {}

    flowReactor = FlowReactor(length=10, diameter=0.8*1e-3)
    flowReactor.setVolumeIn_mL(4.67)

    print(flowReactor)

    for pf in pickleFiles:
        if pf.endswith(".pkl"):
            fullPath = os.path.join(pickleFilePath, pf)
            data = loadFromPickle(fullPath)
            dataSets[pf] = data
            print(f"Loaded data from {pf} with keys: {list(data.keys())} and {len(next(iter(data.values())))} entries.")

            timeVec = np.array([t - data["RelTime"][0] for t in data["RelTime"]])
            temperature = np.array(data["Temp"])
            totalFlowRate = np.array(data["TotalFlowRate"])

            Cin = np.zeros((7, len(data["RelTime"])))

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


            # Plot data
            plotData(data,
                showPlot=True,
                plotDataDictArray=[
                    {
                        "xKey": "RelTime",
                        "xKeyFormatter": lambda tInSeconds: tInSeconds/3600,
                        "yKeys": [key for key in data.keys() if key.startswith("Meas")],
                        "title": f"Measured Concentrations over Time from {pf}",
                        "xLabel": "Time in hours",
                        "yLabel": "Concentration in mol/L",
                    }
                ]
            )

            # Plot outputs
            figure = plt.figure(figsize=(10, 6))
            plt.plot(time, Cout[0, :], label="Cout 2-bromonitrobenze")
            plt.plot(timeVec, data["Meas Bromonitrobenzene (UHPLC)"], label="Measured 2-bromonitrobenze (UHPLC)")
            plt.plot(timeVec, data["Meas Bromonitrobenzene (IR)"], label="Measured 2-bromonitrobenze (IR)")
            
            plt.show()


if __name__ == "__main__":
    main()