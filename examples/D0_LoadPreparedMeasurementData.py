import os

from PyFlow.Measurements import loadFromPickle, plotData

def main():
    pickleFilePath = "C:\\Users\\sebkno\\SynologyDrive\\PhD\\11_Hydrogenation\\B0_BuchwaldHartwigReaction\\PreparedData\\"
    pickleFiles = os.listdir(pickleFilePath)
    dataSets = {}

    for pf in pickleFiles:
        if pf.endswith(".pkl"):
            fullPath = os.path.join(pickleFilePath, pf)
            data = loadFromPickle(fullPath)
            dataSets[pf] = data
            print(f"Loaded data from {pf} with keys: {list(data.keys())} and {len(next(iter(data.values())))} entries.")

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

if __name__ == "__main__":
    main()