import csv
import matplotlib.pyplot as plt

from datetime import datetime

def loadFromFile(filePath, delimiter=',', headerRow=1, startDataRow=2, headerRenameToMapping=None, headerParseFunctionMapping=None):
    data = {}

    with open(filePath, newline="") as f:
        lines = f.readlines()

    # Get header and data lines based on given rows (1-indexed)
    header = lines[headerRow - 1].strip().split(delimiter)
    data_lines = lines[startDataRow - 1:]

    if headerRenameToMapping:
        header = [headerRenameToMapping.get(col, col) for col in header]


    reader = csv.DictReader(data_lines, fieldnames=header, delimiter=delimiter)
    for row in reader:
        for key, value in row.items():
            
            parsed_value = value
            if headerParseFunctionMapping and key in headerParseFunctionMapping:
                try:
                    parsed_value = headerParseFunctionMapping[key](value)
                except Exception as e:
                    print(f"Error parsing value: '{value}' for column: '{key}': {e}")
                    parsed_value = float('nan')

            data.setdefault(key.strip(), []).append(parsed_value)

    return data

def plotData(data, plotDataDictArray, showPlot=True, figure=None):
    if figure is None: figure = plt.figure(figsize=(10, 6))

    for index, plotData in enumerate(plotDataDictArray):
        xKey = plotData.get("xKey")
        yKeys = plotData.get("yKeys", [])
        title = plotData.get("title", "Data over " + xKey)
        xLabel = plotData.get("xLabel", xKey)
        yLabel = plotData.get("yLabel", "Data")

        ax = figure.add_subplot(len(plotDataDictArray), 1, index + 1)
        for yKey in yKeys:  ax.plot(data[xKey], data[yKey], label=yKey)
        ax.set_title(title)
        ax.set_xlabel(xLabel)
        ax.set_ylabel(yLabel)
        ax.legend()

    plt.grid(True)

    if showPlot: plt.show()



if __name__ == "__main__":
    # Example usage
    filePath = "C:\\Users\\sebkno\\SynologyDrive\\PhD\\11_Hydrogenation\\B0_BuchwaldHartwigReaction\\Data\\250610_BH_1.csv"

    measurementData = loadFromFile(
        filePath, 
        startDataRow=5, 
        headerRenameToMapping={
            "HPLC_PUMP_1.F_W": "2-bromonitrobenze",
            "HPLC_PUMP_2.F_W": "Thiophene",
            "HPLC_PUMP_3.F_W": "DBU",
            "HPLC_PUMP_4.F_W": "Xantphos",
            "HPLC_PUMP_5.F_W": "Solvent",
            "THERMOSTAT_3.TB": "Temp"
        },
        headerParseFunctionMapping={
            "Time": lambda d_: datetime.strptime(d_, "%d.%m.%Y %H:%M:%S"),
            "2-bromonitrobenze": float,
            "Thiophene": float,
            "DBU": float,
            "Xantphos": float,
            "Solvent": float,
            "Temp": float
        }
    )
    
    print("Loaded Measurement Data:")
    for entry in measurementData: print(entry)

    plotData(
        measurementData,
        plotDataDictArray=[
            {
                "xKey": "Time",
                "yKeys": ["2-bromonitrobenze", "Thiophene", "DBU", "Xantphos", "Solvent"],
                "title": "Concentration vs Time",
                "xLabel": "Time",
                "yLabel": "Concentration (mol/L)"
            },
            {
                "xKey": "Time",
                "yKeys": ["Temp"],
                "title": "Temperature vs Time",
                "xLabel": "Time",
                "yLabel": "Temperature (°C)"
            }
        ],
    )