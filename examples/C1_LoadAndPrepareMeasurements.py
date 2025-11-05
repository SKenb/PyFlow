from PyFlow.Measurements import deriveNewKeyFromData, loadFromFile, addRelativeTime, mergeDataSetsOnTimeKey, plotData, removeNaNEntries, reduceDataToTimeRange, saveToPickle
from datetime import datetime

import matplotlib.pyplot as plt

def loadIRMeasurements(timeOffsetInSeconds=0):
    filePath = "C:\\Users\\sebkno\\SynologyDrive\\PhD\\11_Hydrogenation\\B0_BuchwaldHartwigReaction\\Data\\Buchwald_Hartwig_all_data_IR_251028.xlsx"
    measurementData = loadFromFile(
        filePath, 
        startDataRow=5, 
        headerRenameToMapping={
            "Sample": "Time",
            "Predicted bromonitrobenzene": "Meas Bromonitrobenzene (IR)",
            "Predicted thiophene": "Meas Thiophene (IR)",
            "Predicted product": "Meas Product (IR)",
        },
        headerParseFunctionMapping={
            "Time": lambda d_: datetime.strptime(d_.replace('_ ', '_').strip('_').replace('.spc#1', '').replace("250", "2025"), "%Y%m%d_%H_%M_%S"),
            "Meas Bromonitrobenzene (IR)": float,
            "Meas Thiophene (IR)": float,
            "Meas Product (IR)": float
        }
    )
    
    measurementData = addRelativeTime(measurementData, "Time", "RelTime", timeOffsetInSeconds=timeOffsetInSeconds)
    return measurementData

def loadUHPLCMeasurements(timeOffsetInSeconds=0):
    filePath = "C:\\Users\\sebkno\\SynologyDrive\\PhD\\11_Hydrogenation\\B0_BuchwaldHartwigReaction\\Data\\Buchwald_Hartwig_all_UHPLC_data.xlsx"
    measurementData = loadFromFile(
        filePath, 
        startDataRow=5, 
        headerRenameToMapping={
            "Timestamp": "Time",
            "Predicted bromonitrobenzene": "Meas Bromonitrobenzene (UHPLC)",
            "Predicted thiophene": "Meas Thiophene (UHPLC)",
            "Predicted product": "Meas Product (UHPLC)",
        },
        headerParseFunctionMapping={
            "Time": lambda d_: datetime.strptime(d_, "%d.%m.%Y %H:%M:%S"),
            "Meas Bromonitrobenzene (UHPLC)": float,
            "Meas Thiophene (UHPLC)": float,
            "Meas Product (UHPLC)": float
        }
    )
    
    measurementData = addRelativeTime(measurementData, "Time", "RelTime", timeOffsetInSeconds=timeOffsetInSeconds)
    measurementData = removeNaNEntries(measurementData, ["RelTime"])
    return measurementData

def loadInputDataFrom_Jun06(): return loadInputDataFrom("C:\\Users\\sebkno\\SynologyDrive\\PhD\\11_Hydrogenation\\B0_BuchwaldHartwigReaction\\Data\\250604_BH_0.csv")
def loadInputDataFrom_Jun10(): return loadInputDataFrom("C:\\Users\\sebkno\\SynologyDrive\\PhD\\11_Hydrogenation\\B0_BuchwaldHartwigReaction\\Data\\250610_BH_1.csv")
def loadInputDataFrom_Jun12(): return loadInputDataFrom("C:\\Users\\sebkno\\SynologyDrive\\PhD\\11_Hydrogenation\\B0_BuchwaldHartwigReaction\\Data\\250612_BH_2.csv")
def loadInputDataFrom_Jun16(): return loadInputDataFrom("C:\\Users\\sebkno\\SynologyDrive\\PhD\\11_Hydrogenation\\B0_BuchwaldHartwigReaction\\Data\\250616_BH_3.csv")
def loadInputDataFrom_Jul01(): return loadInputDataFrom("C:\\Users\\sebkno\\SynologyDrive\\PhD\\11_Hydrogenation\\B0_BuchwaldHartwigReaction\\Data\\250701_BH_4.csv")
def loadInputDataFrom_Jul14(): return loadInputDataFrom("C:\\Users\\sebkno\\SynologyDrive\\PhD\\11_Hydrogenation\\B0_BuchwaldHartwigReaction\\Data\\250714_BH_5.csv")
def loadInputDataFrom_Jul21(): return loadInputDataFrom("C:\\Users\\sebkno\\SynologyDrive\\PhD\\11_Hydrogenation\\B0_BuchwaldHartwigReaction\\Data\\251021_BH_6.csv")

def loadInputDataFrom(filePath):
    inputData = loadFromFile(
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

    inputData = addRelativeTime(inputData, "Time", "RelTime", timeOffsetInSeconds=0)

    return inputData

def deriveInputConcentrations(inputData):
    inputData = deriveNewKeyFromData(inputData, "TotalFlowRate", lambda row: row["2-bromonitrobenze"] + row["Thiophene"] + row["DBU"] + row["Xantphos"] + row["Solvent"])
    inputData = deriveNewKeyFromData(inputData, "Conc 2-bromonitrobenze", lambda row: .7 * row["2-bromonitrobenze"] / row["TotalFlowRate"] if row["TotalFlowRate"] > 1e-3 else 0)
    inputData = deriveNewKeyFromData(inputData, "Conc Thiophene", lambda row: .84 * row["Thiophene"] / row["TotalFlowRate"] if row["TotalFlowRate"] > 1e-3 else 0)
    inputData = deriveNewKeyFromData(inputData, "Conc DBU", lambda row: .98 * row["DBU"] / row["TotalFlowRate"] if row["TotalFlowRate"] > 1e-3 else 0)
    inputData = deriveNewKeyFromData(inputData, "Conc Xantphos", lambda row: 0.275 * row["Xantphos"] / row["TotalFlowRate"] if row["TotalFlowRate"] > 1e-3 else 0)
    return inputData


def main():
    irMeasurements = loadIRMeasurements(timeOffsetInSeconds=-750)
    uhplcMeasurements = loadUHPLCMeasurements(timeOffsetInSeconds=0)  # Adjust time offset as needed

    mergedMeasurements = mergeDataSetsOnTimeKey(irMeasurements, uhplcMeasurements, "RelTime", newCommonTimeVector=irMeasurements["RelTime"])
    mergedMeasurements["Time"] = uhplcMeasurements["Time"]  # Keep original time from UHPLC measurements as timestamp is absolute

    xlimits = []
    for inputDataLoader in [loadInputDataFrom_Jun06, loadInputDataFrom_Jun10, loadInputDataFrom_Jun12, loadInputDataFrom_Jun16, loadInputDataFrom_Jul01, loadInputDataFrom_Jul14, loadInputDataFrom_Jul21]: #loadInputDataFrom_Jun06, loadInputDataFrom_Jun10,
        inputData = inputDataLoader()
        deltaTimeInSeconds = (mergedMeasurements["Time"][0] - inputData["Time"][0]).total_seconds()
        inputData = addRelativeTime(inputData, "Time", "RelTime", timeOffsetInSeconds=-deltaTimeInSeconds)
        inputData = deriveInputConcentrations(inputData)

        print(f"Merging input data from {inputDataLoader.__name__} with time offset {deltaTimeInSeconds} seconds.")
        mergedMeasurements = mergeDataSetsOnTimeKey(mergedMeasurements, inputData, "RelTime", newCommonTimeVector=mergedMeasurements["RelTime"]) 
        mergedMeasurements["Time"] = uhplcMeasurements["Time"]


    xlimits = [[0, 4], [143, 148.1], [192, 194], [288, 290], [647.5, 648.5], [958, 973]]
    for xlimits in xlimits:
        plotData(mergedMeasurements, 
            showPlot=False,
            plotDataDictArray=[
                {
                    "xKey": "RelTime",
                    "xKeyFormatter": lambda tInSeconds: tInSeconds/3600,
                    "yKeys": ["2-bromonitrobenze", "Thiophene", "DBU", "Xantphos", "Solvent", "TotalFlowRate"],
                    "title": "Flow rate over Time",
                    "xLabel": "Time in hours",
                    "yLabel": "Flow rate in mL/min",
                    "xLimit": xlimits,
                },
                {
                    "xKey": "RelTime",
                    "xKeyFormatter": lambda tInSeconds: tInSeconds/3600,
                    "yKeys": ["Meas Bromonitrobenzene (IR)", "Meas Thiophene (IR)", "Meas Product (IR)", "Meas Bromonitrobenzene (UHPLC)", "Meas Thiophene (UHPLC)", "Meas Product (UHPLC)"],
                    "title": "Measured Concentrations over Time",
                    "xLabel": "Time in hours",
                    "yLabel": "Concentration in mol/L",
                    "xLimit": xlimits,
                },
                {
                    "xKey": "RelTime",
                    "xKeyFormatter": lambda tInSeconds: tInSeconds/3600,
                    "yKeys": ["Conc 2-bromonitrobenze", "Conc Thiophene", "Conc DBU", "Conc Xantphos"],
                    "title": "Input Concentrations over Time",
                    "xLabel": "Time in hours",
                    "yLabel": "Concentration in mol/L",
                    "xLimit": xlimits,
                },
                {
                    "xKey": "RelTime",
                    "xKeyFormatter": lambda tInSeconds: tInSeconds/3600,
                    "yKeys": ["Conc 2-bromonitrobenze", "Meas Bromonitrobenzene (IR)", "Meas Bromonitrobenzene (UHPLC)"],
                    "title": "Bromonitrobenzene over Time",
                    "xLabel": "Time in hours",
                    "yLabel": "Concentration in mol/L",
                    "xLimit": xlimits,
                }
            ]
        )

        reducedData = reduceDataToTimeRange(mergedMeasurements, "RelTime", xlimits[0], xlimits[1])
        saveToPickle(reducedData, f"C:\\Users\\sebkno\\SynologyDrive\\PhD\\11_Hydrogenation\\B0_BuchwaldHartwigReaction\\PreparedData\\PreparedMeasurements_{int(xlimits[0])}_{int(xlimits[1])}.pkl")

    plt.show()


if __name__ == "__main__":
    main()