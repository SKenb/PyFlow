from PyFlow.Measurements import filterDataByKey, addTimeOffset, deriveNewKeyFromData, loadFromFile, addRelativeTime, mergeDataSetsOnTimeKey, plotData, removeNaNEntries, reduceDataToTimeRange, saveToPickle
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


def loadInputDataFrom_Jun04(): return loadInputDataFrom("C:\\Users\\sebkno\\SynologyDrive\\PhD\\11_Hydrogenation\\B0_BuchwaldHartwigReaction\\Data\\250604_BH_0.csv")
def loadInputDataFrom_Jun10(): return loadInputDataFrom("C:\\Users\\sebkno\\SynologyDrive\\PhD\\11_Hydrogenation\\B0_BuchwaldHartwigReaction\\Data\\250610_BH_1.csv")
def loadInputDataFrom_Jun12(): return loadInputDataFrom("C:\\Users\\sebkno\\SynologyDrive\\PhD\\11_Hydrogenation\\B0_BuchwaldHartwigReaction\\Data\\250612_BH_2.csv")
def loadInputDataFrom_Jun16(): return loadInputDataFrom("C:\\Users\\sebkno\\SynologyDrive\\PhD\\11_Hydrogenation\\B0_BuchwaldHartwigReaction\\Data\\250616_BH_3.csv")
def loadInputDataFrom_Jul01(): return loadInputDataFrom("C:\\Users\\sebkno\\SynologyDrive\\PhD\\11_Hydrogenation\\B0_BuchwaldHartwigReaction\\Data\\250701_BH_4.csv")
def loadInputDataFrom_Jul14(): return loadInputDataFrom("C:\\Users\\sebkno\\SynologyDrive\\PhD\\11_Hydrogenation\\B0_BuchwaldHartwigReaction\\Data\\250714_BH_5.csv")
def loadInputDataFrom_Oct21(): return loadInputDataFrom("C:\\Users\\sebkno\\SynologyDrive\\PhD\\11_Hydrogenation\\B0_BuchwaldHartwigReaction\\Data\\251021_BH_6.csv")

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
    irMeasurements = loadIRMeasurements(timeOffsetInSeconds=0)
    uhplcMeasurements = loadUHPLCMeasurements(timeOffsetInSeconds=0)

    # filter for 04.06.2025
    filterByDay = lambda day, month: lambda d: d > datetime(2025, month, day, 0, 0, 0) and d < datetime(2025, month, day + 1, 0, 0, 0)

    for dayInfos in [
        {
            "day": 4, "month": 6, "irTimeOffsetInSeconds": -750, 
            "inputDataLoader": loadInputDataFrom_Jun04,
            "xLimits": [0, 3.2], "yLimits": None
        },
        {
            "day": 10, "month": 6, "irTimeOffsetInSeconds": 1094, 
            "inputDataLoader": loadInputDataFrom_Jun10,
            "xLimits": [0, 6], "yLimits": [0, .5]
        },
        {
            "day": 12, "month": 6, "irTimeOffsetInSeconds": -90, 
            "inputDataLoader": loadInputDataFrom_Jun12,
            "xLimits": [0, 4.8], "yLimits": None
        },
        {
            "day": 16, "month": 6, "irTimeOffsetInSeconds": -660, 
            "inputDataLoader": loadInputDataFrom_Jun16,
            "xLimits": [.3, 4.1], "yLimits": [0, .5]
        },
        {
            "day": 1, "month": 7, "irTimeOffsetInSeconds": -415, 
            "inputDataLoader": loadInputDataFrom_Jul01,
            "xLimits": None, "yLimits": None
        },
        {
            "day": 14, "month": 7, "irTimeOffsetInSeconds": 0, 
            "inputDataLoader": loadInputDataFrom_Jul14,
            "xLimits": [2.8, 8], "yLimits": None
        },
        {
            "day": 21, "month": 10, "irTimeOffsetInSeconds": -1000, 
            "inputDataLoader": loadInputDataFrom_Oct21,
            "xLimits": [0.13, 6.8], "yLimits": None
        },
    ]:

        filteredIRMeasurements = filterDataByKey(irMeasurements, "Time", filterByDay(dayInfos["day"], dayInfos["month"]))
        filteredIRMeasurements = addRelativeTime(filteredIRMeasurements, "Time", timeOffsetInSeconds=dayInfos["irTimeOffsetInSeconds"])

        filteredUHPLCMeasurements = filterDataByKey(uhplcMeasurements, "Time", filterByDay(dayInfos["day"], dayInfos["month"]))
        filteredUHPLCMeasurements = addRelativeTime(filteredUHPLCMeasurements, "Time", timeOffsetInSeconds=0)   

        mergedMeasurements = mergeDataSetsOnTimeKey(filteredUHPLCMeasurements, filteredIRMeasurements, "RelTime", newCommonTimeVector=filteredIRMeasurements["RelTime"])

        inputData = dayInfos["inputDataLoader"]()
        deltaTimeInSeconds = (filteredUHPLCMeasurements["Time"][0] - inputData["Time"][0]).total_seconds()
        inputData = addRelativeTime(inputData, "Time", "RelTime", timeOffsetInSeconds=-deltaTimeInSeconds)
        inputData = deriveInputConcentrations(inputData)

        mergedData = mergeDataSetsOnTimeKey(mergedMeasurements, inputData, "RelTime", newCommonTimeVector=mergedMeasurements["RelTime"]) 



    #    plotData(
    #        filteredIRMeasurements,
    #        showPlot=False,
    #        plotDataDictArray=[
    #            {
    #                "xKey": "RelTime",
    #                "xKeyFormatter": lambda tInSeconds: tInSeconds/3600,
    #                "yKeys": ["Meas Bromonitrobenzene (IR)", "Meas Thiophene (IR)", "Meas Product (IR)"],
    #                "title": "IR - Measurement Data over Time",
    #                "xLabel": "Time in hours",
    #                "yLabel": "Concentration in mol/L",
    #            }
    #        ]
    #    )

    #    plotData(
    #        filteredUHPLCMeasurements,
    #        showPlot=False,
    #        plotDataDictArray=[
    #            {
    #                "xKey": "RelTime",
    #                "xKeyFormatter": lambda tInSeconds: tInSeconds/3600,
    #                "yKeys": ["Meas Bromonitrobenzene (UHPLC)", "Meas Thiophene (UHPLC)", "Meas Product (UHPLC)"],
    #                "title": "UHPLC - Measurement Data over Time",
    #                "xLabel": "Time in hours",
    #                "yLabel": "Concentration in mol/L",
    #            }
    #        ]
    #    )
        

    #    plotData(
    #        inputData,
    #        showPlot=False,
    #        plotDataDictArray=[
    #            {
    #                "xKey": "RelTime",
    #                "xKeyFormatter": lambda tInSeconds: tInSeconds/3600,
    #                "yKeys": ["2-bromonitrobenze", "Thiophene", "DBU", "Xantphos", "Solvent"],
    #                "title": "Input Data over Time",
    #                "xLabel": "Time in hours",
    #                "yLabel": "Flow rate in mL/min",
    #            },
    #            {
    #                "xKey": "RelTime",
    #                "xKeyFormatter": lambda tInSeconds: tInSeconds/3600,
    #                "yKeys": ["Temp"],
    #                "title": "Input Temperature over Time",
    #                "xLabel": "Time in hours",
    #                "yLabel": "Temperature in °C",
    #            }
    #        ]
    #    )


        plotData(
            mergedData,
            showPlot=True,
            plotDataDictArray=[
                {
                    "xKey": "RelTime",
                    "xKeyFormatter": lambda tInSeconds: tInSeconds/3600,
                    "yKeys": ["Meas Bromonitrobenzene (IR)", "Meas Thiophene (IR)", "Meas Product (IR)", "Meas Bromonitrobenzene (UHPLC)", "Meas Thiophene (UHPLC)", "Meas Product (UHPLC)"],
                    "title": "Merged measurements over Time",
                    "xLabel": "Time in hours",
                    "yLabel": "Concentration in mol/L",
                    "xLimit": dayInfos["xLimits"],
                    "yLimit": dayInfos["yLimits"],
                },
                {
                    "xKey": "RelTime",
                    "xKeyFormatter": lambda tInSeconds: tInSeconds/3600,
                    "yKeys": ["Conc 2-bromonitrobenze", "Conc Thiophene", "Conc DBU", "Conc Xantphos", "Meas Bromonitrobenzene (IR)", "Meas Thiophene (IR)", "Meas Product (IR)", "Meas Bromonitrobenzene (UHPLC)", "Meas Thiophene (UHPLC)", "Meas Product (UHPLC)"],
                    "title": "Input Concentrations over Time",
                    "xLabel": "Time in hours",
                    "yLabel": "Concentration in mol/L",
                    "xLimit": dayInfos["xLimits"],
                    "yLimit": dayInfos["yLimits"],
                },
                {
                    "xKey": "RelTime",
                    "xKeyFormatter": lambda tInSeconds: tInSeconds/3600,
                    "yKeys": ["Conc 2-bromonitrobenze", "Meas Bromonitrobenzene (IR)", "Meas Bromonitrobenzene (UHPLC)"],
                    "title": "Input Concentrations over Time",
                    "xLabel": "Time in hours",
                    "yLabel": "Concentration in mol/L",
                    "xLimit": dayInfos["xLimits"],
                    "yLimit": dayInfos["yLimits"],
                }
            ]
        )

        reducedData = reduceDataToTimeRange(mergedData, "RelTime", dayInfos["xLimits"][0]*3600, dayInfos["xLimits"][1]*3600) if dayInfos["xLimits"] is not None else mergedData
        saveToPickle(reducedData, f"C:\\Users\\sebkno\\SynologyDrive\\PhD\\11_Hydrogenation\\B0_BuchwaldHartwigReaction\\PreparedData\\PreparedMeasurements_{dayInfos['day']}_{dayInfos['month']}.pkl")

if __name__ == "__main__":
    main()