from PyFlow.Measurements import deriveNewKeyFromData, loadFromFile, addRelativeTime, mergeDataSetsOnTimeKey, plotData, removeNaNEntries
from datetime import datetime

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

def loadInputDataFrom_Jun06():
    filePath = "C:\\Users\\sebkno\\SynologyDrive\\PhD\\11_Hydrogenation\\B0_BuchwaldHartwigReaction\\Data\\250604_BH_0.csv"

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


def main():
    irMeasurements = loadIRMeasurements(timeOffsetInSeconds=-750)
    uhplcMeasurements = loadUHPLCMeasurements(timeOffsetInSeconds=0)  # Adjust time offset as needed

    mergedMeasurements = mergeDataSetsOnTimeKey(irMeasurements, uhplcMeasurements, "RelTime", newCommonTimeVector=irMeasurements["RelTime"])
    mergedMeasurements["Time"] = uhplcMeasurements["Time"]  # Keep original time from UHPLC measurements as timestamp is absolute


    inputData_Jun06 = loadInputDataFrom_Jun06()
    deltaTimeInSeconds = (mergedMeasurements["Time"][0] - inputData_Jun06["Time"][0]).total_seconds()
    inputData_Jun06 = addRelativeTime(inputData_Jun06, "Time", "RelTime", timeOffsetInSeconds=-deltaTimeInSeconds)

    inputData_Jun06 = deriveNewKeyFromData(inputData_Jun06, "TotalFlowRate", lambda row: row["2-bromonitrobenze"] + row["Thiophene"] + row["DBU"] + row["Xantphos"] + row["Solvent"])
    inputData_Jun06 = deriveNewKeyFromData(inputData_Jun06, "Conc 2-bromonitrobenze", lambda row: .7 * row["2-bromonitrobenze"] / row["TotalFlowRate"] if row["TotalFlowRate"] > 1e-3 else 0)
    inputData_Jun06 = deriveNewKeyFromData(inputData_Jun06, "Conc Thiophene", lambda row: .84 * row["Thiophene"] / row["TotalFlowRate"] if row["TotalFlowRate"] > 1e-3 else 0)
    inputData_Jun06 = deriveNewKeyFromData(inputData_Jun06, "Conc DBU", lambda row: .98 * row["DBU"] / row["TotalFlowRate"] if row["TotalFlowRate"] > 1e-3 else 0)
    inputData_Jun06 = deriveNewKeyFromData(inputData_Jun06, "Conc Xantphos", lambda row: 0.275 * row["Xantphos"] / row["TotalFlowRate"] if row["TotalFlowRate"] > 1e-3 else 0)

    mergedMeasurements = mergeDataSetsOnTimeKey(mergedMeasurements, inputData_Jun06, "RelTime", newCommonTimeVector=inputData_Jun06["RelTime"])

    plotData(mergedMeasurements, 
        [
            {
                "xKey": "RelTime",
                "xKeyFormatter": lambda tInSeconds: tInSeconds/3600,
                "yKeys": ["2-bromonitrobenze", "Thiophene", "DBU", "Xantphos", "Solvent", "TotalFlowRate"],
                "title": "Flow rate over Time",
                "xLabel": "Time in hours",
                "yLabel": "Flow rate in mL/min",
                "xLimit": [0, 4],
            },
            {
                "xKey": "RelTime",
                "xKeyFormatter": lambda tInSeconds: tInSeconds/3600,
                "yKeys": ["Meas Bromonitrobenzene (IR)", "Meas Thiophene (IR)", "Meas Product (IR)", "Meas Bromonitrobenzene (UHPLC)", "Meas Thiophene (UHPLC)", "Meas Product (UHPLC)"],
                "title": "Measured Concentrations over Time",
                "xLabel": "Time in hours",
                "yLabel": "Concentration in mol/L",
                "xLimit": [0, 4],
            },
            {
                "xKey": "RelTime",
                "xKeyFormatter": lambda tInSeconds: tInSeconds/3600,
                "yKeys": ["Conc 2-bromonitrobenze", "Conc Thiophene", "Conc DBU", "Conc Xantphos"],
                "title": "Input Concentrations over Time",
                "xLabel": "Time in hours",
                "yLabel": "Concentration in mol/L",
                "xLimit": [0, 4],
            },
            {
                "xKey": "RelTime",
                "xKeyFormatter": lambda tInSeconds: tInSeconds/3600,
                "yKeys": ["Conc 2-bromonitrobenze", "Meas Bromonitrobenzene (IR)", "Meas Bromonitrobenzene (UHPLC)"],
                "title": "Bromonitrobenzene over Time",
                "xLabel": "Time in hours",
                "yLabel": "Concentration in mol/L",
                "xLimit": [0, 4],
            }
        ]
    )


if __name__ == "__main__":
    main()