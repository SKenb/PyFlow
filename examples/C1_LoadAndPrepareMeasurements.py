from PyFlow.Measurements import loadFromFile, addRelativeTime, mergeDataSetsOnTimeKey, plotData, removeNaNEntries
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


def main():
    irMeasurements = loadIRMeasurements(timeOffsetInSeconds=-800)
    uhplcMeasurements = loadUHPLCMeasurements(timeOffsetInSeconds=0)  # Adjust time offset as needed

    mergedMeasurements = mergeDataSetsOnTimeKey(irMeasurements, uhplcMeasurements, "RelTime")
    mergedMeasurements["Time"] = uhplcMeasurements["Time"]  # Keep original time from UHPLC measurements as timestamp is absolute

    plotData(mergedMeasurements, 
        [{
            "xKey": "RelTime",
            "xKeyFormatter": lambda tInSeconds: tInSeconds/3600,
            "yKeys": ["Meas Bromonitrobenzene (IR)", "Meas Thiophene (IR)", "Meas Product (IR)", "Meas Bromonitrobenzene (UHPLC)", "Meas Thiophene (UHPLC)", "Meas Product (UHPLC)"],
            "title": "Measured Concentrations over Time",
            "xLabel": "Time in hours",
            "yLabel": "Concentration in mol/L",
            "xLimit": [0, 10],
        }]
    )


if __name__ == "__main__":
    main()