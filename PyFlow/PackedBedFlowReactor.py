import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

from PyFlow.FlowReactor import FlowReactor

class PackedBedFlowReactor(FlowReactor):
    def __init__(self, length=1.0, diameter=1.0e-3, packedBedInitQuality=1.0, **kwargs):
        self.packedBedInitQuality = packedBedInitQuality
        super().__init__(length=length, diameter=diameter, **kwargs)

        self.lastSimResultForPackedBed = None
        

    def setSpaceSamples(self, setSpaceSamples):
        self.xSamples = setSpaceSamples
        self.packedBedQuality = self.packedBedInitQuality * np.ones(self.xSamples)

    def __str__(self):
        return f"Packed-Bed Flow Reactor (volume={1000*self.getReactorVolumeInLiters()} mL, length={self.length} m, diameter={1000*self.diameter} mm, D={self.dispersionCoefficient})"
    

    def simulateStep(self, startTime, endTime, timeVec, Cin, flowRate, temperature=20, isVolumetricFlowRate=True, method="RK45", rtol=0.000001, atol=0.000001, Cspatial0=None):
        numberOfSpecies = Cin.shape[0]

        if Cspatial0 is None:
            if self.lastSimResult is not None:
                _, prevCspatial = self.lastSimResult
                Cspatial0 = prevCspatial[:, 1:, -1]
            else:
                Cspatial0 = np.zeros((numberOfSpecies, self.xSamples))

        
        if Cspatial0.shape[0] == numberOfSpecies:
            # Need to add packed bed quality to Cspatial0
            Qspatial0 = self.packedBedInitQuality * np.ones(self.xSamples)
            if self.lastSimResultForPackedBed is not None:
                _, Qspatial0 = self.lastSimResultForPackedBed
                Qspatial0 = Qspatial0[:, -1]

            Cspatial0 = np.vstack([Cspatial0, Qspatial0])

            # Add one more Cin which is the placeholder for packed bed quality
            Cin = np.vstack([Cin, np.zeros(len(timeVec))])

        time, Cout, Cspatial = super().simulateStep(startTime, endTime, timeVec, Cin, flowRate, temperature, isVolumetricFlowRate, method, rtol, atol, Cspatial0)
        # store and remove packed bed quality from Cout and Cspatial
        packedBedSimResult = Cspatial[-1, :, :].copy()  # store packed bed quality evolution
        self.lastSimResultForPackedBed = (time, packedBedSimResult[1:, :])  # skip inlet point for better visualization

        Cout = Cout[:-1, :]
        Cspatial = Cspatial[:-1, :, :]
        self.lastSimResult = (time, Cspatial)

        return time, Cout, Cspatial
    

    def simulate(self, timeVec, Cin, **kwargs):
        nSpecies = Cin.shape[0]
        Cspatial0 = kwargs.get("Cspatial0", None)

        # Add one more Cin which is the placeholder for packed bed quality
        Cin = np.vstack([Cin, np.zeros(len(timeVec))])

        # Prepare Cspatial0 with packed bed quality
        Qspatial0 = self.packedBedInitQuality * np.ones(self.xSamples)
        if Cspatial0 is None: Cspatial0 = np.zeros((nSpecies, self.xSamples))

        Cspatial0 = np.vstack([Cspatial0, Qspatial0])
        kwargs["Cspatial0"] = Cspatial0

        startTime = timeVec[0]
        endTime = timeVec[-1]
        time, Cout, Cspatial = super().simulateStep(startTime, endTime, timeVec, Cin, **kwargs)

        # store and remove packed bed quality from Cout and Cspatial
        packedBedSimResult = Cspatial[-1, :, :].copy()  # store packed bed quality evolution
        self.lastSimResultForPackedBed = (time, packedBedSimResult[1:, :])  # skip inlet point for better visualization

        Cout = Cout[:-1, :]
        Cspatial = Cspatial[:-1, :, :]
        self.lastSimResult = (time, Cspatial)

        return time, Cout, Cspatial
    
    def get_dCdt(self, t, Cvec, nSpecies, dx, timeVec, Cin, flowRate, temperature):
        # Cvec includes all species and packed bed quality

        # Csvec = [A1, A2, A3, ... AN, B1, B2, B3, ... BN, ..., Q1, Q2, Q3, ... QN] 
        # A = C1, B = C2, ...
        # A1 at pos 1

        dCdx = np.zeros_like(Cvec)
        dC2d2x = np.zeros_like(Cvec)
        rCi = np.zeros_like(Cvec)

        # It might be, that the ivp_solve also uses time samples in between the time samples of Cin
        input_fcn = interp1d(timeVec, Cin, kind='linear', fill_value="extrapolate")
        temperature_fcn = (lambda _: temperature) if np.isscalar(temperature) else interp1d(timeVec, temperature, kind='linear', fill_value="extrapolate")
        flowRate_fcn = (lambda _: flowRate) if np.isscalar(flowRate) else interp1d(timeVec, flowRate, kind='linear', fill_value="extrapolate")

        for i in range(nSpecies): 
            sIdx = i * self.xSamples
            eIdx = (i + 1) * self.xSamples

            if i < nSpecies - 1:
                # Packed bed quality does not undergo convection/dispersion
                # Convection
                Ci_spatial = Cvec[sIdx:eIdx]
                dCdx[sIdx+1:eIdx] = (Ci_spatial[1:] - Ci_spatial[:-1]) / dx
                dCdx[sIdx] = (Ci_spatial[0] - input_fcn(t)[i]) / dx

                # Dispersion
                dC2d2x[sIdx] = (input_fcn(t)[i] - 2*Ci_spatial[0] + Ci_spatial[1]) / dx**2

                d2Cidx2_center = (Ci_spatial[:-2] - 2*Ci_spatial[1:-1] + Ci_spatial[2:]) / dx**2
                d2Cidx2_end = (Ci_spatial[-2] - Ci_spatial[-1]) / dx**2
                dC2d2x[sIdx+1:eIdx] = np.append(d2Cidx2_center, d2Cidx2_end)

            # Reaction
            if self.reactionNetworkCallback is not None:
                C_Q_spatial = Cvec.reshape(nSpecies, self.xSamples)
                rCi[sIdx:eIdx] = self.reactionNetworkCallback(C_Q_spatial[:-1, :], i if i < nSpecies - 1 else -1, temperature_fcn(t), C_Q_spatial[-1, :])

        return self.dispersionCoefficient * dC2d2x -1 * flowRate_fcn(t) * dCdx + rCi
    
    def plotSpaceTime(self, title="Reactor Time Evolution (3D)", figure=None, rows=1, rowIdx=1, showPlot=False):
        if self.lastSimResult is None: return

        nSpecies = self.getNumberOfSpecies()
        packedBedResultsAvailable = self.lastSimResultForPackedBed is not None
        numberOfCols = nSpecies + (1 if packedBedResultsAvailable else 0)

        if figure is None: figure = plt.figure()
        super().plotSpaceTime(title=title, figure=figure, rows=rows, rowIdx=rowIdx, showPlot=False, _cols=numberOfCols)

        # Plot packed bed quality
        if packedBedResultsAvailable:
            time, packedBedSimResult = self.lastSimResultForPackedBed
            X, Y = np.meshgrid(time, self.getReactorSpaceSamples()[1:])  # skip inlet point for better visualization

            ax = figure.add_subplot(rows, numberOfCols, (((nSpecies + 1) * (rowIdx-1)) + nSpecies+1), projection='3d')
            surf = ax.plot_surface(X, Y, packedBedSimResult, cmap='viridis')
            ax.set_title("Packed Bed Quality")
            ax.set_xlabel("Time t in s")
            ax.set_ylabel("Length x in m")
            ax.set_zlabel("Packed Bed Quality")
            ax.view_init(elev=45, azim=45)
            ax.set_zlim(0, 1.1)


        figure.colorbar(surf, ax=ax, shrink=0.5, aspect=5)
        if showPlot: plt.show()

    def plotOutputConcentration(self, Cin=None, labelCout="Cout", measuredCout=None, measuredCoutTimeVec=None, CinTimeVec=None, additionalTracesToPlot=None, title="Reactor Time Evolution (2D)", figure=None, rows=1, rowIdx=1, showPlot=True):  
        if self.lastSimResult is None: return

        nSpecies = self.getNumberOfSpecies()
        packedBedResultsAvailable = self.lastSimResultForPackedBed is not None
        numberOfCols = nSpecies + (1 if packedBedResultsAvailable else 0)
        
        if figure is None: figure = plt.figure()

        super().plotOutputConcentration(
            Cin=Cin, 
            labelCout=labelCout, 
            measuredCout=measuredCout, 
            measuredCoutTimeVec=measuredCoutTimeVec, 
            CinTimeVec=CinTimeVec, 
            additionalTracesToPlot=additionalTracesToPlot, 
            title=title, 
            figure=figure, 
            rows=rows, 
            rowIdx=rowIdx, 
            showPlot=False, 
            _cols=numberOfCols
        )

        # Plot packed bed quality
        if packedBedResultsAvailable:
            time, packedBedSimResult = self.lastSimResultForPackedBed
            
            # plot of mean and std dev
            ax = figure.add_subplot(rows, numberOfCols, (((nSpecies + 1) * (rowIdx-1)) + nSpecies + 1))
            ax.plot(time, packedBedSimResult.mean(axis=0), label="Avg Packed Bed Quality")
            ax.fill_between(time, packedBedSimResult.min(axis=0), packedBedSimResult.max(axis=0), alpha=0.2)
            ax.set_title("Packed Bed Quality Evolution")
            ax.set_xlabel("Time t in s")
            ax.set_ylabel("Packed Bed Quality")
            ax.set_ylim(0, 1.1 * self.packedBedInitQuality)
            ax.legend()

        if showPlot: plt.show()

    