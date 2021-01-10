# -*- coding: utf-8 -*-
"""
Created on Sun Feb 16 09:44:54 2020

@author: Sathvik

Version 1.0.5
"""

import math

# definig the liquids in shell and tube
benzWeight = 20000
benzUnit = 'kg/hr'
waterWeight = '?'
waterUnit = 'kg/hr'

# temperature of the cold and hot fluid
waterInitialTemp = 80
waterFinalTemp = 60
benzInitialTemp = 25
benzFinalTemp = 65
TempUnit = 'C'

# Molecular Weight
waterMolecularWt = 18
benzMolecularWt = 78


def step1():
    # finding the water and benzene mean Temp using (T1+T2)/2
    global waterTempAvg, benzTempAvg, waterTempAvgK, benzTempAvgK
    waterTempAvg = (waterInitialTemp + waterFinalTemp)/2
    benzTempAvg = (benzInitialTemp + benzFinalTemp)/2
    print('T(Wavg) :', waterTempAvg, 'C')
    print('T(Bavg) :', benzTempAvg, 'C')
    # Converting C to K 
    waterTempAvgK = waterTempAvg + 273
    benzTempAvgK = benzTempAvg + 273
    print('T(Wavg) :', waterTempAvgK, 'K')
    print('T(Bavg) :', benzTempAvgK, 'K')
    

def step2():
    # Density of Water and Benzene 
    global waterTempAvgDensitySI, benzTempAvgDenistySI
    taue = 1 - (waterTempAvgK/647.096)
    print(taue)
    waterTempAvgDensity1 = 17.86
    waterTempAvgDensity2 = (58.60*(pow(taue, 0.35)))
    waterTempAvgDensity3 = -(95.36*(taue**(2/3)))
    waterTempAvgDensity4 = (213.89*taue)
    waterTempAvgDensity5 = -(141.26*(taue**(4/3)))
    print(waterTempAvgDensity1, waterTempAvgDensity2, waterTempAvgDensity3, waterTempAvgDensity4, waterTempAvgDensity5)
    waterTempAvgDensity = waterTempAvgDensity1 + waterTempAvgDensity2 + waterTempAvgDensity3 + waterTempAvgDensity4 + waterTempAvgDensity5
    print(waterTempAvgDensity)
    waterTempAvgDensitySI = waterTempAvgDensity * ((18*10**-3)/10**-3)
    print('Density of Water:', waterTempAvgDensitySI, 'kg/m3')
    c1 = 1.0259
    c2 = 0.266
    c3 = 562.05
    c4 = 0.2839
    benzTempAvgDensity = (c1/(c2**(1+(1-(benzTempAvgK/c3))**c4)))
    print(benzTempAvgDensity)
    benzTempAvgDensitySI = 851.56 # benzTempAvgDensity*(78*10**-3/10**-3)
    print('Density of Benzene :', benzTempAvgDensitySI, 'kg/m3')
    
    
def step3():
    # Heat Capacity of Water and Benzene
    global heatCapacityWaterJperkK, heatCapacityBenzJperkK
    # c1 [0] = water, c2 [1] = benzene
    C1 = [276370, 129440]
    C2 = [-2090.1, -169.5]
    C3 = [8.125, 0.647]
    C4 = [-0.0141, 0]
    C5 = [9.370*10**-6, 0]
    for i in range(0, 1):
        heatCapacityWater = C1[i] + (C2[i]*waterTempAvgK) + (C3[i]*(waterTempAvgK**2)) + (C4[i]*(waterTempAvgK**3)) + (C5[i]*(waterTempAvgK**4))
        heatCapacityWaterJperkK = heatCapacityWater/waterMolecularWt
        print('Specific Heat Capacity of Water :', heatCapacityWaterJperkK, 'J/kK')
    for i in range(1, 2):
        heatCapacityBenz = C1[i] + (C2[i]*benzTempAvgK) + (C3[i]*(benzTempAvgK**2)) + (C4[i]*(benzTempAvgK**3)) + (C5[i]*(benzTempAvgK**4))
        heatCapacityBenzJperkK = heatCapacityBenz/benzMolecularWt
        print('Specific Heat Capacity of Benzene :', heatCapacityBenzJperkK, 'J/kK')
        

def step4():
    # Thermal Conductivity and Viscocity of Water and Benzene
    global thermalConductivityWater, thermalCondutivityBenz
    # C1 [0] = water, C2 [1] = benzene
    C1 = [-0.432, 0.23444]
    C2 = [0.0057255, -0.00030572]
    C3 = [-0.000008078, 0]
    C4 = [1.861*10**-9, 0]
    C5 = [0, 0]
    for i in range(0, 1):
        thermalConductivityWater = C1[i] + (C2[i]*waterTempAvgK) + (C3[i]*(waterTempAvgK**2)) + (C4[i]*(waterTempAvgK**3)) + (C5[i]*(waterTempAvgK**4))
        print('Thermal Conductivity of Water :', thermalConductivityWater, 'W/mK')
    for i in range(1, 2):
        thermalConductivityBenz = C1[i] + (C2[i]*benzTempAvgK) + (C3[i]*(benzTempAvgK**2)) + (C4[i]*(benzTempAvgK**3)) + (C5[i]*(benzTempAvgK**4))
        print('Thrmal Conductivity of Benzene', thermalConductivityBenz, 'W/mK')
    # find the viscocity
        

def step5():
    # Heat Balance and FlowRate
    global benzWeightSI, waterWeightSI, waterUnitSI, totalMassFlowRate
    # water - hot liquid, benzene - cold liquid
    if benzUnit == 'kg/hr':
        benzWeightSI = benzWeight/3600
    print('MassFlow Rate of Benzene :', benzWeightSI, 'kg/s')
    waterUnitSI = 'kg/s'
    waterWeightSI = ((benzWeightSI*heatCapacityBenzJperkK*(benzFinalTemp-benzInitialTemp))/(heatCapacityWaterJperkK*(waterInitialTemp-waterFinalTemp)))
    print('MassFlow Rate of Water :', waterWeightSI, waterUnitSI)
    totalMassFlowRateBenz = benzWeightSI*heatCapacityBenzJperkK*(benzFinalTemp-benzInitialTemp)
    # totalMassFlowRateWater = waterWeightSI*heatCapacityWaterJperkK*(waterInitialTemp-waterFinalTemp)
    totalMassFlowRate = totalMassFlowRateBenz # totalMassFlowRateBenz == totalMassFlowRateWater == 401613.1851
    print('Total MassFlow Rate Q : ', totalMassFlowRate, 'J/s')
    

def step6():
    # LMTD Calculation
    global deltaLMTDCorrectedFt
    deltaT1 = waterInitialTemp-benzFinalTemp
    deltaT2 = waterFinalTemp-benzInitialTemp
    LMTD = (deltaT2-deltaT1)/math.log(deltaT2/deltaT1)
    print(LMTD)
    R = (waterInitialTemp-waterFinalTemp)/(benzFinalTemp-benzInitialTemp)
    S = (benzFinalTemp-benzInitialTemp)/(waterInitialTemp-benzInitialTemp)
    LMTDCorrectionFt = 0.95 # calcultated from the graph
    deltaLMTDCorrectedFt = LMTDCorrectionFt*LMTD
    print(R, S, deltaLMTDCorrectedFt)


def step7():
    # Ud and Atra of the Tube
    # Water - Tube, Benzene - Shell
    global Ud, UdUnit, UdAssumed
    global Rd, RdUnit, At
    Ud = 100 # between 50-150
    UdUnit = 'Btu/Fft2h'
    Rd = 0.003
    RdUnit = 'Fft2h/Btu'
    Rd /= 5.6783
    Ud *= 5.6783
    UdAssumed = Ud
    print('Ud : ', Ud, UdUnit)
    # Overall HeatTransfter Area
    At = totalMassFlowRate/(UdAssumed*deltaLMTDCorrectedFt)
    print('Area of the HeatTransfer : ', At, ' m2')
    

def step8():
    # Tube Dia, Length, Number
    global tubedo, tubedi, tubeUnit, BWGGauge
    global Pt, lengthTube, lengthTubeUnit, numberTubes
    global numberShell, numberPasses
    tubeOD = 1
    tubeID = 0.870
    tubeRawUnit = 'in'
    BWGGauge = 16
    if tubeRawUnit == 'in':
        tubedo = tubeOD*0.0254
        tubedi = tubeID*0.0254
        tubeUnit = 'm'
        print('Outer Dia of Tube : ', tubedo, tubeUnit)
        print('Inner Dia of Tube : ', tubedi, tubeUnit)
    else: 
        tubedo = tubeOD
        tubedi = tubeID
        tubeUnit = tubeRawUnit
    print(tubedo)
    # Assume Triangluar Pitch
    Pt = 1.25*tubedo
    PtUnit = 'm'
    print('Pitch : ', Pt, PtUnit)
    # Standard Lenght = 6, 8, 12, 16, 20, 24 ft
    lengthTubeRaw = 20
    lengthTubeRawUnit = 'ft'
    if lengthTubeRawUnit == 'ft':
        lengthTube = lengthTubeRaw*0.3048
        lengthTubeUnit = 'm'
    else:
        lengthTube = lengthTubeRaw
        lengthTubeUnit = lengthTubeRawUnit
    print('Length of the Tube : ', lengthTube, lengthTubeUnit)
    numberTubes = At/(math.pi*tubedo*lengthTube)
    numberTubes = int(round(numberTubes))
    for i in range(numberTubes, numberTubes+9):
        if numberTubes % 3 != 0:
            numberTubes += 1
    print('Number of Tubes : ', numberTubes)
    # [2 - 4 SHTE], 2 - Shell, 4 - Passes
    numberShell = 2
    numberPasses = 4
    print('Number of Shell : ', numberShell)
    print('Number of Passes : ', numberPasses)
    

def step9():
    # Area of the Tube
    At = (math.pi/4)*(tubedi**2)*(numberTubes/numberPasses)
    print('Area of the Tube : ', At, 'm')
    # Mass Velocity
    Gt = waterWeightSI/At
    print('Mass FlowRate : ', Gt, 'kg/m2s')
    # Raynolds Number
    # Nre = (tubeid*Gt)/ # Find The Viscocity
        

step1()
step2()
step3()
step4()
step5()
step6()
step7()
step8()
step9()