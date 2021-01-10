# -*- coding: utf-8 -*-
"""
Created on Sun Feb 23 19:28:32 2020

@author: Sathvik

Version 2.0.7
"""


import math
import turtle


def Temperature():
    # Temperature Calculation
    # Finding Average and Converting C to K
    global hotLiqInitialTemp, hotLiqFinalTemp, coldLiqInitialTemp, coldLiqFinalTemp, TempRawInput
    global hotLiqAvgTemp, coldLiqAvgTemp, hotLiqAvgTempK, coldLiqAvgTempK
    global hotLiq, coldLiq, FlowRateRawInput
    global hotLiqMolecularWt, coldLiqMolecularWt, hotLiqFlowRate, coldLiqFlowRate
    print('Temperature and Liquid Calculation\n')
    hotLiqFlowRate = input('Flow Rate of Hot Liq : ')
    coldLiqFlowRate = input('Flow Rate of Cold Liq : ')
    FlowRateRawInput = input('Flow Rate. Unit : ')
    hotLiqInitialTemp = float(input('Hot Liquid Inital Temp : '))
    hotLiqFinalTemp = float(input('Hot Liquid Final Temp : '))
    coldLiqInitialTemp = float(input('Cold Liquid Initial Temp : '))
    coldLiqFinalTemp = float(input('Cold Liquid Final Temp : '))
    TempRawInput = input('Temp. Unit : ')
    hotLiqAvgTemp = (hotLiqInitialTemp+hotLiqFinalTemp)/2
    coldLiqAvgTemp = (coldLiqInitialTemp+coldLiqFinalTemp)/2
    print('Hot Liq Avg Temp : ', hotLiqAvgTemp)
    print('Cold Liq Avg Temp : ', coldLiqAvgTemp)
    if TempRawInput == 'C':
        hotLiqAvgTempK = hotLiqAvgTemp + 273
        coldLiqAvgTempK = coldLiqAvgTemp + 273
        print('Hot Liquid Average Temp :', hotLiqAvgTempK, 'K')
        print('Cold Liquid Average Temp :', coldLiqAvgTempK, 'K')
    hotLiq = input('Hot Side Liq : ')
    coldLiq = input('Cold Side Liq : ')
    hotLiqMolecularWt = float(input('Molecular Wt of Hot Liq : '))
    coldLiqMolecularWt  = float(input('Molecular Wt of Cold Liq : '))
    
    
def Density():
    # Hot and Cold Density calculation
    global hotLiqDensity, coldLiqDensity
    print('Density Calculation\n')
    if hotLiq == 'water':
        taue = 1-(hotLiqAvgTemp/647.096)
        print('Hot Liq Taue : ', taue)
    elif coldLiq == 'water':
        taue = 1-(coldLiqAvgTemp/647.096)
        print('Cold Liq Taue : ', taue)
    hotC = []
    coldC = []
    if hotLiq == 'water':
        hotC1 = 17.86
        hotC2 = (58.60*(pow(taue, 0.35)))
        hotC3 = -(95.36*(taue**(2/3)))
        hotC4 = (213.89*taue)
        hotC5 = -(141.26*(taue**(4/3)))
        hotCTotal = hotC1+hotC2+hotC3+hotC4+hotC5
        hotLiqDensity = hotCTotal*((18*10**-3)/10**-3)
    else:
        print('Method 1')
        print('Page: 2-96, Table: 2-32')
        for i in range(4):
            temp = input('Hot C[i] : ')
            hotC.append(float(temp))
        hotLiqTotal = (hotC[0]/(hotC[1]**(1+(1-(hotLiqAvgTempK/hotC[2]))**hotC[3])))
        hotLiqDensity = hotLiqTotal*(hotLiqMolecularWt)
    if coldLiq == 'water':
        coldC1 = 17.86
        coldC2 = (58.60*(pow(taue, 0.35)))
        coldC3 = -(95.36*(taue**(2/3)))
        coldC4 = (213.89*taue)
        coldC5 = -(141.26*(taue**(4/3)))
        coldCTotal = coldC1+coldC2+coldC3+coldC4+coldC5
        coldLiqDensity = coldCTotal*((18*10**-3)/10**-3)
    else:
        print('Method 1')
        print('Page: 2-96, Table: 2-32')
        for i in range(4):
            temp = float(input('Cold C[i] : '))
            coldC.append(float(temp))
        coldLiqTotal = (coldC[0]/(coldC[1]**(1+(1-(coldLiqAvgTempK/coldC[2]))**coldC[3])))
        coldLiqDensity = coldLiqTotal*(coldLiqMolecularWt)
    print('Hot Liquid Density : ', hotLiqDensity)
    print('Cold Liquid Density : ', coldLiqDensity)
    

def SpecificHeat():
    global hotLiqSpecificHeat, coldLiqSpecificHeat
    print('SpecificHeat Calculation\n')
    print('Water Page: 2-413, Table: 2-305')
    print('For Method 1: 2-165 Table: 2-153')
    print('For Method 2: 2-208 Table: 2-184')
    hot = input('for Hot Liq: Method 1 Press 1: Formula, or Method 2: Interploation, Method 3: Direct, Enter 1, 2 or 3: ')
    cold = input('for Hot Liq: Method 1 Press Y: Formula, or Method 2: Interploation, Method 3: Direct, Enter 1, 2 or 3: ')
    hotC = []
    coldC = []
    if hot == '1':
        print('Method 1: Formula: ')
        for i in range(5):
            Cp = input('HotC[i]: ')
            hotC.append(float(Cp))
        hotY = hotC[0]+(hotC[1]*hotLiqAvgTempK)+(hotC[2]*(hotLiqAvgTempK**2))+(hotC[3]*(hotLiqAvgTempK**3))+(hotC[4]*(hotLiqAvgTempK**4))
    elif hot == '2':
        print('Method 2: Interpolation: ')
        hotX = hotLiqAvgTempK
        hotX1 = float(input('Hot X1 Temp : '))
        hotX2 = float(input('Hot X2 Temp : '))
        hotY1 = float(input('Hot Y1 SpecificHeat : '))
        hotY2 = float(input('Hot Y2 SpecificHeat : '))
        hotY = (((hotY2-hotY1)/(hotX2-hotX1))*(hotX-hotX1))+hotY1
    else:
        hotY = float(input('Specific Heat of Hot Liq: '))
    if hotLiq == 'water':
        hotLiqSpecificHeat = (hotY/18)*(10**6)
    else:
        hotLiqSpecificHeat = hotY/hotLiqMolecularWt
    print(hotY, hotLiqSpecificHeat)
    if cold == '1':
        print('Method 1: Formula: ')
        for i in range(5):
            Cp = input('ColdC[i]: ')
            coldC.append(float(Cp))
        coldY = coldC[0]+(coldC[1]*coldLiqAvgTempK)+(coldC[2]*(coldLiqAvgTempK**2))+(coldC[3]*(coldLiqAvgTempK**3))+(coldC[4]*(coldLiqAvgTempK**4))
    elif hot == '2':
        print('Method 2')
        coldX = coldLiqAvgTempK
        coldX1 = float(input('Cold X1 Temp : '))
        coldX2 = float(input('Cold X2 Temp : '))
        coldY1 = float(input('Cold Y1 SpecificHeat : '))
        coldY2 = float(input('Cold Y2 SpecificHeat : '))
        coldY = (((coldY2-coldY1)/(coldX2-coldX1))*(coldX-coldX1))+coldY1
    else:
        coldY = float(input('Specific Heat of Cold Liq: '))
    if coldLiq == 'water':
        coldLiqSpecificHeat = (coldY/18)*(10**6)
    else:
        coldLiqSpecificHeat = coldY/coldLiqMolecularWt
    print(coldY, coldLiqSpecificHeat)


def ThermalConductivity():
    global hotLiqThermalConductivity, coldLiqThermalConductivity
    unit = 'W/mK'
    print('Thermal Conductivity Calculation\n')
    if hotLiq == 'water':
        print('Page: 2-413, Table: 2-305')
        hotX = hotLiqAvgTempK
        hotX1 = float(input('Hot X1 Temp : '))
        hotX2 = float(input('Hot X2 Temp : '))
        hotY1 = float(input('Hot Y1 ThermalConductivity : '))
        hotY2 = float(input('Hot Y2 ThermalConductivity : '))
        hotY = (((hotY2-hotY1)/(hotX2-hotX1))*(hotX-hotX1))+hotY1
        hotLiqThermalConductivity = (hotY*(10**-3))#/18
    else:
        print('Page: 2-439, Table: 2-315')
        hotC1 = float(input('Hot C1 : '))
        hotC2 = float(input('Hot C2 : '))
        hotC3 = float(input('Hot C3 : '))
        hotC4 = float(input('Hot C4 : '))
        hotC5 = float(input('Hot C5 : '))
        hotLiqThermalConductivity = hotC1 + (hotC2*hotLiqAvgTempK) + (hotC3*(hotLiqAvgTempK**2)) + (hotC4*(hotLiqAvgTempK**3)) + (hotC5*(hotLiqAvgTempK**4))
    print('Kh: ', hotLiqThermalConductivity, unit)
    if coldLiq == 'water':
        print('Page: 2-413, Table: 2-305')
        coldX = coldLiqAvgTempK
        coldX1 = float(input('Cold X1 Temp : '))
        coldX2 = float(input('Cold X2 Temp : '))
        coldY1 = float(input('Cold Y1 ThermalConductivity : '))
        coldY2 = float(input('Cold Y2 ThermalConductivity : '))
        coldY = (((coldY2-coldY1)/(coldX2-coldX1))*(coldX-coldX1))+coldY1
        coldLiqThermalConductivity = (coldY*(10**-3))
    else:
        print('Page: 2-439, Table: 2-315')
        coldC1 = float(input('Cold C1 : '))
        coldC2 = float(input('Cold C2 : '))
        coldC3 = float(input('Cold C3 : '))
        coldC4 = float(input('Cold C4 : '))
        coldC5 = float(input('Cold C5 : '))
        coldLiqThermalConductivity = coldC1+(coldC2*coldLiqAvgTempK) + (coldC3*(coldLiqAvgTempK**2)) + (coldC4*(coldLiqAvgTempK**3)) + (coldC5*(coldLiqAvgTempK**4))
    print('Kc: ', coldLiqThermalConductivity, unit)

 
def Viscocity():
    global hotLiqViscocity, coldLiqViscocity
    print('Viscocity Calculation: ')
    if hotLiq == 'water':
        print('Page: 2-413, Table: 2-305')
        hotX = hotLiqAvgTempK
        hotX1 = float(input('Hot X1 Temp : '))
        hotX2 = float(input('Hot X2 Temp : '))
        hotY1 = float(input('Hot Y1 Viscocity : '))
        hotY2 = float(input('Hot Y2 Viscocity : '))
        hotY = (((hotY2-hotY1)/(hotX2-hotX1))*(hotX-hotX1))+hotY1
        hotLiqViscocity = hotY*(10**-6)
    else:
        print('Page: 2-427, Table: 2-313')
        hotC1 = float(input('Hot C1 : '))
        hotC2 = float(input('Hot C2 : '))
        hotC3 = float(input('Hot C3 : '))
        hotC4 = float(input('Hot C4 : '))
        hotC5 = float(input('Hot C5 : '))
        hotLiqViscocity = math.exp(hotC1+(hotC2/hotLiqAvgTempK)+(hotC3*(math.log(hotLiqAvgTempK)))+(hotC4*((hotLiqAvgTempK)**hotC5)))
    print(hotLiqViscocity)
    if coldLiq == 'water':
        print('Page: 2-413, Table: 2-305')
        coldX = coldLiqAvgTempK
        coldX1 = float(input('Cold X1 Temp : '))
        coldX2 = float(input('Cold X2 Temp : '))
        coldY1 = float(input('Cold Y1 Viscocity : '))
        coldY2 = float(input('Cold Y2 Viscocity : '))
        coldY = (((coldY2-coldY1)/(coldX2-coldX1))*(coldX-coldX1))+coldY1
        coldLiqViscocity = (coldY*(10**-6))
    else:
        print('Page: 2-427, Table: 2-313')
        coldC1 = float(input('Cold C1 : '))
        coldC2 = float(input('Cold C2 : '))
        coldC3 = float(input('Cold C3 : '))
        coldC4 = float(input('Cold C4 : '))
        coldC5 = float(input('COld C5 : '))
        coldLiqViscocity = math.exp(coldC1+(coldC2/coldLiqAvgTempK)+(coldC3*(math.log(coldLiqAvgTempK)))+(coldC4*((coldLiqAvgTempK)**coldC5)))
    print(coldLiqViscocity)
    

def HeatBalance():
    global hotLiqFlowRate, coldLiqFlowRate, totalFlowRate
    unit = 'kg/s'
    print('Heat Balance Calculation\n')
    print('Hot Side : ', hotLiq)
    print('Cold side : ', coldLiq)
    coldFlowRateRawInput = FlowRateRawInput
    hotFlowRateRawInput = FlowRateRawInput
    if hotLiqFlowRate == '?':
        coldLiqFlowRate = float(coldLiqFlowRate)
        coldLiqFlowRate = coldLiqFlowRate/3600
        coldFlowRateRawInput = 'kg/s'
        hotLiqFlowRate = ((coldLiqFlowRate*coldLiqSpecificHeat*(coldLiqFinalTemp-coldLiqInitialTemp))/(hotLiqSpecificHeat*(hotLiqInitialTemp-hotLiqFinalTemp)))
    if coldLiqFlowRate == '?':
        hotLiqFlowRate = float(hotLiqFlowRate)
        hotLiqFlowRate = hotLiqFlowRate/3600
        hotFlowRateRawInput = 'kg/s'
        coldLiqFlowRate = ((hotLiqFlowRate*hotLiqSpecificHeat*(hotLiqInitialTemp-hotLiqFinalTemp))/(coldLiqSpecificHeat*(coldLiqFinalTemp-coldLiqInitialTemp)))
    hotLiqFlowRate = float(hotLiqFlowRate)
    coldLiqFlowRate = float(coldLiqFlowRate)
    if hotFlowRateRawInput == 'kg/hr':
        hotLiqFlowRate = hotLiqFlowRate/3600
        hotFlowRateRawInput = 'kg/s'
    if coldFlowRateRawInput == 'kg/hr':
        coldLiqFlowRate = coldLiqFlowRate/3600
        coldFlowRateRawInput = 'kg/s'
    totalFlowRate = hotLiqFlowRate*hotLiqSpecificHeat*(hotLiqInitialTemp-hotLiqFinalTemp)
    print('Hot Liq FlowRate: ', hotLiqFlowRate, hotFlowRateRawInput)
    print('Cold Liq FlowRate: ', coldLiqFlowRate, coldFlowRateRawInput)
    print('Q: ', totalFlowRate, unit)
    
    
def LMTD():
    global LMTDCorrected, FT, AT
    deltaT1 = (hotLiqInitialTemp-coldLiqFinalTemp)
    deltaT2 = (hotLiqFinalTemp-coldLiqInitialTemp)
    LMTD = (deltaT2-deltaT1)/(math.log(deltaT2/deltaT1))
    print('deltaT1: ', deltaT1)
    print('deltaT2: ', deltaT2)
    print(LMTD)
    R = (hotLiqInitialTemp-hotLiqFinalTemp)/(coldLiqFinalTemp-coldLiqInitialTemp)
    S = (coldLiqFinalTemp-coldLiqInitialTemp)/(hotLiqInitialTemp-coldLiqInitialTemp)
    print('R: ', R, '\nS: ', S)
    print('Refer Page 11-6 for LMTD Correction Factor: ')
    FT = float(input('FT form the graph: '))
    LMTDCorrected = LMTD*FT
    print('LMTD Corrected : ', LMTDCorrected)
    print('Page : 11-25 Table: 11-3')
    UD = float(input('Design U: '))
    UD *= 5.6783
    dirt = float(input('Total Dirt Factor: '))
    dirt /= 5.6783
    unit = 'm2K/W'
    print('Design U: ', UD, unit)
    print('Dirt Factor: ', dirt, unit)
    AT = totalFlowRate/(UD*LMTDCorrected)
    print('Overall Heat-Tranfer Area: ', AT, 'm2')
    

def TubeSide():
    global tubeOD, tubeID, tubeRawUnit
    global tubedo, tubedi, tubeUnitSI
    global BWG, pitch, PT, L, LSI, LUnitSI, Nt, tubeLiq
    print('Refer Page 11-42 Table 11-12: ')
    tubeLiq = input('Tube Side Liq: ')
    tubeOD = float(input('Tube Outer Dia: '))
    tubeID = float(input('Tube Inner Dai: '))
    tubeRawUnit = input('Tube Unit: ')
    tubedo = tubeOD
    tubedi = tubeID
    if tubeRawUnit == 'in':
        tubedo *= 0.0254
        tubedi *= 0.0254
        tubeUnitSI = 'm'
    BWG = float(input('BWG Gauge: '))
    print('do: ', tubedo, 'm')
    print('di: ', tubedi, 'm')
    print('Refer Page 11-7 Table 11-5: ')
    pitch = float(input('Pitch Input: '))
    PT = pitch*tubedo
    print('PT: ', PT)
    layout = int(input('Enter the Type of Layout: '))
    L = float(input('Length of the Tube, std[6, 8, 12, 16, 20, 24 ft]: '))
    LUnit = input('Length Unit: ')
    LSI = L
    if LUnit == 'ft':
        LSI *= 0.3048
        LUnit = 'm'
    print('Length of the Tubes: ', LSI, LUnit)
    print('Number of Tubes Calculation: ')
    Nt = int(AT/(math.pi*tubedo*LSI))
    print(Nt)
    while Nt % layout != 0:
        Nt += 1
    print(Nt)
    print('Refer Page: 11-34')
    

def TubeSideHTC():
    global At, Gt, tubeNre, tubeNpr, tubehi, tubehio, Np
    Np = int(input('Number of Passes: '))
    # Area of the Tube
    At = ((math.pi/4)*(tubedi**2)*(Nt/Np))
    print('Area of the Tube: ', At, 'm2')
    # Mass Velocity
    if tubeLiq == hotLiq:
        mass = hotLiqFlowRate
        viscocity = hotLiqViscocity
        specificHeat = hotLiqSpecificHeat
        thermalConductivity = hotLiqThermalConductivity
    else:
        mass = coldLiqFlowRate
        viscocity = coldLiqViscocity
        specificHeat = coldLiqSpecificHeat
        thermalConductivity = coldLiqThermalConductivity
    Gt = mass/At
    print('Gt :', Gt, 'kg/m2s')
    # Nre, Npr Calculation
    tubeNre = (tubedi*Gt)*viscocity
    tubeNpr = (specificHeat*viscocity)/thermalConductivity
    print('Nre :', tubeNre, '\nNpr :', tubeNpr)
    if tubeNre > 10000:
        # Ditus-Bolter Equation
        tubehi = 0.023*(tubeNre**0.8)*(tubeNpr**(1/3))*(thermalConductivity/tubedi)
        print('hi :', tubehi)
        tubehio = tubehi*(tubedi/tubedo)
        print('hio :', tubehio)
    else:
        print('Nre < 10000')
    

def ShellDia():
    global C, D, Ds
    # Find C from
    # Nt = a + b*C + c*C**2 + d*C**3 + e*C**4
    C = float(input('C: '))
    if C <= -24:
        C = -24
    if C >= 24:
        C = 24
    d = tubedo
    D = ((C+36)*d)/0.75
    print('C: ', C, '\nD: ', D, 'm')
    clearence = float(input('Refer Page 11-40 for the Clearence Value: '))
    Ds = D+2*clearence
    print('Ds: ', Ds, 'm')
        
    
def ShellSide():
    global Dotl, Pp, Pn, Nc, Ncw, ls, Sm, Fbp, Stb, Ssb, Sw, Nb, shellLiq, shellNre, b
    Pdash = PT
    Pp = float(input('Refer Page 11-7 for tubedo: '))
    Pn = float(input('Refer Page 11-7 for tubeid: '))
    # Number of Tube Crosses
    # lc/dc = lc
    Dotl = D
    lc = 0.25
    Nc = (Ds*(1-(2*lc)))/Pp
    # Fraction of Total Tubes in Cross-Flow
    lc = float(input('lc: '))
    Fc = ((1/math.pi)*(math.pi*(2*(Ds-(2*lc))/Dotl)*math.sin(math.cosh*((Ds-(2*lc))/Dotl))-(2*(math.cosh((Ds-(2*lc))/Dotl)))))
    print('Fc: ', Fc)
    # Number of Effective
    # Page 11-8
    Ncw = (0.8*lc)/Pp
    # Cross Flow Area
    ls = 0.2*Ds
    Sm = ls*((Ds-Dotl)+((Dotl-tubedo)/Pdash)*(Pdash-tubedo))
    print('Sm :', Sm)
    # Fbp
    Fbp = ((Ds-Dotl)*ls)/Sm
    print(Fbp)
    # Tube to Baffle leakage area for one baffle Stb
    b = 6.223*(10**-4)
    Stb = b*tubedo*Nt(1+Fc)
    print('Tube Baffle Leakage Area: ', Stb, 'm2')
    # Shell to Baffle Leakage Area Ssb
    dauesb = float(input('dauesb: '))
    Ssb = ((Ds*dauesb)/2)*(math.pi-math.cosh(1-((2*ls)/Ds)))
    print('Shell to Baffle Leakage Area Ssb: ', Ssb, 'm2')
    # Area for flow through window Sw
    Swg = ((Ds**2)/4)*((math.cosh(1-(2*(lc/Ds))))-((1-(2*(lc/Ds)))*(math.sqrt(1-((1-2*(lc/Ds))**2)))))
    Swt = (Nt/8)*((1-Fc)*math.pi*(tubedo**2))
    Sw = Swg-Swt
    print('Area for flow through window Sw: ', Sw, 'm2')
    shellLiq = input('Shell Liq: ')
    if shellLiq == hotLiq:
        mass = hotLiqFlowRate
        viscocity = hotLiqViscocity
    else:
        mass = coldLiqFlowRate
        viscocity = coldLiqViscocity
    shellNre = (tubedo*mass)/(viscocity*Sm)
    if shellNre <= 100:
        titab = 2*math.cosh(1-((2*lc)/Ds))
        Dw = ((4*Sw)/((math.pi/2)*Nt*(1-Fc)*tubedo+(Ds*titab)))
        print('Titab: ', titab, '\nDw: ', Dw)
    # Number of Baffles Nb
    Nb = round((L/lc)-1)
    print('Number of Baffles: ', Nb)
    

def ShellSideHTC():
    global Jk, hk, Jc, Jl, Nss, Jb, Jb, shellho
    Jk = float(input('Refer Page 11-9 JkFactor: '))
    if shellLiq == hotLiq:
        mass = hotLiqFlowRate
        viscocity = hotLiqViscocity
        specificHeat = hotLiqSpecificHeat
        thermalConductivity = hotLiqThermalConductivity
    else:
        mass = coldLiqFlowRate
        viscocity = coldLiqViscocity
        specificHeat = coldLiqSpecificHeat
        thermalConductivity = coldLiqThermalConductivity
    shellNre = (tubedo*mass)/(viscocity*Sm)
    print('ShellSide Nre: ', shellNre)
    hk = Jk*specificHeat*(mass/Sm)*((thermalConductivity/(specificHeat*viscocity))**(2/3))
    print('hk :', hk, 'W/m2K')
    Jc = float(input('Refer Table 11-10 for Jc: '))
    print('Ssb+Stb/Sm: ', (Ssb+Stb)/Sm)
    print('Ssb/Ssb+Stb: ', (Ssb/(Ssb+Stb)))
    Jl = float(input('Refer Table 11-11 for Jl: '))
    Nss = float(input('Nss Value: '))
    print('Fbp', Fbp)
    print('Nss/Nc: ', Nss/Nc)
    Jb = float(input('Refer Table 11-12 for Jb: '))
    if shellNre > 100:
        Jr = 1
    else:
        Jr = float(input('Refer Table 11-14 for Jr: '))
    shellho = hk*Jc*Jl*Jb*Jr
    print('Shell Side ho: ', shellho, 'W/m2K')
    

def VelocityCorrectionFactor():
    # Incomplete
    print('Velocity Correction Factor is Incomplete')
    
    
def ShellSidePDrop():
    global fk, deltaPwk, Rl, deltaPs
    fk = float(input('Refer Chart 11-10 for fk: '))
    if shellLiq == hotLiq:
        mass = hotLiqFlowRate
        density = hotLiqDensity
    else:
        mass = coldLiqFlowRate
        density = coldLiqDensity
    deltaPbk = ((b*fk*(mass**2)*Nc)/(hotLiqDensity*(Sm**2)))*((hotLiqViscocity/coldLiqViscocity)**0.14)
    deltaPwk = (b*(mass**2)*(2+(0.6*Ncw)))/(Sm*Sw*density)
    print('Ssb+Stb/Sm: ', (Ssb+Stb)/Sm)
    print('Ssb/Ssb+Stb: ', Ssb/(Ssb+Stb))
    Rl = float(input('Refer Chart 11-16 for Rl: '))
    print('Nss/Nc: ', Nss/Nc)
    print('Fbp: ', Fbp)
    Rb = float(input('Refer Chart 11-17 for Rb: '))
    deltaPs = ((((Nb-1)*deltaPbk*Rb)+(Nb*deltaPwk))*Rl)+(2*deltaPbk*Rb*(1+(Ncw/Nc)))
    deltaPs /= 144
    if deltaPs < 10:
        print('Design is Okay')
    else:
        print('ERROR')
    print('deltaPbk: ', deltaPbk, 'lbf/ft2')
    

def TubeSidePDrop():
    global f, g, deltaPt
    f = float(input('Refer Chart 6-9 for f: '))
    g = 9.8
    if tubeLiq == hotLiq:
        density = hotLiqDensity
    else:
        density = coldLiqDensity
    deltaPt = (((4*f*L*(Gt**2))/(2*g*density*tubedi))+((4*(Gt**2))/2*g*density))*Np
    deltaPt = (deltaPt/4.88)*(1/144)
    print('deltaPt: ', deltaPt)
    if deltaPt < 10:
        print('Design is Okay')
    else:
        print('ERROR')
        

#Temperature()
#Density()
#SpecificHeat()
#ThermalConductivity()
#Viscocity()
#HeatBalance()
#LMTD()
#TubeSide()
#TubeSideHTC()
#ShellDia()
#ShellSide()
#ShellSideHTC()
#VelocityCorrectionFactor()
#ShellSidePDrop()
        

turtle.forward(3.5*100)
turtle.right(90)
turtle.forward(1*100)
turtle.right(90)
turtle.forward(3.5*100)
turtle.right(90)
turtle.forward(1*100)
turtle.backward(0.1*100)
turtle.forward(2*100)