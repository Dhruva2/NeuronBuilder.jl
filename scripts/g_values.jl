"""
Sets of Prinz conductances from 2003 and 2004 papers 
"""

# parameters from Prinz 2003 https://journals.physiology.org/doi/pdf/10.1152/jn.00641.2003, figure 3.
burster_channels = [Prinz.Na(400.0 * Prinz_conv), Prinz.CaS(8.0 * Prinz_conv), Prinz.CaT(0.0 * Prinz_conv), Prinz.H(0.04 * Prinz_conv),
    Prinz.Ka(50.0 * Prinz_conv), Prinz.KCa(20.0 * Prinz_conv), Prinz.Kdr(50.0 * Prinz_conv), Prinz.leak(0.0)]

# parameters from Xolotl, which in turn come from Prinz 2004 "Similar network activity from disparate circuit parameters" 
AB2_ch = [Prinz.Na(100 * Prinz_conv), Prinz.CaS(6 * Prinz_conv), Prinz.CaT(2.5 * Prinz_conv), Prinz.H(0.01 * Prinz_conv),
    Prinz.Ka(50 * Prinz_conv), Prinz.KCa(5 * Prinz_conv), Prinz.Kdr(100 * Prinz_conv), Prinz.leak(0.0)]
LP4_ch = [Prinz.Na(100 * Prinz_conv), Prinz.CaS(4 * Prinz_conv), Prinz.CaT(0.0), Prinz.H(0.05 * Prinz_conv), Prinz.Ka(20 * Prinz_conv),
    Prinz.KCa(0.0), Prinz.Kdr(2.5 * Prinz_conv), Prinz.leak(0.03 * Prinz_conv)]
PY1_ch = [Prinz.Na(100 * Prinz_conv), Prinz.CaS(2 * Prinz_conv), Prinz.CaT(2.4 * Prinz_conv), Prinz.H(0.05 * Prinz_conv), Prinz.Ka(50 * Prinz_conv),
    Prinz.KCa(0.0), Prinz.Kdr(125 * Prinz_conv), Prinz.leak(0.01 * Prinz_conv)]

# parameters from wiktor phillip's Conductor demos
AB_channels = [Prinz.Na(100.0 * Prinz_conv), Prinz.CaS(3.0 * Prinz_conv), Prinz.CaT(1.3 * Prinz_conv), Prinz.H(0.5 * Prinz_conv),
    Prinz.Ka(5.0 * Prinz_conv), Prinz.KCa(10.0 * Prinz_conv), Prinz.Kdr(20.0 * Prinz_conv), Prinz.leak(0.01 * Prinz_conv)]


"""
Sets of Liu conductances from various sources. 
"""

# parameters from Drion et. al. 2015 "Dynamic Input Conductances Shape Neuronal Spiking" https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4586923/
Burster1 = [Liu.Na(700.0 * Liu_conv), Liu.CaS(4.0 * Liu_conv), Liu.CaT(2.0 * Liu_conv), Liu.Ka(50.0 * Liu_conv), Liu.KCa(40.0 * Liu_conv),
    Liu.Kdr(70.0 * Liu_conv), Liu.H(0.03 * Liu_conv), Liu.leak(0.01 * Liu_conv)]

Burster2 = [Liu.Na(700.0 * Liu_conv), Liu.CaS(2.25 * Liu_conv), Liu.CaT(6.25 * Liu_conv), Liu.Ka(85.0 * Liu_conv), Liu.KCa(50.0 * Liu_conv),
    Liu.Kdr(90.0 * Liu_conv), Liu.H(0.03 * Liu_conv), Liu.leak(0.01 * Liu_conv)]

# parameters from O'Leary (2014) 10.1016/j.neuron.2014.04.002  
# Burster. These don't have to multiplied by Liu_conv
AB_channels = [Liu.Na(1831.0), Liu.CaS(27.0), Liu.CaT(23.0), Liu.H(10.1), Liu.Ka(246.0), Liu.KCa(980.0), Liu.Kdr(610.0), Liu.leak(0.99)]


"""
Hodgkin-Huxley model conductances
tonic spiking when an applied current Iapp>2 is used
"""
HH = [Liu.Na(1000), Liu.leak(0.1), Liu.Kdr(200)]

