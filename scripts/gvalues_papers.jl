
````
Sets of Prinz conductances from various sources. 
````

# parameters from Prinz 2003 https://journals.physiology.org/doi/pdf/10.1152/jn.00641.2003, figure 3.
burster_channels = [Prinz.NaV(400.0 * P2NB), Prinz.CaS(8.0 * P2NB), Prinz.CaT(0.0 * P2NB), Prinz.H(0.04 * P2NB),
    Prinz.Ka(50.0 * P2NB), Prinz.KCa(20.0 * P2NB), Prinz.Kdr(50.0 * P2NB), Prinz.Leak(0.0)]

# parameters from Xolotl, which in turn come from Prinz 2004 " Similar network activity from disparate circuit parameters" 
AB2_ch = [Prinz.NaV(100 * P2NB), Prinz.CaS(6 * P2NB), Prinz.CaT(2.5 * P2NB), Prinz.H(0.01 * P2NB),
    Prinz.Ka(50 * P2NB), Prinz.KCa(5 * P2NB), Prinz.Kdr(100 * P2NB), Prinz.Leak(0.0)]
LP4_ch = [Prinz.NaV(100*P2NB), Prinz.CaS(4*P2NB), Prinz.CaT(0.0), Prinz.H(0.05*P2NB), Prinz.Ka(20*P2NB),
     Prinz.KCa(0.0), Prinz.Kdr(2.5*P2NB), Prinz.Leak(0.03*P2NB)]
PY1_ch = [Prinz.NaV(100*P2NB), Prinz.CaS(2*P2NB), Prinz.CaT(2.4*P2NB), Prinz.H(0.05*P2NB), Prinz.Ka(50*P2NB),
     Prinz.KCa(0.0), Prinz.Kdr(125*P2NB), Prinz.Leak(0.01*P2NB)]

# wiktor phillip's conductor
AB_channels = [Prinz.NaV(100.0 * P2NB), Prinz.CaS(3.0 * P2NB), Prinz.CaT(1.3 * P2NB), Prinz.H(0.5 * P2NB), Prinz.Ka(5.0 * P2NB), Prinz.KCa(10.0 * P2NB), Prinz.Kdr(20.0 * P2NB), Prinz.Leak(0.01 * P2NB)]


````
Sets of Liu conductances from various sources. 
````

#From MyModelMenagerie. These don't have to be multiplied by L2NB
stg_liu = [Liu.Na(500.0), Liu.CaS(60.0), Liu.CaT(25), Liu.Ka(250.0), Liu.KCa(40.0), Liu.Kdr(500.0), Liu.H(0.1), Liu.Leak(0.0)]

#dynamic input conductances. Burster
DIC = [Liu.Na(700.0 * L2NB), Liu.CaS(4.0 * L2NB), Liu.CaT(2.0 * L2NB), Liu.Ka(50.0 * L2NB), Liu.KCa(40.0 * L2NB), Liu.Kdr(70.0 * L2NB), Liu.H(0.03 * L2NB), Liu.Leak(0.01 * L2NB)]

#Burster 2
DIC_2 = [Liu.Na(700.0 * L2NB), Liu.CaS(2.25 * L2NB), Liu.CaT(6.25 * L2NB), Liu.Ka(85.0 * L2NB), Liu.KCa(50.0 * L2NB), Liu.Kdr(90.0 * L2NB), Liu.H(0.03 * L2NB), Liu.Leak(0.01 * L2NB)]

# parameters from O'Leary 2014. Burster
#these don't have to be multiplied by L2NB
AB1_channels = [Liu.Na(1831.0), Liu.CaS(27.0), Liu.CaT(23.0), Liu.H(10.1), Liu.Ka(246.0), Liu.KCa(980.0), Liu.Kdr(610.0), Liu.Leak(0.99)]

#high bursts
xolotl_AB_liu = [Liu.Na(100*L2NB), Liu.CaS(3*L2NB), Liu.CaT(1.3*L2NB), Liu.Ka(5*L2NB), Liu.KCa(10*L2NB), Liu.Kdr(20*L2NB), Liu.H(.5*L2NB), Liu.Leak(0.01*L2NB)]

#tonic spiking, these don't have to be multiplied by L2NB
HH = [Liu.Na(1000), Liu.Leak(0.1), Liu.Kdr(200)]

