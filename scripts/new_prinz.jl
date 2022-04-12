const Area = 0.0628 # Prinz/Liu 0.0628 mm2
const Cm = 10.0 # specific capacitance cₘ is a biological constant (around) 10 nF/mm^2

Prinz_conv = Cm / Area

channels = [Prinz.Na(100 * Prinz_conv), Prinz.CaS(6 * Prinz_conv), Prinz.CaT(2.5 * Prinz_conv), Prinz.H(0.01 * Prinz_conv),
    Prinz.Ka(50 * Prinz_conv), Prinz.KCa(5 * Prinz_conv), Prinz.Kdr(100 * Prinz_conv), Prinz.Leak(0.0)]

τCa = 20.0
Ca∞ = 0.05

Dict(:ENa => 50.0, :EH => -20.0, :EK => -80.0, :ELeak => -50.0)


defaults = Dict(Voltage => BasicVoltageDynamics(), Calcium => LiuCalciumDynamics(τCa, Ca∞, 14.96 * 0.0628), Reversal{Calcium} => LiuCaReversalDynamics())

somatic_parameters = Dict(Reversal{Sodium} => -50.0, Reversal{Potassium} => -80.0, Reversal{Leak} => -50.0, Reversal{Proton} => -20.0, Voltage => -60.0, Calcium => 0.05, Reversal{Calcium} => 0.0)

b = BasicNeuron(NoGeometry(Cm), defaults, somatic_parameters, channels, :test_Liu)


neur = b(0)


prob = ODEProblem(neur.sys, [], (0.0, 1000.0), [])

sol = solve(prob, Tsit5())