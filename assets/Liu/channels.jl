"""
remember to add kwargs
"""


Na(x) = BasicComponents.BasicSingleIonChannel(
    :Na, #name
    Sodium(), #ion
    Dict( #dynamics
        UntrackedQuantity(:m) => Liu.channel_dynamics[:Na][:m],
        UntrackedQuantity(:h) => Liu.channel_dynamics[:Na][:h],
        Current{Voltage}() => BasicComponents.basic_mh_current(Sodium(), :m => 3, :h => 1; output=Voltage()),
        Current{Sodium}() => BasicComponents.basic_mh_current(Sodium(), :m => 3, :h => 1; output=Sodium())
    ),
    Dict( #defaults
        Conductance{Sodium}() => x,
        UntrackedQuantity(:m) => calc_defaults_from(Liu.gating_dynamics[:Na][:m∞], [Voltage()]),
        UntrackedQuantity(:h) => calc_defaults_from(Liu.gating_dynamics[:Na][:h∞], [Voltage()])
    )
)

CaS(x) = BasicComponents.BasicSingleIonChannel(
    :CaS,
    Calcium(),
    Dict( #dynamics
        UntrackedQuantity(:m) => Liu.channel_dynamics[:CaS][:m],
        UntrackedQuantity(:h) => Liu.channel_dynamics[:CaS][:h],
        Current{Voltage}() => BasicComponents.basic_mh_current(Calcium(), :m => 3, :h => 1),
        Current{Calcium}() => BasicComponents.basic_mh_current(Calcium(), :m => 3, :h => 1; output=Calcium())
    ),
    Dict( #defaults
        Conductance{Calcium}() => x,
        UntrackedQuantity(:m) => calc_defaults_from(Liu.gating_dynamics[:CaS][:m∞], [Voltage()]),
        UntrackedQuantity(:h) => calc_defaults_from(Liu.gating_dynamics[:CaS][:h∞], [Voltage()])
    )
)

CaT(x) = BasicComponents.BasicSingleIonChannel(
    :CaT,
    Calcium(),
    Dict( #dynamics
        UntrackedQuantity(:m) => Liu.channel_dynamics[:CaT][:m],
        UntrackedQuantity(:h) => Liu.channel_dynamics[:CaT][:h],
        Current{Voltage}() => BasicComponents.basic_mh_current(Calcium(), :m => 3, :h => 1),
        Current{Calcium}() => BasicComponents.basic_mh_current(Calcium(), :m => 3, :h => 1; output=Calcium())
    ),
    Dict( #defaults
        Conductance{Calcium}() => x,
        UntrackedQuantity(:m) => calc_defaults_from(Liu.gating_dynamics[:CaT][:m∞], [Voltage()]),
        UntrackedQuantity(:h) => calc_defaults_from(Liu.gating_dynamics[:CaT][:h∞], [Voltage()])
    )
)

Ka(x) = BasicComponents.BasicSingleIonChannel(
    :Ka,
    Potassium(),
    Dict( #dynamics
        UntrackedQuantity(:m) => Liu.channel_dynamics[:Ka][:m],
        UntrackedQuantity(:h) => Liu.channel_dynamics[:Ka][:h],
        Current{Voltage}() => BasicComponents.basic_mh_current(Potassium(), :m => 3, :h => 1),
        Current{Potassium}() => BasicComponents.basic_mh_current(Potassium(), :m => 3, :h => 1; output=Potassium())
    ),
    Dict( #defaults
        Conductance{Potassium}() => x,
        UntrackedQuantity(:m) => calc_defaults_from(Liu.gating_dynamics[:Ka][:m∞], [Voltage()]),
        UntrackedQuantity(:h) => calc_defaults_from(Liu.gating_dynamics[:Ka][:h∞], [Voltage()])
    )
)



KCa(x) = BasicComponents.BasicMultipleIonChannel(
    :KCa,
    [Calcium(), Potassium()], #sensed
    [Potassium()], #actuated
    Dict( #dynamics
        UntrackedQuantity(:m) => Liu.channel_dynamics[:KCa][:m],
        Current{Voltage}() => BasicComponents.basic_mh_current(Potassium(), :m => 4, :h => 0),
        Current{Potassium}() => BasicComponents.basic_mh_current(Potassium(), :m => 4, :h => 0; output=Potassium())
    ),
    Dict( #defaults
        Conductance{Potassium}() => x,
        UntrackedQuantity(:m) => calc_defaults_from(Liu.gating_dynamics[:KCa][:m∞], [Voltage(), Calcium()]),
    )
)


Kdr(x) = BasicComponents.BasicSingleIonChannel(
    :Kdr,
    Potassium(),
    Dict( #dynamics
        UntrackedQuantity(:m) => Liu.channel_dynamics[:Kdr][:m],
        Current{Voltage}() => BasicComponents.basic_mh_current(Potassium(), :m => 4, :h => 0),
        Current{Potassium}() => BasicComponents.basic_mh_current(Potassium(), :m => 4, :hKdr => 0; output=Potassium())
    ),
    Dict( #defaults
        Conductance{Potassium}() => x,
        UntrackedQuantity(:m) => calc_defaults_from(Liu.gating_dynamics[:Kdr][:m∞], [Voltage()])
    )
)

H(x) = BasicComponents.BasicSingleIonChannel(
    :H,
    Proton(),
    Dict( #dynamics
        UntrackedQuantity(:m) => Liu.channel_dynamics[:H][:m],
        Current{Voltage}() => BasicComponents.basic_mh_current(Proton(), :m => 1, :h => 0),
        Current{Proton}() => BasicComponents.basic_mh_current(Proton(), :m => 1, :h => 0; output=Proton())
    ),
    Dict( #defaults
        Conductance{Proton}() => x,
        UntrackedQuantity(:m) => calc_defaults_from(Liu.gating_dynamics[:H][:m∞], [Voltage()]),
    )
)


Leak(x) = BasicComponents.BasicSingleIonChannel(
    :Leak,
    Voltage(),
    Dict( #dynamics
        Current{Voltage}() => BasicComponents.basic_mh_current(Voltage(), :m => 0, :h => 0),
    ),
    Dict( #defaults
        Conductance{Voltage}() => x
    )
)


