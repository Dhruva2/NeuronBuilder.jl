"""
remember to add kwargs
"""


Na(x) = BasicComponents.BasicSingleIonChannel(
    :Na, #name
    Sodium(),
    Dict( #dynamics
        UntrackedQuantity(:m) => Liu.dynamics(:m, Sodium()),
        UntrackedQuantity(:h) => Liu.dynamics(:h, Sodium()),
        Current{Sodium}() => BasicComponents.basic_mh_current(Sodium(), :m => 3, :h => 1),
        Current{Voltage}() => BasicComponents.basic_mh_current(Sodium(), :m => 3, :h => 1; output=Voltage()),
        Current{Sodium}() => BasicComponents.basic_mh_current(Sodium(), :m => 3, :h => 1; output=Sodium())
    ),
    Dict( #defaults
        Conductance{Sodium}() => x,
        UntrackedQuantity(:m) => 0.00146895,
        UntrackedQuantity(:h) => 0.894999
    )
)


CaS(x) = BasicComponents.BasicSingleIonChannel(
    :CaS,
    Calcium(),
    Dict( #dynamics
        UntrackedQuantity(:mS) => Liu.dynamics(:mS, Calcium()),
        UntrackedQuantity(:hS) => Liu.dynamics(:hS, Calcium()),
        Current{Voltage}() => BasicComponents.basic_mh_current(Calcium(), :mS => 3, :hS => 1),
        Current{Calcium}() => BasicComponents.basic_mh_current(Calcium(), :mS => 3, :hS => 1; output=Calcium())
    ),
    Dict( #defaults
        Conductance{Calcium}() => x,
        UntrackedQuantity(:mS) => 0.0344452,
        UntrackedQuantity(:hS) => 0.5
    )
)

CaT(x) = BasicComponents.BasicSingleIonChannel(
    :CaT,
    Calcium(),
    Dict( #dynamics
        UntrackedQuantity(:mT) => Liu.dynamics(:mT, Calcium()),
        UntrackedQuantity(:hT) => Liu.dynamics(:hT, Calcium()),
        Current{Voltage}() => BasicComponents.basic_mh_current(Calcium(), :mT => 3, :hT => 1),
        Current{Calcium}() => BasicComponents.basic_mh_current(Calcium(), :mT => 3, :hT => 1; output=Calcium())
    ),
    Dict( #defaults
        Conductance{Calcium}() => x,
        UntrackedQuantity(:mT) => 0.0102574,
        UntrackedQuantity(:hT) => 0.993774
    )
)

Ka(x) = BasicComponents.BasicSingleIonChannel(
    :Ka,
    Potassium(),
    Dict( #dynamics
        UntrackedQuantity(:mKa) => Liu.dynamics(:mKa, Potassium()),
        UntrackedQuantity(:hKa) => Liu.dynamics(:hKa, Potassium()),
        Current{Voltage}() => BasicComponents.basic_mh_current(Potassium(), :mKa => 3, :hKa => 1),
        Current{Potassium}() => BasicComponents.basic_mh_current(Potassium(), :mKa => 3, :hKa => 1; output=Potassium())
    ),
    Dict( #defaults
        Conductance{Potassium}() => x,
        UntrackedQuantity(:mKa) => 0.0225301,
        UntrackedQuantity(:hKa) => 0.653091
    )
)



KCa(x) = BasicComponents.BasicMultipleIonChannel(
    :KCa,
    [Calcium(), Potassium()], #sensed
    [Potassium()], #actuated
    Dict( #dynamics
        UntrackedQuantity(:mKCa) => Liu.dynamics(:mKCa, Potassium()),
        Current{Voltage}() => BasicComponents.basic_mh_current(Potassium(), :mKCa => 4, :hKCa => 0),
        Current{Potassium}() => BasicComponents.basic_mh_current(Potassium(), :mKCa => 4, :hKCa => 0; output=Potassium())
    ),
    Dict( #defaults
        Conductance{Potassium}() => x,
        UntrackedQuantity(:mKCa) => 0.00122546,
    )
)


Kdr(x) = BasicComponents.BasicSingleIonChannel(
    :Kdr,
    Potassium(),
    Dict( #dynamics
        UntrackedQuantity(:mKdr) => Liu.dynamics(:mKdr, Potassium()),
        Current{Voltage}() => BasicComponents.basic_mh_current(Potassium(), :mKdr => 4, :hKdr => 0),
        Current{Potassium}() => BasicComponents.basic_mh_current(Potassium(), :mKdr => 4, :hKdr => 0; output=Potassium())
    ),
    Dict( #defaults
        Conductance{Potassium}() => x,
        UntrackedQuantity(:mKdr) => 0.0172529,
    )
)

H(x) = BasicComponents.BasicSingleIonChannel(
    :H,
    Proton(),
    Dict( #dynamics
        UntrackedQuantity(:m) => Liu.dynamics(:m, Proton()),
        Current{Voltage}() => BasicComponents.basic_mh_current(Proton(), :m => 1, :h => 0),
        Current{Proton}() => BasicComponents.basic_mh_current(Proton(), :m => 1, :h => 0; output=Proton())
    ),
    Dict( #defaults
        Conductance{Proton}() => x,
        UntrackedQuantity(:m) => 0.158869
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


