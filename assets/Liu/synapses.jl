chol(x; name=:Chol) = BasicSynapse(
    name,
    Choline(),
    Dict(
        UntrackedQuantity(:s) => Liu.synaptic_channel_dynamics[:s],
        Current{Voltage}() => BasicComponents.basic_mh_current(Choline(), :s => 1, :h => 0; output=Voltage())
    ),
    Dict(
        UntrackedQuantity(:k₋) => 0.025,
        UntrackedQuantity(:Vth) => -35.0,
        UntrackedQuantity(:δ) => 5.0,
        UntrackedQuantity(:s) => 0.0,
        Reversal{Choline}() => -70.0,
        Conductance{Choline}() => x
    )

)