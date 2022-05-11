"""
abstract type Geometry end
"""

struct NoGeometry{C} <: Geometry
    capacitance::C
end
capacitance(g::NoGeometry) = g.capacitance



######### Geometries ###########
# is_geometric(::Compartment{T}) = T
# is_geometric(::Component) = false
"""
Every geometry needs methods for:
    capacitance()
    maybe reversals as dictionary?
    initial_conditions(f::Species)
"""
