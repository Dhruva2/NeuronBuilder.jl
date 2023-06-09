function get_MTK_from(sys, ::Species) #(or anything)
    
end




"""
Template for stuff


struct LiuCalciumDynmamics <: SpeciesDynamics

end


function (l::LiuCalciumDynmamics)(...)
    a,b, = internal_parameters(l)

end




my_synapse = Glut() # with objectid

n = neuron(owned=channels, actuated_by = my_synapse)
n2 = neuron(...)


function connect(n, n2, my_synapse)
    n2sys.Glut = get_name(my_synapse)
    for el in sensed_by(my_synapse)
        nsys.fluence ~ n2sys.el
    end
end




"""

