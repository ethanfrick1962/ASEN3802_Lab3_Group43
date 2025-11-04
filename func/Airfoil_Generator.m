% Name:     AIRFOIL_GENERATOR
% Summary:  Generates 4-digit NACA airfoils.

function [XB, YB] = Airfoil_Generator(MAX_CAMBER, LOC_CAMBER, MAX_THICK)

    arguments (Input)
        MAX_CAMBER
        LOC_CAMBER
        MAX_THICK
    end
    
    arguments (Output)
        XB
        YB
    end

    XB = 0;
    YB = 0;
end