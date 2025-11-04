% Name:     AIRFOIL_GENERATOR
% Summary:  Generates 4-digit NACA airfoils.

function [XB, YB] = Airfoil_Generator(MAX_CAMBER, LOC_CAMBER, MAX_THICK)

    arguments (Input)
        MAX_CAMBER
        LOC_CAMBER
        MAX_THICK
    end
% Yt 
    Yt = t * ( ...
    0.2969 * sqrt(x ./ c) ...
  - 0.1260 * (x ./ c) ...
  - 0.3516 * (x ./ c).^2 ...
  + 0.2843 * (x ./ c).^3 ...
  - 0.1036 * (x ./ c).^4 );

   

  yc = zeros(size(x));

%  Yc with x < p*c and x >= p*c
for i = 1:length(x)
    if x(i) < p * c
        yc(i) = m / (p^2) * (x(i)/c) * (2*p - x(i)/c);
    else
        yc(i) = m / ((1 - p)^2) * ((c - x(i))/c) .* (1 + x(i)/c - 2*p);
    end
end



    arguments (Output)
        XB
        YB
    end

    XB = 0;
    YB = 0;
end