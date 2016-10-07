function F = twoConfinedAverageFrameDiffusionSpeciesModel(...
    parameters,lCum,dT,sigmaNoise,bacArea)
%parameters: D

MSD1 = bacArea/3 * ( 1 - exp( -12 * parameters(2) * dT / bacArea ) ) ...
    - 4/3 * parameters(2) * dT + 4 * sigmaNoise^2;

MSD2 = bacArea/3 * ( 1 - exp( -12 * parameters(3) * dT / bacArea ) ) ...
    - 4/3 * parameters(3) * dT + 4 * sigmaNoise^2;

F = parameters(1) * (1 - exp(-lCum.^2/MSD1)) + (1-parameters(1)) * (1 - exp(-lCum.^2/MSD2));

end