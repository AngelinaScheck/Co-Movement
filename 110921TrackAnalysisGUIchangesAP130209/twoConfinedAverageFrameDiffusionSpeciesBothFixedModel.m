function F = twoConfinedAverageFrameDiffusionSpeciesBothFixedModel(...
    parameters,lCum,dT,sigmaNoise,bacArea,D1fixed,D2fixed)
%parameters: a

MSD1 = bacArea/3 * ( 1 - exp( -12 * D1fixed * dT / bacArea ) ) ...
    - 4/3 * D1fixed * dT + 4 * sigmaNoise^2;

MSD2 = bacArea/3 * ( 1 - exp( -12 * D2fixed * dT / bacArea ) ) ...
    - 4/3 * D2fixed * dT + 4 * sigmaNoise^2;

F = parameters(1) * (1 - exp(-lCum.^2/MSD1)) + (1-parameters(1)) * (1 - exp(-lCum.^2/MSD2));

end