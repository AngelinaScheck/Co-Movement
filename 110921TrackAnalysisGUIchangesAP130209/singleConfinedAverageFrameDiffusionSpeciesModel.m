function F = singleConfinedAverageFrameDiffusionSpeciesModel(...
    parameters,lCum,dT,sigmaNoise,bacArea)
%parameters: D

MSD = bacArea/3 * ( 1 - exp( -12 * parameters * dT / bacArea ) ) ...
    - 4/3 * parameters * dT + 4 * sigmaNoise^2;

F = 1 - exp(-lCum.^2/MSD);

end