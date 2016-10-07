function [ saveTracks, saveLocal ] = ComovCalc( greenLocals, redLocals, greenTracks, redTracks,lobject, distanceThreshold, stepnumber)
%Calculates distance of all tracks red channel, all tracks green channel
%and filters according to stepnumber and distance threshold. Particels that
%fullfill this criterium are considered co-moving. Output co moving tracks,
%percentage co-moving/total number of particles and localizations of
%comoving particles

%initialize output files
saveTracks    = [greenLocals(1:end-5),'_distance',num2str(distanceThreshold),'_steps',num2str(stepnumber),'_startframe', num2str(startframe),'_Comov.tracks'];    
saveLocal  =  [greenLocals(1:end-5),'_distance',num2str(distanceThreshold),'_steps',num2str(stepnumber),'_startframe', num2str(startframe),'_Comov.tracks']; 



fprintf('\nok so far\n');



end

