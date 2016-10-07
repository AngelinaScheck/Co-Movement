function combinedTracks = combineTracks

combinedTracks = [];

[tracksFilename, tracksPathname] = ...
    uigetfile('*.tracks', 'MatFile data:','MultiSelect', 'on');

if ~(isnumeric(tracksFilename)&&tracksFilename==0) %check the user has not pressed cancel

info = whos('tracksFilename');
if strcmp(info.class,'char')
    appData.nFiles = 1;
else
    appData.nFiles = numel(tracksFilename);
end

for ii = 1:appData.nFiles
    
    if appData.nFiles == 1
        loadname  = [tracksPathname tracksFilename];
    else
        loadname = [tracksPathname tracksFilename{1,ii}];
    end
    
    newTracks = importdata(loadname);
    
    if ii>1
    newTracks(:,4) = newTracks(:,4) + max(combinedTracks(:,4));
    end
    
    combinedTracks = [combinedTracks; newTracks];
    
end

[filename, pathname] = uiputfile('.tracks');
save([pathname filename],'combinedTracks');

else
    
    return
    
end

end