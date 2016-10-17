function [ saveTracks, saveLocal ] = ComovCalc( greenLocals, redLocals, greenTracks, redTracks,lobject, distanceThreshold, stepnumber)
%Calculates distance of all tracks red channel, all tracks green channel
%and filters according to stepnumber and distance threshold. Particels that
%fullfill this criterium are considered co-moving. Output co moving tracks,
%percentage co-moving/total number of particles and localizations of
%comoving particles

%initialize output files
saveTracks    = [greenLocals(1:end-5),'_distance',num2str(distanceThreshold),'_steps',num2str(stepnumber),'_Comov.tracks'];    
saveLocal  =  [greenLocals(1:end-5),'_distance',num2str(distanceThreshold),'_steps',num2str(stepnumber),'_Comov.tracks']; 

%Load the tracks, Format: x // y // frame // track-ID
% if tracking in FRET mode, than the tracks are hidden in the struct
data=importdata(greenTracks);
tracksG=data.tracksR; %Tracking is done for both channels independently, but program tracks in red channel, so the red channel values are used, even if tracked for green channel

%Use struct to store filtered tracks
%Track-ID// Number of Steps in Track // x-position // y-position // frames
%// ID of corresponding tracks in other channel (should be one ID only)
filTGreen = struct('id',{},'steps',{},'x',{},'y',{},'frames',{},'partner',{});
filTRed = struct('id',{},'steps',{},'x',{},'y',{},'frames',{},'partner',{});

filtGreenList=[];
mapGreen = containers.Map
mapRed= containers.Map
xy = struct('x', {}, 'y', {});

%Matching IDs will be stored in a table, distances can be checked by looking IDs up in mapRed and mapGreen data structures
idMatch = zeros(0,2); %//greenID //matching red ID



%Filter and Copy to Map
% Aim: loop through every Tracklist only once to save time is not possible
% because of the memory function, use internal functions like find
nMolecules = max(tracksG(:,4)); % number of tracks
for i = 1:nMolecules % loop over tracks
    
    % indexes of track in the tracks file
    index = find(tracksG(:,4)==i);
    steplength=numel(index);
    if steplength > stepnumber % copy tracks with more then stepnumber steps into map filestructure
        x= []; %x positions of the track by time point
        y= []; %y positions of the track by time point
        frame = []; % corresponting framenumbers (should be sorted)
        for j=1:size(index)
            x(end+1) = tracksG(index(j), 1);
            y(end+1) = tracksG(index(j), 2);
            frame(end+1) = tracksG(index(j), 3);
        end
        filTGreen.id=i;
        filtGreen.steps=steplength;
        filtGreen.x=x;
        filtGreen.y=y;
        filtGreen.frame=frame;
		
		xy.x=x;
		xy.y=y;
		mapGreen(id)=xy;
        
        filtGreenList(end+1) = filtGreen;
    end
end

%same for red
data=importdata(redTracks);
tracksR=data.tracksR;
nMoleculesR = max(tracksR(:,4));

for r = 1:nMoleculesR % loop over tracks
    
    % indexes of track in the tracks file
    index = find(tracksR(:,4)==r);
    steplength=numel(index);
    if steplength > stepnumber % compare with step-threshold
        
		
        %Immediatly compare with filtered green Tracks. Save only, if there
        %is a match in the green tracks
        %loop over green tracks
        for g= 1:size(filtGreenList)
            %loop over positions in green track
            for p=1:size(filtGreenList(g).x)
                %Compare distance
                if sqrt((diff(filtGreenList(g).x(p), tracksR(1,r)))^2 + (diff(filtGreenList(g).y(p), tracksR(1,r)))^2 ) <= distanceThreshold
					%copy the redTrack
				
				
                end
           end
        end
        
        
        x= []; %x positions of the track by time point
        y= []; %y positions of the track by time point
        frame = []; % correspontind framenumbers (should be sorted)
        for j=1:size(index)
            x(end+1) = tracksG(index(j), 1);
            y(end+1) = tracksG(index(j), 2);
            frame(end+1) = tracksG(index(j), 3);
        end
        filTGreen.id=i;
        filtGreen.steps=steplength;
        filtGreen.x=x;
        filtGreen.y=y;
        filtGreen.frame=frame;
    end
end






% %Initialize the Counters
% i=1; %counter in loop1
% 
% %Load the tracks
% data=importdata(greenTracks);
% tracks=data.tracksG;
% 
% l=numel(tracks(:,4));
% while i<l
%     steplength=0; %actual number of steps per track
%     id = tracks(i,4);
%     firstOc= i;% first Occurance of an ID
%     
%     if (tracks(i+1, 4) == id)
%         steplength=steplength + 1;
%     else
%         % if the threshold for the steplength is already reached, save the
%         % track
%         if steplength > stepnumber
%             x = []; %x-position of one time point of the track
%             y = []; % y-position
%             frame = []; %frame
%             %loop through track
%             for j=firstOc:1:(firstOc + steplength)
%                 x = [];
%                 y = [];
%                 frame = [];
%                 x(end +1) = tracks(j, 1);
%             end
%         end
%     end
%     
%     fprintf('\nok so far\n');
%     i=i+steplength;
%     
% end




end

