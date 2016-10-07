function tracks = trackWithDummyFRET(localizations,trackParams)
% 2013-02-16 Anne Plochowietz
% tracking of position but link data (Sx,Sy,E,etc.) to mol after tracking
% to achieve smFRET tracking

colNum = length(localizations(1,:));
SumPosBefore = zeros(length(localizations(:,1)),1);
% get pos input ready for tracking
pos(:,1) = localizations(:,2);  %x positions
pos(:,2) = localizations(:,3);  %y positions
pos(:,3) = localizations(:,1);  %frame number or time

% add a dummy track outside the field of view (throughout movie) to ensure continuous spacing in time
dummyTrack(:,1:2) = 10000*ones(max(pos(:,3)),2);
dummyTrack(:,3) = 1:max(pos(:,3)); %continuous frame counter

posAndDummy = [pos; dummyTrack];

[BDummy,IXDummy] = sort(posAndDummy(:,3),1); % sort tracks by ascending frame number
posAndDummy = posAndDummy(IXDummy,1:3); % positions without tracking ordered by frame number

[B,IX] = sort(localizations(:,3),1);  %sort pos by ascending frame number
localizations = localizations(IX,1:colNum);
SumPosBefore(:,1)= localizations(:,2)+localizations(:,3); %save pos sum
SumPosBefore(:,2)=localizations(:,2); %save x position
%track using dummy track data
tracks = track(posAndDummy,trackParams.maxDisp,trackParams);
% delete dummy track
tracks(tracks(:,1) == 10000,:) = [];

%sort out what belongs to what and extend tracks
%has to loop over num of frames (3rd column) where mol are still tracked
count = 0;
for ii =1:max(tracks(:,3))
   indFrame = find(tracks(:,3)==ii);
   %loop over positions within frame
   for jj=1:length(indFrame)
   SumPosAfter = sum(tracks(indFrame(jj),1:2));
   IndSumPosBefore = find(SumPosBefore(:,1)==SumPosAfter);
   if numel(IndSumPosBefore) == 1
   tracks(indFrame(jj),5:colNum+1)=localizations(IndSumPosBefore,4:colNum);
   else
   for kk=1:numel(IndSumPosBefore)    
    indXpos = find(localizations(:,2)==SumPosBefore(IndSumPosBefore(kk),2));
    if numel(indXpos) == 1
    tracks(indFrame(jj),5:colNum+1)=localizations(IndSumPosBefore(kk),4:colNum);
    else
    tracks(indFrame(jj),5:colNum+1)=localizations(IndSumPosBefore(kk),4:colNum);
    count = count + numel(IndSumPosBefore);
    end
   end 
   end    
   end
end    


end