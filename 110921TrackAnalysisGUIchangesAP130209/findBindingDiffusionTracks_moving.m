function bindingTime = findBindingDiffusionTracks_moving(tracks, params, window, minImmobilizedFrames)
% Stephan Uphoff. 14.09.11
% Anne Plochowietz. 22.12.12
% criteria of selection of binding tracks might be subject to change
% at the moment we are looking at unbound-bound-unbound (binding<unbinding)
% scenarios
% implemented D filter (mean within window) which might be subject to
% change as well

%window = 4; %appData.filterDwindow
%minImmobilizedFrames = 8; %appData.minImmobilizedFrames

DThresh = params.DhistThresh; %0.1;

bindingTime = [];
sigmaNoise = params.sigmaNoise;
pixel = params.pixel;
dT = params.dT;
minStep = params.plotTracksMinSteps;
holdFigure = params.holdFigureCheckbox;
if holdFigure
   figureHandle = params.figureHandle;
   figure(figureHandle);
end
hold all

nMolecules = max(tracks(:,4)); % number of tracks

kk = 1;
uu = 1;

for jj = 1:nMolecules % loop over tracks

   % indexes of this track in the tracks data
   xx = find(tracks(:,4)==jj);

   %plot(tracks(xx,1),tracks(xx,2),...
   %        '-','Color',[0.5 0.5 0.5])

   if numel(xx)>minStep % include only tracks with at least minSteps
       uu = uu + 1;
       squaredStepLength = zeros(numel(xx)-1,1);
       u = [];
       v = [];

       % sum all squared displacement
       for ii = 1:numel(xx)-1

           squaredStepLength(ii) = ((tracks(xx(ii+1),1) - tracks(xx(ii),1))^2 +...
               (tracks(xx(ii+1),2) - tracks(xx(ii),2))^2);

           u(ii) = tracks(xx(ii+1),1) - tracks(xx(ii),1);
           v(ii) = tracks(xx(ii+1),2) - tracks(xx(ii),2);

       end

       % D
       D = squaredStepLength * pixel^2/(4*dT)- (sigmaNoise*pixel)^2/dT * ones(size(squaredStepLength));


       movingD = movingFilterAP(D,window); %movingFilter(D,window,'mean'); wrote own moving filter, giving mean D value within window
       %         D(isnan(movingD)) = [];
       movingD(isnan(movingD)) = []; %remove NaN entries at start and end

       % identify immobilized tracks
       immobilizedIndices = movingD < DThresh(1);
       bindingIndices = find(diff(immobilizedIndices)==1);
       unbindingIndices = find(diff(immobilizedIndices)==-1);
                                                                   % criteria of selection of binding tracks
       if ((sum(immobilizedIndices) > minImmobilizedFrames) ...      % number of immobilized time points > minImmobilizedFrames   
                && (numel(bindingIndices)>1) && (numel(unbindingIndices)>1) ... % number of un/bindingIndices > 1
                && (bindingIndices(1) > 1) ...                                   % first time point is not binding
                && (unbindingIndices(1) < numel(D)-1) )                       % last time point is not unbinding
               %&& bindingIndices < unbindingIndices

           %bindingTime(kk) = (unbindingIndices - bindingIndices)*dT;

           figure; plot((1:length(D))*dT,D,'-k','LineWidth',2); hold all; plot((window+1:1:length(movingD)+window)*dT,movingD,'-r'); Title(['black: raw D, red: filtered D, filter interval: ' num2str(2*window+1) ', Dthresh: ' num2str(DThresh) ', MinBind: ' num2str(minImmobilizedFrames)]);
            
           if holdFigure
            figure(figureHandle); hold on;
           else    
            figure; hold on;
           end
           quiver(tracks(xx(1:end-1),1),tracks(xx(1:end-1),2),u',v',0,'LineWidth',2,'Color','r'); hold on;
           plot(tracks(xx,1),tracks(xx,2),'LineWidth',2,'Color','r'); hold on;
           Title(['filter interval: ' num2str(2*window+1) ', Dthresh: ' num2str(DThresh) ', MinBind: ' num2str(minImmobilizedFrames)]); 
           plot(tracks(xx(1),1),tracks(xx(1),2),'bo','MarkerFaceColor','b','MarkerSize',15); hold on;
           plot(tracks(xx(end),1),tracks(xx(end),2),'ko','MarkerFaceColor','k','MarkerSize',15); hold off;

           kk = kk + 1;

           keyboard    
       end 

   end

end

% figure(figureHandle);
% axis image; % same scale on x and y axis.
% hold off;

% figure;
% hist(bindingTime)

nTraces = kk - 1;
nMolecules;
uu;

end


