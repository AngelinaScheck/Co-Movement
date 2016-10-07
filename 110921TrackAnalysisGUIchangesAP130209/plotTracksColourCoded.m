function plotTracksColourCoded(tracks,params)
% function plotTracksColourCoded(tracks,minStep,plotTracksColour) plots all tracks that have
% at least minStep steps. Tracks is the output of the tracking software:
% tracks(:,1) = x coordinates
% tracks(:,2) = y coordinates
% tracks(:,3) = frame
% tracks(:,4) = track id
%
% plotTracksColour defines the one letter colour of the tracks.
% plotTracksColour = 'random' plots tracks in random colours
% plotTracksColour = 'time' plot tracks coloured coded by order of appearance in the
% movie. Tracks early in the movie are plotted blue, later tracks are
% plotted red. 
%
% Stephan Uphoff. 21.09.11
%
% plotTracksColour = 'timetime' plots single tracks with time-dependent
% color
% params.numberOfSubplots: introduce number of subplots to get idea of time-dependence of tracks
% Anne Plochowietz. 13/06/12

numberOfSubplots = params.numberOfSubplots;

if (numberOfSubplots == 1)
    numberOfColSubplots = 1; 
else
    numberOfColSubplots = 2;
end    
minStep = params.plotTracksMinSteps;
plotTracksColour = params.plotTracksColour;
holdFigure = params.holdFigureCheckbox;
if holdFigure
    figureHandle = params.figureHandle;   
end

nMolecules = max(tracks(:,4)); % number of tracks

cmap = colormap(jet); % load colour map

if holdFigure
    figure(figureHandle);
    numberOfSubplots = 1;
else
    figure;
end
hold all

tracksPerSubPlot = ceil(nMolecules/numberOfSubplots);

%initialization
xMinimum = min(tracks(:,1));
xMaximum = max(tracks(:,1));
yMinimum = min(tracks(:,2));
yMaximum = max(tracks(:,2));

kk = 1; %not used at the moment

for jj = 1:nMolecules % loop over tracks
    
    subPlotId = ceil(jj/tracksPerSubPlot);
    
    % indexes of this track in the tracks data
    xx = find(tracks(:,4)==jj);
    
    if numel(xx)>minStep % include only tracks with at least minSteps
        
       
        if strcmp(plotTracksColour, 'time')
        
        % choose the colour by the track id
        colorindex = ceil( length(cmap) * tracks(xx(1),4) /(subPlotId * tracksPerSubPlot) );  %tracks(xx(1),4)/nMolecules

        % plot the track in the chosen colour
        subplot(ceil(numberOfSubplots/2),numberOfColSubplots,subPlotId), plot(tracks(xx,1),tracks(xx,2),...
            '.-','Color',cmap(colorindex,:), 'MarkerSize',3.0), axis([xMinimum xMaximum yMinimum yMaximum]), hold on;%,ylim([yMinimum yMaximum]), hold on;
        %axis image;
        %keyboard;
        elseif strcmp(plotTracksColour,'random') 
        subplot(ceil(numberOfSubplots/2),2,subPlotId), plot(tracks(xx,1),tracks(xx,2),...
            '.-','Color',rand(1,3), 'MarkerSize',3.0), axis([xMinimum xMaximum yMinimum yMaximum]), hold on;            
        axis image;
        
        elseif strcmp(plotTracksColour,'timetime')
        
        %choose colour dependent on track length of same track ID    
            for ll=1:numel(xx) %419 or 628
                colorindex(ll) = ceil(length(cmap) *ll/numel(xx)); %*tracks(xx(ll),4) /(subPlotId*tracksPerSubPlot)
            end
            %keyboard;
            % plot the track in the chosen colour
        subplot(ceil(numberOfSubplots/2),numberOfColSubplots,subPlotId), hold on;
        for mm=1:(numel(xx)-1)
        plot(tracks(xx((mm):(mm+1)),1),tracks(xx((mm):(mm+1)),2),...
            '.-','Color',cmap(colorindex(mm),:), 'MarkerSize',3.0), axis([xMinimum xMaximum yMinimum yMaximum]), hold on;%,ylim([yMinimum yMaximum]), hold on;
        end
        clear colorindex
        else
        subplot(ceil(numberOfSubplots/2),2,subPlotId), plot(tracks(xx,1),tracks(xx,2),...
            '.-','Color',plotTracksColour, 'MarkerSize',3.0), hold on;
        %axis image;
         
        end
        
    end
    kk = kk + 1;
end

if(numberOfSubplots == 1)
axis image; % same scale on x and y axis.
end
hold off;

end