function findBindingDiffusionTracks(tracks, params)
% Stephan Uphoff. 14.09.11
% code modified by Anne Plochowietz 14/06/12 and 22/12/12:
% use to filter tracks for diffusion coefficient range 
% settings 1 (plot binding tracks) and 2 (plot tracks within D interval):
%   only tracks between Dmax and Dmin (Dthresh  = [Dmax Dmin]) are displayed
% settings 3 (Testing mode):
%   only tracks were half of the tracks is larger than Dmin is displayed
% could be manually set after checking the diffusion histogram
% several subplots are not included in this version

sigmaNoise = params.sigmaNoise;
pixel = params.pixel;
dT = params.dT;
DThresh = params.DThresh; %DThesh(1): lower bound, DThresh(2): upper bounds
minStep = params.plotTracksMinSteps;
plotTracksColour = params.plotTracksColour;
holdFigure = params.holdFigureCheckbox;
if holdFigure
    figureHandle = params.figureHandle;    
end

nMolecules = max(tracks(:,4)); % number of tracks

cmap = colormap(jet); % load colour map

D1 = zeros(nMolecules,1);
D2 = zeros(nMolecules,1);

if holdFigure
    figure(figureHandle);
else
    figure;
end
hold all

kk = 1;

for ll = 1:nMolecules

    % indexes of this track in the tracks data
    xx = find(tracks(:,4)==ll);
    
    if(params.PlotTracksSelection == 2)
    plot(tracks(xx,1),tracks(xx,2),...
            '-','Color',[0.5 0.5 0.5]);
    end
    
    if numel(xx)>minStep
        if(params.PlotTracksSelection == 1)
        plot(tracks(xx,1),tracks(xx,2),...
            '-','Color',[0.5 0.5 0.5]);    
        end
    end     
    
end    

%initialization
xMinimum = min(tracks(:,1));
xMaximum = max(tracks(:,1));
yMinimum = min(tracks(:,2));
yMaximum = max(tracks(:,2));
    
for jj = 1:nMolecules % loop over tracks
    
    % indexes of this track in the tracks data
    xx = find(tracks(:,4)==jj);
  
%     if(params.PlotTracksSelection == 2)
%     plot(tracks(xx,1),tracks(xx,2),...
%             '-','Color',[0.5 0.5 0.5]);
%     end
%     
    if numel(xx)>minStep % include only tracks with at least minSteps
        
%         if(params.PlotTracksSelection == 1)
%         plot(tracks(xx,1),tracks(xx,2),...
%             '-','Color',[0.5 0.5 0.5]);    
%         end
%         
        % compare D between two sections of a track
        MSD1 = 0; % reset MSD1
        MSD2 = 0; % reset MSD2
        
        u = [];
        v = [];
        
        % sum all squared displacement
        for ii = 1:floor(numel(xx)/2)-1
            
            MSD1 = MSD1 + ((tracks(xx(ii+1),1) - tracks(xx(ii),1))^2 +...
                (tracks(xx(ii+1),2) - tracks(xx(ii),2))^2);
            
            u(ii) = tracks(xx(ii+1),1) - tracks(xx(ii),1);
            v(ii) = tracks(xx(ii+1),2) - tracks(xx(ii),2);
            
        end
        
        for ll = floor(numel(xx)/2):numel(xx)-1
            
            MSD2 = MSD2 + ((tracks(xx(ll+1),1) - tracks(xx(ll),1))^2 +...
                (tracks(xx(ll+1),2) - tracks(xx(ll),2))^2);
            
            u(ll) = tracks(xx(ll+1),1) - tracks(xx(ll),1);
            v(ll) = tracks(xx(ll+1),2) - tracks(xx(ll),2);
            
        end        
        
        MSD1 = (MSD1/ii) * pixel^2; % mean square displacement
        MSD2 = (MSD2/ll) * pixel^2; % mean square displacement        
 
        % calculate D from MSD
        D1(kk) = MSD1/(4*dT) - (sigmaNoise*pixel)^2/dT;
        D2(kk) = MSD2/(4*dT) - (sigmaNoise*pixel)^2/dT;

        if (params.TracksSelection == 1)   %one part of the track has to fulfil Dcoeff1>DThresh(1) whereas the other has to fulfil Dcoeff2<DThresh(2)
        
            if (D1(kk)<DThresh(1) && D2(kk)>DThresh(2)) ||...
                    (D2(kk)<DThresh(1) && D1(kk)>DThresh(2))  
                if strcmp(plotTracksColour, 'time')    
                    colorindex = ceil( length(cmap) * tracks(xx(1),4) / nMolecules);
                    quiver(tracks(xx(1:end-1),1),tracks(xx(1:end-1),2),u',v',0,'Color',cmap(colorindex,:),'LineWidth',1)
                elseif strcmp(plotTracksColour,'random')
                    quiver(tracks(xx(1:end-1),1),tracks(xx(1:end-1),2),u',v',0,'Color',rand(1,3),'LineWidth',1)
                elseif strcmp(plotTracksColour,'timetime')
                    if (numel(xx) > minStep)    %fixed number has to be inserted here to keep same color-time-scale
                        keyboard;
                    end    
                    for ll=1:numel(xx)
                        colorindex(ll) = ceil(length(cmap) *ll/numel(xx)); %*tracks(xx(ll),4) /(subPlotId*tracksPerSubPlot) %subPlot not included
                    end
                    for mm=1:(numel(xx)-1)
                        quiver(tracks(xx(mm:mm+1),1),tracks(xx(mm:mm+1),2),u(mm:mm+1)',v(mm:mm+1)',0,'Color',cmap(colorindex(mm),:),'LineWidth',1), axis([xMinimum xMaximum yMinimum yMaximum]), hold on;
                    end
                    clear colorindex
                    %keyboard;
                else
                    quiver(tracks(xx(1:end-1),1),tracks(xx(1:end-1),2),u',v',0,'Color',plotTracksColour,'LineWidth',1)
                end    
            end
        
        elseif (params.TracksSelection ==2)     %both parts of the track fulfil: Dcoeffs>DThresh(1), Dcoeffs<DThresh(2)
        
            if (D1(kk)>DThresh(1) && D2(kk)> DThresh(1)) &&...
                    (D2(kk)<DThresh(2) && D1(kk)<DThresh(2))   
                if strcmp(plotTracksColour, 'time')    
                    colorindex = ceil( length(cmap) * tracks(xx(1),4) / nMolecules);
                    quiver(tracks(xx(1:end-1),1),tracks(xx(1:end-1),2),u',v',0,'Color',cmap(colorindex,:),'LineWidth',1)
                elseif strcmp(plotTracksColour,'random')
                    quiver(tracks(xx(1:end-1),1),tracks(xx(1:end-1),2),u',v',0,'Color',rand(1,3),'LineWidth',1)
                elseif strcmp(plotTracksColour,'timetime')
                    if (numel(xx) > minStep)    %fixed number has to be inserted here to keep same color-time-scale
                        keyboard;
                    end    
                    for ll=1:numel(xx)
                        colorindex(ll) = ceil(length(cmap) *ll/numel(xx)); %*tracks(xx(ll),4) /(subPlotId*tracksPerSubPlot) %subPlot not included
                    end
                    for mm=1:(numel(xx)-1)
                        quiver(tracks(xx(mm:mm+1),1),tracks(xx(mm:mm+1),2),u(mm:mm+1)',v(mm:mm+1)',0,'Color',cmap(colorindex(mm),:),'LineWidth',1), axis([xMinimum xMaximum yMinimum yMaximum]), hold on;
                    end
                    clear colorindex
                    %keyboard;
                else
                    quiver(tracks(xx(1:end-1),1),tracks(xx(1:end-1),2),u',v',0,'Color',plotTracksColour,'LineWidth',1)
                end  
            end        
            
         elseif (params.TracksSelection ==3)    %one part of the track has to fulfil Dcoeff>DThresh(1), no upper bound needed (binding and unbinding events might be possible)
        
            if (D1(kk)>DThresh(1)) ||...
               (D2(kk)>DThresh(1))    
                if strcmp(plotTracksColour, 'time')    
                    colorindex = ceil( length(cmap) * tracks(xx(1),4) / nMolecules);
                    quiver(tracks(xx(1:end-1),1),tracks(xx(1:end-1),2),u',v',0,'Color',cmap(colorindex,:),'LineWidth',1)
                    keyboard;
                elseif strcmp(plotTracksColour,'random')
                    quiver(tracks(xx(1:end-1),1),tracks(xx(1:end-1),2),u',v',0,'Color',rand(1,3),'LineWidth',1)
                elseif strcmp(plotTracksColour,'timetime')
                    if (numel(xx) > minStep)    %fixed number has to be inserted here to keep same color-time-scale
                        keyboard;
                    end    
                    for ll=1:numel(xx)
                        colorindex(ll) = ceil(length(cmap) *ll/numel(xx)); %*tracks(xx(ll),4) /(subPlotId*tracksPerSubPlot) %subPlot not included
                    end
                    for mm=1:(numel(xx)-1)
                        quiver(tracks(xx(mm:mm+1),1),tracks(xx(mm:mm+1),2),u(mm:mm+1)',v(mm:mm+1)',0,'Color',cmap(colorindex(mm),:),'LineWidth',1), axis([xMinimum xMaximum yMinimum yMaximum]), hold on;
                    end
                    clear colorindex
                    %keyboard;
                else
                    quiver(tracks(xx(1:end-1),1),tracks(xx(1:end-1),2),u',v',0,'Color',plotTracksColour,'LineWidth',1)
                end  
            end        
                        
            
        end
            
        kk = kk+1;  
        
    end
    
end

D1(kk:end) = [];
D2(kk:end) = [];

axis image; % same scale on x and y axis.
hold off;

end