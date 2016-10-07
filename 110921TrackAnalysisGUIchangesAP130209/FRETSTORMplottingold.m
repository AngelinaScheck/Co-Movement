%% FRET STORM Reconstruction. Stephan Uphoff. 09.07.11
% 1. run localization software i.e. FRETSTORM('*storm.fits',20,7);
% 2. load localization data into workspace
% 3. run this script

%%

AAxpos = data(3:3:end,2);
AAypos = data(3:3:end,3);
E = data(3:3:end,11);
S = data(3:3:end,12);

AA = data(3:3:end,4);
DD = data(1:3:end,4);
DA = data(2:3:end,4);

DDxpos = data(1:3:end,2);
DDypos = data(1:3:end,3);

%%

filter = find(S>0.3 & S<0.5);

figure; hist(E(filter))

%%
figure; plot(E,S,'.'); xlim([0 1]); ylim([0 1]);

%figure; plot(E(filter),S(filter),'.'); xlim([0 1]); ylim([0 1]);

%%

figure;

plot(AAxpos, AAypos, 'rx')

%%
%specific for analysed file
fileName = 'Expt1_hiFRETnorm_EPnorecov_5msalt_ALEX_G1.3mW_R0.65mW_nTIRF_192_260_40.fits';
alternationPeriod = 2;
imageLim = [257 512 1 69; 1 256 1 69];
firstGreenFrame = 2;

tirfIm = TirfImage(fileName,firstGreenFrame, imageLim, alternationPeriod);

green_stack = getGreenStack(tirfIm);
red_stack = getRedStack(tirfIm);

figure;

imageI = getFrame(red_stack,2);
%specific for analysed file
subim=imageI(1:69,1:256);
imavg_R=cast(subim,'double');
imshow(imageI,[min(imavg_R(:)), max(imavg_R(:))]);
set(gca,'YDir','reverse')

%% 
figure; 

%imshow(STD_17_AA) 
imshow(imageI)
set(gca,'YDir','reverse')

hold all


for jj = 15000:length(AAxpos)-1
    
    ii = length(AAxpos) - jj;
    
    if S(ii)>0.35 && S(ii)<0.45
        
        if E(ii)>=0.5 && E(ii)<0.8
            
            colorcode = [1 0 0];
            
            plot(AAxpos(ii),AAypos(ii),'.','Color',colorcode,'MarkerSize',5.0)
            
            %         elseif E(ii)>=0.5 && E(ii)<0.6 % %
            colorcode = [0.5 0.5 0]; % %
            plot(AAxpos(ii),AAypos(ii),'.','Color',colorcode,'MarkerSize',5.0) % %         elseif E(ii)>=0.7 && E(ii)<0.8
            % %             colorcode = [0.8 0.8 0]; % %
            plot(AAxpos(ii),AAypos(ii),'.','Color',colorcode,'MarkerSize',5.0)
            
            
        else
            
            colorcode = [1 1 0];
            
            plot(AAxpos(ii),AAypos(ii),'.','Color',colorcode,'MarkerSize',5.0)
            
        end
        
    end
    
end

hold off

% %% plot autodetect twotone positions too
% 
% for jj = 1:length(intensities)
%     
%     plot(intensities(1,jj).AAadetPos(1),intensities(1,jj).AAadetPos(2),'bo','MarkerSize',10.0)
%     
% end
% 
% hold off

%% prepare data


% bring data into the right format for single particle tracking
pos(:,1) = AAxpos;
pos(:,2) = AAypos;
pos(:,3) = E;
pos(:,4) = S;
pos(:,5) = data(3:3:end,1)/2; %frame (accounting for alternation)

%% tracking

trackWindow = 1.0; % tracking window in pixels

param.mem = 1;
param.dim = 2;
param.good = 1;
param.quiet = 0;

% run tracking function
tracks = track(pos, trackWindow, param);

% number of molecules in field of view
NumTrajectories = max(tracks(:,6));

%%

figure;

%imshow(STD_17_AA)
imshow(imageI)
set(gca,'YDir','reverse')

hold all

for ii = 1:NumTrajectories
    
    xx = find(tracks(:,6)==ii);
    % don't include first and last data point in each track for averaging
    % because they have messed up E and S values
    avPos = mean(tracks(xx(2:end-1),1:4),1);
    
    if length(xx)>10 && avPos(4)>0.3 && avPos(4)<0.5
        
        if avPos(3)>=0.5 && avPos(3)<0.8
            
            colorcode = [1 1 0];
            
            plot(avPos(1),avPos(2),'.','Color',colorcode, 'MarkerSize',5.0)
            
            
        else
            
            colorcode = [1 0 0];
            
            plot(avPos(1),avPos(2),'.','Color',colorcode, 'MarkerSize',5.0)
            
        end
        
    end
    
end

hold off;

%%


% bring data into the right format for single particle tracking
pos(:,1) = AAxpos;
pos(:,2) = AAypos;
pos(:,3) = E;
pos(:,4) = S;
pos(:,5) = DD;
pos(:,6) = DA;
pos(:,7) = AA;
pos(:,8) = data(3:3:end,1)/2; %frame (accounting for alternation)

% tracking

trackWindow = 1.0; % tracking window in pixels

param.mem = 1;
param.dim = 2;
param.good = 1;
param.quiet = 0;

% run tracking function
tracks = track(pos, trackWindow, param);

% number of molecules in field of view
NumTrajectories = max(tracks(:,9));
%%
figure;

%imshow(STD_17)
imshow(imageI)
set(gca,'YDir','reverse')

hold all

for ii = 1:NumTrajectories
    
    xx = find(tracks(:,9)==ii);
    avPos = mean(tracks(xx(2:end-1),1:7),1);
    
    if length(xx)>3 && avPos(4)>0 && avPos(4)<1 && avPos(5)+avPos(6)>1
        
        if avPos(3)>=0.5 && avPos(3)<0.8
            
            colorcode = [1 1 0];
            
            plot(avPos(1),avPos(2),'.','Color',colorcode, 'MarkerSize',5.0)
            
            
        else
            
            colorcode = [1 0 0];
            
            plot(avPos(1),avPos(2),'.','Color',colorcode, 'MarkerSize',5.0)
            
        end
        
    end
    
end

hold off;


%% plot localizations each in a different colour

ColorCoding=rand(NumTrajectories,3);

figure;
hold all;

for ii=1:length(tracks)
    
    %plot(tracks(ii,1), tracks(ii,2), '.', 'Color', ColorCoding(tracks(ii,4),:), 'MarkerSize', 5.0)
    plot(tracks(ii,1),tracks(ii,2),'b.','MarkerSize',5.0)
    % or alternatively all in the same color
    % plot(tracks(ii,1), tracks(ii,2), 'b.', 'MarkerSize', 5.0)
    
end

hold off;

%% plot E map colour coded

figure;

%imshow(STD_17_AA)
imshow(imageI)
set(gca,'YDir','reverse')
%freezeColors;

cmap = colormap(jet); % load colour map

E = zeros(length(tracks),1);
kk = 1;
hold all

for ii = 1:NumTrajectories
    
    xx = find(tracks(:,6)==ii);
    
    if length(xx) > 3
    
    % don't include first and last data point in each track for averaging
    % because they have messed up E and S values
    avPos = mean(tracks(xx(2:end-1),1:4),1);
    
    if  avPos(4)>0.3 && avPos(4)<0.5 && avPos(3)>0 && avPos(1)<235 && avPos(1)>20 && ...
            avPos(2)>15 && avPos(2)<285
        
        
        E(kk) = avPos(3);
        
        % choose the colour according to E value
        colorindex(kk) = ceil( length(cmap(:,1)) * avPos(3) );
        
        % plot the track in the chosen colour
        plot(avPos(1),avPos(2),'.','Color',cmap(colorindex(kk),:), 'MarkerSize',15.0)
        
        kk = kk + 1;
        
    end
    
    end
    
end

% delete unused rows
E(kk:end) = [];

hold off;

figure; hist(E)



