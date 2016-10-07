function FRETSTORMplotting(FLloadname, FRETSTORMFilename, tif, alttype, timeres, plotTracksMinSteps)

% FRET STORM Reconstruction. Stephan Uphoff. 09.07.11
% Edit by RC 31/7/12:
% Integrated into gaussStormGUI, passes in fluorescence movie filename (FL), 
% analysis filename, whether it's tif or not, ALEX or not, and the timeres. 
% Put in sigma limits (to get rid of badly fitted psfs).  
% Can be used with data from ROI selector.
%
% 09/02/13 Anne Plochowietz: 
% Introducing differentiation between filtering, tracking and plotting the analyzed
% data, give output data that can be used by already implemented
% functions, e.g. plot histogram or plot tracks
% image read-in made setup independent
% implemented to draw circles on top of movie/ overlay movie

%08/16 Angelina Scheck: some bugs in ALEX mode fixed

%import data
data = importdata(FRETSTORMFilename);

%check if struct. A struct is produced for whole movie analysis, but just a
%matrix is produced for a ROI
wholemovie = isstruct(data);
if (wholemovie==1)
    colheaders = data.textdata; data = data.data; 
end
    

%Get image limits
if (tif == 0)
    ImageInfo = fits_read_header(FLloadname);
    %Halfdome imageLim
    imageLim = [257 512 1 ImageInfo.NAXIS2; 1 256 1 ImageInfo.NAXIS2];
    %Zugspitze imageLim
    %imageLim = [1 171 1 ImageInfo.NAXIS2; 172 342 1 ImageInfo.NAXIS2];
else
    firstImage = imread(FLloadname,1);
    %Halfdome imageLim
    imageLim = [257 512 1 length(firstImage(:,1)); 1 256 1 length(firstImage(:,1))];
    %Zugspitze imageLim
    %imageLim = [1 171 1 length(firstImage(:,1)); 172 342 1 length(firstImage(:,1))];
end   


%Sigma limits, change if you deem appropriate for your molecules. For now
%these seem reasonable limits.

sigmax_min = 0.1; sigmax_max = 2;
sigmay_min = 0.1; sigmay_max = 2;
    
if (alttype == 1)
    %AS: If acquisition is stopped manually, it might happen, that not
    %equal # of frames are recorded for DD, DA ans AA -> reduce data to
    %last complete set
    rightSize= (fix(length(data(:,1))/2))*2;
    data=data(1:rightSize,:);
    %dig out individual data - CW case    
    DDxpos = data(1:2:end,2);
    DDypos = data(1:2:end,3);    
    DAxpos = data(2:2:end,2);
    DAypos = data(2:2:end,3);
    frames = data(2:2:end,1);

    DD = data(1:2:end,4);
    DA = data(2:2:end,4);

    E = data(2:2:end,11);
    
    LINK = data(2:2:end,12); %AS for combined data

    S_DD_X = data(1:2:end,7);
    S_DD_Y = data(1:2:end,8);

    S_DA_X = data(2:2:end,7);
    S_DA_Y = data(2:2:end,8);
    
    figure;
    subplot(4,1,[1 3])
    plot(DAxpos, DAypos, 'r+'); xlim([0 256]); ylim([0 imageLim(1,4)]); hold on
    subplot(4,1,4)
    hist(E,0:0.02:1); xlim([0 1]); hold on; title('Before filtering');
%     
    figure
    subplot(2,2,1)
    hist(S_DD_X,0:0.05:5); xlim([0 5]); title('S-DD-X'); hold on
    subplot(2,2,2)
    hist(S_DD_Y,0:0.05:5); xlim([0 5]); title('S-DD-Y'); hold on
    subplot(2,2,3)
    hist(S_DA_X,0:0.05:5); xlim([0 5]); title('S-DA-X'); hold on
    subplot(2,2,4)
    hist(S_DA_Y,0:0.05:5); xlim([0 5]); title('S-DA-Y'); hold on

    
    %chop out inner frame due to localisations happening at the edges (see
    framesize = 13;
    indrmv = find(DAxpos < framesize | DAxpos > (imageLim(2,2)-framesize) | DAypos < framesize | DAypos > (imageLim(1,4)-framesize));   
    
    DAxpos(indrmv) = [];
    DAypos(indrmv) = [];
    
    DD(indrmv) = [];
    DA(indrmv) = [];
    
    E(indrmv) = [];
    LINK(indrmv) = []; %AS for combined data
    
    S_DD_X(indrmv) = [];
    S_DD_Y(indrmv) = [];
    
    S_DA_X(indrmv) = [];
    S_DA_Y(indrmv) = [];
    
    frames(indrmv) = [];
    
    figure
    subplot(2,2,1)
    hist(S_DD_X,0:0.05:5); xlim([0 5]); title('S-DD-X'); hold on
    subplot(2,2,2)
    hist(S_DD_Y,0:0.05:5); xlim([0 5]); title('S-DD-Y'); hold on
    subplot(2,2,3)
    hist(S_DA_X,0:0.05:5); xlim([0 5]); title('S-DA-X'); hold on
    subplot(2,2,4)
    hist(S_DA_Y,0:0.05:5); xlim([0 5]); title('S-DA-Y'); hold on
    
    %chop out fitted psfs that don't satisfy sigma criteria
    indrmv_S_DD = find(S_DD_X < sigmax_min | S_DD_X > sigmax_max | S_DD_Y < sigmay_min | S_DD_Y > sigmay_max);
    indrmv_S_DA = find(S_DA_X < sigmax_min | S_DA_X > sigmax_max | S_DA_Y < sigmay_min | S_DA_Y > sigmay_max);
    indrmv_S = [indrmv_S_DD; indrmv_S_DA]; 
    
    DAxpos(indrmv_S) = [];
    DAypos(indrmv_S) = [];
    
    DD(indrmv_S) = [];
    DA(indrmv_S) = [];
    
    E(indrmv_S) = [];
    
    LINK(indrmv_S) = []; %AS for combined data
    
    S_DD_X(indrmv_S) = [];
    S_DD_Y(indrmv_S) = [];
    
    S_DA_X(indrmv_S) = [];
    S_DA_Y(indrmv_S) = [];
    
    frames (indrmv_S) = [];
    %times = frames * timeres;
    
    figure; hist(E,0:0.02:1); xlim([0 1]); title('After filtering');
%     subplot(4,1,[1 3])
%     plot(DAxpos, DAypos, 'r+'); xlim([0 256]); ylim([0 imageLim(1,4)]); hold on
    %efficiency not by roi
    
%      AS: In case data from multiple movies should be combined: save
%     filename of the out file, ROI number, E and S in .comb file for later
%     processing (this part was quickly added, evtually revise later)
    fid = fopen('combinedCells.comb','a');

    for c=1:length(E)
        fprintf(fid,'%s\t',FRETSTORMFilename); %Filename
        fprintf(fid,'%d\t',LINK(c)); %ROI number
        fprintf(fid,'%6.3f\n',E(c)); %E
    end
    
    close(fid);
    
    fprintf('done')

%     %
%     
%         %Bring data into the right format for single particle tracking
%         pos(:,1) = DAxpos;
%         pos(:,2) = DAypos;
%         pos(:,3) = E;
%         pos(:,4) = DD;
%         pos(:,5) = DA;
%         pos(:,6) = frames;
% 
%         %now sort rows based on frame number (needed for track.m)
%         pos = sortrows(pos,6);
%         
%         
%         %Tracking
%         
%         trackWindow = 5.0; % tracking window in pixels
%         
%         param.mem = 3;
%         param.dim = 2;
%         param.good = 1;
%         param.quiet = 0;
%         
%         %Run tracking function
%         tracks = track(pos, trackWindow, param);
%         
%         %don't include first and last tracked point, as it has screwed up
%         %values
%         tracks(1,:) = [];
%         tracks(end,:) = [];
        
%         % number of molecules in field of view
%         NumTrajectories = max(tracks(:,7));
%         savename    = [FRETSTORMFilename '_tracks'];
%         
%         %one off case!! - if only 1 psf in ROI
%         max_time = (size(tracks,1)-1) * timeres;
%         times = 0:timeres:max_time;
        
%         saveFileHeader = 'DAxpos\tDAypos\tE\tDD\tDA\tframe_no\tparticleIP\n';
%         fid = fopen(savename,'w');
%         fprintf(fid, saveFileHeader);
%         fprintf(fid, '%6.3f %6.3f %6.3f %6.3f %6.3f %d %d\n',tracks');
%         fclose(fid);
        
%         %temporary plot - plotting all data obtained as one time trace -
%         %alter later!!  (this was because I isolated 1 psf using the ROI
%         %selector)...
%         figure
%         subplot(2,1,1,'Fontsize',20)
%         plot(times, tracks(:,4),'g-','LineWidth',2); hold on
%         plot(times, tracks(:,5),'r-','LineWidth',2); title('Raw intensities','Fontsize',20);xlabel('Time(s)','Fontsize',20), ylabel('Counts','Fontsize',20);
%         subplot(2,1,2,'Fontsize',20)
%         plot(times, tracks(:,3),'b-','LineWidth',2); title('E_{raw}','Fontsize',20); xlabel('Time (s)','Fontsize',20);
        

%             for jj = 1:NumTrajectories % loop over tracks
% 
%                 indexes of this track in the tracks data
%                 xx = find(tracks(:,7)==jj);
% 
%                 if numel(xx)>plotTracksMinSteps % include only tracks with at least minSteps
%                     
%                     max_tracktime = (length(xx)-1)*timeres;
%                     tracktimes = 0:timeres:max_tracktime;
%                     
%                     form trackname, get start pos X, Y and frame
%                     Xpos = round(tracks(xx(1),1));
%                     Ypos = round(tracks(xx(1),2));
%                     Sframe = tracks(xx(1),6);
%                     trackname = [savename '-X' num2str(Xpos) 'Y' num2str(Ypos) '-frame' num2str(Sframe) '.fig'];
                    
%                     figure
%                     subplot(2,1,1,'Fontsize',20)
%                     plot(tracktimes, tracks(xx,4),'g-','LineWidth',2); hold on
%                     plot(tracktimes, tracks(xx,5),'r-','LineWidth',2); title('Raw intensities','Fontsize',20);xlabel('Time(s)','Fontsize',20), ylabel('Counts','Fontsize',20); title(trackname);
%                     subplot(2,1,2,'Fontsize',20)
%                     plot(tracktimes, tracks(xx,3),'b-','LineWidth',2); title('E_{raw}','Fontsize',20); xlabel('Time (s)','Fontsize',20); ylim([0 1])
% 
%                     %save as fig
%                     saveas(gcf, trackname, 'fig')
                    
% 
%                     if strcmp(plotTracksColour, 'time')
% 
%                     % choose the colour by the track id
%                     colorindex = ceil( length(cmap) * tracks(xx(1),4) / nMolecules );
% 
%                     % plot the track in the chosen colour
%                     plot(tracks(xx,1),tracks(xx,2),...
%                         '.-','Color',cmap(colorindex,:), 'MarkerSize',3.0)
% 
%                     elseif strcmp(plotTracksColour,'random') 
%                     plot(tracks(xx,1),tracks(xx,2),...
%                         '.-','Color',rand(1,3), 'MarkerSize',3.0)            
% 
%                     else
%                      plot(tracks(xx,1),tracks(xx,2),...
%                         '.-','Color',plotTracksColour, 'MarkerSize',3.0)
% 
%                  end
% 
%              end



   
else
    %dig out individual data - ALEX case, not fully tested
    
    %AS: If acquisition is stopped manually, it might happen, that not
    %equal # of frames are recorded for DD, DA ans AA -> reduce data to
    %last complete set
    rightSize= (fix(length(data(:,1))/3))*3;
    data=data(1:rightSize,:);
    
    DDxpos = data(1:3:end,2);
    DDypos = data(1:3:end,3);

    DAxpos = data(2:3:end,2);
    DAypos = data(2:3:end,3);

    AAxpos = data(3:3:end,2);
    AAypos = data(3:3:end,3);

    DD = data(1:3:end,4);
    DA = data(2:3:end,4);
    AA = data(3:3:end,4);

    E = data(3:3:end,11);
    S = data(3:3:end,12);
    
    LINK = data(3:3:end,13); %AS for combined data


    S_DD_X = data(1:3:end,7);
    S_DD_Y = data(1:3:end,8);

    S_DA_X = data(2:3:end,7);
    S_DA_Y = data(2:3:end,8);

    S_AA_X = data(3:3:end,7);
    S_AA_Y = data(3:3:end,8);
    
    frame= data(3:3:end,1); %added by AS
    
    %figure stochimetry
    figure
    subplot(3,2,1)
    hist(S_DD_X,0:0.05:5); xlim([0 5]); title('S-DD-X'); hold on
    subplot(3,2,2)
    hist(S_DD_Y,0:0.05:5); xlim([0 5]); title('S-DD-Y'); hold on
    subplot(3,2,3)
    hist(S_DA_X,0:0.05:5); xlim([0 5]); title('S-DA-X'); hold on
    subplot(3,2,4)
    hist(S_DA_Y,0:0.05:5); xlim([0 5]); title('S-DA-Y'); hold on
    subplot(3,2,5)
    hist(S_AA_X,0:0.05:5); xlim([0 5]); title('S-AA-X'); hold on
    subplot(3,2,6)
    hist(S_AA_Y,0:0.05:5); xlim([0 5]); title('S-AA-Y'); hold on
    
    
    indrmv_S_DD = find(S_DD_X < sigmax_min | S_DD_X > sigmax_max | S_DD_Y < sigmay_min | S_DD_Y > sigmay_max);
    indrmv_S_DA = find(S_DA_X < sigmax_min | S_DA_X > sigmax_max | S_DA_Y < sigmay_min | S_DA_Y > sigmay_max);
    indrmv_S_AA = find(S_AA_X < sigmax_min | S_AA_X > sigmax_max | S_AA_Y < sigmay_min | S_AA_Y > sigmay_max);
    indrmv_S = [indrmv_S_DD; indrmv_S_DA; indrmv_S_AA]; 
    %fprintf(' %d ', max(indrmv_S))
    
    
    AAxpos(indrmv_S) = [];
    AAypos(indrmv_S) = [];
    
    DD(indrmv_S) = [];
    DA(indrmv_S) = [];
    AA(indrmv_S) = [];
    
    E(indrmv_S) = [];
    S(indrmv_S) = [];
    LINK(indrmv_S) = [];
    
    S_DD_X(indrmv_S) = [];
    S_DD_Y(indrmv_S) = [];
    
    S_DA_X(indrmv_S) = [];
    S_DA_Y(indrmv_S) = [];
    
    S_AA_X(indrmv_S) = [];
    S_AA_Y(indrmv_S) = [];
    
    frame (indrmv_S) = []; %added by AS
    
    %plot positions
%     figure;
%     plot(AAxpos, AAypos, 'r.');
    
    filter = find(S>0.3 & S<0.5);
    
    figure; hist(E(filter))
    figure; hist(E)
    figure; plot(E,S,'.'); xlim([0 1]); ylim([0 1]);
    
    figure; plot(E(filter),S(filter),'.'); xlim([0 1]); ylim([0 1]);
    figure; hist(E,0:0.02:1); title('After filtering');
    
    
       %AS: In case data from multiple movies should be combined: save
    %filename of the out file, ROI number, E and S in .comb file for later
    %processing (this part was quickly added, evtually revise later)
    fid = fopen('combinedCells.comb','a');

    for c=1:length(E)
        fprintf(fid,'%s\t',FRETSTORMFilename); %Filename
        fprintf(fid,'%d\t',LINK(c)); %ROI number
        fprintf(fid,'%6.3f\t',E(c)); %E
        fprintf(fid,'%6.3f\n',S(c)); %S
    end
    
    close(fid);
    
%     % prepare data
%    % bring data into the right format for single particle tracking
%     pos(:,1) = AAxpos;
%     pos(:,2) = AAypos;
%     pos(:,3) = E;
%     pos(:,4) = S;
%     size(frame)
%     size(AAxpos)
%     class(frame(1))
%     pos(:,5) = frame/2; %frame (accounting for alternation), corrected by AS
%     
%     %now sort rows based on frame number (needed for track.m), added by AS
%     pos = sortrows(pos,5);
    

    
  %AS: FRET has its own tracking, thats why commented out
    
%     %% tracking
%     
%     trackWindow = 1.0; % tracking window in pixels
%     
%     param.mem = 1;
%     param.dim = 2;
%     param.good = 1;
%     param.quiet = 0;
%     
%     % run tracking function
%     tracks = track(pos, trackWindow, param);
%     
%     %Comment by AS: in CW mode first and last frame not included, should
%     %this be applied to ALEX mode as well?
%     
%     % number of molecules in field of view
%     NumTrajectories = max(tracks(:,6));
%     
%    %AS: here is a bug in generatingtirfIm (commented out, as not needed yet) 
%     
% %     %specific for analysed file - SHOW IMAGE AND THEN..
% %     fileName = 'Expt1_hiFRETnorm_EPnorecov_5msalt_ALEX_G1.3mW_R0.65mW_nTIRF_192_260_40.fits';
% %     alternationPeriod = 2;
% %     imageLim = [257 512 1 69; 1 256 1 69];
% %     firstGreenFrame = 1;
% %     
% %     tirfIm = TirfImage(fileName,firstGreenFrame, imageLim, alternationPeriod);
% %     
% %     green_stack = getGreenStack(tirfIm);
% %     red_stack = getRedStack(tirfIm);
% %     
% %     figure;
% %     
% %     imageI = getFrame(red_stack,2);
% %     %specific for analysed file
% %     subim=imageI(1:69,1:256);
% %     imavg_R=cast(subim,'double');
% %     imshow(imageI,[min(imavg_R(:)), max(imavg_R(:))]);
% %     set(gca,'YDir','reverse')
% %     
% %     %% 
% %     figure; 
% %     
% %     %imshow(STD_17_AA) 
% %     imshow(imageI)
% %     set(gca,'YDir','reverse')
% %     
% %     hold all
% %     
% % %       ..THEN PLOT LOCALISATIONS COLOUR CODED (2 COLOURS) ACCORDING TO S and E ON
% % %       TOP (Double labelled and high FRET or not)
%     
%     for jj = 15000:length(AAxpos)-1
%         
%         ii = length(AAxpos) - jj;
%         
%         if S(ii)>0.35 && S(ii)<0.45
%             
%             if E(ii)>=0.5 && E(ii)<0.8
%                 
%                 colorcode = [1 0 0];
%                 
%                 plot(AAxpos(ii),AAypos(ii),'.','Color',colorcode,'MarkerSize',5.0)
%                 
%                 %         elseif E(ii)>=0.5 && E(ii)<0.6 % %
%                 colorcode = [0.5 0.5 0]; % %
%                 plot(AAxpos(ii),AAypos(ii),'.','Color',colorcode,'MarkerSize',5.0) % %         elseif E(ii)>=0.7 && E(ii)<0.8
%                 % %             colorcode = [0.8 0.8 0]; % %
%                 plot(AAxpos(ii),AAypos(ii),'.','Color',colorcode,'MarkerSize',5.0)
%                 
%                 
%             else
%                 
%                 colorcode = [1 1 0];
%                 
%                 plot(AAxpos(ii),AAypos(ii),'.','Color',colorcode,'MarkerSize',5.0)
%                 
%             end
%             
%         end
%         
%     end
% %     
% %     figure;
% %     
% %     imshow(STD_17_AA)
% %     imshow(imageI)
% %     set(gca,'YDir','reverse')
% %     
% %     hold all
% %     
% % %      HERE UPHOFF IS PLOTTING LONG AND SHORT TRACKS, IE LOCALISED ONES FOR
% % %      HIS POLYMERASE DATA
% %     
% %     for ii = 1:NumTrajectories
% %         
% %         xx = find(tracks(:,6)==ii);
% %         % don't include first and last data point in each track for averaging
% %         % because they have messed up E and S values
% %         avPos = mean(tracks(xx(2:end-1),1:4),1);
% %         
% %         if length(xx)>10 && avPos(4)>0.3 && avPos(4)<0.5
% %             
% %             if avPos(3)>=0.5 && avPos(3)<0.8
% %                 
% %                 colorcode = [1 1 0];
% %                 
% %                 plot(avPos(1),avPos(2),'.','Color',colorcode, 'MarkerSize',5.0)
% %                 
% %                 
% %             else
% %                 
% %                 colorcode = [1 0 0];
% %                 
% %                 plot(avPos(1),avPos(2),'.','Color',colorcode, 'MarkerSize',5.0)
% %                 
% %             end
% %             
% %         end
% %         
% %     end
% %     
% %     hold off;
%     
%     %%
%     
%     

    
end



% %% plot localizations each in a different colour
% 
% ColorCoding=rand(NumTrajectories,3);
% 
% figure;
% hold all;
% 
% for ii=1:length(tracks)
%     
%     %plot(tracks(ii,1), tracks(ii,2), '.', 'Color', ColorCoding(tracks(ii,4),:), 'MarkerSize', 5.0)
%     plot(tracks(ii,1),tracks(ii,2),'b.','MarkerSize',5.0)
%     % or alternatively all in the same color
%     % plot(tracks(ii,1), tracks(ii,2), 'b.', 'MarkerSize', 5.0)
%     
% end
% 
% hold off;
% 
% %% plot E map colour coded
% 
% figure;
% 
% %imshow(STD_17_AA)
% imshow(imageI)
% set(gca,'YDir','reverse')
% %freezeColors;
% 
% cmap = colormap(jet); % load colour map
% 
% E = zeros(length(tracks),1);
% kk = 1;
% hold all
% 
% for ii = 1:NumTrajectories
%     
%     xx = find(tracks(:,6)==ii);
%     
%     if length(xx) > 3
%     
%     % don't include first and last data point in each track for averaging
%     % because they have messed up E and S values
%     avPos = mean(tracks(xx(2:end-1),1:4),1);
%     
%     if  avPos(4)>0.3 && avPos(4)<0.5 && avPos(3)>0 && avPos(1)<235 && avPos(1)>20 && ...
%             avPos(2)>15 && avPos(2)<285
%         
%         
%         E(kk) = avPos(3);
%         
%         % choose the colour according to E value
%         colorindex(kk) = ceil( length(cmap(:,1)) * avPos(3) );
%         
%         % plot the track in the chosen colour
%         plot(avPos(1),avPos(2),'.','Color',cmap(colorindex(kk),:), 'MarkerSize',15.0)
%         
%         kk = kk + 1;
%         
%     end
%     
%     end
%     
% end
% 
% % delete unused rows
% E(kk:end) = [];
% 
% hold off;
% 
% figure; hist(E)



