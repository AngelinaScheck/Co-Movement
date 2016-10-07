function D = histD(tracks,params)
% function D = histD(tracks, params) calculates the diffusion coefficient
% for each track in "tracks" and plots a histogram of D.
%
% Stephan Uphoff. 09.09.11
%
% save D histogram in text file- to combine several data later on
% AP 30.01.12
%
% specify textfile of Dcoeff: number of localized molecules, number of tracked molecules, parameter settings
% AP 03.06.12
%
% plot and save Dapp from single step distance^2 / (4dT)
% introduce MSDall and Dall (single step) and MSD and D (mean values over track)
% AP 22.01.13


totalN = 0;% AS
pixel = params.pixel; % length per pixel
dT = params.dT; % time per frame
sigmaNoise = params.sigmaNoise; % localization noise
DhistMinSteps = params.DhistMinSteps; % minimum number of steps for a track to be analyzed
rangeD = params.rangeD; % D range for the histogram
DallSelect = params.checkboxDallSelect;
nMolecules = max(tracks(:,4));
% for ll= 1:nMolecules
% longestTrack(ll) = length(find(tracks(:,4)==ll));
% end
% LLongestTrack = max(longestTrack);

if DallSelect ==0
    kk = 1;
    MSD = zeros(nMolecules,1);
else
    kk= 0;
    %MSDall = zeros(nMolecules*LLongestTrack,1);
    %Dall = zeros(nMolecules*LLongestTrack,2);
end

for ii = 1:nMolecules
    
    xx = find(tracks(:,4)==ii);
    %Option1: Average D over all MSD of track 
    if (DallSelect ==0)
        if numel(xx) > DhistMinSteps %only tracks longer that DhistMinSteps (same param for Dhist plot)
            totalN = totalN + 1;% AS

            % sum all squared displacement in the track
            for jj = 1:numel(xx)-1

                MSD(kk) = MSD(kk) + ((tracks(xx(jj+1),1) - tracks(xx(jj),1))^2 +...
                    (tracks(xx(jj+1),2) - tracks(xx(jj),2))^2);

            end

            MSD(kk) = MSD(kk)/jj; % mean square displacement
            kk = kk + 1;

        end
    end
    
    %Option2: D for all single SDs (each step)
    if (DallSelect > 0)
        if numel(xx) > DhistMinSteps %only tracks longer that DhistMinSteps (same param for Dhist plot)

            % calculate step wise all MSDs and apparent D coefficients
            for jj = 1:numel(xx)-1

                MSDall(kk+jj,1) = ((tracks(xx(jj+1),1) - tracks(xx(jj),1))^2 +...
                    (tracks(xx(jj+1),2) - tracks(xx(jj),2))^2);
                MSDall(kk+jj,1) = MSDall(kk+jj,1) * pixel^2;
                Dall(kk+jj,1) = MSDall(kk+jj,1)/(4*dT)-sigmaNoise^2*pixel^2/dT; %diffusion coeff with localization noise correction
                %Dall(kk+jj,2) = MSDall(jj,1)/(4*dT); %diffusion coeff without localization noise correction
                                
            end
            %trackLength = numel(xx);
            %trackID = tracks(xx(1),4);
            kk = kk + jj;
        end
    end
end

% calculate D from MSD and correct for localization noise
%Option1
if (DallSelect ==0)
    MSD(kk:end) = []; % delete unused rows
    MSD = MSD * pixel^2; % convert from pixel to length units
    D = MSD/(4*dT) - sigmaNoise^2*pixel^2/dT; %diffusion coeff with localization noise correction
    %D = MSD/(4*dT); %diffusion coeff without localization noise correction

    figure;
    hist(D,rangeD);
    xlim([min(rangeD), max(rangeD)]);
    xlabel('diffusion coefficient [um^2/s]');
    ylabel('histogram count'); drawnow;

    fidDmean = fopen(['Dcoeff_mean' '_thresh' num2str(params.localizationThresh) '_win' num2str(params.trackParams.maxDisp) '_mem' num2str(params.trackParams.mem) '_nMol' num2str(nMolecules) '_MinSteps' num2str(DhistMinSteps) '_nD' num2str(length(D(:,1))) '.txt'],'w');
    for ina1=1:length(D(:,1))
      fprintf(fidDmean,'%6.3f\n',D(ina1,1)); %D coefficients  
      %fprintf(fid,'%6.3f\n',Dxout(1,ina1)); %D range 
    end
    fclose(fidDmean);
    %by AS: display total number of molecules in histogramm
    dim = [.2 .5 .3 .3];
    annotation('textbox',dim, 'String', strcat(num2str(totalN), 'particles in total'),'FitBoxToText','on');
end
%keyboard;
%Option2
if (DallSelect > 0)
    Dall(Dall == 0) = NaN;
    D = Dall;
    figure;
    hist(Dall(:,1),rangeD);
    xlim([min(rangeD), max(rangeD)]);
    xlabel('diffusion coefficient [um^2/s]');
    ylabel('histogram count'); drawnow;

    fidDall = fopen(['Dcoeff_all' '_thresh' num2str(params.localizationThresh) '_win' num2str(params.trackParams.maxDisp) '_mem' num2str(params.trackParams.mem) '_nMol' num2str(nMolecules) '_MinSteps' num2str(DhistMinSteps) '_nDall' num2str(length(Dall(:,1))) '.txt'],'w');
    for ina1=1:length(Dall(:,1))
      fprintf(fidDall,'%6.3f\n',Dall(ina1,1)); %D coefficients  
      %fprintf(fid,'%6.3f\n',Dxout(1,ina1)); %D range 
    end
    fclose(fidDall);
end

end


