function [D locD locFrameD locFrameConfineD] = averageFrameMSD(tracks, params)

pixel = params.pixel;
dT = params.dT;
nMolecules = max(tracks(:,4));
sigmaNoise = params.sigmaNoise;
bacArea = params.bacArea;
maxStep = params.maxStep; % Number of time steps for MSD curve
holdFigure = params.holdFigureCheckbox;
curveColor = params.curveColor;

if holdFigure
    figureHandle = params.figureHandle;    
end

mMSD = zeros(1,maxStep+1);
errorMSD = zeros(1,maxStep+1);
MSD = zeros(nMolecules,1);

for ll = 1:maxStep
    
    for ii = 1:nMolecules
        
        xx = find(tracks(:,4)==ii);
        
        if numel(xx)>ll
            
            SquaredDisplacement = 0;
            
            for jj = 1:numel(xx)-ll
                
                SquaredDisplacement = SquaredDisplacement +...
                (tracks(xx(jj+ll),1) - tracks(xx(jj),1))^2 +...
                    (tracks(xx(jj+ll),2) - tracks(xx(jj),2))^2;
                
            end
            
            MSD(ii) = SquaredDisplacement/jj; % mean squared displacement
            
        else
            
            MSD(ii) = -1;
            
        end
        
    end
    
    MSD = MSD(MSD>=0);
    
    mMSD(ll+1) = mean(MSD);
    
    errorMSD(ll+1) = std(MSD);
    
    MSD = [];
    
end

mMSD = mMSD * pixel^2;
errorMSD = errorMSD * pixel^2;

% simple estimate of D, no corrections
D = mMSD(2)/(4*dT);

% accounting for localization noise
locD =  mMSD(2)/(4*dT) - sigmaNoise^2*pixel^2/dT;

% diffusion coefficient accounting for average frame effect and
% localization noise
locFrameD = 3/8 * mMSD(2)/dT - sigmaNoise^2*pixel^2/dT;

% average diffusion over frame dT offset in MSD curve
mMSD(1) = 4*sigmaNoise^2*pixel^2 - 4/3*locFrameD*dT;

% accounting for confinement, localization accuracy, average frame effect
locFrameConfineD = -bacArea*pixel^2/(12*dT) *...
    log( 1 - 3/(bacArea*pixel^2) *...
    ( mMSD(2) - 4*sigmaNoise^2*pixel^2 + 4/3*D*dT) );

time = (0:maxStep) * dT;

if holdFigure
    figure(figureHandle);
    hold all;
else
    figure;
end
plot(time,mMSD,'Color',curveColor);
%errorbar(time, mMSD, errorMSD,'Color',curveColor);
xlabel('lag time [s]')
ylabel('MSD [um^2]')

fid = fopen(['MSD_thresh15_win5_nMol' num2str(nMolecules) '.txt'],'w');
for ina1=1:length(time)
  fprintf(fid,'%6.3f\t',time(ina1)); %lag time
  fprintf(fid,'%6.3f\t',mMSD(ina1)); %mean MSD
  fprintf(fid,'%6.3f\n',errorMSD(ina1)); %std MSD
end
fclose(fid);

