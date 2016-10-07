function [diffusionFraction, D1, D2] = twoSpeciesMSDThreshold(tracks, params)
% MSD diffusion analysis

pixel = params.pixel;
dT = params.dT;
nMolecules = max(tracks(:,4));
sigmaNoise = params.sigmaNoise;
threshold = params.DhistThresh;
minSteps = params.DhistMinSteps;
rangeD = params.rangeD;
holdFigure = params.holdFigureCheckbox;
if holdFigure
    figureHandle = params.figureHandle;    
end

tracks(:,1:2) = tracks(:,1:2) * pixel;

D = zeros(nMolecules,1);
kk = 1;

if holdFigure
    figure(figureHandle);
else
    figure;
end
hold all
axis image;


for ii = 1:nMolecules
    
    xx = find(tracks(:,4)==ii);
    
    if numel(xx)>minSteps
        
        MSD = 0; % reset MSD
        
        % sum all squared displacement
        for jj = 1:numel(xx)-1
            
            MSD = MSD + ((tracks(xx(jj+1),1) - tracks(xx(jj),1))^2 +...
                (tracks(xx(jj+1),2) - tracks(xx(jj),2))^2);
            
        end
        
        MSD = MSD/jj; % mean square displacement
        
        % calculate D from MSD and correct for localization noise
        D(kk) = MSD/(4*dT) - sigmaNoise^2*pixel^2/dT;
        
        if D(kk) < threshold
            plot(tracks(xx,1)/pixel,tracks(xx,2)/pixel,...
                'r.-', 'MarkerSize',3.0)
        else
            plot(tracks(xx,1)/pixel,tracks(xx,2)/pixel,...
                'b.-', 'MarkerSize',3.0)
        end
        
        % reset MSD and update index
        kk = kk + 1;
        
    end
    
end

D(kk:end) = []; % delete unused rows


Dc = histc(D,rangeD'); %normalized histogram count

figure;
hist(D,rangeD');
hold all
stem(threshold,max(Dc),'marker','none');
xlabel('diffusion coefficient [um^2/s]');
ylabel('histogram count');

D1Indexes = find(D<threshold);
D2Indexes = find(D>threshold);

D1 = mean(D(D1Indexes));
D2 = mean(D(D2Indexes));

diffusionFraction = numel(D1Indexes)/numel(D);

end