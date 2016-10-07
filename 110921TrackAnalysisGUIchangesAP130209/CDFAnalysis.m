function results = CDFAnalysis(tracks, params, initGuess,CDFfittingOption)
% CDFAnalysis
% CDFfittingOption = 1: no fit
% CDFfittingOption = 2: single species fit D
% CDFfittingOption = 3: two species fit a. D1 and D2 fixed
% CDFfittingOption = 4: two species fit a, D1, D2

sigmaNoise = params.sigmaNoise;
dT = params.dT;
pixel = params.pixel;
nMolecules = max(tracks(:,4));
bacArea = params.bacArea;
holdFigure = params.holdFigureCheckbox;
curveColor = params.curveColor;

if holdFigure
    figureHandle = params.figureHandle;    
end

% calculate absolute one-step displacements
distances = zeros(nMolecules,1);
kk = 1;
for ii = 1:nMolecules
    
    xx = find(tracks(:,4)==ii);
    
    if numel(xx)>1
        
        for jj = 1:numel(xx)-1
            
            distances(kk) = sqrt((tracks(xx(jj+1),1) - tracks(xx(jj),1))^2 +...
                (tracks(xx(jj+1),2) - tracks(xx(jj),2))^2);
            
            kk = kk+1;
            
        end
        
    end
    
end
% delete unused rows
distances(kk:end) = [];

% range of cumulative correlation function upper threshold
lCum = min(distances) : 0.1*(min(distances)) : 2*max(distances);

cFraction = zeros(1,length(lCum));

for ii = 1 : length(lCum)
    
    % fraction of measured distances within lCum to total number of steps
    cFraction(ii) = sum(distances<lCum(ii))/numel(distances); % cumulative correlation function value at lCum
    
end

% fit cFraction vs lCum with model

curvefitoptions = optimset( 'lsqcurvefit');
curvefitoptions = optimset( curvefitoptions, 'Display', 'off', 'MaxFunEvals', 2000, 'MaxIter', 2000);


if CDFfittingOption == 1 % no fit
    results = [];
    
if holdFigure
    figure(figureHandle);
    hold all;
else
    figure;
end
    plot(lCum, cFraction,'Color',curveColor)
    xlabel('l [um]');
    ylabel('CDF(l)');
    
    return;
end

if CDFfittingOption == 2 % single species fit
    
    lb = 0;
    ub = inf;
    
    D1free = initGuess(1)/pixel^2;
    
    c = lsqcurvefit(@(parameters,lCum) ...
        singleConfinedAverageFrameDiffusionSpeciesModel(...
        parameters,lCum,dT,sigmaNoise,bacArea),...
        D1free, lCum, cFraction, lb, ub, curvefitoptions);
    
    results(1) = c * pixel^2;
    
    cFit = singleConfinedAverageFrameDiffusionSpeciesModel(...
        c,lCum,dT,sigmaNoise,bacArea);
    
elseif CDFfittingOption == 3 % two species fit a, D1 and D2 fixed
    
    lb = 0;
    ub = 1;
    
    D1fixed = initGuess(2)/pixel^2;
    D2fixed = initGuess(3)/pixel^2;
    
    c = lsqcurvefit(@(parameters,lCum) ...
        twoConfinedAverageFrameDiffusionSpeciesBothFixedModel(...
        parameters,lCum,dT,sigmaNoise,bacArea,D1fixed,D2fixed),...
        initGuess(1), lCum, cFraction, lb, ub, curvefitoptions);
    
    results(1) = c;
    results(2) = D1fixed*pixel^2;
    results(3) = D2fixed*pixel^2;
    
    cFit = twoConfinedAverageFrameDiffusionSpeciesBothFixedModel(...
        c,lCum,dT,sigmaNoise,bacArea,D1fixed,D2fixed);
    
    
elseif CDFfittingOption == 4 % two species fit D1, D2, a
    
    lb = [0 0 0];
    ub = [1 inf inf];
    
    D1free = initGuess(2)/pixel^2;
    D2free = initGuess(3)/pixel^2;
    
    c = lsqcurvefit(@(parameters,lCum) ...
        twoConfinedAverageFrameDiffusionSpeciesModel(...
        parameters,lCum,dT,sigmaNoise,bacArea),...
        [initGuess(1) D1free D2free], lCum, cFraction, lb, ub, curvefitoptions);
    
    results(1) = c(1);
    results(2) = c(2) * pixel^2;
    results(3) = c(3) * pixel^2;
    
    cFit = twoConfinedAverageFrameDiffusionSpeciesModel(...
        c,lCum,dT,sigmaNoise,bacArea);
      
end

if holdFigure
    figure(figureHandle);
    hold all;
else
    figure;
end
plot(lCum, cFraction,'Color',curveColor)
hold all;
plot(lCum, cFit)
legend('data','fit')
xlabel('l [um]');
ylabel('CDF(l)');
hold off;

fid = fopen(['CDF_thresh15_win5_nMol' num2str(nMolecules) '_FitOption' num2str(CDFfittingOption) '.txt'],'w');
for ina1=1:length(lCum)
  fprintf(fid,'%6.3f\t',lCum(ina1)); %lCum
  fprintf(fid,'%6.3f\t',cFraction(ina1)); %cFraction
  fprintf(fid,'%6.3f\n',cFit(ina1)); %cFit 
end
fclose(fid);

end