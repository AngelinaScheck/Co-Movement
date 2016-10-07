function D = singleConfinedAverageFramePICS(macrotracks, params, initGuess)
% single species diffusion PICS

sigmaNoise = params.sigmaNoise;
dT = params.dT;
pixel = params.pixel;
nMolecules = params.nMolecules;
bacArea = params.bacArea;

% calculate absolute one-step displacements
distances = zeros(nMolecules,1);
kk = 1;
for ii = 1:nMolecules
    
    xx = find(macrotracks(:,4)==ii);
    
    if numel(xx)>1
        
        for jj = 1:numel(xx)-1
            
            distances(kk) = sqrt((macrotracks(xx(jj+1),1) - macrotracks(xx(jj),1))^2 +...
                (macrotracks(xx(jj+1),2) - macrotracks(xx(jj),2))^2);

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
lb = 0;
ub = inf;

curvefitoptions = optimset( 'lsqcurvefit');
curvefitoptions = optimset( curvefitoptions, 'Display', 'off', 'MaxFunEvals', 2000, 'MaxIter', 2000);

c = lsqcurvefit(@(parameters,lCum) ...
    singleConfinedAverageFrameDiffusionSpeciesModel(...
    parameters,lCum,dT,sigmaNoise,bacArea),...
    initGuess, lCum, cFraction, lb, ub, curvefitoptions);

D = c * pixel^2;

cFit = singleConfinedAverageFrameDiffusionSpeciesModel(...
    c,lCum,dT,sigmaNoise,bacArea);

figure;
plot(lCum, cFraction,'r')
hold all;
plot(lCum, cFit)
legend('data','fit')
xlabel('l [um]');
ylabel('CDF(l)');
hold off;

end