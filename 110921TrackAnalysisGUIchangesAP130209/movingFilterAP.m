function movingD = movingFilterAP( D, window )
% AP, 22.12.12
% D filtering, such that mean D value within window gets assigned
% other fancy filtering could be tried/research needed
movingD = zeros(size(D));   %start and end stay "0"
Dinterval = zeros(1,2*window+1);    %Def. of window is no. of time points lower/larger than actual time point
                                    %total D values averaged: 2*window+1
for ii=window+1:numel(D)-window
    Dinterval = D(ii-window:ii+window);
    movingD(ii) = mean(Dinterval);
end
movingD(movingD == 0) = [];    %start and end become NaN
end

