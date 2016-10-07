function ROIData = selectROIs(data, alttype)

% load gaussstorm.pos.out into workspace
% pos(:,1) = data.data(:,2);
% pos(:,2) = data.data(:,3);

pos(:,1) = data(:,2);
pos(:,2) = data(:,3);

% figure;
figure(1);
hold on;
plot(pos(:,1),pos(:,2),'.')
hold all;

ROIData = [];

reply = 'Y';
ii = 1;

%added by AS: ALEX mode gives 12 columns in out file (stochiometry), CW
%only 11 
if alttype==1
    lastColumnN=12;
    fprintf('CW')
else
    lastColumnN=13;
    fprintf('ALEX')
end

while reply == 'Y'
    
    p = impoly;
    vertices = getPosition(p);
    
    inIndexes = find( inpolygon(pos(:,1),pos(:,2),vertices(:,1),vertices(:,2)) );
%     inData = data.data(inIndexes,:);
    inData = data(inIndexes,:);
    
    plot(inData(:,2),inData(:,3),'o')
    hold all;
    
    inData(:,lastColumnN) = ii * ones(numel(inIndexes),1); % cell id
    
    ROIData = [ROIData; inData];
    
    reply = input('Do you want more? Y/N: ', 's');
    if isempty(reply)
    reply = 'Y';
    end
    
    ii = ii + 1;
    
end

hold off;

end
