function filterAverParam = findFilterParameters(imagename,tif,findFilterParamsOption)
%small function for testing and finding good parameters for filtering
%calling bpassModifiedAP and findLocalMax
    if (tif == 0)
%load image info
ImageInfo = fits_read_header(imagename);
maxframes = double(ImageInfo.NAXIS3);
ImageXWidth = ImageInfo.NAXIS1;
ImageYLength = ImageInfo.NAXIS2;
    else
        info = imfinfo(imagename);
        %get length of image stack- nFrame
        maxframes = numel(info);
    end    

%number of frames to analyze to find good filter parameters
ImageAna = 5;
%Initialization of filterParam
filterParam = zeros(ImageAna,3);

if(findFilterParamsOption == 1)
%CHOOSE FRAMES RANDOMLY
    
            %get random frames from images
            currentFrame = round(maxframes*rand(1,ImageAna));
              
elseif(findFilterParamsOption == 2)
%TAKE BRIGHTEST FRAMES

            %calculate Sum of intensities for each 50th image
            testSum = zeros(maxframes,1);
            for jj=1:50:maxframes
                    if (tif == 0)
                    testIm = double(flipud(getFrame(imagename,jj)));%for consistency with daostorm
                    else
                    testIm = double(flipud(imread(imagename,jj)));    %load tif-files
                    end  
                    testSum(jj,1) = sum(sum(testIm));
            end
            %sort intensities and store the ImageAna brightest frame indexes in currentFrame
            [testSumSort, SortIndex] = sort(testSum,'descend');
            currentFrame = SortIndex(1:ImageAna)';
            
end
            
for ii=1:ImageAna
%get (ii) random chosen frame from imagestack
    if (tif == 0)
        currentIm = double(flipud(getFrame(imagename,currentFrame(1,ii))));%for consistency with daostorm
    else
        currentIm = double(flipud(imread(imagename,currentFrame(1,ii),'Info',info)));    %load tif-files
    end  
                
quality = 0;

% Make a large figure.
figure('position',[500 300 700 700],'name','Current frame and localized molecules'); 

%  Make subplot to hold plot.
subplot(1,2,1,'position',[0.05 0.1 0.4 0.65]), subplot(1,2,2,'position',[0.55 0.1 0.4 0.65]); drawnow;
%wait til figure is open

% threshold input
uicontrol('Style', 'text', 'String', 'threshold','Position', [100 650 80 20]);
handleThresh = uicontrol('Style', 'edit', 'String', '30','Position', [180 650 30 20]);

% lnoise input
uicontrol('Style', 'text', 'String', 'lnoise','Position', [100 620 80 20]);
handleLnoise = uicontrol('Style', 'edit', 'String', '1','Position', [180 620 30 20]);

% lobject input
uicontrol('Style', 'text', 'String', 'lobject','Position', [100 590 80 20]);
handleLobject = uicontrol('Style', 'edit', 'String', '7','Position', [180 590 30 20]);

% A checkbox if localization quality is ok
uicontrol('Style', 'text', 'String', 'Tick, if localization looks ok:','Position', [300 615 150 15]);
handleRadiobutton = uicontrol('Style','radiobutton','Value', 0 ,...
    'Position',[300 595 15 15]);

% A button to localize the molecules
uicontrol('Style', 'pushbutton', 'String', 'Localize',...
    'Position', [300 640 80 30],'Callback', 'uiresume');

while (quality == 0)

%get filter values and threshold from the edit field
thresh = str2num(get(handleThresh, 'String'));
lnoise = str2num(get(handleLnoise, 'String'));
lobject = str2num(get(handleLobject, 'String'));

%bandpass filter on current frame
filteredIm = bpassModifiedAP(currentIm,lnoise,lobject);
%localize molecules
posFil = findLocalMax(filteredIm,thresh);

%show localization in figure(1)
subplot(1,2,1,'position',[0.05 0.1 0.4 0.65]), imagesc(currentIm), axis equal;
subplot(1,2,2,'position',[0.55 0.1 0.4 0.65]), imagesc(filteredIm), colormap(gray), axis equal, hold on, ...
plot(posFil(:,1),posFil(:,2),'rO'), hold off;

%wait until localize is pressed
uiwait;
%is set to '1', if localization was 'successful'
quality = get(handleRadiobutton,'Value');
end

filterParam(ii,1) = thresh;
filterParam(ii,2) = lnoise;
filterParam(ii,3) = lobject;

close(figure(1));
end
close(figure(1));
%averaged filter parameters
filterAverParam(1,1) = sum(filterParam(:,1))/ImageAna;
filterAverParam(1,2) = sum(filterParam(:,2))/ImageAna;
filterAverParam(1,3) = sum(filterParam(:,3))/ImageAna;


%------------------------------------------------------------------------------
function Frame = getFrame(file, n)

ImageInfo = fits_read_header(file);
NAXIS1 = ImageInfo.NAXIS1;
NAXIS2 = ImageInfo.NAXIS2;
%check whether there is more than one frame
if isfield(ImageInfo,'NAXIS3')
        nFrame= ImageInfo.NAXIS3;
else
        nFrame= 1;
end

lim = [1 NAXIS1 1 NAXIS2];
if n > nFrame || n < 1
	error('Attempted to access out of bounds frame');
end

if nFrame == 1 %ie if there is only one frame
	startVector = double([lim(1), lim(3)]);
	endVector =   double([lim(2), lim(4)]);
else
	startVector = double([lim(1), lim(3), n]);
	endVector =   double([lim(2), lim(4), n]);
end

% the axis definitions in matlab for fits files are 
% 90 degrees to the normal definition
% so rotate it
Frame = rot90(fits_read_image_subset(file, startVector, endVector));

