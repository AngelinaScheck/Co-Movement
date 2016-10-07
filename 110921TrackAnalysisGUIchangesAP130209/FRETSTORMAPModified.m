function savename = FRETSTORM(fileName, thresh, lobject, firstGreenFrame, TFORMpath, alternationPeriod, tif, alttype, startframe, leftRim)

%FRET STORM Analysis

% based on gaussStorm

% generate tform file
% find first green frame
% start with red excitation. localize and measure AA
% use tform to find corresponding positions in DD
% measure DD based on AA localization with tform.
% measure DA based on AA localization.
% calculate E=DA/(DD+DA) and S=AA/(AA+DD+DA) for each localization
% loop over frames
% 
% reconstruct image:
% track localizations
% average AA position, E and S
% plot only localizations with intermediate S
% plot each localization colour coded based on E
%
% example use:
% FRETSTORM('movie.fits', 100, 8, 2, 'registration.tform.mat')
% 
% 15/02/12 Anne Plochowietz: include FRETSTORM functionality into
% gaussStorm GUI (click button and provide required input file)
% 
% Edited 03/06/12 RC: added CW functionality (alttype) and for ALEX, keep molecules
% only if good fit in AA, DD and DA channel (filtering other params happens
% in plotting function).  Added startframe variable.
%
% 02/02/13 Anne Plochowietz
% image limits independent from setup (read-in with movie)
% imageLim is a 2x4 matrix containing the limits of the green and red channels
%	  ie, [ greenXstart, greenXend, greenYstart, greenYend;
%		redXstart, redXend, redYstart, redYend];
% recommended values to setup FRET channels x-dimensions:
% HD: DA: 1-256, DD: 257-512
% ZS: DA: 172-342, DD: 1-171    %BLUE CHANNEL shouldn't be recorded!!!
% !!! Transform filehas to have the same dimensions !!!
%
% now implemented that initial guesses for DDPos are only drop
% when out of windowSize/2 (taken diffusion during alternation period and 
% error through transformMatrix into account)
% filter parameters are set before for-loop
% DDPos guesses are calculated from transformed AAPosFinal/DAPosFinal 
% (after fitting the PSF in the red channel)
%
% 14/02/13 Anne Plochowietz
% set left rim, such that no crap is picked up and then linked
% leftRim parameter used before image is passed to findlocalMax function 
%
% 16/02/13 Anne Plochowietz
% added linkage value in .out file to better link pos/mol in green and red
% channel (col 13 in ALEX mode, col12 in CW mode)
% overcome inconsistency with lobject (now for bpass only),
% bpasskernerldiametre (removed) and windowSize (=lobject)

if (alttype == 2) %ALEX
fprintf('\n\nALEX FRETSTORM analysis for file %s\n',fileName);
savename    = [fileName(1:end-5),'_thresh',num2str(thresh),'_lobject',num2str(lobject),'_startframe', num2str(startframe),'_ALEX_FRET.out'];
%definitions for the output file
saveFileHeader = 'FRAME\tXCENTER\tYCENTER\tBRIGHTNESS\tBG\tI0\tS_X\tS_Y\tTHETA\tECC\tE\tS\tLink\n';
else %CW
fprintf('\n\nCW FRETSTORM analysis for file %s\n',fileName);
savename    = [fileName(1:end-5),'_thresh',num2str(thresh),'_lobject',num2str(lobject),'_startframe', num2str(startframe),'_CW_FRET.out'];    
%definitions for the output file
saveFileHeader = 'FRAME\tXCENTER\tYCENTER\tBRIGHTNESS\tBG\tI0\tS_X\tS_Y\tTHETA\tECC\tE\tLink\n';
end

if (tif == 0)
    ImageInfo = fits_read_header(fileName);
%     imageLim: green channel; red channel
    imageLim = [double(round(ImageInfo.NAXIS1/2+1)) 512 1 ImageInfo.NAXIS2; 1 double(round(ImageInfo.NAXIS1/2)) 1 ImageInfo.NAXIS2];
%     Halfdome imageLim
%     imageLim = [257 512 1 ImageInfo.NAXIS2; 1 256 1 ImageInfo.NAXIS2];
%     Zugspitze imageLim
%     imageLim = [1 171 1 ImageInfo.NAXIS2; 172 342 1 ImageInfo.NAXIS2];
%     check whether there is more than one frame
        if isfield(ImageInfo,'NAXIS3')
                nFrame= ImageInfo.NAXIS3;
        else
                nFrame= 1;
        end
   
else
    info = imfinfo(fileName);
    firstImage = imread(fileName,1,'Info',info);
    [H,W,D] = size(firstImage);
    %imageLim: green channel; red channel
    imageLim = [double(round(W/2+1)) W 1 H; 1 double(round(W/2)) 1 H]; % green channel; red channel
    %Halfdome imageLim
    %imageLim = [257 512 1 H; 1 256 1 H];
    %Zugspitze imageLim
    %imageLim = [1 171 1 H; 172 342 1 H];
    
    nFrame = numel(info); %nFrame = getNumFrames(green_stack);
end    

tirfIm = TirfImage(fileName,firstGreenFrame, imageLim, alternationPeriod);

TFORM=load(TFORMpath);
tform = TFORM.TFORM;

green_stack = getGreenStack(tirfIm);
red_stack = getRedStack(tirfIm);

xsize = double(imageLim(2,2));
ysize = double(imageLim(1,4));

%row definitions, similar to gaussStorm localization routine
rowDef.frameCol = 1;
rowDef.xCol = 2;
rowDef.yCol = 3;

windowSize=lobject;%threshold for disregarding PSF fits

%open the output file for writing
fid = fopen(savename,'w');
fprintf(fid, saveFileHeader);

count = 0;  
%ALEX analysis
if (alttype == 2)
    
    % Increment startframe number by 1 if it doesn't start on green
    if (firstGreenFrame == 2)
        startframe = startframe + 1;
    end
    
    %loop over frames
    for nn = startframe:2:nFrame-1 % for each frame

      if (mod(nn,30) == 0 )
            fprintf(strcat(num2str(nFrame - nn),' frames left.\n'));
      else
        fprintf('.');
      end

      DDIm = double(getFrame( green_stack, nn));
      DAIm = double(getFrame( red_stack, nn));
      AAIm = double(getFrame( red_stack, nn+1));

      AAImFilt = bpass(AAIm,1,lobject);

      % first position guesses in AA
      AAPos=findLocalMax(AAImFilt(:,leftRim:imageLim(2,2)),thresh);  %get rid off left rim after bpass filtering

      % set fitting parameters
      initguess = [];
      AAposThresh = windowSize/2; % if the pos moves more than windowSize/2 from initialPos, drop it
      sigmaLim = [0 inf]; %Sigma can be filtered post-analysis
      DDposThresh = windowSize/2;%windowSize/100;  %takes transform error for DD into account & diffusion
      DAposThresh = windowSize/2; %takes diffusion into account 
      
      %transform green frame to fit red channel
      DDIm_tfm = imtransform(DDIm, tform, 'FillValues', 0, 'Xdata', [1 xsize], 'Ydata', [1 ysize]);
      
      % loop over positions
      for j = 1:size(AAPos,1)
        %do the fit in the red channel
        [AAphot_count, AAa, AAnormChi2, AAposFinal, AAeccentricity] = freeGaussFitEllipse( AAIm, AAPos(j,:), windowSize,AAposThresh, sigmaLim, initguess,'autoDetect');
     
        %do the corresponding DD and DA fits
        [DDphot_count, DDa, DDnormChi2, DDposFinal, DDeccentricity] = freeGaussFitEllipse( DDIm_tfm, AAposFinal, windowSize,DDposThresh, sigmaLim, initguess,'autoDetect');
        [DAphot_count, DAa, DAnormChi2, DAposFinal, DAeccentricity] = freeGaussFitEllipse( DAIm, AAposFinal, windowSize,DAposThresh, sigmaLim, initguess,'autoDetect');

        DDposFinal = tforminv(tform,DDposFinal);    %transform fitted positions in green channel
        
        E = DAphot_count / (DDphot_count + DAphot_count);
        S = (DDphot_count + DAphot_count) / (DDphot_count + DAphot_count + AAphot_count);

        if ~all(AAa == 0 & DDa ==0 & DAa ==0)%if the fit has worked
        count = count +1;
          %for DD
          %need to write following:'FRAME\tXCENTER\tYCENTER\tPHOTONS\tBG\tI0\tS_X\tS_Y\tTHETA\tECC\tE=NaN\tS=NaN\tLink\n';
          fprintf(fid,'%d\t',nn); %FRAME
          fprintf(fid,'%6.3f\t',DDposFinal(1)); %XCENTER
          fprintf(fid,'%6.3f\t',DDposFinal(2)); %YCENTER
          fprintf(fid,'%6.3f\t',DDphot_count); %PHOTONS
          fprintf(fid,'%6.3f\t',DDa(4)); %BG
          fprintf(fid,'%6.3f\t',DDa(1)); %I0
          fprintf(fid,'%6.3f\t',DDa(2)); %S_X
          fprintf(fid,'%6.3f\t',DDa(3)); %S_Y
          fprintf(fid,'%6.3f\t',DDa(7)); %THETA
          fprintf(fid,'%6.3f\t',DDeccentricity); %ECC
          fprintf(fid,'%6.3f\t',NaN); %E
          fprintf(fid,'%6.3f\t',NaN); %S
          fprintf(fid,'%6.3f\n',count); %Link
          
          %for DA
          %need to write following:'FRAME\tXCENTER\tYCENTER\tPHOTONS\tBG\tI0\tS_X\tS_Y\tTHETA\tECC\tE=NaN\tS=NaN\tLink\n';
          fprintf(fid,'%d\t',nn); %FRAME
          fprintf(fid,'%6.3f\t',DAposFinal(1)); %XCENTER
          fprintf(fid,'%6.3f\t',DAposFinal(2)); %YCENTER
          fprintf(fid,'%6.3f\t',DAphot_count); %PHOTONS
          fprintf(fid,'%6.3f\t',DAa(4)); %BG
          fprintf(fid,'%6.3f\t',DAa(1)); %I0
          fprintf(fid,'%6.3f\t',DAa(2)); %S_X
          fprintf(fid,'%6.3f\t',DAa(3)); %S_Y
          fprintf(fid,'%6.3f\t',DAa(7)); %THETA
          fprintf(fid,'%6.3f\t',DAeccentricity); %ECC
          fprintf(fid,'%6.3f\t',NaN); %E
          fprintf(fid,'%6.3f\t',NaN); %S    
          fprintf(fid,'%6.3f\n',count); %Link

          %for AA
          %need to write following:'FRAME\tXCENTER\tYCENTER\tPHOTONS\tBG\tI0\tS_X\tS_Y\tTHETA\tECC\tE\tS\tLink\n';
          fprintf(fid,'%d\t',nn+1); %FRAME
          fprintf(fid,'%6.3f\t',AAposFinal(1)); %XCENTER
          fprintf(fid,'%6.3f\t',AAposFinal(2)); %YCENTER
          fprintf(fid,'%6.3f\t',AAphot_count); %PHOTONS
          fprintf(fid,'%6.3f\t',AAa(4)); %BG
          fprintf(fid,'%6.3f\t',AAa(1)); %I0
          fprintf(fid,'%6.3f\t',AAa(2)); %S_X
          fprintf(fid,'%6.3f\t',AAa(3)); %S_Y
          fprintf(fid,'%6.3f\t',AAa(7)); %THETA
          fprintf(fid,'%6.3f\t',AAeccentricity); %ECC
          fprintf(fid,'%6.3f\t',E); %E
          fprintf(fid,'%6.3f\t',S); %S
          fprintf(fid,'%6.3f\n',count); %Link
        end

      end
    end
    
else
    %CW CASE
    %loop over frames
    for nn = startframe:1:nFrame % for each frame

      if (mod(nn,30) == 0 )
            fprintf(strcat(num2str(nFrame - nn),' frames left.\n'));
      else
        fprintf('.');
      end

      DDIm = double(getFrame( green_stack, nn));
      DAIm = double(getFrame( red_stack, nn));
      
      %transform green frame to fit red
      DDIm_tfm = imtransform(DDIm, tform, 'FillValues', 0, 'Xdata', [1 xsize], 'Ydata', [1 ysize]);

      %sum DA and DD (transformed) to find positions
      DDDAIm_sum = DAIm + DDIm_tfm;

      DDDAIm_sumFilt = bpass(DDDAIm_sum,1,lobject);

      % find position guesses from summed image
      threshCW = thresh * 2; %because we have summed two images.
      DAPos=findLocalMax(DDDAIm_sumFilt(:,leftRim:imageLim(2,2)),threshCW); %get rid off left rim after bpass filtering

      % find corresponding DD position guesses with tform from initial guesses in 
      % red channel
      %DDPos = tforminv(tform,DDDAPos);

      %set fitting parameters
      initguess = [];
      DAposThresh = windowSize/2; % if the pos moves more than windowSize/2 from initialPos, drop it
      sigmaLim = [0 inf]; % Sigma can be filtered post-analysis        
      DDposThresh = windowSize/100;  % only taking error with transform into account and force fit
      
      %transform green frame to fit red channel
      DDIm_tfm = imtransform(DDIm, tform, 'FillValues', 0, 'Xdata', [1 xsize], 'Ydata', [1 ysize]);
    
      % loop over positions
      for j = 1:size(DAPos,1)

        %do the AA fit
        [DAphot_count, DAa, DAnormChi2, DAposFinal, DAeccentricity] = freeGaussFitEllipse( DAIm, DAPos(j,:), windowSize,DAposThresh, sigmaLim, initguess,'autoDetect');
       
        %do the DD fit
        [DDphot_count, DDa, DDnormChi2, DDposFinal, DDeccentricity] = freeGaussFitEllipse( DDIm_tfm, DAposFinal, windowSize,DDposThresh, sigmaLim, initguess,'autoDetect');

        DDposFinal = tforminv(tform,DDposFinal);    %transform fitted positions in green channel
        
        E = DAphot_count / (DDphot_count + DAphot_count);
        

        if ~all(DDa ==0 & DAa ==0)%if the fit has worked
          count = count +1;  
          %for DD
          %need to write following:'FRAME\tXCENTER\tYCENTER\tPHOTONS\tBG\tI0\tS_X\tS_Y\tTHETA\tECC\tE=NaN\tLink\n';
          fprintf(fid,'%d\t',nn); %FRAME
          fprintf(fid,'%6.3f\t',DDposFinal(1)); %XCENTER
          fprintf(fid,'%6.3f\t',DDposFinal(2)); %YCENTER
          fprintf(fid,'%6.3f\t',DDphot_count); %PHOTONS
          fprintf(fid,'%6.3f\t',DDa(4)); %BG
          fprintf(fid,'%6.3f\t',DDa(1)); %I0
          fprintf(fid,'%6.3f\t',DDa(2)); %S_X
          fprintf(fid,'%6.3f\t',DDa(3)); %S_Y
          fprintf(fid,'%6.3f\t',DDa(7)); %THETA
          fprintf(fid,'%6.3f\t',DDeccentricity); %ECC
          fprintf(fid,'%6.3f\t',NaN); %E
          fprintf(fid,'%6.3f\n',count); %Link

          %for DA
          %need to write following:'FRAME\tXCENTER\tYCENTER\tPHOTONS\tBG\tI0\tS_X\tS_Y\tTHETA\tECC\tE\tLink\n';
          fprintf(fid,'%d\t',nn); %FRAME
          fprintf(fid,'%6.3f\t',DAposFinal(1)); %XCENTER
          fprintf(fid,'%6.3f\t',DAposFinal(2)); %YCENTER
          fprintf(fid,'%6.3f\t',DAphot_count); %PHOTONS
          fprintf(fid,'%6.3f\t',DAa(4)); %BG
          fprintf(fid,'%6.3f\t',DAa(1)); %I0
          fprintf(fid,'%6.3f\t',DAa(2)); %S_X
          fprintf(fid,'%6.3f\t',DAa(3)); %S_Y
          fprintf(fid,'%6.3f\t',DAa(7)); %THETA
          fprintf(fid,'%6.3f\t',DAeccentricity); %ECC
          fprintf(fid,'%6.3f\t',E); %E
          fprintf(fid,'%6.3f\n',count); %Link

        end

      end
    end
end
    

fclose(fid);
fprintf('\nDone.\n');


%----------------------------------------------------------------------------------------------
function backspace(nchar)
%function backspace(nchar)

for i = 1:nchar
  fprintf('\b');
end
%------------------------------------------------------------------------------
% function Frame = getFrame(file, n)
% 
% ImageInfo = fits_read_header(file);
% NAXIS1 = ImageInfo.NAXIS1;
% NAXIS2 = ImageInfo.NAXIS2;
% %check whether there is more than one frame
% if isfield(ImageInfo,'NAXIS3')
%         nFrame= ImageInfo.NAXIS3;
% else
%         nFrame= 1;
% end
% 
% lim = [1 NAXIS1 1 NAXIS2];
% if n > nFrame || n < 1
% 	error('Attempted to access out of bounds frame');
% end
% 
% if nFrame == 1 %ie if there is only one frame
% 	startVector = double([lim(1), lim(3)]);
% 	endVector =   double([lim(2), lim(4)]);
% else
% 	startVector = double([lim(1), lim(3), n]);
% 	endVector =   double([lim(2), lim(4), n]);
% end
% 
% % the axis definitions in matlab for fits files are 
% % 90 degrees to the normal definition
% % so rotate it
% Frame = rot90(fits_read_image_subset(file, startVector, endVector));

%----------------------------------------------------------------------------------
function res = bpass(image_array,lnoise,lobject)
% 
% NAME:
%               bpass
% PURPOSE:
%               Implements a real-space bandpass filter that suppresses 
%               pixel noise and long-wavelength image variations while 
%               retaining information of a characteristic size.
% 
% CATEGORY:
%               Image Processing
% CALLING SEQUENCE:
%               res = bpass( image_array, lnoise, lobject )
% INPUTS:
%               image:  The two-dimensional array to be filtered.
%               lnoise: Characteristic lengthscale of noise in pixels.
%                       Additive noise averaged over this length should
%                       vanish. May assume any positive floating value.
%                       May be set to 0 or false, in which case only the
%                       highpass "background subtraction" operation is 
%                       performed.
%               lobject: (optional) Integer length in pixels somewhat 
%                       larger than a typical object. Can also be set to 
%                       0 or false, in which case only the lowpass 
%                       "blurring" operation defined by lnoise is done,
%                       without the background subtraction defined by
%                       lobject.  Defaults to false.
%               threshold: (optional) By default, after the convolution,
%                       any negative pixels are reset to 0.  Threshold
%                       changes the threshhold for setting pixels to
%                       0.  Positive values may be useful for removing
%                       stray noise or small particles.  Alternatively, can
%                       be set to -Inf so that no threshholding is
%                       performed at all. 
%		      SH101104 -set this so is always zero
%
% OUTPUTS:
%               res:    filtered image.
% PROCEDURE:
%               simple convolution yields spatial bandpass filtering.
% NOTES:
% Performs a bandpass by convolving with an appropriate kernel.  You can
% think of this as a two part process.  First, a lowpassed image is
% produced by convolving the original with a gaussian.  Next, a second
% lowpassed image is produced by convolving the original with a boxcar
% function. By subtracting the boxcar version from the gaussian version, we
% are using the boxcar version to perform a highpass.
% 
% original - lowpassed version of original => highpassed version of the
% original
% 
% Performing a lowpass and a highpass results in a bandpassed image.
% 
% Converts input to double.  Be advised that commands like 'image' display 
% double precision arrays differently from UINT8 arrays.

% MODIFICATION HISTORY:
%               Written by David G. Grier, The University of Chicago, 2/93.
%
%               Greatly revised version DGG 5/95.
%
%               Added /field keyword JCC 12/95.
% 
%               Memory optimizations and fixed normalization, DGG 8/99.
%               Converted to Matlab by D.Blair 4/2004-ish
%
%               Fixed some bugs with conv2 to make sure the edges are
%               removed D.B. 6/05
%
%               Removed inadvertent image shift ERD 6/05
% 
%               Added threshold to output.  Now sets all pixels with
%               negative values equal to zero.  Gets rid of ringing which
%               was destroying sub-pixel accuracy, unless window size in
%               cntrd was picked perfectly.  Now centrd gets sub-pixel
%               accuracy much more robustly ERD 8/24/05
%
%               Refactored for clarity and converted all convolutions to
%               use column vector kernels for speed.  Running on my 
%               macbook, the old version took ~1.3 seconds to do
%               bpass(image_array,1,19) on a 1024 x 1024 image; this
%               version takes roughly half that. JWM 6/07
%
%       This code 'bpass.pro' is copyright 1997, John C. Crocker and 
%       David G. Grier.  It should be considered 'freeware'- and may be
%       distributed freely in its original form when properly attributed.  

if nargin < 3, lobject = false; end

normalize = @(x) x/sum(x);

image_array = double(image_array);

if lnoise == 0
  gaussian_kernel = 1;
else      
  gaussian_kernel = normalize(...
    exp(-((-ceil(5*lnoise):ceil(5*lnoise))/(2*lnoise)).^2));
end

if lobject  
  boxcar_kernel = normalize(...
      ones(1,length(-round(lobject):round(lobject))));
end
  
% JWM: Do a 2D convolution with the kernels in two steps each.  It is
% possible to do the convolution in only one step per kernel with 
%
  % gconv = conv2(gaussian_kernel',gaussian_kernel,image_array,'same');
  % bconv = conv2(boxcar_kernel', boxcar_kernel,image_array,'same');
% 
% but for some reason, this is slow.  The whole operation could be reduced
% to a single step using the associative and distributive properties of
% convolution:
%
  % filtered = conv2(image_array,...
  %   gaussian_kernel'*gaussian_kernel - boxcar_kernel'*boxcar_kernel,...
  %   'same');
%
% But this is also comparatively slow (though inexplicably faster than the
% above).  It turns out that convolving with a column vector is faster than
% convolving with a row vector, so instead of transposing the kernel, the
% image is transposed twice.

gconv = conv2(image_array',gaussian_kernel','same');
gconv = conv2(gconv',gaussian_kernel','same');

if lobject
  bconv = conv2(image_array',boxcar_kernel','same');
  bconv = conv2(bconv',boxcar_kernel','same');

  filtered = gconv - bconv;
else
  filtered = gconv;
end

% Zero out the values on the edges to signal that they're not useful.     
lzero = max(lobject,ceil(5*lnoise));

filtered(1:(round(lzero)),:) = 0;
filtered((end - lzero + 1):end,:) = 0;
filtered(:,1:(round(lzero))) = 0;
filtered(:,(end - lzero + 1):end) = 0;

% JWM: I question the value of zeroing out negative pixels.  It's a
% nonlinear operation which could potentially mess up our expectations
% about statistics.  Is there data on 'Now centroid gets subpixel accuracy
% much more robustly'?  To choose which approach to take, uncomment one of
% the following two lines.
% ERD: The negative values shift the peak if the center of the cntrd mask
% is not centered on the particle.
% res = filtered;
%filtered(filtered < threshold) = 0;
filtered(filtered < 0) = 0;
res = filtered;
%--------------------------------------------
%--------------------------------------------
%--------------------------------------------
function out=findLocalMax(im,th)
% finds local maxima in an image to pixel level accuracy.   
%  local maxima condition is >= rather than >
% inspired by Crocker & Griers algorithm, and also Daniel Blair & Eric Dufresne's implementation
%   im = input image for detection
%   th - detection threshold - pixels must be strictly greater than this value
% out : x,y coordinates of local maxima

%identify above threshold pixels
[y,x] = find(im>th);
%delete pixels identified at the boundary of the image
[imsizey, imsizex]= size(im);
edgepixel_idx = find( y==1 | y==imsizey | x==1 | x==imsizex);
y(edgepixel_idx) = [];
x(edgepixel_idx) = [];

%check if its a local maxima
subim = zeros(3,3);
islocalmaxima = zeros(numel(x),1);
for i = 1:numel(x)
  subim = im([y(i)-1:y(i)+1],[x(i)-1:x(i)+1]);
  islocalmaxima(i) = all(subim(2,2)>=subim(:));
end
%assign all above threshold pixels to out initially
out = [x,y];
%delete pixels which are not local maxima
out(find(~islocalmaxima),:) = [];
%-------------------------------------------------------------------------------------------
function [phot_count a normChi2 pos eccentricity ] = freeGaussFitEllipse( im, point_pos, windowSize,posLim, sigmaLim, initguess,varargin)
% function [phot_count pos normChi2 eccentricity a ] = freeGaussFitEllipse( im, point_pos, windowSize,posLim, sigmaLim, initguess)
% fit single 2d gaussian with fixed x, y position. wrapper to gauss fit tools function gaussfit_free_elliptical.cpp
%
% Inputs: 
%   im - fit image should be of type double
%   point_pos - initial estimate of position
%   windowSize  - size of subimage to crop 
%   posLim    - radius to allow shift from initial fit position
%   sigmaLim - [minwidth maxwidth]  - fit limits of psf width
%   initguess - [amplitude widthguess background X_POSim Y_POSim ];
% Outputs:
%   phot_count - volume of psf
%   a -  fit parameters, ie  [A, sigma_x, sigma_y, b, X, Y, theta ];
%    sigma_y is always the width along the major axis and theta is angle from y axis to the major axis
%
% NB ONLY ALLOWS SQUARE SUBIMAGES OTHERWISE RETURN 0
%
% Elliptical gaussian from cpp function:
%    xprime = (X-Xo)*cos(theta) - (Y-Yo)*sin(theta);
%    yprime = (X-Xo)*sin(theta) + (Y-Yo)*cos(theta);
%    e = exp((-(pow(xprime/xdenom,2)))-(pow(yprime/ydenom,2)));
%
%    x[element] = (A * e) + b;
% 
% Twotone TIRF-FRET image analysis software.
% Version 3.1.0 Alpha, released 101115
% Authors: Seamus J Holden, Stephan Uphoff
% Email: s.holden1@physics.ox.ac.uk
% Copyright (C) 2010, Isis Innovation Limited.
% All rights reserved.
% TwoTone is released under an “academic use only�? license; for details please see the accompanying ‘TWOTONE_LICENSE.doc’. Usage of the software requires acceptance of this license
%
n = numel(varargin);
i = 1;
%defaults:
useCPPfit = true;
useAutoDetTol = false;
while i <= n
  if strcmp(varargin{i},'autoDetect')
    useAutoDetTol = true;
    i = i+1;
  elseif strcmp(varargin{i},'useMatlabFit')
    useCPPfit = false;
    i = i+1;
  else
    error('Unrecognised argument');
  end
end

if useAutoDetTol == false
  % default convergence tolerances
  verbose = false;
  epsilon1 =  10^-7; %gradient on lsq -
  epsilon2 =  10^-9; %gradient on fitParam - the most important
  epsilon3 = 0;  %absoluteValueLSQ - problem dependent - usually ~10^3 - NEVER USE IT!
  maxIter  = 100; %how fast it is in the absence of signal
else
  % for autodetection, convergence tolerances can be reduced (sub nanometre localization not required ofr accurate photon count
  % default convergence tolerances
  verbose = false;
  epsilon1 =  10^-2; %gradient on lsq -
  epsilon2 =  10^-2; %gradient on fitParam - the most important
  epsilon3 = 0;  %absoluteValueLSQ - problem dependent - usually ~10^3 - NEVER USE IT!
  maxIter  = 20; %how fast it is in the absence of signal
end

[sizey sizex] = size(im);
X0=point_pos(1);
Y0=point_pos(2);

%round X0, Y0 to use as matrix locations
X0_int = round(X0); 
Y0_int = round(Y0);
windowSize = round(windowSize); %radius should already be an integer anyway

% setup the limits of the cropped image
xstart =  X0_int-windowSize;
xfinish = X0_int+windowSize;
ystart =  Y0_int-windowSize;
yfinish = Y0_int+windowSize;
% check if any of the limits are out of bounds - if so, skip that point
if (xstart<1) || (xstart > sizex) ||  (xfinish<1) || (xfinish > sizex) ...
  || (ystart<1) || (ystart > sizey) ||  (yfinish<1) || (yfinish > sizey) 
  
  a = [0 0 0 0 0 0 0]; 
  normChi2 = NaN;
  %warning('sub im outide limits');
else
  
  %crop to a small area around the point
  fitIm = im( ystart:yfinish, xstart:xfinish);
  [sizeyFit sizexFit] = size(fitIm);
  % set up the point location in the cropped image coords
  X_POSim = X0-xstart+1;
  Y_POSim = Y0-ystart+1;

  %set up the XY Lims
  xLim = [(X_POSim - posLim) (X_POSim + posLim)];
  yLim = [(Y_POSim - posLim) (Y_POSim + posLim)];
  minwidth = sigmaLim(1);
  maxwidth = sigmaLim(2);
  %if an intial guess is not supplied, calculate it
  if all(initguess==0)
    background = min(fitIm(:));
    amplitude = max(fitIm(:));
    widthguess = widthEstimate(fitIm)/2;
    if widthguess < minwidth
      widthguess = minwidth;
    elseif widthguess > maxwidth
      widthguess = maxwidth;
    end
      
    initguess = [amplitude widthguess background X_POSim Y_POSim ];
  end

  %make sure the point position is within the limits of the cropped image and that the image is square
  if ( X_POSim < 1) || ( X_POSim > sizexFit) || ( Y_POSim < 1) || ( Y_POSim > sizeyFit) || (sizexFit ~= sizeyFit)
    a = [0 0 0 0 0 0 0];
  normChi2 = NaN;
  else
    %do the fit 
    if useCPPfit == true %use the faster C++ fitting library
      %do the fit 
      optionVector = [verbose, epsilon1,epsilon2,epsilon3,maxIter];
      [a, normChi2] =freePosEllipseGaussFit_cpp(fitIm,initguess ,xLim,yLim, sigmaLim,optionVector);
    else %use the matlab version
      %do the fit 
      curvefitoptions = optimset( 'lsqcurvefit');
      curvefitoptions = optimset( curvefitoptions,'Jacobian' ,'on','Display', 'off',  'TolX', epsilon2, 'TolFun', epsilon1,'MaxPCGIter',1,'MaxIter',maxIter);
      [a, normChi2] =freePosEllipseGaussFit_matlab(fitIm,initguess ,xLim,yLim, sigmaLim,curvefitoptions);
    end
  end
end

%display( EXITFLAG);
%total photon count I = 2*pi*I0ab
% a: (A, sigma_x, sigma_y, b, Xo, Yo, theta )
I0 = a(1);
stdX = a(2);
stdY = a(3);
BG0 = a(4);
pos = [(a(5)+xstart-1) (a(6)+ystart-1)];
theta = a(7);
eccentricity = sqrt( 1 - (min(stdX,stdY)/max(stdX,stdY) )^2);
phot_count = 2*pi*stdX*stdY*I0;
%----------------------------------------------------------------------------------------------
function widthEst = widthEstimate(m)
%function to better estimate width for initial guess
[sizey sizex] = size(m);
vx = sum(m);
vy = sum(m');

vx = vx.*(vx>0);
vy = vy.*(vy>0);

x = [1:sizex];
y = [1:sizey];

cx = sum(vx.*x)/sum(vx);
cy = sum(vy.*y)/sum(vy);

sx = sqrt(sum(vx.*(abs(x-cx).^2))/sum(vx));
sy = sqrt(sum(vy.*(abs(y-cy).^2))/sum(vy));

widthEst = 0.5*(sx + sy) ;

%----------------------------------------------------------------------------------------------
function [fitParam, normChi2] = freePosEllipseGaussFit_cpp(im,initguess ,xLim,yLim, sigmaLim,optionVector);
%function [fitParam, normChi2] = freePosEllipseGaussFit_cpp(im,initguess ,xLim,yLim, sigmaLim);
%
% ******wrapper for OB gaussFit C++ functions written 091215******
%
% fits point spread function, 
%
% [brightness, fit parameters, lsq] = gaussfit_fixedposition(image, fixed value of Xo, fixed value of Yo, lower bound on
% sigma, upper bound on sigma, lower bound on Xo, upper bound on Xo, lower
% bound on Yo, upper bound on Yo, use supplied initial guess (true or
% false), initial guess as 5 element  - (A,sigma, b, Xo, Yo), options
% vector - (display output on console true or false, epsilon1, 2,3, max
% number of iterations))
%
% fitParam = [A, sigma_x, sigma_y, b, X, Y, theta ];
%
%  initguess = [amplitude widthguess background ];
% initGuess7Vector - (A, sigma_x, sigma_y, b, Xo, Yo, theta )

A0start = initguess(1);
BGstart = initguess(3);
widthStart = initguess(2);
xStart = initguess(4);
yStart = initguess(5);
xMin = xLim(1);
xMax = xLim(2);
yMin = yLim(1);
yMax = yLim(2);
sigmaMin = sigmaLim(1);
sigmaMax = sigmaLim(2);

AMPSCALEFACTOR =max(im(:))/100;
if AMPSCALEFACTOR <= 0 
  AMPSCALEFACTOR = 1;
end
%rescale the variables - to make magnitude of amplitude, background and width similar
im = im./AMPSCALEFACTOR;

useInitGuess = true;
%useInitGuess = false;

A0start = A0start./AMPSCALEFACTOR; %amplitude
BGstart = BGstart./AMPSCALEFACTOR; %backgound
if ( (initguess(2) < sigmaMin) || (initguess(2) > sigmaMax))%if given an out of bounds generate a
  initguess(2) = (sigmaMax+sigmaMin)/2; %sensible one
end
thetaStart = 0; 
initGuess7Vector = [A0start widthStart widthStart BGstart xStart yStart thetaStart];
% A, sigma_x, sigma_y, b, Xo, Yo, theta 


[brightness,fitParam, lsq] =gaussfit_free_elliptical(im, sigmaMin,sigmaMax,sigmaMin,sigmaMax,xMin,xMax,yMin,yMax , useInitGuess, initGuess7Vector, optionVector);

% modify fitParam to fit with "fitParam" output syntax from twotoneMain 
% fitParam = (A, sigma_x, sigma_y, b, Xo, Yo, theta )
if fitParam(1)< 0 %We know that negative amplitude values are patently unphysical so ignore them
  fitParam(1) = 0;
end

fitParam(1) = fitParam(1).*AMPSCALEFACTOR;%amplitude
fitParam(4) = fitParam(4).*AMPSCALEFACTOR;%background
normChi2 = lsq; %this is the lsq normalised per pixel ie lsq_TOT.(Npix)^-1

%----------------------------------------------------------------------------------------------
function [fitParam, normChi2] = freePosEllipseGaussFit_matlab(inputIm,initguess ,xLim,yLim, sigmaLim,curvefitoptions);
% fits point spread function, 
% F = (fitParam(1)*exp(    -(xprime).^2/(2*fitParam(2)^2)+(yprime).^2) /(2*fitParam(3)^2)   ) + fitParam(4))
%        
%           xprime = (X-Xo)*cos(theta) - (Y-Yo)*sin(theta);
%           yprime = (X-Xo)*sin(theta) + (Y-Yo)*cos(theta);
%
% extra fit params fitParam(7) = theta, X0 and Y0 are fitParam(5) and (6)
%

A0start = initguess(1);
BGstart = initguess(3);
widthStart = initguess(2);
xStart = initguess(4);
yStart = initguess(5);
xMin = xLim(1);
xMax = xLim(2);
yMin = yLim(1);
yMax = yLim(2);
sigmaMin = sigmaLim(1);
sigmaMax = sigmaLim(2);

%set up the mesh, size of the input image for use in fitting
[sizey sizex] = size(inputIm);
num_pixels = sizey*sizex;
[X,Y]= meshgrid(1:sizex,1:sizey);
grid = [X Y];

AMPSCALEFACTOR =max(inputIm(:))/100;
if AMPSCALEFACTOR <= 0 
  AMPSCALEFACTOR = 1;
end
%rescale the variables - to make magnitude of amplitude, background and width similar
inputIm = inputIm./AMPSCALEFACTOR;

% initguess input is  [amplitude widthguess background X_POSim Y_POSim ]
% initGuess7Vector output is [amplitude sx sy background X_POSim Y_POSim theta]
A0start = A0start./AMPSCALEFACTOR; %amplitude
BGstart = BGstart./AMPSCALEFACTOR; %backgound
if ( (initguess(2) < sigmaMin) || (initguess(2) > sigmaMax))%if given an out of bounds generate a
  initguess(2) = (sigmaMax+sigmaMin)/2; %sensible one- careful this isnt too close to true val tho
end
thetaStart = 0; 
initGuess7Vector = [A0start widthStart widthStart BGstart xStart yStart thetaStart];
% A, sigma_x, sigma_y, b, Xo, Yo, theta 

% Set fit limits on [amplitude widthx widthy background theta]
% dont set limits on theta but convert it to range 0->2pi afterwards
lb = [0 sigmaMin sigmaMin 0 xMin yMin -inf];
ub = [65535 sigmaMax sigmaMax 65535  xMax yMax inf];
    
%do the fit
try
  [fitParam, res] = ...
    lsqcurvefit(@(x, xdata) gauss2dw(x, xdata), ...
     initGuess7Vector ,grid ,inputIm ,...
      lb,ub,curvefitoptions);    
catch ME
  if strcmp(ME.identifier,'optim:snls:InvalidUserFunction') % supplied absolutely empty image!
    fitParam = [0 0 0 0 0 0 0];
    res = 0;
  else
    rethrow(ME);
  end
end

% modify fitParam to fit with "fitParam" output syntax from twotoneMain 
% fitParam = (A, sigma_x, sigma_y, b, Xo, Yo, theta )
if fitParam(1)< 0 %We know that negative amplitude values are patently unphysical so ignore them
  fitParam(1) = 0;
end

fitParam(1) = fitParam(1).*AMPSCALEFACTOR;%amplitude
fitParam(4) = fitParam(4).*AMPSCALEFACTOR;%background

fitParam(7) = mod(fitParam(7),2*pi);% transform theta onto the interval 0->2*pi
normChi2 = res/num_pixels;     
%-------------------------------------------------------------------------
%---------------------Fitting Subfunctions--------------------------------
%-------------------------------------------------------------------------

function [F J] = gauss2dw(a, data)
% Used by the curve fitter to calculate values for a 2d gaussian
% with the x & y standard deviations equal
% and with fixed positions
% a(1) - A0
% a(2) - sX
% a(3) - sY
% a(4) - B
% a(5) - Xpos
% a(6) - Ypos
% a(7) - theta

%Initialise everything
[sizey sizex] = size(data);
sizex = sizex/2;

F = zeros(sizey, sizex);
X = F;
Y = F;

X = data(:,1:sizex);
Y = data(:,sizex+1:end);

xprime = (X-a(5))*cos(a(7)) - (Y-a(6))*sin(a(7));
yprime = (X-a(5))*sin(a(7)) + (Y-a(6))*cos(a(7));

% Only evaluate the exponential once:
expPart = exp( - ((xprime).^2 /(2*a(2)^2) + (yprime).^2 /(2*a(3)^2) ));

F = a(1)*expPart + a(4);

% compute the jacobian

% initialise everything
n = numel(F);
J = zeros(n,3); % initialise J
Ga1F = zeros(sizey, sizex);% dF/da(1)
Ga2F = Ga1F;% dF/da(2)
Ga3F = Ga1F;% dF/da(3)
Ga4F = Ga1F;% dF/da(4)
Ga5F = Ga1F;% dF/da(7)

% Calculate the grad_a1(F),  grad_a2(F), etc

Ga1F = expPart;

Ga2F = a(1).* expPart .*xprime.^2 .*a(2).^-3;% (A * e) * (pow(xprime,2) * pow(sigma_x,-3)); //dF/dsigma_x
Ga3F = a(1).* expPart .*yprime.^2 .*a(3).^-3;% (A * e) * (pow(yprime,2) * pow(sigma_y,-3)); //dF/dsigma_y

Ga4F = ones(size(X));
%dF/dX0 and dF/dY0 in cpp notation
%   jac[j++] = (A * e) * ( (xprime*cos(theta)*pow(sigma_x,-2)) + (yprime*sin(theta)*pow(sigma_y,-2)) ); //dF/dXo
%           jac[j++] = (A * e) * ( (yprime*cos(theta)*pow(sigma_y,-2)) - (xprime*sin(theta)*pow(sigma_x,-2)) ); //dF/dYo
Ga5F = a(1).* expPart .* ...
      (  xprime.*a(2).^(-2).*cos(a(7)) + yprime.*a(3).^(-2)*sin(a(7)) );
Ga6F = a(1).* expPart .* ...
      ( -xprime.*a(2).^(-2).*sin(a(7)) + yprime.*a(3).^(-2)*cos(a(7)) );

%dF/da(7) in c++ notation:
% (-A * e) *( (-xprime * pow(sigma_x,-2)) * ((X-Xo)*sin(theta) + (Y-Yo)*cos(theta)) + (yprime*pow(sigma_y,-2))*((X-Xo)*cos(theta) - (Y-Yo)*sin(theta)) );
Ga7F = -a(1).* expPart.* ...
      (  (-xprime).*a(2).^(-2).*((X-a(5))*sin(a(7))+ (Y-a(6))*cos(a(7))) ... 
        + (yprime).*a(3).^(-2).*((X-a(5))*cos(a(7))- (Y-a(6))*sin(a(7))) );

% Form the jacobian, see the printed notes on getGaussFit for derivation
J = [Ga1F(:) Ga2F(:) Ga3F(:) Ga4F(:) Ga5F(:) Ga6F(:) Ga7F(:)];



