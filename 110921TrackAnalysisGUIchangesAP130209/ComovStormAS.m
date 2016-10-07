function [savenameGreen, savenameRed] = ComovStormAS(fileName, Greenthreshold, Redthreshold, lobject, firstGreenFrame, TFORMpath, alternationPeriod,tif, alttype, startframe, leftRim)
%Angelina Scheck October 2016, modified from FRETSTORMAPMODIFED. Localize
%in both channels individually, track individually, then link tracks with
%tform --> comovement analysis to check if construct is intact

%output: Colocalization value over time. Number of co-moving tracks,
%number of individual tracks seperated by channel => ratio of intact
%constructs vs deagradation

%Only for CW, error if ALEX (might be updated in later versions)
assert(alttype==1,'comovement analysis only in CW mode');

%save localizations seperated by channel
fprintf('\n\nCW localizations for comovement analysis for file %s\n',fileName);
savenameGreen    = [fileName(1:end-5),'_thresh',num2str(Greenthreshold),'_lobject',num2str(lobject),'_startframe', num2str(startframe),'_CW-green_Comov.out'];    

savenameRed    = [fileName(1:end-5),'_thresh',num2str(Redthreshold),'_lobject',num2str(lobject),'_startframe', num2str(startframe),'_CW-red_Comov.out'];    
%definitions for the output file
saveFileHeader = 'FRAME\tXCENTER\tYCENTER\tBRIGHTNESS\tBG\tI0\tS_X\tS_Y\tTHETA\tECC\tE\tLink\n';


%image properties
if (tif == 0)
    ImageInfo = fits_read_header(fileName);
%     imageLim: green channel; red channel
    imageLim = [double(round(ImageInfo.NAXIS1/2+1)) 512 1 ImageInfo.NAXIS2; 1 double(round(ImageInfo.NAXIS1/2)) 1 ImageInfo.NAXIS2];
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
    
    nFrame = numel(info); %nFrame = getNumFrames(green_stack);
end 

%Translate images for analysis with aid of tform
tirfIm = TirfImage(fileName,firstGreenFrame, imageLim, 1);

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

%open the output files for writing
fid(1) = fopen(savenameGreen,'w');
fprintf(fid(1), saveFileHeader);
fid(2) = fopen(savenameRed,'w');
fprintf(fid(2), saveFileHeader);

countGreen = 0;
countRed = 0;

%Copy green and red localizations for each frame
for nn = startframe:1:nFrame

      %Update status
      if (mod(nn,30) == 0 )
            fprintf(strcat(num2str(nFrame - nn),' frames left.\n'));
      else
        fprintf('.');
      end

      GreenIm = double(getFrame( green_stack, nn));
      RedIm = double(getFrame( red_stack, nn));
      
      %transform green frame to fit red
      GreenIm_tfm = imtransform(GreenIm, tform, 'FillValues', 0, 'Xdata', [1 xsize], 'Ydata', [1 ysize]);

%       %sum DA and DD (transformed) to find positions
%       DDDAIm_sum = DAIm + DDIm_tfm;

      %filter
      GreenIm_tfmFilt = bpass(GreenIm_tfm,1,lobject);
      RedImFilt = bpass(RedIm,1,lobject);
      
      % Localize for both channels
%       threshCW = thresh * 2; %because we have summed two images.

      RedPos=findLocalMax(RedImFilt(:,leftRim:imageLim(2,2)),Redthreshold); %get rid off left rim after bpass filtering
      GreenPos=findLocalMax(GreenIm_tfmFilt(:,leftRim:imageLim(2,2)),Greenthreshold); %get rid off left rim after bpass filtering

      % find corresponding DD position guesses with tform from initial guesses in 
      % red channel
      %DDPos = tforminv(tform,DDDAPos);

      %set fitting parameters, filtering by shape
      initguess = [];
      RedPosThresh = windowSize/2; % if the pos moves more than windowSize/2 from initialPos, drop it
      GreenPosThresh = windowSize/2; % if the pos moves more than windowSize/2 from initialPos, drop it
      sigmaLim = [0 inf]; % Sigma can be filtered post-analysis        
%       DDposThresh = windowSize/100;  % only taking error with transform into account and force fit
      
%       %transform green frame to fit red channel
%       DDIm_tfm = imtransform(DDIm, tform, 'FillValues', 0, 'Xdata', [1 xsize], 'Ydata', [1 ysize]);
%     
% loop over positions for red and green (different sizes)
     for j = 1:size(RedPos,1)

        %do the red fit // WHY NOT CONTINUED WITH FILTERED?
        [Redphot_count, Reda, RednormChi2, RedposFinal, Redeccentricity] = freeGaussFitEllipse( RedIm, RedPos(j,:), windowSize,RedPosThresh, sigmaLim, initguess,'autoDetect');

%         DDposFinal = tforminv(tform,DDposFinal);    %transform fitted positions in green channel
       
        %Write localization file FORMAT NOT CHANGED, SO THAT IT CAN BE READ
        %BY THE OTHER FUNCTIONS IN THE GUI
        if ~all(Reda ==0)%if the fit has worked
          countRed = countRed +1;  
          %for Red
          %need to write following:'FRAME\tXCENTER\tYCENTER\tPHOTONS\tBG\tI0\tS_X\tS_Y\tTHETA\tECC\tE\tLink\n';
          fprintf(fid(2),'%d\t',nn); %FRAME
          fprintf(fid(2),'%6.3f\t',RedposFinal(1)); %XCENTER
          fprintf(fid(2),'%6.3f\t',RedposFinal(2)); %YCENTER
          fprintf(fid(2),'%6.3f\t',Redphot_count); %PHOTONS
          fprintf(fid(2),'%6.3f\t',Reda(4)); %BG
          fprintf(fid(2),'%6.3f\t',Reda(1)); %I0
          fprintf(fid(2),'%6.3f\t',Reda(2)); %S_X
          fprintf(fid(2),'%6.3f\t',Reda(3)); %S_Y
          fprintf(fid(2),'%6.3f\t',Reda(7)); %THETA
          fprintf(fid(2),'%6.3f\t',Redeccentricity); %ECC
          fprintf(fid(2),'%6.3f\t',NaN); %E
          fprintf(fid(2),'%6.3f\n',countRed); %Link

        end
     end
        
     for j = 1:size(GreenPos,1)
       
        %do the green fit
        [Greenphot_count, Greena, GreennormChi2, GreenposFinal, Greeneccentricity] = freeGaussFitEllipse( GreenIm_tfm, GreenPos(j,:), windowSize,GreenPosThresh, sigmaLim, initguess,'autoDetect');

%         DDposFinal = tforminv(tform,DDposFinal);    %transform fitted positions in green channel
       
        %Write localization file FORMAT NOT CHANGED, SO THAT IT CAN BE READ
        %BY THE OTHER FUNCTIONS IN THE GUI
        if ~all(Greena ==0)%if the fit has worked
          countGreen = countGreen +1;  
          
          %need to write following:'FRAME\tXCENTER\tYCENTER\tPHOTONS\tBG\tI0\tS_X\tS_Y\tTHETA\tECC\tE=NaN\tLink\n';
          fprintf(fid(1),'%d\t',nn); %FRAME
          fprintf(fid(1),'%6.3f\t',GreenposFinal(1)); %XCENTER
          fprintf(fid(1),'%6.3f\t',GreenposFinal(2)); %YCENTER
          fprintf(fid(1),'%6.3f\t',Greenphot_count); %PHOTONS
          fprintf(fid(1),'%6.3f\t',Greena(4)); %BG
          fprintf(fid(1),'%6.3f\t',Greena(1)); %I0
          fprintf(fid(1),'%6.3f\t',Greena(2)); %S_X
          fprintf(fid(1),'%6.3f\t',Greena(3)); %S_Y
          fprintf(fid(1),'%6.3f\t',Greena(7)); %THETA
          fprintf(fid(1),'%6.3f\t',Greeneccentricity); %ECC
          fprintf(fid(1),'%6.3f\t',NaN); %E
          fprintf(fid(1),'%6.3f\n',countGreen); %Link

        end
      end
end
    

fclose(fid(1));
fclose(fid(2));
fprintf('\nDone with localizations.\n');







%Now come helper functions directly copied from FRETSTORMAP without review



%----------------------------------------------------------------------------------------------
function backspace(nchar)
%function backspace(nchar)

for i = 1:nchar
  fprintf('\b');
end
%------------------------------------------------------------------------------

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


