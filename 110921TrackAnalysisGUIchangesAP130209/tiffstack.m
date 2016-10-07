function tiffstack(tracks,tif,MinSteps)
%trying to get the path for each molecule but only opening one frame- 
%got it working 28/11/11
%not speeded up up to now
%try O and strait line 06/12/11
%this would remove one while-loop and one if-clause
%therefore get functions Circle and StraitLine implemented
%last modified (color code) on 22/02/12
%Anne Plochowietz

[TIRFFilename, TIRFPathname] = uigetfile(...
   {'*.fits', 'TIRF data: (*.fits)';...
   '*.tif', 'TIRF data: (*.tif)'});

loadname = [TIRFPathname TIRFFilename];

if ~(isnumeric(TIRFFilename)&&TIRFFilename==0) %check the user has not pressed cancel
    if (tif ==0)
ImageInfo = fits_read_header([TIRFPathname TIRFFilename]);
%get image size
imageLim = [1 ImageInfo.NAXIS1 1 ImageInfo.NAXIS2];
%get data of .fits file as image data for each image in the stack        
DataIn = ImageStack(loadname, imageLim);
    else
info = imfinfo(loadname);
%get length of image stack- nFrame
nFrame = numel(info);
I1 = double(imread(loadname, 1, 'Info', info));
[H,W,D] = size(I1);
imageLim = [1 W 1 H];
    end    
%get number of frames in the stack              
maxframes = max(tracks(:,3));
%get total number of tracked molecules
nMolecules = max(tracks(:,4));

%choose directory to save the movie after acquisition
[MovieFilename, MoviePathname] = uiputfile('*.tif', 'Save movie as: (*.tif)');
if ~(isnumeric(MovieFilename)&&MovieFilename==0) %check the user has not pressed cancel
savename = [MoviePathname MovieFilename];

%Initialisation of the colormap for TIRF images using the first image
    if (tif ==0)
I = getFrame(DataIn,1);
    else
I = imread(loadname,1,'Info',info);
    end
map = zeros(max(max(I))+1,3); %colourmatrix
colourrange = single(max(max(I))-min(min(I))); %difference between largeset and smallest intensity value
offset = single(min(min(I))); % smallest colour value
for count = min(min(I)):max(max(I))
    map(count+1,:)=(single(count)-offset)/colourrange;
end
%for cross/tracks colour use ten different colours from jet colormap
cmap = colormap(jet(10));
%Initialization
OfInterest = zeros(numel(tracks(:,4)),12);
%filterParams(1) = 40;%threshold
%filterParams(2) = 1;%lnoise
filterParams(3) = 7;%lobject
windowSize = filterParams(3);
initguess = [];
posThresh = windowSize/2; % if the pos moves more than windowSize/2 from initialPos, drop it
sigmaLim = [0 inf];

counter = 1;
for fra=1:maxframes
    %index in tracks array of tracked particles in frame (fra)
    clear indInFra;
    indInFra = find(tracks(:,3)==fra);
        if (tif == 0)
    %get (fra) image/frame from stack
    I = double(getFrame(DataIn,fra));
        else
    I = double(imread(loadname, fra, 'Info', info));
        end
    %introduced to adjust image contrast per frame independently
    map = zeros(max(max(I))+1,3); %colourmatrix
    colourrange = single(max(max(I))-min(min(I))); %difference between largeset and smallest intensity value
    offset = single(min(min(I))); % smallest colour value
        for count = min(min(I)):max(max(I))
            map(count+1,:)=(single(count)-offset)/colourrange;
        end    
    %flip image vertically, because tracked data is flip so as well
    I = flipud(I);
    %xy = findLocalMax(I, filterParams(1));
    %convert to grayscale image to rgb
    Irgb = ind2rgb(I,map);
      
    for par=1:numel(indInFra)
        %take positions from track
        xpos = round(tracks(indInFra(par),1));
        ypos = round(tracks(indInFra(par),2));
        trackID = round(tracks(indInFra(par),4));
        
        %check,if particle occurs in more than four frames
        clear tracklength;
        tracklength = numel(find(tracks(:,4)==trackID));
        if (tracklength>MinSteps)
        %check, if particle is to close to the edge of the frame
        if (5 < xpos) && (xpos < (imageLim(2)-5)) && (5 < ypos) && (ypos < (imageLim(4)-5))
            %write coloured circles directly into Irgb
            [rowCircle,colCircle] = find(Circle(xpos,ypos));
            %only for a test
            for ina1=1:length(rowCircle)
            Irgb(rowCircle(ina1),colCircle(ina1),1) = cmap(mod(trackID,10)+1,1);
            Irgb(rowCircle(ina1),colCircle(ina1),2) = cmap(mod(trackID,10)+1,2);
            Irgb(rowCircle(ina1),colCircle(ina1),3) = cmap(mod(trackID,10)+1,3);
            end
    %         Irgb(ypos,(xpos-5):(xpos+5),1) = cmap(mod(tracks(trackID,4),10)+1,1);
    %         Irgb(ypos,(xpos-5):(xpos+5),2) = cmap(mod(tracks(trackID,4),10)+1,2);
    %         Irgb(ypos,(xpos-5):(xpos+5),3) = cmap(mod(tracks(trackID,4),10)+1,3);
    %         
    %         Irgb((ypos-5):(ypos+5),xpos,1) = cmap(mod(tracks(trackID,4),10)+1,1);
    %         Irgb((ypos-5):(ypos+5),xpos,2) = cmap(mod(tracks(trackID,4),10)+1,2);
    %         Irgb((ypos-5):(ypos+5),xpos,3) = cmap(mod(tracks(trackID,4),10)+1,3);
            
            %my photon count
            [rowTotalCircle,colTotalCircle] = find(TotalCircle(xpos,ypos));
            my_phot_count = sum(sum(I(rowTotalCircle,colTotalCircle)));
           
            %fit gaussian again for Io, BG and phot_count
            xy = [xpos ypos];
            [phot_count, a, normChi2, posFinal, eccentricity] = freeGaussFitEllipse(I, xy, windowSize,posThresh, sigmaLim, initguess,'autoDetect');
                
            %[phot_count, a, normChi2, posFinal, eccentricity] = freeGaussFitEllipse(I, xy(par,:), windowSize,posThresh, sigmaLim, initguess,'autoDetect');
            OfInterest(counter,:) = [fra, trackID, my_phot_count, phot_count, a(1), a(4), posFinal(1), posFinal(2), a(2), a(3), eccentricity, a(7)];
            counter = counter + 1;
            
            xframe = fra-1;
            while (xframe > 0)
                clear indInXframe;
                indInXframe = find(tracks(:,3)==xframe);
                parX = 1;
                boolean = 0;
                while (parX < (numel(indInXframe)+1))
                    if (tracks(indInXframe(parX),4) == trackID)
                        xold = round(tracks(indInXframe(parX),1));
                        yold = round(tracks(indInXframe(parX),2));
                        %draw line between tracked particles

                        [rowLine,colLine] = find(StraitLine(xold,yold,xpos,ypos));
                        %only for a test
                        for ina2=1:length(rowLine)
                        Irgb(rowLine(ina2),colLine(ina2),1) = cmap(mod(trackID,10)+1,1);
                        Irgb(rowLine(ina2),colLine(ina2),2) = cmap(mod(trackID,10)+1,2);
                        Irgb(rowLine(ina2),colLine(ina2),3) = cmap(mod(trackID,10)+1,3);
                        end
                        boolean = 1;
    %                     xlength = xpos - xold;%get how far they have moved from (fra-1) to (fra)
    %                     ylength = ypos - yold;
    %                     boolean = 1;
    %                     while (abs(xlength) > 0) || (abs(ylength) > 0)
    %                     if (abs(xlength) > abs(ylength))
    %                         xpos = round(xpos-sign(xlength));
    %                         if (xpos > 0) && (ypos > 0)
    %                         Irgb(ypos,xpos,1) = cmap(mod(tracks(trackID,4),10)+1,1);
    %                         Irgb(ypos,xpos,2) = cmap(mod(tracks(trackID,4),10)+1,2);
    %                         Irgb(ypos,xpos,3) = cmap(mod(tracks(trackID,4),10)+1,3);
    %                         end
    %                         xlength = xlength-sign(xlength);
    %                     else
    %                         ypos = round(ypos-sign(ylength));
    %                         if (xpos > 0) && (ypos > 0)
    %                         Irgb(ypos,xpos,1) = cmap(mod(tracks(trackID,4),10)+1,1);
    %                         Irgb(ypos,xpos,2) = cmap(mod(tracks(trackID,4),10)+1,2);
    %                         Irgb(ypos,xpos,3) = cmap(mod(tracks(trackID,4),10)+1,3);
    %                         end
    %                         ylength = ylength-sign(ylength);
    %                     end
    %                     end
                        xpos = xold;
                        ypos = yold;
                        parX = numel(indInXframe)+1;
                        xframe = xframe -1;
                    end
                parX = parX+1;
                end
                if (boolean == 0)
                    xframe = 0;
                end
            end
            
        end
        end
    end
%save images in a 'tif' stack, use WriteMode append
imwrite(Irgb, savename ,'WriteMode','append','Compression','none'); %perhaps add description to get iminfo afterwards

if(mod(fra,100) == 0)
    disp('...................')
end
end
save([savename(1:end-4) '_Idata.mat'],'OfInterest');
disp('all done.')
end
end
%keyboard;
%--------------------------------------------------------------------------
function sparseCircle = Circle(x,y)
%function that returns sparse matrix with ones at positions where rgb
%values should change, use find(Circle(xpos,ypos)) to get the positions of circle,
%go into Irgb with this
%Circle radius = 4 pixels (like in plot command)

jj = [x,x,x+1,x+1,x+2,x+2,x+3,x+3,x+3,x+3,x+4,x+4,x+4 ...
          x-1,x-1,x-2,x-2,x-3,x-3,x-3,x-3,x-4,x-4,x-4];
ii = [y+4,y-4,y+4,y-4,y+3,y-3,y+3,y-3,y+2,y-2,y-1,y,y+1 ...
              y+4,y-4,y+3,y-3,y+3,y-3,y+2,y-2,y-1,y,y+1];

sparseCircle = sparse(ii,jj,1);

%--------------------------------------------------------------------------
function totalCircle = TotalCircle(x,y)
%function that returns sparse matrix with ones at positions where rgb
%values should change, use find(Circle(xpos,ypos)) to get the positions of circle,
%go into Irgb with this
%Circle radius = 4 pixels
%Anne Plochowietz 2012-07-13

jj = [x,x,x,x,x,x,x,x,x, ...
        x+1,x+1,x+1,x+1,x+1,x+1,x+1,x+1,x+1, ...
        x-1,x-1,x-1,x-1,x-1,x-1,x-1,x-1,x-1, ...
        x+2,x+2,x+2,x+2,x+2,x+2,x+2, ...
        x-2,x-2,x-2,x-2,x-2,x-2,x-2, ...
        x+3,x+3,x+3,x+3,x+3,x+3,x+3, ...
        x-3,x-3,x-3,x-3,x-3,x-3,x-3, ...
        x+4,x+4,x+4, ...
        x-4,x-4,x-4];

ii = [y+4,y-4,y+3,y-3,y+2,y-2,y+1,y-1,y, ...
        y+4,y-4,y+3,y-3,y+2,y-2,y+1,y-1,y, ...
        y+4,y-4,y+3,y-3,y+2,y-2,y+1,y-1,y, ...
        y+3,y-3,y+2,y-2,y+1,y-1,y, ...
        y+3,y-3,y+2,y-2,y+1,y-1,y, ...
        y+3,y-3,y+2,y-2,y+1,y-1,y, ...
        y+3,y-3,y+2,y-2,y+1,y-1,y, ...
        y+1,y-1,y, ...
        y+1,y-1,y];

totalCircle = sparse(ii,jj,1);


%--------------------------------------------------------------------------
function sparseStraitLine = StraitLine(xold,yold,xpos,ypos)
%function that returns sparse matrix with ones at positions where rgb
%values should change, use find(StraitLine(xold,yold,xpos,ypos)) to get the positions of strait line,
%go into Irgb with this

m = (ypos-yold)/(xpos-xold);
n = (yold*xpos-ypos*xold)/(xpos-xold);
xlength = xpos-xold;
ylength = ypos-yold;

if(abs(xlength) >= abs(ylength))
x = round(xpos -sign(xlength)*[1:abs(xlength)]);
y = round(m*x + n);
elseif (abs(xlength) < abs(ylength))
y = round(ypos -sign(ylength)*[1:abs(ylength)]);
x = round((y-n)/m);%problem, x = -inf/inf = NaN might appear
    if (isnan(x) ==1)
    x = xpos;    
    end    
end    
sparseStraitLine = sparse(y,x,1);

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