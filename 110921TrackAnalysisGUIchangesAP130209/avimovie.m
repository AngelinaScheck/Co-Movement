function avimovie(tracks,ti,MinSteps)
%use getframe to create .avi movie out of tracked particle
%got it working on 01/12/11 (but figure is still displayed)
%tracks' structure: xpos|ypos|frame|trackID
%code is working, but not speeded up 
[TIRFFilename, TIRFPathname] = uigetfile(...
   {'*.fits', 'TIRF data: (*.fits)';...
   '*.tif', 'TIRF data: (*.tif)'});

loadname = [TIRFPathname TIRFFilename];

if ~(isnumeric(TIRFFilename)&&TIRFFilename==0) %check the user has not pressed cancel
    if (tif == 0)
ImageInfo = fits_read_header([TIRFPathname TIRFFilename]);
%get image size
imageLim = [1 ImageInfo.NAXIS1 1 ImageInfo.NAXIS2];
%get data of .fits file as image data for each image in the stack        
DataIn = ImageStack(loadname, imageLim);
%get first image to calibrate intensity of image stack and to get the image size
I1 = getFrame(DataIn,1);
[H,W,D] = size(I1);
    else
        info = imfinfo(loadname);
        %load current Image
        I1 = double(imread(loadname, 1, 'Info', info));
        [H,W,D] = size(I1);  
    end
%get number of frames in the stack 
maxframes = max(tracks(:,3));
%get total number of tracked molecules
nMolecules = max(tracks(:,4));

%choose directory to save the movie after acquisition
[MovieFilename, MoviePathname] = uiputfile('*.avi', 'Save movie as: (*.avi)');
if ~(isnumeric(MovieFilename)&&MovieFilename==0) %check the user has not pressed cancel
savename = [MoviePathname MovieFilename];


%open avi object
aviobj = avifile(savename,'compression','None');

% prevent the figure window from appearing at all (not reall working)
f = figure('visible','off');
% set figure parameters
set(f, 'paperposition', [0 0 W H]);
set(f, 'papersize', [W H]);

for fra=1:maxframes
    %index in tracks array of tracked particles in frame (fra)
    clear indInFra;
    indInFra = find(tracks(:,3)==fra);
        if (tif==0)
    %get (fra) image/frame from stack
    I = getFrame(DataIn,fra);
        else
    I = double(imread(loadname,fra,'Info',info));
        end
    %flip image vertically, because tracked data is flip so as well
    I = flipud(I);
    imshow(I,[min(min(I1)) max(max(I1))], 'border','tight', 'InitialMagnification',100); hold on;
      
    for par=1:numel(indInFra)
        xpos = round(tracks(indInFra(par),1));
        ypos = round(tracks(indInFra(par),2));
        trackID = round(tracks(indInFra(par),4));
        %check,if particle occurs in more than four frames
        clear tracklength;
        tracklength = length(find(tracks(:,4)==trackID));
        if (tracklength>MinSteps)
        %check, if particle is to close to the edge of the frame
        if (5 < xpos) && (xpos < (W-5)) && (5 < ypos) && (ypos < (H-5))
        %circle localized particle
        plot(xpos,ypos,'ro'); hold on;
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
                    plot(linspace(xold,xpos,10),linspace(yold,ypos,10),'r'); hold on;
                    boolean = 1;
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
%save figure in avi movie
F = getframe(f);
hold off;
aviobj = addframe(aviobj,F);

if(mod(fra,100) == 0)
    disp('...................')
end
end
close(f);
aviobj = close(aviobj);
disp('all done.')
end
end

end