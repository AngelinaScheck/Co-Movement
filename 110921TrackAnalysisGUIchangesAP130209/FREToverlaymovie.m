
function FREToverlaymovie(loadnames, TFormfilename, FirstGreenFrame, tif, alttype, startframe)
% Form overlaid false coloured movies by transforming G frames and
% overlaying with red using using tform from transformMovies.m.  
% RC 19/06/12
% AP 04/02/13
% clearing up loading of images (more general, independent of used setup and FOV)


%load tform
load(TFormfilename);

%make sure loadnames is a cell, normally is unless only 1 file selected
%in which case it becomes just a string
loadnames = cellstr(loadnames);
num_files = size(loadnames,2);

%cycle over files to analyse
for j = 1:num_files
    
    filename = loadnames{j};
    fprintf(['\n\nNow analysing: \n' filename ])
    
    %check image limits dependent on tif or fits
    if (tif == 0)
        ImageInfo = fits_read_header(filename);
        %imageLim
        imageLim = [1 ImageInfo.NAXIS1 1 ImageInfo.NAXIS2];
        %Zugspitze imageLim -- be careful (think about blue channel)
        no_frames = ImageInfo.NAXIS3;
        ysize = double(ImageInfo.NAXIS2);%double(ImageInfo.NAXIS2);
    else
        info = imfinfo(filename);
        firstImage = double(imread(filename,1,'Info',info));
        %imageLim
        [H,W,D] = size(firstImage);
        imageLim = [1 W 1 H];
        %Zugspitze imageLim -- be careful (think about blue channel)
        no_frames = numel(info);
        ysize = double(H);
    end  

    data = ImageStack(filename, imageLim);
    xsize = double(round(imageLim(2)/2));   %235
    xend = imageLim(2);
    %Zugspitze imageLim -- be careful (think about blue channel)
    
    %pre-allocate memory for loop
    Gdata_tfm = zeros(ysize,xsize);
    dummy_blue = zeros(ysize,xsize);
    overlay = zeros(ysize,xsize);
    %overlay_all = zeros(ImageInfo.NAXIS2,xsize, 3, no_frames);
 

    %loop over frames of movie
    if alttype == 1;
        %CW case
        %get image pixel value limits - normalise to 0-1 (in frame loop)
        firstframe = double(getFrame(data,startframe));
        min_I = min(firstframe(:)); max_I = max(firstframe(:));
        
        for i = startframe:no_frames

              if (mod(i,21) == 0 )
                fprintf(strcat(num2str(no_frames - i),' frames left.\n'));
              else
                fprintf('.');
              end

            %get data for each frame from ImageStack & normalise between 0 and
            %255 using values from first frame to analyse
            
            frame = double(getFrame(data,i));
            frame = frame - min_I; frame = frame./(max_I - min_I);

            %split frame data into G/R channels
            Gdata = frame(1:ysize,xsize+1:xend);
            Rdata = frame(1:ysize,1:xsize); 

            %flip data to be in right orientation
            %Gdata = flipud(Gdata); Rdata = flipud(Rdata);
            %tranform G data to fit R data using transform
            Gdata_tfm = imtransform(Gdata, TFORM, 'FillValues', 0, 'Xdata', [1 xsize], 'Ydata', [1 ysize]);

            %make false coloured overlay of transformed G and R data
            overlay = cat(3, Rdata, Gdata_tfm, dummy_blue);
            %save overlays in case needed later
            %overlay_all(:,:,:,i) = overlay;

            %write overlays to tiffstack appending each time
            imwrite(overlay,[filename(1:end-5) '_overlayMovie.tif'],'tif','WriteMode','append','Compression','none');
        end
        fprintf('\nDone.\n');
        
    else
        %ALEX case
        DDdata_tfm = zeros(ysize,xsize);
        ADdata_tfm = zeros(ysize,xsize);
        Aex_overlay = zeros(ysize,xsize);
        Dex_overlay = zeros(ysize,xsize);
        %Dex_overlay_all = zeros(ImageInfo.NAXIS2,xsize, 3, no_frames);
        
        %increase index by 1 if we don't start on a green frame
        if (FirstGreenFrame == 0)
            startframe = startframe + 1;
        end

        firstframe = double(getFrame(data,startframe));
        min_I = min(firstframe(:)); max_I = max(firstframe(:));
        
        %loop over frames, finish on no_frames-1 in case don't end on red
        %(i.e. not an even number). This misses one frame if you do
        %end on red, but it's simpler.     
        for i = startframe:2:(no_frames-1);
            
              if (mod(i,21) == 0 )
                fprintf(strcat(num2str(no_frames - i),' frames left.\n'));
              else
                fprintf('.');
              end
              
            %get data for each frame (G/R) from ImageStack & normalise between 0 and
            %1 using values from first green frame to analyse
            greenframe = double(getFrame(data,i));
            greenframe = greenframe - min_I; greenframe = greenframe./(max_I - min_I);

            redframe = double(getFrame(data,i+1));
            redframe = redframe - min_I; redframe = redframe./(max_I - min_I);
            
            %split frame data into DD/DA/AD/AA channels
            DDdata = greenframe(1:ysize,(xsize+1):xend);
            DAdata = greenframe(1:ysize,1:xsize);
            ADdata = redframe(1:ysize,(xsize+1):xend);
            AAdata = redframe(1:ysize,1:xsize);

            %tranform G data to fit R data using transform
            DDdata_tfm = imtransform(DDdata,TFORM, 'FillValues', 0, 'Xdata', [1 xsize], 'Ydata', [1 ysize]);
            ADdata_tfm = imtransform(ADdata,TFORM, 'FillValues', 0, 'Xdata', [1 xsize], 'Ydata', [1 ysize]);
            
            %make false coloured overlay of transformed G and R data
            Dex_overlay = cat(3, DAdata, DDdata_tfm, dummy_blue);
            Aex_overlay = cat(3, AAdata, ADdata_tfm, dummy_blue);
            
            %save overlays in case needed later
            %overlay_all(:,:,:,i) = overlay;

            %write overlays to tiffstack appending each time
            imwrite(Dex_overlay,[filename(1:end-5) '_overlayMovie.tif'],'tif','WriteMode','append','Compression','none');
            %comment out next line if you just want green frames
            imwrite(Aex_overlay,[filename(1:end-5) '_overlayMovie.tif'],'tif','WriteMode','append','Compression','none');
        end
        fprintf('\nDone.\n');
    end
end
