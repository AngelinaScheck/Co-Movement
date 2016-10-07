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