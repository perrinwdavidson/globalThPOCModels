function coordm = deg2m(coord)
%--------------------------------------------------------------------------
% about:
% convert degrees of latitude and longitude to meters
%
% inputs:
% coord - coordinate array in degrees
%
% outputs:
% coordm - coordinate array in meters
%--------------------------------------------------------------------------

    coordm = deg2km(coord) .* 1000;

end