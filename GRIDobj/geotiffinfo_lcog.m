function info = geotiffinfo_lcog(fileOrURL)
%GEOTIFFINFO Information about GeoTIFF file
%
%   INFO = GEOTIFFINFO(FILENAME) returns a structure whose fields contain
%   image properties and cartographic information about a GeoTIFF file.
%
%   FILENAME is a string scalar or character vector that specifies the name
%   of the GeoTIFF file. Include the folder name in FILENAME or place the
%   file in the current folder or in a folder on the MATLAB path.  If the
%   named file includes the extension '.TIF' or '.TIFF'  (either upper or
%   lower case), you can omit the extension from FILENAME.
%
%   If the named file contains multiple GeoTIFF images, INFO is a structure
%   array with one element for each image in the file.  For example,
%   INFO(3) would contain information about the third image in the file. If
%   multiple images exist in the file, it is assumed that each image has
%   the same cartographic information and the same image width and height.
%
%   INFO = GEOTIFFINFO(URL) reads the GeoTIFF image from a URL. The URL
%   must include the protocol type (e.g., "http://").
%
%   The INFO structure contains these fields:
%
%   Filename       Character vector containing the name of the file or URL.
%
%   FileModDate    Character vector containing the modification date of the
%                  file.
%
%   FileSize       Integer indicating the size of the file in bytes.
%
%   Format         Character vector containing the file format, which 
%                  should always be 'tiff'.
%
%   FormatVersion  Character vector or number specifying the file format
%                  version.
%
%   Height         Integer indicating the height of the image in pixels.
%
%   Width          Integer indicating the width of the image in pixels.
%
%   BitDepth       Integer indicating the number of bits per pixel.  
%
%   ColorType      Character vector indicating the type of image:
%                     'truecolor' for a truecolor (RGB) image,
%                     'grayscale' for a grayscale image,
%                     'indexed'   for an indexed image.
%
%   ModelType      Character vector indicating the type of coordinate 
%                  system used to georeference the image:
%                     'ModelTypeProjected', 'ModelTypeGeographic', 
%                     'ModelTypeGeocentric' or ''.
%
%   PCS            Character vector indicating the projected coordinate
%                  system.
%
%   Projection     Character vector describing the EPSG identifier for the 
%                  underlying projection method.
%
%   MapSys         Character vector indicating the map system, if
%                  applicable:
%                     'STATE_PLANE_27', 'STATE_PLANE_83', 
%                     'UTM_NORTH', 'UTM_SOUTH', or ''.
%
%   Zone           Double indicating the UTM or State Plane Zone number, 
%                  empty ([]) if not applicable or unknown.
%
%   CTProjection   Character vector containing the GeoTIFF identifier for 
%                  the underlying projection method.
%
%   ProjParm       N-by-1 double containing projection parameter values.
%                  The identify of each element is specified by the 
%                  corresponding element of ProjParmId. Lengths are in
%                  meters, angles in decimal degrees.
%
%   ProjParmId     N-by-1 cell array listing the projection parameter 
%                  identifier for each corresponding numerical element 
%                  ProjParm.  The possible values are:
%                     'ProjNatOriginLatGeoKey'   
%                     'ProjNatOriginLongGeoKey'
%                     'ProjFalseEastingGeoKey'   
%                     'ProjFalseNorthingGeoKey'
%                     'ProjFalseOriginLatGeoKey' 
%                     'ProjFalseOriginLongGeoKey'
%                     'ProjCenterLatGeoKey'      
%                     'ProjCenterLongGeoKey'
%                     'ProjAzimuthAngleGeoKey'   
%                     'ProjRectifiedGridAngleGeoKey'
%                     'ProjScaleAtNatOriginGeoKey'
%                     'ProjStdParallel1GeoKey'   
%                     'ProjStdParallel2GeoKey'
%
%   GCS            Character vector indicating the geographic coordinate
%                  system.
%
%   Datum          Character vector indicating the projection datum type,
%                  such as
%                     'North American Datum 1927' or 
%                     'North American Datum 1983'.
%
%   Ellipsoid      Character vector indicating the ellipsoid name as 
%                  defined by the ellipsoid EPSG table.
%
%   SemiMajor      Double indicating the length of the semimajor axis 
%                  of the ellipsoid, in meters.
%
%   SemiMinor      Double indicating the length of the semiminor axis
%                  of the ellipsoid, in meters.
%
%   PM             Character vector indicating the prime meridian location,
%                  for example, 'Greenwich' or 'Paris'.
%
%   PmLongToGreenwich Double indicating the decimal degrees of longitude 
%                  between this prime meridian and Greenwich.  Prime 
%                  meridians to the west of Greenwich are negative.
%
%   UOMLength      Character vector indicating the units of length used in 
%                  the projected coordinate system.
%
%   UOMLengthInMeters Double defining the UOMLength unit in meters.
%
%   UOMAngle       Character vector indicating the angular units used for 
%                  geographic coordinates.
%
%   UOMAngleInDegrees Double defining the UOMAngle unit in degrees.
%
%   TiePoints      Structure containing the image tiepoints.
%                  The structure contains these fields:
%
%                  ImagePoints  Structure containing row and column
%                               coordinates of each image tiepoint. 
%                               The structure contains these fields:
%                                  Row  Double array of size 1-by-N.
%                                  Col  Double array of size 1-by-N. 
%
%                  WorldPoints  Structure containing the x and y
%                               world coordinates of the tiepoints.
%                               The structure contains these fields:
%                                  X    Double array of size 1-by-N.
%                                  Y    Double array of size 1-by-N.
%                 
%   PixelScale     3-by-1 array that specifies the X,Y,Z pixel scale 
%                  values.
%
%   SpatialRef     Map raster reference object if ModelType is 'ModelTypeProjected',
%                  or a geographic raster reference object if ModelType is
%                  'ModelTypeGeographic'.  If ModelType is empty, a warning
%                  is issued and SpatialRef is a map raster reference object.
%                  SpatialRef is empty if ModelType is 'ModelTypeGeocentric',
%                  or if the spatial referencing is ambiguously defined by
%                  the GeoTIFF file, or if ModelType is 'ModelTypeGeographic'
%                  and the geometric transformation type is 'affine'.
%
%   RefMatrix      3-by-2 double referencing matrix that must be 
%                  unambiguously defined by the GeoTIFF file. Otherwise it
%                  is empty.
%
%   BoundingBox    2-by-2 double array that specifies the minimum (row 1)  
%                  and maximum (row 2) values for each dimension of the
%                  image data in the GeoTIFF file.
%
%   CornerCoords   Structure with six fields that contains coordinates of
%                  the outer corners of the GeoTIFF image. Each field is a
%                  1-by-4 double array, or empty ([]) if unknown. The
%                  arrays contain the coordinates of the outer corners of
%                  the corner pixels, starting from the (1,1) corner and
%                  proceeding clockwise:
%
%                  X    Easting coordinates in the projected coordinate 
%                       system. X equals Lon (below) if the model type is
%                       'ModelTypeGeographic'.
%
%                  Y    Northing coordinates in the projected coordinate 
%                       system. Y equals Lat if the model type is 
%                       'ModelTypeGeographic'.
%
%                  Row  Row coordinates of the corner. 
%
%                  Col  Column coordinates of the corner.
%
%                  Lat  Latitudes of the corner. 
%
%                  Lon  Longitudes of the corner. 
%
%   GeoTIFFCodes   Structure containing raw numeric values for those
%                  GeoTIFF fields that are encoded numerically in the file.
%                  These raw values, converted to a character vector
%                  elsewhere in the INFO structure, are provided here for
%                  reference.
%
%                  The GeoTIFFCodes are:
%                     Model
%                     PCS
%                     GCS
%                     UOMLength
%                     UOMAngle
%                     Datum
%                     PM
%                     Ellipsoid
%                     ProjCode
%                     Projection
%                     CTProjection
%                     ProjParmId    
%                     MapSys
%
%                  Each is scalar except for ProjParmId which is a column
%                  vector.
%
%   GeoTIFFTags    Structure containing field names that match the GeoTIFF
%                  tags found in the file. At least one GeoTIFF tag must be
%                  present in the file or an error is issued.
%                  The following fields may be included:
%
%                     ModelPixelScaleTag           1-by-3 double
%                     ModelTiepointTag             1-by-6 double
%                     ModelTransformationTag       1-by-16 double
%                     GeoKeyDirectoryTag           scalar structure
%                     GeoAsciiParamsTag            character vector
%                     GeoDoubleParamsTag           1-by-N double
%                     RPCCoefficientTag            scalar map.geotiff.RPCCoefficientTag
%
%                  The GeoKeyDirectoryTag contains field names that match
%                  the names of the "GeoKeys". For more information about
%                  the "GeoKeys" refer to the GeoTIFF specification at
%                  http://geotiff.maptools.org/spec/geotiff6.html#6.2
%
%                  The RPCCoefficientTag contains properties with names
%                  corresponding to the tag elements listed in the 
%                  "RPCs in GeoTIFF" technical note at:
%                  http://geotiff.maptools.org/rpc_prop.html
%
%   ImageDescription  Character vector describing the image. If no 
%                     description is included in the file, the field is
%                     omitted.
%
%   Example
%   -------
%   info = geotiffinfo('boston.tif')
%
%   See also GEOTIFFWRITE, map.geotiff.RPCCoefficientTag, PROJFWD, PROJINV, READGEORASTER

% Copyright 1996-2020 The MathWorks, Inc.

% 
% Modified for reading Landsat COG data, which may have multiple tiff_infos(see L.279&280 in this code)
% This would require the function constructGeoKeyDirectoryMap at the end.
% by Shi July 14, 2021

% Validate the input. If it's a filename, return the full pathname. If it's
% a URL, download the URL to a temporary file and set isURL to true.
if nargin > 0
    fileOrURL = convertStringsToChars(fileOrURL);
end
[filename, isURL] = internal.map.checkfilename( ...
    fileOrURL, {'tif', 'tiff'}, mfilename, 1, true);

if (isURL)
    % Delete temporary file from Internet download regardless if any
    % errors are issued.
    clean = onCleanup(@() deleteDownload(filename));
end

% Create the information structure from contents in the file.
info = readinfo(filename);

if (isURL)
    % Replace the temporary file name with the URL.
    info.Filename = fileOrURL;
end

%------------------------- info functions ---------------------------------

function info = readinfo(filename)
% Read the GeoTIFF info from the file and return the information structure.

% Obtain the TIFF information structure from the file.
tiff_info = tiffinfo(filename);
tiff_info = tiff_info(1); % force to 1st one (edited by Shi)
% Validate that at least one GeoTIFF tag exists in the file.
map.internal.assert(~isempty(fieldnames(tiff_info(1).GeoTIFFTags)), ...
    'map:geotiff:expectedGeoTIFFTags', filename, 'IMFINFO')
 
% Ensure that all images have the same Height and Width values.
map.internal.assert(isconsistent(tiff_info.Height) && isconsistent(tiff_info.Width), ...
    'map:geotiff:inconsistentImageSizes');

% Obtain the GeoTIFF tags from the TIFF info structure.
tags = tiff_info(1).GeoTIFFTags;
tags = tags(1); % force to 1st (edit by Shi)
% geokey = 'UTM Zone 16 N with WGS84';
% tags.GeoKeyDirectoryTag.GTCitationGeoKey = geokey;
% tags.GeoAsciiParamsTag = [geokey,'|'];

% Obtain height and width variables from the first element.
height = tiff_info(1).Height;
width  = tiff_info(1).Width;

% Obtain the GeoTIFF information structure from the TIFF tags.
gtiff_info = gtiffinfo(tags, height, width);

% Construct a new information structure and copy the input data to it.
info = copyinfo(tiff_info, gtiff_info);

%--------------------------------------------------------------------------

function tiff_info = tiffinfo(filename)
% Obtain subset of TIFF information from filename.

% Turn off TIFF library warnings for this function.
w = warning('off', 'imageio:tiffmexutils:libtiffWarning');
wobj = onCleanup(@()warning(w));

% Open and validate the file.
% Close the TIFF file when the function terminates.
t = Tiff(filename);
tobj = onCleanup(@()close(t));
  
% Obtain the tag names to retrieve from the file.
[tiffTags, optionalTags, geoTiffTags] = getTagNames();

% Obtain the TIFF tag values from the Tiff object and return a structure
% with field names that match the names of the requested tags.
tags = getTiffTags(t, tiffTags);

% Create the TIFF information structure.
tiff_info = tiffstruct(filename, tags);
tiff_info = tiff_info(1); % force to 1st (edited by Shi)

% Add the TIFF tags to the structure.
tiff_info = addTiffTags(tiff_info, tags, optionalTags);

% Add the GeoTIFF tags to the structure.
tiff_info = addGeoTiffTags(tiff_info, tags, geoTiffTags);

%--------------------------------------------------------------------------

function gtiff_info = gtiffinfo(tags, height, width)
% Obtain GeoTIFF information from TIFF GeoTIFF tags.

% Create the GeoTIFF information structure.
gtiff_info = gtiffstruct();

% Obtain the model type from the GeoKeyDirectoryTag, if defined.
gtiff_info = modelinfo(gtiff_info, tags);

% Set UserDefined value.
gtiff_info.UserDefined = 32767;

switch gtiff_info.ModelType
    case 'ModelTypeProjected'
        % Obtain the PCS information from the tags.
        gtiff_info = pcsinfo(gtiff_info, tags);
        
    case 'ModelTypeGeographic'
        % Obtain the GCS information from the tags.
        gtiff_info = gcsinfo(gtiff_info, tags);
        
    case 'ModelTypeGeocentric'
        % Obtain any GCS and PCS information from the tags.
        gtiff_info = gcsinfo(gtiff_info, tags);
        gtiff_info = pcsinfo(gtiff_info, tags);
        
    otherwise
end

% Obtain the spatial referencing information.
gtiff_info = spatialinfo(gtiff_info, tags, height, width);

%--------------------------------------------------------------------------

function gtiff_info = modelinfo(gtiff_info, tags)
% Construct the ModelType string and Model code values from the values in
% the GeoKeyDirectoryTag. If the GTModelTypeTags is contained in the
% GeoTIFFTags structure use that value for the Model code; otherwise,
% set it to []. Use the Model code value to determine the ModelType string.

tagname = 'GeoKeyDirectoryTag';
keyname = 'GTModelTypeGeoKey';

if isfield(tags, tagname) && isfield(tags.(tagname), keyname)
    modelTypeValue = tags.(tagname).(keyname);
    gtiff_info.GeoTIFFCodes.Model = modelTypeValue;
else
    modelTypeValue = '';
    gtiff_info.GeoTIFFCodes.Model = [];
end

switch modelTypeValue
    case 1
        gtiff_info.ModelType = 'ModelTypeProjected';
        
    case 2
        gtiff_info.ModelType = 'ModelTypeGeographic';
        
    case 3
        gtiff_info.ModelType = 'ModelTypeGeocentric';
        
    otherwise
        gtiff_info.ModelType = '';
end

%--------------------------------------------------------------------------

function gtiff_info = pcsinfo(gtiff_info, tags)
% Get the PCS information from the values in the GeoKeyDirectoryTag,
% if it is contained in the TAGS structure. Otherwise, do not modify the
% GTIFF_INFO structure.

tagname = 'GeoKeyDirectoryTag';
keyname = 'ProjectedCSTypeGeoKey';
if isfield(tags, tagname) && isfield(tags.(tagname), keyname)
    
    % Obtain the geokeys from the GeoKeyDirectoryTag.
    geokeys = tags.(tagname);
    
    % Obtain the code  number.
    pcscode = geokeys.(keyname);
    
    if pcscode > 0      
        % Search the projected_crs table for this code.
        result = pcsread(pcscode, gtiff_info.UserDefined);
        
        % Assign the PCS code.
        gtiff_info.GeoTIFFCodes.PCS = pcscode;
        
        % Assign the PCS string.
        gtiff_info.PCS = getStringValue(result, 'name');
        
        % Assign the projcode.
        epsgname = 'conversion_code';
        keyname = 'ProjectionGeoKey';
        if isfield(geokeys, keyname) && ~isempty(geokeys.(keyname))
            default = geokeys.(keyname);
        else
            default = gtiff_info.UserDefined;
        end
        projcode = getDoubleValue(result, epsgname, default);
        gtiff_info.GeoTIFFCodes.ProjCode = projcode;
        
        % Obtain GCS code.
        header_name = 'geodetic_crs_code';
        code = getDoubleValue(result, header_name, gtiff_info.UserDefined);
        
        % Search the geodetic_crs table for this code.
        % If the code is found in the geokeys, then override the projection
        % GCS code.
        gtiff_info = gcsinfo(gtiff_info, tags, code);  

        % Assign the projection values.
        gtiff_info = projinfo(gtiff_info, projcode, geokeys);
        
        % Assign the length unit values.
        keyname = 'ProjLinearUnitsGeoKey';
        axis_result = getFirstAxis(result);
        code = getTagDoubleValue(geokeys, keyname, axis_result, 'uom_code');
        gtiff_info = uomlengthinfo(gtiff_info, code); 
        
    end
end

%--------------------------------------------------------------------------

function gtiff_info = projinfo(gtiff_info, code, geokeys)
% Convert the code to projection values.

% Assign Projection values.
epsgname = 'method_code';
result = tableread('conversion_table', code, epsgname, gtiff_info.UserDefined);
projmethod = getDoubleValue(result, epsgname, gtiff_info.UserDefined);
projname = getStringValue(result, 'name');
gtiff_info.Projection = projname;
gtiff_info.GeoTIFFCodes.Projection = projmethod;

% Assign CTProjection values.
keyname = 'ProjCoordTransGeoKey';
if isfield(geokeys, keyname)
    % Code is found in the geokey. Convert the standard CT code number to a
    % CTProjection name.
    projcode = geokeys.(keyname);
    gtiff_info.CTProjection = codeToCTProjection(projcode);
else
    % Code is found in EPSG table. Convert the code value to a CTProjection
    % name. Using the inverse map, lookup the CTProjection code from the
    % CTProjection name.
    if code ~= gtiff_info.UserDefined
        projcode = projmethod;
        [projname, inverseMap] = codeToCTProjection(projcode);
        projcode = inverseMap(projname);
        gtiff_info.CTProjection = projname;
    else
        % code is user-defined.
        projcode = code;
        gtiff_info.CTProjection = '';
    end
end
gtiff_info.GeoTIFFCodes.CTProjection = projcode;

% Assign UTM, State Plane, and Zone EPSG values.
UTM_North = -9001;
UTM_South = -9002;
State_Plane_27 = -9003;
State_Plane_83 = -9004;
UTM_zone_1N  = 16001;
UTM_zone_60N = 16060;
UTM_zone_1S  = 16101;
UTM_zone_60S = 16160;

if code >= UTM_zone_1N && code  <= UTM_zone_60N
    mapSysCode = UTM_North;
    mapSys = 'UTM_NORTH';
    zone   = code - UTM_zone_1N + 1;
    
elseif code >= UTM_zone_1S && code <= UTM_zone_60S
    mapSysCode = UTM_South;
    mapSys = 'UTM_SOUTH';
    zone   = code - UTM_zone_1S + 1;
    
elseif  code >= 10101 && code <= 15299
    if  rem(code, 100) >= 30 
        mapSysCode = State_Plane_83;
        mapSys = 'STATE_PLANE_83';
        zone   = code - 10000 - 30;
    else
        mapSysCode = State_Plane_27;
        mapSys = 'STATE_PLANE_27';
        zone   = code - 10000;
    end
else
    mapSysCode = gtiff_info.UserDefined;
    zone = gtiff_info.UserDefined;
    mapSys = '';
end

% Assign the MapSys and Zone values.
gtiff_info.GeoTIFFCodes.MapSys = mapSysCode;
gtiff_info.MapSys = mapSys;
gtiff_info.Zone = zone;

% Assign ProjParm info values.
gtiff_info = projparamsinfo(gtiff_info, result, geokeys);

%--------------------------------------------------------------------------

function gtiff_info = projparamsinfo(gtiff_info, result, geokeys)
% Assign the ProjParm info parameters.

[idToNameMap, nameToIdMap] = constructGeoKeyDirectoryMap;
projcode = gtiff_info.GeoTIFFCodes.CTProjection;
[name, inverseMap] = codeToCTProjection(projcode);
id = @(x)(nameToIdMap(lower(x)));

% Construct projParmCodes and EPSG codes values.
[projParmCodes, epsgCodes, EPSG] = constructProjParmCodes(name, id);

% Construct ProjParm values.
projParm = constructProjParm(result, EPSG, epsgCodes);

% Update proj params if the ProjCoordTransGeoKey has been set.
if isfield(geokeys, 'ProjCoordTransGeoKey')
    [projParm, projParmCodes] = updateProjParm( ...
        projParm, projParmCodes, geokeys, name, id);
end

% Translate codes to names.
projParmId = cell(size(projParmCodes));
idToNameMap(0) = 'Unknown-0';
for k=1:numel(projParmCodes)
    projParmId{k} = idToNameMap(projParmCodes(k));
end

% Assign values to gtiff_info.
gtiff_info.GeoTIFFCodes.ProjParmId = projParmCodes;
gtiff_info.ProjParmId = projParmId;
gtiff_info.ProjParm = projParm;

% Adjust projection parameters if using CT_PolarStereographic_B.
if strcmp(name, 'CT_PolarStereographic_B')
    gtiff_info = adjustPolarStereographicBParams(gtiff_info, inverseMap);
end

%--------------------------------------------------------------------------

function gtiff_info = gcsinfo(gtiff_info, tags, code)
% Get the GCS string from the value in the GeoKeyDirectoryTag,
% if it is contained in the GeoTIFFTags structure, otherwise, return empty.

if ~exist('code', 'var')
    code = 0;
end

if isfield(tags,'GeoKeyDirectoryTag')
    % Obtain the geokeys from the GeoKeyDirectoryTag.
    geokeys = tags.GeoKeyDirectoryTag;
    
    if isfield(geokeys, 'GeographicTypeGeoKey')    
        % Obtain the code number for the key.
        code = geokeys.GeographicTypeGeoKey;
    end
else
    geokeys =struct();
end

if code > 0
    % Search the geodetic_crs table for this code.
    result = gcsread(code, gtiff_info.UserDefined);
    
    % Assign the GCS code.
    gtiff_info.GeoTIFFCodes.GCS = code;
    
    % Assign the GCS string.
    gtiff_info.GCS = getStringValue(result, 'name');
    
    % Assign the angle unit values.
    % Note: The GeoTIFF library uses the angle code from the datum if
    % defined rather than the code from the key, if given. For
    % compatibility, use the code from the datum information if provided.
    % In the future, this behavior should change. Remove the four lines
    % below and replace them with this one.
    % code = getTagDoubleValue(geokeys, keyname, result, 'uom_code');
    keyname = 'GeogAngularUnitsGeoKey';
    axis_result = getFirstAxis(result);
    uom_code = getDoubleValue(axis_result, 'uom_code', gtiff_info.UserDefined);
    if uom_code == gtiff_info.UserDefined
        uom_code = getTagDoubleValue(geokeys, keyname, axis_result, 'uom_code');
    end
    gtiff_info = uomangleinfo(gtiff_info, uom_code);

    % Assign the Datum code.
    keyname = 'GeogGeodeticDatumGeoKey';
    code = getTagDoubleValue(geokeys, keyname, result, 'datum_code');
    gtiff_info.GeoTIFFCodes.Datum = code;
    
    % Search the geodetic_datum table for this code.
    epsgname = 'name';
    result = tableread('geodetic_datum', code, epsgname, gtiff_info.UserDefined);
    
    % Assign the Datum string.
    gtiff_info.Datum = getStringValue(result, epsgname);
    
    % Assign the PM code.
    keyname = 'GeogPrimeMeridianGeoKey';
    code = getTagDoubleValue(geokeys, keyname, result, 'prime_meridian_code');
    gtiff_info.GeoTIFFCodes.PM = code;
    
    % Assign the PM values.
    gtiff_info = pminfo(gtiff_info, code);
    
    % Assign the Ellipsoid code.
    keyname = 'GeogEllipsoidGeoKey';
    code = getTagDoubleValue(geokeys, keyname, result, 'ellipsoid_code');
    gtiff_info.GeoTIFFCodes.Ellipsoid = code;
end
    
% Assign the ellipsoid values.
gtiff_info = ellipsoidinfo(gtiff_info, geokeys, code);

%--------------------------------------------------------------------------

function gtiff_info = pminfo(gtiff_info, code)
% Search the prime_meridian table for the code.

% Check to see if the prime meridian code is for Greenwich since it is the
% most common code used.
PM_Greenwich = 8901;
if code == PM_Greenwich
    gtiff_info.PM = 'Greenwich';
    gtiff_info.PMLongToGreenwich = 0;
else
    epsgname = 'name';
    result = tableread('prime_meridian', code, epsgname, gtiff_info.UserDefined);
    gtiff_info.PM = getStringValue(result, epsgname);
    
    epsgname = 'longitude';
    gtiff_info.PMLongToGreenwich = getDoubleValue(result, epsgname, 0);
end
gtiff_info.GeoTIFFCodes.PM = code;

%--------------------------------------------------------------------------

function gtiff_info = uomlengthinfo(gtiff_info, code)
% Search the unit_of_measure table for the code.

% Save time with common values.
switch code
    case 9001
        gtiff_info.UOMLength = 'metre';
        gtiff_info.UOMLengthInMeters = 1.0;
        
    case 9002
        gtiff_info.UOMLength = 'foot';
        gtiff_info.UOMLengthInMeters = .3048;
        
    case 9003
        gtiff_info.UOMLength = 'US survey foot';
        gtiff_info.UOMLengthInMeters = 12.0/39.37;
        
    otherwise
        [name, conv_factor] = uomread(code, gtiff_info.UserDefined);
        gtiff_info.UOMLength = name;
        if isempty(conv_factor)
            dfInMeters = 0.0;
        else
            dfInMeters = conv_factor; 
        end
        gtiff_info.UOMLengthInMeters = dfInMeters;    
end

% Assign the UOMLength code.
gtiff_info.GeoTIFFCodes.UOMLength = code;

%--------------------------------------------------------------------------

function gtiff_info = uomangleinfo(gtiff_info, code)
% Search the unit_of_measure table for the code.

% Save time with common values.
switch code
    case 9101
         gtiff_info.UOMAngle = 'radian';
        gtiff_info.UOMAngleInDegrees = 180.0 / pi;
       
    case {9102 9107 9108 9110 9122}
        gtiff_info.UOMAngle = 'degree';
        gtiff_info.UOMAngleInDegrees = 1.0;
        
    case 9103
        gtiff_info.UOMAngle = 'arc-minute';
        gtiff_info.UOMAngleInDegrees = 1.0 / 60.0;
        
    case 9104
        gtiff_info.UOMAngle = 'arc-second';
        gtiff_info.UOMAngleInDegrees = 1.0 / 3600.0;
        
   case 9105
        gtiff_info.UOMAngle = 'grad';
        gtiff_info.UOMAngleInDegrees = 180.0 / 200.0;
        
   case 9106
        gtiff_info.UOMAngle = 'gon';
        gtiff_info.UOMAngleInDegrees = 180.0 / 200.0;
        
   case 9109
        gtiff_info.UOMAngle = 'microradian';
        gtiff_info.UOMAngleInDegrees = 180.0 / (pi *  1000000.0);
 
    otherwise
        [name, conv_factor] = uomread(code, gtiff_info.UserDefined);
        if ~isempty(conv_factor) && conv_factor ~= 0.0
            angleInRadians = conv_factor;
            angleInDegrees = angleInRadians * 180.0 / pi;
        else
            angleInDegrees = [];
        end
        gtiff_info.UOMAngle = name;
        gtiff_info.UOMAngleInDegrees = angleInDegrees;        
end

% Assign the UOMAngle code.
gtiff_info.GeoTIFFCodes.UOMAngle = code;

%--------------------------------------------------------------------------

function gtiff_info = ellipsoidinfo(gtiff_info, geokeys, code)
% Get the ellipsoid information using the code value. If the code equals
% UserDefined, then set the values to empty. Override SemiMajor and
% SemiMinor values if found in the geokeys.
 
if code == gtiff_info.UserDefined || isempty(code) || code <= 0
    gtiff_info.Ellipsoid = '';  % 'UserDefined'
    gtiff_info.SemiMajor = [];
    gtiff_info.SemiMinor = [];
    
elseif code > 0
    % Obtain the ellipsoid values from the table.
    e = epsgread('ellipsoid', code);
    
    % Assign the Ellipsoid name, SemiMajor, and SemiMinor values.
    gtiff_info.Ellipsoid = e.Ellipsoid;
    gtiff_info.SemiMajor = e.SemiMajor;
    gtiff_info.SemiMinor = e.SemiMinor;
end    

% Check for overridden SemiMajor and SemiMinor parameters. 
keyname = 'GeogSemiMajorAxisGeoKey';     
if isfield(geokeys, keyname)
    gtiff_info.SemiMajor = geokeys.(keyname);
end

keyname = 'GeogSemiMinorAxisGeoKey';
if isfield(geokeys, keyname)
    gtiff_info.SemiMinor = geokeys.(keyname);
end

%--------------------------------------------------------------------------

function gtiff_info = spatialinfo(gtiff_info, tags, height, width)
% Obtain the spatial referencing information.

% Assign the PixelScale value.
if isfield(tags, 'ModelPixelScaleTag') 
    gtiff_info.PixelScale = tags.ModelPixelScaleTag';
else
    gtiff_info.PixelScale = [];
end

% Construct the spatial referencing object.
gtiff_info.SpatialRef = constructSpatialRef( ...
    tags, height, width, gtiff_info.ModelType);

% Construct the referencing matrix.
gtiff_info.RefMatrix = constructRefMatrix(gtiff_info.SpatialRef);

% Construct the TiePoints structure.
gtiff_info.TiePoints = constructTiePoints(tags);

% Construct the BoundingBox array.
gtiff_info.BoundingBox = constructBoundingBox(gtiff_info.SpatialRef);

% Construct the CornerCoords structure.
gtiff_info = constructCornerCoords(gtiff_info, gtiff_info.SpatialRef);

%------------------------- struct functions -------------------------------

function tiff_info = tiffstruct(filename, tags)
% Create the TIFF information structure. tags is a structure array with one
% element per TIFF directory. optionalTags is a cell array of strings.

% Create the information structure.
requiredFields = { ...
    'Filename', ...
    'FileModDate', ...
    'FileSize', ...
    'Format', ...
    'FormatVersion', ...
    'Height', ...
    'Width', ...
    'BitDepth', ...
    'ColorType'};

tiff_info = cell2struct(cell(size(requiredFields)), requiredFields, 2);
numdirs = numel(tags);
tiff_info(numdirs) = tiff_info;

% Obtain information about the file.
d = dir(filename);
fileModDate = d.date;
fileSize = d.bytes;
fileFormat = 'tif';

% Copy data to the tiff_info structure.
[tiff_info.Filename]      = deal(filename);
[tiff_info.FileModDate]   = deal(fileModDate);
[tiff_info.FileSize]      = deal(fileSize);
[tiff_info.Format]        = deal(fileFormat);
[tiff_info.FormatVersion] = deal([]);
[tiff_info.Height]        = deal(tags.ImageLength);
[tiff_info.Width]         = deal(tags.ImageWidth);

%--------------------------------------------------------------------------

function gtiff_info = gtiffstruct()
% Create the GeoTIFF information structure.

codes = struct( ...
    'Model', [], ...
    'PCS', [], ...
    'GCS', [], ...
    'UOMLength', [], ...
    'UOMAngle', [], ...
    'Datum', [], ...
    'PM', [], ...
    'Ellipsoid', [], ...
    'ProjCode', [], ...
    'Projection', [], ...
    'CTProjection', [], ...
    'MapSys', [], ...
    'ProjParmId', []);

gtiff_info = struct( ...
    'ModelType', '', ...
    'PCS', '', ...
    'GCS', '', ...
    'UOMLength', '', ...
    'UOMLengthInMeters', [], ...
    'UOMAngle', '', ...
    'UOMAngleInDegrees', [], ...
    'Datum', '', ...
    'PM', '', ...
    'PMLongToGreenwich', [], ...
    'Ellipsoid', '', ...
    'SemiMajor', [], ...
    'SemiMinor', [], ...
    'Projection', '', ...
    'CTProjection', '', ...
    'ProjParm', [], ...
    'ProjParmId', '', ...
    'PixelScale', [], ...
    'MapSys', '', ...
    'Zone', [], ...
    'CornerCoords', [], ...
    'GeoTIFFCodes', codes);

%--------------------------------------------------------------------------

function info = infostruct(tiff_info)
% Initialize info structure and set it's size to the size of tiff_info.

info = struct( ...
    'Filename', '', ...
    'FileModDate', '', ...
    'FileSize', [], ...
    'Format', '', ...
    'FormatVersion', [], ...
    'Height', [], ...
    'Width', [], ...
    'BitDepth', [], ...
    'ColorType', '', ...
    'ModelType', '', ...
    'PCS', '', ...
    'Projection', '', ...
    'MapSys', '', ...
    'Zone', [], ...
    'CTProjection', '', ...
    'ProjParm', [], ...
    'ProjParmId', [], ...
    'GCS', '', ...
    'Datum', '', ...
    'Ellipsoid', '', ...
    'SemiMajor', [], ...
    'SemiMinor', [], ...
    'PM', '', ...
    'PMLongToGreenwich', [], ...
    'UOMLength', '', ...
    'UOMLengthInMeters', [], ...
    'UOMAngle', '', ...
    'UOMAngleInDegrees', [], ...
    'TiePoints', [], ...
    'PixelScale', [], ...
    'SpatialRef', [], ...
    'RefMatrix',  [], ...
    'BoundingBox', [], ...
    'CornerCoords', [], ...
    'GeoTIFFCodes', [], ...
    'GeoTIFFTags', []);

if isfield(tiff_info,'ImageDescription')
   info.ImageDescription = '';
end
    
info(numel(tiff_info)) = info;

%--------------------------------------------------------------------------

function EPSG = epsgstruct()
% Construct EPSG structure.

EPSG = struct( ...
    'NatOriginLat'         , 8801, ...
    'NatOriginLong'        , 8802, ...
    'NatOriginScaleFactor' , 8805, ...
    'FalseEasting'         , 8806, ...
    'FalseNorthing'        , 8807, ...
    'ProjCenterLat'        , 8811, ...
    'ProjCenterLong'       , 8812, ...
    'Azimuth'              , 8813, ...
    'AngleRectifiedToSkewedGrid' , 8814, ...
    'InitialLineScaleFactor', 8815, ...
    'ProjCenterEasting'    , 8816, ...
    'ProjCenterNorthing'   , 8817, ...
    'PseudoStdParallelScaleFactor' , 8819, ...
    'FalseOriginLat'       , 8821, ...
    'FalseOriginLong'      , 8822, ...
    'StdParallel1Lat'      , 8823, ...
    'StdParallel2Lat'      , 8824, ...
    'FalseOriginEasting'   , 8826, ...
    'FalseOriginNorthing'  , 8827, ...
    'StdParallelLat'       , 8832, ...
    'OriginLong'           , 8833);

%------------------------- read functions ---------------------------------

function result = epsgread(table, code)
% Search the EPSG table for the code value. If found, return a structure,
% otherwise return empty.

% For performance reasons, cache the found results into the persistent
% variable epsgtables. Only store the found results since the entire
% database is over 80 MB when converted to a structure.
persistent epsgtables
if isempty(epsgtables) || ~isstruct(epsgtables)
    epsgtables = struct();
end    
  
% Remove  any '.' from table name.
tablename = strrep(table, '.', '_');

% Assign a logical if reading the ellipsoid table.
readingEllipsoid = strcmp('ellipsoid', tablename);

if readingEllipsoid
    if ~isfield(epsgtables, tablename)
        % The ellipsoid table has never been read. Read the table and cache
        % the found results into the epsgtables structure.
        result = ellipsoid2struct(code);
        
        % Save the ellipsoid structure.
        epsgtables.(tablename) = result;
    else
        % The ellipsoid table has been previously read. Determine if the
        % code is in the cached table. If not, then read the table and
        % append the results.
        S = epsgtables.(tablename);
        result = findtable(S, code);
        if isempty(result)
            result = ellipsoid2struct(code);
            epsgtables.(tablename)(end+1) = result;
        end
    end
elseif ~isfield(epsgtables, tablename)
    % The table has never been read. Read the table and cache the found
    % results into the epsgtables structure.
    result = epsg2struct(table, code);
    if ~isempty(result)
        epsgtables.(tablename) = result;
    end
else
    % The table has been previously read. Determine if the code is in the
    % cached table. If not, then read the table and append the found
    % result.
    S = epsgtables.(tablename);
    result = findtable(S, code);
    if isempty(result)
        S = epsg2struct(table, code);
        result = findtable(S, code);
        if ~isempty(result)
            epsgtables.(tablename)(end+1) = result;
        end
    end   
end

%--------------------------------------------------------------------------

function S = epsg2struct(table, varargin)
% Read the EPSG table and return a structure.

% Read the table.
r = map.internal.epsgread(table, varargin{:});

if isempty(varargin)
    % Output is a cell array. Convert to a structure.
    % Preserve the header.
    header = r(1,:);
    r(1,:) = [];
    
    % Find any empty header names and remove them.
    index = cellfun(@isempty, header);
    header(index) = [];
    r(:,index) = [];
    
    % Convert cell array to structure with the header names as field names.
    S = cell2struct(r,header,2);
else
    S = r;
end

%--------------------------------------------------------------------------

function S = ellipsoid2struct(code)
% Search the EPSG ellipsoid table for the code and return a structure
% containing the code value (as a string), the name, semimajor and
% semiminor values.

strCode = sprintf('%d', code);
try
    % Create a referenceEllipsoid object using the code value.
    e = referenceEllipsoid(code);
    
    % Translate into meters for the GeoTIFF specification.
    e.LengthUnit = 'meter';
    
    % Assign the code as a string value in order to locate it and the
    % Ellipsoid name, SemiMajor, and SemiMinor values.
    S.code = strCode;
    S.Ellipsoid = e.Name;
    S.SemiMajor = e.SemimajorAxis;
    S.SemiMinor = e.SemiminorAxis;    
catch e %#ok<NASGU>
    S.code = strCode;
    S.Ellipsoid = 'Unknown';
    S.SemiMajor = [];
    S.SemiMinor = [];
end

%--------------------------------------------------------------------------

function result = findtable(S, code)
% Given a table structure array, S, search for the code string. If found,
% then return the result (which is a scalar structure). Otherwise, return
% [].

if ~isstruct(S)
    result = [];
else
    if isnumeric(code)
        code = sprintf('%d', code);
    end
    values = {S.code};
    index = strcmp(code, values);
    if any(index)
        result = S(index);
    else
        result = [];
    end
end

%--------------------------------------------------------------------------

function result = tableread(table, code, epsgname, userDefined)
% Search the table for the code.

if code ~= userDefined
    result = epsgread(table, code);
else
    result.(epsgname) = '';   % 'UserDefined';
end

%--------------------------------------------------------------------------

function result = pcsread(code, userDefined)
% Search the projected_crs table for the code.

if code ~= userDefined
    result = epsgread('projected_crs', code);
else
    % result.name = 'UserDefined';
    result.name = '';
end

%--------------------------------------------------------------------------

function result = gcsread(code, userDefined)
% Search the geodetic_crs table for the code.

if code ~= userDefined
    result = epsgread('geodetic_crs', code);
else
    % result.name = 'UserDefined';
    result.name = '';
end

%--------------------------------------------------------------------------

function [name, conv_factor] = uomread(code, userDefined)
% Search the unit_of_measure table for the code. Return [] if the code is
% not found.

epsgname = 'name';
result = tableread('unit_of_measure', code, epsgname, userDefined);
name = getStringValue(result, epsgname);
conv_factor = getDoubleValue(result, 'conv_factor', []);

%--------------------------------------------------------------------------

function axis_result = getFirstAxis(gtif_result)
% Search the axis table for the UOM code.

% For performance reasons, cache the axis table into the persistent
% variable epsgAxisTable.
persistent epsgAxisTable
if isempty(epsgAxisTable) || ~isstruct(epsgAxisTable)
    epsgAxisTable = epsg2struct('axis');
end 

epsgname = 'coordinate_system_code';
if isfield(gtif_result, epsgname)
    S = epsgAxisTable;
    result = S(strcmp({S.(epsgname)},gtif_result.(epsgname)));
    if ~isempty(result)
        axis_result = result(1);
    else
        axis_result = result;
    end
else
    axis_result = [];
end

%------------------------- tag functions ----------------------------------

function [tiffTags, optionalTags, geoTiffTags] = getTagNames()
% Get tag names to retrieve from a file.

% Assign the required and optional tag names.
requiredTags = { ...
    'ImageLength', ...
    'ImageWidth', ...
    'BitsPerSample', ...
    'ColorMap', ...
    'Photometric'};

optionalTags = {'ImageDescription'};

% Obtain a row cell vector of GeoTIFF tag names obtained from the Tiff
% object.
tiffNames = Tiff.getTagNames;
cindex1 = regexp(tiffNames, regexptranslate('wildcard', 'Model*Tag'));
cindex2 = regexp(tiffNames, regexptranslate('wildcard', 'Geo*Tag'));
index1 = ~cellfun(@isempty, cindex1); 
index2 = ~cellfun(@isempty, cindex2); 
geoTiffTags = tiffNames(index1 | index2)';

% Add RPCCoefficientTag.
RPCCoefficientTagName = 'RPCCoefficientTag';
geoTiffTags = [geoTiffTags RPCCoefficientTagName];

% Assign all the tags to tiffTags.
tiffTags = [requiredTags, optionalTags, geoTiffTags];

%--------------------------------------------------------------------------

function tags = getTiffTags(t, tiffTagNames)
% Obtain TIFF tag values from Tiff object t, for the TIFF tags specified
% by the tags names in the cell array tiffTagNames.

% Obtain the TIFF tag values using the Tiff object, t.
tags = getTags(t, tiffTagNames);

% If the RPCCoefficientTag is not supported by the Tiff class, then obtain
% the tag value using imfinfo.
tagName = 'RPCCoefficientTag';
if ~any(strcmp(tagName, Tiff.getTagNames))   
    tags = getRPCCoefficientTag(t.FileName, tags, tagName);
end

%--------------------------------------------------------------------------

function tags = getTags(t, tiffTagNames)
% Obtain TIFF tag values from Tiff object, t.

% Obtain the number of TIFF directories.
numdirs = numberOfDirectories(t);

% Create a structure with the TIFF tags as field names.
tags = cell2struct(cell(size(tiffTagNames)), tiffTagNames, 2);
tags(numdirs) = tags(1);

% Obtain TIFF tag values from file.
for dirnum = 1:numdirs
    t.setDirectory(dirnum);
    for k=1:numel(tiffTagNames)
        tagName = tiffTagNames{k};
        warnObj = turnOffAllWarnings;
        if hasTag(t, tagName)
            value = t.getTag(tagName);
            tags(dirnum).(tagName) = value;
        end
        delete(warnObj)
    end
end

%--------------------------------------------------------------------------

function tf = hasTag(t, tagName)
% Return true if the TIFF file has the tag specified by tagName. Ignore any
% warnings. Return false if an error is issued.
try
    t.getTag(tagName);
    tf = true;
catch 
    tf = false;
end

    
%--------------------------------------------------------------------------

function tags = getRPCCoefficientTag(filename, tags, tagName)
% If present in the file specified by FILENAME, obtain the RPCCoefficient
% tag value. If the tag is found in the file, set its value in the TAGS
% structure using the field name specified by tagName. Otherwise the TAGS
% structure is unchanged.

tagValue = map.geotiff.internal.readRPCCoefficientTag(filename);
index = ~cellfun(@isempty,tagValue);
if any(index)
    [tags(index).(tagName)] = tagValue{index};
end

%--------------------------------------------------------------------------

function warnObj = turnOffAllWarnings
% Turn off all warnings until warnObj is deleted.

warnState = warning('query', 'all');
warnObj = onCleanup(@() warning(warnState));
warning('off', 'all')

%--------------------------------------------------------------------------

function tiff_info = addTiffTags(tiff_info, tags, optionalTags)
% Add the TIFF tags to the information structure.

% Process BitDepth
for k=1:numel(tags)
    tiff_info(k).BitDepth = sum(tags(k).BitsPerSample);
end

% Process ColorType field.
for k = 1:numel(tags)
    if ~isempty(tags(k).ColorMap)
        tiff_info(k).ColorType = 'indexed';
    elseif tags(k).Photometric == Tiff.Photometric.MinIsWhite ...
            || tags(k).Photometric == Tiff.Photometric.MinIsBlack
        tiff_info(k).ColorType = 'grayscale';
    elseif tags(k).Photometric == Tiff.Photometric.RGB
        tiff_info(k).ColorType = 'truecolor';
    else
        tiff_info(k).ColorType = 'unknown';
    end
end

% Add the optionalTags if they are not empty.
for k=1:numel(optionalTags)
    tagname = optionalTags{k};
    if any(~cellfun(@isempty, {tags.(tagname)}))
        [tiff_info.(tagname)] = deal(tags.(tagname));
    end
end

%--------------------------------------------------------------------------

function tiff_info = addGeoTiffTags(tiff_info, tags, geoTiffTagNames)
% Add 'GeoTIFFTags' field to the info structure, tiff_info. Assume that if
% the file contains multiple images, they have the same GeoTIFF tag
% information (consistent image size is enforced later). Translate the Tiff
% tag name, 'GeoASCIIParamsTag' to 'GeoAsciiParamsTag' if present, based on
% the name in the GeoTIFF specification. Translate the GeoKeyDirectoryTag
% to a scalar structure. Translate RPCCoefficientTag to a scalar
% map.geotiff.RPCCoefficientTag, if present.

% Construct a GeoTIFF tags structure.
geoTiffTags = constructGeoTIFFTagsStructure(tags, geoTiffTagNames);

% Translate the first GeoKeyDirectoryTag element to a structure. Copy the
% structure to the other elements since we are assuming that if the file
% contains multiple images, they have the same GeoTIFF tag information.
tagName = 'GeoKeyDirectoryTag';
if isfield(geoTiffTags, tagName)
    [geoTiffTags.(tagName)] = deal(translateGeoKeyDirectoryTag(geoTiffTags(1)));
end

% Translate the Tiff tag name, 'GeoASCIIParamsTag' to 'GeoAsciiParamsTag'
% if present. Ensure field names are in the same order.
GeoASCIIParamsTag = 'GeoASCIIParamsTag';
GeoAsciiParamsTag = 'GeoAsciiParamsTag';
if isfield(geoTiffTags, GeoASCIIParamsTag)
    oldOrder = fieldnames(geoTiffTags);
    newOrder = strrep(oldOrder, GeoASCIIParamsTag, GeoAsciiParamsTag);
    [geoTiffTags.(GeoAsciiParamsTag)] = geoTiffTags.(GeoASCIIParamsTag);
    geoTiffTags = rmfield(geoTiffTags, GeoASCIIParamsTag);
    geoTiffTags = orderfields(geoTiffTags, newOrder);
end

% Translate the RPCCoefficientTag elements to a scalar
% map.geotiff.RPCCoefficientTag, if present. The tag may be different for
% each image in the file.
tagName = 'RPCCoefficientTag';
if isfield(geoTiffTags, tagName)
    for k = 1:length(geoTiffTags)
        tag = geoTiffTags(k).(tagName);
        geoTiffTags(k).(tagName) = translateRPCCoefficientTag(tag);
    end
    if isempty([geoTiffTags.(tagName)])
        geoTiffTags = rmfield(geoTiffTags, tagName);
    end
end

% Assign the values to the output structure. The length of the two
% structures is identical. Use deal since tiff_info does not contain the
% GeoTIFFTags field.
[tiff_info.GeoTIFFTags] = deal(geoTiffTags);

%--------------------------------------------------------------------------

function geoTiffTags = constructGeoTIFFTagsStructure(tags, geoTiffTagNames)
% Construct a GeoTIFF tags structure from the input tags structure.

% Remove fields in tags that are not GeoTIFF tags.
nonGeoTiffTags = setdiff(fieldnames(tags), geoTiffTagNames);
geoTiffTags = rmfield(tags, nonGeoTiffTags);

% Remove empty fields. structfun only accepts scalar structure.
names = fieldnames(geoTiffTags);
index = structfun(@isempty, geoTiffTags(1));
geoTiffTags = rmfield(geoTiffTags, names(index));

%------------------------- get functions ----------------------------------

function rasterInterpretation = getRasterInterpretation(tags)
% Obtain the raster interpretation from the GeoTIFF tags. From the tags, it
% will be 1 for 'RasterPixelIsArea' and 2 for 'RasterPixelIsPoint'. Convert
% the value to a string that is recognized by the RasterReference classes.

% If GTRasterTypeGeoKey ('RasterPixelIsArea' or 'RasterPixelIsPoint') is
% not specified, then default the value to 'RasterPixelIsArea'.
tagname = 'GeoKeyDirectoryTag';
keyname = 'GTRasterTypeGeoKey';
if ~isfield(tags, tagname) ...
        || (isfield(tags, tagname) && ~isfield(tags.(tagname), keyname))
    rasterTypeGeoKey = 1; % RasterPixelIsArea
else
    rasterTypeGeoKey = tags.(tagname).(keyname);
end

% Convert the value to a string.
if rasterTypeGeoKey == 1
    % 'RasterPixelIsArea'
    rasterInterpretation = 'cells';
else
    % 'RasterPixelIsPoint'
    rasterInterpretation = 'postings';
end

%--------------------------------------------------------------------------

function value = getStringValue(S, fname)
% Obtain the string value from the field, FNAME, of the structure S. If
% not found, return the default value, ''.

if ~isempty(S) && isfield(S, fname)
    value = S.(fname);
else
    value = '';
end

%--------------------------------------------------------------------------

function value = getDoubleValue(S, fname, default)
% Obtain the numeric value from the field, FNAME, of the structure S. If
% not found, return the default value, DEFAULT.

if ~isempty(S) && isfield(S, fname) && ~isempty(S.(fname))
    if ~isnumeric(S.(fname))
        value = str2double(S.(fname));
    else
        value = S.(fname);
    end
else
    value = default;
end

%--------------------------------------------------------------------------

function value = getTagDoubleValue(geokeys, keyname, result, epsgname)
% Obtain the keyname value from the geokeys, if present. Otherwise, return
% the value from the result structure.

if isfield(geokeys, keyname)
    value = geokeys.(keyname);
else
    userDefined = 32767;
    value = getDoubleValue(result, epsgname, userDefined);
end

%--------------------------------------------------------------------------

function value = keyget(geokeys, keynames, names, default)
% Return a value from the GEOKEYS structure based on the elements in the
% cell array NAMES. If any element in NAMES matches any element in
% KEYNAMES, then return the first field that matches, if the value is
% numeric. Otherwise, return the DEFAULT value.

index = cellcmp(names, keynames);
if any(index)
    name = names(index);
    value = geokeys.(name{1});
    if ~isnumeric(value)
        value = default;
    end
else
    value = default;
end

%------------------------- construct functions ----------------------------

function [projParmId, anEPSGCodes, EPSG] = constructProjParmCodes(name, id)

projParmId = zeros(7, 1);
anEPSGCodes = zeros(7, 1);
EPSG = epsgstruct();

switch name
    
    case {'CT_CassiniSoldner' ,'CT_NewZealandMapGrid'}
        projParmId(1) = id('ProjNatOriginLatGeoKey');
        projParmId(2) = id('ProjNatOriginLongGeoKey');
        projParmId(6) = id('ProjFalseEastingGeoKey');
        projParmId(7) = id('ProjFalseNorthingGeoKey');
        
        anEPSGCodes(1) = EPSG.NatOriginLat;
        anEPSGCodes(2) = EPSG.NatOriginLong;
        anEPSGCodes(6) = EPSG.FalseEasting;
        anEPSGCodes(7) = EPSG.FalseNorthing;
    
    case 'CT_ObliqueMercator'
        projParmId(1) = id('ProjCenterLatGeoKey');
        projParmId(2) = id('ProjCenterLongGeoKey');
        projParmId(3) = id('ProjAzimuthAngleGeoKey');
        projParmId(4) = id('ProjRectifiedGridAngleGeoKey');
        projParmId(5) = id('ProjScaleAtCenterGeoKey');
        projParmId(6) = id('ProjFalseEastingGeoKey');
        projParmId(7) = id('ProjFalseNorthingGeoKey');
        
        anEPSGCodes(1) = EPSG.ProjCenterLat;
        anEPSGCodes(2) = EPSG.ProjCenterLong;
        anEPSGCodes(3) = EPSG.Azimuth;
        anEPSGCodes(4) = EPSG.AngleRectifiedToSkewedGrid;
        anEPSGCodes(5) = EPSG.InitialLineScaleFactor;
        anEPSGCodes(6) = EPSG.ProjCenterEasting;
        anEPSGCodes(7) = EPSG.ProjCenterNorthing;
            
    case 'CT_ObliqueMercator_Laborde'
        projParmId(1) = id('ProjCenterLatGeoKey');
        projParmId(2) = id('ProjCenterLongGeoKey');
        projParmId(3) = id('ProjAzimuthAngleGeoKey');
        projParmId(5) = id('ProjScaleAtCenterGeoKey');
        projParmId(6) = id('ProjFalseEastingGeoKey');
        projParmId(7) = id('ProjFalseNorthingGeoKey');
        
        anEPSGCodes(1) = EPSG.ProjCenterLat;
        anEPSGCodes(2) = EPSG.ProjCenterLong;
        anEPSGCodes(3) = EPSG.Azimuth;
        anEPSGCodes(5) = EPSG.InitialLineScaleFactor;
        anEPSGCodes(6) = EPSG.ProjCenterEasting;
        anEPSGCodes(7) = EPSG.ProjCenterNorthing;
        
    case { ...
            'CT_LambertConfConic_1SP',  ...
            'CT_Mercator', ...
            'CT_ObliqueStereographic', ...
            'CT_PolarStereographic', ...
            'CT_TransverseMercator', ...
            'CT_TransvMercator_SouthOriented'}
        projParmId(1) = id('ProjNatOriginLatGeoKey');
        projParmId(2) = id('ProjNatOriginLongGeoKey');
        projParmId(5) = id('ProjScaleAtNatOriginGeoKey');
        projParmId(6) = id('ProjFalseEastingGeoKey');
        projParmId(7) = id('ProjFalseNorthingGeoKey');
        
        anEPSGCodes(1) = EPSG.NatOriginLat;
        anEPSGCodes(2) = EPSG.NatOriginLong;
        anEPSGCodes(5) = EPSG.NatOriginScaleFactor;
        anEPSGCodes(6) = EPSG.FalseEasting;
        anEPSGCodes(7) = EPSG.FalseNorthing;

    case 'CT_PolarStereographic_B'
        projParmId(1) = id('ProjNatOriginLatGeoKey');
        projParmId(2) = id('ProjNatOriginLongGeoKey');
        projParmId(5) = id('ProjScaleAtNatOriginGeoKey');
        projParmId(6) = id('ProjFalseEastingGeoKey');
        projParmId(7) = id('ProjFalseNorthingGeoKey');
        
        anEPSGCodes(1) = EPSG.StdParallelLat;
        anEPSGCodes(2) = EPSG.OriginLong;
        anEPSGCodes(5) = EPSG.NatOriginScaleFactor;
        anEPSGCodes(6) = EPSG.FalseEasting;
        anEPSGCodes(7) = EPSG.FalseNorthing;        
    
    case 'CT_LambertConfConic_2SP'
        projParmId(1) = id('ProjFalseOriginLatGeoKey');
        projParmId(2) = id('ProjFalseOriginLongGeoKey');
        projParmId(3) = id('ProjStdParallel1GeoKey');
        projParmId(4) = id('ProjStdParallel2GeoKey');
        projParmId(6) = id('ProjFalseEastingGeoKey');
        projParmId(7) = id('ProjFalseNorthingGeoKey');     
        
        anEPSGCodes(1) = EPSG.FalseOriginLat;
        anEPSGCodes(2) = EPSG.FalseOriginLong;
        anEPSGCodes(3) = EPSG.StdParallel1Lat;
        anEPSGCodes(4) = EPSG.StdParallel2Lat;
        anEPSGCodes(6) = EPSG.FalseOriginEasting;
        anEPSGCodes(7) = EPSG.FalseOriginNorthing;
        
    case 'CT_AlbersEqualArea'
        projParmId(1) = id('ProjStdParallel1GeoKey');
        projParmId(2) = id('ProjStdParallel2GeoKey');
        projParmId(3) = id('ProjNatOriginLatGeoKey');
        projParmId(4) = id('ProjNatOriginLongGeoKey');
        projParmId(6) = id('ProjFalseEastingGeoKey');
        projParmId(7) = id('ProjFalseNorthingGeoKey');
        
        anEPSGCodes(1) = EPSG.StdParallel1Lat;
        anEPSGCodes(2) = EPSG.StdParallel2Lat;
        anEPSGCodes(3) = EPSG.FalseOriginLat;
        anEPSGCodes(4) = EPSG.FalseOriginLong;
        anEPSGCodes(5) = EPSG.FalseOriginEasting;
        anEPSGCodes(7) = EPSG.FalseOriginNorthing;
            
    case {'CT_SwissObliqueCylindrical', 'CT_ObliqueMercator_Rosenmund'}
        projParmId(1) = id('ProjCenterLatGeoKey');
        projParmId(2) = id('ProjCenterLongGeoKey');
        projParmId(6) = id('ProjFalseEastingGeoKey');
        projParmId(7) = id('ProjFalseNorthingGeoKey');
        
    otherwise
end

%--------------------------------------------------------------------------

function projParm = constructProjParm(result, EPSG, anEPSGCodes)
% Assign ProjParm values.

% Initialize values.
n = 7;
[parameter_code, parameter_value, parameter_uom] = ...
    initializeProjParameterArrays(n, result);

adfProjParms = zeros(n, 1);
for k=1:n
    nEPSGCode = anEPSGCodes(k);
    if nEPSGCode == EPSG.AngleRectifiedToSkewedGrid
        adfProjParms(k) = 90;
    elseif nEPSGCode == EPSG.NatOriginScaleFactor ...
            || nEPSGCode == EPSG.InitialLineScaleFactor ...
            || nEPSGCode == EPSG.PseudoStdParallelScaleFactor
        adfProjParms(k) = 1.0;
    else
        adfProjParms(k) = 0;
    end
    
    iEPSG = nEPSGCode == parameter_code;
    nUOM = parameter_uom(iEPSG);
    pszValue = parameter_value(iEPSG);
    
    if isscalar(nUOM)
        if nUOM >= 9100 && nUOM < 9200
            adfProjParms(k) = angleToDegrees( pszValue, nUOM );
        elseif nUOM > 9000 && nUOM < 9100 
            S = uomlengthinfo(struct(), nUOM);
            dfInMeters = S.UOMLengthInMeters;
            adfProjParms(k) = pszValue * dfInMeters;
        else
            adfProjParms(k) = pszValue;
        end
    end
end

% Assign the values.
projParm = adfProjParms;

%--------------------------------------------------------------------------
function [parameter_code, parameter_value, parameter_uom] = ...
    initializeProjParameterArrays(n, result)

parameter_code = zeros(n, 1);
parameter_value = NaN(n, 1);
parameter_uom = zeros(n, 1);
for k=1:n
    epsgname = sprintf('param%d_code', k);
    parameter_code(k) = getDoubleValue(result, epsgname, NaN);
    
    epsgname = sprintf('param%d_value', k);
    if isfield(result,epsgname) && ~isempty(result.(epsgname))
        parameter_value(k) = result.(epsgname);
    end
    
    epsgname = sprintf('param%d_uom_code', k);
    parameter_uom(k) = getDoubleValue(result, epsgname, NaN);
end

%--------------------------------------------------------------------------

function spatialRef = constructSpatialRef(tags, height, width, modelType)
% Construct the spatialRef object from TIFF tags.

useTiepoint = isfield(tags, 'ModelTiepointTag') ...
    && isfield(tags, 'ModelPixelScaleTag') ...
    && numel(tags.ModelPixelScaleTag) >= 2 ...
    && numel(tags.ModelTiepointTag) == 6;

useMatrix = isfield(tags, 'ModelTransformationTag') ...
    && numel(tags.ModelTransformationTag) >= 8;

fileHasSpatialData = (useTiepoint || useMatrix) ...
    && height ~= 0 && width ~= 0;

if fileHasSpatialData
    % The file contains valid spatial referencing data in the tags. Use the
    % values from the tags and the raster size to construct a spatial
    % referencing object.
    rasterSize = [height width];
    rasterInterpretation = getRasterInterpretation(tags);
    spatialRef = constructRasterReference( ...
        useMatrix, tags, rasterSize, rasterInterpretation, modelType);
else
    % The file does not have valid referencing data in the tags or the
    % raster is invalid.
    spatialRef = [];
end

%--------------------------------------------------------------------------

function spatialRef = constructRasterReference( ...
    useMatrix, tags, rasterSize, rasterInterpretation, modelType)
% Construct either a map or geographic raster referencing object based on
% the value of the string, modelType.

% Construct first corners and the Jacobian based on tag values.
if useMatrix
    transform = tags.ModelTransformationTag;
    [x, y, J] = constructFromTransformMatrix(transform);
else
    tiepoint  = tags.ModelTiepointTag;
    scale     = tags.ModelPixelScaleTag;
    [x, y, J] = constructFromSingleTiepoint(tiepoint, scale);
end
   
if isempty(modelType)
    warning(message('map:geotiff:undefinedGTModelTypeGeoKey'));
    modelType = 'ModelTypeProjected';
end

switch modelType
    case 'ModelTypeGeographic'
        isRectilinear = (J(1,2) == 0 && J(2,1) == 0);
        if isRectilinear
            % Set the latitude and longitude corner point from the values
            % in the ModelTiepointTag. Clamp the corner latitude to the
            % interval [-90 90]. The longitude value does not need to be
            % clamped.
            lon = x;
            lat = y;
            if  lat < -90 || lat > 90
                warning(message('map:geotiff:clampingLatitudeCorner', ...
                    sprintf('%.15g', lat)));
                lat = min(90,  lat);
                lat = max(-90, lat);
            end
           
            deltaLatNumerator = J(2,2);
            deltaLonNumerator = J(1,1);
            try
                spatialRef = map.rasterref.internal.constructGeographicRasterReference( ...
                    rasterSize, rasterInterpretation, lat, lon, ...
                    deltaLatNumerator, 1, deltaLonNumerator, 1);
            catch e
                part1 = 'map:spatialref:';
                part2 = 'MATLAB:map:rasterref:GeographicCellsReference:';
                part3 = 'MATLAB:map:rasterref:GeographicPostingsReference:';
                if strncmp(part1, e.identifier, numel(part1)) ...
                        || strncmp(part2, e.identifier, numel(part2)) ...
                        || strncmp(part2, e.identifier, numel(part3))
                    spatialRef = [];
                else
                    rethrow(e);
                end
            end
        else
            % A geographic raster referencing object cannot be created for
            % a (non-rectilinear) affine transformation.
            spatialRef = [];
        end
        
    case 'ModelTypeProjected'
        try
            spatialRef = map.rasterref.internal.constructMapRasterReference( ...
                rasterSize, rasterInterpretation, x, y, J, [1 1; 1 1]);
        catch e
            part = 'map:spatialref:';
            if strncmp(part, e.identifier, numel(part))
                spatialRef = [];
            else
                rethrow(e);
            end
        end
        
    otherwise
        spatialRef = [];
end

%--------------------------------------------------------------------------

function [x, y, J] = constructFromSingleTiepoint(tiepoint, scale)
% Construct the first corners and the Jacobian from a single tiepoint.

dx =  scale(1);
dy = -scale(2);
x = tiepoint(4) - dx * tiepoint(1);
y = tiepoint(5) - dy * tiepoint(2);
J = [dx 0; 0 dy];

%--------------------------------------------------------------------------

function [x, y, J] = constructFromTransformMatrix(transform)
% Construct the first corners and Jacobian from the transformation matrix.

x = transform(4);
y = transform(8);
J = transform([1 2; 5 6]);

%--------------------------------------------------------------------------
    
function gtiff_info = constructCornerCoords(gtiff_info, R)
% Construct the CornerCoords structure.

% Create a CornerCoords structure to hold the outer-edge corners.
S = struct( ...
    'X', [], ...
    'Y', [], ...
    'Row', [], ...
    'Col', [], ...
    'Lat', [], ...
    'Lon', []);
                  
if ~isempty(R)
    yi = R.YIntrinsicLimits([1 1 2 2]);
    xi = R.XIntrinsicLimits([1 2 2 1]);
    S.Row = yi;
    S.Col = xi;
    if strcmp('planar', R.CoordinateSystemType)
        [xw, yw] = R.intrinsicToWorld(xi, yi);
        S.X = xw;
        S.Y = yw;
    else
        [lat, lon] = R.intrinsicToGeographic(xi, yi);
        S.X = lon;
        S.Y = lat;
    end   
    
    % Calculate latitude and longitude values.
    switch gtiff_info.ModelType
        case 'ModelTypeGeographic'
            S.Lat = S.Y;
            S.Lon = S.X;
            
        case 'ModelTypeProjected'
            if supportedCTProjection(gtiff_info.CTProjection)
                [lat, lon] = projinv(gtiff_info, S.X, S.Y);
                S.Lat = lat;
                S.Lon = lon;
            end
                           
        case 'ModelTypeGeocentric'
        otherwise
    end       
end

% Assign S to gtiff_info.
gtiff_info.CornerCoords = S;

%--------------------------------------------------------------------------

function boundingBox = constructBoundingBox(spatialRef)
% Construct the boundingBox array.
                  
if isempty(spatialRef)   
   % The spatial referencing in the GeoTIFF file is not defined, or the 
   % image size is invalid, set boundingBox to [].
   boundingBox = [];   
else
   % Create the boundingBox from the spatial referencing object.
   if strcmp('planar', spatialRef.CoordinateSystemType)
       x = spatialRef.XWorldLimits;
       y = spatialRef.YWorldLimits;
   else
       x = spatialRef.LongitudeLimits;
       y = spatialRef.LatitudeLimits;
   end
   boundingBox = [x' y'];
end

%--------------------------------------------------------------------------

function refMatrix = constructRefMatrix(spatialRef)
% Construct the referencing matrix, refMatrix, from the spatial referencing
% object, spatialRef.

if ~isempty(spatialRef)
    refMatrix = map.internal.referencingMatrix(spatialRef.worldFileMatrix);
else
    refMatrix = [];
end

%--------------------------------------------------------------------------

function TiePoints = constructTiePoints(tags)
% Construct the TiePoints structure based on the values in the GeoTIFF 
% ModelTiepoint tag. tags is a structure containing all of the GeoTIFF
% tags.

% If the ModelTiepointTag is not present, then set it to [].
if ~isfield(tags, 'ModelTiepointTag')
    tiepoints = [];
else
    tiepoints = tags.ModelTiepointTag';
end

% Create an empty TiePoints struct.
imageCoordinates = struct('Row', [], 'Col', []);
mapCoordinates   = struct('X',   [], 'Y',   []);
TiePoints = struct('ImagePoints', imageCoordinates, ...
                   'WorldPoints', mapCoordinates);  

% n is the number of coordinates per tie point 
% (one triplet each of image and world coordinates)
n = 6;
numTiePoints = numel(tiepoints)/n;

if numTiePoints >= 1
    % The gtiff_info.TiePoints field is an array of 
    % GeoTIFF TiePoint coordinate values in the form 
    % (... I, J, K, X, Y, Z, ...); where
    % I, J, K correspond to the image coordinates, and
    % X, Y, Z correspond to the world coordinates.
    %
    % Referenced in the GeoTIFF specification at:
    % http://geotiff.maptools.org/spec/geotiff2.6.html#2.6
    %
    % Create the TiePoints struct with each field name
    % and set the corresponding GTIFF TiePoints values. 
    % Convert to 1-based image-coordinates.

    imageCoordinates.Row = (tiepoints(2:n:end) + .5)';
    imageCoordinates.Col = (tiepoints(1:n:end) + .5)';
    % Z (unused) is defined:
    % imageCoordinates.Z = [tiepoints(3:n:end)]';
    TiePoints.ImagePoints = imageCoordinates;

    mapCoordinates.X = (tiepoints(4:n:end))';
    mapCoordinates.Y = (tiepoints(5:n:end))';
    % Z (unused) is defined:
    % mapCoordinateStruct.Z = [gtiff_info.TiePoints(6:n:end)]';
    TiePoints.WorldPoints = mapCoordinates;
end

%------------------------- copy functions ---------------------------------

function info = copyinfo(tiff_info, gtiff_info)
% Create a new info structure and copy the input data to it.

% Initialize a new structure.
info = infostruct(tiff_info);

% Copy TIFF data.
info = copyTiffFields(info, tiff_info);

% Copy GeoTIFF data.
info = copyGeoTiffFields(info, gtiff_info);

%--------------------------------------------------------------------------

function info = copyTiffFields(info, tiff_info)
% Copy TIFF info fields to info structure.
              
% The tiff_info structure may contain multiple elements.
tiffFields = fieldnames(tiff_info);
for k=1:numel(tiffFields)
   fieldName = tiffFields{k};
    [info.(fieldName)] = deal(tiff_info.(fieldName));
end

% If the ImageDescription field is set, copy it to the info structure.
if isfield(tiff_info,'ImageDescription')
   [info.ImageDescription] = deal(tiff_info.ImageDescription);
end

%--------------------------------------------------------------------------

function info = copyGeoTiffFields(info, gtiff_info)
% Copy GeoTIFF info fields to info structure.

% Set Zone value to [] if undefined.
if gtiff_info.Zone == gtiff_info.UserDefined
   [gtiff_info.Zone] = deal([]);
end

% Remove UserDefined
gtiff_info = rmfield(gtiff_info, 'UserDefined');

% Copy the GeoTIFF info fields. The gtiff_info structure is scalar.
geoTiffFields = fieldnames(gtiff_info);
for k=1:numel(geoTiffFields)
    fieldName = geoTiffFields{k};
    [info.(fieldName)] = deal(gtiff_info.(fieldName));
end

%------------------------- utility functions ------------------------------

function numdirs = numberOfDirectories(t)
% Obtain the number of valid TIFF directories.

numdirs = 1;
nofailures = true;
while ~t.lastDirectory && nofailures
    try
       t.nextDirectory;
       numdirs = t.currentDirectory;
    catch e
        if strcmp(e.identifier, 'MATLAB:imagesci:Tiff:unableToReadDir')
            % The next directory exists but cannot be read. numdirs will
            % not be updated and contains the number of valid directories
            % (one less than this directory number).
            nofailures = false;
        else
            % An unknown error occurred.
            rethrow(e)
        end
    end
end

% Return to the first directory.
t.setDirectory(1);

%--------------------------------------------------------------------------

function GeoKeyDirectoryTag = translateGeoKeyDirectoryTag(tags)
% Translate the TIFF GeoKeyDirectoryTag to a structure that represents
% the GeoTIFF geokey directory.

tag = tags.GeoKeyDirectoryTag;
validateGeoKeyDirectoryTag(tag);

% Construct the key ID to name map.
idMap = constructGeoKeyDirectoryMap;

% Obtain tag ID numbers for special data storage.
geoAscii  = Tiff.TagID.GeoASCIIParamsTag;
geoDouble = Tiff.TagID.GeoDoubleParamsTag;

% Translate the keys.
numKeys = tag(1,end);
for k = 2:numKeys + 1
    id = tag(k,1);
    if idMap.isKey(id)
        key = idMap(id);
    else
        key = 'Unknown';
    end
    
    location = tag(k, 2);
    switch location
        case 0
            value = tag(k,end);
            
        case geoAscii
            value = getGeoASCIIParamsTagValue(tags, tag, k);
            
        case geoDouble
            value = getGeoDoubleParamsTagValue(tags, tag, k);
            
        otherwise
            value = 'Unknown';
    end
    GeoKeyDirectoryTag.(key) = value;
end

%--------------------------------------------------------------------------

function validateGeoKeyDirectoryTag(tag)
map.internal.assert(~isempty(tag), ...
    'map:geotiff:expectedNonEmptyGeoKeyDirectoryTag');

numKeys = tag(1,end);
map.internal.assert(size(tag,2) == 4 && size(tag,1) >= numKeys+1, ...
    'map:geotiff:invalidGeoKeyDirectoryTag', ...
    numKeys+1, size(tag,1), size(tag,2));

%--------------------------------------------------------------------------

function  value = getGeoASCIIParamsTagValue(tags, geoKeyDirectoryTag, k)
% Get GeoASCIIParamsTag value from GeoTIFF tags structure, tags. The
% geoKeyDirectoryTag value provides count and offset values.

asciiTag = 'GeoASCIIParamsTag';
if ~isfield(tags, asciiTag) || isempty(tags.(asciiTag))
    geoAsciiUnknown = 'GeoAsciiParamsTag is not in file.';
    value = geoAsciiUnknown;
else
    count  = geoKeyDirectoryTag(k, 3) - 1; % Remove trailing |
    offset = geoKeyDirectoryTag(k, 4) + 1;
    value = tags.GeoASCIIParamsTag;
    startIndex = min(offset, numel(value));
    endIndex = min(offset+count-1, numel(value));
    value = value(startIndex:endIndex);
end

%--------------------------------------------------------------------------

function  value = getGeoDoubleParamsTagValue(tags, geoKeyDirectoryTag, k)
% Get GeoDoubleParamsTag value from GeoTIFF tags structure, tags. The
% geoKeyDirectoryTag value provides count and offset values.

doubleTag = 'GeoDoubleParamsTag';
if ~isfield(tags, doubleTag) || isempty(tags.(doubleTag))
    geoDoubleUnknown = 'GeoDoubleParamsTag is not in file.';
    value = geoDoubleUnknown;
else
    count  = geoKeyDirectoryTag(k, 3);
    offset = geoKeyDirectoryTag(k, 4) + 1;
    value = tags.GeoDoubleParamsTag;
    startIndex = min(offset, numel(value));
    endIndex = min(offset+count-1, numel(value));
    value = value(startIndex:endIndex);
end

%--------------------------------------------------------------------------

function tag = translateRPCCoefficientTag(value)
% Translate the TIFF RPCCoefficientTag value to a scalar
% map.geotiff.RPCCoefficientTag that represents the GeoTIFF
% RPCCoefficientTag.
% For reference, see: http://geotiff.maptools.org/rpc_prop.html 

expectedLengthOfTag = 92;
actualLengthOfTag = length(value);
if actualLengthOfTag == expectedLengthOfTag
    tag = map.geotiff.RPCCoefficientTag(value);
else
    if actualLengthOfTag ~= 0
        % Invalid tag. Issue a warning rather than an error since this is
        % not a critical tag for GeoTIFF.
        warning(message('map:geotiff:invalidRPCCoefficientTag', actualLengthOfTag));
    end
    tag = [];
end

%--------------------------------------------------------------------------

function tf = isconsistent(varargin)
% True if there are two or more inputs that are numerically equal (see
% ISEQUAL) or if there are less than two inputs.

tf = (nargin < 2) || isequal(varargin{:});

%--------------------------------------------------------------------------

function [name, translationMap] = codeToCTProjection(code)
% Translate a code number to a CTProjection name.
% A code number may be either a GeoTIFF Coordinate Transformation code
% number, which is set via the ProjCoordTransGeoKey, or a code number
% obtained from the EPSG tables.
%
% The GeoTIFF Coordinate Transformation code values range from 1:27.
% These codes are found in the GeoTIFF specification at:
% http://geotiff.maptools.org/spec/geotiff6.html#6.3.3.1
% 
% If the ProjCoordTransGeoKey is not specified, then a code number is
% obtained from the EPSG tables. In this case, the code values range from
% 9801 through 9822.

persistent ctTranslationMap
if isempty(ctTranslationMap) || ~isa(ctTranslationMap, 'containers.Map')
    names = { ...
        'CT_TransverseMercator', ...
        'CT_TransvMercator_Modified_Alaska', ...
        'CT_ObliqueMercator', ...
        'CT_ObliqueMercator_Laborde', ...
        'CT_ObliqueMercator_Rosenmund', ...
        'CT_ObliqueMercator_Spherical', ...
        'CT_Mercator', ...
        'CT_LambertConfConic_2SP', ...
        'CT_LambertConfConic_1SP', ...
        'CT_LambertAzimEqualArea', ...
        'CT_AlbersEqualArea', ...
        'CT_AzimuthalEquidistant', ...
        'CT_EquidistantConic', ...
        'CT_Stereographic', ...
        'CT_PolarStereographic', ...
        'CT_ObliqueStereographic', ...
        'CT_Equirectangular', ...
        'CT_CassiniSoldner', ...
        'CT_Gnomonic', ...
        'CT_MillerCylindrical', ...
        'CT_Orthographic', ...
        'CT_Polyconic', ...
        'CT_Robinson', ...
        'CT_Sinusoidal', ...
        'CT_VanDerGrinten', ...
        'CT_NewZealandMapGrid', ...
        'CT_TransvMercator_SouthOrientated', ...
        'CT_CylindricalEqualArea', ...
        'CT_PolarStereographic_B', ...
        'CT_LambertConfConic_1SP', ...
        'CT_LambertConfConic_2SP', ...
        'CT_LambertConfConic_2SP', ...
        'CT_Mercator', ...
        'CT_Mercator', ...
        'CT_CassiniSoldner', ...
        'CT_TransverseMercator', ...
        'CT_TransvMercator_SouthOriented', ...
        'CT_ObliqueStereographic', ...
        'CT_PolarStereographic', ...
        'CT_NewZealandMapGrid', ...
        'CT_ObliqueMercator', ...
        'CT_ObliqueMercator_Laborde', ...
        'CT_ObliqueMercator_Rosenmund', ...
        'CT_ObliqueMercator', ...
        '', ...
        'CT_AlbersEqualArea', ...
        'CT_PolarStereographic_B', ...
        ''};
    
    % Special considerations:
    % 9816 is Tunisia Mining Grid which has no CT projection, set to ''.
    % 9829 is CT_PolarStereographic_B, which is only used internally.
    % It is translated later back to CT_PolarStereographic.
    values = { ...
        1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19, ...
        20,21,22,23,24,25,26,27,28, 29, ...
        9801, 9802, 9803, 9804, 9805, 9806, 9807, 9808, 9809, ...
        9810, 9811, 9812, 9813, 9814, 9815, 9816, 9822, 9829 ...
        32767};
    ctTranslationMap = containers.Map(values, names);
end

try
    name = ctTranslationMap(code);
    
    % The Map object converts empty strings ('') to 1-by-0 empty. If name
    % is returned empty, then reset it to the empty string ('').
    if isempty(name)
        name = '';
    end
catch e %#ok<NASGU>
    name = '';
end

% This is the location of the element with value 9801 in the values array.
startEPSGIndex = 30;
if nargout > 1
    keys = ctTranslationMap.keys;
    values = ctTranslationMap.values;
    % Remove EPSG codes.
    keys(startEPSGIndex:end-1) = [];
    values(startEPSGIndex:end-1) = [];
    translationMap = containers.Map(values, keys);
end

%--------------------------------------------------------------------------

function tf = supportedCTProjection(name)
% Return true if NAME is a projection string that PROJINV supports.

unsupported = { ...
    'CT_TransvMercator_Modified_Alaska', ...
    'CT_ObliqueMercator_Laborde', ...
    'CT_ObliqueMercator_Rosenmund', ...
    'CT_ObliqueMercator_Spherical', ...
    'CT_NewZealandMapGrid', ...
    'CT_TransvMercator_SouthOrientated', ...
    ' '};
tf = ~isempty(name) && ischar(name) && ~any(strcmp(name, unsupported)) ...
    && length(name) > 3 && strcmp(name(1:3), 'CT_');
    
%--------------------------------------------------------------------------

function [projParm, projParmId] = updateProjParm( ...
    projParm, projParmId, geokeys, name, id)
% Update projection parameters if geokeys have been set.

keynames = fieldnames(geokeys);

names = { ...
    'ProjFalseEastingGeoKey', ...
    'ProjCenterEastingGeoKey', ...
    'ProjFalseOriginEastingGeoKey'};
falseEasting = keyget(geokeys, keynames, names, 0);

names = { ...
    'ProjFalseNorthingGeoKey', ...
    'ProjCenterNorthingGeoKey', ...
    'ProjFalseOriginNorthingGeoKey'};
falseNorthing = keyget(geokeys, keynames, names, 0);

names = { ...
    'ProjNatOriginLongGeoKey', ...
    'ProjFalseOriginLongGeoKey', ...
    'ProjCenterLongGeoKey'};
natOriginLong = keyget(geokeys, keynames, names, 0);

names = { ...
    'ProjNatOriginLatGeoKey', ...
    'ProjFalseOriginLatGeoKey', ...
    'ProjCenterLatGeoKey'};
natOriginLat = keyget(geokeys, keynames, names, 0);

names = {'ProjScaleAtNatOriginGeoKey'};
natOriginScale = keyget(geokeys, keynames, names, 1.0);

switch name
    
    case { ...
            'CT_Stereographic', ...
            'CT_LambertConfConic_1SP',  ...
            'CT_Mercator', ...
            'CT_ObliqueStereographic', ...
            'CT_TransverseMercator', ...
            'CT_TransvMercator_SouthOriented'}
        projParm(1)   = natOriginLat;
        projParm(2)   = natOriginLong;
        projParm(5)   = natOriginScale;
        projParm(6)   = falseEasting;
        projParm(7)   = falseNorthing;
        
        projParmId(1) = id('ProjNatOriginLatGeoKey');
        projParmId(2) = id('ProjNatOriginLongGeoKey');
        projParmId(5) = id('ProjScaleAtNatOriginGeoKey');
        projParmId(6) = id('ProjFalseEastingGeoKey');
        projParmId(7) = id('ProjFalseNorthingGeoKey');
        
    case 'CT_ObliqueMercator'
        name = {'ProjAzimuthAngleGeoKey'};
        azimuth = keyget(geokeys, keynames, name, 0);
        name = {'ProjRectifiedGridAngleGeoKey'};
        rectGridAngle = keyget(geokeys, keynames, name, 90.0);
        names = {'ProjScaleAtNatOriginGeoKey', 'ProjScaleAtCenterGeoKey'};
        natOriginScale = keyget(geokeys, keynames, names, 1.0);
        
        projParmId(1) = id('ProjCenterLatGeoKey');
        projParmId(2) = id('ProjCenterLongGeoKey');
        projParmId(3) = id('ProjAzimuthAngleGeoKey');
        projParmId(4) = id('ProjRectifiedGridAngleGeoKey');
        projParmId(5) = id('ProjScaleAtCenterGeoKey');
        projParmId(6) = id('ProjFalseEastingGeoKey');
        projParmId(7) = id('ProjFalseNorthingGeoKey');
        
        projParm(1)   = natOriginLat;
        projParm(2)   = natOriginLong;
        projParm(3)   = azimuth;
        projParm(4)   = rectGridAngle;
        projParm(5)   = natOriginScale;
        projParm(6)   = falseEasting;
        projParm(7)   = falseNorthing;
        
    case {'CT_CassiniSoldner' ,'CT_Polyconic'}
        names = {'ProjScaleAtNatOriginGeoKey', 'ProjScaleAtCenterGeoKey'};
        natOriginScale = keyget(geokeys, keynames, names, 1.0);
        
        projParmId(1) = id('ProjNatOriginLatGeoKey');
        projParmId(2) = id('ProjNatOriginLongGeoKey');
        projParmId(5) = id('ProjScaleAtNatOriginGeoKey');
        projParmId(6) = id('ProjFalseEastingGeoKey');
        projParmId(7) = id('ProjFalseNorthingGeoKey');
        
        projParm(1)   = natOriginLat;
        projParm(2)   = natOriginLong;
        projParm(5)   = natOriginScale;
        projParm(6)   = falseEasting;
        projParm(7)   = falseNorthing;
        
    case { ...
            'CT_AzimuthalEquidistant', ...
            'CT_MillerCylindrical',  ...
            'CT_Gnomonic', ...
            'CT_LambertAzimEqualArea', ...
            'CT_Orthographic', ...
            'CT_NewZealandMapGrid'}
        
        projParmId(1) = id('ProjCenterLatGeoKey');
        projParmId(2) = id('ProjCenterLongGeoKey');
        projParmId(6) = id('ProjFalseEastingGeoKey');
        projParmId(7) = id('ProjFalseNorthingGeoKey');
        
        projParm(1)   = natOriginLat;
        projParm(2)   = natOriginLong;
        projParm(6)   = falseEasting;
        projParm(7)   = falseNorthing;
        
    case 'CT_Equirectangular'
        
        projParmId(1) = id('ProjCenterLatGeoKey');
        projParmId(2) = id('ProjCenterLongGeoKey');
        projParmId(3) = id('ProjStdParallel1GeoKey');
        projParmId(6) = id('ProjFalseEastingGeoKey');
        projParmId(7) = id('ProjFalseNorthingGeoKey');
        
        names = {'ProjStdParallel1GeoKey'};
        stdParallel1 = keyget(geokeys, keynames, names, 0.0);
        
        projParm(1)   = natOriginLat;
        projParm(2)   = natOriginLong;
        projParm(3)   = stdParallel1;
        projParm(6)   = falseEasting;
        projParm(7)   = falseNorthing;

    case { ...
            'CT_Robinson', ...
            'CT_Sinusoidal',  ...
            'CT_VanDerGrinten'}
        
        projParmId(2) = id('ProjNatOriginLongGeoKey');
        projParmId(6) = id('ProjFalseEastingGeoKey');
        projParmId(7) = id('ProjFalseNorthingGeoKey');
        
        projParm(2)   = natOriginLong;
        projParm(6)   = falseEasting;
        projParm(7)   = falseNorthing;

    case {'CT_PolarStereographic'}
        names = { ...
            'ProjStraightVertPoleLongGeoKey', ...
            'ProjNatOriginLongGeoKey', ...
            'ProjFalseOriginLongGeoKey', ...
            'ProjCenterLongGeoKey'};
        natOriginLong = keyget(geokeys, keynames, names, 0);
        
        names = {'ProjScaleAtNatOriginGeoKey', 'ProjScaleAtCenterGeoKey'};
        natOriginScale = keyget(geokeys, keynames, names, 1.0);
        
        projParmId(1) = id('ProjNatOriginLatGeoKey');
        projParmId(2) = id('ProjStraightVertPoleLongGeoKey');
        projParmId(5) = id('ProjScaleAtNatOriginGeoKey');
        projParmId(6) = id('ProjFalseEastingGeoKey');
        projParmId(7) = id('ProjFalseNorthingGeoKey');
        
        projParm(1)   = natOriginLat;
        projParm(2)   = natOriginLong;
        projParm(5)   = natOriginScale;
        projParm(6)   = falseEasting;
        projParm(7)   = falseNorthing;

    case {'CT_LambertConfConic_2SP'}
        names = {'ProjStdParallel1GeoKey'};
        stdParallel1 = keyget(geokeys, keynames, names, 0.0);
        names = {'ProjStdParallel2GeoKey'};
        stdParallel2 = keyget(geokeys, keynames, names, 0.0);
        
        projParmId(1) = id('ProjFalseOriginLatGeoKey');
        projParmId(2) = id('ProjFalseOriginLongGeoKey');
        projParmId(3) = id('ProjStdParallel1GeoKey');
        projParmId(4) = id('ProjStdParallel2GeoKey');
        projParmId(6) = id('ProjFalseEastingGeoKey');
        projParmId(7) = id('ProjFalseNorthingGeoKey');
        
        projParm(1)   = natOriginLat;
        projParm(2)   = natOriginLong;
        projParm(3)   = stdParallel1;
        projParm(4)   = stdParallel2;
        projParm(6)   = falseEasting;
        projParm(7)   = falseNorthing;
        
    case { ...
            'CT_AlbersEqualArea', ...
            'CT_EquidistantConic'}
        names = {'ProjStdParallel1GeoKey'};
        stdParallel1 = keyget(geokeys, keynames, names, 0.0);
        names = {'ProjStdParallel2GeoKey'};
        stdParallel2 = keyget(geokeys, keynames, names, 0.0);
        
        projParmId(1) = id('ProjStdParallel1GeoKey');
        projParmId(2) = id('ProjStdParallel2GeoKey');
        projParmId(3) = id('ProjNatOriginLatGeoKey');
        projParmId(4) = id('ProjNatOriginLongGeoKey');
        projParmId(6) = id('ProjFalseEastingGeoKey');
        projParmId(7) = id('ProjFalseNorthingGeoKey');
        
        projParm(1)   = stdParallel1;
        projParm(2)   = stdParallel2;
        projParm(3)   = natOriginLat;
        projParm(4)   = natOriginLong;
        projParm(6)   = falseEasting;
        projParm(7)   = falseNorthing;
        
    case {'CT_CylindricalEqualArea'}
        names = {'ProjStdParallel1GeoKey'};
        stdParallel1 = keyget(geokeys, keynames, names, 0.0);
        
        projParmId(1) = id('ProjStdParallel1GeoKey');
        projParmId(2) = id('ProjNatOriginLongGeoKey');
        projParmId(6) = id('ProjFalseEastingGeoKey');
        projParmId(7) = id('ProjFalseNorthingGeoKey');
        
        projParm(1)   = stdParallel1;
        projParm(2)   = natOriginLong;
        projParm(6)   = falseEasting;
        projParm(7)   = falseNorthing;
end

%--------------------------------------------------------------------------

function gtiff_info = adjustPolarStereographicBParams(gtiff_info, inverseMap)
% Adjust the PolarStereographic_B projection parameters to conform to a 
% regular Polar Stereographic projection by computing the scale value and
% changing the CT name.

try
    % Obtain eccentricity.
    ellipsoidCode = gtiff_info.GeoTIFFCodes.Ellipsoid;
    ellipsoid = referenceEllipsoid(ellipsoidCode);
    e = ellipsoid.Eccentricity;
    
    % Calculate scale value.
    phiFd = gtiff_info.ProjParm(1);
    k = computePolarStereographicB_Scale(phiFd, e);
    gtiff_info.ProjParm(5) = k;
    
    % Change name and number to CT_PolarStereographic
    projname = 'CT_PolarStereographic';
    gtiff_info.CTProjection = projname;
    gtiff_info.GeoTIFFCodes.CTProjection = inverseMap(projname);

catch e %#ok<NASGU>
    % Nothing to change. Ellipsoid information is unavailable.
end

%-----------------------------------------------------------------------

function k = computePolarStereographicB_Scale(phiFd, e)
% Compute PolarStereographic_B scale value.

% Compute k using formulas from EPSG.
% Reference: Formulas are obtained from
% http://www.iogp.org/pubs/373-07-2.pdf
% Geomatics
% Guidance Note Number 7, part 2
% Coordinate Conversions and Transformations including Formulas
% Revised - April 2012
% page 67: Polar Stereographic (Variant B)
%   (EPSG dataset coordinate operation method code 9829).

phiF = deg2rad(phiFd);
esinphiF = e * sin(phiF);
v = (1 + esinphiF) / (1 - esinphiF);
if( phiFd >= 0.0)
    % Northern hemisphere.
    tF = tan(pi/4 - phiF/2) * power(v,e/2);
else
    % Southern hemisphere.
    tF = tan(pi/4 + phiF/2) / power(v,e/2);
end

mF = cos(phiF) / sqrt(1 - esinphiF^2);
v1 = power(1+e, 1+e);
v2 = power(1-e, 1-e);
k = mF * sqrt(v1*v2)/(2*tF);

%--------------------------------------------------------------------------

function tf = cellcmp(a,b)
% Compare cell or string arrays. A is either a cell array of strings or a
% string. B is either a cell array of strings or a string. If A is a
% string, then TF is a scalar logical and is true if any member of B
% matches A. If A is a cell array, then TF is the same size as A and is
% true for each member of A that is in B. cellcmp returns the same values
% as ismember but is significantly faster.

if ischar(a) || (ischar(a) && ischar(b))
    tf = any(strcmp(a,b));
else
    tf = false(size(a));
    for k=1:length(a)
        tf(k) = any(strcmp(a{k}, b));
    end
end

%--------------------------------------------------------------------------

function angleDeg = angleToDegrees(angleVal, code)
% Convert the angle, angleVal, to a decimal value, angleDeg, based on the
% value of the EPSG code.

switch code
    case 9110
        % Code 9110 represents a sexagesimal DMS angle.
        angleDeg = sexagesimalAngleToDegrees(angleVal);
        
    case {9105,  9106}
        % grad
        angleDeg = 180 * (angleVal / 200);
        
    case 9101
        % radians
        angleDeg = 180 * (angleVal / pi);
        
    case 9103
        % arc-minute
        angleDeg = angleVal / 60;
        
    case 9104
        % arc-second
        angleDeg = angleVal / 3600;
        
    otherwise
        % decimal degrees
        % Some cases missing but seemingly never used.
        angleDeg = angleVal;
end

%--------------------------------------------------------------------------

function angleDeg = sexagesimalAngleToDegrees(angleStr)
% Convert sexagesmialAngleString to degrees.

% A sexagesimal DMS angle is in the form:
%   DDD.MMSSsss
% From the EPSG unit_of_measure table:
% Pseudo unit. Format:
%    signed degrees - period -
%    minutes (2 digits)
%    integer seconds (2 digits)
%    fraction of seconds (any precision).
% Must include leading zero in minutes and seconds and exclude
% decimal point for seconds.

angleStr = sprintf('%.7f',angleStr);
index = strfind(angleStr, '.');
if isempty(index)
    ddd = angleStr;
    mm = '0';
    ss = '0';
else
    ddd = angleStr(1:index-1);
    angleStr(1:index) = [];
    switch length(angleStr)
        case 0
            mm = '0';
            ss = '0';
        case 1
            mm = [angleStr '0'];
            ss = '0';
        case 2
            mm = angleStr;
            ss = '0';
        case 3
            mm = angleStr(1:2);
            ss = [angleStr(3) '0'];
        case 4
            mm = angleStr(1:2);
            ss = angleStr(3:4);
        otherwise
            mm = angleStr(1:2);
            ss = angleStr(3:4);
            ss_fraction = ['.' angleStr(5:end)];
            ss = [ss ss_fraction];
    end
end
ddd = str2double(ddd);
mm = str2double(mm);
ss = str2double(ss);
if mm > 60 || ss > 60
    angleDeg = 0;
else
    angleDeg = dms2degrees([ddd, mm, ss]);
end

% ---------------- copied by Shi
function [idToNameMap, nameToIdMap] = constructGeoKeyDirectoryMap
%CONSTRUCTGEOKEYDIRECTORYMAP Construct maps for GeoKeyDirectoryTag
%
%   [idToNameMap, nameToIdMap] = constructGeoKeyDirectoryMap maps GeoTIFF
%   key ID (GeoKey) to key name. The function returns one or two maps based
%   on the number of requested outputs. If one output is requested, the
%   output parameter maps the GeoKey ID to its name. If two outputs are
%   requested, the second output parameter maps the key name (in lower
%   case) back to the GeoKey ID. The names and values are obtained from the
%   GeoTIFF specification.

% Copyright 2010-2011 The MathWorks, Inc.

% Note: Added ProjRectifiedGridAngleGeoKey to map to code 3096. This code
% is not listed in the GeoTIFF specification at:
% http://www.remotesensing.org/geotiff/spec/geotiff6.html#6.2.3
% but is included in libgeotiff.

keyName = { ...
    'GTModelTypeGeoKey', ...
    'GTRasterTypeGeoKey', ...
    'GTCitationGeoKey', ...
    'GeographicTypeGeoKey', ...
    'GeogCitationGeoKey', ...           
    'GeogGeodeticDatumGeoKey', ...      
    'GeogPrimeMeridianGeoKey', ...      
    'GeogLinearUnitsGeoKey', ...        
    'GeogLinearUnitSizeGeoKey', ...     
    'GeogAngularUnitsGeoKey', ...       
    'GeogAngularUnitSizeGeoKey', ...    
    'GeogEllipsoidGeoKey', ...          
    'GeogSemiMajorAxisGeoKey', ...      
    'GeogSemiMinorAxisGeoKey', ...     
    'GeogInvFlatteningGeoKey', ...      
    'GeogAzimuthUnitsGeoKey', ...       
    'GeogPrimeMeridianLongGeoKey', ...  
    'ProjectedCSTypeGeoKey', ...          
    'PCSCitationGeoKey', ...              
    'ProjectionGeoKey', ...               
    'ProjCoordTransGeoKey', ...           
    'ProjLinearUnitsGeoKey', ...          
    'ProjLinearUnitSizeGeoKey', ...       
    'ProjStdParallel1GeoKey', ...         
    'ProjStdParallel2GeoKey', ...         
    'ProjNatOriginLongGeoKey', ...       
    'ProjNatOriginLatGeoKey', ...         
    'ProjFalseEastingGeoKey', ...         
    'ProjFalseNorthingGeoKey', ...        
    'ProjFalseOriginLongGeoKey', ...
    'ProjFalseOriginLatGeoKey', ...
    'ProjFalseOriginEastingGeoKey', ...   
    'ProjFalseOriginNorthingGeoKey', ...  
    'ProjCenterLongGeoKey', ...           
    'ProjCenterLatGeoKey', ...            
    'ProjCenterEastingGeoKey', ...        
    'ProjCenterNorthingGeoKey', ...       
    'ProjScaleAtNatOriginGeoKey', ...     
    'ProjScaleAtCenterGeoKey', ...        
    'ProjAzimuthAngleGeoKey', ...         
    'ProjStraightVertPoleLongGeoKey', ...
    'ProjRectifiedGridAngleGeoKey', ...
    'VerticalCSTypeGeoKey', ...           
    'VerticalCitationGeoKey', ...         
    'VerticalDatumGeoKey', ...            
    'VerticalUnitsGeoKey', ...
    'UserDefined'};

keyID = { ...
   1024, ... 
   1025, ...
   1026, ...
   2048, ... 
   2049, ... 
   2050, ... 
   2051, ... 
   2052, ... 
   2053, ...
   2054, ...
   2055, ... 
   2056, ... 
   2057, ... 
   2058, ... 
   2059, ... 
   2060, ... 
   2061, ... 
   3072, ... 
   3073, ... 
   3074, ... 
   3075, ... 
   3076, ...
   3077, ... 
   3078, ... 
   3079, ... 
   3080, ... 
   3081, ... 
   3082, ... 
   3083, ... 
   3084, ... 
   3085, ... 
   3086, ... 
   3087, ...
   3088, ... 
   3089, ... 
   3090, ...
   3091, ... 
   3092, ...
   3093, ... 
   3094, ... 
   3095, ... 
   3096, ...
   4096, ... 
   4097, ... 
   4098, ... 
   4099, ...
   32767};

idToNameMap = containers.Map(keyID, keyName);
if nargout > 1
    nameToIdMap = containers.Map(lower(keyName), keyID);
end
