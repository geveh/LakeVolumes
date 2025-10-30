%% Read data on lakes and tiles
LAKES = readgeotable('glacier_lakes_2020_centroids.shp');


%% Get merit extents
rootdir = 'E:\USER\schwanghart\meritEARTH\merit_tiles';
files = dir(fullfile(rootdir, '**\*.tif'));

for r = 1:numel(files)
    % Get raster information
    info = georasterinfo(fullfile(files(r).folder,files(r).name));

    % Extract limits (needs to be adjusted if raster comes in
    % geographic coordinates
    lon = info.RasterReference.LongitudeLimits;
    lat = info.RasterReference.LatitudeLimits;

    % Create a mapshape
    p = geopolyshape([lat(1) lat(2) lat(2) lat(1) lat(1)],[lon(1) lon(1) lon(2) lon(2) lon(1)]);
    % Add coordinate reference system to the shape
    p.GeographicCRS = info.RasterReference.GeographicCRS;

    % Place mapshape object and file name into a table (this will
    % become a geotable)
    GT = table(p,string(files(r).name),'VariableNames',{'Shape', 'filename'});

    if r == 1
        GTFULL = GT;
    else
        GTFULL = [GTFULL;GT];
    end
end

extents = GTFULL;
clear GTFULL GT lat lon p r 

%% Spatial join of lakes and DEM extents

CIXL = cell(size(extents,1),1);
parfor r = 1:size(extents,1)
    inpoly = isinterior(extents.Shape(r),LAKES.Shape);
    CIXL{r} = find(inpoly);
end

%% Create a neighborhood graph of extents

extents2 = geotable2table(extents,["Lat" "Lon"]);
CP = cellfun(@(x,y) polyshape(x,y),extents2.Lon,extents2.Lat,'UniformOutput',false);
CP = [CP{:}];
neighs = cell(size(CP));

parfor r = 1:size(extents,1)
    if isempty(CIXL{r})
        continue
    end
    I = overlaps(polybuffer(CP(r),0.1),CP);
    I(r) = false;
    neighs{r} = find(I);
end

clear extents2
%% 

CPS = cell(size(extents,1),1);

parfor r = 1:size(extents,1)
    if isempty(CIXL{r})
        CPS{r} = [];
        continue
    end
    
    counter = 0;
    PS = struct();

    % Make sure to find the lake outlet (currently, the location is the
    % centroid (point on surface))
    L = LAKES(CIXL{r},:);

    % Read DEMs
    DEM = GRIDobj(string(files(r).folder) + filesep + string(files(r).name));
    
    % neighboring files
    nfiles = neighs{r};
    if isempty(nfiles)
        % don't do anything
    else
        CDEM = cell(numel(nfiles)+1,1);
        CDEM{1} = DEM;
        for rr = 1:numel(nfiles)
            CDEM{rr+1} = GRIDobj(string(files(nfiles(rr)).folder) + filesep + string(files(nfiles(rr)).name));
        end
        DEM  = mosaic(CDEM{:});
    end

    % Read closed lakes
    ClosedLakes = readgeotable('../GLWD/closed_lakes.shp');
    ext         = getextent(DEM);
    ClosedLakes = geoclip(ClosedLakes.Shape,ext([3 4]),ext([1 2]));
    ClosedLakes = table(ClosedLakes,(1:numel(ClosedLakes))','VariableNames',{'Shape','ID'});

    % reproject to utm
    DEM = reproject2utm(DEM,90);
    
    ClosedLakes = projectshape(ClosedLakes,DEM);
    ClosedLakes = polygon2GRIDobj(DEM,ClosedLakes);

    % Find the minimum in each closed lake
    stats = regionprops(ClosedLakes.Z,DEM.Z,"PixelIdxList");
    for label = 1:numel(stats)
        [~,ix] = min(DEM.Z(stats(label).PixelIdxList));
        DEM.Z(stats(label).PixelIdxList(ix)) = nan;
    end


    % Calculate flow network
    FD = FLOWobj(DEM,'uselibtt',false,'mex',true);

    FD.fastindexing = true;

    % Get lakes
    L = projectshape(L,DEM);
    L = geotable2mapstruct(L);
    IXL = coord2ind(DEM,[L.X]',[L.Y]');
    
    % Go through all lakes and extract their paths
    for rr = 1:numel(IXL)
        counter = counter + 1;
        ID = L(rr).UniqueID;
        [ixpath,d,x,y] = flowpathextract(FD,IXL(rr));
        z = DEM.Z(ixpath);
        [lat,lon] = projinv(parseCRS(DEM),x,y);
        PS(counter).ID = ID;
        PS(counter).d = d(:)';
        PS(counter).lat = lat(:)';
        PS(counter).lon = lon(:)';
        PS(counter).z   = z(:)';
    end
    CPS{r} = PS;
end

CPS = horzcat(CPS{:});
%%
% Open the file for writing
fid = fopen('flowpaths_detailed.txt', 'w');

% Check if file opened successfully
if fid == -1
    error('Could not open file for writing.');
end

% Write each element in the cell array to a new line
for i = 1:numel(CPS)
    fprintf(fid, '%s\n', num2str(CPS(i).ID));
    fprintf(fid, '%s\n', num2str(CPS(i).d/1000));
    fprintf(fid, '%s\n', num2str(CPS(i).lat));
    fprintf(fid, '%s\n', num2str(CPS(i).lon));
    fprintf(fid, '%s\n', num2str(CPS(i).z));
end

% Close the file
fclose(fid);
    