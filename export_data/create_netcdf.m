function create_netcdf

out = load([gdrive filesep 'github' filesep 'GlobalDeltaChange' filesep 'GlobalDeltaData.mat']);




if isfield(out,'delta_name'),
out.delta_name = char(out.delta_name);
end
out.Region_str = char(out.Region_str);
out.Region = double(out.Region);
ftext={'wave_lat','dec deg','Wave data node latitude';...
'wave_lon','dec deg','Wave data node longitude';...
'shelf_depth','m','Continental shelf depth, based on complicated algorithm for finding the shelf break';...
'shelf_slope','','Slope of the continental shelf';...
'shelf_width','','Distance from river mouth to shelf break';...
'shelf_len', 'km','Depth profile of ocean basin w/ multiplier of 20 (indices 1 to 41 indicate distance to a depth of 0 to 40*20=800 m)';...
'shelf_len_lat', 'dec deg','Latitude of ocean basin depth profile';...
'shelf_len_lon','dec deg','Longitude of ocean basin depth profile';...
'Hs','m','Offshore Significant Waveheight';...
'QWave','kg/s','Wave-driven sediment flux';...
'Discharge_tide','m3/s','Tide-driven water discharge';...
'QTide','kg/s','Tide-driven sediment flux';...
'TidalAmp', 'm','Average tidal amplitude';...
'BasinID_ATLAS', '','Drainage Basin ID from the new HydroATLAS dataset';...
'BasinArea','km2','Drainage basin area';...
'BasinID', '','BasinID from HydroSheds';...
'BasinID2', '','Unique BasinID constructed using: 10*BasinID+Continent';...
'Continent', '','ID for different continents';...
'MouthLat', 'dec deg','River mouth latitude';...
'MouthLon', 'dec deg','River mouth longitude';...
'Region','','ID for different regions';...
'Region_str', '','Name of different regions';...
'channel_len','m','Elevation profile of river channel (indices 1,2, 20 indicate distance to 0,1, 19 m of elevation)';...
'channel_len_lat', 'dec deg','Latitude of river channel elevation data';...
'channel_len_lon','dec deg','Latitude of river channel elevation data';...
'ChannelSlope','','Channel slope in delta, from SRTM';...
'Discharge_dist','m3/s','Modern river discharge';...
'Discharge_prist', 'm3/s','Pre-dam river discharge';...
'QRiver_dist','kg/s','Modern river sediment flux';...
'QRiver_prist','kg/s','Pre-dam/pre-land-use change river sediment flux';...
'delta_name','','Name of river'};
    

fna = fieldnames(out);
[~,idx] = ismember(fna,ftext(:,1));
out = rmfield(out,fna(idx==0));

[~,idx] = ismember(fieldnames(out),ftext(:,1));
funits = ftext(idx,2);
fmeta = ftext(idx,3);

write_netcdf('GlobalDeltaData.nc',out,funits,fmeta)



function write_netcdf(ncname,out,funits,fmeta)

fnames = fieldnames(out);

ncid = netcdf.create(ncname,'WRITE');

dim_a = netcdf.defDim(ncid,'n0',size(out.(fnames{1}),1));
dim_b = netcdf.defDim(ncid,'n1',size(out.(fnames{1}),2));


for ii=1:length(fnames),
    
    if isa(out.(fnames{ii}),'integer'),
        out.(fnames{ii}) = int32(out.(fnames{ii}));
        cla = 'NC_INT';
    else,
        cla = class(out.(fnames{ii}));
    end
        
    
    if size(out.(fnames{ii}),1)==10848, 
        dim1 = dim_a;
    else,
        dim1 = netcdf.defDim(ncid,['n' num2str(ii) '_1'],size(out.(fnames{ii}),1));
    end
    
    if size(out.(fnames{ii}),2)==1, 
        dim2 = dim_b;
    else,
        dim2 = netcdf.defDim(ncid,['n' num2str(ii) '_2'],size(out.(fnames{ii}),2));
    end
        
    varid = netcdf.defVar(ncid,fnames{ii},cla,[dim1 dim2]);
        

    
    %ncwrite(ncname,fnames{ii},out.(fnames{ii}))
    netcdf.endDef(ncid);
    netcdf.putVar(ncid,varid,out.(fnames{ii}));
    netcdf.reDef(ncid);
    netcdf.putAtt(ncid,varid,'units',funits{ii});
    netcdf.putAtt(ncid,varid,'about',fmeta{ii});
    
end
    

netcdf.close(ncid);

