

out = load([dropbox filesep 'WorldDeltas' filesep 'scripts' filesep 'GlobalDeltaData.mat']);

fnames = fieldnames(out);

n = length(out.(fnames{1}));

for ii=1:length(fnames),
    
    if length(out.(fnames{ii}))~=n, continue, end
        
    nccreate('GlobalDeltaData.nc',fnames{ii},'Dimensions',{'n',n})
    ncwrite('GlobalDeltaData.nc',fnames{ii},out.(fnames{ii}))
end
    