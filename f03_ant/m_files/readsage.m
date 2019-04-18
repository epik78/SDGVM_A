function yie = readsage(crop)

    %Reads data from SAGE yield database
    dirr{1} = pwd;
    
    %Where to dump
    dirr{2} = '/fastdata-sharc/sm1epk/';
    str{1}=[dirr{2},'HarvestedAreaYield175Crops_NetCDF'];


    %Gets the names of all the crops in the dataset
    cd(dirr{2})
    system('rm -R HarvestedAreaYield175Crops_NetCDF');
    %Sends the zip file with the SAGE data to a folder with more space
    cd /data/sm1epk/crop_sets/SAGE
    system(['cp 175CropsYieldArea_netcdf.zip ',dirr{2}]); 

    cd(dirr{2})
    %Unzip without files
    system('unzip -qq 175CropsYieldArea_netcdf.zip -x "*.nc" "*.png" "*.pdf"');

    %Get the names of all the crop types in the dataset
    cd([dirr{2},'HarvestedAreaYield175Crops_NetCDF']);
    crop_nam=dir;
    crop_nam=crop_nam(3:end);
    %crop_names holds all the crop names in sage
    for i=1:size(crop_nam,1)
        sub=crop_nam(i).name;
        crop_names{i}=strtok(sub,'_');
    end

    cd(dirr{2})
    system('rm -R HarvestedAreaYield175Crops_NetCDF');clear crop_nam sub;

    %Find the index of the specific crop Im looking at
    idx = find(strcmp(lower(crop),crop_names));

    %Get parameters of the specific crop
    cd(dirr{1})
    param=crop_param(idx);

    %Unzip only specific crop data
    cd(dirr{2})
    system(['unzip -qq 175CropsYieldArea_netcdf.zip HarvestedAreaYield175Crops_NetCDF/',crop_names{idx},'_HarvAreaYield_NetCDF/* -d ',dirr{2}]);

    cd([str{1},'/',crop_names{idx},'_HarvAreaYield_NetCDF']);
    fid=netcdf.open([crop_names{idx},'_AreaYieldProduction.nc']);
    var=netcdf.getVar(fid,4);
    
    %Reads yield in metric tons per hectare
    yie = squeeze(var(:,:,2)');
    
    netcdf.close(fid);

    cd(dirr{1})


end