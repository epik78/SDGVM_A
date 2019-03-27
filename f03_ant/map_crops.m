clear;

dirr{1} = pwd;
%Where tha run is and the figures will be printed
dirr{2} = '/fastdata-sharc/sm1epk/SDGVM_runs/fin_outs';
%Where the country masks are
dirr{3} = '/data/sm1epk/crop_sets/coun_masks';
%Where the p7print function is that does the map
dirr{4} = '/data/sm1epk/EW/weath_emp_model/codes/iter_1/';

%Names of crops for file reading
crop_n = {'Maize','Soy','Millet','Swheat','Wwheat','Sbarley','Wbarley','Sorghum'};
%Names of crops for figure titles
str_n = {'Maize','Soy','Millet','Spring Wheat','Winter wheat','Spring Barley','Winter Barley','Sorghum'};
%Variables for file reading
crop_out = {'sow','harv','gdd','seas','ryield'};
%Variables for figure titles
str_out = {'Sowing Day','Harvest Day','Growing Degree Days','Seasonality','RYield'};
%Colorbar limits for each variable
llim(1,:)=[1 365];llim(2,:)=[1 365];llim(3,:)=[1 3000];llim(4,:)=[1 2];llim(5,:)=[1 200];
%Resolution of printed figure
fig_res = 200;
%Resolution of run
ress = 1.;
%lat and lon of grid cell center
latm = 90-ress/2:-ress:-90+ress/2;
lonm = -180+ress/2:ress:180-ress/2;

%Reads mask
mask=zeros(size(latm,2),size(lonm,2));
cd(dirr{3})
for c=1:246
    load(['mask_',num2str(c),'.mat'])
    s2r=imresize(s2r,[size(mask,1) size(mask,2)],'Nearest');
    mask=mask+s2r;
end
mask(mask>0)=1;

%For each crop
for ii = 1:size(crop_n,2)
    %For each variable
    for jj = 1:size(crop_out,2)
        cd(dirr{2})
        [crop_n{ii},'_',crop_out{jj}]
        %Reads file and outputs locs and data for each grid cell
        [locs,dat] = readsdfiles([crop_n{ii},'_',crop_out{jj}],dirr{2});
        
        %Finds lat and lon in pixel space
        locsi = nan(size(locs,1),2);
        for kk = 1:size(locs,1)
            [ia,ib] = min(abs(latm-locs(kk,1)));
            [ic,id] = min(abs(lonm-locs(kk,2)));
            locsi(kk,1) = ib; locsi(kk,2) = id;    
        end
        
        %Gets the year of data
        sub_dat = dat(:,end);

        %Puts the data in the indexes of the pixel space
        sub = nan(size(latm,2),size(lonm,2));
        ind = sub2ind([size(sub,1) size(sub,2)],locsi(:,1),locsi(:,2));
        sub(ind) = sub_dat';
        sub(sub==0) = NaN;

        cd(dirr{4})
        p7print(sub,[str_n{ii},' ',str_out{jj}],[str_n{ii},'_',str_out{jj}],dirr{2},fig_res,'BrBG4',1,0,[llim(jj,1) llim(jj,2)],20,mask,[-90 90],[-180 180],'Countries')
        
    end
end

cd(dirr{1})

function [locs,dat] = readsdfiles(str1,dir_data)

    %Goes in the directory where the SDGVM run is
    cd(dir_data)
    
    %Opens and reads the file
    %fid = fopen([str1,'.dat'],'rt');
    %a = fscanf(fid,'%f');
    %fclose(fid);
    a = dlmread([str1,'.dat']);

    %locs(lat/lon,grid cell)
    locs = a(:,1:2);
    %data(:,grid cell)
    dat = a(:,3:end);   

end

