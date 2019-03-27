clear;
dir1=pwd;

%Directory where the SDGVM run is and where I save figs
dir_data = '/fastdata-sharc/sm1epk/SDGVM_runs/fin_outs';
%Directory where the map maker is
map_d = '/data/sm1epk/EW/weath_emp_model/codes/iter_1';
%Variable I want to read.It does the yearly/monthly/daily outputs
var = 'yearly_lai_Dc_Bl';


%Reads file and outputs locs and data
[locs,dat] = readsdfiles(var,dir_data);


cd(dir1)

%Figures whether its yearly, monthly or daily
strb = extractBefore(var,'_');

switch strb
    case "yearly"
        %Gets the data of the 1st year
        dat=dat(1,:);
    case "monthly"
        %Gets the data of the 6th month of the first year
        dat=squeeze(dat(6,1,:));
    case "daily"
        %Gets the data of day 180 of the first year
        dat=squeeze(dat(180,1,:));
end

%Makes map
makessdmap(locs,dat,var,map_d,dir_data)

cd(dir1)



function [locs,dat] = readsdfiles(str1,dir_data)

    %Goes in the directory where the SDGVM run is
    cd(dir_data)
    %Reads the number of lines in the site info file
    %subtracts one for header and figures the number
    %of grid cells in the run
    [a,b] = unix('wc -l site_info.dat');
    b = extractBefore(b,' ');
    no_grid = str2num(b)-1;

    %Opens and reads the data file
    fid = fopen([str1,'.dat'],'rt');
    a = fscanf(fid,'%f');
    fclose(fid);

    %Figures whether its yearly, monthly or daily
    strb = extractBefore(str1,'_');

    %Whether its yearly,monthly or daily output locs and data
    switch strb
        case "yearly"
            a = reshape(a,[],no_grid);
            locs = a(1:2,:);
            dat = a(3:end,:);   
        case "monthly"
            a = reshape(a,[],no_grid);
            locs = a(1:2,:);
            sub = a(3:end,:);   
            sub = reshape(sub,14,size(sub,1)/14,[]); 
            dat = sub(2:13,:,:);
        case "daily"
            a = reshape(a,[],no_grid);
            locs = a(1:2,:);
            sub = a(3:end,:);   
            sub = reshape(sub,361,[],no_grid);
            sub(1,:,:)=[]; dat=sub;
    end
    
end


function makessdmap(locs,dat,var,map_d,dir_data)

    load coastlines   
    latlim = [-90 90];
    lonlim = [-180 180];

    %Find resolution of data
    sub = diff(locs(1,:));sub(sub==0)=[];res = mode(sub);
    %Initialize map
    mm = nan(180/res,360/res);    
    
    %Finds the map space of the locs
    locz(1,:) = (90-res/2-locs(1,:))/res+1;    
    locz(2,:) = (180+res/2+(locs(2,:)))/res;   
        
    %Assign data to image
    for i=1:size(locz,2) 
        mm(locz(1,i),locz(2,i))=dat(i);
    end
    mm(mm==0) = NaN;

    %Does mask
    mask=zeros(size(mm,1),size(mm,2));
    cd /data/sm1epk/crop_sets/coun_masks
    for c=1:246
        load(['mask_',num2str(c),'.mat'])
        s2r=imresize(s2r,[size(mm,1) size(mm,2)],'Nearest');
        mask=mask+s2r;
    end
    mask(mask>0)=1;

    %Does map
    cd(map_d)
    varr = strrep(var,'_',' ');
    p7print(mm,varr,var,dir_data,400,'BrBG4',0,0,[0 4],20,mask,[-90 90],[-180 180],'Countries')

end
