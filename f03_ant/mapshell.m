clear;dir1=pwd;

%INPUTS
%Directory where the SDGVM run is
dir_data = '/fastdata-sharc/sm1epk/SDGVM_runs/';
%Variable I want to read.It does the yearly/monthly/daily outputs
%Not the default outputs
var = 'daily_gpp';


%Reads file and outputs locs and data
[locs,dat] = readsdfiles(var,dir_data);


cd(dir1)

%Figures whether its yearly, monthly or daily
strb = extractBefore(var,'_');

switch strb
    case "yearly"
        dat=dat(1,:);
    case "monthly"
        dat=squeeze(dat(6,1,:));
    case "daily"
        dat=squeeze(dat(1,1,:));
end

makessdmap(locs,dat,var)

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

    %Opens and reads the file
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


function makessdmap(locs,dat,var)

    load coastlines   
    latlim = [-90 90];
    lonlim = [-180 180];

    sub = diff(locs(1,:));sub(sub==0)=[];res = mode(sub);
    mm = nan(180/res,360/res);    

    locz(1,:) = (90-res/2-(-locs(1,:)))/2+1;    
    locz(2,:) = (180+res/2+(locs(2,:)))/2;     

    rasterSize=[size(mm,1) size(mm,2)];
    
    for i=1:size(locz,2) 
        mm(locz(1,i),locz(2,i))=dat(i);
    end
    mm=flipdim(mm,1);
    axesm('eckert4','Frame','on','MeridianLabel','On','ParallelLabel','On',...
      'PLabelMeridian','West','MLabelParallel',-60,'MlabelLocation',[-90 0 90],...
      'PlabelLocation',[-45 0 45],'Grid','On','MlineLocation',[-90 0 90],...
      'PLineLocation',[-45 0 45],'FontSize',5);
    R=georefcells(latlim,lonlim,rasterSize,'ColumnsStartFrom','north');
    plotm(coastlat,coastlon,'color','black','LineWidth',1)
    h=geoshow(mm,R,'DisplayType','Surface');
    
    var(regexp(var,'_'))=' ';
    title(var,'FontSize',12,'FontWeight','Bold')
end
