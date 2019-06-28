function [locs,mm] = read_sd_outs(str1,dir_data)
    %Goes in the directory where the SDGVM run is
    cd(dir_data)
    
    a = dlmread([str1,'.dat']);

    %locs(lat/lon,grid cell)
    locs = a(:,1:2);
    %data(:,grid cell)
    dat = a(:,3:end);   
  
    %Find resolution of data
    sub = diff(locs(:,1));sub(sub==0)=[];res = mode(sub);
    
    %Initialize map
    mm = nan(180/res,360/res);    
    
    %Finds the map space of the locs
    locz(1,:) = (90-res/2-locs(:,1))/res+1;    
    locz(2,:) = (180+res/2+(locs(:,2)))/res;   
    
    %Assign data to image
    for i=1:size(locz,2)     
        mm(locz(1,i),locz(2,i))=dat(i,end);
    end

    mm(mm==0) = NaN;

end
