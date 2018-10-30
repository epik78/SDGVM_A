clear;dir1=pwd;



%INPUTS
%Directory where the SDGVM run is
dir_data = '/fastdata-sharc/sm1epk/SDGVM_runs/';
%Variable I want to read.It does the yearly/monthly/daily outputs
%Not the default outputs
var = 'daily_gpp';


%Goes in the directory where the SDGVM run is
cd(dir_data)

%Opens and reads the file
fid = fopen([var,'.dat'],'rt');
a = fscanf(fid,'%f');
fclose(fid);

locs = a(1:2,:);

%Gets the first year of the file
yearr=a(3);

ind(1)=4;coun=0;
ch=1;
while ch

        
    if(leapyear(yearr))
        ind(2)=ind(1)+365;
    else
        ind(2)=ind(1)+364;
    end
    
    coun=coun+1;
    dat{coun,1}=yearr;
    dat{coun,2}=a(ind(1):ind(2));
    ind(1)=ind(2)+1+1;
    yearr=yearr+1;
    
    if(ind(2)==size(a,1))
        ch=0;
    end

end

cd(dir1)