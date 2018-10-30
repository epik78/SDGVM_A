clear;
dirr{1}=pwd;

site_num=4;

str{1}='site_';
str{2}='outs';

cd /data/sm1epk/fluxnet/
load('sites_data_all.mat')

for ii = site_num

    A = sites_data{ii};
    %Gets years of fluxnet
    yy = A{5,2};
    %Gets cover of fluxnet
    cov_fl = A{4,2};
    %Gets carbon flux GPP
    fl_fl = A{21,2};

    
    cd(['/fastdata-sharc/sm1epk/',str{1},num2str(ii),'/',str{2}])
  
    %Opens and reads the file
    fid = fopen(['daily_gpp.dat'],'rt');
    a = fscanf(fid,'%f');
    fclose(fid);

    locs = a(1:2,:);

    %Gets the first year of the file
    yearr = a(3);

    ind(1) = 4;coun = 0;
    ch = 1;
    while ch

        
        if(leapyear(yearr))
            ind(2) = ind(1)+365;
        else
            ind(2) = ind(1)+364;
        end
    
        coun = coun+1;
        dat{coun,1} = yearr;
        dat{coun,2} = a(ind(1):ind(2));
        ind(1) = ind(2)+1+1;
        yearr = yearr+1;
    
        if(ind(2)==size(a,1))
            ch = 0;
        end


    end



    coun=0;
    sub_y=cell2mat(dat(:,1));
    [y_ind,sub]=find(sub_y==yy);
    for jj=1:size(y_ind,1)
        flux_sd{jj}=dat{y_ind(jj),2};
        
        flux_fl{jj}=fl_fl{jj};


        
    end






end










cd(dirr{1})

