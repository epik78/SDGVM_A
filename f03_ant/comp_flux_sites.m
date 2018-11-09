clear;
dirr{1}=pwd;

site_num=[1:71,73:114];
site_num=[31 109 46 17 57 52];
str{1}='site_';
str{2}='outs';

flx_i=[16 21];
flx_name_sd={'gpp';'npp'};
cd /data/sm1epk/fluxnet/
load('sites_data_all.mat')

for ii = site_num
    ii
    figure
    A = sites_data{ii};
    
    %Gets years of fluxnet
    yy = A{5,2};

    
    cd(['/fastdata-sharc/sm1epk/',str{1},num2str(ii),'/',str{2}])
  
    %For each flux Im checking
    for jj=1:1%size(flx_i,2)
        %Opens and reads SDGVM flux
        fid = fopen(['daily_',flx_name_sd{jj},'.dat'],'rt');
        a = fscanf(fid,'%f');
        fclose(fid);
        
        %Gets lat/lon
        locs = a(1:2,:);

        %Gets the first year of the run
        yearr = a(3);

        ind(1) = 4;coun = 0;
        ch = 1;
        while ch
            %Checks whether leap year and assigns the end index        
            if(leapyear(yearr))
                ind(2) = ind(1)+365;
            else
                ind(2) = ind(1)+364;
            end
            
            %Adds a year index
            coun = coun+1;
            %Save year
            sd_years(coun) = yearr;
            %Save flux
            sd_run{coun,jj} = a(ind(1):ind(2));
            %Adjust index for next year
            ind(1) = ind(2)+1+1;
            %Sets next year to be read
            yearr = yearr+1;
            %If it reaches the end of the file stop
            if(ind(2)==size(a,1))
                ch = 0;
            end
        end
    end

    for jj = 1:size(flx_i,2)    
        %Gets carbon flux GPP
        fl_flux = A{flx_i(jj),2};
    
        coun = 0;
        %sub_y=cell2mat(dat(:,1));
        [y_ind,sub] = find(sd_years==yy');
        [C,ia,ib] = intersect(sd_years,yy);
        flux_sd=[];flux_fl=[];
        for jj = 1:size(ia,1)
            flux_sd=[flux_sd;sd_run{ia(jj)}];
            flux_fl=[flux_fl;fl_flux{ib(jj)}];
        end

        plot(flux_fl);hold;plot(flux_sd)
        %hold
        %plot(reshape(cell2mat(flux_sd),[],1))
        hold off;


    end
%
%
%
%
%
end
%
%
%
%
%
%
%
%
%
%
cd(dirr{1})

