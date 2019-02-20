function emu_putog(site_num)
    
    dirr{1} = pwd;
    %Directory where the runs are
    dirr{2} = '/fastdata-sharc/sm1epk/';
    
    %Folder name where each site run will be saved
    str{1} = 'site_';
    %Folder where the SDGVM outputs will be saved
    str{2} = '/outs';
    %Folder where the data are collected for emailing
    str{3} = 'phils';

    system(['rm -Rf ',dirr{2},str{3}]);
    %Make directory to gather the data needed for emulator    
    system(['mkdir ',dirr{2},str{3}]);

    file_name{1} = 'daily_gpp';file_name{2} = 'daily_nep';

    %Loads fluxnet data for all sites
    cd /data/sm1epk/fluxnet/    
    load('sites_data_all.mat')
     
    %For each site
    for ii = site_num
        
        system(['mkdir ',dirr{2},str{3},'/',str{1},num2str(ii)]);
        clear obs_flux

        %Reads fluxnet data
        sub = sites_data{ii};
        %Observation years
        obs_y = sub{5,2};
        %obs_flux holds the daily NEE and GPP repsectively
        %for all the years that I have data
        subb = sub{16,2};obs_flux(:,1) = cell2mat(subb(:));
        subb = sub{21,2};obs_flux(:,2) = cell2mat(subb(:));  

        %Prints the observations in file
        cd([dirr{2},str{3},'/',str{1},num2str(ii)])
        fid=fopen('nep_obs.dat','w');
        fprintf(fid,'%8.4f\n',-obs_flux(:,1));
        fclose(fid);
        fid=fopen('gpp_obs.dat','w');
        fprintf(fid,'%8.4f\n',obs_flux(:,2));
        fclose(fid);

        cd([dirr{2},str{1},num2str(ii)])
        sub_dir = dir('run_*');
        %Finds the number of runs for each site
        noruns = size(sub_dir,1);
        

        %For each run
        for jj = 1:noruns

            system(['mkdir ',dirr{2},str{3},'/',str{1},num2str(ii),'/run_',num2str(jj)]);


            %Go into the folder
            %For each file I want to save
            for kk = 1:size(file_name,2)
                cd([dirr{2},str{1},num2str(ii),'/run_',num2str(jj),str{2}])
                %Reads the data for the specific flux
                dat = daily_reads(pwd,file_name{kk});
                mod_y=[dat{:,1}];
                
                %ib holds the indexes of the model years that we also have observations
                [C,ia,ib]=intersect(obs_y,mod_y);                
                dat = dat(ib,:); 

                %Prints the data in file
                cd([dirr{2},str{1},num2str(ii),'/run_',num2str(jj)])
                dat = dat(:,2);dat = cell2mat(dat);
                fid=fopen([file_name{kk},'_n.dat'],'w');
                fprintf(fid,'%8.4f\n',dat);
                fclose(fid);

                system(['cp ./',[file_name{kk},'_n.dat'],' ',[dirr{2},str{3},'/',str{1},num2str(ii),'/run_',num2str(jj),'/',...
                  [file_name{kk},'_n.dat']]]);
                
                system(['cp ./param_values.dat ',[dirr{2},str{3},'/',str{1},num2str(ii),'/run_',num2str(jj)]]);
            end %file name

        end %run

    end %site

    cd(dirr{1})

end

function dat=daily_reads(dir1,file_name)

    
    %Goes in the directory where the SDGVM run is
    cd(dir1)

    %Opens and reads the file
    fid = fopen([file_name,'.dat'],'rt');
    [file_name,'.dat']
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
   
end
