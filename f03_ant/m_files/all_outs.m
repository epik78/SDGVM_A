clear;

dirr{1} = pwd;
%Where the run is
dirr{2} = '/fastdata-sharc/sm1epk/SDGVM_runs/fin_outs3';
%Where the map print is
dirr{3} = '/data/sm1epk/EW/weath_emp_model/codes/iter_1';


fig_id = [2];

crid = [1];
    
countr = {'Argentina','Brazil','Canada','China','Ethiopia','France','Germany',...
  'India','Italy','Kenya','Mexico','Indonesia','Nigeria','Russian Federation',...
  'South Africa','Turkey','United States of America'};
crops = {'Maize','Soybeans','Millet','Wheat','Wheat','Barley','Barley','Sorghum'};
crops_sd = {'Maize','Soy','Millet','Swheat','Wwheat','Sbarley','Wbarley','Sorghum'};
crops_sage = {'maize','soybean','millet','wheat','wheat','barley','barley','sorghum'};


limms(1,1:8) = [13 5 3 8 8 7 7 8];

mask=zeros(180*6,360*6);
cd /data/sm1epk/crop_sets/coun_masks
for c=1:246
    load(['mask_',num2str(c),'.mat'])
    s2r=imresize(s2r,size(mask));
    mask=mask+s2r;
end
mask(mask>0)=1;



switch fig_id

case 1
  
    %Read the shapefile with 246 country borders
    cd /data/sm1epk/crop_sets/TM_WORLD_BORDERS-0.3
    s=shaperead('TM_WORLD_BORDERS-0.3.shp');
    %Correction for country names
    s(175).NAME = 'Russian Federation';
    s(209).NAME = 'United States of America';
    %Gets all country names in s
    countr_s = {s(:).NAME};
   
    %Reads FAO yields.fao_y(country,crop,1961-2017) in hg/ha
    cd(dirr{1}) 
    fao_y = readfao('/data/sm1epk/crop_sets/FAO_Data',countr,crops,{'Yield'});

    for ii = 1:size(crid,2)
        %Holds data of scatter plot
        sc_dat = zeros(size(countr,2),2);

        crops_sd{crid(ii)}

        cd(dirr{1}) 
        %Gets SDGVM yield map 
        [locs,sdgvm_y] = read_crop_outs([crops_sd{crid(ii)},'_ryield'],dirr{2});    

        for jj = 1:size(countr,2)
            
            %Find the index in s for each country Im looking at
            ind_s = find(strcmp(countr_s,countr{jj}));
            %Load mask for particular country
            cd /data/sm1epk/crop_sets/coun_masks
            load(['mask_',num2str(ind_s),'.mat'])
        
            %Mask resize SDGVM resolution
            s2r = imresize(s2r,[size(sdgvm_y,1) size(sdgvm_y,2)]);
            %Multiply by mask
            sdgvm_ym = s2r.*sdgvm_y;
            %Remove zeros and nans
            sdgvm_ym(sdgvm_ym==0)=NaN;sdgvm_ym(isnan(sdgvm_ym))=[];            
            %Converts from g/m2 SDGVM yield to tn/ha
            sdgvm_ym = sdgvm_ym/100;

            cd(dirr{1}) 
            %Gets FAO yield data for the country and crop as read above
            fao_yp = fao_y(jj,ii,end-17);
            %Given in hg/ha and converted here to tn/ha
            fao_yp = fao_yp*100/1000/1000;
            
            sc_dat(jj,1) = fao_yp;
            sc_dat(jj,2) = median(sdgvm_ym);
        end

        figure;hold;
        scatter(sc_dat(:,1),sc_dat(:,2),[],'d','filled','MarkerFaceColor',[.7 0 0],'MarkerEdgeColor',[0 0 0])
        text(sc_dat(:,1),sc_dat(:,2),countr)
         
        hold off
        xlabel('FAO Yield');ylabel('SDGVM Yield')
        title(crops{crid(ii)})
    end


case 2
    
    for ii = 1:size(crid,2)
        cd(dirr{1}) 

        %Gets SDGVM yield map 
        [locs,sdgvm_y] = read_crop_outs([crops_sd{crid(ii)},'_ryield'],dirr{2});    
        %Converts from g/m2 SDGVM yield to tn/ha
        sdgvm_y = sdgvm_y/100;
        mask_m = imresize(mask,size(sdgvm_y),'Nearest');
        cd(dirr{3}) 
        p7print(sdgvm_y,[crops_sd{crid(ii)},' Yield, SDGVM'],[crops_sd{crid(ii)}],'/fastdata-sharc/sm1epk/',100,'RdYlGn6',0,0,[0 limms(1,crid(ii))],20,mask_m,[-90 90],[-180 180],'Countries')

        %Gets SAGE yield map
        cd(dirr{1})
        sage_y = readsage(crops_sage{crid(ii)});
        sage_y(sage_y==0) = NaN;
        mask_m = imresize(mask,size(sage_y),'Nearest');
        cd(dirr{3}) 
        p7print(sage_y,[crops_sd{crid(ii)},' Yield, SAGE'],[crops_sd{crid(ii)}],'/fastdata-sharc/sm1epk/',100,'RdYlGn6',0,0,[0 limms(1,crid(ii))],20,mask_m,[-90 90],[-180 180],'Countries')
         
    end


end

cd(dirr{1})

