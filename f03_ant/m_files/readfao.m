function dat = readfao(dirr,countr,crops,var)

    %Reads FAO crop data
    %Goes to dir with FAO crop data
    cd(dirr)

    %Will hold data for each country,crop,variable and all the years
    dat = nan(size(countr,2),size(crops,2),size(var,2),2017-1961+1);

    %Reads file
    [num,txt,raw] = xlsread('Production_Crops_E_All_Data_NOFLAG.xlsx');

    %Holds countries column
    sub1 = raw(:,2);
    %Holds crops column
    sub2 = raw(:,4);
    %Holds variables column
    sub3 = raw(:,6);

    %Find the rows in the excel file with each country
    for ii = 1:size(dat,1)
        ind1 = strcmp(sub1,countr{ii});ind1=find(ind1);
        
        %Same with crop type
        for jj = 1:size(dat,2)
            ind2 = strcmp(sub2,crops{jj});ind2=find(ind2);
            %Intersection of each country and crop type
            Ca = intersect(ind1,ind2);

            %For each variable
            for kk = 1:size(dat,3)
                ind3 = strcmp(sub3,var{kk});ind3=find(ind3);
                %Interesection of all 3
                Cb = intersect(Ca,ind3);
                %Get data
                if(~isempty(Cb))
                    dat(ii,jj,kk,:) = cell2mat(raw(Cb,8:64));
                else
                    dat(ii,jj,kk,:) = NaN;
                end
            end  
        end
    end
end
