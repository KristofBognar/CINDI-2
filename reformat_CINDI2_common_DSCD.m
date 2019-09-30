% reformat median dSCDs provided by Michel to the format we use for profile
% retrievals

VIS=false;


if VIS
    load original.mat
    data_no2=NO2visMaxdoasMedianDSCDs;
    data_o4=O4visMaxdoasMedianDSCDs;
else
    load original_UV.mat
    data_no2=NO2uvMaxdoasMedianDSCDs;
    data_o4=O4uvMaxdoasMedianDSCDs;
end

    %% no2
    % create table
    no2=table;

    % save DOY (jan 1 00:00 = 0) and time
    no2.DOY=data_no2(:,1);
    no2.Fractionaltime=data_no2(:,2);

    % create DDMMYY and HHMMSS fields
    no2.DateDDMMYYYY=datestr(yeartime(2016)+no2.DOY,'dd/mm/yyyy');
    % no2.DateDDMMYYYY=datetime(yeartime(2016)+no2.DOY,'ConvertFrom','datenum');
    no2.Timehhmmss=datestr(yeartime(2016)+no2.DOY,'HH:MM:SS');

    % pointing info
    no2.SZA=data_no2(:,3);
    no2.SolarAzimuthAngle=data_no2(:,4);
    no2.Elevviewingangle=data_no2(:,5);
    no2.Azimviewingangle=data_no2(:,6);

    % dSCDs
    no2.NO2dSCD=data_no2(:,7)*1e15;
    no2.NO2dSCDerr=data_no2(:,8)*1e15;

    no2.O4dSCD=data_o4(:,7)*1e42;
    no2.O4dSCDerr=data_o4(:,8)*1e42;

    % remove unnecessary zenith measurements
    time=datetime(yeartime(2016)+no2.DOY,'ConvertFrom','datenum');
    ind=find(minute(time)<10 & hour(time)~=11 & hour(time)>5 & hour(time)<17);
    ind2=find(hour(time)==11);
    ind=unique(sort([ind;ind2]));

    no2=no2(ind,:);
    time=time(ind,:);

    % merge noon zenith measurements
    for i=[258,266]
        ind=find(hour(time)==11 & minute(time)>=30 & minute(time)<=40 &...
                 no2.Elevviewingangle==90 & floor(no2.DOY)==i);
        no2(ind,1:8)=repmat(no2(ind(6),1:8),length(ind),1);
        no2(ind(2:end),:)=[];
        time(ind(2:end))=[];

        no2.NO2dSCD(ind(1))=0;
        no2.NO2dSCDerr(ind(1))=0;
        no2.O4dSCD(ind(1))=0;
        no2.O4dSCDerr(ind(1))=0;

    end


    % % % remove noon zenith measurements
    % % ind=find(hour(time)==11 & no2.Elevviewingangle==90 & floor(no2.DOY)==258);
    % % no2(ind(3:end-1),:)=[];
    % % time(ind(3:end-1),:)=[];
    % % ind=find(hour(time)==11 & no2.Elevviewingangle==90 & floor(no2.DOY)==266);
    % % no2(ind(3:end-1),:)=[];
    % % time(ind(3:end-1),:)=[];

    % % % remove dSCD values for ZS measurements
    % % ind=find(no2.Elevviewingangle==90);
    % % for i=ind'
    % %     no2.NO2dSCD(i)=0; 
    % %     no2.NO2dSCDerr(i)=0; 
    % %     no2.O4dSCD(i)=0; 
    % %     no2.O4dSCDerr(i)=0; 
    % % end

    % wite to csv file
    ind=find(floor(no2.DOY)==258 & no2.Azimviewingangle==287);
    if VIS
        writetable(no2(ind,:),'commonDSCD_2016_09_15_v2.csv');
    else
        writetable(no2(ind,:),'commonDSCD_2016_09_15_UV.csv');
    end
    
    ind=find(floor(no2.DOY)==266 & no2.Azimviewingangle==287);
    if VIS 
        writetable(no2(ind,:),'commonDSCD_2016_09_23_v2.csv');
    else
        writetable(no2(ind,:),'commonDSCD_2016_09_23_UV.csv');
    end

