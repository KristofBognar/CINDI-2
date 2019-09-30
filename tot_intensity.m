% calculate total intensities fo HsS and AS

lampscan=true;

if lampscan

    % for lamp scan on Sept. 10 (CINDI-2)
    elevlist=0.2:-0.05:-0.8;
    elevlist=elevlist';


%     inarr=load('PC1H_2016_263.hs');
    inarr=load('PC1H_2016_254.hs');
    st=30;
    fin=50;

    flux=sum(inarr(st:fin,1:2048),2); % sum all counts across the CCD
    elev=inarr(st:fin,2049); 
    az=inarr(st:fin,2050); 
    ft=inarr(st:fin,2051); 
    
    
    plot(elev,flux, 'bx')

%     ind=find(az~=287 & az~=290 & elev<=-0.2 & elev>=-1.0);
%     flux=flux(ind);
%     elev=elev(ind);
%     az=az(ind);
%     ft=ft(ind);
% 
%     ii=size(elevlist,1);
%     azlist=[];
%     elevplot=[];
%     fluxplot=[];
% 
%     for i=1:size(elev,1)-ii+1
%     %     round(elev(i:ii+i-1)-elevlist,4)
%         if ~any(round(elev(i:ii+i-1)-elevlist,4))
%             elevplot=[elevplot, elev(i:ii+i-1)];
%             fluxplot=[fluxplot, flux(i:ii+i-1)];
%             azlist=[azlist,az(i+1)];
%         end
% 
%     end
% 
%     fluxplot=normc(fluxplot);
% 
%     legendcell=cell(1,size(azlist,2));
%     for i=1:size(azlist,2)
%         plot(elevplot(:,i), fluxplot(:,i)), hold on
%         legendcell{i}=num2str(azlist(i));
%     end
% 
%     legend(legendcell,'location','northwest')
%     xlabel('Elevation angle')
%     ylabel('Normalized total flux')
%     
    % azind=find(az==288.6);
    % az=az(azind);
    % elev=elev(azind);
    % flux=flux(azind)
    
%     figure
%     for i=1:size(fluxplot,1)
%         plot(azlist(2:4),fluxplot(i,2:4)), hold on
%     end

else

    % calculate flux difference between each level for horizon scan
    
    inarr=load('PC1H_2016_273.hs');

    flux=sum(inarr(:,1:2048),2); % sum all counts across the CCD
    elev=inarr(:,2049); 
    az=inarr(:,2050); 
    ft=inarr(:,2051); 
    
    
    delta_flux=diff(flux);
    [maxval,ind]=max(delta_flux(1:end-1));
    horizon_elev=elev(ind:ind+1)
%     figure(1)
%     plot(elev(1:end-2), delta_flux(1:end-1), 'bx');
%     xlabel('Elevation angle')
%     ylabel('\Delta flux')
    
    % for curve fitting
    
%     ind=find(elev <= 2 & elev >= -2);
    ind=1:27; % for one scan only...
    elev2_temp=elev(ind);
    delta2=diff(flux(ind));
    
    elev2=elev2_temp(1:end-1)+diff(elev2_temp)./2;

    [f,g]=fit(elev2,delta2,'gauss1')
    
    figure
    plot(elev(ind),flux(ind),'bx')
    xlabel('Elevation angle')
    ylabel('Flux')
    
    
    figure
    plot(f,elev2,delta2)
    xlabel('Elevation angle')
    ylabel('\Delta flux')
    
%     clearvars

    
end

% % To see how QDOAS DSCDs compare with different settings
% % first group: no differential XS
% % second group: no2 and no2a differential XS
% % third group: no2 diffXS, no2a orthogonalized to no2
% file_nm= 'diff_vs_not.ASC';
% fid = fopen(file_nm, 'r');
% fgetl(fid);
% qdoas_raw = (fscanf(fid,'%f', [60,inf]))';
% fclose(fid);
% % 
% % plot(qdoas_raw(:,4),qdoas_raw(:,11),'kx'), hold on
% % plot(qdoas_raw(:,4),qdoas_raw(:,27),'ro'), hold on
% % plot(qdoas_raw(:,4),qdoas_raw(:,43),'b+'), hold on
% 
% arr=zeros(27,3);
% 
% arr(:,1)=qdoas_raw(:,11);
% arr(:,2)=qdoas_raw(:,27);
% arr(:,3)=qdoas_raw(:,43);
% 
% arrdiff=diff(arr,1,2);
% 
% plot(abs((arr(:,1)./arr(:,3))-1)*100, 'bx')
% 
