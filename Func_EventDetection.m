
function [PeaksData, PercAct, Good_TC ]= Func_EventDetection(T_Cells_DF_nosub,t , mskNrn, im_stack, T_date, T_num)

            amp = 1; % min Amplitude of Peak in %
            prom = 2 ;% min Prominance of Peak in %
            minWid = .05; % min Width of peaks in sec 
            maxWid = 200; % min Width of peaks in sec 
            Col = [.5 .5 .5];
            PeaksData=[];
            m=0;n=1;j=1;Good_idx=[];
    for i=1:size(T_Cells_DF_nosub,1)
            CC = smoothdata(T_Cells_DF_nosub(i,:), 'movmean' , 20);
            [pks,locs,w,p] = findpeaks(CC , t);
%             [pks,locs,w,p] =findpeaks(CC, t,'MinPeakHeight',amp,'MinPeakProminence', prom , 'MaxPeakWidth' , maxWid , 'Annotate','extents');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
         Peaks = find(p>prom & w>minWid & w<maxWid) ; %pks>amp & 
         locs_par = locs(Peaks); 

         if ~isempty(locs_par)
            plot(t, CC , 'color', Col, 'LineWidth' , 2); hold on; 
%             title(num2str(i))
            for ii=1:length(locs_par), idx_par(ii) = find(t==locs_par(ii));end
            pks_par = pks(Peaks); 
            w_par = w(Peaks);
            p_par = p(Peaks);
            for iii=1:length(locs_par),  plot([locs_par(iii) locs_par(iii)], [-6 14] , 'r'); hold all;end            
            m = n + length(locs_par) ;
            PeaksData(n:m-1) = locs_par ;
            n = m ;
            if p_par > .5
            Good_idx(j) = i ; j=j+1;
            end
%          ginput(1)
         end
    end

Good_TC = T_Cells_DF_nosub(Good_idx,:);
Good_msk = mskNrn(Good_idx,:,:);
%

figure(11);clf;
set(gcf, 'position', [350 100 900 600]);
subplot(211);
    for i=1:size(Good_TC,1)
      C = smooth2d(  Good_TC(i,:) , 5 );
    plot(t, C ,'color' , Col,  'LineWidth' , 3 );hold on;
    end
    axis tight
    yLim=get(gca, 'ylim');
    for i=1:length(PeaksData)
    plot([PeaksData(i) PeaksData(i)] ,  yLim , 'r'); hold on;   
    end
xlim([0 t(end)]); box off;
ylabel('\DeltaF/F, %');
axis tight

if length(PeaksData)>=3
subplot(212);
Matt_MEA_nhist(PeaksData,  'minbins',  length(PeaksData), 'noerror', 'color' , 'b');
xlim([0 t(end)]);
end
%%
figure(22);clf;
set(gcf, 'position', [50 100 900 600]);
            imshow(im_stack); hold on;
            for i=1:size(mskNrn,1)
                ms = squeeze( mskNrn(i,:,:) );
                imcontour(ms , 'r');
            end
            for i=1:size(Good_msk,1)
                ms = squeeze( Good_msk(i,:,:) );
                imcontour(ms , 'g');
            end          
            PercAct = size(Good_msk,1)/size(mskNrn,1)*100;
title([num2str(size(Good_msk,1)) '/' num2str(size(mskNrn,1)) '... ' num2str(floor(PercAct)) '%' ]);
axis image  ; axis off; hold on;
%
%%
figure(33); clf;
Blou = [0 1 0.6];
set(gcf, 'position', [300 50 500 900]);
subplot(3,1,1:2); n=0;
    for i=1:size(Good_TC,1)
      C = smooth2d(  Good_TC(i,:) , 5 );
    plot(t, C+n ,'color' , Col,  'LineWidth' , 1);hold on;n=n+10;
    end
    axis tight
    yLim=get(gca, 'ylim');
    for i=1:length(PeaksData)
    plot([PeaksData(i) PeaksData(i)] ,  yLim , 'color' , Blou); hold on;   
    end
    n=0;
    for i=1:size(Good_TC,1)
      C = smooth2d(  Good_TC(i,:) , 5 );
    plot(t, C+n ,'color' , 'k',  'LineWidth' , 1);hold on;n=n+10;
    end
    box off
xlabel('time,s');
ylabel('\DeltaF/F, %');
title([T_date ' run' T_num ' - nROI=' num2str(size(Good_TC,1)) ]);
subplot(313);
if length(PeaksData)>1
Matt_MEA_nhist(PeaksData,  'binfactor', 20,  length(PeaksData), 'noerror', 'color' , Blou);hold on;
% Matt_MEA_nhist(PeaksData, 'binfactor',20);
xlim([0 t(end)]);
end