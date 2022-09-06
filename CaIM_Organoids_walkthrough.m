close all
clear all

filename = 'D:\DATA\CaIM\Organoids_SCZ\2021\m_08\d_100821\10082021_sz10_d240_06.sif';

    trial_name= filename(47:end-4);
    T_date = filename(47:54);
    T_num = filename(end-5:end-4);
    disp([trial_name])
    
%% MAKING Stack of Tiff images extracted from the .sif file
[tiff_stack_1, tiff_stack_2, nx, ny, info, framenumber]= Func_Tiff_Stack_Organoids (filename);
if isempty (tiff_stack_2), clear tiff_stack_2 ; end
%%  GETTING some data used for the next steps2
if framenumber>1000    % to make sure to have the actual long recording                           
    frameRate = info.framerate;
    Period = 1/frameRate ;
    expTime = info.exposureTime ;
    t = linspace(0 , framenumber*Period  , framenumber);         
    framenumber_1 = size(tiff_stack_1,3); 

%% MAKING Reference Image
close all
    im_stack = max(tiff_stack_1,[],3);
    im_stack = (im_stack-min(im_stack(:)))/(max(im_stack(:))-min(im_stack(:))); 
    im_stack = imadjust(im_stack);
    im_stack = imadjust(im_stack);
figure(1);
imshow(im_stack, 'border', 'tight');  axis image  ; axis off;

%% EXTRACT ROI THROUGH SEGMENTATION
close all
fprintf('\n');
fprintf('......Segmentation........');
    [mskNrn, msknum] = Func_SegmentationOrganoids(im_stack);
    if msknum==0
        T_Cells_F_nosub = 0;
        T_Cells_DF_nosub= 0;
        t= 0;
        im_stack= 0;
        mskNrn= 0;
        msknum = 0;
        nx= 0;
        ny= 0;
        info= 0;
        framenumber  = 0;   
        return; 
    end
    
fprintf('DONE');
%% Time course From  selected ROIs
[T_Cells_F_nosub, T_Cells_DF_nosub ]= Func_TimeCourse(msknum, t,tiff_stack_2,tiff_stack_1,mskNrn);
           
%% FIGURES
if msknum~=0

figure(3); clf;
set(gcf, 'position', [400 50 900 800]);
            imshow(im_stack); hold on;
            for i=1:size(mskNrn,1)
                ms = squeeze( mskNrn(i,:,:) );
                imcontour(ms , 'g');
%                 drawnow;title([num2str(i)]);pause(0.2);
            end
axis image  ; axis off; hold on;

figure(4); clf;
set(gcf, 'position', [400 50 500 300]);
plot(t, T_Cells_DF_nosub);hold on;
xlabel('time,s');
ylabel('\DeltaF/F, %');
title([T_date ' ' T_num ' - nROI=' num2str(msknum) ]);

figure(5); clf;
set(gcf, 'position', [900 50 500 900]); n=0;
    for i=1:size(T_Cells_DF_nosub,1)
    plot(t, T_Cells_DF_nosub(i,:)+n);hold on;n=n+10;
    end
xlabel('time,s');
ylabel('\DeltaF/F, %');
title([T_date ' ' T_num ' - nROI=' num2str(msknum) ]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Detecction of Calcium Event and getting Percentage of Active ROIs ... Also making the figures from data
Func_EventDetection(T_Cells_DF_nosub,t , mskNrn, im_stack, T_date, T_num) ;

%% SAving the files

clearvars tiff_stack_1 tiff_stack_2 

mkdir([DATA_fldr d_fldr(dte).name '\MatFile\' ]);
save([DATA_fldr d_fldr(dte).name '\MatFile\' T_date '-' T_num  '.mat']);
mkdir([DATA_fldr d_fldr(dte).name '\Traces\' ]);
str = [DATA_fldr d_fldr(dte).name '\Traces\'  T_date '-' T_num   '_ROI.tif'];
saveas(figure(3), str);
str = [DATA_fldr d_fldr(dte).name '\Traces\'  T_date '-' T_num   '_TC.tif'];
saveas(figure(5), str);
str = [DATA_fldr d_fldr(dte).name '\Traces\'  T_date '-' T_num   '_ROI_active.tif'];
saveas(figure(22), str);
str = [DATA_fldr d_fldr(dte).name '\Traces\'  T_date '-' T_num   '_TC_active_a.tif'];
saveas(figure(11), str);
str = [DATA_fldr d_fldr(dte).name '\Traces\'  T_date '-' T_num   '_TC_active_b.tif'];
saveas(figure(33), str);

fprintf('\n'); 
disp([DATA_fldr d_fldr(dte).name '\' TIF_files(counter).name '... DONE......']);

end
end


return



%% VISUALIZE STACK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [IM, rect] = imcrop(im_ref);
% IM2 = im_ref(rect(1):rect(4), rect(2):rect(3));
% imshow(IM2);
            Sfld = [DATA_fldr d_fldr(dte).name '\vid\' trial_name '\'] ;
            mkdir(Sfld)
figure(1);
% if ssimval<0.1
crange = [min(tiff_stack_1(:)) max(tiff_stack_1(:))*.9];
figure(1);clf;
    for i=1:10:length(tiff_stack_1)
        im = tiff_stack_1(:,:,i);
%         i= imcrop(im, rect);
%         im = smooth2d(im, 5);
        imagesc(im, crange); %
        colormap jet; axis image  ; axis off; colorbar; %set(gca,'visible','off')
%         title([num2str(i) '/' num2str(length(tiff_stack_1))]);        
                    sec =   round(((i)*Period)-Period) ;
                    minutes = 0 ;
                    if sec > 59 , minutes = floor(sec/60) ; sec = round(sec-(minutes*60)) ; end
                    if sec/10 < 1 , sec = ['0' num2str(sec)]; else sec = num2str(sec); end
                    if minutes/10 < 1 , minutes = ['0' num2str(minutes)]; else minutes = num2str(minutes); end
                    timestmp =[ minutes 'min ' sec 'sec'  ];
                
                text(10,495, timestmp , 'color' , 'y', 'FontWeight' , 'bold' , 'FontSize' , 15);
                drawnow;
            saveas(figure(1), [Sfld '\frame_' num2str(i) '.tif']);
    end
 
    if exist ('tiff_stack_2')
    for ii=1:10:length(tiff_stack_2)-100
        im = tiff_stack_2(:,:,ii);
%         im = imcrop(im, rect);
        imagesc(im, crange); 
        colormap jet; axis image  ; axis off; colorbar        
                    sec =   round(((i+ii)*Period)-Period) ;
                    if sec > 59 , minutes = floor(sec/60) ; sec = round(sec-(minutes*60)) ; end
                    if sec/10 < 1 , sec = ['0' num2str(sec)]; else sec = num2str(sec); end
                    if minutes/10 < 1 , minutes = ['0' num2str(minutes)]; else minutes = num2str(minutes); end
                    timestmp =[ minutes 'min ' sec 'sec'  ];                
                text(10,495, timestmp , 'color' , 'y', 'FontWeight' , 'bold', 'FontSize' , 15);
%                 title([num2str(i) ' / ' num2str(framenumber)]);    
%                 h = suptitle(trial_name);
%                 h.Interpreter  = 'none';
                drawnow;
            saveas(figure(1), [Sfld '\frame_' num2str(i+ii) '.tif']);
    end
    end
% end



