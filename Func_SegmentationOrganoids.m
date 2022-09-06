
function [mskNrn, msknum, MASK ]= Func_SegmentationOrganoids(im_stack)

figure(5);clf;
i_eq = im_stack ;
    level = 0.45 ;
    min_size = 10 ;
    max_size = 100 ;
    
%extract ROIs
%             if method == 1
                %WATERSHED
%                 set(handles.text1,'String', 'Measuring ROIS (watershed)..');
%                 pause(0.01);
%                 BW3 = ROIExtractionMain.watershedSegmentation(i_eq, min_size , level);
%             else
%                 %EXTEND CENTER PIXEL
%                 set(handles.text1,'String', 'Measuring ROIS (extend center)..');
%                 pause(0.01);
%                 BW3 = ROIExtractionMain.extendPixelSegmentation(i_eq);
%             end
%                 BW3 = imextendedmax(i_eq, 0.13);
                BW3 = imextendedmax(i_eq, 0.13);

            %segment into individual areas which get a label and border
            labeledImage = bwlabel(BW3, 8);
            
            %make color for plot
            coloredLabels = label2rgb (labeledImage, 'hsv', 'k', 'shuffle');
            
            %overlay on plot
            imshow(i_eq); hold on;
            ih = imshow( coloredLabels);
            set( ih, 'AlphaData', 0.4 );
 %%
            % Extract all areas from above
%             ROIareas = regionprops(labeledImage, i_eq, 'all');
            ROIareas = regionprops(labeledImage, BW3, 'all');
            
            %get borders of all ROIs
            boundaries = bwboundaries(BW3);
            
            %remove ROIS with certain characteristics
            %-too large
            i_roi = 1;
            counter = 1;
            
            clear boundaries_selected ROIareas_selected
            while i_roi < size(ROIareas, 1)+1
                if ROIareas(i_roi).Area < max_size 
                    if ROIareas(i_roi).Area > min_size
                    %copy only smaller
                    ROIareas_selected(counter,1) = ROIareas(i_roi);
                    boundaries_selected (counter,1) = boundaries(i_roi);
                    counter = counter + 1;
                    end
                end
                i_roi = i_roi+1;
            end
            
            %stop if no rois found
            if ~exist('ROIareas_selected')
                disp('STOPPED, No ROIs could be extracted, try different settings');
                mskNrn=0; msknum=0; 
                return
            end
            
            ROIareas = ROIareas_selected;
            boundaries = boundaries_selected;
            
            %count of rois
            numberOfRois = size(ROIareas, 1)
            %plot outlines of ROIS over original image
%             subaxis(1,4,3,'sh',0,'sv',0,'MR',0, 'ML', 0, 'MT', 0, 'MB', 0, 'PaddingTop', 0, 'PaddingBottom', 0, 'PaddingLeft', 0, 'PaddingRight', 0);
            imshow(i_eq, 'border', 'tight');
            
            RoiImages = gcf;
            axis image;
            hold on;
            mskNrn=[];msknum=1;msk_tmp=[];
            numberOfBoundaries = size(boundaries, 1);
            MASK = zeros(size(im_stack));
            for k = 1 : numberOfBoundaries
                thisBoundary = boundaries{k};
                plot(thisBoundary(:,2), thisBoundary(:,1), 'g', 'LineWidth', 1);
                            msk_tmp(:,:) = poly2mask(thisBoundary(:,2),thisBoundary(:,1),size(im_stack,1),size(im_stack,2));
                            msknum = msknum + 1; mskNrn(msknum,:,:) = msk_tmp;
                  index = (msk_tmp~=0);
                  MASK(index) = msk_tmp(index);
            end
            hold off;
            
            suptitle([' nROI = ' num2str(msknum) ]);

% figure(3); clf;
% set(gcf, 'position', [400 50 900 800]);
%             imshow(A); hold on;
%             for i=1:size(mskNrn,1)
%                 ms = squeeze( mskNrn(i,:,:) );
%                 imcontour(ms , 'g');
% %                 drawnow;title([num2str(i)]);pause(0.2);
%             end
% axis image  ; axis off; hold on;
%%
%                   index = (MASK~=0);
%                   MASK(index) = msk_tmp(index);
% 
% coloredLabels = label2rgb (MASK, 'parula', 'k', 'shuffle');
% figure(1);clf,
%             imshow(i_eq); hold on;
%             ih = imshow( coloredLabels);
%             set( ih, 'AlphaData', 0.4 );
% % imagesc(msk_tmp);
% % imagesc(MASK);
% axis image  ; axis off; hold on;



