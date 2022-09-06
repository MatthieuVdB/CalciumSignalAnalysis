
function [T_Cells_F_nosub, T_Cells_DF_nosub ]= Func_TimeCourse(msknum, t,tiff_stack_2,tiff_stack_1,mskNrn,framenumber)

T_Cells_DF_nosub= zeros(msknum-1 , length(t)) ; 
T_Cells_F_nosub= zeros(msknum-1 , length(t)) ; 
fprintf('\n');
str = ['......making TCs, progress...%% '] ;       
fprintf(str);                     
    if exist('tiff_stack_2')
        framenumber_1 = size(tiff_stack_1,3);
        framenumber_2 = size(tiff_stack_2,3);
                                 for kk= 1:msknum
                                    msk=find(mskNrn(kk,:,:) > 0);
                                    yMsk= zeros(1,length(t));
                                    dff_tmp = zeros(1,length(t));
                                    for ii=1:framenumber_1
                                        foo = tiff_stack_1(:,:,ii); 
                                        yMsk(ii)=nanmean(foo(msk)); 
                                    end
                                    for ii=framenumber_1+1:framenumber_1+framenumber_2
                                        foo = tiff_stack_2(:,:,ii-framenumber_1); 
                                        yMsk(ii)=nanmean(foo(msk)); 
                                    end
                                    T_Cells_F_nosub(kk,:)=yMsk;
                                    minF = nanmean(yMsk(2:length(t)/5)); %disp([num2str(minF)]);
                                    dff_tmp = ((yMsk-minF) / minF )*100;
                                    T_Cells_DF_nosub(kk,:)=dff_tmp;                                    
                                    if kk>1
                                        for j=0:log10(kk-1), fprintf('\b'); end % delete previous counter display
                                    end
                                    fprintf('%d', kk);
                                end  
    else
                                for kk=1:msknum
                                    msk=find(mskNrn(kk,:,:) > 0);
                                    yMsk= zeros(1,length(t));
                                    for ii=1:framenumber
                                        foo = tiff_stack_1(:,:,ii); 
                                        yMsk(ii)=mean(foo(msk)); 
                                    end
                                    T_Cells_F_nosub(kk,:)=yMsk;
                                    minF = nanmean(yMsk(20:framenumber/2)); %disp([num2str(minF)]);
                                    T_Cells_DF_nosub(kk,:)=((yMsk-minF) / minF )*100;
                                    if kk>1
                                        for j=0:log10(kk-1), fprintf('\b'); end % delete previous counter display
                                    end
                                      fprintf('%d', kk);
                                end
                                   
    end
 % VIZ
fprintf('...DONE');
                figure(6); clf;
                set(gcf, 'position', [400 50 1100 700]);
                subplot(211);
                plot(t, T_Cells_F_nosub);hold on; title('F');
                subplot(212);
                plot(t, T_Cells_DF_nosub);hold on;title('deltaF/F');                    
%          