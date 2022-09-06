function [tiff_stack_1, tiff_stack_2, nx, ny, info, framenumber] = Func_Tiff_Stack_Organoids(filename)

                    framenumber=1;
                    info=func_ext_uli_sifread2(filename, framenumber); %getting info from sif file
                    Tifs=[filename(1:end-4) '.tif']; 
                    disp('Getting Tif images from Sif file');
                    tiff_info = imfinfo(Tifs); % return tiff structure, one element per image
                   
                        tiff_stack_1 = imread(Tifs, 1) ; % read in first image
                        tiff_stack_2 = [];
                        nx = size(tiff_stack_1, 1);ny = size(tiff_stack_1, 2);
                    framenumber=length(tiff_info);
                    str = ['......making Tiff stack - n = ' num2str(framenumber) '   progression........' ] ;
                    disp(str);  fprintf('\n');
                        if framenumber > 2000
                        tiff_stack_1 = zeros( nx, ny, 2500);
                        tiff_stack_2 = zeros( nx, ny, framenumber-2500);
                            for ii = 1 : 2500
                                tiff_stack_1(:,:, ii) = double(imread(Tifs, ii, 'info' , tiff_info));
%                                 disp([ num2str(ii) ' / ' num2str(framenumber) ]);  
                                II=round(ii/framenumber*100)  ;
                                if II>1 && II<10, fprintf('\b');fprintf('\b');fprintf('\b');
                                fprintf('%d %%', II);   % delete previous counter display
                                elseif II>=10, fprintf('\b');fprintf('\b'); fprintf('\b');fprintf('\b');
                                fprintf('%d %%', II);
                                end                           
                            end
                            for ii = 1 : framenumber-2500
                                tiff_stack_2(:,:, ii) = double(imread(Tifs, 2500+ii, 'info' , tiff_info));
%                                 disp([ num2str(2500+ii) ' / ' num2str(framenumber) ]);
                                II=round((ii+2500)/framenumber*100)  ;
                                if II>1 && II<10, fprintf('\b');fprintf('\b');fprintf('\b');
                                fprintf('%d %%', II);   % delete previous counter display
                                elseif II>=10, fprintf('\b');fprintf('\b'); fprintf('\b');fprintf('\b');
                                fprintf('%d %%', II);
                                end    
                            end
                            
                        else
                        tiff_stack_1 = zeros( nx, ny, framenumber);
                            for ii = 1: framenumber
                                tiff_stack_1(:,:, ii) =double(imread(Tifs, ii, 'info' , tiff_info));   
                                if ii>1,  for j=0:log10(ii-1), fprintf('\b'); end , end% delete previous counter display
                                fprintf('%d %%', ii);
                            end
                        end
                        
              
                            
end