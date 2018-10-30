function [gs_for, min_xy, max_xy] = CEC_zslice_thresh(slice)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    I_eq = adapthisteq(slice);
    bslice = imfilter(I_eq,fspecial('disk',3),'same');
%    bslice = medfilt2(I_eq, [5 5],'symmetric');
%         ws_bslice = watershed(imcomplement(bslice))>0;
% keyboard
   
    min_xy = double(ordfilt2(bslice,1,ones(49,49),'symmetric'))+1;
    max_xy = double(ordfilt2(bslice,49.^2,ones(49,49),'symmetric'));
    
    gs_thresh = bslice>(max_xy-(min_xy-1))/2+(min_xy-1);
    gs_background = min_xy+35>max_xy;
    
    gs_init_for = (gs_thresh & ~gs_background);
    % | bwareaopen(edge(bslice,'sobel'),5);
%     imopen
      gs_for = imopen(gs_init_for,strel('disk',3,0));
      gs_for2 = imopen(gs_for,strel('disk',6,0));
% keyboard
    
    
    

% outputArg2 = inputArg2;
end

