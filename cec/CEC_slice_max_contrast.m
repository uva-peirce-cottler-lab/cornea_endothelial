function [bw_for,im_elev] = CEC_slice_max_contrast(bw_xyz,xyz)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


    [bw_for, im_min] = zslice_map(xyz, @(x) CEC_zslice_thresh(x));
    
    % Find top surface
    
    % Recheck all z slices based on contract between 2
    
    % Elevation image
end

