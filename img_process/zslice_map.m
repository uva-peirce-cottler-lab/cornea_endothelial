function out_arg_cell = zslice_map(img, fun)

% proc_img = zeros(size(img));
X=cell(size(img,3),nargout(fun));

% keyboard
for z=1:size(img,3)
    [X{z,:}] = fun(img(:,:,z));

end

out_arg_cell = cell(1,nargout(fun));

for c=1:nargout(fun)
    
   out_arg_cell{c} = cat(3, X{:,c}); 
end

%    proc_img(:,:,z) = 
end


