function proc_img = subimage_map(img, sdim,fun)


for r=1:size(img,1)/sdim(1)
   for c = 1:size(img,2)/sdim(2) 
       r_ind = sdim(1)*(r-1)+1:sdim(1)*r;
       c_ind = sdim(2)*(c-1)+1:sdim(2)*c;
       temp_img = fun(img(r_ind,c_ind,:));
       
        if r==1 && c==1
           if size(temp_img,3)==1
                proc_img = zeros(size(img,1),size(img,2));
           else
              proc_img = zeros(size(img,1),size(img,2),size(img,3)); 
           end
        end
        proc_img(r_ind,c_ind,:) = temp_img;
   end
end


end


