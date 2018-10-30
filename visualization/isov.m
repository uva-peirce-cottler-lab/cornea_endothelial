function isov(zimg)

rzimg = false([512 512 size(zimg,3)]);
for z=1:size(zimg,3)
    rzimg(:,:,z) = imresize(zimg(:,:,z),[512 512]);
end

isosurface(rzimg, 0.5);

end