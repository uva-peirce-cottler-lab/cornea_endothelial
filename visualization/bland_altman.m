function bland_altman(y1,y2)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

plot((y1+y2./2), y1-y2,'bo','markersize',2);

std_error = std(y1-y2);
mean_error = mean(y1-y2);
keyboard
% plot(

end

