function bland_altman(y1,y2,x_text, y_text)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

plot((y1+y2./2), y1-y2,'bo','markersize',2);
hold on
plot([xlim' xlim'], repmat([mean(y)-1.96*std(y1-y21) ...
    mean(y)+1.96*std(y1-y2)],[2 1]),'g');
plot(xlim, [mean(y) mean(y)],'b'); hold off

xlabel(['Mean of ' x_text])
ylabel(['\delta of ' y_text])

end

