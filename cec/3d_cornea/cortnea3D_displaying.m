if ~exist('fvB','var')
	file = uigetfile('*.mat','Pick the isosurfaceData file') ;
	load(file(1:end-4))
end
clearvars -except fvB fvM fvT
figure('units','normalized','position',[0.05,0.1,0.9,0.8])

colorBottom = [0,0.184,0.655] ;
% colorBottom = [0,0,1] ;
colorMiddle    = [0,160,230]/255 ;
% colorMiddle = [0,0,0.7] ;
colorTop    = [0,255,255]/255 ;
% colorTop    = [0,0,0.5] ;


zLayers = 8 ; % [1,2,3,4,5,6,7,8]

pB = gobjects(1,512) ;
pM = gobjects(1,512) ;
pT = gobjects(1,512) ;

for i = 1:(zLayers*64)%64%384%length(fv)
	pB(i) = patch(fvB(i),'Visible','off'); 
	pM(i) = patch(fvM(i),'Visible','off'); 
	pT(i) = patch(fvT(i),'Visible','off'); 

	pB(i).EdgeColor = 'none';
	pM(i).EdgeColor = 'none';
	pT(i).EdgeColor = 'none';

	pB(i).FaceColor = colorBottom ;
	pM(i).FaceColor = colorMiddle ;
	pT(i).FaceColor = colorTop ;

end

set(pB,'Visible','on')
set(pM,'Visible','on')
set(pT,'Visible','on')

daspect([1,1,1]);
view(3);
axis tight;
camlight;
lighting gouraud
