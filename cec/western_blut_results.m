

f2_myh11_norm_by_viniculin = [0.5683537 0.186279171;...
    0.23616734 0.235015603;...
    0.277432854 0.314049816];

f2_cd31_norm_by_viniculin=[2.379404878 0.031601968;3.431349307 0.063448445;1.808083407 0.012809718];

f7_1a_hsp90_norm_viniculin = [0.69564845 0.233310859;0.728569211 0.016485812;...
0.845500337 0.18357265;0.860928932 0.305148261;...
0.341430991 0.417587567;0.1071465 0.436071723];


f7_2_asma_norm_viniculin = [2.259249694 1.896825594;2.142224833 0.193872912;...
3.314895346 2.540249974;2.92800762 2.691814914;2.722224422 3.276322581;...
    1.935349279 2.190928193];

plot_2sample(f2_myh11_norm_by_viniculin(:,1),f2_myh11_norm_by_viniculin(:,2),...
    'Myh11 Expression',{'Sclera','Cornea'});
[h,p]=ttest2(f2_myh11_norm_by_viniculin(:,1),f2_myh11_norm_by_viniculin(:,2))
set(gcf,'position',[200 200 120 125])

plot_2sample(f2_cd31_norm_by_viniculin(:,1),f2_cd31_norm_by_viniculin(:,2),...
    'CD31 Expression',{'Sclera','Cornea'});
[h,p]=ttest2(f2_cd31_norm_by_viniculin(:,1),f2_cd31_norm_by_viniculin(:,2))
set(gcf,'position',[200 200 120 125])

plot_2sample(f7_1a_hsp90_norm_viniculin(:,1),f7_1a_hsp90_norm_viniculin(:,2),...
    'HSP90 Expression',{sprintf('Ctrl'),sprintf('CEC')});
set(gcf,'position',[ 680   652   170   150])
[h,p]=ttest2(f7_1a_hsp90_norm_viniculin(:,1),f7_1a_hsp90_norm_viniculin(:,2))

plot_2sample(f7_2_asma_norm_viniculin(:,1),f7_2_asma_norm_viniculin(:,2),...
    'aSMA Expression',{'Ctrl','CEC'});
set(gcf,'position',[ 680   652   170   150])
[h,p]=ttest(f7_2_asma_norm_viniculin(:,1),f7_2_asma_norm_viniculin(:,2))


