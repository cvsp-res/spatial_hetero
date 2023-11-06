

%% 
%%%% SECTION 0 : HEADING AND INITIATION
% MATLAB file name :  Qunatification Analysis for C1-C1(R) = 23 ; C2-C2(Y) = 92 ; C1-C2 = 92
% The parameter set details --> CM : [ GR , RD , C-Co1 , C-C , C-L , MMP ] = [ 2 , .75 , 46 , 46 , 4 , .01 ]

close all; 
clear all; 
clc;

% Defining the outputs of the model.
Number_of_disconnected_cell_clusters = [];          % Matrix for 1 st o/p
Area_of_largest_cell_cluster = [];                  % Matrix for 2 nd o/p
Total_area = [];                                    % Matrix for 3 rd o/p
Area_of_min_enclosing_cir = [];                     % Matrix for 4rth o/p          
Radius_of_min_enclosing_cir = [];                   % Matrix for 5 th o/p 


% Defining the variables and parameters
n = 50 ;                 % Number of replicates
mcs = 7;                 % { 1:0, 2:200, 3:400, 4:600, 5 : 800, ...  mcs:(mcs+1)*200 }  --> ( key : value )
format shortG




%%
%%% SECTION 1 : Image analysis and output calculation
 
% 1.1 : Heterogeneous case with various spatial configurations
ds1 = imageDatastore('D:\cvsp\1) Mar_22\1) C-C\3.3) R_halved_Y_doubled_92\2.1) Output_50','IncludeSubfolders',1,'LabelSource','foldernames');
% Give the path of folder that contains image screenshots to do image processing and quantification analysis
ds1.Files;                           % Datastore helps to store the file names and their paths at one place so that we can call them for analysis.
IM1 = reshape(ds1.Files,[9,n+2,6]);  % Re-shaping the Image matrix to access the paths of images as required.

% Image Analysis of the heterogeneous cases
for i = 1 : 6                          % i --> index running from 1 to 6 as we analyse 6 spatial arrangements including one swap case.
    for j = 1 : n
        I = imread(cell2mat(IM1(mcs,j,i)));                     % reading the images one by one
        %figure, imshow(I);
        [r,g,b] = imsplit(I);                                   % obtaining the content of red, gree, blue in the image
        r = im2double(r);                                   
        g = im2double(g);
        b = im2double(b);
        a1 = ( r == 1 | (r+g) == 2 ) & ((r+g+b) ~= 3);          % This logic expression basically act as mask / filter
        %figure, imshow(a1);
        a = regionprops(a1,'area');
        area_1 = cat(1,a.Area);
        S = sum(area_1);
        c = regionprops(a1,'centroid');
        e = regionprops(a1,'Extrema');
        centroids = cat(1,c.Centroid);
        extrema = cat(1,e.Extrema);
        %imshow(a1);
        pointsX = extrema(:,1);   % x - coordinates of extrema
        pointsY = extrema(:,2);   % y - coordinates of extrema
    
        [iX, iY, NoOfColorChannel] = size(a1);         % optimization to find the centre and fit the min. enclosing cir.
        textPosition = [10 10];
        P0 = [iX/2 iY/2];
        [P, radiusMin] = fminsearch(@(P) radiusFromPoint(P,pointsX,pointsY), P0);
        [xx,yy]= meshgrid(1:iX, 1:iY);
        mask = hypot(P(1)-xx, P(2)-yy) <=radiusMin ;
        %figure, imshow(mask);
        circleArea = length(nonzeros(mask));
        Area_of_min_enclosing_cir(i,j) = circleArea;
        Radius_of_min_enclosing_cir(i,j) = radiusMin ;

        %surface
        %Image morphological operations
        se = strel('diamond',1.0);
        s1 = imdilate(a1,se);
        s2 = imfill(s1,'holes');  % Used to calculate Area of largest cell cluster
        %figure, imshow(s2);
        s22 = bwareaopen(s2,20);  
        %figure, imshow(s22);
        s3 = bwperim(s2,8);       
        %figure, imshow(s3);
        s4 = bwpropfilt(s3,'perimeter',1);
        %figure, imshow(s4)
        e2= regionprops(s4,'Extrema');
        extrema2= cat(1,e2.Extrema);
        %imshow(s3);  
        T(i) = sum(s3(:));      
        mmfill = imfill(s4,'holes');
        mmarea= length(nonzeros(mmfill));
        Tmm(i) = mmarea; %perim of mainmass
        mmoctagon= roipoly(s4,extrema2(:,1),extrema2(:,2));
        octarea= length(nonzeros(mmoctagon));
        MIrat(i)= (mmarea)/octarea;
        [labeled_image, N] = bwlabel(s3);

        Number_of_disconnected_cell_clusters(i,j) = N ;
        A = regionprops('table',s22,'area');                %%% s2 or s22 should work. 
        AREA = cat(1,A.Area);
        Area_of_largest_cell_cluster(i,j) = max(AREA);
        Total_area(i,j) = sum(AREA);


        RGB = I;
        drawCircle = insertShape(RGB,'circle',[P(1) P(2) radiusMin],'LineWidth',2);
        temp_img = imoverlay(drawCircle,s3,'white');
        text_str = ['radius:' num2str(radiusMin)  '; circle area:' num2str(circleArea) '; cell area:' num2str(Total_area(i,j))  '; dconn_obj:' num2str(Number_of_disconnected_cell_clusters(i,j))];
        finalImg = insertText(temp_img,textPosition,text_str,'FontSize',10,'BoxColor','white','BoxOpacity',0.2,'TextColor','white');
        %figure, imshow(finalImg);

    end
        
end


% 1.2 : Homogeneous case --> First extreme case 
ds2 = imageDatastore('D:\cvsp\1) Mar_22\1) C-C\4) All Identical\2) cc3d_output\cancol2_PS_C_C_23','IncludeSubfolders',1,'LabelSource','foldernames');
ds2.Files;
IM2 = reshape(ds2.Files,[9,n+1]);
for j = 1 : n 
    I = imread(cell2mat(IM2(mcs,j)));                           % reading the images one by one
    %figure, imshow(I);
    [r,g,b] = imsplit(I);                                       % obtaining the content of red, gree, blue in the image
    r = im2double(r);                                           % This logic expression basically act as mask / filter
    g = im2double(g);
    b = im2double(b);
    a1 = (r-g) > .4 & (r-b) > .4;
    %figure, imshow(a1);
    a = regionprops(a1,'area');
    area_1 = cat(1,a.Area);
    S = sum(area_1);
    c = regionprops(a1,'centroid');
    e = regionprops(a1,'Extrema');
    centroids = cat(1,c.Centroid);
    extrema = cat(1,e.Extrema);
    %imshow(image);
    pointsX = extrema(:,1);   % x - coordinates of extrema
    pointsY = extrema(:,2);   % y - coordinates of extrema
    
    [iX, iY, NoOfColorChannel] = size(a1);             % optimization to find the centre and fit the min. enclosing cir.
    textPosition = [10 10];
    P0 = [iX/2 iY/2];
    [P, radiusMin] = fminsearch(@(P) radiusFromPoint(P,pointsX,pointsY), P0);
    [xx,yy]= meshgrid(1:iX, 1:iY);
    mask = hypot(P(1)-xx, P(2)-yy) <=radiusMin ;
    %figure, imshow(mask);
    circleArea = length(nonzeros(mask));
    Area_of_min_enclosing_cir(7,j) = circleArea;
    Radius_of_min_enclosing_cir(7,j) = radiusMin ;

    %surface
    se = strel('diamond',1.0);
    s1 = imdilate(a1,se);
    s2 = imfill(s1,'holes');  % Used to calculate Area of Largest Cluster
    %figure, imshow(s2);
    s22 = bwareaopen(s2,20);  
    %figure, imshow(s22);
    s3 = bwperim(s2,8);       
    %figure, imshow(s3);
    s4 = bwpropfilt(s3,'perimeter',1);
    e2= regionprops(s4,'Extrema');
    extrema2= cat(1,e2.Extrema);
    %imshow(s3);  
    T(i) = sum(s3(:));      
    mmfill = imfill(s4,'holes');
    mmarea= length(nonzeros(mmfill));
    Tmm(i) = mmarea; %perim of mainmass
    mmoctagon= roipoly(s4,extrema2(:,1),extrema2(:,2));
    octarea= length(nonzeros(mmoctagon));
    MIrat(i)= (mmarea)/octarea;
    [labeled_image, N] = bwlabel(s3);
    Number_of_disconnected_cell_clusters(7,j) = N ;
    A = regionprops('table',s22,'area');        
    AREA = cat(1,A.Area);
    Area_of_largest_cell_cluster(7,j) = max(AREA);
    Total_area(7,j) = sum(AREA);

    RGB = I;
    drawCircle = insertShape(RGB,'circle',[P(1) P(2) radiusMin],'LineWidth',2);
    temp_img = imoverlay(drawCircle,s3,'white');
    text_str = ['radius:' num2str(radiusMin)  '; circle area:' num2str(circleArea) '; cell area:' num2str(Total_area(7,j))  '; dconn_obj:' num2str(Number_of_disconnected_cell_clusters(i,j))];
    finalImg = insertText(temp_img,textPosition,text_str,'FontSize',10,'BoxColor','white','BoxOpacity',0.2,'TextColor','white');
    %figure, imshow(finalImg);
    
end


% 1.3 : Homogeneous case --> Second extreme case 
ds3 = imageDatastore('D:\cvsp\1) Mar_22\1) C-C\4) All Identical\2) cc3d_output\cancol2_PS_C_C_92','IncludeSubfolders',1,'LabelSource','foldernames');
ds3.Files;
IM3 = reshape(ds3.Files,[9,n+1]);
for j = 1 : n 
    I = imread(cell2mat(IM3(mcs,j)));                               % reading the images one by one
    %figure, imshow(I);
    [r,g,b] = imsplit(I);                                           % obtaining the content of red, gree, blue in the image
    r = im2double(r);
    g = im2double(g);
    b = im2double(b);
    a1 = (r-g) > .4 & (r-b) > .4;                                   % This logic expression basically act as mask / filter
    %figure, imshow(a1);
    a = regionprops(a1,'area');
    area_1 = cat(1,a.Area);
    S = sum(area_1);
    c = regionprops(a1,'centroid');
    e = regionprops(a1,'Extrema');
    centroids = cat(1,c.Centroid);
    extrema = cat(1,e.Extrema);
    %imshow(image);
    pointsX = extrema(:,1);   % x - coordinates of extrema
    pointsY = extrema(:,2);   % y - coordinates of extrema
    
    [iX, iY, NoOfColorChannel] = size(a1);             % optimization to find the centre and fit the min. enclosing cir.
    textPosition = [10 10];
    P0 = [iX/2 iY/2];
    [P, radiusMin] = fminsearch(@(P) radiusFromPoint(P,pointsX,pointsY), P0);
    [xx,yy]= meshgrid(1:iX, 1:iY);
    mask = hypot(P(1)-xx, P(2)-yy) <=radiusMin ;
    %figure, imshow(mask);
    circleArea = length(nonzeros(mask));
    Area_of_min_enclosing_cir(8,j) = circleArea;
    Radius_of_min_enclosing_cir(8,j) = radiusMin ;

    %surface
    se = strel('diamond',1.0);
    s1 = imdilate(a1,se);
    s2 = imfill(s1,'holes');  % Used to calculate Area of largest cell cluster
    %figure, imshow(s2);
    s22 = bwareaopen(s2,20);  
    %figure, imshow(s22);
    s3 = bwperim(s2,8);      
    %figure, imshow(s3);
    s4 = bwpropfilt(s3,'perimeter',1);
    e2= regionprops(s4,'Extrema');
    extrema2= cat(1,e2.Extrema);
    %imshow(s3);  
    T(i) = sum(s3(:));     
    mmfill = imfill(s4,'holes');
    mmarea= length(nonzeros(mmfill));
    Tmm(i) = mmarea; %perim of mainmass
    mmoctagon= roipoly(s4,extrema2(:,1),extrema2(:,2));
    octarea= length(nonzeros(mmoctagon));
    MIrat(i)= (mmarea)/octarea;
    [labeled_image, N] = bwlabel(s3);
    Number_of_disconnected_cell_clusters(8,j) = N ;
    A = regionprops('table',s2,'area');        
    AREA = cat(1,A.Area);
    Area_of_largest_cell_cluster(8,j) = max(AREA);
    Total_area(8,j) = sum(AREA) ;


    RGB = I;
    drawCircle = insertShape(RGB,'circle',[P(1) P(2) radiusMin],'LineWidth',2);
    temp_img = imoverlay(drawCircle,s3,'white');
    text_str = ['radius:' num2str(radiusMin)  '; circle area:' num2str(circleArea) '; cell area:' num2str(Total_area(8,j))  '; dconn_obj:' num2str(Number_of_disconnected_cell_clusters(i,j))];
    finalImg = insertText(temp_img,textPosition,text_str,'FontSize',10,'BoxColor','white','BoxOpacity',0.2,'TextColor','white');
    %figure, imshow(finalImg);
end





%%
%%%% SECTION 2 : Outputs, Plots and Statistical Analysis

% Outputs
(Number_of_disconnected_cell_clusters);
(Area_of_largest_cell_cluster);
(Total_area);
(Area_of_min_enclosing_cir);
(Radius_of_min_enclosing_cir);


% Section 5 : Simplifying the representation in plots
NDC = Number_of_disconnected_cell_clusters;
ALC = Area_of_largest_cell_cluster;
TA  = Total_area ;
AMC = Area_of_min_enclosing_cir;
RMC = Radius_of_min_enclosing_cir ;


% Outputs you need to consider for performing the statistical analysis in prism !
op_1 = NDC';
op_2 = ALC';
op_3 = TA';
op_4 = AMC';
op_5 = RMC';

mkdir 1200_mcs\

names = {'Low HI - I','Low HI - II','Low HI < HI < Mid HI','Mid HI - I','Mid HI - II','High HI','C-C_high_adhesion','C-C_low_adhesion'};
 
% You need to convert these arrays/matrices into table and try to save them as .csv or .xlsx files instead of txt.

save('1200_mcs\1) NDC.txt','op_1','-ascii','-tabs');
save('1200_mcs\2) ALC.txt','op_2','-ascii','-tabs');
save('1200_mcs\3) TA.txt','op_3','-ascii','-tabs');
save('1200_mcs\3) AMC.txt','op_4','-ascii','-tabs');
save('1200_mcs\4) RMC.txt','op_5','-ascii','-tabs');



%  BOX PLOTS OF INDIVIDUAL  OUTPUTS !

fig_1 = figure, boxplot(op_1,'Notch','on','Labels',names);
hold on
plot(1:length(op_1(1,:)),mean(op_1),'Marker','_','MarkerSize',16,'Color','k','LineStyle','none','LineWidth',1);  % x-axis is the intergers of position
hold off
title('Number of disconnected cell clusters');
ax = gca, ax.FontSize = 10 ;
ax.FontName = 'Cambria Math';
fig_1.Position(3:4) = [600,300];
saveas(fig_1,'1)ndc.png');

fig_2 = figure, boxplot(op_2,'Notch','on','Labels',names);
hold on
plot(1:length(op_2(1,:)),mean(op_2),'Marker','_','MarkerSize',16,'Color','k','LineStyle','none','LineWidth',1);  % x-axis is the intergers of position
hold off
title('Area of largest cell cluster');
ax = gca, ax.FontSize = 10 ;
ax.FontName = 'Cambria Math';
fig_2.Position(3:4) = [600,300];
saveas(fig_2,'2)alc.png');

fig_3 = figure, boxplot(op_3,'Notch','on','Labels',names);
hold on
plot(1:length(op_3(1,:)),mean(op_3), 'Marker','_','MarkerSize',16,'Color','k','LineStyle','none','LineWidth',1);  % x-axis is the intergers of position
hold off
title('Total area');
ax = gca, ax.FontSize = 10 ;
ax.FontName = 'Cambria Math';
fig_3.Position(3:4) = [600,300];
saveas(fig_3,'3)ta.png');

fig_4 = figure, boxplot(op_4,'Notch','on','Labels',names);
hold on
plot(1:length(op_4(1,:)),mean(op_4), 'Marker','_','MarkerSize',16,'Color','k','LineStyle','none','LineWidth',1);  % x-axis is the intergers of position
hold off
title('Area of minimum enclosing cirlce');
ax = gca, ax.FontSize = 10 ;
ax.FontName = 'Cambria Math';
fig_4.Position(3:4) = [600,300];
saveas(fig_4,'4)amc.png');

fig_5 = figure, boxplot(op_5,'Notch','on','Labels',names);
hold on
plot(1:length(op_5(1,:)),mean(op_5), 'Marker','_','MarkerSize',16,'Color','k','LineStyle','none','LineWidth',1);  % x-axis is the intergers of position
hold off
title('Radius of minimum enclosing circle');
ax = gca, ax.FontSize = 10 ;
ax.FontName = 'Cambria Math';
fig_5.Position(3:4) = [600,300];
saveas(fig_5,'5)rmc.png');

Mean_NDC = round(mean(NDC,2));      % since ndc has to be an integer, rounding it to an integer
Mean_ALC = mean(ALC,2);
Mean_TA = mean(TA,2);
Mean_AMC = mean(AMC,2);
Mean_RMC = mean(RMC,2);


Std_NDC = std(NDC,0,2);
Std_ALC = std(ALC,0,2);
Std_TA  = std(TA,0,2);
Std_AMC = std(AMC,0,2);
Std_RMC = std(RMC,0,2);


Std_err_NDC = Std_NDC/sqrt(n);
Std_err_ALC = Std_ALC/sqrt(n);
Std_err_TA  = Std_TA/sqrt(n);
Std_err_AMC = Std_AMC/sqrt(n);
Std_err_RMC = Std_RMC/sqrt(n);


Median_NDC = median(NDC,2);
Median_ALC = median(ALC,2);
Median_TA  = median(TA,2);
Median_AMC = median(AMC,2);
Median_RMC = median(RMC,2);



%%%%%%%%%%%     ONE TABLE CONTAINING MEANS, MEDIANS, STD, STD.ERR FOR ALL CASES.
op_1_table = array2table([Mean_NDC';Median_NDC';Std_NDC';Std_err_NDC'],'VariableNames',names,'RowNames',{'Mean','Median','Std. Dev.','Std. Error'})
writetable(op_1_table,'1200_mcs\op_1_table.xlsx','WriteRowNames',true,'Range','A4:H10'); 
op_2_table = array2table([Mean_ALC';Median_ALC';Std_ALC';Std_err_ALC'],'VariableNames',names,'RowNames',{'Mean','Median','Std. Dev.','Std. Error'})
writetable(op_2_table,'1200_mcs\op_2_table.xlsx','WriteRowNames',true,'Range','A4:H10'); 
op_3_table = array2table([Mean_TA';Median_TA';Std_TA';Std_err_TA'],'VariableNames',names,'RowNames',{'Mean','Median','Std. Dev.','Std. Error'})
writetable(op_3_table,'1200_mcs\op_3_table.xlsx','WriteRowNames',true,'Range','A4:H10'); 
op_4_table = array2table([Mean_AMC';Median_AMC';Std_AMC';Std_err_AMC'],'VariableNames',names,'RowNames',{'Mean','Median','Std. Dev.','Std. Error'})
writetable(op_4_table,'1200_mcs\op_4_table.xlsx','WriteRowNames',true,'Range','A4:H10'); 
op_5_table = array2table([Mean_RMC';Median_RMC';Std_RMC';Std_err_RMC'],'VariableNames',names,'RowNames',{'Mean','Median','Std. Dev.','Std. Error'})
writetable(op_5_table,'1200_mcs\op_5_table.xlsx','WriteRowNames',true,'Range','A4:H10'); 


% %  Phenotypic space plots

g1 = phase_plot_mean_se(Mean_ALC,Mean_NDC,Std_err_ALC,Std_err_NDC,'Area of largest cell cluster','Number of disconnected cell clusters  vs  Area of largest cell cluster');
%g1.Position(3:4) = [1000,600];
saveas(g1,'6)ndc-alc.png');
g2 = phase_plot_mean_se(Mean_TA,Mean_NDC,Std_err_TA,Std_err_NDC,'Total area','Number of disconnected cell clusters  vs  Total area');
%g2.Position(3:4) = [1000,600];
saveas(g2,'7)ndc-ta.png');
g3 = phase_plot_mean_se(Mean_AMC,Mean_NDC,Std_err_AMC,Std_err_NDC,'Area of minimum enclosing circle','Number of disconnected cell clusters  vs  Area of minimum enclosing circle');
saveas(g3,'8)ndc-amc.png');
g4 = phase_plot_mean_se(Mean_RMC,Mean_NDC,Std_err_RMC,Std_err_NDC,'Radius of minimum enclosing circle','Number of disconnected cell clusters  vs  Radius of minimum enclosing circle');
saveas(g4,'9)ndc-rmc.png'); 

% save('600_mcs\anova_1.mat','p_anova','tbl_anova','stats_anova','-v7.3')
% save('600_mcs\anova_2.mat','p_anova2','tbl_anova2','stats_anova2','-v7.3')


% Statistical Analysis perfomed on outputs

[pa_1,ta_1,sa_1,c_1,m_1,h_1] = stat_analysis(op_1);
[pa_2,ta_2,sa_2,c_2,m_2,h_2] = stat_analysis(op_2);
[pa_3,ta_3,sa_3,c_3,m_3,h_3] = stat_analysis(op_3);
[pa_4,ta_4,sa_4,c_4,m_4,h_4] = stat_analysis(op_4);
[pa_5,ta_5,sa_5,c_5,m_5,h_5] = stat_analysis(op_5);

writecell(ta_1,'1200_mcs\op_1_table.xlsx','Range','A14:G18');
writecell(ta_2,'1200_mcs\op_2_table.xlsx','Range','A14:G18');
writecell(ta_3,'1200_mcs\op_3_table.xlsx','Range','A14:G18');
writecell(ta_4,'1200_mcs\op_4_table.xlsx','Range','A14:G18');
writecell(ta_5,'1200_mcs\op_5_table.xlsx','Range','A14:G18');


% NDC_p_matrix = p_mat(c_1)
% ALC_p_matrix = p_mat(c_2)
% TA_p_matrix  = p_mat(c_3)
% AMC_p_matrix = p_mat(c_4)
% RMC_p_matrix = p_mat(c_5)

%op_1_p_matrix = array2table(NDC_p_matrix,'VariableNames',names,'RowNames',names);
%writetable(op_1_p_martix,'1200_mcs\op_1_table.xlsx','Range','B24:G30');
% writematrix(NDC_p_matrix,'1200_mcs\op_1_table.xlsx','Range','B24:I30');
% writematrix(ALC_p_matrix,'1200_mcs\op_2_table.xlsx','Range','B24:I30');
% writematrix(TA_p_matrix,'1200_mcs\op_3_table.xlsx','Range','B24:I30');
% writematrix(AMC_p_matrix,'1200_mcs\op_4_table.xlsx','Range','B24:I30');
% writematrix(RMC_p_matrix,'1200_mcs\op_5_table.xlsx','Range','B24:I30');


% Tukey -  Krammer statistical test inferences
% NDC_tk_inf = NDC_p_matrix <= .05;
% ALC_tk_inf = ALC_p_matrix <= .05;
% TA_tk_inf  = TA_p_matrix <= .05;
% AMC_tk_inf = AMC_p_matrix <= .05;
% RMC_tk_inf = RMC_p_matrix <= .05;

% [p_anova,tbl_anova,stats_anova] = anova1(op_2,names);
% % Calculate and Display Std. Dev. and Variance also
% [c,m,h,gnames] = multcompare(stats_anova,"CType","tukey-kramer","Display","off","alpha",.05); 
% %[c,m,h,gnames] = multcompare(stats_anova,"CriticalValueType","tukey-kramer","Display","on");  




%%
%%%% SECTION 3 : Functions used in the matlab code

function [sdist,spos] = sdistFromPoint(sC1,sC2,Xmm,Ymm)
    sdistSq = ( Xmm - sC1 ).^2 + ( Ymm - sC2 ).^2;
    [sdist,spos] = min(sqrt(sdistSq));
end


function r = radiusFromPoint(P,pointsX,pointsY)
    px = P(1);
    py = P(2);
    distSq = ( pointsX - px ).^2 + ( pointsY - py ).^2;
    r = sqrt(max(distSq));
end


function [p_anova,tbl_anova,stats_anova,c,m,h] = stat_analysis(matrix)
    names = {'Low HI - I','Low HI - II','Low HI < HI < Mid HI','Mid HI - I','Mid HI - II','High HI','C-C_low_adhesion','C-C_high_adhesion'};
    [p_anova,tbl_anova,stats_anova] = anova1(matrix,names);
    % Calculate and Display Std. Dev. and Variance also
    [c,m,h,gnames] = multcompare(stats_anova,"CType","tukey-kramer","Display","on","Alpha",0.05); % You can add Alpha parameter as well !
    %[c,m,h,gnames] = multcompare(stats_anova,"CriticalValueType","tukey-kramer","Display","on");        
end


function [g] = phase_plot_mean_se(vec_1,vec_2,vec_3,vec_4,xl,ti)  % (x,y,x_err,y_err)
    g = figure, errorbar(vec_1(1),vec_2(1),vec_4(1),vec_4(1),vec_3(1),vec_3(1),'o','Color','#5f4b8b','LineWidth',1.5);  % UV
    hold on
    errorbar(vec_1(2),vec_2(2),vec_4(2),vec_4(2),vec_3(2),vec_3(2),'ob','LineWidth',1.5);  % Blue
    errorbar(vec_1(3),vec_2(3),vec_4(3),vec_4(3),vec_3(3),vec_3(3),'o','Color','#B88A00','LineWidth',1.5);   % dark yellow
    errorbar(vec_1(4),vec_2(4),vec_4(4),vec_4(4),vec_3(4),vec_3(4),'o','Color','#FFA500','LineWidth',1.5);   % orange-yellow
    errorbar(vec_1(5),vec_2(5),vec_4(5),vec_4(5),vec_3(5),vec_3(5),'o','Color','#04878F','LineWidth',1.5);  % swap - green
    errorbar(vec_1(6),vec_2(6),vec_4(6),vec_4(6),vec_3(6),vec_3(6),'or','LineWidth',1.5);  % Chk
    errorbar(vec_1(7),vec_2(7),vec_4(7),vec_4(7),vec_3(7),vec_3(7),'ok','LineWidth',1.5);  % 
    errorbar(vec_1(8),vec_2(8),vec_4(8),vec_4(8),vec_3(8),vec_3(8),'o','Color','#00377A','LineWidth',1.5);  % 
    %xlim([5000 10500]);
    ylim([0 5]);
    xlabel(xl);
    ylabel('Number of disconnected cell clusters');
    %title(ti);
    legend({'1) Low HI - I','2) Low HI - II','3) Low HI < HI < Mid HI','4) Mid HI - I','5) Mid HI - II','6) High HI',...
        ['7) Homogeneous case of C-C adhesion' newline 'double of model adhesion'],...
        ['8) Homogeneous case of C-C adhesion' newline 'half of model adhesion']},'location','northwest');
    hold off
    ax = gca, ax.FontSize = 12 ;
    set(ax,'box','off');
    ax.FontName = 'Cambria Math';
    g.Position = [100,100,800,400];
end

% 
% function P_Matrix = p_mat(mat)
% 
%     P_Matrix(1,:) = [1        mat(1,6)  mat(2,6)   mat(3,6)  mat(4,6)  mat(5,6)  mat(6,6)  mat(7,6)];
%     P_Matrix(2,:) = [mat(1,6) 1         mat(8,6)   mat(9,6)  mat(10,6) mat(11,6) mat(12,6) mat(13,6)];
%     P_Matrix(3,:) = [mat(2,6) mat(8,6)  1          mat(14,6) mat(15,6) mat(16,6) mat(17,6) mat(18,6)];
%     P_Matrix(4,:) = [mat(3,6) mat(9,6)  mat(14,6)  1         mat(19,6) mat(20,6) mat(21,6) mat(22,6)];
%     P_Matrix(5,:) = [mat(4,6) mat(10,6) mat(15,6)  mat(19,6) 1         mat(23,6) mat(24,6) mat(25,6)];
%     P_Matrix(6,:) = [mat(5,6) mat(11,6) mat(16,6)  mat(20,6) mat(23,6) 1         mat(26,6) mat(27,6)];
%     P_Matrix(7,:) = [mat(6,6) mat(12,6) mat(17,6)  mat(21,6) mat(24,6) mat(26,6) 1         mat(28,6)];
%     P_Matrix(8,:) = [mat(7,6) mat(13,6) mat(18,6)  mat(22,6) mat(25,6) mat(27,6) mat(28,6) 1     ];
%     P_Matrix ;
% 
% end



