function [] = radianceplot(ctr_freqs, rad_patt, angs, figureTitle)
% RADIANCEPLOT
% This utility plot the radiance pattern in polar coordinates at different
% frequency
%
% Musical Acoustics Course
% Riccardo R. De Lucia
% 2018
% Mirco Pezzoli
% 2019-21
figure
tiledlayout('flow')
for f = 1:length(ctr_freqs)
      
    rad_patt_f = mean(rad_patt( (ctr_freqs(f)+1)-50:(ctr_freqs(f)+1)+50 , :) , 1);
    rad_patt_f = rad_patt_f/max(rad_patt_f);
    rad_patt_f_dB = db(rad_patt_f,'voltage')+30;
    %%%%
    rad_patt_f_dB(rad_patt_f_dB<0) = 0;
    %%%%
    rad_patt_f_dB = [rad_patt_f_dB rad_patt_f_dB(1)];
    
    theta = circshift(angs,-6,2);
    theta = [theta theta(1)];
%     figure, clf
    nexttile
    mypolarplot(theta(:),rad_patt_f_dB(:),30,[10 20],{'-20dB','-10dB'},[0:45:315],90,{'b'},{'-'},{2},[10 10],'full',90)
%     text(0,40,[figureTitle num2str(ctr_freqs(f)) ' Hz'],'fontsize',20,'HorizontalAlignment','center');
    axis([-35 35 -35 40])
    title([figureTitle num2str(ctr_freqs(f)) ' Hz'])
    
%     pause
end


end