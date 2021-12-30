function mypolarplot(theta,rho,max_rho,rho_ticks,rho_ticks_label,theta_ticks,rotation,color,linestyle,linewidth,fontsize,range,pos_rho_ticks)
% MYPOLARPLOT
% This utility helps in plotting radiance patterns in polar plots in a
% intuitive fashion
%
% Musical Acoustics Course
% Riccardo R. De Lucia
% 2018
% Mirco Pezzoli
% 2019-21

if nargin<12
    range='full';
    pos_rho_ticks = 90;
end


hold on;
rho_ticks = sort(rho_ticks);
nRhoTicks = length(rho_ticks);
if strcmp(range,'full')
    viscircles(repmat([0 0],nRhoTicks,1),rho_ticks(1:end),'EdgeColor','k','linewidth',1,'linestyle',':');
    viscircles([0 0],max_rho,'EdgeColor','k','linewidth',1,'linestyle','-');
elseif strcmp(range,'half_up')
    viscircles_halfUp(repmat([0 0],nRhoTicks,1),rho_ticks(1:end),'EdgeColor','k','linewidth',1,'linestyle',':');
    viscircles_halfUp([0 0],max_rho,'EdgeColor','k','linewidth',1,'linestyle','-');
elseif strcmp(range,'half_right');
    viscircles_halfRight(repmat([0 0],nRhoTicks,1),rho_ticks(1:end),'EdgeColor','k','linewidth',1,'linestyle',':');
    viscircles_halfRight([0 0],max_rho,'EdgeColor','k','linewidth',1,'linestyle','-');
end

lx = max_rho * cosd(theta_ticks + rotation);
ly = max_rho * sind(theta_ticks + rotation);
lx_t = (max_rho + max_rho*0.12) * cosd(theta_ticks + rotation);
ly_t = (max_rho + max_rho*0.12) * sind(theta_ticks + rotation);
for i = 1:length(theta_ticks)
    plot([0 lx(i)],[0 ly(i)],'k:');
    text(lx_t(i),ly_t(i),[num2str(theta_ticks(i)) '^\circ'],...
         'HorizontalAlignment','center',...
         'FontSize',fontsize(1));
end

for i = 1:length(rho_ticks)
    x = rho_ticks(i) * cosd(pos_rho_ticks);
    y = rho_ticks(i) * sind(pos_rho_ticks);
    text(x,y,rho_ticks_label{i},'fontsize',fontsize(2),...
         'HorizontalAlignment','center',...
         'VerticalAlignment','baseline');
end

pts_x = rho .* cosd(theta);
pts_y = rho .* sind(theta);

for i = 1:size(pts_x,2)
    plot(pts_x(:,i),pts_y(:,i),'color',color{i},'linestyle',linestyle{i},'linewidth',linewidth{i});
end


axis off
axis equal
axis([-38.036809815950917  38.036809815950917 -30 30]);