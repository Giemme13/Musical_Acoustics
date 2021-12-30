function ph = grpdelay2phase(grd)
% GRPDELAY2PHASE
% This function get the phase from the group delay
%
% Musical Acoustics Course
% Riccardo R. De Lucia
% 2018
% Mirco Pezzoli
% 2019-20

ph = -cumsum(grd);
ph = 2*pi*ph/length(grd);
