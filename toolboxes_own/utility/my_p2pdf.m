function my_p2pdf( fig_, cm_width, cm_height, name, path )
%MY_P2PDF prints figure to pdf on desktop.
%
% USAGE:
%  [FCST_COVM] = MY_P2PDF(FIG_,CM_WIDTH,CM_HEIGHT,NAME,PATH)
%
% INPUTS:
%   FIG_      - figure handle
%   CM_WIDTH  - width in cm
%   CM_HEIGHT - height in cm
%   NAME      - [Optional] Name of pdf (default: print)
%   PATH      - [Optional] Path of pdf (default: Desktop)
%
% COMMENTS:
%
%  See also 
%

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 19.02.2017

if isempty(name)
    name='print2pdf';
end
if isempty(path)
    [~, userdir] = system('echo %USERPROFILE%');
    path=strcat(userdir,'\Desktop');
end

set(fig_,'PaperUnits','centimeters','PaperPosition',[0 0 cm_width cm_height],'PaperSize',[cm_width cm_height])
set(fig_.CurrentAxes,'LooseInset',get(gca,'TightInset'))
print(fig_,strcat(path,'\',name),'-dpdf','r1200')

end


