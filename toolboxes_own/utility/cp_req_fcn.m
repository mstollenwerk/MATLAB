function cp_req_fcn( parentFcn )
%CP_REQ_FCN creates a folder with all user-made functions that parentFcn
%depends on. Location: desktop
%
% USAGE:
%  cp_req_fcn('fcn')
%
% INPUTS:
%   PARENT_FCN   - parent function
%
% COMMENTS:
%
%  See also 
%

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 08.08.2017
% 05.09.2017: changed description

parentpath = which(parentFcn);
FileList = matlab.codetools.requiredFilesAndProducts(parentFcn);
basedir = winqueryreg('HKEY_CURRENT_USER', 'Software\Microsoft\Windows\CurrentVersion\Explorer\Shell Folders', 'Desktop');
newfolder = [basedir '\reqFcns'];
mkdir(newfolder)
copyfile(parentpath,newfolder)
FileList = arrayfun(@(x) x,FileList);
cellfun(@(x) copyfile(x,newfolder),FileList);
end
