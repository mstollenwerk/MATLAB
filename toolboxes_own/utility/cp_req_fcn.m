function cp_req_fcn( parentpath, destpath )
%CP_REQ_FCN creates a folder with all user-made functions that parentFcn
%depends on. Location: desktop
%
% USAGE:
%  cp_req_fcn('fcn')
%
% INPUTS:
%   PARENT_FCN - parent function
%
% COMMENTS:
%
%  See also 
%

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 08.08.2017
% 05.09.2017: changed description
% 10.09.2021: added destination folder input

% parentpath = which(parentFcn);
FileList = matlab.codetools.requiredFilesAndProducts(parentpath);
% basedir = winqueryreg('HKEY_CURRENT_USER', 'Software\Microsoft\Windows\CurrentVersion\Explorer\Shell Folders', 'Desktop');
% dest_fldr = [basedir '\reqFcns'];
mkdir(destpath)
copyfile(parentpath,destpath)
FileList = arrayfun(@(x) x,FileList);
cellfun(@(x) copyfile(x,destpath),FileList);
end
