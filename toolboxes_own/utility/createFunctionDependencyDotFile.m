function createFunctionDependencyDotFile(calls)
%CREATEFUNCTIONDEPENDENCYDOTFILE Create a GraphViz DOT diagram file from function call list
%
% Calls (cellstr) is an n-by-2 cell array in format {caller,callee;...}.
%
% Example:
% calls = { 'foo','X'; 'bar','Y'; 'foo','Z'; 'foo','bar'; 'bar','bar'};
% createFunctionDependencyDotFile(calls)

baseName = 'functionCalls';
dotFile = [baseName '.dot'];
fid = fopen(dotFile, 'w');
fprintf(fid, 'digraph G {\n');
for i = 1:size(calls,1)
    [parent,child] = calls{i,:};
    fprintf(fid, '   "%s" -> "%s"\n', parent, child);
end
fprintf(fid, '}\n');
fclose(fid);

% Render to image
imageFile = [baseName '.png'];
% Assumes the GraphViz bin dir is on the path; if not, use full path to dot.exe
cmd = sprintf('dot -Tpng -Gsize="2,2" "%s" -o"%s"', dotFile, imageFile);
system(cmd);
fprintf('Wrote to %s\n', imageFile);