function fork_files(files_, rootdir, destdir)

for ii = 1:length(files_)
    cpfile = strcat(rootdir,files_{ii});
    backslashes = strfind(files_{ii},'\');
    dest = strcat(destdir,files_{ii}(1:backslashes(end)));
    mkdir(dest)
    copyfile(cpfile, dest)
end

end