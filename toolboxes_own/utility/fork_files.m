function fork_files(files_, rootdir, destdir)

for ii = 1:length(files_)
    cpdir = strcat(rootdir,files_{ii});
    dest = strcat(destdir,files_{ii});
    mkdir(dest)
    copyfile(cpdir, dest)
end

end