function quantquote_raw_to_matlab(path_raw_data, filename_matlab_data)
%QUANTQUOTE_RAW_TO_MATLAB transforms rc data calculated by R in csv format
%to MatLab data.

cd(path_raw_data)

gunzip("vech_rc.csv.gz")
vech_rc = readmatrix("vech_rc.csv");
delete vech_rc.csv

rc = ivech_(vech_rc);

dates = readmatrix("dates.csv", 'OutputType', 'datetime');

symbols  = readmatrix("symbols.csv", 'OutputType', 'char');

log_ret_otc = readmatrix("log_ret_otc.csv");

dt_5min = readtable("log_ret_5min.csv", 'Range', 'A:A');
dt_5min = datetime(dt_5min{:,1}, 'InputFormat', 'uuuu-MM-dd''T''HH:mm:ss''Z');

log_ret_5min = readmatrix("log_ret_5min.csv", 'Range', 'B:end');

save(filename_matlab_data, 'dates', 'symbols', 'vech_rc', 'rc', 'log_ret_otc', 'dt_5min', 'log_ret_5min')

end
