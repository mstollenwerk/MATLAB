function [dates1,dates2]=createDates(ini,fin)
%==========================================================================
% Create vector of dates given an initial and an ending point.
%
%INPUT:
% ini (string) the initial point (EX: '1-Jan-2007')
% fin (string) the final   point (EX: '31-Jan-2010')
%
%OUTPUT:
% dates1 (string) column vector of days excluding holidays
% dates2 (numeric) column vector of days excluding holidays
%==========================================================================
t1 = datenum(ini);
t2 = datenum(fin);
full_dates=t1:t2;
working_days = isbusday(full_dates, [], [1 0 0 0 0 0 1]);
dates1=full_dates(working_days)';
dates2=datestr(dates1);
end
