function [SAT,PRES] = cmgread(filename,td)
%cmgread function reads data from an cmg .out file
%   inputs: filename: a string represnting the name of the .out file'
%           td: time duration of the cmg simulation eg: 365 days


fileID=fopen(filename);
file=textscan(fileID,'%s','delimiter','\n');

for time=1:td
    timestampstart='Time =';
    timestampnum=num2str(time);
    timestampend='  ';
    timestamp=strcat(timestampstart,{' '},timestampnum,{' '});
    %ts='Time = 1                     **********************************************************************         2019 JAN.  2'
    k=~cellfun(@isempty,regexp(file{1},timestamp));
    loc(:,time)=find(k);
    
end

%% PRESSURE
for time=1:td
    clear sl sl2 sl3 sl4
    %read pressure
    index=loc(1,time);
    index=index+5;
    count=1;
    for idk=index:index+39
        sl(:,count)=split(file{1}{idk});
        for i=3:16
            sl2(i-2,count)=str2double(sl(i,count));
        end
        count=count+1;
    end
    clear sl
    count=1;
    for idk2=idk+3:idk+42
        sl(:,count)=split(file{1}{idk2});
        for i=3:16
            sl3(i-2,count)=str2double(sl(i,count));
        end
        count=count+1;
    end
    
    clear sl
    count=1;
    for idk3=idk2+3:idk2+42
        sl(:,count)=split(file{1}{idk3});
        for i=3:14
            sl4(i-2,count)=str2double(sl(i,count));
        end
        count=count+1;
    end
    clear sl
    PRES(:,1:14,time)=sl2';
    PRES(:,15:28,time)=sl3';
    PRES(:,29:40,time)=sl4';
end
%% SATURATION
for time=1:td
    clear sl sl2 sl3 sl4
    %read pressure
    index=loc(2,time);
    index=index+5;
    count=1;
    for idk=index:index+39
        sl(:,count)=split(file{1}{idk});
        for i=3:16
            sl2(i-2,count)=str2double(sl(i,count));
        end
        count=count+1;
    end
    clear sl
    count=1;
    for idk2=idk+3:idk+42
        sl(:,count)=split(file{1}{idk2});
        for i=3:16
            sl3(i-2,count)=str2double(sl(i,count));
        end
        count=count+1;
    end
    
    clear sl
    count=1;
    for idk3=idk2+3:idk2+42
        sl(:,count)=split(file{1}{idk3});
        for i=3:14
            sl4(i-2,count)=str2double(sl(i,count));
        end
        count=count+1;
    end
    clear sl
    SAT(:,1:14,time)=sl2';
    SAT(:,15:28,time)=sl3';
    SAT(:,29:40,time)=sl4';
end
fclose(fileID);
end

