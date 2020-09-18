function [sat_f,pres_f] = SWINTcmg(sw_b,pres,tstamp,en)
%SWINT Generates Forecasts for saturation and pressure evolution by
%using your input water saturation (sw), pressure( kPa), and time (stamp),
%and ensemble size(default use 1).

%create a text file to copy input sw,sg, and pres (prepping data for .dat
%file)

tstamp1=num2str(tstamp); %marks your files with the tstamp input
enstamp=num2str(en); %Originally this function was created as a subroutine to study ensemble-based data assimilation techniques. (sorry! if this isn't the ideal code)
swbtext=strcat('swb',tstamp1,enstamp);
sgtext=strcat('sg',tstamp1,enstamp);
prestext=strcat('pres',tstamp1,enstamp);
swbtext=strcat(swbtext,'.txt');
sgtext=strcat(sgtext,'.txt');
prestext=strcat(prestext,'.txt');
[x,y]=size(sw_b);
fid=fopen(swbtext,'w');
fid2=fopen(sgtext,'w');
fid3=fopen(prestext,'w');
tab=' ';
for j=1:x
   for i=1:y
      fprintf(fid,'%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c',num2str(sw_b(i,j),'%2.2f'),tab) ;
      fprintf(fid2,'%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c',num2str(1-sw_b(i,j),'%2.2f'),tab);
      fprintf(fid3,'%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c',num2str(pres(i,j),'%2.2f'),tab);
   end
   fprintf(fid,'\n');
   fprintf(fid2,'\n');
   fprintf(fid3,'\n');
end


%% Function to generate the SWINT indexing for .dat file
cmgdat=strcat('DAcmgfile',tstamp1,enstamp);
fileID2=strcat(cmgdat,'.out');
cmgdat=strcat(cmgdat,'.dat');
[~]=cmgcall(cmgdat,swbtext,sgtext,prestext,tstamp);%Writes out a .dat file using a predefined reservoir. (see
command1=['./RunSim.sh gem 2019.10 ',cmgdat]; %this command should be changed based on version you use
system(command1);
command1=['cd ..'];
system(command1);
[sat_f pres_f]=cmgread(fileID2,1);
fclose(fid);
fclose(fid2);
fclose(fid3);
end
