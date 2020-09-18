fid1 = fopen('sandtank2D.dat');
   tline = fgets(fid1);
   fid2 = fopen('sandtankMATLAB2.txt','w');
   Character2addfront = 'fprintf(fid,''\n'; % here put the character to be added
   Character2addend = ''');';
   while ischar(tline)
      disp(tline)
      fprintf(fid2,'%s %s %s \n',Character2addfront,tline);
      tline = fgets(fid1);
   end
   

   fclose(fid1);
   fclose(fid2);
   