function fsdt_write_csv(path,headers,data)
if size(data,2)~=numel(headers), error('CSV header/data size mismatch.'); end
fid=fopen(path,'w'); if fid<0, error('Cannot write %s.',path); end
cleanup=onCleanup(@() fclose(fid));
for j=1:numel(headers)
    if j>1, fprintf(fid,','); end
    fprintf(fid,'%s',headers{j});
end
fprintf(fid,'\n');
fmt=[repmat('%.15g,',1,size(data,2)-1) '%.15g\n'];
fprintf(fid,fmt,data');
end
