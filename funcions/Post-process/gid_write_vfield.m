function [ ] = gid_write_vfield(fid,nameres,time,vfield)

npnod = size(vfield,1);
s =['Result' ' "' nameres '" ' '"time"' '%12.5d' ' Vector OnNodes \n'];
fprintf(fid,s,time);
fprintf(fid,['Values \n']);
for i = 1 : npnod
   fprintf(fid,['%6.0i %13.5d %13.5d %13.5d \n'],i,vfield(i,:));
end
fprintf(fid,['End Values \n']);
fprintf(fid,'# \n');

end

