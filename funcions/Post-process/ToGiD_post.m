function ToGiD_post(file_name,iter,nnode,ngaus,disp)

gtype = 'Linear';
% Escribe el fichero de resultados
res_file = strcat(file_name,'_',num2str(iter),'.flavia.res');
fid = fopen(res_file,'w');
job=1;
gid_write_headerpost(fid,gtype,ngaus,job)


% DISPLACEMENT
nameres = 'DISPLACEMENT';
gid_write_vfield(fid,nameres,iter,disp);

fclose(fid);

end
