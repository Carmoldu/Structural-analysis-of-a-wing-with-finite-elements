function ToGID (file_name,istep,xnodes,conectivities,nnode,nelem,npnod)
% Construcción de malla postproceso

gtype = 'Linear';
msh_file = strcat(file_name,'_',num2str(istep),'.flavia.msh');
fid = fopen(msh_file,'w');

fprintf(fid,['MESH "WORKPIECE" dimension %2.0f  Elemtype %s   Nnode %2.0f \n \n'],3,gtype,nnode);
fprintf(fid,['Coordinates \n']);
for i = 1 : npnod
    fprintf(fid,['%6.0f %12.5d %12.5d %12.5d \n'],i,xnodes(i,:));
end
fprintf(fid,['End Coordinates \n \n']);
fprintf(fid,['Elements \n']);

fprintf(fid,['%6.0f %6.0f %6.0f \n'],[1:nelem;conectivities']);

fprintf(fid,['End elements \n \n']);


fclose(fid);

end
