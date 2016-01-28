clear
% MAIN PROGRAM: EXAMPLE OF POSTPROCESSING 
% - 3D TRUSS STRUCTURE: POSTPROCESS
% - ELEMENT: BAR
% - RESULTS: DISPLACEMENT AND AXIAL STRESS
fname = 'TASK_B';

% list of nodal coordinates
coord = [ 
0 0 0
10 0 0
10 0 10
10 10 10
10 10 0             
];         

% Table of connectivities
element.connec = [
1 2 
2 3
3 4 
4 5
];

% POSTPROCESS INFORMATION
% element type (Linear means two nodes)
element.etype = 'Linear';

% NODAL DISPLACEMENTS
post.disp = [
-0.5 0.5 -0.5
0.5 -0.5 0.5
-0.5 0.5 -0.5
0.5 -0.5 0.5
-0.5 0.5 -0.5
];

% ELEMENT AXIAL STRESS
post.stress = [           
-2.00475e-001 
-3.25376e-001 
-3.25704e-001 
-2.42826e-001 ]; 

%===============================================
postinfo(fname,coord,element,post)


