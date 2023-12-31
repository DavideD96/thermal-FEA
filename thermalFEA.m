magnification = 1000;

substrate_length = 5;
substrate_width = 4;
substrate_thickness = 0.7;

CAF_length = 0.5;
CAF_width = 0.5;
CAF_thickness = 0.05;

electrode_thickness = 0.16;
electrode_length = 1.5;
electrode_width = 1;

CAF_highedge = 0.5*(substrate_width+CAF_width);
CAF_lowedge = 0.5*(substrate_width-CAF_width);

Elsx_edge = 0.5*(substrate_length-CAF_length);
Eldx_edge = 0.5*(substrate_length+CAF_length);

El_upedge = 0.5*(substrate_width-electrode_width)+electrode_width;
El_dwedge = 0.5*(substrate_width-electrode_width);

Elsx_sxedge = 0.5*(substrate_length-CAF_length)-electrode_length;
Eldx_dxedge = substrate_length-Elsx_sxedge;

SublowerLeft  = [0 , 0];
SublowerRight = [substrate_length , 0];
SubupperRight = [substrate_length , substrate_width];
SubupperLeft =  [0 , substrate_width];

ElsxlowerLeft = [Elsx_sxedge , El_dwedge];
ElsxlowerRight = [Elsx_edge , El_dwedge];
ElsxupperRight = [Elsx_edge , El_upedge];
ElsxupperLeft =  [Elsx_sxedge , El_upedge];

EldxlowerLeft = [Eldx_edge , El_dwedge];
EldxlowerRight = [Eldx_dxedge , El_dwedge];
EldxupperRight = [Eldx_dxedge , El_upedge];
EldxupperLeft =  [Eldx_edge , El_upedge];

CAFlowerLeft = [Elsx_edge, CAF_lowedge];
CAFlowerRight = [Eldx_edge, CAF_lowedge];
CAFupperRight = [Eldx_edge, CAF_highedge];
CAFupperLeft = [Elsx_edge, CAF_highedge];

Sub = [3;4;SublowerLeft(1);SublowerRight(1);SubupperRight(1);SubupperLeft(1); ...
    SublowerLeft(2);SublowerRight(2);SubupperRight(2);SubupperLeft(2)];
Elsx = [3;4;ElsxlowerLeft(1);ElsxlowerRight(1);ElsxupperRight(1);ElsxupperLeft(1); ...
    ElsxlowerLeft(2);ElsxlowerRight(2);ElsxupperRight(2);ElsxupperLeft(2)];
Eldx = [3;4;EldxlowerLeft(1);EldxlowerRight(1);EldxupperRight(1);EldxupperLeft(1); ...
    EldxlowerLeft(2);EldxlowerRight(2);EldxupperRight(2);EldxupperLeft(2)];
CAF = [3;4;CAFlowerLeft(1);CAFlowerRight(1);CAFupperRight(1);CAFupperLeft(1); ...
    CAFlowerLeft(2);CAFlowerRight(2);CAFupperRight(2);CAFupperLeft(2)];

%C1 = [1;0;0;0.4];
%C1 = [C1;zeros(length(R1) - length(C1),1)];

gd = Sub;
sf = 'Sub';
ns = char('Sub')'; %OCCHIO ALL'APICE
gSub = decsg(gd,sf,ns);



gd = [Sub, Elsx, Eldx, CAF];
sf = 'Sub+Elsx+Eldx+CAF';
ns = char('Sub','Elsx','Eldx','CAF')'; %OCCHIO ALL'APICE
g = decsg(gd,sf,ns);

pdegplot(g,'FaceLabels','on');
%[1,2],[-2,0.5])

%g = addFace(g,[9 5 3 7]);
%g = addFace(g,[1 3 8 4]);
%g = addFace(g,[2 7 10 6]);
tmodel = createpde('thermal', 'steadystate');
g = geometryFromEdges(tmodel,g);

figure;
g = extrude(g,substrate_thickness); % estrudo substrato
tmodel.Geometry = g;
%pdegplot(tmodel,'CellLabels','on','EdgeLabel','on','FaceLabels','on');


g = extrude(g,[6,7,8],CAF_thickness); % estrudo elettrodi e CAF di CAF_thickness

%pdegplot(tmodel,'CellLabels','on','EdgeLabel','on','FaceLabels','on');
g = extrude(g,[27,29],electrode_thickness-CAF_thickness); %estrudo elettrodi
%g = extrude(g,6,0.5);
tmodel.Geometry = g;
pdegplot(tmodel,'CellLabels','on','EdgeLabel','off','FaceLabels','on');

generateMesh(tmodel, 'Hmax',0.15);
pdemesh(tmodel);

%ATTENZIONE: FORSE MATLAB USA GLI inc INVECE CHE I metri. Convertire le
%conducibilità termiche!
thermalProperties(tmodel, "Cell",[1 2 3 4], "ThermalConductivity",1*magnification); %substrato vetro
thermalProperties(tmodel, "Cell",[5 8 7 9], "ThermalConductivity",[300/(magnification^2);300/(magnification^2); 300]); %elettrodo dx
%thermalProperties(tmodel, "Cell",[13 14 15 16 17 18], "ThermalConductivity",314*magnification); %elettrodo sx
thermalProperties(tmodel, "Cell",6, "ThermalConductivity",[100/(magnification^2); 100/(magnification^2); 100]); %CAF

internalHeatSource(tmodel, 50, "Cell",6);
thermalBC(tmodel, 'Face',[1,2,3,4],'Temperature',300); %sotto
thermalBC(tmodel, 'Face',[9,10,24,25],'Temperature',300); %lati
%thermalBC(tmodel, 'Face',[57,43,46,30],'Temperature',300); %lati elettrodi

result = solve(tmodel);
T = result.Temperature;
figure
pdeplot3D(tmodel,ColorMapData=T)