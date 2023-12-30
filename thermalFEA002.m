magnification = 1000;
stretching = 1000;

%use meters, then multiply for magnification and stretching

substrate_length = 0.005*magnification; %0.005 [m]
substrate_width = 0.004*magnification;
substrate_thickness = 0.0007*magnification;

CAF_length = 0.0005*magnification;
CAF_width = 0.0005*magnification;
CAF_thickness = 25e-9*magnification*stretching;

electrode_thickness = 100e-9*magnification*stretching;
electrode_length = 0.0015*magnification;
electrode_width = 0.001*magnification;

%NOTA IMPORTANTE!
%quando magnifichi tutto il sistema (tutte e tre le direzioni), le tre componenti della conducibilità
%termica si modificano allo stesso modo: per un ingrandimento di un fattore
%C, esse si ridurranno di un fattore C. Ad esempio 300 W/(m*K) corrispondono
%a 0.3 W/(mm*K).
%
%Se invece ingrandisci solo una direzione (stretching), tutte e tre si modificano, ma in
%modo diverso: le due componenti trasverse rispetto allo stretching si
%riducono di C, quella parallela aumenta di C.

kvetro = 1/magnification;
kelectr = [300/stretching; 300/stretching; 300*stretching]/magnification;
kcaf = [100/stretching; 100/stretching; 100*stretching]/magnification;

q = 5; %internal heat source = rho*j^2 = (1e-6 Ohm/m) * (1e+9 A/m^2)^2 / ( magnification^3 * stretching )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
figure
pdegplot(tmodel,'CellLabels','on','EdgeLabel','off','FaceLabels','on');

figure
generateMesh(tmodel, 'Hmax',0.1);
pdemesh(tmodel);

%ATTENZIONE: FORSE MATLAB USA GLI inc INVECE CHE I metri. Convertire le
%conducibilità termiche!

thermalProperties(tmodel, "Cell",[1 2 3 4], "ThermalConductivity",kvetro); %substrato vetro
thermalProperties(tmodel, "Cell",[5 8 7 9], "ThermalConductivity",kelectr); %elettrodo dx
%thermalProperties(tmodel, "Cell",[13 14 15 16 17 18], "ThermalConductivity",314*magnification); %elettrodo sx
thermalProperties(tmodel, "Cell",6, "ThermalConductivity",kcaf); %CAF

internalHeatSource(tmodel, q, "Cell",6);
thermalBC(tmodel, 'Face',[1,2,3,4],'Temperature',300); %sotto
thermalBC(tmodel, 'Face',[9,10,24,25],'Temperature',300); %lati
%thermalBC(tmodel, 'Face',[57,43,46,30],'Temperature',300); %lati elettrodi

%thermalIC(tmodel,0); non necessario, siamo nel caso stazionario
tmodel.StefanBoltzmannConstant = 5.670367e-8; %SI units

result = solve(tmodel);
%result = solvepde(tmodel);
T = result.Temperature;
figure
pdeplot3D(tmodel,ColorMapData=T);

%quiver plot
[qx,qy,qz] = evaluateHeatFlux(result);

figure
pdeplot3D(result.Mesh,FlowData=[qx qy qz])

%slice plot
x=tmodel.Mesh.Nodes(1,:)';       %  x-coordinate of nodes
y=tmodel.Mesh.Nodes(2,:)';       %  y-coordinate of nodes
z=tmodel.Mesh.Nodes(3,:)';       %  z-coordinate of nodes
u= T(1,:)';        %  temperature of nodes


% figure
% slice(x,y,z,T,[1,0.5,1],[],[],'nearest')

% [X,Y,Z] = meshgrid(0:substrate_length/100:substrate_length,0:substrate_width/100:substrate_width,0:0.01:substrate_thickness);
% V = interpolateSolution(result,X,Y,Z);
% V = reshape(V,size(X));
% figure
% colormap jet
% contourslice(X,Y,Z,V,[],[],0:5:60)
% xlabel("x")
% ylabel("y")
% zlabel("z")
% colorbar
% view(-11,14)
% axis equal