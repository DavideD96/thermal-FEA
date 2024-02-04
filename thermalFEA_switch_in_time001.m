magnification = 1000;
stretching = 5000;

times = linspace(0.001,700.501,2);

%use meters, then multiply for magnification and stretching

substrate_length = 0.002*magnification; %0.005 [m]
substrate_width = 0.002*magnification;
substrate_thickness = 0.0007*magnification;

CAF_radius = 0.000125*magnification;
CAF_thickness = 22e-9*magnification*stretching;

h = 0.15; %mesh parameter

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
kcaf = [100/stretching; 100/stretching; 100*stretching]/magnification;

cCAF = 129E+3; %J/(kg*K)
rhoCAF = ((19.2E+3/magnification^3)/stretching)*0.7; %gold: 19.2 g/cm^3 = 0.0192 kg/cm^3 = 19200 kg/m^3 but CAF has high percentage of vacuum => x 0.7

cVetro = 840E+3; %J/(kg*K)
rhoVetro = 2.5E+3/magnification^3; %(kg/m^3)/magnification^3

q = 7.4e+12/( magnification^3 * stretching ); %internal heat source = rho*j^2 = (1e-6 Ohm/m) * (1e+9 A/m^2)^2 / ( magnification^3 * stretching )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SublowerLeft  = [0 , 0];
SublowerRight = [substrate_length , 0];
SubupperRight = [substrate_length , substrate_width];
SubupperLeft =  [0 , substrate_width];

Sub = [3;4;SublowerLeft(1);SublowerRight(1);SubupperRight(1);SubupperLeft(1); ...
    SublowerLeft(2);SublowerRight(2);SubupperRight(2);SubupperLeft(2)];
CAF = [1;substrate_length/2;substrate_length/2;CAF_radius];
CAF = [CAF;zeros(length(Sub) - length(CAF),1)];

%C1 = [1;0;0;0.4];
%C1 = [C1;zeros(length(R1) - length(C1),1)];

gd = Sub;
sf = 'Sub';
ns = char('Sub')'; %OCCHIO ALL'APICE
gSub = decsg(gd,sf,ns);

gd = [Sub, CAF];

sf = 'Sub+CAF';
ns = char('Sub','CAF')'; %OCCHIO ALL'APICE
g = decsg(gd,sf,ns);

figure
pdegplot(g,'FaceLabels','on');
%[1,2],[-2,0.5])

%g = addFace(g,[9 5 3 7]);
%g = addFace(g,[1 3 8 4]);
%g = addFace(g,[2 7 10 6]);
tmodel = createpde('thermal', 'transient');
g = geometryFromEdges(tmodel,g);

figure;
g = extrude(g,substrate_thickness); % estrudo substrato
tmodel.Geometry = g;
%pdegplot(tmodel,'CellLabels','on','EdgeLabel','on','FaceLabels','on');


g = extrude(g,4,CAF_thickness); % estrudo CAF di CAF_thickness

%pdegplot(tmodel,'CellLabels','on','EdgeLabel','on','FaceLabels','on');

tmodel.Geometry = g;
figure
pdegplot(tmodel,'CellLabels','on','EdgeLabel','off','FaceLabels','on');

figure
generateMesh(tmodel, 'Hmax',h);
pdemesh(tmodel);

%ATTENZIONE: FORSE MATLAB USA GLI inc INVECE CHE I metri. Convertire le
%conducibilità termiche!

thermalProperties(tmodel, "Cell",[1 2], "ThermalConductivity",kvetro,"MassDensity",rhoVetro,"SpecificHeat",cVetro); %substrato vetro
thermalProperties(tmodel, "Cell",3, "ThermalConductivity",kcaf,"MassDensity",rhoCAF,"SpecificHeat",cCAF); %CAF

%condizoini iniziali
thermalIC(tmodel,300);

internalHeatSource(tmodel, q, "Cell",3);
thermalBC(tmodel, 'Face',[1,2],'Temperature',300); %sotto
thermalBC(tmodel, 'Face',[5,6,7,8],'Temperature',300); %lati
%thermalBC(tmodel, 'Face',[57,43,46,30],'Temperature',300); %lati elettrodi

%thermalIC(tmodel,0); non necessario, siamo nel caso stazionario
tmodel.StefanBoltzmannConstant = 5.670367e-8; %SI units

result = solve(tmodel,times);
%result = solvepde(tmodel);
T1 = result.Temperature(:,1);
T2 = result.Temperature(:,2);
figure
pdeplot3D(tmodel,ColorMapData=T1);
figure
pdeplot3D(tmodel,ColorMapData=T2);

%results al tempo finale
% results2 = results;
% results2.Temperature = results2.Temperature(:,2);

[X,Y,Z] = meshgrid(0:0.05:substrate_length,0:0.05:substrate_width,0:0.05:substrate_thickness);
%V = interpolateSolution(result,X,Y,Z);
T_3D = interpolateTemperature(result,X,Y,Z, 2); %evaluate at timestep 2
T_3D = reshape(T_3D,size(X));

figure
slice(X,Y,Z,T_3D,[],substrate_width/2,[substrate_thickness*0.5, substrate_thickness*1],'nearest');
colorbar;


% figure
% colormap jet
% contourslice(X,Y,Z,T_3D,[],[],0:0.02:substrate_thickness)
% xlabel("x")
% ylabel("y")
% zlabel("z")
% colorbar

% quiver plot
% [qx,qy,qz] = evaluateHeatFlux(result);
% 
% figure
% pdeplot3D(result.Mesh,FlowData=[qx qy qz])
% 
% slice plot
% x=tmodel.Mesh.Nodes(1,:)';       %  x-coordinate of nodes
% y=tmodel.Mesh.Nodes(2,:)';       %  y-coordinate of nodes
% z=tmodel.Mesh.Nodes(3,:)';       %  z-coordinate of nodes
% u= T(1,:)';        %  temperature of nodes
% 
% 
% figure
% slice(x,y,z,T,[1,0.5,1],[],[],'nearest')
% 
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