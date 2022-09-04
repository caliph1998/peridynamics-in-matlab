close all
clear all
%% PARAMETERS REQUIRING USER INPUT
radiusMajor = 0.05; % Unit: m central hole major radius
radiusMinor = 0.05; % Unit: m central hole minor radius
mpLinearDensity = 200; %material points per m
TimeInterval = 6000; %TimeInterval: Number of time iTimeIntervalervals



%% GEOMETRICAL PARAMETERS

length = 0.5; % Unit: m Plate length
width = 0.5; % Unit: m Plate width
ellipseCurvature = radiusMinor ^ 2 / radiusMajor;
numOfDivX = mpLinearDensity * length; %NumofDiv_x: Number of divisions in x direction
numOfDivY = numOfDivX;
dx = 1 / mpLinearDensity; %Incremental distance between a material points pair
thick = dx; % Unit: m thickness of the plate
area = dx * dx; % minimum rectanble area containing one material point
volume = area * thick; % unit: m^3 %minimum volume of a cube containing one material point
InitialTotalNumMatPoint = numOfDivX*numOfDivY; %Initial estimation of number of material points needed
%% MECHANICAL PROPERTIES

ELASTIC_MODULUS = 200e9; % Unit: N/m^2
POISSON_RATIO = 0.3;  
APPLIED_PRESSURE = 500e7; % Unit: N/m^2
SHEAR_MODULUS = ELASTIC_MODULUS / (2 * (1 + POISSON_RATIO)); 
BULK_MODULUS = ELASTIC_MODULUS / (2 * (1 - POISSON_RATIO)); 
%% PERIDYNAMICS PARAMETERS

delta = 3.015 * dx; % Unit: m horizon size: a peridynamics (PD) parameter
volumeCorrectionRadius = dx / 2; % volume correction related number
alpha=0.5*(BULK_MODULUS-2*SHEAR_MODULUS); % a PD parameter used ONLY in state-based PD
bcd = 2/(pi*thick*delta^3); % a PD parameter referenced in literature as d
bcs = 6*SHEAR_MODULUS/(pi*thick*delta^4); % a PD parameter referenced in literature as b
neighborsPerNode = floor(pi * delta ^ 2 / area);
totalNumOfBonds = InitialTotalNumMatPoint * neighborsPerNode * 100;
nodefam = zeros(totalNumOfBonds,3); % Total array allocated to storing the neighbors of every material point and their bondforces
%% OTHER PARAMTETERS

%TODO: TimeInterval variable should be automatically approximated by the
%model.
warning('TODO: TimeInterval variable should be automatically approximated by the model.');
dt = 1; % time unit (s)

%% DENSE MATERIAL POINT AREA PARAMTERS 

tip_radius = radiusMajor;
%dx_tip = tip_radius / tipNumOfDiv * 40; %Recommended material point sizes at the crack tips
dxDenseArea = dx / 2;
%dx_tip = 6.25e-04;
%% DEFINING THE TIPS AND THE ELLIPSE HOLE REGIONS USING CLASS ELLIPSECLASS

% Constructor: EllipseClass(x_origin, y_origin_, major radius, minor radius);
center_hole = EllipseClass(0, 0, radiusMajor, radiusMinor / 4);
left_tip = EllipseClass(0, 0, radiusMajor * 2 * 0 , radiusMajor * 2 * 0 ); %if major and minor radii are equal => circle.
right_tip = EllipseClass(radiusMajor, 0, 0, 0);
%% COORDINATE GENERATION FOR EACH MATERIAL POINT

%TODO:Build an interface for desired inputs.
warning('TODO:Build an interface for desired inputs.');
[coord, path_horizontal, numOfRemoteMPs, totalNumMatPoint] = generateBarDense(numOfDivX, length, width, dx, dxDenseArea, InitialTotalNumMatPoint, center_hole);
%[coord, path_horizontal, numOfRemoteMPs, totalNumMatPoint] = generateCircleDense(NumofDiv_x, length, width, dx, dxDense, InitialTotalNumMatPoint, center_hole, left_tip, right_tip);


%% HORIZON SIZE ASSIGNEMNTS

Deltas = zeros(totalNumMatPoint, 1); %an array holding all material points horizons
BCS = zeros(totalNumMatPoint, 1); % an array holding bcs of all material points
BCD = zeros(totalNumMatPoint, 1); % an array holding bcd of al material points

for i = 1: numOfRemoteMPs
    Deltas(i, 1) = delta;
end

for i = numOfRemoteMPs+1: totalNumMatPoint
    Deltas(i, 1) = 3.015 * dxDenseArea;
end

for i = 1:totalNumMatPoint
    BCS(i, 1) = 6*SHEAR_MODULUS/(pi*thick*Deltas(i,1)^4);
    BCD(i, 1) = 2/(pi*thick*Deltas(i,1)^3);
end

%% DEFINITION OF NEEDED ARRAYS (COORD_EXCESS IS TRIMMED AT THIS STAGE)

path_horizontal = SortPath('hor',path_horizontal, coord); %Sorting horizontal path from head to tail of the path
numfam = zeros(totalNumMatPoint,1); %numfam: Number of family nodes
pointfam = zeros(totalNumMatPoint,1); %pointfam: Pointer
PDforce = zeros(totalNumMatPoint,2);%PDforce: Peridynamic force
BodyForce = zeros(totalNumMatPoint,2);%Body force
PDforceold = zeros(totalNumMatPoint,2);
PD_SED_distorsion = zeros(totalNumMatPoint,2);
SurCorrFactor_dilatation = zeros(totalNumMatPoint,2); % PD surface correction factor for dilatation
SurCorrFactor_distorsion = zeros(totalNumMatPoint,2); % PD surface correction factor old for dilatation
disp = zeros(totalNumMatPoint,2); % displacement
total_disp = zeros(totalNumMatPoint,2); %The total sum of displacement for each material point
vel = zeros(totalNumMatPoint,2); % velocity
velhalfold = zeros(totalNumMatPoint,2);% velocity of half old
velhalf = zeros(totalNumMatPoint,2); % velocity of half
massvec = zeros(totalNumMatPoint,2);% mass vector
PD_SED_dilatation = zeros(totalNumMatPoint,2); % Peridynamic strain energy density for dilatation
PD_SED_dilatation_Fixed = zeros(totalNumMatPoint,2); % Fixed Peridynamic strain energy density for dilatation
Check_time = zeros(TimeInterval,1);
Steady_check_x = zeros(TimeInterval,1);
Steady_check_y = zeros(TimeInterval,1);
fprintf('Large arrays assigned. \n');
%path_edge = SortPath('edge',path_edge, coord);
%path_circular = SortPath('circle',path_circular, coord);



%% NEIGHBOR MEMBER SEARCH ALGORITHM
tic
for i = 1:totalNumMatPoint
    
    delta = Deltas(i,1);
    if (i == 1) 
    pointfam(i,1) = 1;
    else
    pointfam(i,1) = pointfam(i-1,1) + numfam(i-1,1);
    end

    for j = 1:totalNumMatPoint
        RelativePosition_Vector = sqrt((coord(j,1) - coord(i,1))^2 + (coord(j,2) - coord(i,2))^2);
        if(i~=j)
            volCorr_radius = Deltas(j,1) / 3.015 / 2;
            if(RelativePosition_Vector <= delta )%+ volCorr_radius)
            numfam(i,1) = numfam(i,1) + 1;
            nodefam(pointfam(i,1)+numfam(i,1)-1,1) = j;
            end
        end
    end
end
toc
fprintf('\nNeighbors assigned. \n');
lastBondIter = pointfam(i,1) + numfam(i,1) - 1; % Used to trim nodefam and to define bondForces array
%nodefam = nodefam(1:lastBondIter, 1); %Trimming the nodefame to the last non-zero value.
bondForces = zeros(lastBondIter, 1);

%% SURFACE CURRACTION FACTORING

%Surface correction factor calculation - start
for i = 1:totalNumMatPoint
disp(i,1) = 0.001*coord(i,1);
disp(i,2) = 0.0;
end

%%Surface correction factor calculation - function%%
[PD_SED_distorsion(:,1),SurCorrFactor_distorsion(:,1),PD_SED_dilatation(:,1), SurCorrFactor_dilatation(:,1)] = Calculate_SurCorrection(Deltas,BCS,BCD,disp,totalNumMatPoint,numfam,nodefam,pointfam,coord);
fprintf('\nPD_SEDs  in x direction assigned\n');
    
for i = 1:totalNumMatPoint
    disp(i,1) = 0;
    disp(i,2) = 0.001 * coord(i,2);
end
[PD_SED_distorsion(:,2), SurCorrFactor_distorsion(:,2),PD_SED_dilatation(:,2), SurCorrFactor_dilatation(:,2)] = Calculate_SurCorrection(Deltas,BCS,BCD,disp,totalNumMatPoint,numfam,nodefam,pointfam,coord);
%Surface correction factor calculation - end
fprintf('\nPD_SEDs  in y direction assigned\n');
    
%initial displacemeTimeInterval
for i = 1:totalNumMatPoint
disp(i,1) = 0;
disp(i,2) = 0; 
end


%Stable mass vector computation
for i = 1:totalNumMatPoint
massvec(i,1) = 1.25 * dt ^ 2 * (pi * (Deltas(i,1))^2 * thick)  * BCS(i,1) / (Deltas(i,1) / 3.015);
massvec(i,2) = 1.25 * dt ^ 2 * (pi * (Deltas(i,1))^2 * thick) * BCS(i,1) / (Deltas(i,1) / 3.015);
end

%% APPLYING EXTRENAL FORCES

for i = 1:totalNumMatPoint
    if (coord(i,2) == min(coord(:,2))) %applying force to the lower edge
        BodyForce(i,2) = (-1) * APPLIED_PRESSURE/dx;
    elseif(coord(i,2) == max(coord(:,2)))%applying force to the upper edge
        BodyForce(i,2) =  APPLIED_PRESSURE/dx; 
    end
end
%{
%Applied loading - Left
for i = 1:NumofDiv_y
BodyForce(i,1) = -1 * Applied_pressure/dx;% * (coord(i,2) * 2 + 0.5);
end


%Applied loading - Right
for i = (nnum-NumofDiv_y+1):nnum
BodyForce(i,1) = Applied_pressure/dx;% * (coord(i,2) * 2 + 0.5);
end
%}
%{
for i = 1:TotalNumMatPoint
    if (coord(i,1) == min(coord(:,1)))
        BodyForce(i,1) = -1 * Applied_pressure/dx;
    elseif(coord(i,1) == max(coord(:,1)))
        BodyForce(i,1) = Applied_pressure/dx;
    end
end
%}
%{
%Applied loading - Left
for i = 1:NumofDiv_y
BodyForce(i,1) = -1 * Applied_pressure/dx;
end

%Applied loading - Right
for i = (nnum-NumofDiv_y+1):nnum
BodyForce(i,1) = Applied_pressure/dx;
end
%}

%% TIME INTEGRATION STEP

%testNode = 555;
testCoordinates = [tip_radius * 6 / 5, 0];
testNode = get_closest_point(testCoordinates(1,1), testCoordinates(1,2), coord);


%One of these lists is redundant. The reason is that for every direct
%force, there is one reactionary force in equal magnitude in the opposite
%direction in the other list with even the same array index.
directForces = zeros(totalNumOfBonds,4);
reactForces = zeros(totalNumOfBonds,4);
iterDirect = 1;
iterReact = 1;

tic
for tt = 1:TimeInterval
time = tt
    
    for i = 1:totalNumMatPoint
    PD_SED_dilatation_Fixed(i,1) = 0;
        for j = 1:numfam(i,1)
         cnode = nodefam(pointfam(i,1)+j-1,1);
         RelativePosition_Vector = sqrt((coord(cnode,1) - coord(i,1))^2 + (coord(cnode,2) - coord(i,2))^2);
         RelativeDisp_Vector=sqrt((coord(cnode,1)+disp(cnode,1)-coord(i,1)-disp(i,1))^2+(coord(cnode,2)+disp(cnode,2)-coord(i,2)-disp(i,2))^2);
         Stretch = (RelativeDisp_Vector - RelativePosition_Vector) / RelativePosition_Vector;
         AbsoluteValue_x_y = RelativeDisp_Vector * RelativePosition_Vector;
         Coeff_x = (coord(cnode,1) + disp(cnode,1) - coord(i,1) - disp(i,1)) * (coord(cnode,1) - coord(i,1));
         Coeff_y = (coord(cnode,2) + disp(cnode,2) - coord(i,2) - disp(i,2)) * (coord(cnode,2) - coord(i,2));
         Directional_cosine = (Coeff_x + Coeff_y) / AbsoluteValue_x_y;
         volumeCorrectionRadius = Deltas(cnode,1) / 3.015 / 2;
            if (RelativePosition_Vector <= Deltas(i,1)-volumeCorrectionRadius) 
            fac = 1;
            elseif (RelativePosition_Vector <= Deltas(i,1)+volumeCorrectionRadius && RelativePosition_Vector > Deltas(i,1)-volumeCorrectionRadius )
            fac = (Deltas(i,1)+volumeCorrectionRadius-RelativePosition_Vector)/(2*volumeCorrectionRadius);
            else
            fac = 0;
            end

            if (abs(coord(cnode,2) - coord(i,2)) <= 1e-10)
            theta = 0;
            elseif (abs(coord(cnode,1) - coord(i,1)) <= 1e-10)
            theta = 90*pi/180;
            else
            theta = atan(abs(coord(cnode,2) - coord(i,2)) / abs(coord(cnode,1) - coord(i,1)));
            end

            SurCorrFactor_x = (SurCorrFactor_dilatation(i,1) + SurCorrFactor_dilatation(cnode,1))/2;
            SurCorrFactor_y = (SurCorrFactor_dilatation(i,2) + SurCorrFactor_dilatation(cnode,2))/2;
            SurCorrFactor_Arbitrary_dilatation = 1 / (((cos(theta))^2/(SurCorrFactor_x)^2) + ((sin(theta))^2 / (SurCorrFactor_y)^2));
            SurCorrFactor_Arbitrary_dilatation = sqrt(SurCorrFactor_Arbitrary_dilatation);

            %Critic stretch if statement is deactivated as it is unused
            %if (failm(i,j)==1) 
            volume = (Deltas(cnode,1) / 3.015) ^ 2 * thick;
            bcs = BCS(i,1);
            bcd = BCD(i,1);
            
            PD_SED_dilatation_Fixed(i,1) = PD_SED_dilatation_Fixed(i,1) + BCD(i,1) * Deltas(i,1) * Stretch * Directional_cosine * volume * SurCorrFactor_Arbitrary_dilatation * fac;                          
            %else
            %PD_SED_dilatation_Fixed(i,1) = 0;
            %end                                             
        end
    end

PDforce = zeros(totalNumMatPoint, 2);
nodefam(:,2) = 0;
nodefam(:,3) = 0;
for i = 1:totalNumMatPoint
    for j = 1:numfam(i,1)
    cnode = nodefam(pointfam(i,1)+j-1,1); %cnode: the current neighbor node/material-point.
    RelativePosition_Vector = sqrt((coord(cnode,1) - coord(i,1))^2 + (coord(cnode,2) - coord(i,2))^2); %the initial distance
    RelativeDisp_Vector=sqrt((coord(cnode,1)+disp(cnode,1)-coord(i,1)-disp(i,1))^2+(coord(cnode,2)+disp(cnode,2)-coord(i,2)-disp(i,2))^2);
    Stretch = (RelativeDisp_Vector - RelativePosition_Vector) / RelativePosition_Vector;
    AbsoluteValue_x_y = RelativeDisp_Vector * RelativePosition_Vector;
    Coeff_x = (coord(cnode,1) + disp(cnode,1) - coord(i,1) - disp(i,1)) * (coord(cnode,1) - coord(i,1)); %(C'x - Cx)(Cx-Ax)
    Coeff_y = (coord(cnode,2) + disp(cnode,2) - coord(i,2) - disp(i,2)) * (coord(cnode,2) - coord(i,2)); %(C'y - Cy)(Cy-Ay)
        Directional_cosine = (Coeff_x + Coeff_y) / AbsoluteValue_x_y;%Some weird constant that will be used in bond_constant calculations.
        volumeCorrectionRadius = Deltas(cnode,1) / 3.015 / 2;
        if (RelativePosition_Vector <= Deltas(i,1)-volumeCorrectionRadius) %if all the way inside the horizon
         fac = 1;
        elseif (RelativePosition_Vector <= Deltas(i,1)+volumeCorrectionRadius) %if partially inside the A horizon
         fac = (Deltas(i,1)+volumeCorrectionRadius-RelativePosition_Vector)/(2*volumeCorrectionRadius); %VolCorr_radius = dx / 2
        else
         fac = 0; %Unnecassary else since it will never happen
         error("fac2");
        end
        if (abs(coord(cnode,2) - coord(i,2)) <= 1.0e-10) %if both C and A were on a horizontal line
         theta = 0;
        elseif (abs(coord(cnode,1) - coord(i,1)) <= 1.0e-10) %if both C and A were on a vertical line
         theta = 90 * pi / 180;
        else
         theta = atan(abs(coord(cnode,2) - coord(i,2)) / abs(coord(cnode,1) - coord(i,1))); %theta: C and A angle
        end

        SurCorrFactor_x = (SurCorrFactor_distorsion(i,1) + SurCorrFactor_distorsion(cnode,1)) / 2;
        SurCorrFactor_y = (SurCorrFactor_distorsion(i,2) + SurCorrFactor_distorsion(cnode,2)) / 2;
        SurCorrFactor_Arbitrary_distorsion = 1 / (((cos(theta))^2/(SurCorrFactor_x)^2) + ((sin(theta))^2/(SurCorrFactor_y)^2));
        SurCorrFactor_Arbitrary_distorsion = sqrt(SurCorrFactor_Arbitrary_distorsion);
        SurCorrFactor_x = (SurCorrFactor_dilatation(i,1) + SurCorrFactor_dilatation(cnode,1)) / 2;
        SurCorrFactor_y = (SurCorrFactor_dilatation(i,2) + SurCorrFactor_dilatation(cnode,2)) / 2;
        SurCorrFactor_Arbitrary_dilatation = 1 / (((cos(theta))^2 / (SurCorrFactor_x)^2) + ((sin(theta))^2 / (SurCorrFactor_y)^2));
        SurCorrFactor_Arbitrary_dilatation = sqrt(SurCorrFactor_Arbitrary_dilatation);
        %Note: same variable SurCorrFactor_x and SurCorrFactor_y are used
        %to hold two different concepts as variables.
% bcd = d in Madenci; bcs: b in Madenci; PD_SED_dilatation_Fixed: thetha_k
% in madenci; SurCorrFactor_Arbitrary_dilation: Gd in Madenci;
% SurCorrFactor_Arbitrary_distorction: Gb in Madenci; Directional_cosine:
% Arrowhead_k_j in Madenci; alpha: a in madenci
%So here it seems t_kj and t_jk are summed up to form bonForce_const
%instead of using their difference.
%Anyway, in order to apply the dual horizon, you need to use t_kj of all
%the neighbors BUT t_jk of only those who are in your dual horizon (i.e.
%recognize you as well).
%For points with smaller horizon, this doesn't make a difference, but it
%does make a difference for the points with the larger horizon at the
%vicinity of the smaller horzion points.
        volume = (Deltas(cnode,1) / 3.015) ^ 2 * thick;
        bondForce_const = (2 * BCD(i,1)*Deltas(i,1) * alpha / RelativePosition_Vector * Directional_cosine * PD_SED_dilatation_Fixed(i,1)  + ...
                      2 * BCS(i,1)*Deltas(i,1) * Stretch * SurCorrFactor_Arbitrary_distorsion) * volume * fac / RelativeDisp_Vector;
        %The point of the two lines below is to split the force vector into its components.
        %But it needs to be divided by RelativeDisp_Vector which is done in
        %the above line which is unclear and confusing.
        directForce_x =  (coord(cnode,1) + disp(cnode,1) - coord(i,1) - disp(i,1)) * bondForce_const ;           
        directForce_y =  (coord(cnode,2) + disp(cnode,2) - coord(i,2) - disp(i,2)) * bondForce_const ;           
        
        PDforce(i,1) = PDforce(i,1) + directForce_x;     
        PDforce(i,2) = PDforce(i,2) + directForce_y;
        
        
        reactionaryForce_x = directForce_x / (volume * fac);
        reactionaryForce_y = directForce_y / (volume * fac);
        
        
        volume = (Deltas(i,1) / 3.015) ^ 2 * thick;
        reactionaryForce_x = reactionaryForce_x * (volume * fac);
        reactionaryForce_y = reactionaryForce_y * (volume * fac);
        


        
        
        reactionaryForce_x = (-1) * reactionaryForce_x;
        reactionaryForce_y = (-1) * reactionaryForce_y;
        PDforce(cnode,1) = PDforce(cnode,1) + reactionaryForce_x;     
        PDforce(cnode,2) = PDforce(cnode,2) + reactionaryForce_y;
        
        

        
 
       if (tt == TimeInterval)
        directForces(iterDirect, 1) = i;
        directForces(iterDirect, 2) = cnode;
        directForces(iterDirect, 3) = directForce_x;
        directForces(iterDirect, 4) = directForce_y;
        iterDirect = iterDirect + 1;
        
        reactForces(iterReact, 1) = cnode;
        reactForces(iterReact, 2) = i;
        reactForces(iterReact, 3) = reactionaryForce_x;
        reactForces(iterReact, 4) = reactionaryForce_y;
        iterReact = iterReact + 1;
        nodefam(pointfam(i,1)+j-1, 2) = nodefam(pointfam(i,1)+j-1, 2) +  directForce_x;
        nodefam(pointfam(i,1)+j-1, 3) = nodefam(pointfam(i,1)+j-1, 3) +  directForce_y;
        
           for iter = pointfam(cnode, 1): pointfam(cnode,1) + numfam(cnode, 1) - 1
                neighb = nodefam(iter, 1);
              if (neighb ~= i)
                 continue;
              end
         i_index = iter;
         nodefam(i_index, 2) = nodefam(i_index, 2) + reactionaryForce_x;
         nodefam(i_index, 3) = nodefam(i_index, 3) + reactionaryForce_y;
           end
       end    
           
    end
end
%% Cutting the force lists to remove excess rows
directForces = directForces(1:iterDirect-1,:);
reactForces = reactForces(1:iterReact-1,:);

%%
%%% Getting steady-state solutions through Adaptive Dynamic relaxation %%%
cn1 = 0;
cn2 = 0;
for i = 1:totalNumMatPoint
    if (velhalfold(i,1)~=0)
    cn1 = cn1 - disp(i,1) * disp(i,1) * (PDforce(i,1) / massvec(i,1) - PDforceold(i,1) / massvec(i,1)) / (dt * velhalfold(i,1));
    end
    
    if (velhalfold(i,2)~=0) 
    cn1 = cn1 - disp(i,2) * disp(i,2) * (PDforce(i,2) / massvec(i,2) - PDforceold(i,2) / massvec(i,2)) / (dt * velhalfold(i,2));
    end
    
    cn2 = cn2 + disp(i,1) * disp(i,1);
    cn2 = cn2 + disp(i,2) * disp(i,2);
end

    if (cn2~=0)
        if ((cn1 / cn2) > 0) 
        cn = 2 * sqrt(cn1 / cn2);
        else
        cn = 0;
        end
    else
    cn = 0;
    end

    %if (cn > 2)
    %cn = 1.9;
    %end

    for i = 1:totalNumMatPoint
        if (tt == 1)
        velhalf(i,1) = 1 * dt / massvec(i,1) * (PDforce(i,1) + BodyForce(i,1)) / 2;		
        velhalf(i,2) = 1 * dt / massvec(i,2) * (PDforce(i,2) + BodyForce(i,2)) / 2;
        else	
        velhalf(i,1) = ((2 - cn * dt) * velhalfold(i,1) + 2 * dt / massvec(i,1) * (PDforce(i,1) + BodyForce(i,1))) / (2 + cn * dt);
        velhalf(i,2) = ((2 - cn * dt) * velhalfold(i,2) + 2 * dt / massvec(i,2) * (PDforce(i,2) + BodyForce(i,2))) / (2 + cn * dt);
        end
        
        %%%Deactivated velocity calculation that had no use%%%
        vel(i,1) = 0.5 * (velhalfold(i,1) + velhalf(i,1));
        vel(i,2) = 0.5 * (velhalfold(i,2) + velhalf(i,2));
        disp(i,1) = disp(i,1) + velhalf(i,1) * dt;
        disp(i,2) = disp(i,2) + velhalf(i,2) * dt;
        total_disp(i,1) = total_disp(i,1) + disp(i,1);
        total_disp(i,2) = total_disp(i,2) + disp(i,2);       
        velhalfold(i,1) = velhalf(i,1);
        velhalfold(i,2) = velhalf(i,2);
        PDforceold(i,1) = PDforce(i,1);
        PDforceold(i,2) = PDforce(i,2);
    end

Check_time(tt,1)= tt;
Steady_check_x(tt,1) = disp(testNode,1);
Steady_check_y(tt,1) = disp(testNode,2);
end
toc
%time interval iteration ends

% 50*40 +(25+12) = 2037
%%
%Dongjun_hole_stress = CalculateStressforPoint(coord,totalNumMatPoint,numfam,nodefam, pointfam, Deltas);
unpunched_d = coord(path_horizontal(1),1) - coord(path_horizontal(end),1); %Distance from the crack tip to the plate edge
unpunched_d = unpunched_d * (-1);
normal_path_horizontal = (coord(path_horizontal,1) - coord(path_horizontal(1),1)) / (unpunched_d); %normalized path distance
%%
%Diffenece with stressE: only neighbors are considered if below.
stressE = zeros(totalNumMatPoint, 1);
for i = 1 : size(directForces, 1)
    if (coord(directForces(i,2), 2) > coord(directForces(i,1), 2))
        stressE(directForces(i,1),1) = stressE(directForces(i,1),1) + directForces(i,4);
    end
end
for i = 1 : size(reactForces, 1)
    if (coord(reactForces(i,2), 2) > coord(reactForces(i,1), 2))
        stressE(reactForces(i,1), 1) = stressE(reactForces(i,1),1) + reactForces(i,4);
    end
end


for i = 1 : totalNumMatPoint
    for l = pointfam(i, 1) : pointfam(i, 1) + numfam(i, 1) - 1
        j = nodefam(l,1);
        if (coord(j,2) < coord(i,2))
            if (coord(j,1) == coord(i,1))
                for k = 1 : size(directForces, 1)
                    if(directForces(k,1) == j && coord(directForces(k,2), 2) > coord(i, 2))
                        stressE(i, 1) = stressE(i,1) + directForces(k,4);
                    end
                    if(reactForces(k,1) == j && coord(reactForces(k,2), 2) > coord(i,2))
                        stressE(i, 1) = stressE(i,1) + reactForces(k,4);
                    end
                end
            %note to future self: remember to use get_closest_point for
            %stress calculation in stressE and stressD  even.
            elseif (Deltas(j,1) > Deltas(i,1) && coord(j,1) < coord(i,1))
                rPortion = (coord(i,1) - coord(j,1)) / (Deltas(j,1) / 3.015);
                lPortion = 1 - rPortion;
                rMP = get_closest_point(coord(j,1) + Deltas(j,1) / 3.015, coord(j,2), coord);
                for k = 1 : size(directForces, 1) %could use switch here
                    if(directForces(k,1) == j && coord(directForces(k,2), 2) > coord(i, 2))
                        stressE(i, 1) = stressE(i,1) + directForces(k,4) * lPortion;
                    end
%                     if(reactForces(k,1) == j && coord(reactForces(k,2), 2) > coord(i,2))
%                         stressE(i, 1) = stressE(i,1) + reactForces(k,4) * lPortion;
%                     end
                    if(directForces(k,2) == j && coord(reactForces(k,1), 2) > coord(i,2))
                        stressE(i,1) = stressE(i,1) - directForces(k,4) * lPortion;
                    end
                    if(directForces(k,1)                  == rMP && coord(directForces(k,2), 2) > coord(i, 2))
                        stressE(i, 1) = stressE(i,1) + directForces(k,4) * rPortion;
                    end
%                     if(reactForces(k,1) == rMP && coord(reactForces(k,2), 2) > coord(i,2))
%                         stressE(i, 1) = stressE(i,1) + reactForces(k,4) * rPortion;
%                     end
                    if(directForces(k,2) == rMP && coord(directForces(k,1), 2) > coord(i,2))
                        stressE(i,1) = stressE(i,1) - directForces(k,4) * rPortion;
                    end
                end
            end
        end
    end
end

for node = 1: totalNumMatPoint
    T = Deltas(node,1) / 3.015;
    stressE(node,1) = stressE(node,1) * T / APPLIED_PRESSURE;
end
%% GETTING THE SAVING ADDRESS FROM THE USER

while 1 == 1
    address = 0;
    prompt = 'Where to save files? Enter the address. Press Enter if no save is needed:\n';
    addressOld = input(prompt, 's');
    if (isempty(addressOld) == 0 && isfolder(addressOld) == 0)
     fprintf('\ninvalid address. This folder does not exist. Please try again\n');
     continue;
    end
    prompt = 'For security, please re-enter the address. Press Enter if no save is needed:\n';
    addressRepeat = input(prompt, 's');
    if (strcmp(addressOld, addressRepeat) == 0)
        fprintf('\nRepeated address does not match with the previous one. Try again\n');
        continue;
    end
    address = addressRepeat;
    break;
end
%% SAVING THE WORKSPACE

if (isempty(address) == 0)
    save(strcat(address, '\workspace'));
end
%% DEFORMED VS UNDEFORMED FIGURE

figure(1)
hold on
figureNameDeformation = '\deformation.jpeg';
h1=plot(coord(:,1),coord(:,2),'.r');
h2=plot(coord(:,1)+disp(:,1),coord(:,2)+disp(:,2),'.b');
h3=plot(coord(testNode,1),coord(testNode,2),'ro','MarkerSize',3.3,'MarkerFaceColor','r');
h4=plot(coord(testNode,1)+disp(testNode,1),coord(testNode,2)+disp(testNode,2),'bo','MarkerSize',3.3,'MarkerFaceColor','b');
legend([h1 h2 h3 h4],{'Undeformed state','Deformed state','Check point in the undeformed state',...
'Check point in the deformed state'});
xlim([-0.4 0.4])
ylim([-0.4 0.4])
xlabel('x axis [m]');
ylabel('y axis[m]');


saveAndPrintFigure(gcf, address, figureNameDeformation);
%% DISPLACEMENT FIELD
figure(2)
figureNameDisp = '\dispXXYY.jpeg';
sz = 10;
subplot(1,2,1);
%plotting the absolute values of displacements%
scatter(coord(:,1), coord(:,2), sz, -abs(disp(:,1)), 'filled');
xlabel('x');
ylabel('y');
title(['U11 in ', num2str(numOfDivX), ' * ', num2str(numOfDivY)]);
colorbar('southoutside');
colormap('jet');

subplot(1,2,2);
scatter(coord(:,1), coord(:,2), sz, abs(disp(:,2)), 'filled');
xlabel('x');
ylabel('y');
title(['U22 in ', num2str(numOfDivX), ' * ', num2str(numOfDivY)]);
colorbar('southoutside');
colormap('jet');

saveAndPrintFigure(gcf, address, figureNameDisp);
% %% STRESS FIELD
% figure(3)
% %figureNameStress = '\stressXXYY(INCORRECT).jpeg';
% sz = 10;
% subplot(1,2,1);
% scatter(coord(:,1), coord(:,2), sz, (Dongjun_hole_stress(:,1)), 'filled');
% xlabel('x');
% ylabel('y');
% colorbar('southoutside');
% colormap('jet');
% title(['S11 in ', num2str(NumofDiv_x), ' * ', num2str(NumofDiv_y)]);
% subplot(1,2,2);
% scatter(coord(:,1), coord(:,2), sz, (Dongjun_hole_stress(:,2)), 'filled');
% xlabel('x');
% ylabel('y');
% colorbar('southoutside');
% colormap('jet');
% title(['S22 in ', num2str(NumofDiv_x), ' * ', num2str(NumofDiv_y)]);
% 
% %saveAndPrintFigure(gcf, address, figureNameStress);
% %% NORMAL STRESS FIELD
% figure(4)
% figureNameStress = '\stressXXYYNormal(INCORRECT).jpeg';
% sz = 10;
% subplot(1,2,1);
% scatter(coord(:,1), coord(:,2), sz, (Dongjun_hole_stress(:,1)) / Applied_pressure, 'filled');
% xlabel('x');
% ylabel('y');
% colorbar('southoutside');
% colormap('volacno');
% title(['S11 in ', num2str(NumofDiv_x), ' * ', num2str(NumofDiv_y)]);
% subplot(1,2,2);
% scatter(coord(:,1), coord(:,2), sz, (Dongjun_hole_stress(:,2)) / Applied_pressure, 'filled');
% xlabel('x');
% ylabel('y');
% colorbar('southoutside');
% colormap('jet');
% title(['S22 in ', num2str(NumofDiv_x), ' * ', num2str(NumofDiv_y)]);
% 
% saveAndPrintFigure(gcf, address, figureNameStress);
% %% DONGJUN: STRESS vs NORMALIZED DISTANCE
% figure(5)
% %figureNameStressPathNorm = '\stressPathNorm.jpeg';
% subplot(1,2,1);
% ssx = ExtractPathData(path_horizontal, Dongjun_hole_stress, 1);
% plot(normal_path_horizontal, ssx);
% xlabel('Material Points');
% ylabel('S');
% title('S11 in the horizontal edge');
% subplot(1,2,2);
% ssy = ExtractPathData(path_horizontal, Dongjun_hole_stress, 2);
% plot(normal_path_horizontal, ssy);
% xlabel('Material Points');
% ylabel('S');
% title('S22 in the horizontal edge');
% %saveAndPrintFigure(gcf, address, figureNameStressPathNorm);
% %% DONGJUN: STRESS vs MATERIAL POINTS
% figure(6)
% %figureNameStressPath = '\stressPath.jpeg';
% subplot(1,2,1);
% plot(ExtractPathData(path_horizontal, Dongjun_hole_stress, 1));
% xlabel('Material Points');
% ylabel('S');
% title('S11 in the horizontal edge');
% subplot(1,2,2);
% plot(ExtractPathData(path_horizontal, Dongjun_hole_stress, 2));
% xlabel('Material Points');
% ylabel('S');
% title('S22 in the horizontal edge');
% %saveAndPrintFigure(gcf, address, figureNameStressPath);
%% DISPLACEMENT vs MATERIAL POINTS
figure(7)
figureNameDispPath = '\dispPath.jpeg';
subplot(1,2,1);
plot((ExtractPathData(path_horizontal, disp, 1)));
xlabel('Material Points');
ylabel('U');
title('U11 in the horizontal edge');
subplot(1,2,2);
plot(abs(ExtractPathData(path_horizontal, disp, 2)));
xlabel('Material Points');
ylabel('U');
title('U22 in the horizontal edge');

saveAndPrintFigure(gcf, address, figureNameDispPath);
%% STEADY STATE CHECK
figure(8)
figureNameStability = '\stability.jpeg';
hold on
h1=plot(Check_time(:,1),Steady_check_x(:,1),'.k');
h2=plot(Check_time(:,1),Steady_check_y(:,1),'.g');
legend([h1 h2],{'Displacement of x direction at blue point','Displacement of y direction at blue point'});
%ylim([-0.01 0.01])
title({'Steady state checking'});
xlabel('Time');
ylabel('Displacement [m]');

saveAndPrintFigure(gcf, address, figureNameStability);
%% HORIZONTAL PATH DEMONSTRATION

figure(9)
figureNamePath = '\pathDemonstration.jpeg';
hold on
scatter(coord(:,1), coord(:,2), '.g');
sz = 10;
scatter(ExtractPathData(path_horizontal, coord, 1), ExtractPathData(path_horizontal, coord, 2), sz + 10, '.b');
hold off

saveAndPrintFigure(gcf, address, figureNamePath);
%% DEMONSTRATION OF NEIGHBORS

node = get_closest_point(0.12,-0.135,coord);
sz = 50;
figure(10)
figureNameNeighbors = '\arrangement.jpeg';
hold on
xlabel('x');
ylabel('y');
title(['Neighbors of Material Point = ', num2str(node), ' in MP density of ' ,num2str(numOfDivX), '*', num2str(numOfDivX)]);
scatter(coord(:,1), coord(:,2), '.g');
for j = 1:numfam(node,1)
   cnode = nodefam(pointfam(node,1)+j-1,1);
   scatter(coord(cnode, 1), coord(cnode, 2), sz, '.r');
end

scatter(coord(node, 1), coord(node, 2), sz, '.b');
hold off

saveAndPrintFigure(gcf, address, figureNameNeighbors);
%%
% stressD = zeros(totalNumMatPoint, 1);
% stressA = zeros(totalNumMatPoint, 1);
% for testNode = 1:totalNumMatPoint
% fyA = 0;
% fyD = 0;
% for i = pointfam(testNode, 1) : pointfam(testNode, 1) + numfam(testNode, 1) - 1
%     neighbor = nodefam(i, 1);
%     if (coord(neighbor,2) > coord(testNode,2))
%         fyA = fyA + (nodefam(i,3));
%         fyD = fyD + (nodefam(i,3));
% 
%     elseif (coord(neighbor,1) == coord(testNode,1))
%         for j = pointfam(neighbor,1) : pointfam(neighbor, 1) + numfam(neighbor, 1) - 1
%             neighborOfNeighbor = nodefam(j, 1);
%             if (isBetween(coord, neighbor, neighborOfNeighbor, testNode) == 1)
%             %if (coord(neighborOfNeighbor,2) > coord(testNode,2))
%                 fyA = fyA + (nodefam(j,3));
%             end
%             if (coord(neighborOfNeighbor,2) > coord(testNode,2))
%                 fyD = fyD + (nodefam(j,3));
%             end
%         end
%     end
% end
% T = Deltas(testNode,1) / 3.015;
% stressA(testNode,1) = fyA * T / Applied_pressure;
% stressD(testNode,1) = fyD * T / Applied_pressure;
% end
%%
% All the forces in the y direction applied to All the material points with
% x= x0 and y <= y0 is sigmaY for x0,y0
% Wrong: all NEIGHBOR forces to all material points with x = x0 and y <=y0
% Wrong: all forces to all NEIGHBOR material points to MP x = x0 when y <= y0;
% i.e. don't bother with neighbor stuff. just sum the forces up.
% we expect the problem of no MP being x=x0 below you to still persist
% here, causing zebra-pattern stress field at the boundaries between the
% dense and the remote areas for the dense MPs.

% stressC = zeros(totalNumMatPoint, 1);
% for i = 1 : size(directForces, 1)
%     if (coord(directForces(i,2), 2) > coord(directForces(i,1), 2))
%         stressC(directForces(i,1),1) = stressC(directForces(i,1),1) + directForces(i,4);
%     end
% end
% for i = 1 : size(reactForces, 1)
%     if (coord(reactForces(i,2), 2) > coord(reactForces(i,1), 2))
%         stressC(reactForces(i,1), 1) = stressC(reactForces(i,1),1) + reactForces(i,4);
%     end
% end
% 
% 
% for i = 1 : totalNumMatPoint
%     for j = 1: totalNumMatPoint
%         if (coord(j,2) < coord(i,2))
%             if(coord(j,1) == coord(i,1))
%                 for k = 1 : size(directForces, 1)
%                     if(directForces(k,1) == j && coord(directForces(k,2), 2) > coord(i, 2))
%                         stressC(i, 1) = stressC(i,1) + directForces(k,4);
%                     end
%                     if(reactForces(k,1) == j && coord(reactForces(k,2), 2) > coord(i,2))
%                         stressC(i, 1) = stressC(i,1) + reactForces(k,4);
%                     end
%                 end
%             elseif (Deltas(j,1) > Deltas(i,1) && coord(j,1) < coord(i,1))
%                     rPortion = (coord(i,1) - coord(j,1)) / (Deltas(j,1) / 3.015);
%                     lPortion = 1 - rPortion;
%                     rMP = get_closest_point(coord(j,1) + Deltas(j,1) / 3.015, coord(j,2), coord);
%                     for k = 1 : size(directForces, 1) %could use switch here
%                         if(directForces(k,1) == j && coord(directForces(k,2), 2) > coord(i, 2))
%                             stressC(i, 1) = stressC(i,1) + directForces(k,4) * lPortion;
%                         end
%                         if(reactForces(k,1) == j && coord(reactForces(k,2), 2) > coord(i,2))
%                             stressC(i, 1) = stressC(i,1) + reactForces(k,4) * lPortion;
%                         end
%                         if(directForces(k,1) == rMP && coord(directForces(k,2), 2) > coord(i, 2))
%                             stressC(i, 1) = stressC(i,1) + directForces(k,4) * rPortion;
%                         end
%                         if(reactForces(k,1) == rMP && coord(reactForces(k,2), 2) > coord(i,2))
%                             stressC(i, 1) = stressC(i,1) + reactForces(k,4) * rPortion;
%                         end
%                     end
%             end
%         end
%     end
% end
% 
% for node = 1: totalNumMatPoint
%     T = Deltas(node,1) / 3.015;
%     stressC(node,1) = stressC(node,1) * T / Applied_pressure;
% end

    
%%
figure(997)
figureNameStress = '\stressXXYY.jpeg';

% subplot(1,3,1);
% sz = 10;
% scatter(coord(:,1), coord(:,2), sz, (stressD(:,1)), 'filled');
% xlabel('x');
% ylabel('y');
% colorbar('southoutside');
% colormap('jet');
% title(['S22 StressD in ', num2str(NumofDiv_x), ' * ', num2str(NumofDiv_y)]);
% subplot(1,3,2);
% sz = 10;
% scatter(coord(:,1), coord(:,2), sz, (stressC(:,1)), 'filled');
% xlabel('x');
% ylabel('y');
% colorbar('southoutside');
% colormap('jet');
% title(['S22 StressC in ', num2str(NumofDiv_x), ' * ', num2str(NumofDiv_y)]);
% subplot(1,3,3);
sz = 5;
scatter(coord(:,1), coord(:,2), sz, (stressE(:,1)), 'filled');
pbaspect([1 1 1])
xlabel('x');
ylabel('y');
colorbar('southoutside');
colormap('jet');
title(['S22 StressE in ', num2str(numOfDivX), ' * ', num2str(numOfDivY)]);

saveAndPrintFigure(gcf, address, figureNameStress);
%%
%DONGJUN: STRESS vs MATERIAL POINTS
figure(6)
figureNameStressPathNorm = '\stressPathNormYY.jpeg';
stressOnPath = ExtractPathData(path_horizontal, stressE, 1);
plot(normal_path_horizontal, stressOnPath);
xlabel('Material Points');
ylabel('S');
title('S22 stressE');
saveAndPrintFigure(gcf, address, figureNameStressPathNorm);
%%




















































function [PD_SED_distorsion, SurCorrFactor_distorsion, PD_SED_dilatation, SurCorrFactor_dilatation] = Calculate_SurCorrection(Deltas,BCS,BCD,disp,TotalNumMatPoint,numfam,nodefam,pointfam,coord)
        Elastic_Modulus = 200e9; % unit: N/m^2
        % Classical strain energy for dilatation
        SED_analytical_dilatation = 0.001;
        % PD strain energy for distorsion=
        Possion_ratio = 0.3;  
        Shear_Modulus = Elastic_Modulus / (2 * (1 + Possion_ratio)); 
        % bulk modulus % unit: N/m^2
        Bulk_Modulus = Elastic_Modulus / (2 * (1 - Possion_ratio)); 
        % PD material parameter
        alpha=0.5*(Bulk_Modulus-2*Shear_Modulus);
        SED_analytical_distorsion = 1 / (2*(1 - Possion_ratio*Possion_ratio)) * (Elastic_Modulus) * (0.001)^2 - alpha * (0.001)^2; 
        
        PD_SED_distorsion = zeros(TotalNumMatPoint,1);
        SurCorrFactor_distorsion = zeros(TotalNumMatPoint,1);
        PD_SED_dilatation = zeros(TotalNumMatPoint,1);
        SurCorrFactor_dilatation = zeros(TotalNumMatPoint,1);
        thick = Deltas(1,1) / 3.015;
        
       for i = 1:TotalNumMatPoint
        for j = 1:numfam(i,1)
        cnode = nodefam(pointfam(i,1)+j-1,1);
        RelativePosition_Vector = sqrt((coord(cnode,1) - coord(i,1))^2 + (coord(cnode,2) - coord(i,2))^2);
        RelativeDisp_Vector=sqrt((coord(cnode,1)+disp(cnode,1)-coord(i,1)-disp(i,1))^2+(coord(cnode,2)+disp(cnode,2)-coord(i,2)-disp(i,2))^2);
        Stretch = (RelativeDisp_Vector - RelativePosition_Vector) / RelativePosition_Vector;
        AbsoluteValue_x_y = RelativeDisp_Vector * RelativePosition_Vector;
        Coeff_x = (coord(cnode,1) + disp(cnode,1) - coord(i,1) - disp(i,1)) * (coord(cnode,1) - coord(i,1));
        Coeff_y = (coord(cnode,2) + disp(cnode,2) - coord(i,2) - disp(i,2)) * (coord(cnode,2) - coord(i,2));
        Directional_cosine = (Coeff_x + Coeff_y) / AbsoluteValue_x_y;
        VolCorr_radius = Deltas(cnode,1) / 3.015 / 2;
            if (RelativePosition_Vector <= Deltas(i,1)-VolCorr_radius)
            fac = 1;
            elseif (RelativePosition_Vector <= Deltas(i,1)+VolCorr_radius)
            fac = (Deltas(i,1)+VolCorr_radius-RelativePosition_Vector)/(2*VolCorr_radius);
            else
            fac = 0;
            error("error for material point rel and del as follows %f %f", RelativePosition_Vector, Deltas(i,1)-VolCorr_radius);
            end
        Volume = (Deltas(cnode,1) / 3.015) ^ 2 * thick;
        PD_SED_distorsion(i,1) = PD_SED_distorsion(i,1) + BCS(i,1) * Deltas(i,1) * (Stretch^2) * (RelativePosition_Vector) * Volume * fac;
        PD_SED_dilatation(i,1) = PD_SED_dilatation(i,1) + BCD(i,1) * Deltas(i,1) * Stretch * Directional_cosine * Volume * fac;
        end
        SurCorrFactor_distorsion(i,1) = SED_analytical_distorsion / PD_SED_distorsion(i,1);
        SurCorrFactor_dilatation(i,1) = SED_analytical_dilatation / PD_SED_dilatation(i,1);
       end
end

function [sorted] = SortPath(p_name, list, coord)
    
    if (strcmp(p_name,'hor')) %Checks if p_name is equal to 'hor'
        corr_coord = zeros(length(list), 1); % Creat a list of zeros with the length of path(list) size.
        for i = 1: length(list)
            corr_coord(i,1) = coord(list(i), 1); %Derive the coordinates of all the points in the path into corr_coord
        end
    [~,sortIdx] = unique(corr_coord, 'last');
    sorted = list(sortIdx);
    elseif(strcmp(p_name, 'ver') || strcmp(p_name,'edge'))
        corr_coord = zeros(length(list), 1);
        for i = 1: length(list)
            corr_coord(i,1) = coord(list(i), 2);
        end
     [~,sortIdx] = unique(corr_coord);
     sorted = list(sortIdx);
    end
        
end
function [extracted_data] = ExtractPathData(path, data, direction)
    extracted_data = [];
    for i = 1:max(size(path))
        extracted_data(end+1) = data(path(i), direction);
    end
end

function [averaged_data] = AverageDispPath(path, data, direction)
    extracted_data = [];
    for i = 1:max(size(path))
        extracted_data(end+1) = data(path(i), direction);
    end
end

%%%The function below calculates the stresses in the x and y direction of
%%%any point on the plate and returns a stress array.
function [stress] = CalculateStressforPoint(coord,TotalNumMatPoint,numfam,nodefam, pointfam, Deltas)
    stress = zeros(TotalNumMatPoint, 2); %Predifining the size of the stress array to speed up the calculations.

    for i = 1: TotalNumMatPoint
        f_x = 0; %Sum of all forces in X direction
        f_y = 0; %Sum of all forces in Y direction
        satisfy_x = 0;
        satisfy_y = 0;

        localThickness = Deltas(1,1) / 3.015;
        point_x = coord(i,1); %% Extract coordinates of point i
        point_y = coord(i,2);
        for j = 1: TotalNumMatPoint %% Iterating every other material point for conditions
            if (coord(j,2) == point_y && coord(j,1) < point_x + localThickness / 2) %% coord(j,2) == point_y &&
                    satisfy_x = satisfy_x + 1;
                    for k = 1:numfam(j,1) %%Iterate for all the neighbors of point j that satisfy the conditions.
                        cnode = nodefam(pointfam(j,1)+k-1,1);
                        if (coord(cnode,1) > point_x + localThickness / 2) % Sum bondforce of the neighbors that satisfy this condition
                            f_x = f_x + nodefam(pointfam(j,1)+k-1, 2); %%% Requires Saving BondForce_x in nodefam(:,2)
                        end
                    end
            end
            %Same below for stress in y direction
            if (coord(j,1) == point_x && coord(j,2) < point_y + localThickness / 2) %%%%coord(j,1) == point_x &&
                    satisfy_y = satisfy_y + 1;
                    for k = 1:numfam(j,1) %%Iterate for all the neighbors of point j that satisfy the conditions.
                        cnode = nodefam(pointfam(j,1)+k-1,1);
                        if (coord(cnode,2) > point_y + localThickness / 2) % Sum bondforce of the neighbors that satisfy this condition
                            f_y = f_y + nodefam(pointfam(j,1)+k-1, 3); %%% Requires Saving BondForce_y in nodefam(:,3)
                        end
                    end
            end
            stress_x = f_x * l;
             
            stress_y = f_y * l;
        end
        %Saving the stress for point i before moving on to the next point
        %on the plate
        stress(i,1) = stress_x;
        stress(i,2) = stress_y;
    end
end

function [coord, path, numOfRemoteMPs, totalNumMatPoint] = generateBarDense(NumofDiv_x, length, width, dx, dxDense, matPointNum, center_hole)

NumofDiv_y = NumofDiv_x;
path = [];
coord_excess = zeros(matPointNum, 2);
counter = 0; % a counter used in defining coordinates of each material point
local_start_yy = (-1) * width / 4;
%TODO: Fix the dependency of the width of the local area with its density
%in generateBarDense.
warning('TODO: Fix the dependency of the width of the local area with its density in generateBarDense.');
for i = 1:NumofDiv_x
    for j = 1:NumofDiv_y
      coordx = -1/2*length + (dx/2) + (i - 1)*dx;
      coordy = -1/2*width + (dx/2) + (j - 1)*dx;
      
      %Applying the hole in the plate (can be deactivated by commenting the
      %if statement below%
      
      if (coordy > local_start_yy && coordy < (-1) * local_start_yy)
          continue
          %nullpoint(nnum,1) = 0;
      end
      
      %Defining paths of material points similar to the path feature in
      %ABAQUS%%%%%%%(ONLY THE HORIZONTAL PATH WORKS)

      counter = counter + 1;
      coord_excess(counter,1) = coordx; %A coord_excess is defined since initially a larger than needed array size has to be used.
      coord_excess(counter,2) = coordy; %coord-excess is trimmed later as the end elements are empty.
       
    end
end

coord_excess = coord_excess(1:counter, :); %coord_excess is trimmed here once
numOfRemoteMPs = counter;

for i = 1:NumofDiv_x * 2 - 1
    for j = 1:NumofDiv_y
        
      coordx = -1/2*length + dx / 2 +  (i - 1)*dxDense;
      coordy = local_start_yy + (dxDense/2) + (j - 1)*dxDense;
      
     if (center_hole.inEllipse(coordx, coordy))
          continue
          %nullpoint(nnum,1) = 0;
     end
     
     if ( abs(coordy) <= dxDense && coordx >= 0 && coordy >= 0)
          path(end+1) = counter + 1;
     end
      counter = counter + 1;
      coord_excess(counter,1) = coordx; %A coord_excess is defined since initially a larger than needed array size has to be used.
      coord_excess(counter,2) = coordy; %coord-excess is trimmed later as the end elements are empty.
       
    end
end

coord = coord_excess(1:counter, :); %coord: Material point locations
totalNumMatPoint = size(coord, 1);

end

function [coord, path, numOfRemoteMPs, totalNumMatPoint] = generateCircleDense(NumofDiv_x, length, width, dx, dxDense, matPointNum, center_hole, left_tip, right_tip)
NumofDiv_y = NumofDiv_x;
path = [];
coord_excess = zeros(matPointNum, 2);
counter = 0; % a counter used in defining coordinates of each material point

for i = 1:NumofDiv_x
    for j = 1:NumofDiv_y
      coordx = -1/2*length + (dx/2) + (i - 1)*dx;
      coordy = -1/2*width + (dx/2) + (j - 1)*dx;
      
      %Applying the hole in the plate (can be deactivated by commenting the
      %if statement below%
      
      if (center_hole.inEllipse(coordx, coordy) || left_tip.inEllipse(coordx, coordy) || right_tip.inEllipse(coordx, coordy))
          continue
          %nullpoint(nnum,1) = 0;
      end
      
      %Defining paths of material points similar to the path feature in
      %ABAQUS%%%%%%%(ONLY THE HORIZONTAL PATH WORKS)
           
     if ( abs(coordy) <= dx && coordx >= 0 && coordy > 0)
          path(end+1) = counter + 1;
     end

      counter = counter + 1;
      coord_excess(counter,1) = coordx; %A coord_excess is defined since initially a larger than needed array size has to be used.
      coord_excess(counter,2) = coordy; %coord-excess is trimmed later as the end elements are empty.
       
    end
end
coord = coord_excess(1:counter, :); %coord_excess is trimmed here once
numOfRemoteMPs = counter;

seedLeft = get_circle(left_tip.x_center, left_tip.y_center, left_tip.radius_major, dxDense, center_hole);
seedRight = get_circle(right_tip.x_center, right_tip.y_center, right_tip.radius_major, dxDense, center_hole);
coord = [coord; seedLeft; seedRight];
totalNumMatPoint = size(coord, 1);


for i = 1: numOfRemoteMPs
    
    coordx = coord(i, 1);
    coordy = coord(i, 2);
    if ( abs(coordy) <= dx && coordx >= 0 && coordy > 0)
          path(end+1) = i;
    end
end

for i = numOfRemoteMPs + 1: totalNumMatPoint
    coordx = coord(i, 1);
    coordy = coord(i, 2);
    if (abs(coordy) <= dxDense && coordx >= 0 && coordy > 0)
        path(end+1) =i;
    end
end

end

function saveAndPrintFigure(gcf, address, figureName)
    if (isempty(address) == 0)
        saveas(gcf, strcat(address, figureName));
    end
end

