%% Notes
%Note that it is assumed that there are two sheets per unit cell, and in the POSCAR file, the last half of the coordinates correspond to the top sheet.
% Assumed that the unit cell is defined as square or rectangular.
clear; clc; warning off;
    
%% Define Parameters
%For making many structures, fill out the "vect" vector.
    % Column 1 = "name". This can be anything, but it must be unique among the rows
    % Column 2 = Rotation angle in degrees.
    % Column 3 = Translation in x- direction. In units of fraction of lattice vector, a. (e.g., 0.5 is a translation of a/2).
    % Column 4 = Radius of disk. Units of angstrom.
vect = [
100	60	0	25
101	50	0.1	25
102	40	0.25	25
103	20	0.333	25
104	30	0.5	25
105	30	0	20
106	60	0.1	20
107	20	0.25	20
108	50	0.333	20
109	40	0.5	20
110	50	0	15
111	30	0.1	15
112	60	0.25	15
113	40	0.333	15
114	20	0.5	15
115	20	0	10
116	40	0.1	10
117	30	0.25	10
118	60	0.333	10
119	50	0.5	10
120	40	0	5
121	20	0.1	5
122	6	0.25	5
123	30	0.333	5
124	60	0.5	5
];


for mm = 1:size((vect),1)
    %% Extracting values from vect
    Name = num2str(vect(mm,1)); %used for appending output filename
    angle = vect(mm,2); % angle of rotation for top sheet
    Translation = vect(mm,3); % fraction of lattice vector for translation
    R = vect(mm,4); % radius of nanosheet
    
        
    %% Reading POSCAR
    clear dummy;
    FileIn = 'POSCAR_50x50';
    fid = fopen(FileIn);
    count = 1;
    tline = fgetl(fid);
    data{count} = {tline};
    while ischar(tline)
        count = count +1;
        tline = fgetl(fid);
        data{count} = {tline};
    end
    fclose(fid);
    for ii = 8:(length(data)-1)
        dummy(ii-7,:) = data{ii};
    end
    for ii = 1:length(dummy)
        coordinates(ii,:) = str2num(dummy{ii});
    end
    BottomSheet = coordinates(1:(length(coordinates)/2),:);
    TopSheet = coordinates((length(coordinates)/2 + 1):end,:);
    

    %% Applying rotation
    RotationMatrix = @(theta) [cosd(theta) -sind(theta) 0;    sind(theta) cosd(theta) 0;    0 0 1];
    
    %Defining the center of the sheet. Assumed to be half of x- and y- lattice vectors. (works only for square or rectangular unit cells.)
    LatticeVectors = [str2num(data{3}{1}); str2num(data{4}{1}); str2num(data{5}{1})];
    xcenter = LatticeVectors(1,1)./2;
    ycenter = LatticeVectors(2,2)./2;
    
    %Translating the center of the top sheet to the origin
    TopSheet(:,1) = TopSheet(:,1) - xcenter;
    TopSheet(:,2) = TopSheet(:,2) - ycenter;
    
    %Rotating about the origin
    RotatedCoordinates = [RotationMatrix(angle)*TopSheet']';
    
    %Translating the center of the top sheet back to the center of the unit cell
    RotatedCoordinates(:,1) = RotatedCoordinates(:,1) + xcenter;
    RotatedCoordinates(:,2) = RotatedCoordinates(:,2) + ycenter;
    
    %% Applying translation of top sheet
    a = 3.8919225;
    T = Translation*a;
    RotatedCoordinates(:,1) = RotatedCoordinates(:,1) + T;
    
    %% Concatenating position matrix
    AllCoordinates = [BottomSheet; RotatedCoordinates];
    
    %% Getting rid of all atoms at distances greater than "R"
    AtomDistances = sqrt((AllCoordinates(:,1)-xcenter).^2 + (AllCoordinates(:,2)-ycenter).^2);
    AllCoordinates(AtomDistances>R,3) = NaN;
    
    %% Removing Atoms at distances > D
    Atom_Amounts = str2num(data{6}{1});
    dummy = [];
    dummy(1) = 0;
    for ii = 1:length(Atom_Amounts);
        dummy(ii+1) = dummy(ii)+Atom_Amounts(ii);
    end
    for ii = 1:(length(dummy)-1)
        new_Atom_Amounts(ii) = sum(~isnan(AllCoordinates((dummy(ii)+1):dummy(ii+1),3)));
    end
    data{6}{1} = num2str(new_Atom_Amounts);
    AllCoordinates(isnan(AllCoordinates(:,3)),:) = [];
    
    
    %% Translating sheets closer to the origin, with a 2.5 angstrom offset
    AllCoordinates(:,1)  = AllCoordinates(:,1) - min(AllCoordinates(:,1)) + 2.5;
    AllCoordinates(:,2)  = AllCoordinates(:,2) - min(AllCoordinates(:,2)) + 2.5;
    OrigionalLength = length(AllCoordinates);
    
    
    %% Trimmig unit cell to adjust for nearst sheet distance to be > 10 angstrom
    a_vector = str2num(data{3}{1});
    b_vector = str2num(data{4}{1});
    
    a_vector(1) = max(AllCoordinates(:,1)) + 15;
    b_vector(2) = max(AllCoordinates(:,2)) + 15;
    
    data{3}{1} = num2str(a_vector);
    data{4}{1} = num2str(b_vector);
    
    
    %% Adding hydrogens to saturate all Si atoms
    xcenter = (max(AllCoordinates(:,1)) -   min(AllCoordinates(:,1)))/2 + min(AllCoordinates(:,1));
    ycenter = (max(AllCoordinates(:,2)) -  min(AllCoordinates(:,2)))/2 + min(AllCoordinates(:,2));
    
    NewCoord = []; count = 1; AtomCount = [];
    %For each atom, find the number of nearest neighbors.
    for ii = 1:length(AllCoordinates)
        Pos1 = AllCoordinates(ii,:);
        AtomCount(ii) = 0;
        for jj = 1:length(AllCoordinates)
            Pos2 = AllCoordinates(jj,:);
            distance = sqrt((Pos1(1) - Pos2(1))^2 + (Pos1(2) - Pos2(2))^2 + (Pos1(3) - Pos2(3))^2);
            if distance < 2.4 & distance > 0.1
                if abs(Pos1(3) - Pos2(3)) < 2 % Intersheet H-H is > 2, so this avoids them.
                    test(ii) = distance;
                    AtomCount(ii) = AtomCount(ii)+ 1;
                end
            end
        end
        
        % If the iith atom has 2 nearest neighbors, add two hydrogens at positions defined by equations E1, E2, E3, and E4.
        if AtomCount(ii) == 2 % note that V = [xh1 yh1 xh2 yh2]
            xa = AllCoordinates(ii,1);
            ya = AllCoordinates(ii,2);
            dac = sqrt((xa - xcenter).^2 + (ya - ycenter).^2); % distance between iith atom and the center of the sheet
            E1 = @(V) (-1.5 + sqrt((xa - V(1)).^2 + (ya - V(2)).^2)).^2; % Distance between 1st H and Si = 1.5 angstrom
            E2 = @(V) (-1.5 + sqrt((xa - V(3)).^2 + (ya - V(4)).^2)).^2; % Distance between 2nd H and Si = 1.5 angstrom
            E3 = @(V) (-2.6 + sqrt((V(1) - V(3))^2 + (V(2) - V(4))^2)).^2; % Distance between 1st and 2nd H is 2.6 angstroms (120 degrees apart)
            E4 = @(V) (sqrt((V(1) - xcenter)^2 + (V(2) - ycenter)^2) - sqrt((V(3) - xcenter)^2 + (V(4) - ycenter)^2)).^2; % Distances between both hydrogens and the center of the sheet are equal.
            E5 = @(V) (sqrt((xcenter - V(1)).^2 + (ycenter - V(2)).^2) - sqrt(dac.^2 + 1.5^2 - 2*(1.5)*dac*cosd(120))).^2; % Angle of H1-Si-center = 120 degrees.
            
            rff = @(V) [E1(V); E2(V); E3(V); E4(V); E5(V)];
            out = fminsearch(@(V) norm(rff(V)),[xa ya xa ya]);
            NewCoord(count,:) = [out(1) out(2) (AllCoordinates(ii,3))]; count = count+1; % Adding new atom to vector called "NewCoord"
            NewCoord(count,:) = [out(3) out(4) AllCoordinates(ii,3)]; count = count+1; % Adding new atom to vector called "NewCoord"
        end
        
        % If the iith atom has 3 nearest neighbors, add one hydrogen at the position defined by "rff"
        if AtomCount(ii) == 3 % Note that y = mx + b. We are adding hydrogen radially out from the Si, at a distance of 1.5 angstrom.
            %y = mx +  b is defined by the line connecting the center of the sheet to the iith atom.
            xa = AllCoordinates(ii,1);
            ya = AllCoordinates(ii,2);
            m = (ya - ycenter)./(xa - xcenter);
            b = ya - m*xa;
            yH = @(xH) m*xH + b;
            rff = @(xH) 1.5 + sqrt((xcenter - xa).^2 + (ycenter - ya).^2) - sqrt((xcenter - xH).^2 + (ycenter - yH(xH)).^2);
            xH = fzero(rff,xa);
            yH = yH(xH);
            NewCoord(count,:) = [xH yH AllCoordinates(ii,3)]; count = count+1 ;
        end
    end
    
    %% Appending additional atoms
    AllCoordinates = [AllCoordinates; NewCoord];
    data{1}{1} = [data{1}{1} 'H'];
    data{6}{1} = [data{6}{1} ' ' num2str(length(AllCoordinates) - OrigionalLength)];
    
    
    %% Writing File
    new_poscar = {data{1}{1}; data{2}{1}; data{3}{1}; data{4}{1}; data{5}{1}; data{6}{1}; data{7}{1}};
    filePh = fopen(['POSCAR_' Name],'w');
    fprintf(filePh,'%s\n',new_poscar{:});
    fclose(filePh);
    dlmwrite(['POSCAR_' Name],AllCoordinates,'-append','delimiter',' ', 'precision', 16)
    Name
end