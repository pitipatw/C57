
close all
clear
clc

fprintf("LOOP\n")
set_ran = 0;
set_rmin = [sqrt(3) 2 2.5]; %3 
set_cnt_interval = [5 10 20 40]; %4
set_cnt_step = [1 0.5];%2
set_penal0 = [1 2 3]; %3
%total = 3*4*2*3 =  problems

%%test set 
set_rmin = [sqrt(3)]; %3 
set_cnt_interval = [5]; %4
set_cnt_step = [1];%2
set_penal0 = [1]; %3

for rmin = set_rmin
    for cnt_interval = set_cnt_interval 
        for cnt_step = set_cnt_step
            for penal0 = set_penal0 
close all

disp("Checking save directory (Results)...")
if ~exist('Results')
    disp("      Results directory could not be found")
    disp("      Results directory created")
    mkdir('Results')
else
    disp("      Results directory is already exist")
    disp("      Doing nothing")
end
disp("Entering the main script...")


% addpath(genpath('Functions'))
% addpath(genpath('LoadCases'))
addpath(genpath('MMA'))
addpath(genpath('stenglib-master'))

disp("Related functions and paths added")

try
    parpool("Processes")
catch
    disp("Parallel already exist, moving on!")
end

% generating unique random value.
todayDate =  datetime; 
% Pid = datenum(todayDate); This is unreadable by a person :(
Pid = replace(string(todayDate), ":", "-");
filename = ['Results//vol0.3_output-' replace(num2str(Pid),' ','-')];

nelx = 1;
nely = 90;
nelz = 30;


nely = 10;
nelz = 3*nely;


showPlot = true ; 
volfrac = 0.5;
penal = penal0;
% rmin = rmin; %already defined.
ft = 1;
ftBC = 'notN';
eta = 0.5;
beta = 0.001;
divis = 0.1;
move = 0.05; %0.1 or 0.05, do this with MMA instead.
maxit = 400;
mode = "MMA";
% mode = "OC";
problemCase = 1 ; %ContinuumCantilever.
symmetry = "" ; %axes to reflect.
%should make a report on the variables printed into file log
fprintf( ['Run script with\n' ...
    '   nelx, nely, nelz : %i, %i, %i\n' ...
    '   volfrac: %g\n' ...
    '   penal: %d\n' ...
    '   rmin: %g\n' ...
    '   ft: %d, ftBC:%s\n' ...
    '   eta: %g\n' ...
    '   beta: %g\n' ...
    '   divis: %g\n'...
    '   move: %g\n' ...
    '   maxit:%d\n' ...
    '   cnt_interval:%d\n'...
    '   cnt_step:%d\n'...
    '   mode:%s\n'], ...
      nelx,nely,nelz,...
      volfrac, penal ,rmin,...
      ft, ftBC, eta, beta,divis, move,...
      maxit,...
      cnt_interval, cnt_step,...
      mode);

%update the filename to contain input parameters.
filename = strcat(filename ,'_' ,string(nelx), '-', string(nely), '-' ,string(nelz) ,...
    '-', string(volfrac) ,'-' ,string(nelz) ,'-', string(rmin), '-' ,string(divis),...
    '_', string(cnt_interval),'-', string(cnt_step)) ;

tic
[xPhys,nodalCoord, elemCoord,F, fixedDof] = topOrtho3Dloop(nelx,nely,nelz,volfrac,penal,rmin,ft,ftBC,eta,beta,divis,move,maxit,...
    mode,problemCase, showPlot,filename, cnt_interval, cnt_step);
toc

% xPhys
isoplot(xPhys,nelx,nely,nelz, "end")


% use saving scheme from pMatlab documentation
% filename = ['output.' num2str(Pid) '.mat'];
% save(filename);

% filename = ['output.' num2str(Pid) '.mat'];
% save(filename, 'variable1', 'variable2', ..., 'variableN');

% generating unique random value.
% todayDate =  datetime; 
% Pid = datenum(todayDate);
% filename = ['Results//output.' num2str(Pid) '.mat'];
save(strcat(filename , '.mat'));
filename

figname = strcat(replace(filename, "output", "fig"), ".png"); 
saveas(figure(2), figname); 



%update the plot with a checkmark

set_ran = set_ran + 1 ;
fprintf("######\nSet ran: %d\n######\n", set_ran)
            end
        end
    end
end
% state = sum( (Pid + 1)*clock*100. );
% rng(state);

% [X,Y,Z] = meshgrid(1:nelx ,1:nely, 1:nelz)
% scatter3(X,Y,Z, reshape(xPhys, 1, 40,25))

%should make a report on the variables printed into file log



function [xPhys,nodalCoord, elemCoord,F, fixedDof] = topOrtho3Dloop(nelx,nely,nelz,volfrac,penal,rmin,ft,ftBC,eta,beta,divis,move,maxit, mode, ...
    problemCase, showPlot, filename, cnt_interval, cnt_step)

%% PRE. 1) MATERIAL AND CONTINUATION PARAMETERS
fprintf("#######\nOpt:" + mode+"\n#######\n")

Ec = 10000;                                                                % Young modulus of solid
Emin = 1e-9;                                                               % Young modulus of "void"
Et = divis*Ec;
nuC = 0.3;  %default is 0.3                                                % Poisson ratio
nuT = divis*nuC ;

% penalCnt = { 1, 1, 25, 0.25 };                                           % continuation scheme on penal
% betaCnt  = { 1, 1, 25,    2 };

% penalCnt = { 1, 3, 40, 0.5 };                                              % continuation scheme on penal
penalCnt = { 1, 3, cnt_interval, cnt_step };                                              % continuation scheme on penal
betaCnt  = { 1, 1, 200,    2 };   % continuation scheme on beta

if ftBC == 'N', bcF = 'symmetric'; else, bcF = 0; end                      % filter BC selector
%% PRE. 2) DISCRETIZATION FEATURES
%input: nelx, nely, nelz
%output: nEl, nodeNrs, elemNrs, cMat, nDof

% FROM 125.
nEl = nelx*nely*nelz; %number of elements
nDof = ( 1 + nely ) * ( 1 + nelz ) * ( 1 + nelx ) * 3;                     % total number of DOFs        #3D#

nodeNrs = int32( reshape( 1 : ( 1 + nelx ) * ( 1 + nely ) * ( 1 + nelz ), ...
    1 + nely, 1 + nelz, 1 + nelx ) ); %3D matrix representing node number and their positions.
elemNrs = int32( reshape( 1 : ( nelx ) * ( nely ) * ( nelz ), ...
    nely, nelz, nelx ) ); %3D matrix representing node number and their positions.
cVec = reshape( 3 * nodeNrs( 1 : nely, 1 : nelz, 1 : nelx ) + 1, nEl, 1 ); % Vector of the next Dof of the node.   #3D#
cMat = cVec+int32( [0,1,2,3*(nely+1)*(nelz+1)+[0,1,2,-3,-2,-1],-3,-2,-1,3*(nely+...
    1)+[0,1,2],3*(nely+1)*(nelz+2)+[0,1,2,-3,-2,-1],3*(nely+1)+[-3,-2,-1]]);% connectivity matrix         #3D#
%%

lcDof = 3 * nodeNrs( 1 , nelz+1, :)-1 ; 

fixedNodes1 = nodeNrs(nely+1, 1, :); 
fixedNodes2 = nodeNrs(1:nely+1, nelz+1, :); %symmetry part
fixedNodes3 = nodeNrs(:,:,:) ; %restrain out of plane.
% fixedNodes2 = [];


fixedDof1 = [3*fixedNodes1-1];
fixedDof2 = [3*fixedNodes2] ;
fixedDof3 = [3*fixedNodes3-2];

fixedDof12 = union(fixedDof1, fixedDof2);
fixedDof= union(fixedDof12, fixedDof3); 

[ pasS, pasV ] = deal( [], [] );                                           % passive solid and void elements
act = setdiff( ( 1 : nEl )', union( pasS, pasV ) );                        % set of active d.v.

F = fsparse( lcDof(:), 1, -1, [ nDof, 1 ] );             % define load vector
freeDof = setdiff( 1 : nDof, fixedDof )';
%%
% NOT FROM 125.
% [nodeNrs, elemNrs, nodalCoord, cArr, elemNodalCoord, cMat] = mesh3DSolid(nelx + 1, nely + 1, nelz + 1);

%note the shape 1 by nEl of cVec, this is a column vector
%the int32 ... is a row vector
% adding them together permutes the result.
%the int32... is the relative position of the cube dof (24 positions)
%the result follows this orientation
%
%     3----7
% 3  /|   /|   <-- this is a hexahedron.
%   4 2  8 6
%   |/   |/
%   1----5

%iK and jK maps the local coordinate (order by the elements to the global
%coordinate K.
[ sI, sII ] = deal( [ ] ); %sI, sII maps rows and columns of upper tri matrix.
for j = 1 : 24
    sI = cat( 2, sI, j : 24 );
    sII = cat( 2, sII, repmat( j, 1, 24 - j + 1 ) );
end
[ iK , jK ] = deal( cMat( :,  sI )', cMat( :, sII )' );

%% PRE. 3) LOADS, SUPPORTS AND PASSIVE DOMAINS

[nodalCoord, elemCoord] = showProblemDef2(nodeNrs, elemNrs, F, fixedDof(:), false);
% showProblemDef(nodalCoord, F, fixedDof, freeDof)

% *** From top125
% [ pasS, pasV ] = deal( [], [] );                                           % passive solid and void elements
% [ pasS, pasV ] = deal( elemNrs(floor(nely/2), 1:nelz, :), [] );                                           % passive solid and void elements
% act = setdiff( ( 1 : nEl )', union( pasS, pasV ) );                        % set of active d.v.

%% PRE. 4) DEFINE IMPLICIT FUNCTIONS
prj = @(v,eta,beta) (tanh(beta*eta)+tanh(beta*(v(:)-eta)))./...
    (tanh(beta*eta)+tanh(beta*(1-eta)));                                   % projection
deta = @(v,eta,beta) - beta * csch( beta ) .* sech( beta * ( v( : ) - eta ) ).^2 .* ...
    sinh( v( : ) * beta ) .* sinh( ( 1 - v( : ) ) * beta );                % projection eta-derivative
dprj = @(v,eta,beta) beta*(1-tanh(beta*(v-eta)).^2)./(tanh(beta*eta)+tanh(beta*(1-eta)));% proj. x-derivative
cnt = @(v,vCnt,l) v+(l>=vCnt{1}).*(v<vCnt{2}).*(mod(l,vCnt{3})==0).*vCnt{4};
%% PRE. 5) PREPARE FILTER
[dy,dz,dx]=meshgrid(-ceil(rmin)+1:ceil(rmin)-1,...
    -ceil(rmin)+1:ceil(rmin)-1,-ceil(rmin)+1:ceil(rmin)-1 );
h = max( 0, rmin - sqrt( dx.^2 + dy.^2 + dz.^2 ) );                        % conv. kernel                #3D#
Hs = imfilter( ones( nely, nelz, nelx ), h, bcF );                         % matrix of weights (filter)  #3D#
dHs = Hs;
% ------------------------ PRE. 6) ALLOCATE AND INITIALIZE OTHER PARAMETERS
[ x, dsK, dV ] = deal( zeros( nEl, 1 ) );                                  % initialize vectors
dV( act, 1 ) = 1/nEl/volfrac;                                              % derivative of volume
x( act ) = ( volfrac*( nEl - length(pasV) ) - length(pasS) )/length( act );% volume fraction on active set
x( pasS ) = 1;                                                             % set x = 1 on pasS set
[ xPhys, xOld, ch, loop, U ] = deal( x, 1, 1, 0, zeros( nDof, 1 ) );       % old x, x change, it. counter, U

%% precalculate stiffness matrix
[gp, gw] = gaussQD(8);
G1 = Ec/nuC/(Ec/(1+nuC)/(1-2*nuC));
G2 = G1 ;
G3 = G1 ;
D0 = Ec/(1+nuC)/(1-2*nuC)*[1+nuC  nuC   nuC   0  0  0
    nuC  1+nuC nuC   0  0  0
    nuC  nuC   1+nuC 0  0  0
    0    0     0     G1 0  0
    0    0     0     0  G2 0
    0    0     0     0  0  G3];

KEi = zeros(24,24); % This is a dummy, as usually KEi is different for each element i.

for i = 1:size(gp,1)
    r = gp(i, 1);
    s = gp(i, 2);
    t = gp(i, 3);
    B = strainDisp(r, s, t);
    % J = jacMap(r, s, t, xE, yE, zE); %determinant should be irrerevant of xE,yE, and zE.
    KEi = KEi + gw(i)*B'*D0*B/8;
end

KE = repmat(KEi, 1, 1 , nEl); %KE is the same for the whole thing.
D  = repmat(D0, 1, 1, nEl);
    U = FEA3D(KE, Emin, xPhys, penal, iK, jK, F, freeDof, U, nEl, nDof);
    c = F'*U;

    fprintf("Initial c: %d\n", c)


conv_count_hist = [];
conv_hist = [];
old_penal = 0 ; %for tracking change in penal
%% INITIALIZE MMA OPTIMIZER
m = 1  ;  % The number of general constraints
n = nEl; % The number of design variables.
xmin = zeros(n,1);
xmax = ones(n,1);
xold1 = x(:);                       % xval, one iteration ago (provided that iter>1).
xold2 = x(:);                       % xval, two iterations ago (provided that iter>2).
low   = ones(n,1);                  % Column vector with the lower asymptotes from the previous iteration (provided that iter>1).
upp   = ones(n,1);                  % Column vector with the upper asymptotes from the previous iteration (provided that iter>1).
% move  = 0.05;

a0    = 1;                % The constants a_0 in the term a_0*z.
a     = zeros(m,1);       % Column vector with the constants a_i in the terms a_i*z.
c_MMA = 10000*ones(m,1);  % Column vector with the constants c_i in the terms c_i*y_i.
d     = zeros(m,1);       % Column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.

%% Initilize a plot that will hold everything from now.


%plots
%5 is the theta convergence
nexttile(5); title("Theta1");

%1 is the total stress
nexttile(1); title("Total Stress")
%6 7 8 are the stresses.
nexttile(2); title("S1");
nexttile(3); title("S2");
nexttile(4); title("S3");

%result (isoplot) is (11, [2,2])
nexttile(7,[2,2])
title("Result 0")
%number of iterations is 13
nexttile(6)
title("Conv Count History")

nexttile(10)
title("Conv Value")


% ================================================= START OPTIMIZATION LOOP
while (ch > 5e-6 && loop < maxit) || penal < 3
    % update iter. counter
    % cont = cont + 1;

    %% Calculate Strain-Displacement Matrix? Local stiffness?
    %for every element

    theta1_hist = [];
    theta2_hist = [];
    E_hist      = [];
    % conv_hist   = [];


    %this is a loop.
    [KE , D, conv_count_hist]= setKE3D(Et,Ec,nelx,nely,nelz,Emin,xPhys,penal,nEl,...
        iK,jK,F,freeDof,cMat,KE,U,nuC,...
        gp,gw,...
        theta1_hist, theta2_hist,E_hist,D,conv_count_hist,...
        elemCoord, loop, filename);


    nexttile(6)
    plot(1:length(conv_count_hist), conv_count_hist)


    % ----------- RL. 1) COMPUTE PHYSICAL DENSITY FIELD (AND ETA IF PROJECT.)
    xTilde = imfilter( reshape( x, nely, nelz, nelx ), h, bcF ) ./ Hs;       % filtered field              #3D#
    xPhys( act ) = xTilde( act );                                            % reshape to column vector
    if ft > 1                              % compute optimal eta* with Newton
        f = ( mean( prj( xPhys, eta, beta ) ) - volfrac )  * (ft == 3);      % function (volume)
        while abs( f ) > 1e-6           % Newton process for finding opt. eta
            eta = eta - f / mean( deta( xPhys, eta, beta ) );
            f = mean( prj( xPhys, eta, beta ) ) - volfrac;
        end
        dHs = Hs ./ reshape( dprj( xPhys, eta, beta ), nely, nelz, nelx );   % sensitivity modification    #3D#
        xPhys = prj( xPhys, eta, beta );                                     % projected (physical) field
    end

    ch = norm( xPhys - xOld,2 ) ./ nEl; %should be no norm, or norm 2
    conv_hist = [conv_hist ch];
    
    xOld = xPhys;

    %% RL. 2) SETUP AND SOLVE EQUILIBRIUM EQUATIONS

    % sK = ( Emin + xPhys.^penal * ( E0 - Emin ) ); %scale for each Ke0
    %Ke0 must be a matrix that has all of the Ke -> the one that is the
    %output of setKE3D.
    dsK( act ) = -penal * ( 1 - Emin ) * xPhys( act ) .^ ( penal - 1 );

    U = FEA3D(KE, Emin, xPhys, penal, iK, jK, F, freeDof, U, nEl, nDof);
    c = F'*U;

    % ------------------------------------------ RL. 3) COMPUTE SENSITIVITIES
    %have to fix the Ke0 situation for derivatives;
    dc = dsK .* sum( squeeze(pagemtimes(reshape(U(cMat)',1,24,[]), KE))' .* U( cMat ), 2) ;
    %  dc = dsK .* sum( ( U( cMat ) * Ke0 ) .* U( cMat ), 2 );                  % derivative of compliance
    dc  = imfilter( reshape( dc, nely, nelz, nelx ) ./ dHs, h, bcF );         % filter objective sens.      #3D#
    dV0 = imfilter( reshape( dV, nely, nelz, nelx ) ./ dHs, h, bcF );        % filter compliance sens.     #3D#
    % ----------------- RL. 4.1) UPDATE DESIGN VARIABLES AND APPLY CONTINUATION
    if mode == "OC"
        xT = x( act );
        [ xU, xL ] = deal( xT + move, xT - move );                               % current upper and lower bound
        ocP = xT .* sqrt( - dc( act ) ./ dV0( act ) ); %this line causes ComplexN.                      % constant part in resizing rule
        l = [ 0, mean( ocP ) / volfrac ];                                        % initial estimate for LM
        while ( l( 2 ) - l( 1 ) ) / ( l( 2 ) + l( 1 ) ) > 1e-4                   % OC resizing rule
            lmid = 0.5 * ( l( 1 ) + l( 2 ) );
            x( act ) = max( max( min( min( ocP / lmid, xU ), 1 ), xL ), 0 );
            if mean( x ) > volfrac, l( 1 ) = lmid; else, l( 2 ) = lmid; end
        end

        %OC report
        fprintf( 'It.:%5i C:%6.5e V:%7.3f ch.:%0.2e penal:%7.2f beta:%7.1f eta:%7.2f lm:%0.2e \n', ...
            loop, F'*U, mean(xPhys(:)), ch, penal, beta, eta, lmid );
    elseif mode == "MMA"
        % ----------------- RL. 4.2) MMA

        xval  = x(:);                                 % design variables
        f0val = c;                                    % objective function
        df0dx = reshape(dc, nEl,1);                   % obj gradient

        dv_c = dV0; %transfer

        fval_c  = sum(xPhys(:))/(volfrac*nEl) - 1;  % calculate contin vol
        dfdx_c   = dv_c(:)';        % normalize cont vol sens

        fval = [fval_c];
        dfdx = [dfdx_c];

        [xmmap, ~, ~, ~, ~, ~, ~, ~, ~, low,upp] = ...
            mmasub(m, n, loop, xval, xmin, xmax, xold1, xold2, ...
            f0val,df0dx,fval,dfdx,low,upp,a0,a,c_MMA,d,0,move);

        % Update MMA Variables with Normal Weighting Projection
        xnew     = xmmap;
        % change = max(abs(xnew(:)-x(:))); currently using norm1.
        xold2    = xold1(:);
        xold1    = x(:);
        x        = xnew;

        %MMA report
        fprintf( 'It.:%5i C:%6.5e V:%7.3f ch.:%0.2e penal:%7.2f\n', ...
            loop, F'*U, mean(xPhys(:)), ch, penal);
    else
        disp("Invalid Mode")
    end
    %%
    old_penal = penal; 
    [penal,beta] = deal(cnt(penal,penalCnt,loop), cnt(beta,betaCnt,loop));   % apply conitnuation on parameters

    % -------------------------- RL. 5) PRINT CURRENT RESULTS AND PLOT DESIGN
    if showPlot
        isoplot(xPhys,nelx,nely,nelz, loop)
    end
    
    loop = loop + 1;

    nexttile(10)
    hold on
    set(gca, 'YScale', 'log')
    xlabel("Iterations")
    ylabel("conv")
    plot(1:loop, conv_hist)
    if old_penal ~= penal
        xline(loop)
    end
    hold off

end
if showPlot
    isoplot(xPhys,nelx,nely,nelz, loop)

    % figure(6)
    % patch(isosurface(isovals, .01),'FaceColor','b','EdgeColor','none');
    % patch(isocaps(isovals, .01),'FaceColor','r','EdgeColor','none');
    % drawnow; view( [ 145, 25 ] ); axis equal tight off;
    % isoplot(xPhys, nelx, nely,nelz)
end

end


%local functions

function [gp, gw] = gaussQD(nnpe)
switch nnpe
    case 8
        pt = 1/sqrt(3);
        gp = [...
            -pt, -pt, -pt;...
            -pt,  pt, -pt;...
            -pt,  pt,  pt;...
            -pt, -pt,  pt;...
             pt, -pt, -pt;...
             pt,  pt, -pt;...
             pt,  pt,  pt;...
             pt, -pt,  pt;...
            ];
        gw = [1,1,1,1,1,1,1,1];
    otherwise
        error('Gauss Quadrature not yet defined')
end
end




function [nodalCoord, elemCoord] = showProblemDef2(nodeNrs,elemNrs, ...
    fNatural, fixedDof, showNodes)
%SHOWPROBLEMDEF
%Plot the problem definition

%fix dimension, in case it's wrong.
if size(fixedDof,1) == 1 
    fixedDof = fixedDof';
end

nEl = size(elemNrs(:),1) ;
nNode = size(nodeNrs(:),1) ;

nDofs = nNode*3; 
allDofs = 1:nDofs ; 
xDofs = 1:3:nDofs ;
yDofs = 2:3:nDofs ;
zDofs = 3:3:nDofs ;


elemMaxX = size(elemNrs,3);
elemMaxY = size(elemNrs,1);
elemMaxZ = size(elemNrs,2);
X = zeros(1,nEl);
Y = zeros(1,nEl);
Z = zeros(1,nEl);
for i = 1:nEl
    [r,c,v] = ind2sub(size(elemNrs),find(elemNrs == i));
    xIdx = elemMaxX - v + 1  ;
    yIdx = r;
    zIdx = c ;
    X(i) = xIdx;
    Y(i) = yIdx;
    Z(i) = zIdx;
end
elemCoord = [X' Y' Z'];

nodeMaxX = size(nodeNrs,3);
% nodeMaxY = size(nodeNrs,1);
% nodeMaxZ = size(nodeNrs,2);
% subplot(2,3,1)

%plot all of the nodes
% axis([0 max(nodeNrs,1) 0 max(nodalCoord(:,2)) 0 max(nodalCoord(:,3))])
X = zeros(1,nNode);
Y = zeros(1,nNode);
Z = zeros(1,nNode);
for i = 1:nNode
    [r,c,v] = ind2sub(size(nodeNrs),find(nodeNrs == i));
    xIdx = nodeMaxX - v + 1;
    yIdx = r;
    zIdx = c;
    X(i) = xIdx;
    Y(i) = yIdx;
    Z(i) = zIdx;
end

nodalCoord = [X' Y' Z'];

%% PLOTTING TIME!
% figure %new unique figure.
figure(2)
hold on
t = tiledlayout(3,4);
title(t,"Dashboard (3D)")
%element plot
% subplot(1,2,1) % node and element plotted separately.
% hold on
% set ( gca, 'YDir', 'reverse' )
% title("Element plot")

% scatter3(elemCoord(:,1), elemCoord(:,2), elemCoord(:,3))

% ePlotIdx = 1:elemMaxY:(elemMaxY)*(elemMaxZ);
% txt = string(ePlotIdx);

% text(elemCoord(ePlotIdx,1), elemCoord(ePlotIdx,2), elemCoord(ePlotIdx,3),txt)

% xlabel("X")
% ylabel("Y")
% zlabel("Z")
% axis equal
% axis tight
% view([45,45])

%node plot
nexttile(9)
hold on
title("Node plot")
scatter3(nodalCoord(:,1), nodalCoord(:,2), nodalCoord(:,3),'x')

%% Supports 
%xRes
xRes = logical(sum(fixedDof == xDofs,1)') ; 
xPts = nodalCoord(xRes,1);
yPts = nodalCoord(xRes,2); 
zPts = nodalCoord(xRes,3);
scatter3(xPts, yPts, zPts,20, 'o', 'r')

%yRes
yRes = logical(sum(fixedDof == yDofs,1)') ; 
xPts = nodalCoord(yRes,1);
yPts = nodalCoord(yRes,2); 
zPts = nodalCoord(yRes,3);
scatter3(xPts, yPts, zPts,30,'o', 'g') 

%zRes
zRes = logical(sum(fixedDof == zDofs,1)') ; 
xPts = nodalCoord(zRes,1);
yPts = nodalCoord(zRes,2); 
zPts = nodalCoord(zRes,3);
scatter3(xPts, yPts, zPts,40,'o', 'b') 

%% Forces
xF = fNatural(xDofs,1) ; 
yF = fNatural(yDofs,1) ; 
zF = fNatural(zDofs,1) ; 
zerosF = zeros(size(xF));

quiver3(nodalCoord(:,1), nodalCoord(:,2), nodalCoord(:,3), xF,zerosF,zerosF,2,'r')
quiver3(nodalCoord(:,1), nodalCoord(:,2), nodalCoord(:,3), zerosF, -yF,zerosF,2,'g')
quiver3(nodalCoord(:,1), nodalCoord(:,2), nodalCoord(:,3), zerosF ,zerosF, zF,2,'b')




xlabel("X")
ylabel("Y")
zlabel("Z")
% legend("Nodes","x-res","y-res","z-res","x-Force","y-Force","z-Force")
axis tight
axis equal
view([45,45])


if showNodes 
    text(nodalCoord(:,1), nodalCoord(:,2), nodalCoord(:,3),string(nodeNrs(:)))
end
end

function isoplot(xPhys,nelx,nely,nelz,loop)

isovals = shiftdim(reshape(xPhys,nely,nelz,nelx),2);
%isoval now is x, y, z

thin = "nah";
% isovals = reshape(xPhys,nelz,nely,nelx);
if size(isovals,3) == 1 
    %repeat on the second 3rd dim
    isovals(:,:,2) = isovals(:,:,1);
    nelx = nelx+ 1;
    thin = "x" ;
end

if size(isovals,2) == 1 
    %repeat on the second 3rd dim
    isovals(:,2,:) = isovals(:,1,:);
    nely = nely+ 1;
    thin = "y" ;
end


if size(isovals,1) == 1 
    %repeat on the second 3rd dim
    isovals(2,:,:) = isovals(1,:,:);
    nelz = nelz+ 1;
    thin = "z" ; 
end
    
isovals = cat(2, isovals, flip(isovals, 2)); 

isovals = smooth3(isovals,'box',1);
% figure(100)
% cla();
nexttile(7,[2,2])
title("Result " + string(loop))
hold on


patch(isosurface(isovals,0.9),'FaceColor','k','EdgeColor','k','FaceAlpha',1,'visible','on');
patch(isocaps(isovals,0.9),'FaceColor','k','EdgeColor','k','FaceAlpha',1,'visible','on');
patch(isosurface(isovals,0.8),'FaceColor','b','EdgeColor','c','FaceAlpha',0.3,'visible','on');
patch(isocaps(isovals,0.8),'FaceColor','c','EdgeColor','c','FaceAlpha',0.3,'visible','on');
patch(isosurface(isovals,0.5),'FaceColor','r','EdgeColor','none','FaceAlpha',0.2,'visible','on');
patch(isocaps(isovals,0.5),'FaceColor','m','EdgeColor','none','FaceAlpha',0.2,'visible','on');

% set ( gca, 'XDir', 'reverse' )
xlabel("Z")
ylabel("Y") 
zlabel("X")
% patch(isosurface(isovals,0.3),'FaceColor','r','EdgeColor','none','FaceAlpha',0.2,'visible','on');
% patch(isocaps(isovals,0.3),'FaceColor','m','EdgeColor','none','FaceAlpha',0.2,'visible','on');
drawnow; %title('Ortho3D result'); 
axis([0 2*nelz+0.5 0 nely+0.5 0 nelx+0.5]); 
axis equal;

xticks(0:5:2*nelz) 
yticks(0:5:nely) 
zticks(0:5:nelx) 
grid on

% view([-37.5, 30])
% view([0, 0])
if thin == "x"
view([0 90])
elseif thin == "y"
    view([0 0])
elseif thin == "z" 
    view([90 0])
end


end



function [U] = FEA3D(KE,Emin,xPhys,penal,...
    iK,jK,F,freeDof,U,nEl, nDof)

%called in set_KE and inside the while loop in top88

% TRUSS K MATRIX
% computation of the system stiffness matrix and convolution parameters,
% done in separate function for clarity

% make xPhys o.o.p array and penalize
% xPhysPen = reshape( (Emin+ xPhys(:).^eta*(1-Emin) ) ,1,1,nele_c);
%mult KE by penalized xPhys

% SimpKE = pagemtimes( reshape(KE,64,1,nele_c) , xPhysPen); 
xPhysPen = (Emin+ xPhys(:).^penal*(1-Emin));

%upper triangular matrix of a 24 by 24 matrix is (24)(24+1)/2 = 300
% SimpKE = zeros(24,24,nele_c);
SimpKE = zeros(300,1, nEl); 
 
% make sK
for e = 1:nEl %go through every element, assigned the value of density (0 and 1, hopefully).
    tmp = triu(KE(:,:,e)*xPhysPen(e)); 
    tmpt = tmp';
    m = tril(true(size(tmpt)));
    SimpKE(:,:,e) = tmpt(m).' ;
    % SimpKE(:,:,e) = KE(:,:,e)*xPhysPen(e) ;

    %this might not help much, but it is good for now.
    %the better way is to initialize D0 upstream as upper triangular matrix
    %already.
end

sK = reshape( SimpKE ,300*nEl,1);
% sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);

Iar = sort( [ iK( : ), jK( : ) ], 2, 'descend' );
K = fsparse( Iar(:, 1), Iar(:, 2) , sK, [nDof, nDof] ) ;
% K = K_c;

% L = chol( K( free, free ), 'lower' );
% U( free ) = L' \ ( L \ F( free ) );  

% COMBINE TRUSS AND CONTINUUM K MATRICES
% solving the partitioned matrix
L = chol( K( freeDof, freeDof ), 'lower' );
U( freeDof ) = L' \ ( L \ F( freeDof ) );

% figure(4)
% colormap(gray); imagesc(-xPhys); axis equal; axis tight; axis off;pause(1e-6); 
% % 
% clf;
% colormap(gray); axis equal; 
% for ely = 1:nely 
%     for elx = 1:nelx 
%         n1 = (nely+1)*(elx-1)+ely; 
%         n2 = (nely+1)* elx +ely; 
%         Ue = 1*U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1); 
%         ly = ely-1; lx = elx-1; 
%         xx = [Ue(1,1)+lx Ue(3,1)+lx+1 Ue(5,1)+lx+1 Ue(7,1)+lx ]'; 
%         yy = [-Ue(2,1)-ly -Ue(4,1)-ly -Ue(6,1)-ly-1 -Ue(8,1)-ly-1]'; 
%         patch(xx,yy,-xPhys(ely,elx))
%     end 
% end
% drawnow; 


end

function   [KE,D, conv_count_hist] = ...
    setKE3D(Et,Ec,nelx,nely,nelz,Emin,xPhys,penal,nEl,...
    iK,jK,F,freedofs,cMat,KE,U,nuC,...
    gp, gw,...
    theta1_hist, theta2_hist,E_hist,D,conv_count_hist,...
    elemCoord, loop, filename)
% Assembling the stiffness matrix from the orthotopic material properties.

% Set parameters
it_max = 1000;
iter = 0;
it_stop = 10;

conv_hist = [];

% if loop == 1 %only the first time.
%     divis_local = max(0.3, Et/Ec);
% else
%     divis_local = Et/Ec;
% end

divis_local = max(0.3 ,Et/Ec);

diff  = 0.05; %step size for the continuation scheme for divis.

Et0 = Ec*divis_local;
nuT = nuC*divis_local; % nu concrete tensile is "divis" times that of nu concrete compression

nDof = 3*(nelx+1)*(nely+1)*(nelz+1) ;
% [gp, gw] = gaussQD(8);

% thresh   = 0.001;
thresh   = 0.1 ;
% Bcell    = {B1; B2; B3; B4};

%pre-allocation output
S1all    = zeros(nelx*nely*nelz,8);
S2all    = zeros(nelx*nely*nelz,8);
S3all    = zeros(nelx*nely*nelz,8);

all_tht  = zeros(8,2,nelx*nely*nelz);
all_l    = zeros(8,3,nelx*nely*nelz);
all_m    = zeros(8,3,nelx*nely*nelz);
all_n    = zeros(8,3,nelx*nely*nelz);

tht1_avg = zeros(nelx*nely*nelz,1);
tht2_avg = zeros(nelx*nely*nelz,1);
l_avg    = zeros(3,nelx*nely*nelz); % nelem by (l1,l2,l3).
m_avg    = zeros(3,nelx*nely*nelz); % nelem by (m1,m2,m3).
n_avg    = zeros(3,nelx*nely*nelz); % nelem by (n1,n2,n3).

% initial FEA
[U] = FEA3D(KE,Emin,xPhys,penal,iK,jK,F,freedofs,U,nEl, nDof);
c = F'*U;

while iter < it_max

    iter        = iter+1;   % max iterations of inner loop
    c0          = c;        %put c into c0 (aka, c_old)
    tht1_avg0    = tht1_avg;  % previous iteration's avg theta, 1x1xnele
    tht2_avg0    = tht2_avg;  % previous iteration's avg theta, 1x1xnele
    fprintf("HelloParfor")
    parfor e = 1:nEl

        KE_matrix = zeros(24,24,8);  % 24x24x (8 gaussian points)
        D_matrix  = zeros(6,6,8);    % 6x6x8 (8 gaussian points)
        displacement = U(cMat(e,:)); % 24x1 displacement of the current element.

        for g = 1:8

            r = gp(g, 1);
            s = gp(g, 2);
            t = gp(g, 3);
            B = strainDisp(r,s,t); % 6x24
            strain = B*displacement; % 6x1 (6x24 x 24x1)

            %stressVoigt = [xx yy zz yz zx xy]
            %this is 6x6 x 6x1 = 6x1
            stressVoigt = D(:,:,e)*strain;  %Please check Voight notaion.

            % loop over gauss points ->  8 points for 3D.
            [D0, stressVoigtPr, coordSys] = matPropCal(stressVoigt, Ec, Et0, nuC, nuT); %change this later

            % D0 = matIsotropic3D(Ec, nuC);
            % stressVoigtPr = diag(stressVoigt) ;
            % coordSys = ones(3,3);

            S1all(e,g) = stressVoigtPr(1,1);
            S2all(e,g) = stressVoigtPr(2,2);
            S3all(e,g) = stressVoigtPr(3,3);

            l1 = coordSys(1, 1);
            m1 = coordSys(2, 1);
            n1 = coordSys(3, 1);
            l2 = coordSys(1, 2);
            m2 = coordSys(2, 2);
            n2 = coordSys(3, 2);
            l3 = coordSys(1, 3);
            m3 = coordSys(2, 3);
            n3 = coordSys(3, 3);

            [l1,m1,n1] = uniqueCosines(l1,m1,n1) ;
            [l2,m2,n2] = uniqueCosines(l2,m2,n2) ;
            [l3,m3,n3] = uniqueCosines(l3,m3,n3) ;

            theta1 = atan(l1/m1); %theta
            theta2 = asin(n1); %phi

            %save them.
            all_tht(g,:,e) = [theta1 theta2];

            all_l(g,:,e) = [l1 l2 l3];
            all_m(g,:,e) = [m1 m2 m3];
            all_n(g,:,e) = [n1 n2 n3];

            % for a more general case : J = jacMap(r, s, t,xE, yE, zE) ;
            D_matrix(:,:,g)= D0;
            KE_matrix(:,:,g) = gw(g)*B'*D0*B ;
        end

        KE_eighth = KE_matrix/8;
        KE(:,:,e) = sum(KE_eighth,3); %sum along the 3 (third) axis.
        D_eighth  = D_matrix/8;
        D(:,:,e)  = sum(D_eighth,3); %sum along the 3 (third) axis.

    end
    
    fprintf("hello FEA")
    [U] = FEA3D(KE,Emin,xPhys,penal,iK,jK,F,freedofs,U,nEl, nDof);
    c = F'*U;

    %check for convergence between old and new c.
    conv = abs((c-c0)/c0);

    %this part should not be divided by 8? they are not from the
    %integration.
    tht1_avg = sum(all_tht(:,1,:))/8;
    tht2_avg = sum(all_tht(:,2,:))/8;

    l_avg =  squeeze(sum(all_l, 1))/8;
    m_avg =  squeeze(sum(all_m, 1))/8;
    n_avg =  squeeze(sum(all_n, 1))/8;

    S1_avg = sum(S1all,2)/8;
    S2_avg = sum(S2all,2)/8;
    S3_avg = sum(S3all,2)/8;

    delta1 = mean( abs(tht1_avg(:) - tht1_avg0(:)).*xPhys(:) );
    delta2 = mean( abs(tht2_avg(:) - tht2_avg0(:)).*xPhys(:) );

    %keep track of variables.
    theta1_hist = [theta1_hist delta1];
    theta2_hist = [theta2_hist delta2];
    E_hist     = [E_hist Et0/Ec];
    conv_hist  = [conv_hist 10*conv];

    %visualize the tracked variables.
    % plotTheta3D(tht1_avg, tht2_avg,S1_avg,S2_avg, S3_avg,...
    %     nelx,nely,nelz,...
    %     theta1_hist,theta2_hist,...
    %     E_hist,thresh,E_ct,E_cc,divis,conv_hist)

    plotTheta3Dnew(l_avg, m_avg, n_avg, elemCoord,...
        S1_avg, S2_avg, S3_avg,...
        theta1_hist,theta2_hist,...
        E_hist,thresh,Et,Ec,divis_local, conv_hist,...
        xPhys, loop, filename);

    thetaSeparated(l_avg, m_avg, n_avg, elemCoord,...
        S1_avg, S2_avg, S3_avg,...
        xPhys, loop, filename);
    % theta1_hist,theta2_hist,...
    % E_hist,thresh,Et,Ec,divis_local, conv_hist,...

    %check convergence
    if conv < thresh
        if round(Et0,3) == round(Et,3)
            break
        end
        divis_local = divis_local-diff;
        Et0 = max(Ec*divis_local,Et);
        divis_local = Et0/Ec;
        nuT = nuC*Et0/Ec;
        % iter = 0;
    end

    if iter == it_max
        disp("exit by itMax")
        break
    end

    % emergency exit if stiffness and compliance aren't converging
    % if iter == it_max-it_stop
    %     if round(E_ct0) == round(E_ct)
    %         break
    %     end
    %     divis = divis-diff;
    %     E_ct0 = max(E_cc*divis,E_ct);
    %     nu_ct = nu_cc*E_ct0/E_cc;
    %     iter = 0;
    % end

end

conv_count_hist = [conv_count_hist iter];

end


%local functions

function [D0, stressVoigtPr, coordSys] = matPropCal(stressVoigt, Ec, Et, nuC, nuT)
inputSize = size(stressVoigt, 1);

% Calculate principal stresses (in Voigt) and transformation matrix
[stressVoigtPr, Q, coordSys] = transMat(stressVoigt);
% [stressVoigtPr, Q, coordSys] = findPrStress(stressVoigt);

% stressVoigtPr
% Gconst = Ec/nuC; %this ratio is constant.
tmpG        = ((1+nuC)/Ec)+ ((1+nuT)/Et);
GconstMixed      = 1/tmpG;
GconstT = Et/2/(1+nuT);
GconstC = Ec/2/(1+nuC);
switch inputSize
    case 3
        if (stressVoigtPr(1) <= 0) && (stressVoigtPr(2) <= 0)
            D = matIsotropic2D(Ec, nuC);
        elseif (stressVoigtPr(1) > 0) && (stressVoigtPr(2) > 0)
            D = matIsotropic2D(Et, nuT);
        elseif size(stressVoigtPr(stressVoigtPr(1:2) <= 0), 1) == 1 %one axis in compression
            D = matOrtho2D(Ec, Et, nuC, nuT);
            els
            error('no material property defined for this stress state.')
        end
    case 6
        if (stressVoigtPr(1,1) <= 0) && (stressVoigtPr(2,2) <= 0) && (stressVoigtPr(3,3) <= 0)
            %all compression
            G = GconstC/(Ec/(1 + nuC)/(1 - 2*nuC));
            D = Ec/(1 + nuC)/(1 - 2*nuC)* ...
                [1+nuC,   nuC,   nuC, 0,0,0; ...
                nuC, 1+nuC,   nuC, 0,0,0; ...
                nuC,   nuC, 1+nuC, 0,0,0; ...
                0,0,0,G,0,0; ...
                0,0,0,0,G,0; ...
                0,0,0,0,0,G; ...
                ];

        elseif (stressVoigtPr(1,1) > 0) && (stressVoigtPr(2,2) > 0) && (stressVoigtPr(3,3) > 0)
            %all tension
            G = GconstT/(Et/(1 + nuT)/(1 - 2*nuT));
            D = Et/(1 + nuT)/(1 - 2*nuT)* ...
                [1+nuT,   nuT,   nuT, 0,0,0; ...
                nuT, 1+nuT,   nuT, 0,0,0; ...
                nuT,   nuT, 1+nuT, 0,0,0; ...
                0,0,0,G,0,0; ...
                0,0,0,0,G,0; ...
                0,0,0,0,0,G; ...
                ];

        elseif (stressVoigtPr(1,1) <= 0) && (stressVoigtPr(2,2) > 0) && (stressVoigtPr(3,3) > 0) %one axis in compression
            %1 compression

            G12 = GconstMixed*(1 - nuT - 2*nuC*nuT) ;
            G13 = GconstMixed*(1 - nuT - 2*nuC*nuT) ;
            G23 = GconstT*(1 - nuT - 2*nuC*nuT) ;
            %could use symmetry on this
            D = 1/(1 - nuT - 2*nuC*nuT)*...
                [Ec*(1-nuT), Ec*nuT                     , Ec*nuT                    ,0,0,0; ...
                Et*nuC    , Et*(1-nuT*nuC)/(1+nuT)     , Et*nuT*(1 + nuC)/(1 + nuT) ,0,0,0; ...
                Et*nuC    , Et*nuT*(1 + nuC)/(1 + nuT) , Et*(1-nuT*nuC)/(1+nuT)     ,0,0,0; ...
                0,0,0,G23,0,0; ...
                0,0,0,0,G13,0; ...
                0,0,0,0,0,G12; ...
                ];

        elseif (stressVoigtPr(1,1) <= 0) && (stressVoigtPr(2,2) <= 0) && (stressVoigtPr(3,3) > 0)%two axis in compression
            % 2 compression
            G12 = GconstC*(1 - nuC - 2*nuC*nuT) ;
            G13 = GconstMixed*(1 - nuC - 2*nuC*nuT) ;
            G23 = GconstMixed*(1 - nuC - 2*nuC*nuT) ;

            %could use symmetry on this
            D = 1/(1 - nuC - 2*nuC*nuT)*...
                [Ec*(1 - nuC*nuT)/(1 + nuC), Ec*nuC*(1 + nuT)/(1 + nuC), Ec*nuT      , 0,0,0; ...
                Ec*nuC*(1 + nuT)/(1 + nuC) , Ec*(1-nuC*nuT)/(1+nuC)    , Ec*nuT      , 0,0,0; ...
                Et*nuC                     , Et*nuC                    , Et*(1 - nuC), 0,0,0; ...
                0,0,0,G23,0,0; ...
                0,0,0,0,G13,0; ...
                0,0,0,0,0,G12; ...
                ];
        else
            stressVoigtPr(1,1)
            stressVoigtPr(2,2)
            stressVoigtPr(3,3)
            error('no material property defined for this stress state.')
        end
        %             This is ideal case, where we can generalize constitutive
        %             matrix in 3D.
        %             D = matOrtho3D(E1, E2, E3, nu1, nu2, nu3);
end

D0=Q*D*Q';
end


function [stressVoigtPr, Q, coordSys] = transMat(stressVoigt)

inputSize = size(stressVoigt, 1);

switch inputSize
    case 3
        stressVoigt = [stressVoigt(1), stressVoigt(2), 0,...
            0, 0, stressVoigt(3)];
    case 6
        % do nothing
    otherwise
        error('Engineering stress input wrong size')
end

[stressVoigtPr, coordSys] = voigtPrincipal(stressVoigt);

switch inputSize
    case 3
        l1 = coordSys(1, 1);
        m1 = coordSys(2, 1);
        l2 = coordSys(1, 2);
        m2 = coordSys(2, 2);
        Q = [...
            l1^2,  l2^2,     2*l1*l2;...
            m1^2,  m2^2,     2*m1*m2;...
            l1*m1, l2*m2, l1*m2+l2*m1];
    case 6

        l1 = coordSys(1, 1);
        m1 = coordSys(2, 1);
        n1 = coordSys(3, 1);
        l2 = coordSys(1, 2);
        m2 = coordSys(2, 2);
        n2 = coordSys(3, 2);
        l3 = coordSys(1, 3);
        m3 = coordSys(2, 3);
        n3 = coordSys(3, 3);

        %fix direction, ensures positive quadrant orientation.
        [l1,m1,n1] = uniqueCosines(l1,m1,n1) ;
        [l2,m2,n2] = uniqueCosines(l2,m2,n2) ;
        [l3,m3,n3] = uniqueCosines(l3,m3,n3) ;

        coordSys2 = [l1 m1 n1; l2 m2 n2; l3 m3 n3] ;

        test_coordSys(coordSys2);

        Q = [...
            l1^2,  l2^2,  l3^2,     2*l2*l3,     2*l1*l3,     2*l1*l2;...
            m1^2,  m2^2,  m3^2,     2*m2*m3,     2*m1*m3,     2*m1*m2;...
            n1^2,  n2^2,  n3^2,     2*n2*n3,     2*n1*n3,     2*n1*n2;...
            m1*n1, m2*n2, m3*n3, m2*n3+m3*n2, m1*n3+m3*n1, m1*n2+m2*n1;...
            l1*n1, l2*n2, l3*n3, l2*n3+l3*n2, l1*n3+l3*n1, l1*n2+l2*n1;...
            l1*m1, l2*m2, l3*m3, l2*m3+l3*m2, l1*m3+l3*m1, l1*m2+l2*m1];

end
end

function test_coordSys(coordSys)



l1 = coordSys(1, 1);
m1 = coordSys(2, 1);
n1 = coordSys(3, 1);
l2 = coordSys(1, 2);
m2 = coordSys(2, 2);
n2 = coordSys(3, 2);
l3 = coordSys(1, 3);
m3 = coordSys(2, 3);
n3 = coordSys(3, 3);

u1 = l1^2 + m1^2 + n1^2;
u2 = l2^2 + m2^2 + n2^2;
u3 = l3^2 + m3^2 + n3^2;


if abs(u1-1) > 1e-3
    disp("1 fails")
end
if abs(u2-1) > 1e-3
    disp("2 fails")
end
if abs(u3-1) > 1e-3
    disp("3 fails")
end

%check orthogonality
v1 = [l1 m1 n1];
v2 = [l2 m2 n2];
v3 = [l3 m3 n3] ;

if dot(v1,v2) >1e-3
    disp("1 2 fails")
end

if dot(v1,v3) >1e-3
    disp("1 3 fails")
end

if dot(v2,v3) >1e-3
    disp("2 3 fails")
end



end


function [stressVoigtPr, coordSys] = voigtPrincipal(stressVoigt)
    % Convert Voigt stresses to tensor stresses
    stressTensor = [...
        stressVoigt(1), stressVoigt(6), stressVoigt(5);...
        stressVoigt(6), stressVoigt(2), stressVoigt(4);...
        stressVoigt(5), stressVoigt(4), stressVoigt(3)];
    
    [coordSys, stressVoigtPr] = eig(stressTensor);
    
    % coordSys(:, 2) = dot(coordSys(:, 2), cross(coordSys(:, 3), coordSys(:, 1))) ...
    %     * coordSys(:, 2); % Ensure output is right-handed system
    % coordSys(:, 3) = dot(coordSys(:, 3), cross(coordSys(:, 1), coordSys(:, 2))) ...
    %      * coordSys(:, 3); % Ensure output is right-handed system
end


function [newl, newm, newn] = uniqueCosines(l,m,n)

if exist('n', 'var')

    if n < 0
        newl = -l ; newm = -m ; newn = -n;
    elseif n < 1e-3
        if m < 0
            newl = -l ; newm = -m ; newn = -n;
        elseif m < 1e-3
            if l < 0
                newl = -l ; newm = -m ; newn = -n;
            else
                newl = l ; newm = m ; newn = n;
            end
        else
            newl = l ; newm = m ; newn = n;
        end
    else
        newl = l ; newm = m ; newn = n;
    end
else

    if m < 0
        newl = -l ; newm = -m ;
    elseif m < 1e-3
        if l < 0
            newl = -l ; newm = -m ;
        else
            newl = l ; newm = m ;
        end
    else
        newl = l ; newm = m ;
    end
end

end

function plotTheta3Dnew(l_avg,m_avg,n_avg,elemCoord, S1,S2,S3,...
                   theta1_hist,theta2_hist,...
                   E_hist,threshold,E_ct,E_cc,divis,conv_hist,...
                   xPhys, loop, filename)
%function plotTheta3Dnew(l_avg,m_avg,n_avg,elemNodalCoord, S1,S2,S3, xPhys)
%This function finds the highest absolute value of stresss at each
%element and shows them as vectors.

%%
% f2 = figure(2);
% clf
% tiledlayout(2)
% f2.Position = p2;

% subplot(2,1,1)
nexttile(5)
plot(1,1)
hold on
% clf
%fprintf(string(divis))
title(['Theta1 at ', num2str(divis)])
xaxis = 1:length(theta1_hist);
plot(xaxis,theta1_hist,'Color',"#A2142F",'LineWidth',2)
plot(xaxis,E_hist    ,'Color',"#EDB120",'LineWidth',2)
plot(xaxis,conv_hist ,'Color',"#7E2F8E",'LineWidth',2)
yline(threshold*10    ,'Color',"#000000",'LineWidth',1.5,'LineStyle',':')
yline(E_ct/E_cc      ,'Color',"#D95319",'LineWidth',1.0,'LineStyle','--')
axis([1 inf  0 5])
% legend('Theta1','Trusses /10','Concrete','Unorm x100')
% legend('Theta1','E_hist','Convergence','Threshold')
hold off

pause(0.001)

%%element location (must be the same as the actual element, the order of
%l, m, and n are the same as the element anyway, therefore, plotting them
%with the elemNodalCoord should be ok. They also have the same order as
%xPhys, therefore, this should work too.
X = elemCoord(:,1)';
Y = elemCoord(:,2)';
Z = elemCoord(:,3)';

%A vector indicating maximum `size` of stress on each element.
D = abs([S1 S2 S3]); %concatinating the stresses from 3 directions
maxCheck = max( D, [], 2); %find the maximum of each row.

%check if it's the maximum among 3 axis 
%Maximum Check
mchk1 = (maxCheck == abs(S1))';
mchk2 = (maxCheck == abs(S2))';
mchk3 = (maxCheck == abs(S3))';

%check if it's 0 or not. we won't plot tick if it's 0
%Zero Check
zchk1 = (abs(S1) < 1e-3)';
zchk2 = (abs(S2) < 1e-3)';
zchk3 = (abs(S3) < 1e-3)';

%check if it's tension (positive) or not (compression, negative) 
% Tension Check
tchk1 = (S1 > 0)';
tchk2 = (S2 > 0)';
tchk3 = (S3 > 0)';

%if they are tensile -> blue, otherwise (compression) red.

%for visual, only show xPhys that are more than 0.5 (the same as in
%isoplot.m
% xchk = (xPhys >= 0.9)' ;
xchk = (xPhys >= 0.0)' ;

%final chk parameters
chk1 = mchk1 & ~zchk1 ;
chk2 = mchk2 & ~zchk2 ;
chk3 = mchk3 & ~zchk3 ;

chk1pos = chk1 & tchk1 ;
chk1neg = ~chk1pos ;

chk2pos = chk2 & tchk2 ;
chk2neg = ~chk2pos ;

chk3pos = chk3 & tchk3 ;
chk3neg = ~chk3pos ;

nexttile(1)
plot(1,1)
hold on
title("Total Stress Plot")
%plot in this order: 
% first Principal Direction (S1 -> S2 -> S3) 
% highest value (LineWidth 2 -> .1) 
% color ('b' -> 'r') 

%% S1
%higher
quiver3(X,Y,Z,  xchk.*chk1pos.*l_avg(1,:), xchk.*chk1pos.*m_avg(1,:), xchk.*chk1pos.*n_avg(1,:),'off', 'Color', 'b', 'LineWidth',2); %Highest value, non-zero, positive
quiver3(X,Y,Z,  xchk.*chk1neg.*l_avg(1,:), xchk.*chk1neg.*m_avg(1,:), xchk.*chk1neg.*n_avg(1,:),'off', 'Color', 'r', 'LineWidth',2); %Highest value, non-zero, negative

%lower
quiver3(X,Y,Z,  xchk.*~chk1pos.*l_avg(1,:), xchk.*~chk1pos.*m_avg(1,:), xchk.*~chk1pos.*n_avg(1,:),'off', 'Color', 'r', 'LineWidth',.01); 
quiver3(X,Y,Z,  xchk.*~chk1neg.*l_avg(1,:), xchk.*~chk1neg.*m_avg(1,:), xchk.*~chk1neg.*n_avg(1,:),'off', 'Color', 'b', 'LineWidth',.01); 

%% S2
quiver3(X,Y,Z,  xchk.*chk2pos.*l_avg(2,:), xchk.*chk2pos.*m_avg(2,:), xchk.*chk2pos.*n_avg(2,:),'off', 'Color', 'b', 'LineWidth',2); %Highest value, non-zero, positive
quiver3(X,Y,Z,  xchk.*chk2neg.*l_avg(2,:), xchk.*chk2neg.*m_avg(2,:), xchk.*chk2neg.*n_avg(2,:),'off', 'Color', 'r', 'LineWidth',2); %Highest value, non-zero, negative

quiver3(X,Y,Z,  xchk.*~chk2pos.*l_avg(2,:), xchk.*~chk2pos.*m_avg(2,:), xchk.*~chk2pos.*n_avg(2,:),'off', 'Color', 'r', 'LineWidth',.01); 
quiver3(X,Y,Z,  xchk.*~chk2neg.*l_avg(2,:), xchk.*~chk2neg.*m_avg(2,:), xchk.*~chk2neg.*n_avg(2,:),'off', 'Color', 'b', 'LineWidth',.01); 

%% S3
quiver3(X,Y,Z,  xchk.*chk3pos.*l_avg(3,:), xchk.*chk3pos.*m_avg(3,:), xchk.*chk3pos.*n_avg(3,:),'off', 'Color', 'b', 'LineWidth',2); %Highest value, non-zero, positive
quiver3(X,Y,Z,  xchk.*chk3neg.*l_avg(3,:), xchk.*chk3neg.*m_avg(3,:), xchk.*chk3neg.*n_avg(3,:),'off', 'Color', 'r', 'LineWidth',2); %Highest value, non-zero, negative

quiver3(X,Y,Z,  xchk.*~chk3pos.*l_avg(3,:), xchk.*~chk3pos.*m_avg(3,:), xchk.*~chk3pos.*n_avg(3,:),'off', 'Color', 'r', 'LineWidth',.01); 
quiver3(X,Y,Z,  xchk.*~chk3neg.*l_avg(3,:), xchk.*~chk3neg.*m_avg(3,:), xchk.*~chk3neg.*n_avg(3,:),'off', 'Color', 'b', 'LineWidth',.01);

axis equal; %colorbar
% view([-37.5, 30])
axis([0 max(X) 0 max(Y) 0  max(Z)])
view([90,0])
% view([45,45])
hold off


%this is for debugging the xPhys location.
% figure(99)
% scatter3(X,Y,Z, 100*(xPhys.^4));

pause(0.001)
end



function [outStat] = thetaSeparated(l_avg,m_avg,n_avg,elemCoord, S1,S2,S3,...
                   xPhys, loop, filename)
                   % theta1_hist,theta2_hist,...
                   % E_hist,threshold,E_ct,E_cc,divis, conv_hist,...
%THETASEPARATED
%(l_avg,m_avg,n_avg,elemNodalCoord, S1,S2,S3, xPhys)
%This function finds the highest absolute value of stresss at each
%element and shows them as vectors,
%separated by the 3 perpendicular angles.

%%element location (must be the same as the actual element, the order of
%l, m, and n are the same as the element anyway, therefore, plotting them
%with the elemNodalCoord should be ok. They also have the same order as
%xPhys, therefore, this should work too.
l_avg(1,:) = (abs(S1).*l_avg(1,:)')';
l_avg(2,:) = (abs(S2).*l_avg(2,:)')';
l_avg(3,:) = (abs(S3).*l_avg(3,:)')';

m_avg(1,:) = -(abs(S1).*m_avg(1,:)')';
m_avg(2,:) = -(abs(S2).*m_avg(2,:)')';
m_avg(3,:) = -(abs(S3).*m_avg(3,:)')';

n_avg(1,:) = (abs(S1).*n_avg(1,:)')';
n_avg(2,:) = (abs(S2).*n_avg(2,:)')';
n_avg(3,:) = (abs(S3).*n_avg(3,:)')';

X = elemCoord(:,1)';
Y = elemCoord(:,2)';
Z = elemCoord(:,3)';
% scale down the l,m,n by S1, S2, S3
% l_avg = (abs(S1).*l_avg')';
% m_avg = (abs(S2).*m_avg')';
% n_avg = (abs(S3).*n_avg')';

%A vector indicating maximum `size` of stress on each element.
D = abs([S1 S2 S3]); %concatinating the stresses from 3 directions
maxCheck = max( D, [], 2); %find the maximum of each row.

%check if it's the maximum among 3 axis 
%Maximum Check
mchk1 = (maxCheck == abs(S1))';
mchk2 = (maxCheck == abs(S2))';
mchk3 = (maxCheck == abs(S3))';

%check if it's 0 or not. we won't plot thick if it's 0
%Zero Check
zchk1 = (abs(S1) < 1e-1)';
zchk2 = (abs(S2) < 1e-1)';
zchk3 = (abs(S3) < 1e-1)';

%check if it's tension (positive) or not (compression, negative) 
% Tension Check
tchk1 = (S1 > 0)';
tchk2 = (S2 > 0)';
tchk3 = (S3 > 0)';
%if they are tensile -> blue, otherwise (compression) red.
%for visual, only show xPhys that are more than 0.5 (the same as in
%isoplot.m
xchk = (xPhys >= 0.8)' ;

%final chk parameters
chk1 = mchk1 & ~zchk1 ;
chk2 = mchk2 & ~zchk2 ;
chk3 = mchk3 & ~zchk3 ;

%%
chk1pos = chk1 & tchk1 ;
chk1neg = ~chk1pos ;
chk2pos = chk2 & tchk2 ;
chk2neg = ~chk2pos ;
chk3pos = chk3 & tchk3 ;
chk3neg = ~chk3pos ;
% figure(31)
% clf

%plot in this order: 
% first Principal Direction (S1 -> S2 -> S3) 
% highest value (LineWidth 2 -> .1) 
% color ('b' -> 'r') 

% figure(31)
% clf

%% S1
nexttile(2)
plot(1,1)
hold on
% title("S1")
%higher
quiver3(X,Y,Z,  xchk.*chk1pos.*l_avg(1,:), xchk.*chk1pos.*m_avg(1,:), xchk.*chk1pos.*n_avg(1,:),0.5,'Color', 'b', 'LineWidth',2); %Highest value, non-zero, positive
quiver3(X,Y,Z,  xchk.*chk1neg.*l_avg(1,:), xchk.*chk1neg.*m_avg(1,:), xchk.*chk1neg.*n_avg(1,:),0.5,'Color', 'r', 'LineWidth',2); %Highest value, non-zero, negative

%lower
quiver3(X,Y,Z,  xchk.*~chk1pos.*l_avg(1,:), xchk.*~chk1pos.*m_avg(1,:), xchk.*~chk1pos.*n_avg(1,:),0.5,'Color', 'r', 'LineWidth',.1); 
quiver3(X,Y,Z,  xchk.*~chk1neg.*l_avg(1,:), xchk.*~chk1neg.*m_avg(1,:), xchk.*~chk1neg.*n_avg(1,:),0.5,'Color', 'b', 'LineWidth',.1); 

axis equal; %colorbar
axis([0 max(X) 0 max(Y) 0 max(Z)])
view([90, 0])
hold off

%% S2
nexttile(3)
plot(1,1)
hold on
% title("S2")
quiver3(X,Y,Z,  xchk.*chk2pos.*l_avg(2,:), xchk.*chk2pos.*m_avg(2,:), xchk.*chk2pos.*n_avg(2,:), 0.5,'Color', 'b', 'LineWidth',2); %Highest value, non-zero, positive
quiver3(X,Y,Z,  xchk.*chk2neg.*l_avg(2,:), xchk.*chk2neg.*m_avg(2,:), xchk.*chk2neg.*n_avg(2,:), 0.5,'Color', 'r', 'LineWidth',2); %Highest value, non-zero, negative

quiver3(X,Y,Z,  xchk.*~chk2pos.*l_avg(2,:), xchk.*~chk2pos.*m_avg(2,:), xchk.*~chk2pos.*n_avg(2,:),0.5,'Color', 'r', 'LineWidth',.1); 
quiver3(X,Y,Z,  xchk.*~chk2neg.*l_avg(2,:), xchk.*~chk2neg.*m_avg(2,:), xchk.*~chk2neg.*n_avg(2,:),0.5,'Color', 'b', 'LineWidth',.1); 

axis equal; %colorbar
axis([0 max(X) 0 max(Y) 0 max(Z)])
view([90, 0])
hold off

%% S3
nexttile(4)
plot(1,1)
hold on
% title("S3")
quiver3(X,Y,Z,  xchk.*chk3pos.*l_avg(3,:), xchk.*chk3pos.*m_avg(3,:), xchk.*chk3pos.*n_avg(3,:),0.5,'Color', 'b', 'LineWidth',2); %Highest value, non-zero, positive
quiver3(X,Y,Z,  xchk.*chk3neg.*l_avg(3,:), xchk.*chk3neg.*m_avg(3,:), xchk.*chk3neg.*n_avg(3,:),0.5,'Color', 'r', 'LineWidth',2); %Highest value, non-zero, negative

quiver3(X,Y,Z,  xchk.*~chk3pos.*l_avg(3,:), xchk.*~chk3pos.*m_avg(3,:), xchk.*~chk3pos.*n_avg(3,:),0.5,'Color', 'r', 'LineWidth',.1); 
quiver3(X,Y,Z,  xchk.*~chk3neg.*l_avg(3,:), xchk.*~chk3neg.*m_avg(3,:), xchk.*~chk3neg.*n_avg(3,:),0.5,'Color', 'b', 'LineWidth',.1);

axis equal; %colorbar
axis([0 max(X) 0 max(Y) 0 max(Z)])
view([90, 0])
hold off

%this is for debugging the xPhys location.
% figure(99)
% scatter3(X,Y,Z, 100*(xPhys.^4));


if loop == 0
    frame1 = getframe(2);
      im = frame2im(frame1);
      [imind,cm] = rgb2ind(im,256);
      imwrite(imind,cm,strcat(filename ,'_separated.gif'),'gif', 'Loopcount',inf);
else
    frame1 = getframe(2);
      im = frame2im(frame1);
      [imind,cm] = rgb2ind(im,256);
      imwrite(imind,cm,strcat(filename ,'_separated.gif'),'gif','WriteMode','append');
end

end

function B = strainDisp(r, s, t)
G = shapeFuncGradient(r, s, t)';
B = sparse(6, 24);
% Axial strain
B(1, 1:3:24) = G(:, 1); %x
B(2, 2:3:24) = G(:, 2); %y
B(3, 3:3:24) = G(:, 3); %z
% Shear strain
B(4, 2:3:24) = G(:, 3); %yz
B(4, 3:3:24) = G(:, 2);
B(5, 1:3:24) = G(:, 3); %zx
B(5, 3:3:24) = G(:, 1);
B(6, 1:3:24) = G(:, 2); %xy
B(6, 2:3:24) = G(:, 1);
end


function G = shapeFuncGradient(r, s, t)
%Gradient of the shape function:
%1/8*(1 + cci)(1 + nni)(1 + cc) 
    G =1/8*...
        [...
        -(1-s)*(1-t), -(1-t)*(1-r), -(1-s)*(1-r);
         (1-s)*(1-t), -(1-t)*(1+r), -(1-s)*(1+r);
         (1+s)*(1-t),  (1-t)*(1+r), -(1+s)*(1+r);
        -(1+s)*(1-t),  (1-t)*(1-r), -(1+s)*(1-r);
        -(1-s)*(1+t), -(1+t)*(1-r),  (1-s)*(1-r);
         (1-s)*(1+t), -(1+t)*(1+r),  (1-s)*(1+r);
         (1+s)*(1+t),  (1+t)*(1+r),  (1+s)*(1+r);
        -(1+s)*(1+t),  (1+t)*(1-r),  (1+s)*(1-r)]';

end