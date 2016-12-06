%Joel Cheverie
%1002924393

%% File: dinoTestF
%% A3 2016 handout code
%% Estimate F matrix from synthetic corresponding points.
%%
%% ADJ

clear
close all
FALSE = 1 == 0;
TRUE = ~FALSE;

global matlabVisRoot

%% We need to ensure the path is set for the iseToolbox.
if length(matlabVisRoot)==0
  dir = pwd;
  cd '../../matlab/'   %% CHANGE THIS to your startup directory
  startup;
  cd(dir);
end

reconRoot = '.';  %% CHANGE THIS to the directory you installed A4 handout/
addpath([reconRoot '/utils']);

% Random number generator seed:
seed = round(sum(1000*clock));
rand('state', seed);
seed0 = seed;
% We also need to start randn. Use a seedn generated from seed:
seedn = round(rand(1,1) * 1.0e+6);
randn('state', seedn);

nTrial = 10;  %% Number of ransac trials to use


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up cameras
%% The cameras are automatically rotated by projectDino to fixate
%% on the mean of the 3D points.  We do not always want to allow
%% the cameras to move in this fashion (eg. when we change sclZ).
%% So we will compute the rotations for the left and right cameras
%% once and for all, and then use these.
f = 100; % focal length
dLeft = [-50, 0, -150]';  % Center of projection for left camera
dRight = [50, 0, -150]';  % Center of projection for right camera
%% Compute camera rotations to fixate on Dino's center.
[pLeft polys MintLeft MextLeft] = projectDino(f, dLeft, [], 1.0);
Rleft = MextLeft(:, 1:3);
[pRight polys MintRight MextRight] = projectDino(f, dRight, [], 1.0);
Rright = MextRight(:, 1:3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Generate data...
sclZ = 1.0;
%% Dino left image data
[pLeft polys MintLeft MextLeft] = projectDino(f, dLeft, Rleft, sclZ);

%% Dino right image data
[pRight polys MintRight MextRight] = projectDino(f, dRight, Rright, sclZ);


%% Show left and right images
hFig = figure(1); clf; 
plot(pLeft(1,:), pLeft(2,:), '.b');
axis xy; axis equal;
xlabel('X'); ylabel('Y');
axis([-150 150 -100 100]);
title('Left image of Dino');
pause(0.1);

hFig = figure(2); clf; 
plot(pRight(1,:), pRight(2,:), '.b');
axis xy; axis equal;
axis([-150 150 -100 100]);
xlabel('X'); ylabel('Y');
title('Right image of Dino');
pause(0.1);


%% Build correspondence data
clear imPts;

imPts = cat(3, pLeft, pRight);
nPts = size(imPts,2);
if size(imPts,1) == 2
  imPts = cat(1, imPts, ones(1,nPts,2));
end



%% RANSAC for F
seeds = {};
sigma = 2.0; rho = 2;
for kTrial = 1: nTrial
  %% Test out F matrix on a random sample of 8 points
  idTest = randperm(nPts);
  nTest = min(8, nPts);
  idTest = idTest(1:nTest);

  %% Solve for F matrix on the random sample
  [F Sa Sf] = linEstF(imPts(:,idTest,1), imPts(:,idTest,2),1);
  
  %% Compute perpendicular error of all points to epipolar lines
  perpErrL = zeros(1,nPts);
  for k = 1:nPts
    lk = imPts(:,k,2)' * F';
    perpErrL(k) = (lk* imPts(:,k,1))/norm(lk(1:2));
  end
  
  %% Detect inliers
  idInlier = abs(perpErrL) < rho*sigma;
  
  %% Count inliers
  nInlier = sum(idInlier);
  if nInlier > 20
    %% Store sets of sampled points with at least 20 inliers
    seed = struct;
    seed.id = idTest;
    seed.idInlier = idInlier;
    seed.nInlier = nInlier;
    seed.F = F;
    
    kSeed = length(seeds)+1
    seeds{kSeed} = seed;
  end
end 
%% Done RANSAC trials

%% Extract best solution
nInliers = zeros(1, length(seeds));
for ks = 1:length(seeds)
  nInliers(ks) = seeds{ks}.nInlier;
end 
[nM ks] = max(nInliers);
nInliers(ks)

%%  Refine estimate of F using all inliers.
F = seeds{ks}.F;
idInlier = seeds{ks}.idInlier;

idInlierOld = idInlier;
sum(idInlier)
%% Do at most 10 iterations attempting to entrain as many points as possible.
for kIt = 1: 10
  %% Fit F using all current inliers
  [F Sa Sf] = linEstF(imPts(:,idInlier,1), imPts(:,idInlier,2),1);
  
  %% Compute perpendicular error to epipolar lines
  perpErrL = zeros(1,nPts);
  for k = 1:nPts
    lk = imPts(:,k,2)' * F';
    perpErrL(k) = (lk* imPts(:,k,1))/norm(lk(1:2));
  end
  idInlier = abs(perpErrL) < rho*sigma;
  nInlier = sum(idInlier)
  
  %% If we have the same set of inliers as the previous iteration then stop.
  if all(idInlier == idInlierOld)
    break;
  end
  idInlierOld = idInlier;
end
  
%%%%%%%%%%%%%%%%%%%%% Plot results
nTest = 64;  %% Number of epipolar lines to plot
nCol = 16;   %% Number of different colours to use.
col = hsv(nCol);  %% Colour map.

%% Random sample the lines to plot
idLines = 1:floor(nPts/nTest):nPts;  
idLines = idLines(1:nTest);

%% Show left image
hFig = figure(1);
clf; hold on;
% Plot all interest point locations as blue .'s
plot(imPts(1,:,1), imPts(2,:,1), '.b');
axis([-150 150 -100 100]); axis xy; axis equal;
ax = axis;
cropBox = [ax(1) ax(3) ax(2) ax(4)];
title('Epipolar Lines in Left Image');
for kl = 1:length(idLines)
  % Plot interest point location corresponding to epipolar line as a "o"
  % in the same colour as the epipolar line.
  k = idLines(kl);
  plot(imPts(1,k,1), imPts(2,k,1), 'o', 'Color', col(mod(k,nCol)+1,:));
  % Plot epipolar line.
  lk = imPts(:,k,2)' * F';
  epk = cropLineInBox(lk(1:2), lk(3), cropBox); 
  set(line(epk(:,1), epk(:,2)), 'Color', col(mod(k,nCol)+1,:));
end

%% Show right image
hFig = figure(2);
clf; hold on;
% Plot all interest point locations as blue .'s
plot(imPts(1,:,2), imPts(2,:,2), '.b');
axis([-150 150 -100 100]); axis xy; axis equal;
ax = axis;
cropBox = [ax(1) ax(3) ax(2) ax(4)];
title('Epipolar Lines in Right Image');
perpErrR = [];
for kl = 1:length(idLines)
  % Plot interest point location corresponding to epipolar line as a "o"
  % in the same colour as the epipolar line.
  k = idLines(kl);
  plot(imPts(1,k,2), imPts(2,k,2), 'o', 'Color', col(mod(k,nCol)+1,:));
  % Plot epipolar line.
  lk = imPts(:,k,1)' * F;
  epk = cropLineInBox(lk(1:2), lk(3), cropBox); 
  set(line(epk(:,1), epk(:,2)), 'Color', col(mod(k,nCol)+1,:));
  perpErrR = [perpErrR (lk * imPts(:,k,2)/norm(lk(1:2)))];
end

%% Compute perpendicular distance to epipolar lines in left and right images.
perpErrL = [];
for k = 1:nPts
  lk = imPts(:,k,2)' * F';
  perpErrL = [perpErrL (lk * imPts(:,k,1))/norm(lk(1:2))];
end
perpErrR = [];
for k = 1:nPts
  lk = imPts(:,k,1)' * F;
  perpErrR = [perpErrR (lk * imPts(:,k,2)/norm(lk(1:2)))];
end

%% Plot a histogram of the perpendicular distances
err = [perpErrL perpErrR];
err = min(err, 10);
err = max(err, -10);
[n b] = histo(err, 64);
figure(4); clf;
plot(b,n);
title('Distance to epipolar line');
xlabel('Error in pixels');
ylabel('Frequency');

%% Count inliers
idL = abs(perpErrL)< rho*sigma;
idR = abs(perpErrR) < rho*sigma;
idInlier = idL & idR;
sum(idInlier)
sum(idInlier)/nPts

%% save 'Fcorr' F Sa Sf idInlier nInliers

%Compute F_0

%QA1
%Step 1 center the world coordinates onto MextLeft.
r_l = MextLeft(:,1:3);
t_l = MextLeft(:,4);
upper_left = r_l';
upper_right = -1 * transpose(r_l)*t_l;
bottom = [0;0;0;1]';
H = horzcat(upper_left, upper_right);
H = vertcat(H,bottom);
relative = MextRight*H;
%Step 2 compute E
R = relative(:,1:3);
t = relative(:,4);
t_x = zeros(3,3);
t_x(1,2)=-t(3);
t_x(1,3)=t(2);
t_x(2,1)= t(3);
t_x(2,3) = -t(1);
t_x(3,1)= -t(2);
t_x(3,2) = t(1);
E = t_x * R;
%Step 3 Compute F
F_0 = transpose(inv(MintRight))*E * inv(MintLeft);

%QA2
x = -150:15:150;
y = -100:10:100;
pts=zeros(3,length(x)*length(y));
count = 1;
for i=1:length(x)
    for j=1:length(y)
        pts(:,count) = [x(i);y(j);1];
        count = count+1;
    end
end


idLines = 1:floor(nPts/nTest):nPts;  
idLines = idLines(1:nTest);
%% Show left image
error = 0;
hFig = figure(5);
clf; hold on;
% Plot all interest point locations as blue .'s
%plot(pts(1,:), pts(2,:), '.b');
axis([-150 150 -100 100]); axis xy; axis equal;
ax = axis;
cropBox = [ax(1) ax(3) ax(2) ax(4)];
title('Epipolar Lines in Left Image For F_0');
errors=zeros(length(x)*length(y));
for k = 1:size(pts,2)
  % Plot interest point location corresponding to epipolar line as a "o"
  % in the same colour as the epipolar line.
  %k = idLines(kl);
  %plot(pts(1,k), pts(2,k), 'o', 'Color', col(mod(k,nCol)+1,:));
  % Plot epipolar line.
  lk = pts(:,k)' * F';
  %F_0's line.
  lk_f0 = pts(:,k)' * F_0;
  epk = cropLineInBox(lk_f0(1:2), lk_f0(3), cropBox); 
  %Ingnore if outside of the box
  if(isnan(epk(1)))
      continue
      disp('Nope')
  end
  %Otherwise consider both endpoints of the line.
  ep_1 = epk(1,:);
  ep_1 = [ep_1 1];
  ep_2 = epk(2,:);
  ep_2 = [ep_2 1];
  perpErr1=[abs(lk*ep_1')/norm(lk(1:2))];
  perpErr2=[abs(lk*ep_2')/norm(lk(1:2))];
  perpErr = max(perpErr1, perpErr2);
  errors(k)=perpErr;
  error = max(perpErr, error);
  set(line(epk(:,1), epk(:,2)), 'Color', col(mod(k,nCol)+1,:));
end

error

%% Show left image
hFig = figure(6);
clf; hold on;
% Plot all interest point locations as blue .'s
%plot(pts(1,:), pts(2,:), '.b');
axis([-150 150 -100 100]); axis xy; axis equal;
ax = axis;
cropBox = [ax(1) ax(3) ax(2) ax(4)];
title('Epipolar Lines in Left Image For F');
errors=zeros(length(x)*length(y));
for k = 1:size(pts,2)
  % Plot interest point location corresponding to epipolar line as a "o"
  % in the same colour as the epipolar line.
  %k = idLines(kl);
  %plot(pts(1,k), pts(2,k), 'o', 'Color', col(mod(k,nCol)+1,:));
  % Plot epipolar line.
  lk = pts(:,k)' * F';
  %F_0's line.
  lk_f0 = pts(:,k)' * F_0;
  epk = cropLineInBox(lk(1:2), lk(3), cropBox); 
  %Ingnore if outside of the box
  if(isnan(epk(1)))
      continue
      disp('Nope')
  end
  set(line(epk(:,1), epk(:,2)), 'Color', col(mod(k,nCol)+1,:));
end