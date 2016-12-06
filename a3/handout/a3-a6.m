%Joel Cheverie
%1002924393
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
dLeft = [-5, 0, -150]';  % Center of projection for left camera
dRight = [5, 0, -150]';  % Center of projection for right camera
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

%QA3

sigma_n = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1];
num_sigma_trials = 100;

%Make a grid of points.
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

%Now we do the trials.
median_errors = zeros(length(sigma_n), 1);
for i = 1:length(sigma_n)
    errors = zeros(num_sigma_trials,1);
    for j=1:num_sigma_trials
        %Add noise.
        imPts(1,:,1) = imPts(1,:,1) + normrnd(0,sigma_n(i));
        imPts(2,:,1) = imPts(2,:,1) + normrnd(0,sigma_n(i));
        imPts(1,:,2) = imPts(1,:,2) + normrnd(0,sigma_n(i));
        imPts(2,:,2) = imPts(2,:,2) + normrnd(0,sigma_n(i));
          %% Test out F matrix on a random sample of 8 points
          idTest = randperm(nPts);
          nTest = min(8, nPts);
          idTest = idTest(1:nTest);

          %% Solve for F matrix on the random sample
          [F Sa Sf] = linEstF(imPts(:,idTest,1), imPts(:,idTest,2),1);


        %%%%%%%%%%%%%%%%%%%%% Plot results
        nTest = 64;  %% Number of epipolar lines to plot
        nCol = 16;   %% Number of different colours to use.
        col = hsv(nCol);  %% Colour map.

        idLines = 1:floor(nPts/nTest):nPts;  
        idLines = idLines(1:nTest);
        %% Show left image
        
        error = 0;
        %hFig = figure(5);
        clf; hold on;
        % Plot all interest point locations as blue .'s
        plot(imPts(1,:,1), imPts(2,:,1), '.b');
        axis([-150 150 -100 100]); axis xy; axis equal;
        ax = axis;
        cropBox = [ax(1) ax(3) ax(2) ax(4)];
        title('Epipolar Lines in Left Image For F_0');
        for k = 1:size(pts,2)
          % Plot interest point location corresponding to epipolar line as a "o"
          % in the same colour as the epipolar line.
          %k = idLines(kl);
          %plot(imPts(1,k,1), imPts(2,k,1), 'o', 'Color', col(mod(k,nCol)+1,:));
          % Plot epipolar line.
          lk = pts(:,k)' * F';
          %F_0's line.
          lk_f0 = pts(:,k)' * F_0;
          epk = cropLineInBox(lk_f0(1:2), lk_f0(3), cropBox); 
          %Ingnore if outside of the box
          if(isnan(epk(1)))
              continue
          end
          %Otherwise consider both endpoints of the line.
          ep_1 = epk(1,:);
          ep_1 = [ep_1 1];
          ep_2 = epk(2,:);
          ep_2 = [ep_2 1];
          perpErr1=[abs(lk*ep_1')/norm(lk(1:2))];
          perpErr2=[abs(lk*ep_2')/norm(lk(1:2))];
          perpErr = max(perpErr1, perpErr2);
          error = max(perpErr, error);
          set(line(epk(:,1), epk(:,2)), 'Color', col(mod(k,nCol)+1,:));
        end

        errors(j)=error;
    end
    median_errors(i) = median(errors);
end
figure(3)
plot(sigma_n, median_errors)
title('Median Errors Vs. Sigma_N Plot')
xlabel('Sigma_N')
ylabel('Median Error')

figure(4)
semilogx(sigma_n, median_errors)
title('Median Errors Vs. Log(Sigma_N) Plot')
xlabel('Log(Sigma_N)')
ylabel('Median Errors')


%% Show left image
hFig = figure(5);
clf; hold on;
% Plot all interest point locations as blue .'s
%plot(pts(1,:), pts(2,:), '.b');
axis([-150 150 -100 100]); axis xy; axis equal;
ax = axis;
cropBox = [ax(1) ax(3) ax(2) ax(4)];
title('Epipolar Lines in Left Image For F_0');
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
  perpErr1
  perpErr2=[abs(lk*ep_2')/norm(lk(1:2))];
  perpErr2
  perpErr = max(perpErr1, perpErr2);
  perpErr
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


