% Joel Cheverie
% 1002924393
close all;
fclose all;
clear;
clc;

% Set up paths and constants.

% Need to change this to refer to your startup.m 
run '../../../matlab/startup.m';

% Need to change this to be the install directory for
% the A2 handout code (with subdirectories circleFit, and util,
% and sibling directory ../data).
hmdir = '../cellFinder/';

% Link the other handout directories for A2
addpath([hmdir,'circleFit/']);
addpath([hmdir 'util/']);
addpath([hmdir '../data/']);
addpath(['../../../matlab/utvisToolbox/edge/']); % for cannyEdgels.m

cd([hmdir]);

%close all
clear

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters you need to set.
%40 seems to work best.
sigmaGM = 40; % Geman-McLure estimator scale parameter

% Additional parameters (you do not need to set)
numGuesses = 15;   % Number of circles proposals, if possible
maxNumTrials = 10; % Max number of propose-fit-check loops
minNumPts = 10;    % Min number of data edgels to bother fitting

% Display parameters
showEdgelResults = false;  % Show edgel norm and orientation images
demoRobustConv = false;  % Display robust fitting on each proposal.
interactive = true;  % Set to false to run in batch

% Names of test images for which we have stored background results.
% Here y,x = 0/1,0/1 denotes which quadrant of the original image.
fname = {'cellGFP2a_00', 'cellGFP2a_10', 'cellGFP2a_01', 'cellGFP2a_11'};
nFname = length(fname);  % Number of test images


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop over cell images

for iFname = 2:nFname-1  % Use 1 image for now, change to 1:nFname to process all

  fileName = fname{iFname};
  ['Running file: ' fileName]

  % Load image and connected component labels.
  im = pgmRead([fileName '.pgm']);
  ccomps = pgmRead([fileName '_conCom.pgm']);
  % Make component labels start at 0 instead of 1.
  ccomps = ccomps-1;
  imageBox = [1 1; size(ccomps,2) size(ccomps,1)];
  
  % Generate mask of detected cell blobs.
  cind = (ccomps > 0);

  % Show image and detected cell blobs
  figure(1); clf; 
  displayImage(im);
  sizeIm = size(im);
  figure(2); clf; 
  displayImage(double(cind));
  pause(1);

  % Echo number of components
  Nccs = max(max(ccomps))

  %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Extract edgels from the mask image cind

  sigma = 1.25;  % For the Gaussian used in the Canny edge detector.
  minStrength = 20;  % Minimum edgel strength
  [eIm nIm dIm] = cannyEdgels(double(cind)*255, sigma, minStrength);
  % eIm is 0/1 with 1 marking positions of edgels,
  % nIm is the strength of the edge (2-norm of gradient)
  % dIm * pi = theta, the edgel angle, where edgel
  %   normal = [cos(theta), sin(theta)].

  % Use log edge strength image (add one to avoid divide by zero).
  lognIm = log(nIm+1);

  maskIm = nIm > minStrength;
  dIm(isnan(dIm)) = -1;
  
  if showEdgelResults
    % Show Canny results
    figure(4);  displayImage(maskIm.*lognIm);  
    title('log gradient magnitude');
    figure(5); clf; displayImage(double(maskIm .* eIm));  
    title('edgel locations');
    figure(6); clf; showIm(maskIm .* dIm, [0 2]); 
    title('edgel orientation');
    if interactive
      fprintf(2, 'Hit enter to continue...');
      pause;
      fprintf(2, '\n');
    end

    % Close figures with results of edge detection
    close 4 5 6
  end


  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Let's fit the circles to find the disk-shaped cells

  circleEstimates = zeros(3, maxNumTrials, Nccs);
  nFound = zeros(Nccs,1);

  % You can switch between looping through all connected
  % components versus interactively selecting components
  % by switching the comments on the following two loop headers.
  
  mouseSelect = false;
  for ccInd = 1:Nccs
  
  %mouseSelect = true;
  %while mouseSelect
    
    if mouseSelect
      % Choose a connected component by selecting it with a mouse
      % Show segments in a random colour image.
      kMap = 64;
      cmap = [0 0 0; jet(kMap)];
      c = mod(randperm(Nccs)-1, kMap) + 1;
      qqc = zeros(size(ccomps));
      qqc(cind) = c(ccomps(cind));

      figure(2);
      displayImage(qqc);
      colormap(cmap);
      pause(0.1);

      figure(2); 
      fprintf('Select a component...');
      x = round(ginput(1));
      fprintf('\n');
      if size(x,1) > 0 && ptInCropBox(x, imageBox);
        ccInd = ccomps(x(2), x(1));
      else
        break;
      end
      if ccInd == 0
        break;
      end
    end
    
    % Connected component number ccInd has been selected.
    nFound(ccInd) = 0;
    cc = (ccomps == ccInd);
    
    % Build a crop box around the selected component, for display.
    [pp qq] = find(cc > 0.5);
    ccBox = [min(qq), min(pp); max(qq), max(pp)] + 32*[-1 -1; 1 1];
    % to map from window coords to original image coords (for display fns)
    cropOffset = [1; 1] - ccBox(1,:)';
    
    
    % Get the edgels and orientations within the crop box
    m = (nIm > minStrength) .* cc;
    % log nIm is for output visualization
    lognImCrop = cropImage(lognIm.*m, ccBox, 1);
    [p theta] = getEdgelData(eIm.*m, dIm.*m, ccBox);
    normals = [cos(theta), sin(theta)];
    Npts = size(p,1);
    
    % Search for plausible circles (Guess, fit, check)
    goodCircle = false;
    trialInd = 0;
    while ((trialInd < maxNumTrials) & (Npts > minNumPts))
      
      % As circles are found, edgels are removed from the data
      % set.  As a result, p, theta, and normals become shorter.
      
      % Get proposals for possible circles
      circles = getProposals(normals, p, numGuesses);
      
      if demoRobustConv
        % Show fitting behaviour on each initial guess.
        for kCirc = 1:size(circles,1)
          % Robustly fit a circle based on the initial guess
          % provided by each proposal.
          figure(13); clf;
          displayImage(lognImCrop);
          hold on;
          title('Convergence Demo (Guess (g) and Robust fit (r))');
          xGuess = circles(kCirc, 1:2)';
          rGuess = circles(kCirc, 3);
          plotCircle(xGuess+cropOffset, rGuess, 'g');
          plot(p(:,1)+cropOffset(1), p(:,2)+cropOffset(2), 'ob');
          [x0, r, w, maxW] = fitCircleRobust(p, xGuess, rGuess,...
                                             normals, sigmaGM);
          if abs(r) < max(imageBox(:))
            ptLabels = floor(w * 4/maxW);   % Label weights by quartile
            displayFittedPoints(p, ptLabels, cropOffset);
            plotCircle(x0+cropOffset, r, 'r');
          end
          fprintf(2, 'Hit enter to continue...');
          pause;
          fprintf(2, '\n');
        end
      end

      % Select the best proposal.
      circle = bestProposal(circles, sigmaGM, normals, p);
      
      if size(circle,1) > 0
        % Robustly fit a circle based on the initial guess
        % provided by the best proposal.
        [x0, r, w, maxW] = fitCircleRobust(p, circle(1:2)', circle(3), ...
                                           normals, sigmaGM);
        
        % Rescale weights to have a max of  one.
        w = w/maxW;  
        maxW = 1;
        
        % Test for accepting robust circle fit
        ptLabels = floor(w * 4);   % Label weights by quartile
        goodCircle = isGoodCircle(x0, r, w, ...
                          circleEstimates(:,:, ccInd), nFound(ccInd));
      else
        goodCircle = false;
      end
      
      % Display results and remove fitted edgels (if any) from data set
      if goodCircle
        nFound(ccInd) = nFound(ccInd) + 1;
        circleEstimates(:,nFound(ccInd),ccInd) = [x0; r];
        
        if (nFound(ccInd) == 1)
          figure(11); clf; displayImage(lognImCrop); hold on;
          title('Robust fit');
        end
        figure(12); clf;
        displayImage(lognImCrop);  title('Fit Weights'); hold on;
        displayFittedPoints(p, ptLabels, cropOffset);
        figure(11);
        plotCircle(x0+cropOffset, r, 'r');
        figure(1); hold on;
        plotCircle(x0, r, 'r');
        pause(1);
        
        % Remove edgels in newly fit circle.
        pInd = (ptLabels > 0);
        Npts = sum(~pInd);
        if (Npts > minNumPts)
          p = p(~pInd, :);
          theta = theta(~pInd, :);
          normals = normals(~pInd,:);
        end
      
      end % good circle found
        
      trialInd = trialInd + 1;
    end % Guess, Fit, Check while loop
    
    pause(1);
  end  % Looping over components (either while mouseSelect or over 1:Nccs)

  %%
  % Display all the fitted circles

  figure(1);
  clf;
  displayImage(im);
  hold on;

  for ccInd = 1:Nccs
    for circleInd = 1:nFound(ccInd)
      x0 = circleEstimates(1:2,circleInd,ccInd);
      r = circleEstimates(3,circleInd,ccInd);  
      plotCircle(x0, r, 'r');
    end
  end

  % For a postscript copy of figure 1 uncomment
  % the following line.
  % print('-depsc', fileName);

  % Pause before the next test image.
  if ~mouseSelect & interactive
    fprintf(2, 'Hit enter to continue...');
    pause;
    fprintf(2, '\n');
  end

end  % for loop over each data set.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

