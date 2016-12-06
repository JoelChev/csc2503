% Joel Cheverie
% 1002924393
function [L] = fitChromeSphere(chromeDir, nDir, chatty)
  % [L] = fitChromeSphere(chromeDir, nDir, chatty)
  % Input:
  %  chromeDir (string) -- directory containing chrome images.
  %  nDir -- number of different light source images.
  %  chatty -- true to show results images. 
  % Return:
  %  L is a 3 x nDir image of light source directions.

  % Since we are looking down the z-axis, the direction
  % of the light source from the surface should have
  % a negative z-component, i.e., the light sources
  % are behind the camera.
    
  if ~exist('chatty', 'var')
    chatty = false;
  end
    
  mask = ppmRead([chromeDir, 'chrome.mask.ppm']);
  mask = mask(:,:,1) / 255.0;

  for n=1:nDir
    fname = [chromeDir,'chrome.',num2str(n-1),'.ppm'];
    im = ppmRead(fname);
    imData(:,:,n) = im(:,:,1);           % red channel
  end

  % Need to determine centre and radius of 
  % the sphere. This is possible by finding the 
  % brightest pixels, then averaging the location
  % to recover x and y. 
  % The pixel that is farthest right less the 
  % pixel that is farthest left gives the diameter.
  % Divide this by 2 to get the radius.
  
  [column, row] = find(mask>0);
  centre = [mean(row), mean(column)];
  diameter = double((max(column) - min(column)));
  radius = double(diameter/2);

  
  % Now determine the direct specular part of the Phong model.
  d_e = [0 0 -1.0];
  L=zeros(3, nDir);
  for i=1:nDir
  	% First we find the brightest spot in the image.
    % This is where the most reflection is coming from.
    [bright_x, bright_y] = find(imData(:,:,i) == 255);
    reflection_point= [mean(bright_y), mean(bright_x)];
    % determine the normal at bright_x and bright_y
    norm_x = reflection_point(1) - centre(1);
    norm_y = reflection_point(2) - centre(2);
    norm_z = -sqrt(radius^2 - norm_x^2 - norm_y^2);
    unnormalized_normal = [norm_x, norm_y, norm_z];
    normal = unnormalized_normal/norm(unnormalized_normal);
    % As we are assuming
    % perfect specular reflection,
    % we can recover d_i, aka L.
    d_i =2*dot(normal,d_e)*normal - d_e;
    L(:,i)=d_i;
      if (chatty)
          figure(1); clf;
          showIm(imData(:,:,n) .* mask);
      end
  end
  return;

