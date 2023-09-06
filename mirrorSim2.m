%this produces all of the values for a spherical mirror, using the constraints of the mirror for maximum values
%truth matrix tells what values of mirrors with dimensions a & b work, 1 means it fits requirements and 2 means it doesn't
%all matrixes are based on the values that would come from a in the first row and b in the first column

increment = 0.005; %1/2 cm precision

%mirror data
aMax = 0.915; % max diameter is 1.83m so max aperature is 0.915m
bMax = 0.915; % maximum height of mirror is 3.05m, but if b > a mirror starts to curve inwards
rMax = 6.10; % r = 2f and maximum f is 3.15 meters as that is the launch profile


%satellite data
height = 222000; %satellite is at 222km
resolutionMax = 2.5; %max resolution is 2.5 m, but can be less

%Glass values, change for different glass
indexRefraction = 1.47;
massDensity = 2230;
MatCost = 1.41;
SurfCost = 1.9;

%launch properties
costKiloLaunch = 5000; %price varies, but good for choosing final mirror

%ccd values
powerMin = 1.07*10^-9; %wattage needed to be recieved
pixelLength = 5.5*10^-6; %meters
CCD_Cost =  5050;


%calculating a, b, radius, and area. If radius is greater than max, remove it from truth table
b = [0:increment:bMax]; %create a vector of b values
a = [0:increment:aMax]; %create a vector of a values

truth = ones(length(b),length(a));
radius = zeros(length(b),length(a));
area = zeros(length(b),length(a));

for i=2:length(b)
  for j=2:length(a)
    radius(i,j) =(a(j)^2 + b(i)^2)*((2*b(i))^(-1));
    area(i,j) = (a(j)^2 + b(i)^2)*pi;
    if radius(i,j) > rMax
      truth(i,j) = 0;
    endif
    truth(1,1) = 0;
    truth(1,j) = 0;
    truth (j,1) = 0;
  endfor
endfor


%calculating focal length, image distance, and magnification
focalLength = ones(length(b),length(a));
imageDistance = zeros(length(b),length(a));
magnification = zeros(length(b),length(a));

for i=2:length(b)
  for j=2:length(a)
    focalLength(i,j) = radius(i,j)/2;
    imageDistance(i,j) = radius(i,j)/2;
    magnification(i,j) = -imageDistance(i,j)/height;
  endfor
endfor


%calculating wattage recieved by mirror for CCD
Wattage = zeros(length(b),length(a)); %wattage matrix

for i=2:length(b)
  for j=2:length(a)
    Wattage(i,j) = ((area(i,j)*100)/(pi * height*height));
    if(Wattage(i,j) < powerMin);
      truth(i,j) = 0;
    endif
  endfor
endfor


%calulating resolution and removing values that do not fit
resolution = zeros(length(b),length(a));

for j=2:length(a)
    resolution(:,j) = (115.8/(2*1000*a(j)))*(1/3600)*(pi/180)*height;
    if(resolution(1,j) > resolutionMax);
      truth(:,j) = 0;
    endif
endfor


%calculating ccd resolution
CCDResolution = zeros(length(b),length(a));

for i=2:length(b)
  for j=2:length(a)
    CCDResolution(i,j) = -3*pixelLength * magnification(i,j)^-1;
    if(CCDResolution(i,j) > resolutionMax);
      truth(i,j) = 0;
    endif
  endfor
endfor


%calculate glass values

thickness = zeros(length(b),length(a));
weight = zeros(length(b),length(a));
GlassCost = zeros(length(b),length(a));

for i=2:length(b)
  for j=2:length(a)
    thickness(i,j) = (0.0951*a(j)*massDensity)*10^(-4);
    weight(i,j) = 9.81*(massDensity*thickness(i,j)*area(i,j));
    GlassCost(i,j) = MatCost*area(i,j)*thickness(i,j) + (110/indexRefraction)*SurfCost*area(i,j);
  endfor
endfor


%calculate cost
cost = GlassCost + CCD_Cost*ones(length(b),length(a)) + weight*costKiloLaunch;

costTruth = cost.*truth;

for i=1:length(b)
  for j=1:length(a)
    if costTruth(i,j) == 0
      costTruth(i,j) = inf;
    endif
  endfor
endfor


%find cheapest mirror
[minimumCost,I] = min(costTruth);
[minimumCost,J] = min(minimumCost);


cheapest = [I(J),J]
cheapestA = a(J);
cheapestB = b(I(J));
cheapestCost = minimumCost
cheapestRadius = radius(I(J),J);
cheapestThickness = thickness(I(J),J);
cheapestImageDistance = imageDistance(I(J),J);
cheapestFocalLength = focalLength(I(J),J);
cheapestWeight = weight(I(J),J);


cheapestMirror = [cheapestA  cheapestRadius cheapestB cheapestThickness cheapestImageDistance cheapestFocalLength cheapestCost cheapestWeight]

