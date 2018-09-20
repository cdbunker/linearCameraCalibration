%This script uses linear camera calibration to calibrate a camera
%given n 3D to 2D point correspondeces

%Read an image and keep only the calibration grid
im1 = imread('rig1.jpg');
im1 = im1(:,1:530,:);

%Set up matrix with 3D points

points3d = 5*[0 0 0; 5 0 0; 0 5 0; 0 0 5; 5 5 0; 5 0 5; 0 5 5;
    3 6 0; 6 3 0; 3 0 6; 1 5 0; 5 0 1; 0 1 5]'; %5 cm grid spacing

one = ones(1,length(points3d));
points3dh = [points3d; one]; %homogoneous coordinates

%select corresponding 2D points (done using ginput)
x = [289.5099;
  219.4109;
  365.5495;
  283.5693;
  297.2327;
  211.6881;
  359.6089;
  342.3812;
  249.1139;
  240.7970;
  351.8861;
  218.8168;
  297.8267];
y = [353.0743
  375.0545;
  363.7673;
  269.3119;
  389.3119;
  282.9752;
  273.4703;
  381.5891;
  388.7178;
  258.6188;
  369.1139;
  356.6386;
  270.5000];

%Build P matrix
P = [];
zero_vector = [0, 0, 0, 0];
for i=1:length(points3d)
    point = points3dh(:,i);
    p1 = [point', zero_vector, -x(i)*point'];
    p2 = [zero_vector, point', -y(i)*point'];
    P = [P;p1;p2];
end

%eigenvector corresponding to smallest eigenvalue is our solution
[V, D] = eig(P'*P);
M = V(:,1);
M = reshape(M,4,3)';

%reproject 3D points and calculate error
reprojections = M*points3dh;
reprojections = bsxfun(@rdivide,reprojections,reprojections(3,:));
xyh = [x'; y'; one];
error = sum(sum((reprojections-xyh).^2))/length(xyh);
fprintf('\nMean Squared Error: %f\n\n', error)

%Plot original and reprojected points
figure
imshow(im1);
hold on
plot(x,y, 'r+','markers',12)
plot(reprojections(1,:),reprojections(2,:), 'go','markers',12)
title('True Point Locations: +; Reprojected Point Locations: O')

%Plot the i,j,k vectors
%i basis
figure
imshow(im1)

p1 = [0,0,0,1]';
p1_projected = M*p1;
p1_projected = bsxfun(@rdivide,p1_projected,p1_projected(3,:));
p2 = [50,0,0,1]';
p2_projected = M*p2;
p2_projected = bsxfun(@rdivide,p2_projected,p2_projected(3,:));
line([p1_projected(1), p2_projected(1)],[p1_projected(2), p2_projected(2)])

%j basis
p1 = [0,0,0,1]';
p1_projected = M*p1;
p1_projected = bsxfun(@rdivide,p1_projected,p1_projected(3,:));
p2 = [0,50,0,1]';
p2_projected = M*p2;
p2_projected = bsxfun(@rdivide,p2_projected,p2_projected(3,:));
line([p1_projected(1), p2_projected(1)],[p1_projected(2), p2_projected(2)])

%k basis
p1 = [0,0,0,1]';
p1_projected = M*p1;
p1_projected = bsxfun(@rdivide,p1_projected,p1_projected(3,:));
p2 = [0,0,50,1]';
p2_projected = M*p2;
p2_projected = bsxfun(@rdivide,p2_projected,p2_projected(3,:));
line([p1_projected(1), p2_projected(1)],[p1_projected(2), p2_projected(2)])

%Reprojected grid lines and plot on xy plane
for i=1:10
    j=10;
    p1 = [5*i,0,0,1]';
    p1_projected = M*p1;
    p1_projected = bsxfun(@rdivide,p1_projected,p1_projected(3,:));
    p2 = [5*i,5*j,0,1]';
    p2_projected = M*p2;
    p2_projected = bsxfun(@rdivide,p2_projected,p2_projected(3,:));
    hold on
    plot([p1_projected(1), p2_projected(1)],[p1_projected(2), p2_projected(2)], 'g')
end
for j=1:10
    i=10;
    p1 = [0,5*j,0,1]';
    p1_projected = M*p1;
    p1_projected = bsxfun(@rdivide,p1_projected,p1_projected(3,:));
    p2 = [5*i,5*j,0,1]';
    p2_projected = M*p2;
    p2_projected = bsxfun(@rdivide,p2_projected,p2_projected(3,:));
    hold on
    plot([p1_projected(1), p2_projected(1)],[p1_projected(2), p2_projected(2)], 'g')
end

%Reprojected grid lines and plot on xz plane
for i=1:10
    j=10;
    p1 = [5*i,0,0,1]';
    p1_projected = M*p1;
    p1_projected = bsxfun(@rdivide,p1_projected,p1_projected(3,:));
    p2 = [5*i,0,5*j,1]';
    p2_projected = M*p2;
    p2_projected = bsxfun(@rdivide,p2_projected,p2_projected(3,:));
    hold on
    plot([p1_projected(1), p2_projected(1)],[p1_projected(2), p2_projected(2)], 'r')
end
for j=1:10
    i=10;
    p1 = [0,0,5*j,1]';
    p1_projected = M*p1;
    p1_projected = bsxfun(@rdivide,p1_projected,p1_projected(3,:));
    p2 = [5*i,0,5*j,1]';
    p2_projected = M*p2;
    p2_projected = bsxfun(@rdivide,p2_projected,p2_projected(3,:));
    hold on
    plot([p1_projected(1), p2_projected(1)],[p1_projected(2), p2_projected(2)], 'r')
end

%Reprojected grid lines and plot on yz plane
for i=1:10
    j=10;
    p1 = [0,5*i,0,1]';
    p1_projected = M*p1;
    p1_projected = bsxfun(@rdivide,p1_projected,p1_projected(3,:));
    p2 = [0,5*i,5*j,1]';
    p2_projected = M*p2;
    p2_projected = bsxfun(@rdivide,p2_projected,p2_projected(3,:));
    hold on
    plot([p1_projected(1), p2_projected(1)],[p1_projected(2), p2_projected(2)], 'b')
end
for j=1:10
    i=10;
    p1 = [0,0,5*j,1]';
    p1_projected = M*p1;
    p1_projected = bsxfun(@rdivide,p1_projected,p1_projected(3,:));
    p2 = [0,5*i,5*j,1]';
    p2_projected = M*p2;
    p2_projected = bsxfun(@rdivide,p2_projected,p2_projected(3,:));
    hold on
    plot([p1_projected(1), p2_projected(1)],[p1_projected(2), p2_projected(2)], 'b')
end

title('Grid Lines Reprojected onto Image')