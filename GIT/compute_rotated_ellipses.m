 function [rot_ellipse_upper_x,rot_ellipse_upper_y,rot_ellipse_bottom_x,rot_ellipse_bottom_y,mask] = compute_rotated_ellipses(I,centroid_upper_x_mean,centroid_upper_y_mean,centroid_bottom_x_mean,centroid_bottom_y_mean)

%   Synopsis
%       [rot_ellipse_upper_x,rot_ellipse_upper_y,rot_ellipse_bottom_x,rot_ellipse_bottom_y] = compute_rotated_ellipses(I,centroid_upper_x_mean,centroid_upper_y_mean,centroid_bottom_x_mean,centroid_bottom_y_mean)
%   Description
%        The function returns coordinates of ellipses aproximating upper
%        and bottom vertebrae. Firstly it computes y-radius of required
%        ellipses, X-radius is given by hand. Then it computes reguired angle of rotation between vertical  
%        orientated ellipse and required oriented ellipse. After that,
%        coordinates of rotated ellipses are returned.
%   Inputs 
%         - I                         image
%         - centroid_upper_x_mean     x-coordinate of upper vertebrae centroid
%         - centroid_upper_y_mean     y-coordinate of upper vertebrae centroid
%         - centroid_bottom_x_mean    x-coordinate of bottom vertebrae centroid
%         - centroid_bottom_y_mean    y-coordinate of bottom vertebrae centroid
%   Outputs
%        Coordinates of ellipses aproximating upper and bottom vertebrae.



% compute y_radius
canny = edge(I, 'Canny');
m = (centroid_bottom_y_mean-centroid_upper_y_mean)/(centroid_bottom_x_mean-centroid_upper_x_mean);
N = 200000 ; 
x = linspace(1,size(I,2),N);  
y = centroid_upper_y_mean+m*(x-centroid_upper_x_mean);  
x = round(x);
y = round(y);
x= x(y>0 & y < size(canny,1));
y=y(y>0 & y < size(canny,1));
pom = [];
j = 1;
for i=1:length(x)
    if canny(y(i),x(i)) == 1
        pom(j,1) = x(i);
        pom(j,2) = y(i);
        j = j+1;
    end
end
pom = unique(pom, 'rows');
if size(pom,1) > 2
    [~,ind1] = min(pom(:,2));
    [~,ind2] = max(pom(:,2));
    pom = [pom(ind1,:);pom(ind2,:)];
end
if size(pom,1) == 0 
    Distance = size(canny,1);
elseif size(pom,1) == 1 & pom(1,2) > (size(canny,1) - size(canny,1)/3)
    x = x(y == 1 );
    x = x(1,1);
    Distance =  sqrt((pom(1,1)-x)^2 + (pom(1,2)-1)^2);  
elseif size(pom,1) == 1 & pom(1,2) < size(canny,1)/3
    y = max(y);
    Distance = sqrt((pom(1,1)-x(y==910))^2 + (pom(1,2)-y)^2);  
elseif size(pom,1) == 2
    Distance = sqrt((pom(1,1)-pom(2,1))^2 + (pom(1,2)-pom(2,2))^2);
    if Distance < size(canny,1)/1.4
        Distance = size(canny,1);
    end
end

% Distance
% pom
% figure
% imshow(canny)
% hold on
% plot(x,y,'.r')
% plot(pom(:,1),pom(:,2),'gx')
y_radius = Distance/4 + 20;
x_radius = 300; %x_radius is given by hand


% compute vertiacal orientated ellipses
ellipse_upper_x = zeros(1000,1);
ellipse_upper_y = zeros(1000,1);
ellipse_bottom_x = zeros(1000,1);
ellipse_bottom_y = zeros(1000,1);
theta = linspace(0,2*pi,1000);
for k = 1:1000
        ellipse_upper_x(k) = x_radius * cos(theta(k));
        ellipse_upper_y(k) = y_radius * sin(theta(k));
        ellipse_bottom_x(k) = x_radius * cos(theta(k));
        ellipse_bottom_y(k) = y_radius * sin(theta(k));
end

% compute the required angle of rotation
x1 = [centroid_upper_x_mean, centroid_bottom_x_mean];
y1 = [centroid_upper_y_mean, centroid_bottom_y_mean];
x2 = [centroid_bottom_x_mean,  centroid_bottom_x_mean];
y2 = [centroid_bottom_y_mean, centroid_bottom_y_mean - 10];
difference = (atan((y1(2)-y1(1))/(x1(2)-x1(1))) - atan((y2(2)-y2(1))/(x2(2)-x2(1))));
R  = [cos(difference) -sin(difference); ...
      sin(difference)  cos(difference)];
rCoords_upper = R*[ellipse_upper_x' ; ellipse_upper_y'];   
rCoords_bottom = R*[ellipse_bottom_x' ; ellipse_bottom_y'];   
% compute required rotated ellipses
rot_ellipse_upper_x = rCoords_upper(1,:)' + centroid_upper_x_mean;
rot_ellipse_upper_y = rCoords_upper(2,:)' + centroid_upper_y_mean;  
rot_ellipse_bottom_x = rCoords_bottom(1,:)' + centroid_bottom_x_mean;
rot_ellipse_bottom_y = rCoords_bottom(2,:)' + centroid_bottom_y_mean;  

% display the result
% figure
% imshow(I)
% hold on
% plot(centroid_upper_x_mean,centroid_upper_y_mean,'b*')
% plot(centroid_bottom_x_mean,centroid_bottom_y_mean,'r*')
% plot(ellipse_upper_x+centroid_upper_x_mean,ellipse_upper_y+centroid_upper_y_mean,'r','LineWidth',2);
% plot(ellipse_bottom_x+centroid_bottom_x_mean,ellipse_bottom_y+centroid_bottom_y_mean,'r','LineWidth',2);
% plot(rot_ellipse_upper_x,rot_ellipse_upper_y,'b','LineWidth',2);
% plot(rot_ellipse_bottom_x,rot_ellipse_bottom_y,'b','LineWidth',2);
% hold off

% mask_bottom = poly2mask(rot_ellipse_bottom_x,rot_ellipse_bottom_y,size(I,1),size(I,2));
% mask_upper = poly2mask(rot_ellipse_upper_x,rot_ellipse_upper_y,size(I,1),size(I,2));
% mask = mask_bottom+mask_upper;
mask=0;

end