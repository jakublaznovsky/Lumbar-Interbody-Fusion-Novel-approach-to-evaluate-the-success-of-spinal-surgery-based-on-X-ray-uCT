function [mask,mask_ind] = find_roi(I,intersection,centroid_upper_x_mean,centroid_upper_y_mean,centroid_bottom_x_mean,centroid_bottom_y_mean)

%%   Synopsis
%       [mask,mask_ind] = find_roi(I,intersection,centroid_upper_x_mean,centroid_upper_y_mean,centroid_bottom_x_mean,centroid_bottom_y_mean)
%   Description
%        Returns 2D mask based on points of intersection of ellipses. To
%        change size of mask, change variables length_of_mask and width_of_mask.
%   Inputs 
%         - I                          image
%         - intersection               coordinates of intersection of two ellipses 
%         - centroid_upper_x_mean      x-coordinate of centroid of upper ellipse
%         - centroid_upper_y_mean      y-coordinate of centroid of upper ellipse
%         - centroid_bottom_x_mean     x-coordinate of centroid of bottom ellipse
%         - centroid_bottom_y_mean     y-coordinate of centroid of bottom ellipse       
%   Outputs
%        Mask is just for later displaying of result. Mask_ind defines ROI.



%%
% distance to side points 
% %d = whole distance, d2= distance from A to C
d=sqrt((intersection(1,1) - intersection(2,1))^2 + (intersection(1,2) - intersection(2,2))^2);
length_of_mask = 1600; %1600?
d2 = (length_of_mask - sqrt((intersection(1,1) - intersection(2,1))^2 + (intersection(1,2) - intersection(2,2))^2))/2;
      %900??? 1000??? pro 9A9 1600
x_c = intersection(2,1) + (d2*(intersection(2,1)-intersection(1,1))/d);
y_c = intersection(2,2) + (d2*(intersection(2,2)-intersection(1,2))/d);
% hold on;plot(x_c,y_c,'gx')

x_d = intersection(1,1) + (d2*(intersection(1,1)-intersection(2,1))/d);
y_d = intersection(1,2) + (d2*(intersection(1,2)-intersection(2,2)))/d;

% hold on;plot(x_d,y_d,'rx')


%distance to corners
d = sqrt((centroid_upper_x_mean - centroid_bottom_x_mean)^2 + (centroid_upper_y_mean - centroid_bottom_y_mean)^2);
width_of_mask = 500;%
d2 = width_of_mask/2;
     %250??? 500??? pro 9A9 800
x_e = centroid_bottom_x_mean - (d2*(centroid_bottom_x_mean-centroid_upper_x_mean)/d);
y_e = centroid_bottom_y_mean - (d2*(centroid_bottom_y_mean-centroid_upper_y_mean)/d);
u = [centroid_bottom_x_mean-x_e,centroid_bottom_y_mean-y_e];
P0 = [x_d,y_d];
P1 = [P0+u;P0-u];
P0 = [x_c,y_c];
u = [centroid_bottom_x_mean-x_e,centroid_bottom_y_mean-y_e];
P1 = floor([P1;P0-u;P0+u]);

pgon = polyshape(P1(:,1),P1(:,2),'Simplify',false);
hold on;plot(pgon)
disp(P1)

mask = poly2mask(P1(:,1),P1(:,2),size(I,1),size(I,2));

[maximum_x,maximum_pos] = max(P1(:,1));
maximum_y = P1(maximum_pos,2);
[minimum_x,minimum_pos] = min(P1(:,1));
minimum_y = P1(minimum_pos,2);



if minimum_x < 1
    minimum_x = 1;
elseif minimum_y < 1
    minimum_y = 1;
elseif maximum_x > size(I,2)
    maximum_x = size(I,2);
elseif maximum_y > size(I,1)
    maximum_y = size(I,1);
end

mask_ind = [maximum_x maximum_y;minimum_x minimum_y];

