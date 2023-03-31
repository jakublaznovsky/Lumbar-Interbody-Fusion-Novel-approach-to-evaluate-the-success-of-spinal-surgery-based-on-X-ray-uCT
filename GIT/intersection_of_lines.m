function [x1,y1,x2,y2,x_intersect, y_intersect] = intersection_of_lines(centroid_upper_x_mean,centroid_bottom_x_mean,centroid_upper_y_mean,centroid_bottom_y_mean,intersection)
%   Synopsis
%       [x1,y1,x2,y2,x_intersect, y_intersect] = intersection_of_lines(centroid_upper_x_mean,centroid_bottom_x_mean,centroid_upper_y_mean,centroid_bottom_y_mean,intersection)
%   Description
%        Returns coordinates of intersection of two lines defined by centroids of ellipses and intersections of ellipses. 
%   Inputs 
%         - centroid_upper_x_mean      x-coordinate of centroid of upper ellipse
%         - centroid_bottom_x_mean     x-coordinate of centroid of bottom ellipse
%         - centroid_upper_y_mean      y-coordinate of centroid of upper ellipse
%         - centroid_bottom_y_mean     y-coordinate of centroid of bottom ellipse
%         - intersection               coordinates of intersection of two ellipses
%   Outputs
%         - x1,y1,x2,y2                points for later displaying lines
%         - x_intersect                x-coordinate of intersection of lines 
%         - y_intersect                y-coordinate of intersection of lines 


% point for line1
x1  = [centroid_upper_x_mean centroid_bottom_x_mean];
y1  = [centroid_upper_y_mean centroid_bottom_y_mean];
% points for line2
x2 = intersection(:,1);
y2 = intersection(:,2);
% fit linear polynomial
p1 = polyfit(x1,y1,1);
p2 = polyfit(x2,y2,1);
% calculate intersection
x_intersect = fzero(@(x) polyval(p1-p2,x),3);
y_intersect = polyval(p1,x_intersect);
% display result
% line(x1,y1);
% hold on;
% line(x2,y2);
% plot(x_intersect,y_intersect,'r*')


x12 = centroid_upper_x_mean-x_intersect; y12 = centroid_upper_y_mean-y_intersect;
x32 = intersection(2,1)-x_intersect; y32 = intersection(2,2)-y_intersect;
ang = atan2(abs(x12*y32-x32*y12),x12*x32+y12*y32);
ang = ang*180/pi;
if (ang < 35) | ((180 - ang) < 35)
%     disp('Detection of vertebraes was not succesfull.')
    error('Detection of vertebraes was not succesfull. Computation must be moved to another image')
end
end