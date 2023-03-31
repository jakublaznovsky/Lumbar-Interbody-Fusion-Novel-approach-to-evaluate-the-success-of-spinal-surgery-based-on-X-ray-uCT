function [intersection] = intersection_of_ellipses(ellipse_upper_x,ellipse_bottom_x,ellipse_upper_y,ellipse_bottom_y)

%   Synopsis
%       [intersection] = intersection_of_ellipses(ellipse_upper_x,ellipse_bottom_x,ellipse_upper_y,ellipse_bottom_y)
%   Description
%        Returns coordinates of intersections of two ellipses. 
%   Inputs 
%         - ellipse_upper_x      x-coordinates of upper ellipse
%         - ellipse_bottom_x     x-coordinates of bottom ellipse
%         - ellipse_upper_y      y-coordinates of upper ellipse
%         - ellipse_bottom_y     y-coordinates of bottom ellipse
%   Outputs
%        Intersections of two ellipses. In first column there are
%        x-coordinates of intersections. In second column there are
%        y-coordinates.


% compute intersections
Ax = ellipse_upper_x;
Bx = ellipse_bottom_x;
Ay = ellipse_upper_y;
By = ellipse_bottom_y;

TMP = bsxfun(@(x,y) abs(x-y), Ax(:), reshape(Bx,1,[]));
[Dx, idxB] = min(TMP,[],2) ;
Resultx = Bx(idxB);
Dx = Dx';
Dy = zeros(1,length(Ay));

for i=1:length(Ax)
    Dy(1,i) = abs(Ay(i) - By(idxB(i)));
end

res = Dx + Dy;
[~, index] = sort(res);
res_Ax = Ax(sort(index(1:3)));
res_Ay = Ay(sort(index(1:3)));
res_Bx = Bx(idxB(sort(index(1:3))));
res_By = By(idxB(sort(index(1:3))));

Ax = ellipse_bottom_x;
Bx = ellipse_upper_x;
Ay = ellipse_bottom_y;
By = ellipse_upper_y;

TMP = bsxfun(@(x,y) abs(x-y), Ax(:), reshape(Bx,1,[]));
[Dx, idxB] = min(TMP,[],2) ;
Resultx = Bx(idxB);
Dx = Dx';
Dy = zeros(1,length(Ay));

for i=1:length(Ax)
    Dy(1,i) = abs(Ay(i) - By(idxB(i)));
end

res = Dx + Dy;
[~, index] = sort(res);
res_Ax = [res_Ax;Ax(sort(index(1:3)))];
res_Ay = [res_Ay;Ay(sort(index(1:3)))];
res_Bx = [res_Bx;Bx(idxB(sort(index(1:3))))];
res_By = [res_By;By(idxB(sort(index(1:3))))];
result = [[res_Ax;res_Bx],[res_Ay;res_By]];
[~,intersection] = kmeans(result,2);
end
