function [centroid_upper_x_mean,centroid_upper_y_mean,centroid_bottom_x_mean,centroid_bottom_y_mean] = find_quadrilaterals(I,xc,yc)

%   Synopsis
%       [centroid_upper_x_mean,centroid_upper_y_mean,centroid_bottom_x_mean,centroid_bottom_y_mean] = find_quadrilaterals(I,xc,yc)
%   Description
%        Returns coordinates of two centroids computed from centroids of quadrilaterals approximating upper and
%        bottom vertebrae.
%   Inputs 
%         - I       image
%         - xc      x-coordinates of detected corners
%         - yc      y-coordinates of detected corners
%   Outputs
%        Coordinates of two centroids.


% find approximating quadrilaterals
pom_x = zeros(1,4);
pom_y = zeros(1,4);
quadrilateral_x = [];
quadrilateral_y = [];

for i=1:length(xc)
    for j=1:length(xc)
        dist = sqrt((xc(i) - xc(j))^2 + (yc(i) - yc(j))^2);
        if dist <= 520 && dist >= 370                                       %>=400
            pom_x(1) = xc(i);
            pom_x(2) = xc(j);
            pom_y(1) = yc(i);
            pom_y(2) = yc(j);
            for k=1:length(xc)
                dist = sqrt((xc(j) - xc(k))^2 + (yc(j) - yc(k))^2);
                if dist <= 870 && dist >= 630
                    pom_x(3) = xc(k);
                    pom_y(3) = yc(k);
                    x12 = pom_x(1)-pom_x(2); y12 = pom_y(1)-pom_y(2);
                    x32 = pom_x(3)-pom_x(2); y32 = pom_y(3)-pom_y(2);
                    ang = atan2(abs(x12*y32-x32*y12),x12*x32+y12*y32); 
                    ang = ang*180/pi;
                    if ang >= 70 && ang <= 110
                        for l=1:length(xc)
                            dist = sqrt((xc(k) - xc(l))^2 + (yc(k) - yc(l))^2);
                            if dist <= 600 && dist >= 400
                                pom_x(4) = xc(l);
                                pom_y(4) = yc(l);
                                dist = sqrt((xc(i) - xc(l))^2 + (yc(i) - yc(l))^2);
                                x12 = pom_x(1)-pom_x(4); y12 = pom_y(1)-pom_y(4);
                                x32 = pom_x(3)-pom_x(4); y32 = pom_y(3)-pom_y(4);
                                ang = atan2(abs(x12*y32-x32*y12),x12*x32+y12*y32);   
                                ang = ang*180/pi;
                                if dist <= 900 && dist >= 620 && ang >= 70 && ang <= 110
                                    quadrilateral_x = [quadrilateral_x;pom_x];
                                    quadrilateral_y = [quadrilateral_y;pom_y];
                                   
                                end
                            end
                        end
                    end
                    
                end              
            end
        end
    end
    pom_x = zeros(1,4);
    pom_y = zeros(1,4);
end

% remove duplicate coordinates from quadrilaterals
cond = 1;
i = 1;
j = 1;
pom = sort(quadrilateral_x,2);
while cond == 1
    if j > size(quadrilateral_x,1)
        i = i+1;
        j = i;
    end
    if i ~= j & isequal(pom(i,:),pom(j,:))
        pom(j,:) = 0;
        i = i+1;
        j = j+1;
    else
        j = j+1;
    end
    if i > size(quadrilateral_x,1)
        cond = 0;
    end
end
quadrilateral_x = quadrilateral_x(find(~all(pom==0,2)),:);
quadrilateral_y = quadrilateral_y(find(~all(pom==0,2)),:);

%%compute 2 centroids of the appropriate quadrilaterals
pgon_x = [];
pgon_y = [];
centroid_upper_x = [];
centroid_upper_y = [];
centroid_bottom_x = [];
centroid_bottom_y = [];
for i=1:size(quadrilateral_x,1)
    pom = polyshape(quadrilateral_x(i,:),quadrilateral_y(i,:));
    if pom.NumRegions == 1 & area(pom) > 65000
        pgon_x = [pgon_x; pom.Vertices(:,1)'];
        pgon_y = [pgon_y; pom.Vertices(:,2)'];
%     center of quadrilateral
        if any(pom.Vertices(:,2)' < size(I,1)/5)  %the center of upper vertebrae is computing from quadrilaterals which have at least one of the vertices < size(BW,1)/7
            [x,y] = centroid(pom);
            centroid_upper_x = [centroid_upper_x;x];
            centroid_upper_y = [centroid_upper_y;y];      
        elseif any(pom.Vertices(:,2)' > (size(I,1) - size(I,1)/5)-30) %the center of bottom vertebrae is computing from quadrilaterals which have at least one of the vertices > (size(BW,1) - size(BW,1)/7)-30)
            [x,y] = centroid(pom);
            centroid_bottom_x = [centroid_bottom_x;x];
            centroid_bottom_y = [centroid_bottom_y;y];
        end
    end   
end
%% output
centroid_upper_x_mean = sum(centroid_upper_x)/size(centroid_upper_x,1);
centroid_upper_y_mean = sum(centroid_upper_y)/size(centroid_upper_y,1);
centroid_bottom_x_mean = sum(centroid_bottom_x)/size(centroid_bottom_x,1);
centroid_bottom_y_mean = sum(centroid_bottom_y)/size(centroid_bottom_y,1);

end