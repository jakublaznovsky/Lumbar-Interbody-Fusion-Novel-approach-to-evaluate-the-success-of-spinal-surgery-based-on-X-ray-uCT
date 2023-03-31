%% This is the main script, for computation of the Lumbar Interbody Fusion ration. The script is based on the data from X-ray micro Computed Tomography and porcine vertebrae.
clc;clear all;close all
warning off;

%% LOAD ALL SLICES FROM .VOL FILE:
fid = fopen('21A8.vol');
voxels = fread(fid, '*uint16');  %read as 16 bit unsigned. I'm assuming they're unsigned
fclose(fid);
voxels = reshape(voxels, [1554, 1393, 2020]);  %Reshape the data, according to the .pcr file
array=imrotate3(voxels, 90,[1 0 0]); %Rotate the data, according to  top-cranialis and bottom-caudalis (ventralis, dorsalis) orientation
array(:,:,1:700)=[];
array(:,:,650:end)=[];

%% BINARIZATION

pom = array(:,:,(floor(size(array,3)/2)));  %7A8,9A8
thresh = graythresh(array(:,:,200))%(pom);
binar = imbinarize(array,thresh); %array_bin=binar


%% LOCATE VERTEBRAE IN 2D

% ind = round(size(array,3)/2-40); %Find a cross-section suitable for the detection of the LIF area
% ind = [ind,ind]; %in order to go through both directiona
% i = 1;
% while 1;
%     I = binar(:,:,ind(i)); % Binarize the image
%     I_orig = array(:,:,ind(i)); 
%     [I2,xc,yc,imgfill] = detect_concave_and_convex_points(I); %Detection of convex and concave points
%     [centroid_upper_x_mean,centroid_upper_y_mean,centroid_bottom_x_mean,centroid_bottom_y_mean] = find_quadrilaterals(I_orig,xc,yc); %approximate the vertebrae by quadrilaterals
%     if any(isnan([centroid_upper_x_mean,centroid_bottom_x_mean])) %Find centroids of both vertebrae
%         disp("At least one vertebrae wasn't detected.") %Case, if the centroid was not detected
%         if i == 2 || ind(1) == ind(2)
%             ind = [ind(1)+100,ind(2)-100]
%             i = 1;
%         else
%             i = 2;
%         end
%     else
%         break
%     end
%     if ind(2) < 1 || ind(1) > size(array,3)
%         disp("Dection of vertebraes wasn't sucessfull.")
%         break
%     end
% end
% 
% % compute the reguired rotation of ellipses 
% [rot_ellipse_upper_x,rot_ellipse_upper_y,rot_ellipse_bottom_x,rot_ellipse_bottom_y,mask] = compute_rotated_ellipses(imgfill,centroid_upper_x_mean,centroid_upper_y_mean,centroid_bottom_x_mean,centroid_bottom_y_mean);
% % compute the intersection of ellipses
% intersection = intersection_of_ellipses(rot_ellipse_upper_x,rot_ellipse_bottom_x,rot_ellipse_upper_y,rot_ellipse_bottom_y);
% % compute the intersections of two lines
% [x1,y1,x2,y2,x_intersect, y_intersect] = intersection_of_lines(centroid_upper_x_mean,centroid_bottom_x_mean,centroid_upper_y_mean,centroid_bottom_y_mean,intersection);
% % display results
% figure
% imshow(I_orig,[])
% hold on
% plot(centroid_upper_x_mean,centroid_upper_y_mean,'b*')
% plot(centroid_bottom_x_mean,centroid_bottom_y_mean,'b*')
% plot(rot_ellipse_upper_x,rot_ellipse_upper_y,'b','LineWidth',2);
% plot(rot_ellipse_bottom_x,rot_ellipse_bottom_y,'b','LineWidth',2);
% plot(intersection(:,1),intersection(:,2),'rx','MarkerSize',15,'LineWidth',3)
% line(x1,y1,'LineWidth',2);
% line(x2,y2,'LineWidth',2);
% plot(x_intersect,y_intersect,'r*')
% hold off

%% LOCATE VERTEBRAE IN 2D

ind = round(size(array,3)/2);
ind = [ind,ind]; %in order to go through both direction
i = 1;
while 1;
    I = binar(:,:,ind(i)); % Binarize the image
    I(:,1:round(size(I,2)/3))=0; %Crop the spinous process
    I_orig = array(:,:,ind(i)); %Select one crossśection from the original data
    [I2,xc,yc,imgfill] = detect_concave_and_convex_points(I); %Detection of convex and concave points
    [centroid_upper_x_mean,centroid_upper_y_mean,centroid_bottom_x_mean,centroid_bottom_y_mean] = find_quadrilaterals(I_orig,xc,yc); %approximate the vertebrae by quadrilaterals
    if any(isnan([centroid_upper_x_mean,centroid_bottom_x_mean])) %Find centroids of both vertebrae
        disp("At least one vertebrae wasn't detected.") %Case, if the centroid was not detected
        if i == 2 || ind(1) == ind(2)
            ind = [ind(1)+100,ind(2)-100]
            i = 1;
        else
            i = 2;
        end
    else
        % compute the reguired rotation of ellipses 
        [rot_ellipse_upper_x,rot_ellipse_upper_y,rot_ellipse_bottom_x,rot_ellipse_bottom_y,mask] = compute_rotated_ellipses(imgfill,centroid_upper_x_mean,centroid_upper_y_mean,centroid_bottom_x_mean,centroid_bottom_y_mean);
        if sum(sum(I2(mask==1)))*100/sum(sum(mask)) < 40 %area of vertebrae in ellipses is less than 40%
            disp("At least one vertebrae wasn't detected.")
            if i == 2 || ind(1) == ind(2)
                ind = [ind(1)+100,ind(2)-100]
                i = 1;
            else
                i = 2;
            end
        else
            break
        end
    end
    if ind(2) < 1 || ind(1) > size(array,3)
        disp("Dection of vertebraes wasn't sucessfull.")
        break
    end
end

% compute the intersection of ellipses
intersection = intersection_of_ellipses(rot_ellipse_upper_x,rot_ellipse_bottom_x,rot_ellipse_upper_y,rot_ellipse_bottom_y);
% compute the intersections of two lines
[x1,y1,x2,y2,x_intersect, y_intersect] = intersection_of_lines(centroid_upper_x_mean,centroid_bottom_x_mean,centroid_upper_y_mean,centroid_bottom_y_mean,intersection);
% display results
figure
imshow(I,[])
hold on
plot(centroid_upper_x_mean,centroid_upper_y_mean,'rp','LineWidth',7)
plot(centroid_bottom_x_mean,centroid_bottom_y_mean,'rp','LineWidth',7)
plot(rot_ellipse_upper_x,rot_ellipse_upper_y,'b','LineWidth',4);
plot(rot_ellipse_bottom_x,rot_ellipse_bottom_y,'b','LineWidth',4);
plot(intersection(:,1),intersection(:,2),'rx','MarkerSize',15,'LineWidth',5)
line(x1,y1,'Color','b','LineWidth',4);
line(x2,y2,'Color','c','LineWidth',4);
plot(x_intersect,y_intersect,'r*','LineWidth',4)
hold off

%% FIND MASK
[mask,mask_ind] = find_roi(I,intersection,centroid_upper_x_mean,centroid_upper_y_mean,centroid_bottom_x_mean,centroid_bottom_y_mean); %Find ROI according to the detected vertebrae
figure;imshow(mask)
I(mask==0)=0;
figure;imshow(I)


[Ld] = apply_watershed(array,mask,mask_ind); %Divide both vertebrae in order to locate the LIF area
[final_stack,area,fusion] = area_ratio(Ld); %Compute the area of LIP and area of vertebral bodies

