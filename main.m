%%
clc;clear all;close all
warning off;
tic

%% LOAD ALL SLICES FROM .VOL FILE:
fid = fopen('9A8.vol');
voxels = fread(fid, '*uint16');  %read as 16 bit unsigned. I'm assuming they're unsigned
fclose(fid);
voxels = reshape(voxels, [2014, 2014, 2024]);  %you may need to permute the numbers here
%and you may need to permute the dimensions afterward
array=imrotate3(voxels,90,[1 0 0]);

%% LOAD ALL SLICES AS TIFF STACK:
% file_folder = fullfile('Q:\Vojtova\IGA_pig_spine\24_10_17_impl_6_vzorku_8tyd\12A8\Image_analysis2');
% % file_folder = fullfile('U:\obratle\image_analyses_part');
% files = dir(fullfile(file_folder,'*.tif'));
% file_names = {files.name};
% t = Tiff(fullfile(file_folder,file_names{1}));
% I = read(t);
% % t = Tiff('8B16_1216.tif','r');
% % I = read(t);
% array = zeros(size(I,1),size(I,2),length(file_names),class(I));
% for i=1:length(file_names)
%     t = Tiff(fullfile(file_folder,file_names{i}));
%      t=rgb2gray(imread(t.FileName)); %optional
% %      array(:,:,i) = read(t); 
%      array(:,:,i) = t;
% end

% ULOZENI SOUBORU INPUT_CROPPED
%% Loading

% file_folder = dir('Q:\Vojtova\IGA_pig_spine\24_10_17_impl_6_vzorku_8tyd\8A8\Images_analysis2');
% CT=zeros(2130,2793,2857);
% for i=3:length(file_folder) 
%     nazev_DICOM=file_folder(i,1).name; 
%     pom=imread([nazev_DICOM]); 
% 
%     CT(:,:,i)=rgb2gray(pom);
% end
%% if images are exported in different plane, therefore:
% array = permute(array,[1 3 2]);

%% BINARIZATION

pom = array(:,:,(floor(size(array,3)/2)));  %7A8,9A8
thresh = graythresh(array(:,:,200))%(pom);
binar = imbinarize(array,thresh); %array_bin=binar
%     pom1 = histeq(pom1);           % jen pro 8A8???
%     pom1 =  im2bw(pom1,0.85);


% pomoci ekvalizace
% array1 = histeq(array);
% pom = array(:,:,(floor(size(array,3)/2)));
% thresh = graythresh(pom);
% binar = imbinarize(array1,thresh); %thresh=0.98 - nesmysl

%% LOCATE VERTEBRAES IN 2D

ind = round(size(array,3)/2);
ind = [ind,ind]; %in order to go through both directiona
i = 1;
while 1;
    I = binar(:,:,ind(i));
    I_orig = array(:,:,ind(i));
    [I2,xc,yc,imgfill] = detect_concave_and_convex_points(I);
    [centroid_upper_x_mean,centroid_upper_y_mean,centroid_bottom_x_mean,centroid_bottom_y_mean] = find_quadrilaterals(I_orig,xc,yc);
    if any(isnan([centroid_upper_x_mean,centroid_bottom_x_mean]))
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
    if ind(2) < 1 || ind(1) > size(array,3)
        disp("Dection of vertebraes wasn't sucessfull.")
        break
    end
end

% compute the reguired rotation of ellipses 
[rot_ellipse_upper_x,rot_ellipse_upper_y,rot_ellipse_bottom_x,rot_ellipse_bottom_y,mask] = compute_rotated_ellipses(imgfill,centroid_upper_x_mean,centroid_upper_y_mean,centroid_bottom_x_mean,centroid_bottom_y_mean);
% compute the intersection of ellipses
intersection = intersection_of_ellipses(rot_ellipse_upper_x,rot_ellipse_bottom_x,rot_ellipse_upper_y,rot_ellipse_bottom_y);
% compute the intersections of two lines
[x1,y1,x2,y2,x_intersect, y_intersect] = intersection_of_lines(centroid_upper_x_mean,centroid_bottom_x_mean,centroid_upper_y_mean,centroid_bottom_y_mean,intersection);
% display results
figure
imshow(I_orig,[])
hold on
plot(centroid_upper_x_mean,centroid_upper_y_mean,'b*')
plot(centroid_bottom_x_mean,centroid_bottom_y_mean,'b*')
plot(rot_ellipse_upper_x,rot_ellipse_upper_y,'b','LineWidth',2);
plot(rot_ellipse_bottom_x,rot_ellipse_bottom_y,'b','LineWidth',2);
plot(intersection(:,1),intersection(:,2),'rx','MarkerSize',15,'LineWidth',3)
line(x1,y1,'LineWidth',2);
line(x2,y2,'LineWidth',2);
plot(x_intersect,y_intersect,'r*')
hold off

%% NEW location of vertebrae(checking area of vertebrae in ellipses)

ind = round(size(array,3)/2);
ind = [ind,ind]; %in order to go through both direction
i = 1;
while 1;
    I = binar(:,:,ind(i));
    I_orig = array(:,:,ind(i));
    [I2,xc,yc,imgfill] = detect_concave_and_convex_points(I);
    [centroid_upper_x_mean,centroid_upper_y_mean,centroid_bottom_x_mean,centroid_bottom_y_mean] = find_quadrilaterals(I_orig,xc,yc);
    if any(isnan([centroid_upper_x_mean,centroid_bottom_x_mean]))
        disp("At least one vertebrae wasn't detected.")
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
[mask,mask_ind] = find_roi(I,intersection,centroid_upper_x_mean,centroid_upper_y_mean,centroid_bottom_x_mean,centroid_bottom_y_mean);
figure;imshow(mask)
I(mask==0)=0;
figure;imshow(I)


[Ld] = apply_watershed_v2(array,mask,mask_ind);
[final_stack,area,fusion] = vypocet_plochy_v9(Ld)

