function [Ld] = apply_watershed_v1(array,mask,mask_ind)
%%   Synopsis
%       [area,array_after_watershed] = apply_watershed(array_bin)
%   Description
%        After appropriate preproccesing of volume, it returns area of joint is micromillimeters besides vertebraes and
%        array with drawn joint.
%   Inputs 
%         - binar                              volume after binarization
%         - mask                               logical mask for ROI
%         - mask_ind                           indices of mask
%   Outputs
%        Area in micromillimeters besides vertebraes. Array_after_watershed for displaying result.

%% apply mask on binar

for i=1:size(array,3)
    pom1 = array(:,:,i);
    pom1(mask == 0) = 0;
    array_mask(:,:,i) = pom1(min(mask_ind(:,2)):max(mask_ind(:,2)),min(mask_ind(:,1)):max(mask_ind(:,1)));
end


%% Binarization
p=imadjustn(array_mask);
thresh = graythresh(p);
binar = imbinarize(p,thresh);

%% preproces binary array

pom1=binar;
pom2=bwareaopen(pom1,200);
pom3=imgradient3(pom2,'central');
pom4=pom3~=0;
pom5 = imfill(pom4,'holes');

pom6=imgradient3(pom5,'central');
pom7=pom6~=0;
pom8 = imfill(pom7,'holes');

pom9=pom1+pom8;
pom10=pom9~=0;
pom11 = imfill(pom10,'holes');

pom12=bwareaopen(pom11,2000);
%pom13=imclose(pom12,strel('cube',3));
pom13=pom12;
pom14=imfill(pom13,'holes');

pom15=pom14;
pom15(1,:,:)=1;pom15(end,:,:)=1;
pom15(:,:,1)=1;pom15(:,:,end)=1;
pom16=imfill(pom15,'holes');

        
preprocessed=pom16;

preprocessed_distance=imclose(preprocessed,strel('cube',5));      
preprocessed_distance=imfill(preprocessed_distance,'holes');

disp('preproces binary array done')





volume = sum(sum(sum(preprocessed)));
local_minima = 140;   
con = 1;
iter = 0;
while con == 1;
    iter = iter + 1
    D = -bwdist(~preprocessed_distance);
    watershed_mask1 = imextendedmin(D,local_minima);  
    D = imimposemin(D,watershed_mask1); %adjustment of the distance map
    D(~preprocessed_distance) = Inf; %adjustment of the distance map
    Ld = watershed(D); %Watershed separation
    Ld(~preprocessed) = 0; %Set to zero, all 
    unique(Ld)
    regions = regionprops3(Ld);
    maximum = max(regions.Volume);
    if maximum > volume/2 + 24000000   %if volume separeted approximately into half
        local_minima = local_minima - 20;%30
    else
        con = 0;
    end
end
%% automated watershed


%reduction of grafts - parts smaller than selected threshold
regions = regionprops3(Ld);
for i=1:size(regions,1)
    if regions.Volume(i) < 4000000
        Ld(Ld==i)=0;
    end
end
unique_values = unique(Ld);


for i=1:length(unique_values)
    Ld(Ld == unique_values(i)) = i-1;
end
unique_values = unique(Ld);

disp('automated watershed done')

%% remove processes
Ld1=Ld;

Ld(1,:,:)=[];Ld(end,:,:)=[];
Ld(:,:,1)=[];Ld(:,:,end)=[];

CC = bwconncomp(Ld);
        numPixels = cellfun(@numel,CC.PixelIdxList);
        [s,idx]=sort( numPixels, 'descend' );

for i=3:length(s) % first two biggest volumes belong to vertebraes
    idx2=idx(i);
    pixels=CC.PixelIdxList(idx2);
    Ld(pixels{1})=0;
end
       













