function [final_stack,area,fusion] = area_ratio(Ld)

%load('7A8_Ld_new_methodology.mat');
%Ld=flip(Ld,1);
 Ld_cropped=Ld;

% Ld_normal=Ld;
%% normalizace na 0=pozadí, 1,2=obratle (fce bwlabel?)
% Ld_normal=zeros(size(Ld_cropped));
% 
% for i=1:size(Ld_cropped,1)
%     for j=1:size(Ld_cropped,2)
%         for k=1:size(Ld_cropped,3)
%             if Ld_cropped(i,j,k)==1
%                 Ld_normal(i,j,k)=2;
%             elseif Ld_cropped(i,j,k)==2
%                 Ld_normal(i,j,k)=1;                
% %                 elseif Ld_cropped(i,j,k)==3
% %                 Ld_normal(i,j,k)=1;
%             end
%         end
%     end
% end

%% dataset cropping
% Ld_cropped=zeros(size(Ld_normal,1),size(Ld_normal,2));
% 
% for i=1:size(Ld_normal,3)
%   if length(unique(Ld_normal(:,:,i)))>=3 %zajištění že jsou v řezu oba obratle
%         Ld_cropped(:,:,end+1)=Ld_normal(:,:,i);        
%   end
% end
% Ld_cropped(:,:,1)=[]; 

%% fusion detection
p=1
for i=2:size(Ld_cropped,1)-1
    for j=2:size(Ld_cropped,2)-1
        for k=2:size(Ld_cropped,3)-1
            if sum(sum(Ld_cropped(i-1:i+1,j-1:j+1,k)==1))>=1
                obratel1=1;
            else
                obratel1=0;
            end
            if sum(sum(Ld_cropped(i-1:i+1,j-1:j+1,k)==2))>=1
                obratel2=1;
            else
                obratel2=0;
            end
            if Ld_cropped(i,j,k)==0&&obratel1==1&&obratel2==1
                Ld_cropped(i,j,k)=3;
            end
            
        end
    end
    p=p+1
end

       
        
%% replacing the connection point with zeros and vertebrae detection
final_stack=zeros(size(Ld_cropped,1),size(Ld_cropped,2));

obratle_bez_fuze=zeros(size(Ld_cropped));
fuze=zeros(size(Ld_cropped));
obratel1=zeros(size(Ld_cropped));
obratel2=zeros(size(Ld_cropped));

for i=1:size(obratle_bez_fuze,1)
    for j=1:size(obratle_bez_fuze,2)
        for k=1:size(obratle_bez_fuze,3)
            if Ld_cropped(i,j,k)~=3 %extraction of both vertebrae into a new variable - without fusion sites. The fusion sites themselves go into another variable
                obratle_bez_fuze(i,j,k)=Ld_cropped(i,j,k);
            else
                fuze(i,j,k)=3; %preservation of fusion sites for future use
            end
            if Ld_cropped(i,j,k)==1
                obratel1(i,j,k)=1;%obratle_bez_fuze(i,j,k); %vertebral extraction 1
            elseif Ld_cropped(i,j,k)==2
                obratel2(i,j,k)=1;%obratle_bez_fuze(i,j,k); %vertebral extraction 2
            end
        end
    end
end

obratle=obratel1+obratel2;
%% Vertebral surface detection

SE = strel('cube',60);
obr_closed=imclose(obratle_bez_fuze,SE);
m1=[0 0 0; 0 1 0; 0 0 0];
m2=[0 1 0; 1 1 1; 0 1 0];
mask=zeros(3,3,3);
mask(:,:,1)=m1;mask(:,:,2)=m2;mask(:,:,3)=m1;
obr_closed_boundary=obr_closed;

p=1;
for i=2:size(obr_closed,1)-1
    for j=2:size(obr_closed,2)-1
        for k=2:size(obr_closed,3)-1
            pom=obr_closed(i-1:i+1,j-1:j+1,k-1:k+1);
            if sum(unique(pom(mask==1)))==3
                obr_closed_boundary(i,j,k)=3; %marks adjecent px with 3
            end
        end
    end
    p=p+1
end

%remove false boundary pixels of dilated mask
obr_closed_boundary(obr_closed_boundary~=3)=0;
CC=bwconncomp(obr_closed_boundary);
        numPixels = cellfun(@numel,CC.PixelIdxList);
        [s,idx]=sort( numPixels, 'descend' );

for i=2:length(s) % first two biggest volumes belong to vertebraes
    idx2=idx(i);
    pixels=CC.PixelIdxList(idx2);
    obr_closed_boundary(pixels{1})=0;
end
       

%extraction of individual vertebraes
obratel1_hranice=obratel1;
obratel2_hranice=obratel2;


for k=2:size(obr_closed_boundary,3)-1

[rows, columns] = size(obr_closed_boundary(:,:,k));
upperRows=zeros(size(obr_closed_boundary,1),1);
for column = 1 : columns
    thisColumn = obr_closed_boundary(:,column,k);
    row = find(thisColumn==3, 1, 'first'); %indicates the column on which the boundary occurs
    if ~isempty(row)
        upperRows(column) = row; 
    end
end
first=find(upperRows,1,'first');%finding the column in which to start looking for the boundary
last=find(upperRows,1,'last');%finding the column in which to start looking for the border from the other side
coef=round((last-first)/10);
first=first+coef;
last=last-coef;

no=nnz(upperRows)-coef; %number of boundary pixels
% marking adjecent pixels with 2

pp=1;
firsta=first;
while pp==1
    obratel1_row=find(obratel1(:,firsta,k),1,'first'); %finding the line in which to start looking for the border
    if isempty(obratel1_row)
        firsta=firsta+1;
    else
        obr1_hranice = bwtraceboundary(obratel1(:,:,k),[obratel1_row firsta],'E',8,no,'clockwise'); %find boundary pixels
        if obr1_hranice(1,:)==obr1_hranice(end,:)
            firsta=firsta+1;
        else
            pp=0;
        end
    end
end
for i=1:length(obr1_hranice)
    obratel1_hranice(obr1_hranice(i,1),obr1_hranice(i,2),k)=2;
end

pp=1;
firstb=first;
while pp==1
    obratel2_row=find(obratel2(:,firstb,k),1,'last'); %finding a line in which to start looking for the border from the other side
    if isempty(obratel2_row)
        firstb=firstb+1;
    else
        obr2_hranice = bwtraceboundary(obratel2(:,:,k),[obratel2_row firstb],'E',8,no,'counterclockwise'); %find boundary pixels
        if obr2_hranice(1,:)==obr2_hranice(end,:)
            firstb=firstb+1;
        else
            pp=0;
        end
    end
end
for i=1:length(obr2_hranice)
    obratel2_hranice(obr2_hranice(i,1),obr2_hranice(i,2),k)=2;
end

pp=1;
lasta=last;
while pp==1
    obratel1_row2=find(obratel1(:,lasta,k),1,'first'); %finding the line in which to start looking for the boundary
    if isempty(obratel1_row2)
        lasta=lasta-1;
    else
        obr1_hranice2 = bwtraceboundary(obratel1(:,:,k),[obratel1_row2 lasta],'W',8,no,'counterclockwise'); %find boundary pixels from opposite direction
        if obr1_hranice2(1,:)==obr1_hranice2(end,:)
            lasta=lasta-1;
        else
            pp=0;
        end
    end
end
    for i=1:length(obr1_hranice2)
        obratel1_hranice(obr1_hranice2(i,1),obr1_hranice2(i,2),k)=2;
    end

pp=1;
lastb=last;
while pp==1
    obratel2_row2=find(obratel2(:,lastb,k),1,'last'); %finding the line in which to start looking for the border from the other side
    if isempty(obratel2_row2)
        lastb=lastb-1;
    else        
        obr2_hranice2 = bwtraceboundary(obratel2(:,:,k),[obratel2_row2 lastb],'W',8,no,'clockwise'); %find boundary pixels from opposite direction
        if obr2_hranice2(1,:)==obr2_hranice2(end,:)
            lastb=lastb-1;
        else
            pp=0;
        end
    end
end
for i=1:length(obr2_hranice2)
    obratel2_hranice(obr2_hranice2(i,1),obr2_hranice2(i,2),k)=2;
end

end


obratle_hranice=obratel1_hranice+obratel2_hranice; %marked borders on both vertebrae

%odstraneni hranicnich px ze samostatnych objektu (artefakty)
obratle_hranice_pom=obratle_hranice;
obratle_hranice_pom(obratle_hranice_pom~=2)=0;
CC=bwconncomp(obratle_hranice_pom);
        numPixels = cellfun(@numel,CC.PixelIdxList);
        [s,idx]=sort( numPixels, 'descend' );

for i=3:length(s) % first two biggest volumes belong to vertebraes
    idx2=idx(i);
    pixels=CC.PixelIdxList(idx2);
    obratle_hranice(pixels{1})=0;
end


%creation of final stack
final_stack=obratle_hranice + fuze;

%%
% plot(obr1_hranice2(:,2),obr1_hranice2(:,1),'r','LineWidth',2)
% plot(obr2_hranice2(:,2),obr2_hranice2(:,1),'r','LineWidth',2)

figure;imshow3Dfull(final_stack);

area=sum(sum(sum(final_stack==2)))
area_half=sum(sum(sum(final_stack==2)))/2
fusion=sum(sum(sum(final_stack==3)))



end









