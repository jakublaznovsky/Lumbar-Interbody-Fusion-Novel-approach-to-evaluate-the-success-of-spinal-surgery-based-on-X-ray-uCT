function [I2,xc,yc,imgfill] = detect_concave_and_convex_points(imgbw)

% Image Binarization by otsu's method
% level = graythresh(I);
% imgbw =  im2bw(I,level); %level? 0.3
% I = histeq(I);
% imgbw =  im2bw(I,0.85);
figure
imshow(imgbw)
title('BW');
imgbw = imclose(imgbw, strel('diamond',5));                                  %10
% parameters for css method
C = 1.5;
T_angle = 162;
sig = 3;
H = 0.35;
L = 0;
Endpoint = 1;
Gap_size = 1;

k = 10; %kth adjucnet points to the corner point                             20

% preprocessing to smooth the objects boundarries
I2 = imgbw;
% se = strel('disk',4); 
% se1 = strel('disk',1);
% I2 = imerode(I2,se);
% I2 = imdilate(I2,se1);

I2 = imclose(I2,strel('disk',13));
I2 = bwareaopen(I2, 300000);
figure
imshow(I2)
title('after processing');

% remove vertebrae process???
% CC = bwconncomp(I2);
% if CC.NumObjects ==2
%     I2(CC.PixelIdxList{1,2})=0;
% end

% extracting the image edge by canny
BW = edge(I2,'canny',[L,H]);
%extract curves from binary edge map

[L,W]=size(BW);
BW1=zeros(L+2*Gap_size,W+2*Gap_size);
BW_edge=zeros(L,W);
BW1(Gap_size+1:Gap_size+L,Gap_size+1:Gap_size+W)=BW;
[r,c]=find(BW1==1);
cur_num=0;
while size(r,1)>0
    point=[r(1),c(1)];
    cur=point;
    BW1(point(1),point(2))=0;
    [I,J]=find(BW1(point(1)-Gap_size:point(1)+Gap_size,point(2)-Gap_size:point(2)+Gap_size)==1);
    while size(I,1)>0
        dist=(I-Gap_size-1).^2+(J-Gap_size-1).^2;
        [min_dist,index]=min(dist);
        point=point+[I(index),J(index)]-Gap_size-1;
        cur=[cur;point];
        BW1(point(1),point(2))=0;
        [I,J]=find(BW1(point(1)-Gap_size:point(1)+Gap_size,point(2)-Gap_size:point(2)+Gap_size)==1);
    end
    
    % Extract edge towards another direction
    point=[r(1),c(1)];
    BW1(point(1),point(2))=0;
    [I,J]=find(BW1(point(1)-Gap_size:point(1)+Gap_size,point(2)-Gap_size:point(2)+Gap_size)==1);
    while size(I,1)>0
        dist=(I-Gap_size-1).^2+(J-Gap_size-1).^2;
        [min_dist,index]=min(dist);
        point=point+[I(index),J(index)]-Gap_size-1;
        cur=[point;cur];
        BW1(point(1),point(2))=0;
        [I,J]=find(BW1(point(1)-Gap_size:point(1)+Gap_size,point(2)-Gap_size:point(2)+Gap_size)==1);
    end
    
    if size(cur,1)>(size(BW,1)+size(BW,2))/65
        cur_num=cur_num+1;
        curve{cur_num}=cur-Gap_size;
    end
    [r,c]=find(BW1==1);
    
end

for i=1:cur_num
    curve_start(i,:)=curve{i}(1,:);
    curve_end(i,:)=curve{i}(size(curve{i},1),:);
    if (curve_start(i,1)-curve_end(i,1))^2+...
            (curve_start(i,2)-curve_end(i,2))^2<=32
        curve_mode(i,:)='loop';
    else
        curve_mode(i,:)='line';
    end
    
    BW_edge(curve{i}(:,1)+(curve{i}(:,2)-1)*L)=1;
end

%get corners
corner_num=0;
cout=[];

GaussianDieOff = .0001;
pw = 1:30;
ssq = sig*sig;
width = max(find(exp(-(pw.*pw)/(2*ssq))>GaussianDieOff));
if isempty(width)
    width = 1;
end
t = (-width:width);
gau = exp(-(t.*t)/(2*ssq))/(2*pi*ssq);
gau=gau/sum(gau);

for i=1:cur_num; %%only loop
    x=curve{i}(:,1);
    y=curve{i}(:,2);
    W=width;
    L=size(x,1);
    if L>W
        
        % Calculate curvature
        if curve_mode(i,:)=='loop'
            x1=[x(L-W+1:L);x;x(1:W)];
            y1=[y(L-W+1:L);y;y(1:W)];
        else
            x1=[ones(W,1)*2*x(1)-x(W+1:-1:2);x;ones(W,1)*2*x(L)-x(L-1:-1:L-W)];
            y1=[ones(W,1)*2*y(1)-y(W+1:-1:2);y;ones(W,1)*2*y(L)-y(L-1:-1:L-W)];
        end
        
        xx=conv(x1,gau);
        xx=xx(W+1:L+3*W);
        yy=conv(y1,gau);
        yy=yy(W+1:L+3*W);
        Xu=[xx(2)-xx(1) ; (xx(3:L+2*W)-xx(1:L+2*W-2))/2 ; xx(L+2*W)-xx(L+2*W-1)];
        Yu=[yy(2)-yy(1) ; (yy(3:L+2*W)-yy(1:L+2*W-2))/2 ; yy(L+2*W)-yy(L+2*W-1)];
        Xuu=[Xu(2)-Xu(1) ; (Xu(3:L+2*W)-Xu(1:L+2*W-2))/2 ; Xu(L+2*W)-Xu(L+2*W-1)];
        Yuu=[Yu(2)-Yu(1) ; (Yu(3:L+2*W)-Yu(1:L+2*W-2))/2 ; Yu(L+2*W)-Yu(L+2*W-1)];
        K=abs((Xu.*Yuu-Xuu.*Yu)./((Xu.*Xu+Yu.*Yu).^1.5));
        K=ceil(K*100)/100;
        
        % Find curvature local maxima as corner candidates
        extremum=[];
        N=size(K,1);
        n=0;
        Search=1;
        
        for j=1:N-1
            if (K(j+1)-K(j))*Search>0
                n=n+1;
                extremum(n)=j;  % In extremum, odd points is minima and even points is maxima
                Search=-Search;
            end
        end
        if mod(size(extremum,2),2)==0
            n=n+1;
            extremum(n)=N;
        end
        
        n=size(extremum,2);
        flag=ones(size(extremum));
        
        % Compare with adaptive local threshold to remove round corners
        for j=2:2:n
            %I=find(K(extremum(j-1):extremum(j+1))==max(K(extremum(j-1):extremum(j+1))));
            %extremum(j)=extremum(j-1)+round(mean(I))-1; % Regard middle point of plateaus as maxima
            
            [x,index1]=min(K(extremum(j):-1:extremum(j-1)));
            [x,index2]=min(K(extremum(j):extremum(j+1)));
            ROS=K(extremum(j)-index1+1:extremum(j)+index2-1);
            K_thre(j)=C*mean(ROS);
            if K(extremum(j))<K_thre(j)
                flag(j)=0;
            end
        end
        extremum=extremum(2:2:n);
        flag=flag(2:2:n);
        extremum=extremum(find(flag==1));
        
        % Check corner angle to remove false corners due to boundary noise and trivial details
        flag=0;
        smoothed_curve=[xx,yy];
        while sum(flag==0)>0
            n=size(extremum,2);
            flag=ones(size(extremum));
            for j=1:n
                if j==1 & j==n
                    ang=mia_curve_tangent(smoothed_curve(1:L+2*W,:),extremum(j));
                elseif j==1
                    ang=mia_curve_tangent(smoothed_curve(1:extremum(j+1),:),extremum(j));
                elseif j==n
                    ang=mia_curve_tangent(smoothed_curve(extremum(j-1):L+2*W,:),extremum(j)-extremum(j-1)+1);
                else
                    ang=mia_curve_tangent(smoothed_curve(extremum(j-1):extremum(j+1),:),extremum(j)-extremum(j-1)+1);
                end
                if ang>T_angle & ang<(360-T_angle)
                    flag(j)=0;
                end
            end
            
            if size(extremum,2)==0
                extremum=[];
            else
                extremum=extremum(find(flag~=0));
            end
        end
        
        extremum=extremum-W;
        extremum=extremum(find(extremum>0 & extremum<=L));
        idx{1,i} =  extremum; % keep extermum for each curve
        n=size(extremum,2);
        for j=1:n
            corner_num=corner_num+1;
            cout(corner_num,:)=curve{i}(extremum(j),:);
            
        end
    end
end


% Add Endpoints
if Endpoint
    for i=1:cur_num
        if size(curve{i},1)>0 & curve_mode(i,:)=='line'
            
            % Start point compare with detected corners
            compare_corner=cout-ones(size(cout,1),1)*curve_start(i,:);
            compare_corner=compare_corner.^2;
            compare_corner=compare_corner(:,1)+compare_corner(:,2);
            if min(compare_corner)>25       % Add end points far from detected corners
                corner_num=corner_num+1;
                cout(corner_num,:)=curve_start(i,:);
            end
            
            % End point compare with detected corners
            compare_corner=cout-ones(size(cout,1),1)*curve_end(i,:);
            compare_corner=compare_corner.^2;
            compare_corner=compare_corner(:,1)+compare_corner(:,2);
            if min(compare_corner)>25
                corner_num=corner_num+1;
                cout(corner_num,:)=curve_end(i,:);
            end
        end
    end
end

% edit concave and convex points 
vis=1;
[idxpp,xc,yc,imgfill] = convex_or_concave(idx,curve,I2,k,vis);

end
