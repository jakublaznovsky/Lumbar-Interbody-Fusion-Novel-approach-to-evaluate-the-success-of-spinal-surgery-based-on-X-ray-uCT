function [idxpp,xc,yc,imgfill] = convex_or_concave(idx,curve,mask,k,vis)

%   Synopsis
%       [idxpp,xc,yc,imgfill] = convex_or_concave(idx,curve,mask,k,vis)
%   Description
%        Returns the convex and concave points in objects boundary.The corner 
%        point pi is qualified as concave if the line
%        connecting p i−k to p i+k does not reside inside 
%        contour segmentation and segment grouping
%   Inputs 
%         - idx      index of detected corener in the input image.    
%         - curve   objects boundareis
%         - mask    binar mask of image
%         - k       kth adjacent contour points
%         - vis     0 or 1, for visualizing the dtected concave points
%   Outputs
%        - idxpp    index of detected concave in the input image.
%        - xc       x-coordinate of detected corners
%        - yc       y-coordinate of detected corners 
%        - imgfill  preproccessed image
%   Authors
%          Sahar Zafari <sahar.zafari(at)lut(dot)fi>       MODIFIED
%
%   Changes
%       14/01/2016  First Edition


%preprocessing of image
    se = strel('disk',1); 
    mask = imdilate(mask,se);
%finding convex and concave corners
    [nmrows, nmcols] = size(mask);
    l = 0;
    for i = 1:size(idx,2)
        if length(idx{i})== 0
            idxpp{i}=[];
        end
         for j = 1:length(idx{i})
            idxi =idx{i}(j);
            x = curve{i}(idxi,2);
            y = curve{i}(idxi,1);
            xyp = curve{i}(max(1,mod(idxi+k,size(curve{i},1))),:);
            xym = curve{i}(max(1,mod(end+(idxi-k),size(curve{i},1))),:);
            [psubx,psuby]=mia_bresenham(xyp(2),xyp(1),xym(2),xym(1));
            pidx = sub2ind([nmrows, nmcols], psuby, psubx);
            th = 0.90; 
            ratfg = sum(mask(pidx))/length(pidx);
            res = ratfg > th;
          if res    % if res == 1, then it is a convex corner
              l = l+1;
              x_convex(l) = x;
              y_convex(l) = y;
              idxpp{i}(j) = idxi;
          else
              l = l+1;
              x_concave(l) = x;
              y_concave(l) = y;
              idxpp{i}(j) = idxi;
          end

        end
    end
    for i=1:length(idxpp)
      if  ~isempty(idxpp{i})
       idxpp{i}(idxpp{i}==0)=[];
      end
    end
%imfill the holes but leave the one from vertebrae process
    imgfill = imfill(mask,'holes');
    inv = ~mask;
    inv = imclearborder(inv);
    CC = bwconncomp(inv);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [biggest,idx] = max(numPixels);
    if 10000<biggest & biggest<26000   
        imgfill(CC.PixelIdxList{idx}) = 0;
    end

    x_concave = x_concave(x_concave ~= 0);
    y_concave = y_concave(y_concave ~= 0);
    x_convex = x_convex(x_convex ~= 0);
    y_convex = y_convex(y_convex ~= 0);
    
%remove detected concave corners inside the object
    for i=1:length(x_concave)  
        sub_matrix = imgfill(y_concave(i)-1:y_concave(i)+1,x_concave(i)-1:x_concave(i)+1);
        if sum(sum(sub_matrix)) == 9
            x_concave(i) = 0;
            y_concave(i) = 0;
        end
    end
%remove detected convex corners inside the object    
    for i=1:length(x_convex)  
        sub_matrix = imgfill(y_convex(i)-1:y_convex(i)+1,x_convex(i)-1:x_convex(i)+1);
        if sum(sum(sub_matrix)) == 9
            x_convex(i) = 0;
            y_convex(i) = 0;
        end
    end
    
%reduce the number of detected concaive points which are close
% together 15pixels
    for i=1:length(x_concave)
        for j=1:length(x_concave)
            if (i ~= j) && (sqrt((x_concave(i)-x_concave(j))^2 + (y_concave(i)-y_concave(j))^2) < 15) && (x_concave(i)~=0 && y_concave(i)~=0)
                x_concave(j) = 0;
                y_concave(j) = 0;
                
                  
            end
        end
    end
    
    for i=1:length(x_convex)
        for j=1:length(x_convex)
            if (i ~= j) && (sqrt((x_convex(i)-x_convex(j))^2 + (y_convex(i)-y_convex(j))^2) < 15) && (x_convex(i)~=0 && y_convex(i)~=0)
                x_convex(j) = 0;
                y_convex(j) = 0;
            end
        end
    end

    x_concave = x_concave(x_concave ~= 0);
    y_concave = y_concave(y_concave ~= 0);
    x_convex = x_convex(x_convex ~= 0);
    y_convex = y_convex(y_convex ~= 0);
%display the result
     if vis == 1
        figure;
        imshow(imgfill); hold on
        plot(x_concave,y_concave,'sg','markerfacecolor','g','LineWidth',5); 
        plot(x_convex,y_convex,'sr','markerfacecolor','r','LineWidth',5);
        hold on
        title('Detected Concave(green) and Convex(red) points')
     end

     xc = [x_concave,x_convex];
     yc = [y_concave,y_convex];

    end
