%Code of retina layer segmentation based on graph theory using Dijkstra's
%algorirhm, coded by Yicheng Zhang
warning off;
close all;
clear;
clc;

origImg = imread('\exampleOCTimage0001.tif');

origImg = imresize(origImg,0.2);

origImg = im2double(origImg);

%filter image
 img = imfilter(origImg,fspecial('gaussian',[5 5],1),'replicate');
% img = imfilter(origImg, fspecial('average', [3,3]));
figure, imshow(img);


% % image flattening 
% [rows, cols] = size(img);
% 
% brightestRow = zeros(1, cols);
% 
% for i = 1 : cols
%     brightestRow(i) = max(img(:,i));
% end
% 
% brightestValue = max(brightestRow);
% 
% [pilotX, pilotY] = find(brightestRow == brightestValue);
% 
% meanValue = mean(brightestRow);
% 
% adjustBrightestRow = zeros(1, cols);
% 
% % for i = 1 : cols 
% %     if  abs(meanValue - brightestRow(i)) > 5
% %         adjustBrightestRow(i) = meanValue;
% %     else
% %         adjustBrightestRow(i) = brightestRow(i);
% %     end
% % end
% 
% shiftOffset = brightestRow - meanValue;
% shiftOffset = ceil(shiftOffset);

% img = imsharpen(img);



%pad image with vertical column on both sides
sizeOfOrigImg = size(img);

% paddedImg = zeros([sizeOfOrigImg(1) sizeOfOrigImg(2)+2]);
% paddedImg(:,2:1+sizeOfOrigImg(2)) = img;
numOfRows = sizeOfOrigImg(1);
numOfCols = sizeOfOrigImg(2);
rowStart = numOfRows * 0.4;
rowEnd = numOfRows *  0.8;

paddedImg = zeros([rowEnd - rowStart + 1, sizeOfOrigImg(2)+2]);
paddedImg(:, 2:1+sizeOfOrigImg(2)) = img(rowStart:rowEnd,:);

sizeOfPaddedImg = size(paddedImg);

%vertical gradient image
gradientImg = zeros(sizeOfPaddedImg);

for x = 1:sizeOfPaddedImg(2)
    gradientImg(:, x) = -1 * gradient(paddedImg(:, x), 2);
end

gradientImg = (gradientImg-min(gradientImg(:)))/(max(gradientImg(:))-min(gradientImg(:)));
figure, imshow(gradientImg);

%inverse of the vertical gradient image
invertGradientImg = gradientImg * -1 + 1; 
figure, imshow(invertGradientImg);


%minimum weight
minWeight = 10^-5;

%array to store gradient image edge weights
gradientEdgeWeightArray = nan([numel(paddedImg),8]);

%array to store inverse gradient image edge weigths
inverseGradientEdgeWeightArray = nan([numel(paddedImg),8]);

%array to store starting node positions
startNodePositionsArray = nan([numel(paddedImg),8]);

%array to store ending node positions
endNodePositionsArray = nan([numel(paddedImg),8]);

%offset array defined to limit the path search direction
searchDirectionOffsetArray = [1, 1, 1, 0, 0, -1, -1, -1;
    1, 0, -1, 1, -1, 1, 0, -1];
            
%fill in the arrays with edge weights
sizeOfGradEdgeWArray = size(gradientEdgeWeightArray);

for index = 1 : sizeOfGradEdgeWArray(1) * sizeOfGradEdgeWArray(2) 
    [x, y] = ind2sub(sizeOfGradEdgeWArray, index);    
    [xx,yy] = ind2sub(sizeOfPaddedImg, x);    
    X = xx + searchDirectionOffsetArray(1, y);
    Y = yy + searchDirectionOffsetArray(2, y);
    
    if X >= 1 && X <= sizeOfPaddedImg(1) && Y >= 1 && Y <= sizeOfPaddedImg(2)
         if Y == 1 || Y == sizeOfPaddedImg(2);
            gradientEdgeWeightArray(x, y) = minWeight;
            inverseGradientEdgeWeightArray(x, y) = minWeight;
         else
            gradientEdgeWeightArray(x, y) = 2 - gradientImg(xx, yy) - gradientImg(X, Y) + minWeight;
            inverseGradientEdgeWeightArray(x, y) = 2 - invertGradientImg(xx,yy) - invertGradientImg(X,Y) + minWeight;
         end
        startNodePositionsArray(x, y) = sub2ind(sizeOfPaddedImg, xx, yy);
        endNodePositionsArray(x, y) = sub2ind(sizeOfPaddedImg, X, Y);
    end
    
end

%assemble the adjacency matrix
assemblyPositionMat = ~isnan(gradientEdgeWeightArray) & ~isnan(startNodePositionsArray) & ~isnan(endNodePositionsArray) & ~isnan(inverseGradientEdgeWeightArray);
gradientEdgeWeightArray = gradientEdgeWeightArray(assemblyPositionMat);
inverseGradientEdgeWeightArray = inverseGradientEdgeWeightArray(assemblyPositionMat);
startNodePositionsArray = startNodePositionsArray(assemblyPositionMat);
endNodePositionsArray = endNodePositionsArray(assemblyPositionMat);

%sparse adjacency matrices of gradient images
gradienAdjMatrix = sparse(startNodePositionsArray, endNodePositionsArray, gradientEdgeWeightArray, numel(paddedImg), numel(paddedImg));

%sparse adjacency matrices of inverse gradient images
inverseGradientAdjMatrix = sparse(startNodePositionsArray, endNodePositionsArray, inverseGradientEdgeWeightArray, numel(paddedImg), numel(paddedImg));

%get shortest path going from light to dark
[~, gradientPath{1}] = graphshortestpath( gradienAdjMatrix, 1, numel(paddedImg) );
[gradientPathY, gradientPathX] = ind2sub(sizeOfPaddedImg,gradientPath{1});

%delete first and last few points that is by the image borders
gradientPathX = gradientPathX(gradient(gradientPathY)~=0);
gradientPathY = gradientPathY(gradient(gradientPathY)~=0);

%get shortest path going from dark to light
[~, inverseGradPath{1}] = graphshortestpath( inverseGradientAdjMatrix, 1, numel(paddedImg) );
[inverseGradPathY, inverseGradPathX] = ind2sub(sizeOfPaddedImg, inverseGradPath{1});

%delete first and last few points which are image borders
inverseGradPathX = inverseGradPathX(gradient(inverseGradPathY)~=0);
inverseGradPathY = inverseGradPathY(gradient(inverseGradPathY)~=0);


%show result
figure, imshow( origImg(rowStart:rowEnd,:));
axis image;
colormap('gray'); 
hold on;
plot(inverseGradPathX, inverseGradPathY, 'g-', 'linewidth', 1); 
plot(gradientPathX, gradientPathY, 'r-', 'linewidth', 1);











