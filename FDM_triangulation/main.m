  	
clc;
clear;
fprintf('\nFDM for triangulation tested on the dataset Alcatraz Water Tower\n\n');

fprintf('Loading dataset ...\n\n');

data_path = fullfile('../','dataset','vis(2).mat');
load(data_path);

M = size(PArray,1)/3;
N = size(Img,2);

meantri = mean(sum(isfinite(Img)))/2;
fprintf('# scene points: %d\n', N);  
fprintf('# views: %d\n', M);
fprintf('# average L: %.2f\n', meantri);
fprintf('\nFDM running ... \n');
fprintf('(expect to finish in about 2 minutes)\n');
fprintf('\nProgress: \n');

X = zeros(3,N);
timedes = 0;
blksize = 1000;
numblk = ceil(N/blksize);
%for blk = 1:3
for blk = 1:numblk

    pts = blksize*(blk-1)+1 : blksize*blk;
    pts(pts>N) = [];

    tic
    [ X_2v ] = f_2views_fea( PArray, Img(:,pts) );
    inittime = toc;

    tic
    %[ X_blk, mres_blk, nitr_blk, ~] = f_tri_L2_FDM(PArray, Img(:,pts), X_2v, 1e-2);
    [ X_blk, mres_blk, nitr_blk, ~] = f_tri_L2_FDM(PArray, Img(:,pts), X_2v, 1e-2);
    timedes = timedes + toc;

    X(:,pts) = X_blk;
    fprintf('          %.2f%%,\n', 100*blk/numblk);
end

fprintf('\nDone\n');

plot3(X(1,:),X(2,:),X(3,:),'b.', 'markersize', 2);
fprintf('Total time: %.2fs.\n', timedes);
fprintf('Average time: %.4fs / triangulation instance.\n', timedes/N);
