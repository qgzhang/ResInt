
clc;
clear;
fprintf('\nRes-Int for KRot tested on the dataset Alcatraz Water Tower\n\n');

try
    fprintf('Checking mex ...\n');
    c_tri_single_L2(1, 1, 1, 1, 1, 1, 1,1);
    fprintf('Mex check passed.\n');
catch exception
    fprintf('Compiling source code...\n');
    mex c_tri_single_L2.c
end


fprintf('\nLoading dataset ...\n\n');

data_path = fullfile('../','dataset','vis(2).mat');
load(data_path);

M = size(PArray,1)/3;
N = size(Img,2);
meantri = mean(sum(isfinite(Img)))/2;
meantrans = mean(sum(isfinite(Img),2))/2;
obs = sum(sum(isfinite(Img),2))/2;
fprintf('# scene points M: %d\n', N);  
fprintf('# images L:       %d\n', M);
fprintf('# observations:   %d\n', obs);
fprintf('  average L:      %.2f\n', meantri);
fprintf('  average M:      %.2f\n\n', meantrans);



%% Initialisation
initpath = fullfile('../','dataset','X0T0V(2).mat');
if exist(initpath,'file')
    load(initpath);
    u = cell(1,M);                   %p = cell(1,M);
    for i = 1:M
        u{i} = Img(2*(i-1)+1:2*i,:); %p{i} = P(3*(i-1)+1:3*i,:);
    end
else
    u = cell(1,M);                   %p = cell(1,M);
    for i = 1:M
        u{i} = Img(2*(i-1)+1:2*i,:); %p{i} = P(3*(i-1)+1:3*i,:);
    end
    tic
    ub = 100;
    lb = 0;
    [Uf,Pf,dual,s] = krot_feas(u,A,(ub-lb)/2,1,1000); 
    X0 = Uf(1:3,:);
    T0 = Pf(:,4);    T0 = reshape(T0, 3, M);
    init_time = toc;
    save(initpath, 'X0','T0','init_time');
end
Imgv = cell(1,N);
for i = 1:N
    Imgv{i} = Img(:,i); 
end
            
clear PArray;
clear XArray;
clear Img;
        
% % Res-Int
% fprintf('\nRes-Int(sequential) running ... \n');
% fprintf('(expect to finish in about 7 minutes)\n\n\n');
% parallel = 0;
% tic
%     [X, T, mres, nitr] = f_krot_L2(R, Imgv, A, u, 1e-3, X0, T0, 1e-2,parallel,1);
% time_seq = toc;
% fprintf('----------------------Res_Int(sequential) Done--------------------------\n');
% fprintf('Total runtime: %12.4fs', time_seq);
% fprintf('\n\n');



try
    parpool();
catch
end

fprintf('\nRes-Int(parallel) running ... \n');
fprintf('(expect to finish in about 3 minutes (depending on available number of workers))\n');


parallel = 1;
tic
    [X, T, mres, nitr] = f_krot_L2(R, Imgv, A, u, 1e-3, X0, T0, 1e-2,parallel,1);
time_par = toc;
fprintf('----------------------Res_Int(parallel) Done--------------------------\n');
fprintf('Total runtime: %12.4fs', time_par);
fprintf('\n\n');
try
    poolobj = gcp('nocreate');
    delete(poolobj);
catch
end        
        
%% plot
K_path = fullfile('../','dataset','KK.mat');
load(K_path);
T_recon = zeros(3,M);
X(:,X(3,:)>400)=[];
PArray_recon = [R, T(:)];
for i = 1:M
    P = KK \ PArray_recon(3*(i-1)+1:3*i,:);
    T_recon(:,i) = - P(:,1:3) \ P(:,4);
end
figure; 
plot3(X(1,:), X(2,:), X(3,:), 'b.', 'markersize', 3); hold on;
plot3(T_recon(1,:), T_recon(2,:), T_recon(3,:), 'r+', 'markersize', 6); 
axis equal;
axis off;
grid off;
view([ 90.5, -66.0]);


