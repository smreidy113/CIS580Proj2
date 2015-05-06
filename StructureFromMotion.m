%function StructureFromMotion

close all;

% K = [568.996140852 0 643.21055941;
%      0 568.988362396 477.982801038;
%      0 0 1];
% nImages = 6;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% % Load images
% for iImage = 1 : nImages
%     str = sprintf('image%07d.bmp', iImage);
%     im{iImage} = imread(str);
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% % Load matching
% Mx = []; My = []; M = [];
% for iImage = 1 : nImages-1;
%     str = sprintf('matching%d.txt', iImage);
%     [mx, my, m] = LoadMatching(str, iImage, nImages);
%     Mx = [Mx;mx];
%     My = [My;my];
%     M = [M;m];
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% % Initialize 3D points, reconstruction index, and visibility matrix
% X3D = zeros(size(M,1), 3);
% ReconX = zeros(size(M,1),1);
% V = zeros(size(M,1), nImages);
% 
% c = 1;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% % Exclude outliers using F matrix
% for iImage1 = 1 : nImages-1
%     for iImage2 = iImage1+1 : nImages
%         idx1 = find(M(:,iImage1)==1);
%         idx2 = find(M(:,iImage2)==1);
%         idx = intersect(idx1, idx2);
%         
%         subplot(4,3,c);
%         hold on;
%         str = sprintf('image%07d.bmp', iImage1);
%         C = imread(str);
%         l = floor(length(C(1,:))/3);
%         image(0, 0, C);
%         str = sprintf('image%07d.bmp', iImage2);
%         C = imread(str);
%         image(l, 0, C);
%         set(gca,'YDir','reverse');
%         xlim([0 l*2])
%         ylim([0 floor(length(C(:,1)))])
%         
%         x1 = [Mx(idx,iImage1) My(idx,iImage1)];
%         x2 = [Mx(idx,iImage2) My(idx,iImage2)];
%         if size(x1,1) < 8
%             continue;
%         end
%         
%         prevSize = length(x1);
%         
%         for i=1:length(x1)
%             plot([x1(i,1) x2(i,1)+l], [x1(i,2) x2(i,2)], 'r--');
%         end
%         
%         [x1, x2, inlier] = GetInliersRANSAC(x1, x2);
%         M(idx(~inlier),iImage1) = 0;
%         
%         perc = length(x1)/prevSize;
%         
%         text(0,floor(length(C(:,1))),sprintf('Filtered: %3f of original', perc));
%         text(l/2,-30,sprintf('Image %d', iImage1));
%         text(3*l/2,-30,sprintf('Image %d', iImage2));
%         
%         for i=1:length(x1)
%             plot([x1(i,1) x2(i,1)+l], [x1(i,2) x2(i,2)], 'b-');
%         end
%         
%         hold off;
%         
%         c = c + 1;
%         
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Set initial two frames
initialframe1 = 2;
initialframe2 = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Get point index for two frames
idx1 = find(M(:,initialframe1)==1);
idx2 = find(M(:,initialframe2)==1);
idx = intersect(idx1, idx2);

x1 = [Mx(idx,initialframe1) My(idx,initialframe1)];
x2 = [Mx(idx,initialframe2) My(idx,initialframe2)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Get fundamental matrix and essential mtraix
F = EstimateFundamentalMatrix(x1, x2)
E = EssentialMatrixFromFundamentalMatrix(F,K)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Camera pose estimation
[Cset Rset] = ExtractCameraPose(E);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Triangulation and pose disambiguation
for i = 1 : 4
    Xset{i} = LinearTriangulation(K, zeros(3,1), eye(3), Cset{i}, Rset{i}, x1, x2);
    
end
[C R X] = DisambiguateCameraPose(Cset, Rset, Xset);

idx2rem = find(X(:,3) < 0 | abs(X(:,1)) > 100 | abs(X(:,2)) > 100 | abs(X(:,3)) > 100);

X(idx2rem,:) = [];
idx(idx2rem,:) = [];
x1(idx2rem,:) = [];
x2(idx2rem,:) = [];

figure;
hold on;
plot3(X(:,1),X(:,3),X(:,2),'b.');
Xpose1 = R*(C' + [0 0 20])';
plot3([C(1) Xpose1(1)], [C(3) Xpose1(3)], [C(2) Xpose1(2)], 'b-');
plot3(C(1), C(2), C(3), 'bo');
xlabel('X');
ylabel('Z');
zlabel('Y');
grid on;
view(45, 45);
%hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Nonlinear triangulation
disp('Nonlinear triangulation');
X = NonlinearTriangulation(K, zeros(3,1), eye(3), C, R, x1, x2, X);

idx2rem = find(X(:,3) < 0 | abs(X(:,1)) > 100 | abs(X(:,2)) > 100 | abs(X(:,3)) > 100);

X(idx2rem,:) = [];
idx(idx2rem,:) = [];
x1(idx2rem,:) = [];
x2(idx2rem,:) = [];

plot3(X(:,1),X(:,3),X(:,2),'r.');

hold off;

figure;

P = K*R*[eye(3) -C];

for i=1:length(X(:,1))
    
    Xhom = [X(i,:) 1]';
    
    x1rep(i,:) = [P(1,:)*Xhom/(P(3,:)*Xhom),P(2,:)*Xhom/(P(3,:)*Xhom)];
    
    if norm(x1rep(i,:) - x2(i,:)) > 5
        to_rem = [to_rem i];
    else
        %to_rem(i) = 0;
    end
end

to_rem(to_rem == 0) = [];

X(to_rem',:) = [];
idx(to_rem',:) = [];
x1(to_rem',:) = [];
x2(to_rem',:) = [];

hold on;
%plot(x1(:,1),x1(:,2),'b.');
plot(x2(:,1),x2(:,2),'g.');
plot(x1rep(:,1),x1rep(:,2),'r.');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Set reconstructed frame
r_idx = [initialframe1, initialframe2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Set camera pose
Cr_set{1} = zeros(3,1);
Rr_set{1} = eye(3,3);
Cr_set{2} = C;
Rr_set{2} = R;

% return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Set points and visibility matrix
X3D(idx,:) = X;
ReconX(idx) = 1;
V(idx, initialframe1) = 1;
V(idx, initialframe2) = 1;

Cr_set = {};
Rr_set = {};
r_idx = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Add images
for iImage = 1 : nImages
    if ~isempty(find(r_idx==iImage))
        continue;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
    % Get 2D-3D correspondences
    idx1 = find(ReconX==1);
    idx2 = find(M(:,iImage)==1);
    idx = intersect(idx1, idx2);
    if length(idx) < 6
        continue;
    end
    
    X = X3D(idx,:);
    x = [Mx(idx,iImage) My(idx,iImage)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
    % Run PnP
    [C, R] = PnPRANSAC(X, x, K);
C
R
    disp('Nonlinear PnP');
    [C, R] = NonlinearPnP(X, x, K, C, R);
C
R
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
    % Set camera poses and reconstructed frame index
    Cr_set{end+1} = C;
    Rr_set{end+1} = R;
    r_idx(end+1) = iImage;
    V(idx, iImage) = 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
    % Triangulation for additional points
    disp('Adding more points');
    for iImage1 = 1 : length(r_idx)-1
        idx1 = find(ReconX~=1);
        idx2 = find(M(:,r_idx(iImage1))==1);
        idx3 = find(M(:,iImage)==1);
        idx = intersect(idx1, idx2);
        idx = intersect(idx, idx3);
        x1 = [Mx(idx,r_idx(iImage1)) My(idx,r_idx(iImage1))];
        x2 = [Mx(idx,iImage) My(idx,iImage)];
        X = LinearTriangulation(K, Cr_set{iImage1}, Rr_set{iImage1}, C, R, x1, x2);
        X = NonlinearTriangulation(K, Cr_set{iImage1}, Rr_set{iImage1}, C, R, x1, x2, X);
        
        X3D(idx,:) = X;
        ReconX(idx) = 1;

        V(idx,r_idx(iImage1)) = 1;
        V(idx,iImage) = 1;

    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
    % Set visibiltiy and measurements for bundle adjustment
    V_bundle = V(:,r_idx);
    Mx_bundle = Mx(:,r_idx);
    My_bundle = My(:,r_idx);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
    % Run bundle adjustment
    disp('Bundle adjustment');
    [Cset Rset, X] = BundleAdjustment(K, Cr_set, Rr_set, X3D, ReconX, V_bundle, Mx_bundle, My_bundle);
end