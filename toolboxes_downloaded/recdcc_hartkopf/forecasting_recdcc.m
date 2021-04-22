function forecasting_recdcc

% load('simdata_RecDCC.mat'); 
load('../data/JP60FFXdata_9factor_p60.mat');

[T,~,p] = dim{:};

opt = optimoptions('fmincon', 'GradObj', 'off', 'Display', 'none', 'TolX', 1e-5, 'TolFun', 1e-5);

startDate = 998;
estimFreq = 10;


M_t = zeros(p,p,T);
L_t = zeros(p,p,T);
thetav_t = zeros(2,p,T);
thetav_opt = zeros(2,p);
thetac_t = zeros(2,1,T);
H = zeros(p,startDate-1);
H_ = zeros(p,startDate);
S_t = zeros(p,p,T);


load('results_JP60FFX_p60_outsample_full_h1.mat', 'thetav_t', 'thetac_t');

for tt = startDate:T
    disp(tt);
    thetav_0 = thetav_t(:,:,tt);
    thetac_0 = thetac_t(:,:,tt);

    if ~mod((tt-startDate),estimFreq) || (tt-startDate)==0
        %% Step 1
        [M, L, cstar] = step1(cr(:,:,tt-startDate+1:tt-1));

        %% Step 2
        disp('Estimation : Idiosyncratics');
        for ii = 1:p
            thetav_opt(:,ii) = fmincon(@step2, thetav_0(:,ii), [1 1],[1],[],[],[0; 0],[],[], opt, cstar(ii,ii,:), startDate-1);
            [~, H(ii,:)] = step2(thetav_opt(:,ii), cstar(ii,ii,:), startDate-1);
        end
        thetav_0 = thetav_opt;

        %% Step 3
        disp('Estimation : Correlations');
        thetac_opt = fmincon(@step3, thetac_0, [1 1],[1],[],[],[0; 0],[],[], opt, cstar, H, p, startDate-1);
        thetac_0 = thetac_opt;
%         [~, R] = step3(thetac_opt, cstar, H, p, startDate-1);
    end
    
    cstar = calc_cstar(cr(:,:,tt-startDate+1:tt-1), L);
    for ii = 1:p
        H_(ii,:) = calc_h(thetav_opt(:,ii), cstar(ii,ii,:), startDate);
    end
    R_ = calc_R(thetac_opt, cstar, H_, p, startDate);

    %% Calculate S
    S = calc_S(H_, R_, L);

    S_t(:,:,tt) = S(:,:,end);

    M_t(:,:,tt) = M;
    L_t(:,:,tt) = L;
    thetav_t(:,:,tt) = thetav_opt;
    thetac_t(:,:,tt) = thetac_opt;

end

save results_JP60FFX_p60_outsample_full_h1 S_t thetac_t thetav_t M_t L_t cr

end % end of main


function [M, L, cstar] = step1(c)
%%
M = mean(c,3);
L = chol(M)';
cstar = zeros(size(c));
for tt = 1:size(c,3)
    cstar(:,:,tt) = L\c(:,:,tt)/L';
end

end


function [llf, h] = step2(par, cstar, T)
%%
gam = par(1);
del = par(2);
cons = 1 - gam - del;


h = zeros(1,T);
h(1) = cstar(:,:,1);
llf = 0;
for tt = 2:T
    h(tt) = cons + gam*cstar(:,:,tt-1) + del*h(tt-1);
    llf = llf - 0.5*log(h(tt)) - 0.5*cstar(:,:,tt)/h(tt);
end

llf = llf/(T-1);
llf = - llf;

end


function [llf, R] = step3(par, cstar, H, p, T)
%%
al = par(1);
be = par(2);
A = sqrt(al)*eye(p);
B = sqrt(be)*eye(p);
cons = eye(p) - A*A' - B*B';

llf = 10e8;
if cons(1,1) <= 0
    return
end


Ip = eye(p);
llf = 0;
Q = eye(p);
R = zeros(size(cstar));
R(:,:,1) = eye(p);
cq = eye(p);
for tt = 2:T
    D = diag(sqrt(H(:,tt)));
    DcsD = D\cstar(:,:,tt)/D;
    
    Q = cons + A*cq*A' + B*Q*B';
    diagQ = diag(sqrt(diag(Q)));
    cq = diagQ\DcsD/diagQ;
    R(:,:,tt) = diagQ\Q/diagQ;
    
    llf = llf - 0.5*log(det(R(:,:,tt))) - 0.5*trace((inv(R(:,:,tt))-Ip)*DcsD);
end

llf = llf/(T-1);
llf = - llf;

end


function cstar = calc_cstar(c, L)
%%
cstar = zeros(size(c));
for tt = 1:size(c,3)
    cstar(:,:,tt) = L\c(:,:,tt)/L';
end

end


function h = calc_h(par, cstar, T)
%%
gam = par(1);
del = par(2);
cons = 1 - gam - del;


h = zeros(1,T);
h(1) = cstar(:,:,1);
for tt = 2:T
    h(tt) = cons + gam*cstar(:,:,tt-1) + del*h(tt-1);
end

end


function R = calc_R(par, cstar, H, p, T)
%%
al = par(1);
be = par(2);
A = sqrt(al)*eye(p);
B = sqrt(be)*eye(p);
cons = eye(p) - A*A' - B*B';


Q = eye(p);
R = zeros(p,p,T);
R(:,:,1) = eye(p);
cq = eye(p);
DcsD = cstar(:,:,1);
for tt = 2:T
    Q = cons + A*cq*A' + B*Q*B';
    diagQ = diag(sqrt(diag(Q)));
    cq = diagQ\DcsD/diagQ;    
    R(:,:,tt) = diagQ\Q/diagQ;

    if tt < T
        D = diag(sqrt(H(:,tt)));
        DcsD = D\cstar(:,:,tt)/D;
    end
end

end


function S = calc_S(Hii, R, L)
%%
S = zeros(size(R));
T = size(R,3);

for tt = 1:T
    D = diag(sqrt(Hii(:,tt)));
    H = D*R(:,:,tt)*D;
    S(:,:,tt) = L*H*L';
end

end
%% end of file