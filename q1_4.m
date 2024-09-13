close all
clear all
clc
format long

xi = @(x) (exp(x)-1);

D = 1;
D2 = 1;
alpha = [1.999 1.95 1 -1 -100]';
eps = 0.01;

xbar = zeros(size(alpha)); % Initialize xbar array
abar = zeros(size(alpha)); % Initialize Abar array
a = zeros(size(alpha)); % Initialize A array
pi_r = zeros(size(alpha)); % Initialize pi_r array
pi_r2 = zeros(size(alpha)); % Initialize pi_r2 array
varsplit = zeros(size(alpha)); % Initialize varsplit array
exkurtsplit = zeros(size(alpha)); % Initialize exkurtsplit array

for i = 1:numel(alpha)
    xbar(i) = 0.2*sqrt(0.02*(4-alpha(i))/(2-alpha(i)));
    abar(i) = (2-alpha(i))^2.5/(((4-alpha(i))^(1.5))*(0.02^1.5)*0.2);
    a(i) = ((2-alpha(i))^2*xbar(i)^alpha(i))/((4-alpha(i))*0.02);
    pi_r(i) = a(i)*xbar(i)^(-alpha(i))/-alpha(i);
    if pi_r(i) < 0
        pi_r(i) = Inf;
    end
    pi_r2(i) = abar(i)/(xbar(i)^(-1-alpha(i)))*(integral(@(x) (x).^(-1-alpha(i)), 0.5*xbar(i), xbar(i))); 
    var(i) = a(i)*(xbar(i)^(2-alpha(i)))/(2-alpha(i));
    exkurt(i) = a(i)*(xbar(i)^(4-alpha(i)))/(4-alpha(i));
    m2(i) = abar(i)/(xbar(i)^(-1-alpha(i)))*(integral(@(x) (x).^(1-alpha(i)), 0.5*xbar(i), 1*xbar(i))); 
    m4(i) = abar(i)/(xbar(i)^(-1-alpha(i)))*(integral(@(x) (x).^(3-alpha(i)), 0.5*xbar(i), 1*xbar(i))); 
    varsplit(i) = round(1 - (m2(i))/var(i),4);
    exkurtsplit(i) = round(1 - (m4(i))/exkurt(i),4);
end

alpha = [1.999 1.95 1 -1 -100]';
n = length(alpha);
p = zeros (1, n) ;
elogS = zeros(size(alpha)); % Initialize elogS array
vlogS = zeros(size(alpha)); % Initialize vlogS array
exkurtlogS = zeros(size(alpha)); % Initialize exkurtlogS array
min_ent_val = zeros(size(alpha)); % Initialize min_ent_val array

cgvf = zeros(1, n);
variation_future = zeros(1, n);
denom3 = zeros(1, n);
num4 = zeros(1, n);
beta2 = zeros(1, n);
hedging_coeff = zeros(1, n);
beta_4b = zeros(1, n);
stddev = zeros(1, n);
std_dev_hedging_error_to_maturity = zeros(1, n);

for i=1:n  
    kappa = @(v) TRUSTb(@(x) exp(v*x)-1,v,v^2,alpha(i));
    p(i) = fzero(@(p) kappa(1+p) - kappa(p), [-10 0]);
    elogS (i) = kappa (1); 
    vlogS(i) = sqrt(kappa(2)-2*kappa(1));
    num(i) = kappa(4)-4*kappa(3)+6*kappa(2)-4*kappa(1);
    denom(i) = (kappa(2)-2*kappa(1))^2;
    exkurtlogS(i) = num(i)/denom(i);
    
    %Q1
    ex4_1(i) = kappa(1)-0.1-0.5*0.25^2; % Exercise 4.1 the difference is due to the variation in jumps 
    
    %Q2a
    num1(i) = exp(0.5 * (kappa(1.5)-kappa(1)-kappa(0.5))) - 1;
    denom1(i) = sqrt(exp(0.5^2*(kappa(2)-2*kappa(1)))-1) * sqrt(exp(kappa(1) - 2 * kappa(0.5))-1);
    corr(i) = num1(i) / denom1(i);

    %Q2b
    esscher_eq = @(y) kappa(1 - y) - kappa(-y);
    y(i) = fzero(esscher_eq, 0); 
    esscher(i) = kappa(1 - y(i)) - kappa(-y(i));

    %Q2d
    varoptimal(i) = elogS(i)/vlogS(i)^2;
 
    %Q2c
    minent_eq = @(z) TRUSTb( @(x) (exp(x)-1).*exp(-z.*(exp(x)-1)),1,-2*z+1,alpha(i) );
    lambda(i) = fzero(minent_eq,0);
    minent(i) = TRUSTb( @(x) (exp(x)-1).*exp(-lambda(i).*(exp(x)-1)),1,-2*lambda(i)+1,alpha(i) );

    %Q2e
    slc(i) =  TRUSTb( @(x) -x.*(1-varoptimal(i).*(exp(x)-1)),-1,2*varoptimal(i),alpha(i) );

    %Q2f
    num2(i) = TRUSTb(@(x) x.*(exp(x)-1),0,2,alpha(i));
    denom2(i) = TRUSTb(@(x)(exp(x)-1).^2,0,2,alpha(i));
    beta(i) = num2(i)/denom2(i);
    
    %Q2g
    bvv(i) = TRUSTb( @(x) x.^2,0,2,alpha(i));
    e0(i) = sqrt(bvv(i) - beta(i)^2*vlogS(i)^2);

    %Q3
    vf(i) =  TRUSTb( @(x) (exp(x)-1).^2.*(1-varoptimal(i).*(exp(x)-1)),0,2,alpha(i) );

    num3(i) = TRUSTb(@(x) (exp(x)-1).^3,0,0,alpha(i));
    beta1(i) = num3(i)/denom2(i);

    bvv1(i) = TRUSTb( @(x) (exp(x)-1).^4,0,0,alpha(i));
    e1(i) = sqrt(bvv1(i) - beta1(i)^2*vlogS(i)^2);

    %Q4
    cgvf(i) = TRUSTb(@(x) x.^2 .* (1 - lambda(i).* (exp(x) - 1)), 0, 2, alpha(i));
    denom3(i) = TRUSTb(@(x) (exp(x) - 1).^2, 0, 2, alpha(i));
    num4(i) = TRUSTb(@(x) x.^2 .* (exp(x) - 1), 0, 0, alpha(i));
    beta2(i) = num4(i) / denom3(i);
    beta_comparison(i) = 1 - beta(i);
    stddev(i) = sqrt(bvv1(i) - beta2(i)^2 * vlogS(i));
end

ELogS = round(elogS,8);
SDLogS = round(vlogS,7);
ExKurtLogS = round(exkurtlogS,6);
Xbar = arrayfun(@(x) sprintf('%.1f%%', x * 100), xbar, 'UniformOutput', false);
A = arrayfun(@(x) sprintf('%.2e', x), a, 'UniformOutput', false) ;
Abar = arrayfun(@(x) sprintf('%.2e', x), abar, 'UniformOutput', false);
Pi_r = round(pi_r,1);
Pi_r2 = [round(pi_r2(1,1),7) round(pi_r2(2,1),3) round(pi_r2(3,1),1) round(pi_r2(4,1),1) round(pi_r2(5,1),1)]';
Var_split = arrayfun(@(x) sprintf('%.1f%%', x * 100), varsplit, 'UniformOutput', false);
Exkurt_split = arrayfun(@(x) sprintf('%.1f%%', x * 100), exkurtsplit, 'UniformOutput', false);

% Question 1 - Exercise 4.1
disp('Q1:');
arrayfun(@(x) sprintf('%.2e', x), ex4_1, 'UniformOutput', false)'

%Table 3.1 - Question 2
Table_1 = table(alpha, Xbar, A, Abar, Pi_r, Pi_r2,Var_split,Exkurt_split)

%Table 4.1 - Question 2
Table_2 = table(alpha, Xbar, Abar, ELogS, SDLogS,ExKurtLogS)

Corr = arrayfun(@(x) sprintf('%.5f%', x), corr, 'UniformOutput', false)';
Ess = arrayfun(@(x) sprintf('%.5f%', x), y, 'UniformOutput', false)';
Ent = arrayfun(@(x) sprintf('%.5f%', x), lambda, 'UniformOutput', false)';
Opt = arrayfun(@(x) sprintf('%.5f%', x), varoptimal, 'UniformOutput', false)';
column = {'alpha','Correlation', 'Esscher parameter', 'Minimal Entropy Parameter', 'Variance Optimal Parameter'};

% Table showing results for  Q2a)-Q2d)
Table_3  = table(alpha,Corr,Ess,Ent,Opt,'VariableNames',column)

SLC = arrayfun(@(x) sprintf('%.6f%', x), slc, 'UniformOutput', false)';
Beta = arrayfun(@(x) sprintf('%.3f%%', x * 100), beta, 'UniformOutput', false)';
Err = arrayfun(@(x) sprintf('%.1f bp%', x*10000), e0, 'UniformOutput', false)';
column1 = {'Alpha','Short log Contract Parameter', 'Hedging Coefficient','Std Dev of Hedging error to maturity'};

% Table showing results for  Q2e)-Q2g)
Table_4  = table(alpha,SLC,Beta,Err,'VariableNames',column1)

VF = arrayfun(@(x) sprintf('%.6f%', x), vf, 'UniformOutput', false)';
Beta1 = arrayfun(@(x) sprintf('%.3f bp%', x * 10000), beta1, 'UniformOutput', false)';
Err1 = arrayfun(@(x) sprintf('%.2f bp%', x*10000), e1, 'UniformOutput', false)';
column2 = {'alpha','Variance Future Parameter', 'Hedging Coefficient','Std Dev of Hedging error to maturity'};

% Table showing results for  Q2e)-Q2g)
Table_5  = table(alpha,VF,Beta1,Err1,'VariableNames',column2)

Beta2 = arrayfun(@(x) sprintf('%.5f%%', x*100), beta2, 'UniformOutput', false)';
Err2 = arrayfun(@(x) sprintf('%.10f%', x), stddev, 'UniformOutput', false)';
Err3 = arrayfun(@(x) sprintf('%.1f bp%', x*10000), stddev, 'UniformOutput', false)';
CGVF = arrayfun(@(x) sprintf('%.6f%', x), cgvf, 'UniformOutput', false)';
column3 = {'alpha','CG Variance Future', 'CG Hedging coefficient','CG Hedging error '};
column4 = {'alpha','Hedging error SLC', 'Hedging error VF','Hedging error CG'};

Table_6  = table(alpha,CGVF,Beta2,Err2,'VariableNames',column3)
Table_7  = table(alpha,Err,Err1,Err3,'VariableNames',column4)


function b = TRUSTb(xi,D,D2,alpha)
eps = 0.01;
xbar = 0.2*sqrt(0.02*(4-alpha)/(2-alpha));
Abar = (2-alpha)^2.5/(4-alpha)^1.5/0.02^1.5/0.2/2;
b = 0.1*D + 0.5*D2*(0.15^2+2*Abar*eps^(2-alpha)/(2-alpha)*xbar^(1+alpha))...
    + Abar*(integral(@(x) (xi(x)-D*x).*(abs(x/xbar)).^(-alpha-1),-xbar,-eps)...
    + integral(@(x) (xi(x)-D*x).*(abs(x/xbar)).^(-alpha-1),eps,xbar));
end