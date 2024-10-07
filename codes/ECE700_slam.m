
% indexing help:
% pos = 2*T*(M+1)*(i-1) where i stands for robot i
% x((pos+1):(pos+2*T*M)) are the states x_i for robot i
% x((pos+2*T*M+1):(pos+2*T*(M+1))) are the states y_i for robot i
% x((pos+1):(pos+2*T*(M+1))) are the states x_i and y_i for robot i

T = 5;
z = zeros(2*T,5);
m = zeros(2*T,6);

M = 6;
n = 5;

z(:,1) = [4;
          1;
          3.9;
          1.1;
          3.8;
          1.2;
          3.7;
          1.3;
          3.6;
          1.4];

z(:,2) = [6;
          2.5;
          6.1;
          2.6;
          6.2;
          2.7;
          6.3;
          2.8;
          6.4;
          2.9];

z(:,3) = [2;
          4.5;
          2.2;
          4.6;
          2.4;
          4.7;
          2.5;
          4.9;
          2.6;
          5.1];

z(:,4) = [4;
          8.5;
          4.1;
          8.3;
          4.2;
          8.1;
          4.3;
          8.0;
          4.4;
          7.9];

z(:,5) = [7;
          7.5;
          7.1;
          7.4;
          7.2;
          7.3;
          7.3;
          7.1;
          7.4;
          7.0];

m(:,1) = [2;
          2.5;
          2;
          2.5;
          2;
          2.5;
          2;
          2.5;
          2;
          2.5];

m(:,2) = [7.5;
          2.5;
          7.6;
          2.5;
          7.4;
          2.5;
          7.6;
          2.5;
          7.4;
          2.5];

m(:,3) = [2.5;
          6;
          2.5;
          6.2;
          2.5;
          5.8;
          2.5;
          6.1;
          2.5;
          5.9];

m(:,4) = [7.5;
          6;
          7.5;
          6;
          7.5;
          6;
          7.5;
          6;
          7.5;
          6];

m(:,5) = [5;
          3;
          5.1;
          2.8;
          5.1;
          2.7;
          5.2;
          2.6;
          5.4;
          2.4];

m(:,6) = [6;
          9;
          5.8;
          8.7;
          5.6;
          8.4;
          5.2;
          8.1;
          5.0;
          7.8];

ztrue = z;
mtrue = m;

N = cell(1,5);

N{1} = [1 2 3 5];
N{2} = [2 4 5];
N{3} = [1 3 6];
N{4} = [3 5];
N{5} = [3 4 5];

deltaZ = [ -0.0986    0.1158    0.2441    0.0990    0.0984;
           0.0981    0.0922    0.0902   -0.1509   -0.1336;
           -0.0817    0.1180    0.2244    0.0872    0.0812;
           0.1033    0.1178    0.1164   -0.1823   -0.0925;
           -0.0456    0.0344    0.0685    0.0981    0.0702;
           0.1094    0.0602    0.2119   -0.1607   -0.1708;
           -0.0459    0.0568    0.0774    0.0705    0.0808;
           0.0783    0.1121    0.2455   -0.0816   -0.0457];

deltaM{1} = [-1.9520    3.4745   -1.4476    1.0278;
             1.5869    1.4339    5.0330    1.9722;
             -1.9215    3.6682   -1.2746    1.1552;
             1.3186    1.4159    5.1532    1.6795;
             -1.7917    3.6069   -1.2422    1.2920;
             1.3188    1.2645    4.6026    1.5205;
             -1.7113    3.9389   -1.2644    1.4524;
             1.1426    1.2311    4.7814    1.3159;
             -1.4988    3.8324   -1.1379    1.8039;
             0.9820    1.0787    4.4718    1.0662];

deltaM{2} = [1.4893    1.5524   -1.0206;
             -0.0067    3.4887    0.4816;
             1.4414    1.3919   -1.0680;
             -0.1693    3.4345    0.2390;
             1.2155    1.3278   -1.0780;
             -0.2125    3.2440   -0.0045;
             1.3252    1.1234   -1.0489;
             -0.3446    3.1451   -0.2437;
             1.0954    1.0292   -0.9793;
             -0.3939    3.1030   -0.4826];

deltaM{3} = [0.0175    0.4851    3.9215;
             -2.0365    1.3384    4.4761;
             -0.1837    0.2457    3.5331;
             -2.1257    1.5287    4.1015;
             -0.4448    0.0493    3.2427;
             -2.2602    1.0893    3.7202;
             -0.4481   -0.0163    2.6650;
             -2.4423    1.2972    3.1185;
             -0.6086   -0.1286    2.4730;
             -2.6604    0.7875    2.8025];

deltaM{4} = [-1.4940    1.0048;
             -2.5495   -5.4752;
             -1.5401    1.0541;
             -2.1296   -5.4515;
             -1.7235    0.8716;
             -2.2557   -5.3595;
             -1.8693    0.9087;
             -1.9978   -5.4253;
             -1.8790    0.9403;
             -1.9800   -5.4677];

deltaM{5} = [-4.5177    0.4532   -2.0389;
             -1.4977   -1.5635   -4.4842;
             -4.6396    0.4249   -1.9297;
             -1.2775   -1.2605   -4.5799;
             -4.6914    0.3364   -2.0535;
             -1.5031   -1.3387   -4.6803;
             -4.7400    0.2418   -2.0669;
             -0.9599   -1.1564   -4.3931;
             -4.8473    0.0288   -1.9729;
             -1.1374   -0.9641   -4.6770];

V = diag(-ones(T,1)) + diag(ones(T-1,1),1);
V = V(1:end-1,:);
V = kron(V,eye(2));


A = cell(1,5);

for i=1:5
    A{i} = zeros(2*(6*T+T));
    A{i}((2*T*6+1):end,(2*T*6+1):end) = V'*V + length(N{i})*eye(2*T);
    for j=1:length(N{i})
        index = N{i}(j);
        pos = (2*T)*(index-1)+1;
        A{i}(pos:(pos+(2*T-1)),pos:(pos+(2*T-1))) = eye(2*T);
        A{i}(pos:(pos+(2*T-1)),(2*T*6+1):end) = -eye(2*T);
        A{i}((2*T*6+1):end,pos:(pos+(2*T-1))) = -eye(2*T);
    end
end

b = cell(1,5);

for i=1:5
    b{i} = zeros(2*(6*T+T),1);
    b{i}((2*T*6+1):end) = - V'*deltaZ(:,i);
    for j=1:length(N{i})
        index = N{i}(j);
        pos = (2*T)*(index-1)+1;
        b{i}(pos:(pos+(2*T-1))) = - deltaM{i}(:,j);
        b{i}((2*T*6+1):end) = b{i}((2*T*6+1):end) + deltaM{i}(:,j);
    end
end

c = cell(1,5);

for i=1:5
    c{i} = (1/2)*norm(deltaZ(:,i),2)^2;
    for j=1:length(N{i})
        c{i} = c{i} + (1/2)*norm(deltaM{i}(:,j),2)^2;
    end
end

% stack mi and zi, each mi is 6*2*5, zi is 5*2
T = 200;
x_FBS = zeros(70,T,5);
y_FBS = zeros(70,T,5);
iter = zeros(70,T);
diff_FBS = zeros(70,T,5);
diff_FBS_norm = zeros(1,T);
alpha = 0.2;
%FBS
for j=1:T
    for i=1:5
        y_FBS(:,j,i) = (eye(70)-alpha*A{i})*x_FBS(:,j,i)-alpha*b{i};
        iter(1:60,j) = iter(1:60,j) + y_FBS(1:60,j,i);
    end
    for i=1:5
        x_FBS(:,j+1,i) = [0.2*iter(1:60, j);% mean(m)
                          y_FBS(61:70, j, i)];
        diff_FBS(:,j,i) = x_FBS(:,j+1,i)-x_FBS(:,j,i);
    end
end
alpha = 10;
for j=1:T
    stat = [diff_FBS(:,j,1);
            diff_FBS(:,j,2);
            diff_FBS(:,j,3);
            diff_FBS(:,j,4);
            diff_FBS(:,j,5)];
    diff_FBS_norm(j) = norm(stat);
    if diff_FBS_norm(j) <= 0.00001 % thershold
    disp(['when k=', num2str(j), ', the norm of FBS is less than 0.00001.'])
    break
    end
end

% DRS
x_DRS = zeros(70,T,5);
x_half = zeros(70,T,5);
y_DRS = zeros(70,T,5);
iter = zeros(60,T);
diff_DRS = zeros(70,T,5);
diff_DRS_norm = zeros(1,T);
for j=1:T
    for i=1:5
        iter(:,j) = iter(:,j) + y_DRS(1:60,j,i);
    end
    for i=1:5
        x_half(:,j,i) = [0.2*iter(1:60, j);% mean(m)
                         y_DRS(61:70, j, i)];
        x_DRS(:,j+1,i) = (eye(70)+alpha*A{i})\(2*x_half(:,j,i)-y_DRS(:,j,i)-alpha*b{i});
        y_DRS(:,j+1,i) = y_DRS(:,j,i) + x_DRS(:,j+1,i) - x_half(:,j,i);
        diff_DRS(:,j,i) = x_DRS(:,j+1,i)-x_DRS(:,j,i);
    end
end
for j=1:T
   stat = [diff_DRS(:,j,1);
            diff_DRS(:,j,2);
            diff_DRS(:,j,3);
            diff_DRS(:,j,4);
            diff_DRS(:,j,5)];
    diff_DRS_norm(j) = norm(stat);
    if diff_DRS_norm(j) <= 0.00001 % thershold
    disp(['when k=', num2str(j), ', the norm of DRS is less than 0.00001.'])
    break
    end
end

%plot convergence
plot(1:T, diff_FBS_norm, 1:T, diff_DRS_norm, 'linewidth',2);
legend("FBS","DRS");
xlabel('$k$','interpreter','latex');
ylabel('$\|x^{k+1}-x^k\|_2$','interpreter','latex');

%obtain the m* and z*
mznew_FBS = x_FBS(:,T,1);
mznew_DRS = x_DRS(:,T,1);
for i=1:4
mznew_FBS = [mznew_FBS x_FBS(:,T,i+1)];
mznew_DRS = [mznew_DRS x_DRS(:,T,i+1)];
end

%get x,y of mi and zi from new result
mxnew_FBS = zeros(30,1);
mynew_FBS = zeros(30,1);
mxnew_DRS = zeros(30,1);
mynew_DRS = zeros(30,1);
for i =1:30
    mxnew_FBS(i,1) = mznew_FBS(2*i-1,1);
    mynew_FBS(i,1) = mznew_FBS(2*i,1);
    mxnew_DRS(i,1) = mznew_DRS(2*i-1,1);
    mynew_DRS(i,1) = mznew_DRS(2*i,1);
end
zxnew_FBS = zeros(5,5);
zynew_FBS = zeros(5,5);
zxnew_DRS = zeros(5,5);
zynew_DRS = zeros(5,5);
for i =31:35
    zxnew_FBS(i-30,:) = mznew_FBS(2*i-1,:);
    zynew_FBS(i-30,:) = mznew_FBS(2*i,:);
    zxnew_DRS(i-30,:) = mznew_DRS(2*i-1,:);
    zynew_DRS(i-30,:) = mznew_DRS(2*i,:);
end

%get x,y of mtrue and ztrue
mtruestack = [m(:,1); m(:,2); m(:,3);m(:,4);m(:,5);m(:,6)];
mtruex = zeros(30,1);
mtruey = zeros(30,1);
ztruex = zeros(5,5);
ztruey = zeros(5,5);
for i =1:30
    mtruex(i,1) = mtruestack(2*i-1,1);
    mtruey(i,1) = mtruestack(2*i,1);
end
for i=1:5
    ztruex(i,:) = ztrue(2*i-1,:);
    ztruey(i,:) = ztrue(2*i,:);
end
ztruexcol = [ztruex(:,1);
                  ztruex(:,2);
                  ztruex(:,3);
                  ztruex(:,4);
                  ztruex(:,5)];
ztrueycol = [ztruey(:,1);
                  ztruey(:,2);
                  ztruey(:,3);
                  ztruey(:,4);
                  ztruey(:,5)];

%compute "centre of mass"
Num = 5*5+30;
comx_true = (sum(ztruex,"all")+sum(mtruex,"all"))/Num;
comy_true = (sum(ztruey,"all")+sum(mtruey,"all"))/Num;

comx_FBS = (sum(mxnew_FBS,"all")+sum(zxnew_FBS,"all"))/Num;
comy_FBS = (sum(mynew_FBS,"all")+sum(zynew_FBS,"all"))/Num;

comx_DRS = (sum(mxnew_DRS,"all")+sum(zxnew_DRS,"all"))/Num;
comy_DRS = (sum(mynew_DRS,"all")+sum(zynew_DRS,"all"))/Num;

shiftx_FBS = comx_true-comx_FBS;
shifty_FBS = comy_true-comy_FBS;

shiftx_DRS = comx_true-comx_DRS;
shifty_DRS = comy_true-comy_DRS;

mxshift_FBS = mxnew_FBS+shiftx_FBS*ones(30,1);
myshift_FBS = mynew_FBS+shifty_FBS*ones(30,1);
zxshift_FBS = zxnew_FBS+shiftx_FBS*ones(5,5);
zyshift_FBS = zynew_FBS+shifty_FBS*ones(5,5);
zxshiftcol_FBS = [zxshift_FBS(:,1);
                  zxshift_FBS(:,2);
                  zxshift_FBS(:,3);
                  zxshift_FBS(:,4);
                  zxshift_FBS(:,5)];
zyshiftcol_FBS = [zyshift_FBS(:,1);
                  zyshift_FBS(:,2);
                  zyshift_FBS(:,3);
                  zyshift_FBS(:,4);
                  zyshift_FBS(:,5)];
mxshift_DRS = mxnew_DRS+shiftx_DRS*ones(30,1);
myshift_DRS = mynew_DRS+shifty_DRS*ones(30,1);
zxshift_DRS = zxnew_DRS+shiftx_DRS*ones(5,5);
zyshift_DRS = zynew_DRS+shifty_DRS*ones(5,5);
zxshiftcol_DRS = [zxshift_DRS(:,1);
                  zxshift_DRS(:,2);
                  zxshift_DRS(:,3);
                  zxshift_DRS(:,4);
                  zxshift_DRS(:,5)];
zyshiftcol_DRS = [zyshift_DRS(:,1);
                  zyshift_DRS(:,2);
                  zyshift_DRS(:,3);
                  zyshift_DRS(:,4);
                  zyshift_DRS(:,5)];

%scatter FBS
% h1=scatter(mxshift_FBS,myshift_FBS);%,zxshiftcol_FBS,zyshiftcol_FBS
% xlabel('$x$','interpreter','latex');
% ylabel('$y$','interpreter','latex');
% title('Measurement of obstacles and robots applying FBS');
% hold on;
% h2=scatter(zxshiftcol_FBS,zyshiftcol_FBS);
% h3=scatter(mtruex,mtruey,'*');
% h4=scatter(ztruexcol,ztrueycol,'*');
% legend([h1,h2,h3,h4],"obst_m","rob_m","obst_t","rob_t");
% hold off;

%scatter DRS
% h1=scatter(mxshift_DRS,myshift_DRS);
% xlabel('$x$','interpreter','latex');
% ylabel('$y$','interpreter','latex');
% title('Measurement of obstacles and robots applying DRS');
% hold on;
% h2=scatter(zxshiftcol_DRS,zyshiftcol_DRS);
% h3=scatter(mtruex,mtruey,'*');
% h4=scatter(ztruexcol,ztrueycol,'*');
% legend([h1,h2,h3,h4],"obst_m","rob_m","obst_t","rob_t");
% hold off;
