%% data
data = zeros(6,10,9);

data(:,:,1) = [208    71   244   202   173   180   177   195   181    30;
            231   139   124   245   193     8    81   203   192   127;
            32   244   204   167   189    71   242    48    70   245;
            233   246    36     9   100    12     9   125   173    87;
            161    40   108   217   167    25   112   114   167   149;
            25   248   234   238    44   210    97   165    41    57];

data(:,:,2) = [192   140   208   157   234    19   145    79   176    39;
            65    35    62   121    73    14   120   135   191   211;
            129    38   237    90   193   135     3    42   115   137;
            178    66    89   212   192   199    86   154    21   254;
            227   214    50   149    97   238    41    67    58    20;
            245    65    64   140   145    33   203   167   233   113];

data(:,:,3) = [27    22    46   140   102   106    86    62   147    11;
            245   102    67    37    19    13   230   103    15    43;
            1    66    37   218    61   230    94    25    60   166;
            198   204    35   159    31   241    28    34    90   187;
            208   110   222    89    47   125   199   240   209   165;
            222   232   148   131    61   125    99   244     4   115];

data(:,:,4) = [139    94   124   208    89    53    58   110    66    57;
            76   160   111   203   239    77    44    47   104    30;
            190   199   114   164   223   120    58   231   152    76;
            48    21    78    97   140    59   111   250    67    81;
            175   237   130   207   159   215    79   112   154   108;
            47   198   130   136   150    50   235    28   181   130];


data(:,:,5) = [22   125   133    94    25    27   227   128   205   125;
            67   148    59   252    67   167    85   122   147    43;
            204    61   125    10    86   126   178   231    47   250;
            7   117   159   226   173   199    50   156    61   182;
            237   246   173   233    35   182     8   158   226   128;
            186   139   101   203   184   230   190   219     7   120];

data(:,:,6) = [15   209   248    21    15    74    95    13   107   178;
            174   208   165    34   102   110    51   188   251   170;
            11   184   204    44   134     4   125    69    77    45;
            18    38   116   100   106   251    87   108   179    33;
            133   168   110   212   167    43   243   140   170   255;
            25   132   210   205   160    27   235   240   137    44];

data(:,:,7) = [8   117    49    98   210   231   108   153    18   183;
            143   250   109   149   251   224    80   120    81   247;
            225    40   123    64   186   209    41   177   135   135;
            171   218    31    74    88    66    46   178   167    83;
            49   164   150   157   149   152   108   163   104    27;
            94    96    58    68    27     6    24     9   209   156];

data(:,:,8) = [199   112   163   177    88   234    82   121    62    23;
            108   134   244    17   199     0   200    39   234   147;
            23   117    61    65   172   118   120    87    69   174;
            68   223   172    57     2   108     9   155   195   139;
            39   132    74   170   154   118    45    49    48   109;
            72   241   171   215    99   196   184   188    73   164];

data(:,:,9) = [165    60   196    65    81   139    56    93    49   220;
            173    30    89   156    30   165    27   195    35   124
            162   155   169   148   240   139    28   160   178   100
            241   115   106   138   165   184    16   197    24   171
            53   117   215   222   122   133   103   238   134   189
            181   169   212    68   163   253   114   248   135   133];

labels = [-1     1     1    -1     1    -1    -1    -1     1;
          1    -1    -1     1     1     1    -1     1    -1;
          1    -1    -1    -1    -1     1    -1    -1     1;
          1    -1    -1     1     1    -1    -1    -1    -1;
          -1    -1    -1     1    -1    -1     1     1     1;
          -1    -1    -1    -1     1    -1    -1    -1     1;
          -1    -1    -1    -1    -1     1    -1    -1    -1;
          -1    -1     1    -1     1     1    -1    -1     1;
          -1     1    -1     1    -1     1    -1    -1    -1;
          -1    -1    -1    -1     1    -1    -1     1    -1];
%% algorithm
%initialization
data = data/255;% normalizing data!! Answer: we don't want the gradient become too large or too small.

T = 4000;% set iteration time, T_0.5 = 10000, T_10 = 600
omega = zeros(6,T);%
b = zeros(1,T);
x = [omega; b];
gf = zeros(7,T);
norm_gf = zeros(1,T);
f_value = zeros(1,T);
alpha = 0.5;% change alpha to 0.5, 10

% the following code is used for problem c
% for k =1:T
%     omega = x(1:6,k);
%     b = x(7,k);
%     [gf(:,k),f_value(1,k)] = twofunction(data,labels,omega,b);
%     x(:,k+1) = x(:,k) - alpha*gf(:,k);% update
%     norm_gf(k)=norm(gf(:,k),2);% norm of gradient f
% if norm_gf(k) <= 0.01
%     disp(['when k=', num2str(k), ', the norm of gradient is less than 0.01.'])
%     break
% end
% end
% %plot gradient
% % plot(1:T,norm_gf);
% % xlabel('time');
% % ylabel('\nablaf');
% %plot f(x)
% plot(1:T,f_value,'r');
% xlabel('time');
% ylabel('f(x)');
% % test classification
% A = test(data,omega,b); % calculate the value of omega'x+b
% B = A-labels;

%prombelm d, please comment the above code out EXCEPT initialization
rho=0.0001;%1,0.0001
for k =1:T
    omega = x(1:6,k);
    b = x(7,k);
    [gf(:,k),~] = twofunction(data,labels,omega,b);
    gf(:,k) = gf(:,k)+rho*[omega;0];
    [~,f_value(1,k)] = twofunction(data,labels,omega,b);
    f_value(1,k) = f_value(1,k)+rho*(omega'*omega)/2;
    x(:,k+1) = x(:,k) - alpha*gf(:,k);% update
    norm_gf(k)=norm(gf(:,k),2);% norm of gradient f
if norm_gf(k) <= 0.01
    disp(['when k=', num2str(k), ', the norm of gradient is less than 0.01.'])
    break
end
end
%plot gradient
% plot(1:T,norm_gf);
% xlabel('time');
% ylabel('\nablaf');
%plot f(x)
plot(1:T,f_value,'r');
xlabel('time');
ylabel('f(x)');
% test classification
A = test(data,omega,b); % calculate the value of omega'x+b
B = A-labels;
%% Functions

% Function for caculating gradient of f
function [grad_f,f_value] = twofunction(data,labels,omega,b)
grad_f_i = zeros(7,10,9);
grad_f_i_sum = zeros(7,9);
f_i_value = zeros(1,10,9);
f_i_value_sum = zeros(1,9);
for i=1:9
    for j=1:10
    grad_f_i(:,j,i) = [-labels(j,i)*data(:,j,i)/(exp(labels(j,i)*(omega'*data(:,j,i)+b))+1);
                   -labels(j,i)/(exp(labels(j,i)*(omega' *data(:,j,i)+b))+1)];
    f_i_value(1,j,i) = log(1+exp(-labels(j,i)*(omega' *data(:,j,i)+b)));
    end
    grad_f_i_sum(:,i) = sum(grad_f_i(:,:,i),2);% rowsum
    f_i_value_sum(:,i) = sum(f_i_value(:,:,i),2);% rowsum
end
grad_f = sum(grad_f_i_sum,2);%rowsum
f_value = sum(f_i_value_sum,2);%rowsum
end

% test the classifier
function A = test(data,omega,b)
A = zeros(10,9);
for i=1:9
    for j=1:10
if omega'*data(:,j,i)+b>0
A(j,i) = 1;
else 
    A(j,i) = -1;
end
    end
end
end
