function Tm = GTVR_Tm(B,L,h0,Tm0,h,doy,c)

%获取目标点经纬度所在格网点四个角点经纬度
B1 = floor(B/2.5)*2.5;
B2 = ceil(B/2.5)*2.5;
L1 = floor(L/2.5)*2.5;
L2 = ceil(L/2.5)*2.5;

if B1==B2&&L1==L2
    
    a = c(c(:,1)==B&c(:,2)==L,:);
    
elseif B1==B2
    
    a1 = c(c(:,1)==B1&c(:,2)==L1,:);
    a2 = c(c(:,1)==B1&c(:,2)==L2,:);
    a = ((L2-L)*a1+(L-L1)*a2)/(L2-L1);    
    
elseif L1==L2
    
    a1 = c(c(:,1)==B1&c(:,2)==L1,:);
    a2 = c(c(:,1)==B2&c(:,2)==L1,:);
    a = ((B2-B)*a1+(B-B1)*a2)/(B2-B1);   
    
else
    %得到四个格网点的系数
    a1 = c(c(:,1)==B1&c(:,2)==L1,:);
    a2 = c(c(:,1)==B1&c(:,2)==L2,:);
    a3 = c(c(:,1)==B2&c(:,2)==L1,:);
    a4 = c(c(:,1)==B2&c(:,2)==L2,:);

    %双线性内插得到目标点系数
    a = ((B2-B)*(L2-L)*a1+(B-B1)*(L-L1)*a4+(B-B1)*(L2-L)*a3+(B2-B)*(L-L1)*a2)/((B2-B1)*(L2-L1));    
end

X = a(:,3:9)';

%计算设计矩阵
dh = h-h0;
len = length(h(:,1));
i = 1;
I = repmat(i,len,1);

C1 = cos(2*pi*doy/365.25);
S1 = sin(2*pi*doy/365.25);
C2 = cos(4*pi*doy/365.25);
S2 = sin(4*pi*doy/365.25);
C3 = cos(2*pi*h/19000)-cos(2*pi*h0/19000);
S3 = sin(2*pi*h/19000)-sin(2*pi*h0/19000);

H = [I.*dh,C1.*dh,S1.*dh,C2.*dh,S2.*dh,C3,S3];

Tm = H*X+Tm0;
