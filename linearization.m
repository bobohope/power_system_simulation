function [Asys, Bsys] = linearization(Pg,Qg,Vg,thg,Vbus,theta,PV,PQ,Pl,Ql)
%
%   
disp('Starting to read all the data');
%% pre-define all variables
global sys_base mac_con exc_con pss_con Ybus nbus

%================= machine data=========================
len_mac=size(mac_con,1);
mac_indx = mac_con(:,1);
mac_bus_indx = mac_con(:,2);
BM=mac_con(:,3)/sys_base; % uniform base conversion matrix
%xls=mac_con(:,4)./BM; % leakage reactance,                           unused
Ra=mac_con(:,5)./BM; % resistance,                                   unused
xd=mac_con(:,6)./BM; % d-axis synchronous reactance
xdd=mac_con(:,7)./BM; % d-axis transient reactance
%xddd=mac_con(:,8)./BM; % d-axis subtransient reactance
Td0d=mac_con(:,9); % d-axis time constant
%Td0dd=mac_con(:,10); % d-axis subtransient time constant             unused 
xq=mac_con(:,11)./BM; % q-axis reactance
xqd=mac_con(:,12)./BM; % q-axis transient reactance
%xqdd=mac_con(:,13)./BM; % q-axis subtransient reactance
Tq0d=mac_con(:,14); % q-axis time constant
%Tq0dd=mac_con(:,15); % q-axis subtransient time constant             unused
H=mac_con(:,16).*BM; % Inertial constant
D=mac_con(:,17).*BM; % damping
%D=[0;0;0]; 
f=60; % 60Hz base frequency
wB=2*pi*f; %base frequency in rad/s
M=2.*H/wB;
N_Machine=size(mac_con,1);

%===========machine end=================

%===========exciter data================
Tr=ones(N_Machine,1);
KA=zeros(N_Machine,1);
Ta=ones(N_Machine,1);
Ke=ones(N_Machine,1);
Te=ones(N_Machine,1);
Kf=zeros(N_Machine,1);
Tf=ones(N_Machine,1);
E1=zeros(N_Machine,1);
E2=zeros(N_Machine,1);
SEE1=zeros(N_Machine,1);
SEE2=zeros(N_Machine,1);
len_exc=size(exc_con,1);
Exc_indx = exc_con(:,2);
for i=1:len_exc
    Exc_m_indx=exc_con(i,2);%present machine index
    Tr(Exc_m_indx)=exc_con(i,3);
    Ta(Exc_m_indx)=exc_con(i,5);
    Ke(Exc_m_indx)=exc_con(i,10);
    Te(Exc_m_indx)=exc_con(i,11);
    Kf(Exc_m_indx)=exc_con(i,16);
    Tf(Exc_m_indx)=exc_con(i,17);
    KA(Exc_m_indx)=exc_con(i,4);
    E1(Exc_m_indx)=exc_con(i,12);
    E2(Exc_m_indx)=exc_con(i,14);
    SEE1(Exc_m_indx)=exc_con(i,13);
    SEE2(Exc_m_indx)=exc_con(i,15);
end
%================exciter end=================

%================pss1A=======================
Ks=zeros(N_Machine,1);
Tw=ones(N_Machine,1);
T1=ones(N_Machine,1);
T2=ones(N_Machine,1);
T3=ones(N_Machine,1);
T4=ones(N_Machine,1);
Vs_max=zeros(N_Machine,1);
Vs_min=zeros(N_Machine,1);
len_pss=size(pss_con,1);
for i=1:len_pss
    Pss_m_indx=pss_con(i,2);%present machine index
    Ks(Pss_m_indx)=pss_con(i,3);%pssgain
    Tw(Pss_m_indx)=pss_con(i,4);%washout time constant
    T1(Pss_m_indx)=pss_con(i,5);%first lead time constant
    T2(Pss_m_indx)=pss_con(i,6);%first lag time constant
    T3(Pss_m_indx)=pss_con(i,7);%second lead time constant
    T4(Pss_m_indx)=pss_con(i,8);%second lag time constant  
    Vs_max(Pss_m_indx)=pss_con(i,9);%maximum output limit
    Vs_min(Pss_m_indx)=pss_con(i,10);%minimum output limit
end
%==============pss end==========================
disp('Finish reading all data!');
disp('---------------------------');
disp('Start to calculate initial values');
%% ==============calculate initial condition============
% 1. calculate complex injection current for all generator bus (page 190)
j=sqrt(-1);
voltage = Vg.*(cos(thg) + j*sin(thg)); %Vg in complex form
current = conj((Pg+j*Qg)./voltage); % Ig in complex form 
% 2. compute delta0 (voltage angle reference) for generator
delta0 = angle(voltage + (Ra+j.*xq).*current); % in rad
% 3. compute w0
w0 = wB;
% 4. compute Id and Iq
id0 = abs(current) .* sin(delta0 - angle(current));
iq0 = abs(current) .* cos(delta0 - angle(current));
vd0 = abs(voltage) .* (sin(delta0 - angle(voltage)));
vq0 = abs(voltage) .* cos(delta0 - angle(voltage));
% 5. compute Eqd0,Edd0
Edd0 = (xq-xqd) .* iq0; 
assert(max(abs(vd0+Ra.*id0-xqd.*iq0 - Edd0))<10e-10);
Eqd0 = abs(voltage).*cos(delta0 - angle(voltage))+xdd.*id0;
assert(max(abs(vq0+Ra.*iq0+xdd.*id0 - Eqd0))<10e-10);
% 6. compute Efd
Efd0=Eqd0+(xd-xdd).*id0;
% compute Vref
Vref0 = Efd0; % no exciter
Vref0(Exc_indx,:) = Efd0(Exc_indx,:)./KA(Exc_indx,:) + Vg(Exc_indx,:);% with exciter
% compute Tm
Tm = Eqd0.*iq0+(xq-xdd).*id0.*iq0;
%============end initial condition calculation=======
disp('All the initial values has been calculated');
disp('------------------------------------');
%% =========get state space representation============


A1 = [];
B1 = [];
B2 = [];
E_1 = [];
C1=[];
D1=[];
D2=[];
C2=[];
for k=1:len_mac
    i=mac_indx(k);
    if find(Exc_indx==mac_indx(k)) %if the exciter exist, 7 states
        
        % calculate S_E(Efd)  and its derivative--page 70, equation 4.22
        Bxi = log(SEE1(i)/SEE2(i))/(E1(i)-E2(i));
        Axi = SEE1(i)/exp(Bxi*E1(i));
%         Axi = 0.0039; 
%         Bxi = 1.555;
        fsi = -(Ke(i)+Efd0(i)*Axi*Bxi*exp(Bxi*Efd0(i)) + Axi*exp(Bxi*Efd0(i)))/Te(i); 
        % get sub matrices
        % A1i, B1i, B2i, E1i (eqn 8.27), page 220
        A1i = [
            0 1          0            0            0                        0        0;          %delta 
            0 -D(i)/M(i) -iq0(i)/M(i) -id0(i)/M(i) 0                        0        0;          %w
            0 0          -1/Td0d(i)    0           1/Td0d(i)                0        0;          %Eqd
            0 0          0            -1/Tq0d(i)   0                        0        0;          %Edd
            0 0          0            0            fsi                      1/Te(i)  0;          %Efd
            0 0          0            0            -KA(i)*Kf(i)/Ta(i)/Tf(i) -1/Ta(i) KA(i)/Ta(i);%VR
            0 0          0            0            Kf(i)/(Tf(i))^2          0        -1/Tf(i);   %Rf
            ];
        B1i = [   %id iq
            0                                      0;                                      
            (iq0(i)*(xdd(i)-xqd(i))-Edd0(i))/M(i) (id0(i)*(xdd(i)-xqd(i))-Eqd0(i))/M(i);
            -(xd(i)-xdd(i))/Td0d(i)                0;
            0                                     (xq(i)-xqd(i))/Tq0d(i);
            0                                      0;
            0                                      0;
            0                                      0;
            ];
        B2i = [       % theta, V
            0    0;
            0    0;
            0    0;
            0    0;
            0    0;
            0    -KA(i)/Ta(i);
            0    0;
            ];
        E1i = [       % Tm Vref
            0         0;
            1/M(i)    0;
            0         0;
            0         0;
            0         0;
            0         KA(i)/Ta(i);
            0         0;
            ];
        
        % C1i, D1i, D2i (eqn 8.32)
        C1i =[
            -Vg(i)*cos(delta0(i)-thg(i)) 0 0 1 0 0 0;
            Vg(i)*sin(delta0(i)-thg(i))  0 1 0 0 0 0; 
            ];
        D1i =[
            -Ra(i)  xqd(i);
            -xdd(i) -Ra(i);
            ];
        D2i=[
            Vg(i)*cos(delta0(i)-thg(i))  -sin(delta0(i)-thg(i));
            -Vg(i)*sin(delta0(i)-thg(i)) -cos(delta0(i)-thg(i));
            ];
        % C2
        C2i= [id0(i)*Vg(i)*cos(delta0(i)-thg(i))-iq0(i)*Vg(i)*sin(delta0(i)-thg(i)) 0 0 0 0 0 0;
             -iq0(i)*Vg(i)*cos(delta0(i)-thg(i))-id0(i)*Vg(i)*sin(delta0(i)-thg(i)) 0 0 0 0 0 0];
        A1 =blkdiag(A1, A1i);
        B1 =blkdiag(B1, B1i);
        B2 =blkdiag(B2, B2i);
        E_1 =blkdiag(E_1, E1i);
        C1 =blkdiag(C1, C1i);
        D1 =blkdiag(D1, D1i);
        D2 =blkdiag(D2, D2i);
        C2 =blkdiag(C2, C2i);
    else %if the exciter do not exist, 4 states
        % only the first 4 equations
        % A1i, B1i, B2i, E1i (eqn 8.27), page 220
        A1i = [
            0 1          0            0            ;          %delta 
            0 -D(i)/M(i) -iq0(i)/M(i) -id0(i)/M(i) ;          %w
            0 0          -1/Td0d(i)    0           ;          %Eqd
            0 0          0            -1/Tq0d(i)   ;          %Edd
            ];
        B1i = [%Id, Iq
             0                                      0;                                       
            (iq0(i)*(xdd(i)-xqd(i))-Edd0(i))/M(i) (id0(i)*(xdd(i)-xqd(i))-Eqd0(i))/M(i);
            -(xd(i)-xdd(i))/Td0d(i)                0;
            0                                     (xq(i)-xqd(i))/Tq0d(i);
            ];
        B2i = zeros(4,2);
        E1i =  [  % Tm Vref
            0         0;
            1/M(i)    0;
            0         1/Td0d(i);
            0         0;
            ];
        
        % C1i, D1i, D2i (eqn 8.32)
        C1i =[
            -Vg(i)*cos(delta0(i)-thg(i)) 0 0 1;
            Vg(i)*sin(delta0(i)-thg(i))  0 1 0; 
            ];
        D1i =[
            -Ra(i)  xqd(i);
            -xdd(i) -Ra(i);
            ];
        D2i=[
            Vg(i)*cos(delta0(i)-thg(i))  -sin(delta0(i)-thg(i));
            -Vg(i)*sin(delta0(i)-thg(i)) -cos(delta0(i)-thg(i));
            ];
        % C2
        C2i= [id0(i)*Vg(i)*cos(delta0(i)-thg(i))-iq0(i)*Vg(i)*sin(delta0(i)-thg(i)) 0 0 0;
             -iq0(i)*Vg(i)*cos(delta0(i)-thg(i))-id0(i)*Vg(i)*sin(delta0(i)-thg(i)) 0 0 0];
        A1 =blkdiag(A1, A1i);
        B1 =blkdiag(B1, B1i);
        B2 =blkdiag(B2, B2i);
        E_1 =blkdiag(E_1, E1i);
        C1 =blkdiag(C1, C1i);
        D1 =blkdiag(D1, D1i);
        D2 =blkdiag(D2, D2i);
        C2 =blkdiag(C2, C2i);
    end
end

D3=[];
D4 = zeros(2*len_mac,2*len_mac);
D5 = zeros(2*len_mac,2*(nbus-len_mac));
%  D3, D4, D5 (eqn 8.36)
for l=1:len_mac
    p=mac_indx(l); % index of generator, usually 1,2,3
    m=mac_bus_indx(l); % index of generator's bus
    % D3
    D3i =[Vg(p)*sin(delta0(p)-thg(p))  Vg(p)*cos(delta0(p)-thg(p));
        Vg(p)*cos(delta0(p)-thg(p))  -Vg(p)*sin(delta0(p)-thg(p))]; % Id Iq
    D3 =blkdiag(D3, D3i);
    
    
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % D4 m x 2len_mac
    for k=1:len_mac
        q = mac_indx(k);
        n = mac_bus_indx(k);
        if p==q  
            % d_thetai
            D4(2*l-1,2*k-1)= -id0(p)*Vg(p)*cos(delta0(p)-thg(p))+iq0(p)*Vg(p)*sin(delta0(p)-thg(p));
            D4(2*l,2*k-1) = iq0(p)*Vg(p)*cos(delta0(p)-thg(p))+id0(p)*Vg(p)*sin(delta0(p)-thg(p));
            for kk = 1:nbus
                if kk ~= m
                D4(2*l-1,2*k-1)= D4(2*l-1,2*k-1) + Vg(p)* Vbus(kk)*abs(Ybus(m,kk))*sin(thg(p)-theta(kk)-angle(Ybus(m,kk)));
                D4(2*l,2*k-1)= D4(2*l,2*k-1) - Vg(p)* Vbus(kk)*abs(Ybus(m,kk))*cos(thg(p)-theta(kk)-angle(Ybus(m,kk)));
                end
            end
            % d_Vi
            D4(2*l-1,2*k)= id0(p)*sin(delta0(p)-thg(p))+iq0(p)*cos(delta0(p)-thg(p))...
                       - Vg(p)*abs(Ybus(m,m))*cos(thg(p)-thg(p)-angle(Ybus(m,m)));
            D4(2*l,2*k)= -iq0(p)*sin(delta0(p)-thg(p))+id0(p)*cos(delta0(p)-thg(p))...
                       - Vg(p)*abs(Ybus(m,m))*sin(thg(p)-thg(p)-angle(Ybus(m,m)));
            for kk = 1:nbus
                D4(2*l-1,2*k)=D4(2*l-1,2*k) - Vbus(kk)*abs(Ybus(m,kk))*cos(thg(p)-theta(kk)-angle(Ybus(m,kk)));
                D4(2*l,2*k)=D4(2*l,2*k)- Vbus(kk)*abs(Ybus(m,kk))*sin(thg(p)-theta(kk)-angle(Ybus(m,kk)));
            end
        else
            D4(2*l-1,2*k-1) = -Vg(p)*Vg(q)*abs(Ybus(m,n))*sin(thg(p)-thg(q)-angle(Ybus(m,n)));
            D4(2*l,2*k-1)= Vg(p)*Vg(q)*abs(Ybus(m,n))*cos(thg(p)-thg(q)-angle(Ybus(m,n)));
            D4(2*l-1,2*k) = -Vg(p)*abs(Ybus(m,n))*cos(thg(p)-thg(q)-angle(Ybus(m,n)));
            D4(2*l,2*k) =  -Vg(p)*abs(Ybus(m,n))*sin(thg(p)-thg(q)-angle(Ybus(m,n)));
        end
    
    end
    %D5 m x 2(nbus-len_mac)
    for kp=1:(nbus-len_mac)
        kk = PQ(kp);
        D5(2*l-1,2*kp-1) = -Vg(p)*Vbus(kk)*abs(Ybus(m,kk))*sin(thg(p)-theta(kk)-angle(Ybus(m,kk)));
        D5(2*l,2*kp-1) = Vg(p)*Vbus(kk)*abs(Ybus(m,kk))*cos(thg(p)-theta(kk)-angle(Ybus(m,kk)));
        D5(2*l-1,2*kp) = -Vg(p)*abs(Ybus(m,kk))*cos(thg(p)-theta(kk)-angle(Ybus(m,kk)));
        D5(2*l,2*kp) = -Vg(p)*abs(Ybus(m,kk))*sin(thg(p)-theta(kk)-angle(Ybus(m,kk)));
    end
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end
D6=zeros(2*(nbus-len_mac),2*len_mac);
D7=zeros(2*(nbus-len_mac),2*(nbus-len_mac));
% D6, D7 (eqn 8.40)
for l=1:(nbus-len_mac)
    % D6 2npq x 2npv
    m = PQ(l); % index of PQ bus
    for kg = 1:len_mac
        n = PV(kg);
        D6(2*l-1,2*kg-1)=-Vbus(m)*Vbus(n)*abs(Ybus(m,n))*sin(theta(m)-theta(n)-angle(Ybus(m,n)));
        D6(2*l,2*kg-1)= Vbus(m)*Vbus(n)*abs(Ybus(m,n))*cos(theta(m)-theta(n)-angle(Ybus(m,n)));
        D6(2*l-1,2*kg)=-Vbus(m)*abs(Ybus(m,n))*cos(theta(m)-theta(n)-angle(Ybus(m,n)));
        D6(2*l,2*kg)=-Vbus(m)*abs(Ybus(m,n))*sin(theta(m)-theta(n)-angle(Ybus(m,n)));
    end
    % D7 2npq x 2npq
    for k = 1:(nbus-len_mac)
        n = PQ(k);
        if m == n
            for kk = 1:nbus
              %d_theta  
            if kk ~= m
                D7(2*l-1,2*k-1)= D7(2*l-1,2*k-1)+ Vbus(m)*Vbus(kk)*abs(Ybus(m,kk))*sin(theta(m)-theta(kk)-angle(Ybus(m,kk)));
                D7(2*l,2*k-1)= D7(2*l,2*k-1)-Vbus(m)*Vbus(kk)*abs(Ybus(m,kk))*cos(theta(m)-theta(kk)-angle(Ybus(m,kk)));
            end
            end
            D7(2*l-1,2*k)= -Vbus(m)*abs(Ybus(m,n))*cos(theta(m)-theta(n)-angle(Ybus(m,n)));
            D7(2*l,2*k)= -Vbus(m)*abs(Ybus(m,n))*sin(theta(m)-theta(n)-angle(Ybus(m,n)));  
              for kk = 1:nbus
                  D7(2*l-1,2*k)= D7(2*l-1,2*k) - Vbus(kk)*abs(Ybus(m,kk))*cos(theta(m)-theta(kk)-angle(Ybus(m,kk)));
                  D7(2*l,2*k)= D7(2*l,2*k) -Vbus(kk)*abs(Ybus(m,kk))*sin(theta(m)-theta(kk)-angle(Ybus(m,kk)));
              end
        else
            D7(2*l-1,2*k-1)= -Vbus(m)*Vbus(n)*abs(Ybus(m,n))* sin(theta(m)-theta(n)-angle(Ybus(m,n)));
            D7(2*l,2*k-1) =  Vbus(m)*Vbus(n)*abs(Ybus(m,n))* cos(theta(m)-theta(n)-angle(Ybus(m,n)));
            D7(2*l-1,2*k)= -Vbus(m)*abs(Ybus(m,n))*cos(theta(m)-theta(n)-angle(Ybus(m,n)));
            D7(2*l,2*k)= -Vbus(m)*abs(Ybus(m,n))*sin(theta(m)-theta(n)-angle(Ybus(m,n)));
        end
    end
end
% Get Asys
Asys = A1-B1*inv(D1)*C1-(B2-B1*inv(D1)*D2)*inv(D4-D3*inv(D1)*D2-D5*inv(D7)*D6)*(C2-D3*inv(D1)*C1);
% get B matrix for different setup
Bsys= E_1;



end

