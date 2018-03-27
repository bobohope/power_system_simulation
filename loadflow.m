function [Vg,thg,Pgen,Qgen,Pl,Ql,pv,pq,V,theta]=loadflow;
%  Newton Raphson method for loadflow
%   Detailed explanation goes here

%% all variables
global nbus sys_base Ybus bus line
type = bus(:,8);
V = bus(:,2);
Pg = bus(:,4);
Pl = bus(:,6);
Ql = bus(:,7);
% pre allocate variables
theta = zeros(nbus,1);
Qg = zeros(nbus,1);

P = Pg-Pl;
Q = Qg-Ql;
Psp = P;
Qsp = Q;
G = real(Ybus);
B = imag(Ybus);
pv = find(type == 2 | type == 1);
pq = find(type == 3);
npv = length(pv);
npq = length(pq);
err = 1; % error
Iter = 1; % iteration counter
disp('Starting loadflow');
%% start the while loop
while (err>1e-12)
% calculate P and Q using V, theta, G and B
P = zeros(nbus,1);
Q = zeros(nbus,1);
for i = 1:nbus
    for k = 1:nbus
        P(i)= P(i) + V(i)*V(k)*(G(i,k)*cos(theta(i)-theta(k))+B(i,k)*sin(theta(i)-theta(k)));
        Q(i)= Q(i) + V(i)*V(k)*(G(i,k)*sin(theta(i)-theta(k))-B(i,k)*cos(theta(i)-theta(k)));
    end
end
% calculat the change from specified value
dPa = Psp-P;
dQa = Qsp-Q;
k = 1; % iterator for load bus
dQ = zeros(npq,1);
for i = 1:nbus
    if type(i) ==3
        dQ(k,1) = dQa(i);
        k = k+1;
    end
end
dP = dPa(type~=1); %dP for slack bus is 0...
M = [dP; dQ];   %mismatch vector

% Jacobian
% J1 = dP/d theta
J1 = zeros(nbus,nbus);
for i=1:nbus
    for k=1:nbus
        if k== i
            for n=1:nbus
                J1(i,k)=J1(i,k)+V(i)*V(n)*(-G(i,n)*sin(theta(i)-theta(n))+B(i,n)*cos(theta(i)-theta(n)));
            end
            J1(i,k)=J1(i,k)-V(i)^2*B(i,i);
        else
            J1(i,k)= V(i)*V(k)*(G(i,k)*sin(theta(i)-theta(k))-B(i,k)*cos(theta(i)-theta(k)));
        end
    end
end
J1=J1(type~=1,type~=1);
% J2 = dP/ dV , we only need PQ bus here because dP/dV for PV bus is 0
J2 = zeros(nbus,npq);
for i=1:nbus
    for k=1:npq
        n=pq(k); % get the index of pq bus
        if n == i
            for np=1:nbus % this n is inside this scope only
                J2(i,k)=J2(i,k)+V(np)*(G(i,np)*cos(theta(i)-theta(np))+B(i,np)*sin(theta(i)-theta(np)));
            end
            J2(i,k)=J2(i,k)+V(i)*G(i,i);
        else
            J2(i,k)=V(i)*(G(i,n)*cos(theta(i)-theta(n))+B(i,n)*sin(theta(i)-theta(n)));
        end
    end
end
J2 = J2(type~=1,:);
% J3 = dQ/ d theta, we only need PQ bus here
J3 = zeros(npq,nbus);
for i = 1:npq
    m = pq(i);% get the index of pq bus
    for k = 1:nbus
        if k==m
            for n=1:nbus
                J3(i,k)=J3(i,k)+V(m)*V(n)*(G(m,n)*cos(theta(m)-theta(n))+B(m,n)*sin(theta(m)-theta(n)));
            end
            J3(i,k)=J3(i,k)-V(m)^2*G(m,m);
        else
            J3(i,k)= V(m)*V(k)*(-G(m,k)*cos(theta(m)-theta(k))-B(m,k)*sin(theta(m)-theta(k)));
        end
    end
end
J3 = J3(:,type~=1);
% J4 = dQ/ dV
J4 = zeros(npq,npq);
for i= 1:npq
    m=pq(i);
    for k = 1:npq
        n=pq(k);
        if m==n
            for n = 1:nbus % this n is local scope
                J4(i,k)=J4(i,k)+V(n)*(G(m,n)*sin(theta(m)-theta(n))-B(m,n)*cos(theta(m)-theta(n)));
            end
            J4(i,k)=J4(i,k) - V(m)*B(m,m);
        else
            J4(i,k)=V(m)*(G(m,n)*sin(theta(m)-theta(n))-B(m,n)*cos(theta(m)-theta(n)));
        end
    end
end
J = [J1 J2;J3 J4];
% calculate mismatch
X = J\M; % the correction term
dtheta = X (1:nbus-1);
dV = X(nbus:end);
% update the state vector 
theta(type~=1) = dtheta + theta(type~=1);
V(pq) = dV+V(pq);
err= max(abs(M));
Iter = Iter+1;
end
fprintf('Finish! The number of iteration of loadflow program is %d.\n',Iter);
disp('------------------------------------');
%% get all outputs
Vg(1:npv,1)=V(pv); %terminal voltage of generator unchange
thg(1:npv,1)=theta(pv); %angle of generator, this will change
Pgen(1:npv,1) = P(pv)+Pl(pv); % P for generator
Qgen(1:npv,1) = Q(pv)+Pl(pv); % Q for generator
end

