function Y = ybus
%  Calculate the system ybus
%  input data is line and bus data
global line nbus
i=sqrt(-1);
%% read data
FromIndex = line(:,1);
ToIndex = line(:,2);
R = line(:,3);
X = line(:,4);
B = line(:,5); % each end assign B/2 !!!\
tap = line(:,6); 
tap(tap == 0) = 1;
Z = R + i*X; % impedance of that line
y = 1./Z; % admintance of that line
nline = length(FromIndex); % number of lines

%% form Ybus
Y=zeros(nbus,nbus);
% off diagonal elements
for k=1:nline
    Y(FromIndex(k),ToIndex(k))= Y(FromIndex(k),ToIndex(k))-y(k)/tap(k);
    Y(ToIndex(k),FromIndex(k))= Y(FromIndex(k),ToIndex(k));
end
% diagonal elements
for m=1:nbus
    for n=1:nline
        if FromIndex(n)==m
            Y(m,m)=Y(m,m)+y(n)/(tap(n)^2) + i*B(n)/2;
        elseif ToIndex(n)==m
            Y(m,m)=Y(m,m)+y(n)+i*B(n)/2;
        end
    end
end


end

