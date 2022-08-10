%% Single CV Step 
function [MO,I_ref,V_ref] = Profile_CV_Steps(t)
j = 0;
%% Set the steps 
% j = j + 1;
% MO(j) = 1;
% I_ref(j) = 0.001; % Current Value
% V_ref(j) = 0;
% Time(j)  = 0.001; % [s]

j = j + 1;
MO(j)    = 2;
I_ref(j) = 0; 
V_ref(j) = 1;
Time(j)  = 20; % [s]

j = j + 1;
MO(j)    = 2;
I_ref(j) = 0; % Current Value
V_ref(j) = 2;
Time(j)  = 20; % [s]

j = j + 1;
MO(j)    = 2;
I_ref(j) = 0; % Current Value
V_ref(j) = 3;
Time(j)  = 20; % [s]
% 
%% Make Reference Time Vector
time_vec = zeros(1,length(MO));
for i = 1:length(MO)
    if i == 1
        time_vec(i) = time_vec(i)+Time(i);
    else
        time_vec(i) = time_vec(i-1)+Time(i);
    end
end

%% Find which section we are in
idx = find(t<time_vec);

%% Set All outputs
if isempty(idx)
    MO    = MO(end);
    I_ref = I_ref(end);
    V_ref = V_ref(end);
else
    MO    = MO(idx(1));
    I_ref = I_ref(idx(1));
    V_ref = V_ref(idx(1));
end

end