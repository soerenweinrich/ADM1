function [x_in,q_in]= ADM1_interp_input(input,t)

% -------------------------------------------------------------------------
%
% Interpolate input values
% 
% Function corresponds to:
% q_in = interp1(input(:,1),input(:,2),t,'previous','extrap');
% x_in = interp1(input(:,1),input(:,3:end),t,'previous','extrap');
%
% -------------------------------------------------------------------------

% Get number of input sample points
num_t = size(input,1);
% Get number model components
num_input = size(input,2)-2;
% Find the last (previous) input value at t
for i=num_t:-1:1 
    if t>=double(input(i,1))
        q_in=input(i,2);
        x_in(1:num_input)=input(i,3:num_input+2);
        break; 
	elseif i==1
        q_in=0;
        x_in(1:num_input)=0; 
   end
end

