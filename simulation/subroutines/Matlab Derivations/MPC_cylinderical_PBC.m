%% Abaqus MPC derivations

clear all
close all

%% Symbolic structures
NDEP = 6;
MDOF = 6;
N = 3;

% Positions
xa = sym('xa_',[1,MDOF],'Real');
xb = sym('xb_',[1,MDOF],'Real');

% Displacements
ua = sym('ua_',[1,MDOF],'Real');
ub = sym('ub_',[1,MDOF],'Real');
urp = sym('urp_',[1,MDOF],'Real');


% Rotations in teh form
%     theta_B = sqrt(UR1_B^2 + UR2_B^2 + UR3_B^2);
%     omega_B = [UR1_B/theta_B; UR2_B/theta_B; UR3_B/theta_B];
theta_B = sqrt((ub(4))^2+(ub(5))^2+(ub(6))^2);
omega_B = [(ub(4))/theta_B; (ub(5))/theta_B; (ub(6))/theta_B];
Z = [0; 0; 1];

% Node A
Xa_R = sqrt(xa(1)^2+xa(2)^2);
Xa_Phi = atan2(xa(2),xa(1));
Xa_Z = xa(3);
Ua_R = sqrt((xa(1)+ua(1))^2+(xa(2)+ua(2))^2)-Xa_R;
Ua_Phi = atan2((xa(2)+ua(2)),(xa(1)+ua(1)))-Xa_Phi;
Ua_Z = ua(3);

% Node B
Xb_R = sqrt(xb(1)^2+xb(2)^2);
Xb_Phi = atan2(xb(2),xb(1));
Xb_Z = xb(3);
Ub_R = sqrt((xb(1)+ub(1))^2+(xb(2)+ub(2))^2) - Xb_R;
Ub_Phi = atan2((xb(2)+ub(2)),(xb(1)+ub(1))) - Xb_Phi;
Ub_Z = ub(3);

% Apply rotation
[theta_A,omega_A] = CompoundAxAng( ...
   theta_B,omega_B, ...
    Xa_Phi-Xb_Phi,Z);
UR_A = theta_A*omega_A;

f1 = ua(1)-sqrt(ub(1)*ub(1)+ub(2)*ub(2))*cos(atan2(ub(2),ub(1))+((Xa_Phi+Ua_Phi)-(Xb_Phi+Ub_Phi)));
f2 = ua(2)-sqrt(ub(1)*ub(1)+ub(2)*ub(2))*sin(atan2(ub(2),ub(1))+((Xa_Phi+Ua_Phi)-(Xb_Phi+Ub_Phi)));
f3 = Ua_Z-Ub_Z-(Xa_Z-Xb_Z)*urp(3);
f4 = ua(4)-UR_A(1);
f5 = ua(5)-UR_A(2);
f6 = ua(6)-UR_A(3);
f = [f1,f2,f3,f4,f5,f6];
u = [transpose(ua),transpose(ub),transpose(urp)];

%% Derivations
NDEP = length(f);
MDOF = length(f);
N = size(u,2);

for n = 1:size(u,2)
    for i = 1:length(f)
        for j = 1:length(f)
            A(i,j,n) = diff(f(i),u(j,n));
        end
    end
end
A_inv = inv(A(:,:,1));

%% Write MPC_File subroutine file automatically

fid = fopen('MPC_equations.f','w');
sol1 = solve(f(1),ua(1));
fprintf(fid,fortranFormat(...
sprintf('UE(%d)',1), -(f(1)-ua(1)),0,false));
sol2 = solve(f(1),ua(2));
fprintf(fid,fortranFormat(...
sprintf('UE(%d)',2), -(f(2)-ua(2)),0,false));
fprintf(fid,fortranFormat(...
sprintf('UE(%d)',3),solve(f(3),ua(3)),0,false));
fprintf(fid,fortranFormat(...
sprintf('UE(%d)',4),UR_A(1),0,false));
fprintf(fid,fortranFormat(...
sprintf('UE(%d)',5),UR_A(2),0,false));
fprintf(fid,fortranFormat(...
sprintf('UE(%d)',6),UR_A(3),0,false));

for n = 1:N
    for i = 1:MDOF
        for j = 1:MDOF
            if A(i,j,n) ~= 0
                fprintf(fid,fortranFormat(...
                sprintf('A(%d,%d,%d)',i,j,n),A(i,j,n),0,false));
            end
        end
    end
end

% Specify the file name
filename = 'MPC_temp.f';

% Open the file for reading
fid = fopen(filename, 'r');
if fid == -1
    error('File cannot be opened');
end

% Read the entire file into a character array
file_contents = fread(fid, '*char')';
fclose(fid);

% Open the file for writing
fid = fopen(filename, 'w');
if fid == -1
    error('File cannot be opened for writing');
end

% Write the modified contents back to the file
fwrite(fid, file_contents);
fclose(fid);

disp('File has been updated.');


%% Helper functions:
%
% Source:
%
% https://math.stackexchange.com/questions/382760/composition-of-two-axis-angle-rotations
%
% https://i.stack.imgur.com/nIVLC.jpg
% 

function [gamma,N] = CompoundAxAng(alpha,L, beta,M)
    %CompoundAxAng Compounds two rotations given by axis-angle pairs
    %   alpha - original angle of rotation about axis L (scalar, in radians)
    %   L - original axis of rotation (3x1 vector)
    %   beta - angle of additional rotation about axis M (scalar, in radians)
    %   M - axis about which to apply the additional rotation (3x1 vector)

    % Normalize the axes to ensure they are unit vectors
    L = L / norm(L);
    M = M / norm(M);

    % Rodrigues' rotation formula to rotate L around M by beta
    cosBeta = cos(beta);
    sinBeta = sin(beta);
    MdotL = dot(M, L);
    MCrossL = cross(M, L);

    % Compute the new axis N after rotation
    N = L * cosBeta + MCrossL * sinBeta + M * MdotL * (1 - cosBeta);

    % The new angle remains the original alpha for the simplicity of the example,
    % since compounding angles generally needs more context about how rotations combine.
    gamma = alpha;  % Returning the original rotation angle unchanged

end

% Fortran syntax writing functions
function out = fortranComment(str,numIndent)
    out = ['C     ',repelem('  ',numIndent),str,'\n'];
end

function out = fortranFormat(varName,expr,numIndent,increment)
    % Convert symbolic expression to string and clean up formatting
    % Perform replacements for ua_ variables
    str = erase(fortran(expr),{'      t0 = ','     &',newline,char(13)});
    for i = 1:6
        old_str = sprintf('ua_%d', i);
        new_str = sprintf('U(%d,1)', i);
        str = strrep(str, old_str, new_str);
    end

    % Perform replacements for xa_ variables
    for i = 1:6
        old_str = sprintf('xa_%d', i);
        new_str = sprintf('X(%d,1)', i);
        str = strrep(str, old_str, new_str);
    end

    % Perform replacements for ua_ variables
    for i = 1:6
        old_str = sprintf('ub_%d', i);
        new_str = sprintf('U(%d,2)', i);
        str = strrep(str, old_str, new_str);
    end

    % Perform replacements for xa_ variables
    for i = 1:6
        old_str = sprintf('xb_%d', i);
        new_str = sprintf('X(%d,2)', i);
        str = strrep(str, old_str, new_str);
    end
   
    % Perform replacements for xa_ variables
    for i = 1:6
        old_str = sprintf('urp_%d', i);
        new_str = sprintf('U(%d,3)', i);
        str = strrep(str, old_str, new_str);
    end
    firstLine = true;
    newString = '';
    while ~isempty(str)
        if firstLine
            if increment
                startString = [...
                    '      ',repelem('  ',numIndent),varName,' = ',...
                    varName,' + '];
            else
                startString = [...
                    '      ',repelem('  ',numIndent),varName,' = '];
            end
            firstLine = false;
        else
            startString = ['     *',repelem('  ',numIndent+1)];
        end
        maxWidth = 72 - length(startString);
        
        if length(str) > maxWidth
            ind = maxWidth;
            while true
                if any(isstrprop(str(ind:ind+1),'punct')) && ...
                        ~all(str(ind:ind+1)=='*') && ...
                        ~any(str(ind:ind+1)=='.')
                    break;
                else
                    ind = ind-1;
                end
            end
            newString = sprintf('%s%s%s\n',...
                newString,startString,str(1:ind));
            str = str(ind+1:end);
        else
            newString = sprintf('%s%s%s\n',...
                newString,startString,str);
            str = '';
        end
    end
    out = newString;
end
