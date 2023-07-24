% ============================================================
% experimental_attitude.m
% 
% Experimental code to try and get attitude dynamics into
% MPC-MHE. Currently closely followign this paper.
% https://arxiv.org/pdf/2211.11653.pdf
% 
% Author: An Dang
%
% ============================================================

%{
    NOTES:

    Dynamics are the same as HCW but it now also adds Rotation to the thrusters

    In other words:
    ---------------
    statedot = A*state + Rot*B*u

    where Rot is SO(3)


    using MATLAB's zoh will be very useful.

    seems like having attitude be supported should stay within 3d.

    QUESTION: decouple translation and rotation??

    seems like the rotation is based on our pose relative to the target spacecraft.

    in this case, we jus include a rotation matrix as usual.

    when thruster misalignment happens, we induce a change in the rotation matrix
%}



disp("Hello World\n")