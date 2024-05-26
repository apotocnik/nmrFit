%--------------------------------------------------------------------------
% sim_lib - fit function library
%
% Author: Anton Potocnik, F5, IJS
% Date:   28.01.2009 - 07.07.2014
% Arguments:
%       [simstr params] = sim_lib(name, number)
% Input:    name    ... simulation name ('Powder', 'Quadrupole',
% 'Quadrupole with jumping')
% Output:   funstr  ... string of simulation function
%           params  ... a cell array of parameters
%--------------------------------------------------------------------------
% NOTE!
% All simulations return spectra and frequency space: [spc, f] !!!
% Parameter order must be the same as function parametere order!!!
function [simstr params startVal] = nmrsim_lib(name)

switch name

    case 'Lorentzian'
        simstr = ['Lorentzian(fmax,NOP,A,LB,fc)'];
        params = {'fmax','NOP','A','LB','fc'};
        startVal = [600 1024 1 2 200];
        
    case 'Gaussian'
        simstr = ['Gaussian(fmax,NOP,A,LB,fc)'];
        params = {'fmax','NOP','A','LB','fc'};
        startVal = [600 1024 1 2 200];
        
    case 'Gaussian3'
        simstr = ['Gaussian3(fmax,NOP,A1,LB1,fc1,A2,LB2,fc2,A3,LB3,fc3)'];
        params = {'fmax','NOP','A1','LB1','fc1','A2','LB2','fc2','A3','LB3','fc3',};
        startVal = [600 1024 1 2 180 1 2 200 1 2 210];
    case 'Gaussian4'
        simstr = ['Gaussian4(fmax,NOP,A1,LB1,fc1,A2,LB2,fc2,A3,LB3,fc3,A4,LB4,fc4)'];
        params = {'fmax','NOP','A1','LB1','fc1','A2','LB2','fc2','A3','LB3','fc3','A4','LB4','fc4'};
        startVal = [600 1024 1 2 180 1 2 200 1 2 210 1 2 210];
        
    case 'Gaussian4link'
        simstr = ['Gaussian4link(fmax,NOP,A1,LB1,fc1,A2A1,LB2,fc2,A3A1,LB3,fc3,A4A1,LB4,fc4)'];
        params = {'fmax','NOP','A1','LB1','fc1','A2A1','LB2','fc2','A3A1','LB3','fc3','A4A1','LB4','fc4'};
        startVal = [600 1024 1 2 180 1 2 200 1 2 210 1 2 210];
        
    case 'Gaussian4link2'
        simstr = ['Gaussian4link2(fmax,NOP,A1,LB1,fc1,A2A1,LB2LB1,fc2,A3A1,LB3LB1,fc3,A4A1,LB4LB1,fc4)'];
        params = {'fmax','NOP','A1','LB1','fc1','A2A1','LB2LB1','fc2','A3A1','LB3LB1','fc3','A4A1','LB4LB1','fc4'};
        startVal = [600 1024 1 2 180 1 2 200 1 2 210 1 2 210];
        
    case 'Pwd Kanisotropy'
        simstr = ['nmrpowderF(N,fmax,NOP,A,LB,Kxx,Kyy,Kzz)'];
        params = {'N','fmax','NOP','A','LB','Kxx','Kyy','Kzz'};
        startVal = [10 600 1024 1 2 23 23 43];
        
    case 'Pwd KanisotropySym'
        simstr = ['nmrpowderFsym(N,fmax,NOP,A,LB,Kiso,Kax,Kasym)'];
        params = {'N','fmax','NOP','A','LB','Kiso','Kax','Kasym'};
        startVal = [10 600 1024 1 2 23 23 43];
        
   case 'Pwd KanisotropySymStrain'
        simstr = ['nmrpowderFsymStrain(N,fmax,NOP,A,LB,Kiso,Kax,Kasym,KaxS)'];
        params = {'N','fmax','NOP','A','LB','Kiso','Kax','Kasym','KaxS'};
        startVal = [10 600 1024 1 2 23 23 43,5];
        
    case 'Pwd 2KanisotropySym'
        simstr = ['nmrpowderFsym2(N,fmax,NOP,A1,LB1,K1iso,K1ax,K1asym,A2_A1,LB2,K2iso,K2ax,K2asym)'];
        params = {'N','fmax','NOP','A1','LB1','K1iso','K1ax','K1asym','A2_A1','LB2','K2iso','K2ax','K2asym'};
        startVal = [10 600 1024 1 2 198 100 30 1.5 2 170 20 0];
        
        
%     case 'Pwd Kanisotropy + Lorentz'
%         simstr = ['simAnisotropyLor1(N,DW,TD,A,LB,Kxx,Kyy,Kzz,B,Klor,W)'];
%         params = {'N','DW (us)','TD','A','LB (kHz)','Kxx (kHz)','Kyy (kHz)','Kzz (kHz)', 'B (ratio)', 'Klor', 'W (kHz)'};
%         startVal = [50 0.5 4096 1 2 12 -8 -8 1 0 2];

    case 'Pwd Quadrupole + Kanisotropy'
        simstr = ['qpole(N,fmax,NOP,Spin,A,LB,Kiso,Kani,Kasy,niQ,eta,DniQ,offset)'];
        params = {'N','fmax','NOP','Spin','A','LB','Kiso','Kani','Kasy','niQ','eta','DniQ','offset'};
        startVal = [20 1000 1024 1.5 1 10 0 0 0 1000 0 0 0];
        

    otherwise
        simstr = {'Lorentzian', ...
                  'Gaussian', ...
                  'Gaussian3', ...
                  'Gaussian4', ...
                  'Gaussian4link', ...
                  'Gaussian4link2', ...
                  'Pwd Kanisotropy', ...
                  'Pwd KanisotropySym', ...
                  'Pwd KanisotropySymStrain', ...
                  'Pwd 2KanisotropySym', ...
                  'Pwd Quadrupole + Kanisotropy', ...
                  };
        params = {};
        
end

