clear all; close all

SetDefaultFigure

%
% ORTHOGONAL CORRELATION

% The disturbance can be:
% - coloured noise - the parameters are:
%        - s (standard deviation of the Gaussian white noise at the input of the noise filter)
%        - omega_g (the bandwith of the noise filter)
% - harmonic distrurbance - the parameters are:
%        - A_sin (amplitude of the disturbance)
%        - omega_g (angular frequency of the disturbance)
% - drift - the parameters are:
%        - Drift - the inclination of the ramp function Drift*t
%

% which disturbances to include (0 means they are not present)
s = 0;
A_sin = 0;
Drift =0;

% Parameters of the disturbance

% Do we want to remove the drift (detrend)
flag_detrend  = 0;

% Signal length
n=5000;

% Process model (continuous transfer function numerator and denominator)
b=0.02;
a=[1 0.05 0.02];

% Number and values of test frequencies
N = 50;
f=2*logspace(-3,-1,N);


ID_Ortogonalna_Korelacija_Izracun
