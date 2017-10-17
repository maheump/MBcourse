function [ PofA, PofB ] = LML_LocalMatchingLaw( pa, pb )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if numel(pa) ~= numel(pb)
    error('pa and pb must have the same size');
end

Va = pa;
Vb = pb;
V = pa + pb;

PofA = Va ./ V;
PofB = Vb ./ V;

end