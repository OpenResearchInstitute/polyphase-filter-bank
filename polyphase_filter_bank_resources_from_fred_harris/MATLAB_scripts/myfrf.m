function [DH,DW] = myfrf(N, F, GF, W, A, diff_flag)
%---function [DH,DW] = remezfrf(N, F, GF, W, A, diff_flag)
%REMEZFRF Frequency Response Function for REMEZ.
%  REMEZ(N,F,A, ...) or
%  REMEZ(N,F,{'remezfrf',A}, ...) designs a linear-phase FIR filter
%  using REMEZ.
%
%  The symmetry option SYM defaults to 'even' if unspecified in the
%  call to REMEZ.
%
%  See also REMEZ.

%   Copyright 1988-2002 The MathWorks, Inc.
% $Revision: 1.6 $
%  modified by Karl Moerder SDSU, copyright 2003

%  [DH,DW]=REMEZFRF(N,F,GF,W,A,diff_flag)
%      N: filter order (length minus one)
%      F: vector of band edges
%     GF: vector of interpolated grid frequencies
%      W: vector of weights, one per band
%      A: vector of amplitudes of desired frequency response at band edges F
% diff_flag: ==1 for differentiator (1/f) weights, ==0 otherwise
%
%     DH: vector of desired filter response (mag & phase)
%     DW: vector of weights (positive)
%
% NOTE: DH(GF) and DW(GF) are specified as functions of frequency

% Support query by REMEZ for the default symmetry option:
if nargin==2,
  % Return symmetry default:
  if strcmp(N,'defaults'),
    DH = 'even';   % can be 'even' or 'odd'
    return
  end
end

if nargin < 6
    diff_flag = 0;
else%+++
    error('Differentiator option is not allowed with myfrf.')%+++
end

% Prevent discontinuities in desired function
for k=2:2:length(F)-2
    if F(k) == F(k+1)
        error('Adjacent bands not allowed.')
    end
end
if length(F) ~= length(A)
    error('Frequency and amplitude vectors must be the same length.')
end

nbands = length(A)/2;
l = 1;
while (l+1)/2 <= nbands
    sel = find( GF>=F(l) & GF<=F(l+1) );
    % desired magnitude is line connecting A(l) to A(l+1)
    if F(l+1)~=F(l)   % 
        slope=(A(l+1)-A(l))/(F(l+1)-F(l));
        DH(sel) = polyval([slope A(l)-slope*F(l)],GF(sel));
    else   % zero bandwidth band 
        DH(sel) = (A(l)+A(l+1))/2;  
    end
%---DW(sel) = W((l+1)/2) ./ (1 +(diff_flag & A(l+1) >= .0001)*(GF(sel)/2 - 1));
    if A(l+1) > 0.0001;%+++ check for passband or stopband
        DW(sel) = W((l+1)/2);%+++ passband
    else;%+++
        DW(sel) = W((l+1)/2) * (GF(sel)/GF(sel(1)));%+++ stopband
    end;%+++
    l = l + 2;
end