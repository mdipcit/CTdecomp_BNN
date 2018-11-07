% Author: Shunsuke Ono (ono@sp.ce.titech.ac.jp)

function [out] = EvalImgQuality(u, u_org, Qtype, varargin)

if ~isempty(varargin)
    dynamic = varargin{1};
else
    dynamic = 1;
end
if strcmp(Qtype, 'PSNR')
    MSE = sum(sum(sum((u - u_org).^2)));
    MSE= MSE/(numel(u));
    out = 10 * log10(dynamic^2/MSE);
elseif strcmp(Qtype, 'SSIM')
    out = ssim_index((u/dynamic)*255, (u_org/dynamic)*255);
end