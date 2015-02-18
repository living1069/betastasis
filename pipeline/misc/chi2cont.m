function [p,X2] = chi2cont(x,varargin)
% chi2cont chi-square test of contingency table
%   p = chi2cont(x) performs a chi-square test on the data in the
%   m-by-n contingency table x. The null hypothesis is that there is no difference
%   in the row variable distribution ('outcomes') between the columns 
%   ('treatments'). The result of the test is returned in h. h=1 indicates
%   a rejection of the null hypothesis at the 5% significance level.h=0
%   indicates that the null hypothesis can not be rejected at the 5%
%   significance level.
%
%   Reference http://www.psychstat.missouristate.edu/introbook/sbk28m.htm
%   Mark Snaterse, January 22 2014

e = sum(x,2)*sum(x)/sum(x(:));
X2 = (x-e).^2./e;
X2 = sum(X2(:));
df = prod(size(x)-[1 1]);
p = 1-chi2cdf(X2,df);
%if isnan(p), x, e, X2, end

