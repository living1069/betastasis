function pval = logrank(varargin)

pval = nan;

if length(varargin) == 4
	Am(:, 1) = varargin{1};
	Am(:, 2) = varargin{2};
	Bm(:, 1) = varargin{3};
	Bm(:, 2) = varargin{4};
elseif length(varargin) == 2
	A = varargin{1};
	B = varargin{2};
	
	if isfield(A, 'Meta') && isfield(A.Meta, 'Patient') && ...
		isfield(B, 'Meta') && isfield(B.Meta, 'Patient')
		A = A.Meta.Patient;
		B = B.Meta.Patient;
	elseif isfield(A, 'Patient') && isfield(B, 'Patient')
		A = A.Patient;
		B = B.Patient;
	else
		error 'No survival time metadata available.';
	end
	
	Am(:, 1) = A.SurvivalTime;
	Am(:, 2) = A.Censored;
	Bm(:, 1) = B.SurvivalTime;
	Bm(:, 2) = B.Censored;
else
	error 'Invalid amount or aguments.';
end

Am = Am(~any(isnan(Am), 2), :);
Bm = Bm(~any(isnan(Bm), 2), :);

if numel(Am) == 0 || numel(Bm) == 0
	fprintf('Survival data was not available for both groups.\n');
	return;
end

if max(Am(:, 1)) > 365 || max(Bm(:, 1)) > 365
	fprintf(['Data contains extended survival times. ' ...
		'Showing the time axis in months...\n']);
	Am(:, 1) = Am(:, 1) / 30;
	Bm(:, 1) = Bm(:, 1) / 30;
end

pval = logrank_giuseppe(Am, Bm);
return;




function p = logrank_giuseppe(varargin)
% LOGRANK Comparing survival curves of two groups using the log rank test
% Comparison of two survival curves can be done using a statistical
% hypothesis test called the log rank test. It is used to test the null
% hypothesis that there is no difference between the population survival
% curves (i.e. the probability of an event occurring at any time point is
% the same for each population). This function use the Kaplan-Meier
% procedure to estimate the survival function, so it is mandatory to download
% KMPLOT (http://www.mathworks.com/matlabcentral/fileexchange/22293).
%
% Syntax: 	logrank(x1,x2,alpha)
%      
%     Inputs:
%           X1 and X2 (mandatory)- Nx2 data matrix:
%                     (X:,1) = survival time of the i-th subject
%                     (X:,2) = censored flag 
%                             (0 if not censored; 1 if censored)
%           note that if X is a vector, all the flags of the second column
%           will be set to 0 (all data are not censored).
%           ALPHA (optional) - significance level (default 0.05) 
%     Outputs:
%           Kaplan-Meier plot
%           Log-rank statistics
%
%      Example: 
%           load logrankdata x1 x2
%           logrank(x1,x2)
%
%LOG-RANK TEST FOR KAPLAN-MEIER SURVIVAL FUNCTIONS
%
%--------------------------------------------------------------------------------
%UL				S.E.			z				p-value			alpha
%--------------------------------------------------------------------------------
%6.57226		2.80788			2.16258			0.01529			0.050
%--------------------------------------------------------------------------------
%		The survival functions are statistically different
%
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo-edta@poste.it
%
% To cite this file, this would be an appropriate format:
% Cardillo G. (2008). LogRank: Comparing survival curves of two groups
% using the log rank test
% http://www.mathworks.com/matlabcentral/fileexchange/22317

%Input Error handling
args=cell(varargin);
nu=numel(args);
if isempty(nu) || nu<2
    error('Warning: Data vectors are required')
elseif nu>3
    error('Warning: Max three input data are required')
end
default.values = {[],[],0.05};
default.values(1:nu) = args;
[x1 x2 alpha] = deal(default.values{:});
if ~all(isfinite(x1(:))) || ~all(isnumeric(x1(:))) ...
        || ~all(isfinite(x2(:))) || ~all(isnumeric(x2(:)))
    error('Warning: all X1 and X2 values must be numeric and finite')
end
if isvector(x1) 
    x1(:,2)=0;
else
    if ~isequal(size(x1,2),2)
        error('LOGRANK requires Nx2 matrix data.');
    end
    if ~all(x1(:,2)==0 | x1(:,2)==1)
        error('Warning: all X1(:,2) values must be 0 or 1')
    end
end
if isvector(x2) 
    x2(:,2)=0;
else
    if ~isequal(size(x2,2),2)
        error('LOGRANK requires Nx2 matrix data.');
    end
    if ~all(x2(:,2)==0 | x2(:,2)==1)
        error('Warning: all X2(:,2) values must be 0 or 1')
    end
end
if nu>2
    if ~isscalar(alpha) || ~isnumeric(alpha) || ~isfinite(alpha) || isempty(alpha)
        error('Warning: it is required a numeric, finite and scalar ALPHA value.');
    end
    if alpha <= 0 || alpha >= 1 %check if alpha is between 0 and 1
        error('Warning: ALPHA must be comprised between 0 and 1.')
    end
end
clear args default nu

%recall KMPLOT function to construct tables of data (table1 and table2),
%tables of censored data (table12 and table 22), Kaplan-Meier variables
%(t1, t2, T1 and T2) and Kaplan-Meier graphical data for censored data 
%(xcg and ycg).
try
    [table1 table12 t1 T1 xcg1 ycg1]=kmplot(x1,0.05,0);
    [table2 table22 t2 T2 xcg2 ycg2]=kmplot(x2,0.05,0);
catch ME
    disp('Download KMPLOT: http://www.mathworks.com/matlabcentral/fileexchange/22293')
    rethrow(ME);
end

T1 = 100 * T1;
T2 = 100 * T2;
ycg1 = 100 * ycg1;
ycg2 = 100 * ycg2;

%plot both Kaplan-Meier curves
%clf
hold on
S1=stairs(t1,T1,'b', 'LineWidth', 1); %Kaplan-Meier curve for treatment 1

S2=stairs(t2,T2,'r', 'LineWidth', 1); %Kaplan-Meier curve for treatment 2
if ~isempty(table12)
    S3=plot(xcg1,ycg1,'k+'); %Censored data for treatment 1 (if there are)
else
    S3=[];
end

if ~isempty(table22)
    S3=plot(xcg2,ycg2,'k+'); %Censored data for treatment 2 (if there are)
end
hold off
%set the axis properly
xmax=max([t1;t2])+1;
axis([0 xmax 0 100]);
axis square
%add labels and legend
title('Kaplan-Meier estimate of survival functions')
ylabel('Estimated survival functions')
xlabel('Time')
if isempty(S3)
    legend([S1 S2],'Treatment 1','Treatment 2')
else
    legend([S1 S2 S3],'Treatment 1','Treatment 2','Censored')
end
clear S1 S2 S3 xmax T1 T2 t1 t2 xcg1 ycg1 xcg2 ycg2

%Full-blown LOGRANK procedure
%Merge the first columns of Table1 and Table2 (time intervals)
%and pick-up unique values
A=unique([table1(:,1);table2(:,1)]);
table=zeros(length(A),9); %matrix preallocation
%Out in the first column the time intervals
table(:,1)=A; 
%Put in the columns 2 and 3 and in the proper rows the deaths and alive
%taken from table1 columns 2 and 3
[c ia ib]=intersect(table1(:,1),A);
table(ib,2:3)=table1(ia,2:3);
%Put in the columns 4 and 5 and in the proper rows the deaths and alive
%taken from table2 columns 2 and 3
[c ia ib]=intersect(table2(:,1),A);
table(ib,4:5)=table2(ia,2:3);
%remove the rows where there arent't deaths in both treatments
table((table(:,2)==0 & table(:,4)==0),:)=[];
clear A c ia ib table1 table2
%fill the "pigeon-holes"
c=find(table(:,3)==0); %find the "pigeon-holes" of treatment 1
for I=1:length(c)
    if c(I)~=1
        %find the first interval time before the hole where there is almost 1
        %death
        J=find(table(1:c(I)-1,3)>0,1,'last');
        table(c(I),3)=table(J,3)-table(J,2);
        if ~isempty(table12)
        %find eventually censored data
            K=find((table12(:,1)<table(c(I),1) & table12(:,1)>=table(J,1)),1,'last');
            %Put in the hole how many subject were alive before the interval time
            %of the hole
            if ~isempty(K)
                table(c(I),3)=table(c(I),3)-sum(table12(K,2));
            end
        end
    else
        table(1,3)=length(x1);
    end
end
%Do the same for tratment 2
c=find(table(:,5)==0);
for I=1:length(c)
    if c(I)~=1
        J=find(table(1:c(I)-1,5)>0,1,'last');
        table(c(I),5)=table(J,5)-table(J,4);
        if ~isempty(table22)
            K=find((table22(:,1)<table(c(I),1) & table22(:,1)>=table(J,1)),1,'last');
            if ~isempty(K)
                table(c(I),5)=table(c(I),5)-sum(table22(K,2));
            end
        end
    else
        table(1,5)=length(x2);
    end
end
clear c I J K table12 table22

%Fill the table and compute the statistic variable
%Compute the total deaths and alive before the i-th time interval
table(:,6:7)=[sum(table(:,[2 4]),2) sum(table(:,[3 5]),2)];
%Compute the difference between observed deaths for treatment 1 and
%expected deaths in the hyphthesis that the treatments are similar
table(:,8)=table(:,2)-table(:,3).*table(:,6)./table(:,7);
%Log-rank statistic is the sum of column 8 values
UL=abs(sum(table(:,8)));
%Compute the contribute to the standard error
table(:,9)=prod(table(:,[3 5 6]),2).*(table(:,7)-table(:,6)) ...
    ./(table(:,7).^2.*(table(:,7)-ones(size(table,1),1)));
%find if there is some NaN (i.e. 0/0)
loc=isnan(table(:,9));
if any(loc)
    table(loc,9)=0;
end
SUL=sqrt(sum(table(:,9))); %Compute the totale standard error
z=abs((UL-0.5)/SUL); %normalized UL with Yates'es correction
p=1-0.5*erfc(-z/realsqrt(2)); %p-value

%display results
disp('LOG-RANK TEST FOR KAPLAN-MEIER SURVIVAL FUNCTIONS')
disp(' ')
tr=repmat('-',1,80);
disp(tr)
fprintf('UL\t\t\tS.E.\t\tz\t\tp-value\t\talpha\n')
disp(tr)
fprintf('%0.5f\t\t\t%0.5f\t\t%0.5f\t\t%0.5f\t\t%0.3f\n',UL,SUL,z,p,alpha)
disp(tr)
if p<alpha
    fprintf('\t\tThe survival functions are statistically different\n')
 else
    fprintf('\t\tThe survival functions are not statistically different\n')
end
return;
