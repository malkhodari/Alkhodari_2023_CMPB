function f_hrv=ada_f(hrv)               % f_hrv is filtered hrv by adaptive filtering method 

%    Function: Use adaptive method to remove abnormal RR intervals, especially the patients have VPC(Ventricular premature complexes)
%    Reference: 1. ' Renormalised Entropy: A New Method of Non-Linear Dynamics for the Analysis of Heart Rate Variability',N Wessel, A Voss, Computers in Cardiology 1994
%               2. ' Nonlinear analysis of complex phenomena in cardiological data, N Wessel, A Voss, Herzschr Elektrophys 11:159-173(2000)'
%
%    Steps: 1. The adaptive percent-filter: (1) Use tachogram X(n) to calculate t(n),Ua(n),Ca(n),;a(n)    ( ; is lampda )
%   											  t(n):  the binomial filtered series
%   										      Ua(n): the adaptive mean
%   											  Ca(n): the adaptive standard deviation
%   											  ;a(n): the adaptive second moment
%    										  
%    											  t(n)=( X(n-3)+6X(n-2)+15X(n-1)+20X(n)+15X(n+1)+6X(n+2)+X(n+3) )/64
%                                                 Ua(n)=Ua(n-1)-c*( Ua(n-1)-t(n-1) ) ;  c: cotrolling coeffient(here is 0.05)
%                                                 Ca(n)=sqrt(Ua(n)^2-;a(n)) ; 
%    											  ;a(n)=;a(n-1)-c*(;a(n-1)-t(n-1)^2);  
%
%                                           (2) The RR-interval X(n) is classified as abnormal if the following inequalities satisfied;
%										 		  |X(n)-X(n-1)|>p*X(n-1)/100+Cf*Ca_mean ;  p: propotional limit(here is 10), Cf: 3.0,  Ca_mean: averaged Ca(n)
%												  |X(n)-Xlv|>p*Xlv/100+Cf*Ca_mean
%                                                 
%                                                 Values recognised as not normal are replaced with a random number from [Ua(n)-0.5*Ca(n),Ua(n)+0.5*Ca(n)]
%                                                 to avoid false decreased variabilities.
% 
%          2. The adaptive controlling filter: From the resulting time series X(n)% of the adaptive percent-filter, the binomial filtered series 
%             and the respective adaptive mean value and standard deviation are again calculated.
% 
%             The value X(n)% is considered to be not normal if |X(n)%-Ua(n)|>Cf1*Ca(n)+Cb;  
%             Cf1: filter coefficient(here Cf1=3.0); Cb is a basic variability(for HRV here Cb=20ms)
% 
%             Abnormal values are replaced with the respective values of the binomial filtered series
% 

X=hrv;
n=length(X);

% Calculate t(n) %
% When i =1,2,3 or i=n-2,n-1,n, it can't have i-3,i-2,i-1,i+3,i+2,i+1. So use t(i)to replace t(i-3),t(i-2),i.e. So does the other situations.
% For exmaple, when i=1, it can't have i-3,i-2,i-1, so use t(1) to replace the first 3 items of the equation ***

for N=1:2

    
   t(1)=( X(1)+6*X(1)+15*X(1)+20*X(1)+15*X(1+1)+6*X(1+2)+X(1+3) )/64;                                                                  
   t(2)=( X(2)+6*X(2)+15*X(2-1)+20*X(2)+15*X(2+1)+6*X(2+2)+X(2+3) )/64;
   t(3)=( X(3)+6*X(3-1)+15*X(3-2)+20*X(3)+15*X(3+1)+6*X(3+2)+X(3+3) )/64;

   for i=4:n-3
       t(i)=( X(i-3)+6*X(i-2)+15*X(i-1)+20*X(i)+15*X(i+1)+6*X(i+2)+X(i+3) )/64;    
   end

   t(n-2)=( X(n-2-3)+6*X(n-2-2)+15*X(n-2-1)+20*X(n-2)+15*X(n-2+1)+6*X(n-2+2)+X(n-2) )/64;
   t(n-1)=( X(n-1-3)+6*X(n-1-2)+15*X(n-1-1)+20*X(n-1)+15*X(n-1+1)+6*X(n-1)+X(n-1) )/64;
   t(n)=( X(n-3)+6*X(n-2)+15*X(n-1)+20*X(n)+15*X(n)+6*X(n)+X(n) )/64;

   % Calculate <a(n), in this code, use Ua(n) to replace <a(n), lampda for ;, sigma for C, cause matlab can't recognize them. 

   Ua(1)=mean(X);lampda(1)=Ua(1)^2; 

   c=0.05;  % cotrolling coeffient
  
   % lampda(i)=

   for i=2:n
       Ua(i)=Ua(i-1)-c*( Ua(i-1)-t(i-1) );  
	   lampda(i)=lampda(i-1)-c*( lampda(i-1)-t(i-1)^2 );
	   
   end
   
   for i=1:n
       sigma(i)=abs(sqrt( Ua(i)^2-lampda(i) ));
   end
   
   % Classify normal and abmormal RR-intervals 

   Cf=3;p=10;
   sigma_mean=mean(sigma);
   last=Ua(1);
   for i=2:n
	   if abs(X(i)-X(i-1))>p*X(i-1)/100+Cf*sigma_mean & abs(X(i)-last)>p*last/100+Cf*sigma_mean 
	      X(i)=rand(1)*sigma(i)+(Ua(i)-0.5*sigma(i));  % produce a random number between Ua(i)-0.5*sigma(i) and Ua(i)+0.5*sigma(i)
      else
          last=X(i); 
      end
   end
% After the above loop, the X(n) is changed into X(n)%, but in the program X(n)% is still denoted by X(n)
% and use X(n)% to recalculate t(n), Ua, lampda and sigma, so I use 2 times loop. 'for N=1:2'
end

 
% The adaptive controlling procedure %

Cf1=3; sigmaB=0.02;   % sigmaB is 200ms %

f_hrv=X;

for i=1:n
	if abs(X(i)-Ua(i))>(Cf1*sigma(i)+sigmaB)
		f_hrv(i)=t(i);
	end
end


