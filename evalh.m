function [f] = evalh(n,x,ind)

global problem
	
%  ----------------------------------------------------------------------

%  AP1: Exemple 1 of "A modified Quasi-Newton method for vector optimization problem"

if ( strcmp(problem,'AP1' )  )
		
	if ( ind == 1 ) 
		f = 0.25d0 * ( ( x(1) - 1.0d0 ) ^ 4 + 2.0d0 * ( x(2) - 2.0d0 ) ^ 4 );
		
		return
	end
	
	if ( ind == 2 ) 
		f = exp( ( x(1) + x(2) ) / 2.0d0 ) + x(1) ^ 2 + x(2) ^ 2;
		
		return
	end	
	
	if ( ind == 3 ) 
		f = 1.0d0/6.0d0 * ( exp( - x(1) ) + 2.0d0 * exp( - x(2) ) );
		
		return
	end	
		
end		

%  ----------------------------------------------------------------------

%  AP2: Exemple 2 of "A modified Quasi-Newton method for vector optimization problem"

if ( strcmp(problem,'AP2' )  )
		
	if ( ind == 1 ) 
		f = x(1) ^ 2 - 4.0d0;
		
		return
	end
	
	if ( ind == 2 ) 
		f = ( x(1) - 1.0d0 ) ^ 2;
		
		return
	end	
		
end		

%  ----------------------------------------------------------------------

%  AP3: Exemple 3 of "A modified Quasi-Newton method for vector optimization problem"

if ( strcmp(problem,'AP3' )  )
		
	if ( ind == 1 ) 
		f = 0.25d0 * ( ( x(1) - 1.0d0 ) ^ 4 + 2.0d0 * ( x(2) - 2.0d0 ) ^ 4 );
		
		return
	end
	
	if ( ind == 2 ) 
		f = ( x(2) - x(1) ^ 2 ) ^ 2 + ( 1.0d0 - x(1) ) ^ 2;
		
		return
	end	
		
end	

%  ----------------------------------------------------------------------

%  AP4: Exemple 4 of "A modified Quasi-Newton method for vector optimization problem"

if ( strcmp(problem,'AP4' )  )
		
	if ( ind == 1 ) 
		f = 1.0d0/9.0d0 * ( ( x(1) - 1.0d0 ) ^ 4 + 2.0d0 * ( x(2) - 2.0d0 ) ^ 4 ...
		+ 3.0d0 * ( x(3) - 3.0d0 ) ^ 4 );
		
		return
	end
	
	if ( ind == 2 ) 
		f = exp( ( x(1) + x(2) + x(3) ) / 3.0d0 ) + x(1) ^ 2 + x(2) ^ 2 + x(3) ^ 2;
		
		return
	end	
	
	if ( ind == 3 ) 
		f = 1.0d0/1.2d1 * ( 3.0d0 * exp( -x(1) ) + 4.0d0 * exp( - x(2) ) + 3.0d0 * exp( -x(3) ) );
		
		return
	end	
		
end		

%  ----------------------------------------------------------------------	

%   BK1
%   A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit

if ( strcmp(problem,'BK1' )  )
	
	if ( ind == 1 ) 
		f = x(1) ^ 2 + x(2) ^ 2;
		
		return
	end
	
	if ( ind == 2 ) 
		f = ( x(1) - 5.0d0 ) ^ 2 + ( x(2) - 5.0d0 ) ^ 2;
		
		return	
	end
				
end			


%  ----------------------------------------------------------------------	

%   DD1

if ( strcmp(problem,'DD1' )  )
	
	if ( ind == 1 ) 
		f = x(1) ^ 2 + x(2) ^ 2 + x(3) ^ 2 + x(4) ^ 2 + x(5) ^ 2;
		
		return
	end
	
	if ( ind == 2 ) 
		f = 3.0d0 * x(1) + 2.0d0 * x(2) - x(3) / 3.0d0 + 1.0d-2 * ( x(4) - x(5) ) ^ 3;
		
		return	
	end
		
end

%  ----------------------------------------------------------------------

%   DGO1

if ( strcmp(problem,'DGO1' )  )
		
	if ( ind == 1 ) 
		f =  sin( x(1) );
		
		return
	end
	
	if ( ind == 2 ) 
		f =  sin( x(1) + 0.7d0 );
		
		return	
	end
end		

%  ----------------------------------------------------------------------

%   DGO2

if ( strcmp(problem,'DGO2' )  )
		
	if ( ind == 1 ) 
		f = x(1) ^ 2;
		
		return
	end
	
	if ( ind == 2 ) 
		f =  9.0d0 - sqrt( 8.1d1 - x(1) ^ 2 );
		
		return	
	end
end		

	

%  ----------------------------------------------------------------------

%   DTLZ1

if ( strcmp(problem,'DTLZ1' )  )

	%k = dimk;
	%m = dimm;
    
    k = 5;
    m = 3;

	faux = 0.0d0;
	for i = m:n
		faux = faux + ( x(i) - 0.5d0 ) ^ 2 - cos( 2.0d1 * pi * ( x(i) - 0.5d0 ) );
	end
	faux = 1.0d2 * ( k + faux );
	
	f = 0.5d0 * ( 1.0d0 + faux );				
	for i = 1:m-ind
		f = f * x(i);
	end
	
	if ( ind > 1 ) 
        f = f * ( 1.0d0 - x(m-ind+1) );
    end
			
	
	return	
	
end		

%  ----------------------------------------------------------------------

%   DTLZ2

if ( strcmp(problem,'DTLZ2' )  )

	%k = dimk;
	%m = dimm;
    
    k = 5;
    m = 3;

	faux = 0.0d0;
	for i = m:n
		faux = faux + ( x(i) - 0.5d0 ) ^ 2;
	end
	
	f = 1.0d0 + faux;			
	for i = 1:m-ind
		f = f * cos( x(i) * pi / 2.0d0 );
	end
	
	if ( ind > 1 ) 
        f = f * sin( x(m-ind+1) * pi / 2.0d0 );
    end
			
	
	return	
	
end

%  ----------------------------------------------------------------------

%   DTLZ3

if ( strcmp(problem,'DTLZ3' )  )

	%k = dimk;
	%m = dimm;
    
    k = 5;
    m = 3;

	faux = 0.0d0;
	for i = m:n
		faux = faux + ( x(i) - 0.5d0 ) ^ 2 - cos( 2.0d1 * pi * ( x(i) - 0.5d0 ) );
	end
	faux = 1.0d2 * ( k + faux );
	
	f = 1.0d0 + faux;					
	for i = 1:m-ind
		f = f * cos( x(i) * pi / 2.0d0 );
	end
	
	if ( ind > 1 ) 
        f = f * sin( x(m-ind+1) * pi / 2.0d0 );
    end
			
	
	return	
	
end	

%  ----------------------------------------------------------------------

%   DTLZ4

if ( strcmp(problem,'DTLZ4' )  )

	%k = dimk;
	%m = dimm;
    
    alpha = 2.0d0;
    
    k = 5;
    m = 3;

	faux = 0.0d0;
	for i = m:n
		faux = faux + ( x(i) - 0.5d0 ) ^ 2;
	end
	
	f = 1.0d0 + faux;		
	for i = 1:m-ind
		f = f * cos( x(i)^alpha * pi / 2.0d0 );
	end
	
	if ( ind > 1 ) 
        f = f * sin( x(m-ind+1)^alpha * pi / 2.0d0 );
    end
	
	return	
	
end

%  ----------------------------------------------------------------------

%   FA1

if ( strcmp(problem,'FA1' )  )
		
	if ( ind == 1 ) 
		f = ( 1.0d0 - exp( -4.0d0 * x(1) ) ) / ( 1.0d0 - exp( -4.0d0 ) );
		
		return
	end
	
	if ( ind == 2 ) 
		faux = ( 1.0d0 - exp( -4.0d0 * x(1) ) ) / ( 1.0d0 - exp( -4.0d0 ) );
		f = ( x(2) + 1.0d0 ) * ( 1.0d0 - ( faux / ( x(2) + 1.0d0  ) ) ^ 0.5d0 );
		
		return	
	end
	
	if ( ind == 3 ) 
		faux = ( 1.0d0 - exp( -4.0d0 * x(1) ) ) / ( 1.0d0 - exp( -4.0d0 ) );
		f = ( x(3) + 1.0d0 ) * ( 1.0d0 - ( faux / ( x(3) + 1.0d0  ) ) ^ 0.1d0 );
		
		return	
	end
end		

%  ----------------------------------------------------------------------

%   Far1

if ( strcmp(problem,'Far1' )  )
		
	if ( ind == 1 ) 
		f =  - 2.0d0 * exp( 1.5d1 * ( - ( x(1) - 0.1d0 ) ^ 2 - x(2) ^ 2 ) ) ...
			 - exp( 2.0d1 * ( - ( x(1) - 0.6d0 ) ^ 2 - ( x(2) - 0.6d0 ) ^ 2 ) ) ...
			 + exp( 2.0d1 * ( - ( x(1) + 0.6d0 ) ^ 2 - ( x(2) - 0.6d0 ) ^ 2 ) ) ...
			 + exp( 2.0d1 * ( - ( x(1) - 0.6d0 ) ^ 2 - ( x(2) + 0.6d0 ) ^ 2 ) ) ...
			 + exp( 2.0d1 * ( - ( x(1) + 0.6d0 ) ^ 2 - ( x(2) + 0.6d0 ) ^ 2 ) ) ;
		
		return
	end
	
	if ( ind == 2 ) 
		f =  2.0d0 * exp( 2.0d1 * ( - x(1) ^ 2 - x(2) ^ 2 ) ) ...
			 + exp( 2.0d1 * ( - ( x(1) - 0.4d0 ) ^ 2 - ( x(2) - 0.6d0 ) ^ 2 ) ) ...
			 - exp( 2.0d1 * ( - ( x(1) + 0.5d0 ) ^ 2 - ( x(2) - 0.7d0 ) ^ 2 ) ) ...
			 - exp( 2.0d1 * ( - ( x(1) - 0.5d0 ) ^ 2 - ( x(2) + 0.7d0 ) ^ 2 ) ) ...
			 + exp( 2.0d1 * ( - ( x(1) + 0.4d0 ) ^ 2 - ( x(2) + 0.8d0 ) ^ 2 ) ) ;
		
		return	
	end
		
end

%  ----------------------------------------------------------------------

%  FDS
%  NEWTONS METHOD FOR MULTIOBJECTIVE OPTIMIZATION
		
if ( strcmp(problem,'FDS' )  )
		
	if ( ind == 1 ) 
		f = 0.0d0;
		for i = 1:n
			f = f + i * ( x(i) - i ) ^ 4;
		end	
		f = f / ( n ^ 2 );
		
		return
	end
	
	if ( ind == 2 ) 
		f = exp( sum(x)/n ) + norm(x) ^ 2;
		
		return
	end	
	
	if ( ind == 3 ) 
		f = 0.0d0;
		for i = 1:n
			f = f + i * ( n - i + 1.0d0 ) * exp( - x(i) );
		end	
		f = f / ( n * ( n + 1.0d0 ) ) ;
		
		return
	end	
		
end	

%  ----------------------------------------------------------------------

%  FF1 
%  C. M. Fonseca and P. J. Fleming: ???An overview of evolutionary algorithms in multiobjective optimization

if ( strcmp(problem,'FF1' )  )

	if ( ind == 1 ) 
		f = 1.0d0 - exp( - ( x(1) - 1.0d0 ) ^ 2 - ( x(2) + 1.0d0 ) ^ 2 );
		
		return
	end
	
	if ( ind == 2 ) 
		f = 1.0d0 - exp( - ( x(1) + 1.0d0 ) ^ 2 - ( x(2) - 1.0d0 ) ^ 2 );
		
		return
	end	
		
end		

%  ----------------------------------------------------------------------

%   Hil1

if ( strcmp(problem,'Hil1' )  )
	a = 2.0d0 * pi / 3.6d2 * ( 4.5d1 + 4.0d1 * sin( 2.0d0 * pi * x(1) ) ...
		+ 2.5d1 * sin( 2.0d0 * pi * x(2) ) );
	b = 1.0d0 + 0.5d0 * cos( 2.0d0 * pi * x(1) );
	
	if ( ind == 1 ) 
		f = cos( a ) * b;
		
		return
	end
	
	if ( ind == 2 ) 
		f = sin( a ) * b;
		
		return	
	end
		
end			

%  ----------------------------------------------------------------------

%   IKK1

if ( strcmp(problem,'IKK1' )  )
	
	if ( ind == 1 ) 
		f = x(1) ^ 2;
		
		return
	end
	
	if ( ind == 2 ) 
		f = ( x(1) - 2.0d1 ) ^ 2;
		
		return	
	end
	
	if ( ind == 3 ) 
		f = x(2) ^ 2;
		
		return	
	end
		
end			

%  ----------------------------------------------------------------------

%   IM1

if ( strcmp(problem,'IM1' )  )
	
	if ( ind == 1 ) 
		f = 2.0d0 * sqrt( x(1) );
		
		return
	end
	
	if ( ind == 2 ) 
		f = x(1) * ( 1.0d0 - x(2) ) + 5.0d0;
		
		return	
	end
		
end

%  ----------------------------------------------------------------------	

%  JOS1 
%  Dynamic Weighted Aggregation for Evolutionary Multi-Objetive Optimization: Why fores It Work and How?

if ( strcmp(problem,'JOS1' )  )

	if ( ind == 1 ) 
		f = 0.0d0;
		for i = 1:n
				f = f + x(i) ^ 2;
		end
		f = f / n;
		
		return
	end
	
	if ( ind == 2 ) 
		f = 0.0d0;
		for i = 1:n
				f = f + ( x(i) - 2.0d0 ) ^ 2  ;
		end
		f = f / n;
		
		return
	end

end

%  ----------------------------------------------------------------------	

%  JOS4
%  Dynamic Weighted Aggregation for Evolutionary Multi-Objetive Optimization: Why fores It Work and How?

if ( strcmp(problem,'JOS4' )  )

	if ( ind == 1 ) 
	
		f = x(1);
		
		return
	end
	
	if ( ind == 2 ) 
	
		faux = 1.0d0 + 9.0d0 * sum(x(2:n)) / ( n - 1 );
		
		f = faux * ( 1.0d0 - ( x(1) / faux ) ^ 0.25d0 - ( x(1) / faux ) ^ 4.0d0 );
		
		
		return
	end

end	

%  ----------------------------------------------------------------------

%  	KW2

if ( strcmp(problem,'KW2' )  )
		
	if ( ind == 1 ) 
		f = - 3.0d0 * ( 1.0d0 - x(1) )^2 * exp( -x(1)^2 - ( x(2) + 1.0d0 ) ^ 2 ) ...
			+ 1.0d1 * ( x(1) / 5.0d0 - x(1)^3 - x(2)^5 ) * exp( - x(1)^2 - x(2)^2 ) ...
			+ 3.0d0 * exp( -( x(1) + 2.0d0 )^2 - x(2)^2 ) - 0.5d0 * ( 2.0d0 * x(1) + x(2) );
		
		return
	end
		
	if ( ind == 2 ) 
		f = - 3.0d0 * ( 1.0d0 + x(2) )^2 * exp( -x(2)^2 - ( 1.0d0 - x(1) ) ^ 2 ) ...
			+ 1.0d1 * ( - x(2) / 5.0d0 + x(2)^3 + x(1)^5 ) * exp( - x(1)^2 - x(2)^2 ) ...
			+ 3.0d0 * exp( -( 2.0d0 - x(2) )^2 - x(1)^2 ) ;
		
		return	
	end
		
end		

%  ----------------------------------------------------------------------

%  	LE1

if ( strcmp(problem,'LE1' )  )
										
	if ( ind == 1 ) 
		f = ( x(1) ^ 2 + x(2) ^ 2 ) ^ 0.125d0;
		
		return
	end
	
	if ( ind == 2 ) 
		f = ( ( x(1) - 0.5d0 ) ^ 2 + ( x(2) - 0.5d0 ) ^ 2 ) ^ 0.25d0;
		
		return
	end
		
end		

%  ----------------------------------------------------------------------

%  Lov1  

if ( strcmp(problem,'Lov1' )  )

	if ( ind == 1 ) 
		f = - ( -1.05d0 * x(1) ^ 2 - 0.98d0 * x(2) ^ 2 );
		
		return
	end
	
	if ( ind == 2 ) 
		f = - ( -0.99d0 * ( x(1) - 3.0d0 ) ^ 2 - 1.03d0 * ( x(2) - 2.5d0 ) ^ 2 );
		
		return
	end	
		
end		

%  ----------------------------------------------------------------------

%  Lov2

if ( strcmp(problem,'Lov2' )  )

	if ( ind == 1 ) 
		f = x(2);
		
		return
	end
	
	if ( ind == 2 ) 
		f = ( x(2) - x(1) ^ 3 ) / ( x(1) + 1.0d0 );
		f = - f;
		
		return
	end	
		
end		

%  ----------------------------------------------------------------------

%  Lov3  

if ( strcmp(problem,'Lov3' )  )
		
	if ( ind == 1 ) 
			f = - ( - x(1) ^ 2 - x(2) ^ 2 );
			
		return
	end
	
	if ( ind == 2 ) 
			f = - ( - ( x(1) - 6.0d0 ) ^ 2 + ( x(2) + 0.3d0 ) ^ 2 );
			
		return
	end	
	
end		

%  ----------------------------------------------------------------------

%  Lov4  

if ( strcmp(problem,'Lov4' )  )

	if ( ind == 1 ) 
		f = - x(1) ^ 2 - x(2) ^ 2 - 4.0d0 * ( exp( - ( x(1) + 2.0d0 ) ^ 2 - x(2) ^ 2 ) + ...
			exp( - ( x(1) - 2.0d0 ) ^ 2 - x(2) ^ 2 ) );
		f = - f;
		
		return
	end
	
	if ( ind == 2 ) 
		f = - ( x(1) - 6.0d0 ) ^ 2 - ( x(2) + 0.5d0 ) ^ 2;
		f = - f;
		
		return
	end	
		
end	

%  ----------------------------------------------------------------------

%  Lov5

if ( strcmp(problem,'Lov5' )  )

	MM = [ -1.0d0  , -0.03d0,  0.011d0; ...
	       -0.03d0 , -1.0d0 ,  0.07d0 ; ...
	        0.011d0,  0.07d0, -1.01d0 ];
	
	p = [ x(1)  ; x(2) - 0.15d0;  x(3)];
	a = 0.35d0;
	
	A1 = sqrt( 2.0d0 * pi / a ) * exp( dot( p, MM*p ) / a ^ 2 );
	
	p = [ x(1)  ; x(2) + 1.1d0 ;  0.5d0 * x(3)];
	a = 3.0d0;
	
	A2 = sqrt( 2.0d0 * pi / a ) * exp( dot( p, MM*p ) / a ^ 2 );
	
	faux = A1 + A2;

	if ( ind == 1 ) 
		f = sqrt(2.0d0)/2.0d0 * ( x(1) + faux );
		f = - f;
		
		return
	end
	
	if ( ind == 2 ) 
		f = sqrt(2.0d0)/2.0d0 * ( - x(1) + faux );
		f = - f;
		
		return
	end	
		
end		

%  ----------------------------------------------------------------------

%  Lov6

if ( strcmp(problem,'Lov6' )  )

	if ( ind == 1 ) 
		f = x(1);
		
	end
	
	if ( ind == 2 ) 
		f = 1.0d0 - sqrt( x(1) ) - x(1) * sin( 1.0d1 * pi * x(1) ) ...
		+ x(2) ^ 2 + x(3) ^ 2 + x(4) ^ 2 + x(5) ^ 2 + x(6) ^ 2;
		
		return
	end	
		
end	

%  ----------------------------------------------------------------------

%  	LTDZ
% 	Combining convergence and diversity in evolutionary multiobjective optimization

if ( strcmp(problem,'LTDZ' )  )
										
	if ( ind == 1 ) 
		f = 3.0d0 - ( 1.0d0 + x(3) ) * cos( x(1) * pi / 2.0d0 ) * cos( x(2) * pi / 2.0d0 );
		f = - f;
		
		return
	end
	
	if ( ind == 2 ) 
		f = 3.0d0 - ( 1.0d0 + x(3) ) * cos( x(1) * pi / 2.0d0 ) * sin( x(2) * pi / 2.0d0 );
		f = - f;
		
		return
	end
	
	if ( ind == 3 ) 
		f = 3.0d0 - ( 1.0d0 + x(3) ) * cos( x(1) * pi / 2.0d0 ) * sin( x(1) * pi / 2.0d0 );
		f = - f;
		
		return
	end
		
end						

%  ----------------------------------------------------------------------

%  	MGH9 

if ( strcmp(problem,'MGH9' )  )

	if ( ind == 1 || ind == 15 ) 
		y = 9.0d-4;
	elseif ( ind == 2 || ind == 14 )  		
		y = 4.4d-3;
	elseif ( ind == 3 || ind == 13 )  		
		y = 1.75d-2;
	elseif ( ind == 4 || ind == 12 )  		
		y = 5.4d-2;
	elseif ( ind == 5 || ind == 11 )  		
		y = 1.295d-1;		
	elseif ( ind == 6 || ind == 10 )  		
		y = 2.42d-1;
	elseif ( ind == 7 || ind == 9  )  		
		y = 3.521d-1;
	elseif ( ind == 8 )  		
		y = 3.989d-1;
	end
	
	t = ( 8.0d0 - ind ) / 2.0d0;
			
	f = x(1) * exp( - x(2) * ( t - x(3) ) ^ 2 / 2.0d0 ) - y	;		
	
	
	return
	
end	

%  ----------------------------------------------------------------------

%  	MGH16 

if ( strcmp(problem,'MGH16' )  )
	
	t = ind / 5.0d0;
			
	f = ( x(1) + t * x(2) - exp(t) ) ^ 2 + ( x(3) + x(4) * sin(t) - cos(t) ) ^ 2;
	
	
	return
		
end				

%  ----------------------------------------------------------------------

%  	MGH26 

if ( strcmp(problem,'MGH26' )  )
							
	t = 0.0d0;
	for i = 1:n
		t = t + cos(x(i));
	end
	
	f = ( n - t + ind * ( 1.0d0 - cos(x(ind)) ) - sin(x(ind)) ) ^ 2;
	
	
	return
		
end				

%  ----------------------------------------------------------------------

%  	MGH33

if ( strcmp(problem,'MGH33' )  )
							
	faux = 0.0d0;					
	for i = 1:n
		faux = faux + i * x(i);
	end
	
	f = ( ind * faux - 1.0d0 ) ^ 2	;		
	
	
	return
		
end	

%  ----------------------------------------------------------------------

%  	MHHM2

if ( strcmp(problem,'MHHM2' )  )
										
	if ( ind == 1 ) 
		f = ( x(1) - 0.8d0 ) ^ 2 + ( x(2) - 0.6d0 ) ^ 2;
		
		return
	end
	
	if ( ind == 2 ) 
		f = ( x(1) - 0.85d0 ) ^ 2 + ( x(2) - 0.7d0 ) ^ 2;
		
		return
	end
		
	if ( ind == 3 ) 
		f = ( x(1) - 0.9d0 ) ^ 2 + ( x(2) - 0.6d0 ) ^ 2;
		
		return
	end
			
end		

%  ----------------------------------------------------------------------

%   MLF1

if ( strcmp(problem,'MLF1' )  )
		
	if ( ind == 1 ) 
		f = ( 1.0d0 + x(1) / 2.0d1 ) * sin( x(1) );
		
		return
	end
	
	if ( ind == 2 ) 
		f = ( 1.0d0 + x(1) / 2.0d1 ) * cos( x(1) );
		
		return	
	end
		
end	

%  ----------------------------------------------------------------------

%   MLF2

if ( strcmp(problem,'MLF2' )  )
		
	if ( ind == 1 ) 
		f = 5.0d0 - ( ( x(1) ^ 2 + x(2) - 1.1d1 ) ^ 2 ...
			+ ( x(1) + x(2) ^ 2 - 7.0d0 ) ^ 2 ) / 2.0d2;
		f = - f;
		
		return
	end
	
	if ( ind == 2 ) 
		f = 5.0d0 - ( ( 4.0d0 * x(1) ^ 2 + 2.0d0 * x(2) - 1.1d1 ) ^ 2 ...
			+ ( 2.0d0 * x(1) + 4.0d0 *  x(2) ^ 2 - 7.0d0 ) ^ 2 ) / 2.0d2;
		f = - f;
		
		return	
	end
		
end	

%  ----------------------------------------------------------------------

%   MMR1 		
%   Box-constrained multi-objective optimization: A gradient-like method without ??????a priori?????? scalarization	

% 	if ( strcmp(problem,'MMR1' )  )
		
% 		if ( ind == 1 ) 
% 			f = 1.0d0 + x(1) ^ 2
% 			
% 			return
% 		end
	
% 		if ( ind == 2 ) 
% 			f = 2.0d0 - 0.8d0 * exp( - ( ( x(2) - 0.6d0 ) / 0.4d0 ) ^ 2 ) -...
% 			 exp( - ( ( x(2) - 0.2d0 ) / 0.04d0 ) ^ 2 )
% 			f = f / ( 1.0d0 + x(1) ^ 2 )
% 			
% 			return	
% 		end	
		
% 	end		

if ( strcmp(problem,'MMR1' )  )
		
	if ( ind == 1 ) 
		f = x(1);
		
		return
	end
	
	if ( ind == 2 ) 
		f = 2.0d0 - 0.8d0 * exp( - ( ( x(2) - 0.6d0 ) / 0.4d0 ) ^ 2 ) -...
		 exp( - ( ( x(2) - 0.2d0 ) / 0.04d0 ) ^ 2 );
		f = f / x(1);
		
		return	
	end	
		
end

%  ----------------------------------------------------------------------

%   MMR2
%   Box-constrained multi-objective optimization: A gradient-like method without ??????a priori?????? scalarization	

if ( strcmp(problem,'MMR2' )  )
		
	if ( ind == 1 ) 
		f = x(1);
		
		return
	end
	
	if ( ind == 2 ) 			
		faux = x(1) / ( 1.0d0 + 1.0d1 * x(2) );
		f = 1.0d0 - faux ^ 2 - faux * sin( 8.0d0 * pi * x(1) );
		f = f * ( 1.0d0 + 1.0d1 * x(2) );
		
		return	
	end	
		
end		
	
%  ----------------------------------------------------------------------

%   MMR3
%   Box-constrained multi-objective optimization: A gradient-like method without ??????a priori?????? scalarization	

if ( strcmp(problem,'MMR3' )  )
		
	if ( ind == 1 ) 
		f = x(1) ^ 3;
		
		return
	end
	
	if ( ind == 2 ) 
		f = ( x(2) - x(1) ) ^ 3;
		
		return	
	end	
		
end		

%  ----------------------------------------------------------------------

%   MMR4
%   Box-constrained multi-objective optimization: A gradient-like method without ??????a priori?????? scalarization	

if ( strcmp(problem,'MMR4' )  )
		
	if ( ind == 1 ) 
		f = x(1) - 2.0d0 * x(2) - x(3) - 3.6d1 / ( 2.0d0 * x(1) + x(2) + 2.0d0 * x(3) + 1.0d0 );
		
		return
	end
	
	if ( ind == 2 ) 
		f = - 3.0d0 * x(1) + x(2) - x(3);
		
		return	
	end	
		
end						

%  ----------------------------------------------------------------------

%   MOP 2	

if ( strcmp(problem,'MOP2' )  )
		
	if ( ind == 1 ) 
		f = 0.0d0;
		for i = 1:n
			f = f + ( x(i) - 1.0d0 / ( sqrt( real(n) ) ) ) ^ 2;
		end
		
		f = 1.0d0 - exp ( - f )	;
		
		return
	end
	
	if ( ind == 2 ) 
		f = 0.0d0;
		for i = 1:n
			f = f + ( x(i) + 1.0d0 / ( sqrt( real (n) ) ) ) ^ 2;
		end
		
		f = 1.0d0 - exp ( - f )	;
		
		return	
	end	
		
end	

%  ----------------------------------------------------------------------

%   MOP 3	

if ( strcmp(problem,'MOP3' )  )
		
	if ( ind == 1 ) 
		A1 = 0.5d0 * sin(1.0d0) - 2.0d0 * cos(1.0d0) + sin(2.0d0) - 1.5d0 * cos(2.0d0) ;
		A2 = 1.5d0 * sin(1.0d0) - cos(1.0d0) + 2.0d0 * sin(2.0d0) - 0.5d0 * cos(2.0d0);
		B1 = 0.5d0 * sin(x(1)) - 2.0d0 * cos(x(1)) + sin(x(2)) - 1.5d0 * cos(x(2)) ;
		B2 = 1.5d0 * sin(x(1)) - cos(x(1)) + 2.0d0 * sin(x(2)) - 0.5d0 * cos(x(2));	
		f = - 1.0d0 - ( A1 - B1 ) ^ 2 - ( A2 - B2 ) ^ 2;
		f = - f
		
		return
	end
	
	if ( ind == 2 ) 
		f = - ( x(1) + 3.0d0 ) ^ 2 - ( x(2) + 1.0d0 ) ^ 2;
		f = - f;
		
		return	
	end	
		
end			

%  ----------------------------------------------------------------------

%   MOP 5	

if ( strcmp(problem,'MOP5' )  )
		
	if ( ind == 1 ) 
		f = 0.5d0 * ( x(1) ^ 2 + x(2) ^ 2 ) + sin( x(1) ^ 2 + x(2) ^ 2 );
		
		return
	end
	
	if ( ind == 2 ) 
		f = ( 3.0d0 * x(1) - 2.0d0 * x(2) + 4.0d0 ) ^ 2 / 8.0d0 + ...
			( x(1) - x(2) + 1.0d0 ) ^ 2 / 2.7d1 + 1.5d1;
		
		return	
	end
	
	if ( ind == 3 ) 
		f = 1.0d0 / ( x(1) ^ 2 + x(2) ^ 2 + 1.0d0 ) - 1.1d0 * exp( - x(1) ^ 2 - x(2) ^ 2 );
		
		return
	end	
		
end	

%  ----------------------------------------------------------------------

%   MOP6

if ( strcmp(problem,'MOP6' )  )
		
	if ( ind == 1 ) 
		f = x(1);
		
		return
	end
	
	if ( ind == 2 ) 
		a = 1.0d0 + 1.0d1 * x(2) ;
		t = x(1) / a;
		f = a * ( 1.0d0 - t ^ 2 - t * sin( 8.0d0 * pi * x(1) ) );
		
		return	
	end
		
end	


%  ----------------------------------------------------------------------

%   MOP 7

if ( strcmp(problem,'MOP7' )  )
		
	if ( ind == 1 ) 
		f =  ( x(1) - 2.0d0 ) ^ 2 / 2.0d0 ...
		+ ( x(2) + 1.0d0 ) ^ 2 / 1.3d1 + 3.0d0;
		
		return
	end
	
	if ( ind == 2 ) 
		f =  ( x(1) + x(2) - 3.0d0 ) ^ 2 / 3.6d1 ...
		+ ( - x(1) + x(2) + 2.0d0 ) ^ 2 / 8.0d0 - 1.7d1;
		
		return	
	end
	
	if ( ind == 3 ) 
		f =  ( x(1) + 2.0d0 * x(2) - 1.0d0 ) ^ 2 / 1.75d2 ...
		+ ( - x(1) + 2.0d0 * x(2) ) ^ 2 / 1.7d1 - 1.3d1;
		
		return
	end	
			
end			

%  ----------------------------------------------------------------------

%  	PNR

if ( strcmp(problem,'PNR' )  )
										
	if ( ind == 1 ) 
		f = x(1) ^ 4 + x(2) ^ 4 - x(1) ^ 2 + x(2) ^ 2 - 1.0d1 * x(1) * x(2) + 2.0d1;
		
		return
	end
	
	if ( ind == 2 ) 
		f = x(1) ^ 2 + x(2) ^ 2;
		
		return
	end
		
end

%  ----------------------------------------------------------------------

%   QV1

if ( strcmp(problem,'QV1' )  )
		
	if ( ind == 1 ) 
		f = 0.0d0;
		for i = 1:n
			f = f + x(i) ^ 2 - 1.0d1 * cos( 2.0d0 * pi * x(i) ) + 1.0d1;
		end
	
		f = ( f / n ) ^ 0.25 ;
		
		return
	end
	
	if ( ind == 2 ) 
		f = 0.0d0;
		for i = 1:n
			f = f + ( x(i) - 1.5d0 ) ^ 2 - 1.0d1 * cos( 2.0d0 * pi * ( x(i) -1.5d0 ) ) + 1.0d1;
		end
	
		f = ( f / n ) ^ 0.25;
		
		return	
	end
		
end	

%  ----------------------------------------------------------------------

%  SD

if ( strcmp(problem,'SD' )  )
		
	if ( ind == 1 ) 
		f = 2.0d0 * x(1) + sqrt(2.0d0) * ( x(2) + x(3) ) + x(4);
		
		return
	end
	
	if ( ind == 2 ) 
		f = 2.0d0 / x(1) + 2.0d0 * sqrt(2.0d0) / x(2) + 2.0d0 * sqrt(2.0d0) / x(3) + 2.0d0 / x(4);
		
		return		
	end	
		
end		

%  ----------------------------------------------------------------------

%  SLCDT1

if ( strcmp(problem,'SLCDT1' )  )
		
	if ( ind == 1 ) 
		f = 0.5d0 * ( sqrt( 1.0d0 + ( x(1) + x(2) ) ^ 2 ) + ...
				sqrt( 1.0d0 + ( x(1) - x(2) ) ^ 2 ) + x(1) - x(2) ) +  ...
				0.85d0 * exp( - ( x(1) + x(2) ) ^ 2 )	;	
		
		return
	end
	
	if ( ind == 2 ) 
		f = 0.5d0 * ( sqrt( 1.0d0 + ( x(1) + x(2) ) ^ 2 ) + ...
				sqrt( 1.0d0 + ( x(1) - x(2) ) ^ 2 ) - x(1) + x(2) ) +  ...
				0.85d0 * exp( - ( x(1) + x(2) ) ^ 2 )	;
		
		return		
	end	
		
end	

%  ----------------------------------------------------------------------

%  	SLCDT2
% 	Convergence of stochastic search algorithms to finite size pareto set approximations 

if ( strcmp(problem,'SLCDT2' )  )
									
	if ( ind == 1 ) 
		f = ( x(1) - 1.0d0 ) ^ 4;
		for i = 2:n
			f = f + ( x(i) - 1.0d0 ) ^ 2;
		end
		
		return
	end
	
	if ( ind == 2 ) 
		f = ( x(2) + 1.0d0 ) ^ 4;
		for i = 1:n
			if ( i ~= 2 ) 
                f = f + ( x(i) + 1.0d0 ) ^ 2;
            end
		end
		
		return
	end
	
	if ( ind == 3 ) 
		f = ( x(3) - 1.0d0 ) ^ 4;
		for i = 1:n
			if ( i ~= 3 ) 
                f = f + ( x(i) - ( - 1.0d0 ) ^ (i+1) ) ^ 2;
            end
		end
		
		return
	end
		
end		
					
%  ----------------------------------------------------------------------

%   SP1

if ( strcmp(problem,'SP1' )  )
		
	if ( ind == 1 ) 
		f = ( x(1) - 1.0d0 ) ^ 2 + ( x(1) - x(2) ) ^ 2;
		
		return
	end
	
	if ( ind == 2 ) 
		f = ( x(2) - 3.0d0 ) ^ 2 + ( x(1) - x(2) ) ^ 2;
		
		return	
	end
		
end				
	
%  ----------------------------------------------------------------------

%   SSFYY2

if ( strcmp(problem,'SSFYY2' )  )
	
	if ( ind == 1 ) 
		f = 1.0d1 + x(1) ^ 2 - 1.0d1 * cos( x(1) * pi / 2.0d0 );
		
		return
	end
	
	if ( ind == 2 ) 
		f = ( x(1) - 4.0d0 ) ^ 2;
		
		return	
	end
		
end			
			
%  ----------------------------------------------------------------------

%   SK1

if ( strcmp(problem,'SK1' )  )
		
	if ( ind == 1 ) 
		f = - x(1) ^ 4 - 3.0d0 * x(1) ^ 3 + 1.0d1 * x(1) ^ 2 + 1.0d1 * x(1) + 1.0d1;
		f = - f;
		
		return
	end
	
	if ( ind == 2 ) 
		f = - 0.5d0 * x(1) ^ 4 + 2.0d0 * x(1) ^ 3 + 1.0d1 * x(1) ^ 2 - 1.0d1 * x(1) + 5.0d0;
		f = - f;
		
		return	
	end
		
end	

%  ----------------------------------------------------------------------

%   SK2

if ( strcmp(problem,'SK2' )  )
		
	if ( ind == 1 ) 
		f = - ( x(1) - 2.0d0 ) ^ 2 - ( x(2) + 3.0d0 ) ^ 2 ...
		- ( x(3) - 5.0d0 ) ^ 2 - ( x(4) - 4.0d0 ) ^ 2 + 5.0d0;
		f = - f;
		
		return
	end
	
	if ( ind == 2 ) 
		f = ( sin( x(1) ) + sin( x(2) ) + sin( x(3) ) + sin( x(4) ) ) / ...
		( 1.0d0 + ( x(1) ^ 2 + x(2) ^ 2 + x(3) ^ 2 + x(4) ^ 2 ) / 1.0d2 );
		f = - f;
		
		return	
	end
		
end		


%  ----------------------------------------------------------------------

%   TKLY1

if ( strcmp(problem,'TKLY1' )  )
		
	if ( ind == 1 ) 
		f = x(1);
		
		return
	end
	
	if ( ind == 2 ) 
	
		A1 = ( 2.0d0 - exp( - ( ( x(2) - 0.1d0 ) / 4.0d-3 ) ^ 2 ) ...
			         - 0.8d0 * exp( - ( ( x(2) - 0.9d0 ) / 4.0d-1 ) ^ 2 ) );
		A2 = ( 2.0d0 - exp( - ( ( x(3) - 0.1d0 ) / 4.0d-3 ) ^ 2 ) ...
			         - 0.8d0 * exp( - ( ( x(3) - 0.9d0 ) / 4.0d-1 ) ^ 2 ) );
		A3 = ( 2.0d0 - exp( - ( ( x(4) - 0.1d0 ) / 4.0d-3 ) ^ 2 ) ...
			         - 0.8d0 * exp( - ( ( x(4) - 0.9d0 ) / 4.0d-1 ) ^ 2 ) );
		f = A1 * A2 * A3 / x(1);
		
		return	
	end
		
end	


%  ----------------------------------------------------------------------

%  	Toi4

if ( strcmp(problem,'Toi4' )  )
		
	if ( ind == 1 ) 
		f = x(1) ^ 2 + x(2) ^ 2 + 1.0d0;
		
		return
	end
	
	if ( ind == 2 ) 
		f = 0.5d0 * ( ( x(1) - x(2) ) ^ 2 + ( x(3) - x(4) ) ^ 2  ) + 1.0d0;
		
		return	
	end
		
end					

%  ----------------------------------------------------------------------

%  	Toi8

if ( strcmp(problem,'Toi8' )  )
		
	if ( ind == 1 ) 
		f = ( 2.0d0 * x(1) - 1.0d0 ) ^ 2 ;
		
		return
	end
	
	if ( ind ~= 1 ) 
		f = ind * ( 2.0d0 * x(ind-1) - x(ind) ) ^ 2;
		
		return	
	end
		
end		

%  ----------------------------------------------------------------------

%  	Toi9

if ( strcmp(problem,'Toi9' )  )
		
	if ( ind == 1 ) 
		f = ( 2.0d0 * x(1) - 1.0d0 ) ^ 2 + x(2) ^ 2;
		
		return
	end
	
	if ( ind > 1 && ind < n ) 
		f = ind * ( 2.0d0 * x(ind-1) - x(ind) ) ^ 2  ...
		   - ( ind - 1.0d0 ) * x(ind-1)^2 + ind * x(ind) ^ 2;
		
		return	
	end
	
	if ( ind == n ) 
		f = n * ( 2.0d0 * x(n-1) - x(n) ) ^ 2 - ( n - 1.0d0 ) * x(n-1)^2;
		
		return
	end
		
end		

%  ----------------------------------------------------------------------

%  	Toi10 (Rosenbrock)

if ( strcmp(problem,'Toi10' )  )
		
	f = 1.0d2 * ( x(ind+1) - x(ind) ^ 2 ) ^ 2 + ( x(ind+1) - 1.0d0 ) ^ 2;
	
		return
		
end	

%  ----------------------------------------------------------------------

%   VU1

if ( strcmp(problem,'VU1' )  )
		
	if ( ind == 1 ) 
		f = 1.0d0 / ( x(1) ^ 2 + x(2) ^ 2 + 1.0d0 );
		
		return
	end
	
	if ( ind == 2 ) 
		f = x(1) ^ 2 + 3.0d0 * x(2) ^ 2 + 1.0d0;
		
		return	
	end
		
end		

%  ----------------------------------------------------------------------

%   VU2

if ( strcmp(problem,'VU2' )  )
		
	if ( ind == 1 ) 
		f = x(1) + x(2) + 1.0d0;
		
		return
	end
	
	if ( ind == 2 ) 
		f = x(1) ^ 2 + 2.0d0 * x(2) - 1.0d0;
		
		return	
	end
		
end	

%  ----------------------------------------------------------------------

%   ZDT1

if ( strcmp(problem,'ZDT1' )  )
		
	if ( ind == 1 ) 
		f = x(1);
		
		return
	end
	
	if ( ind == 2 ) 
		faux = 1.0d0 + 9.0d0 * sum(x(2:n)) / ( n - 1 );
		
		f = faux * ( 1.0d0 - sqrt( x(1) / faux ) );
		
		return	
	end
		
end			

%  ----------------------------------------------------------------------

%   ZDT2

if ( strcmp(problem,'ZDT2' )  )
		
	if ( ind == 1 ) 
		f = x(1);
		
		return
	end
	
	if ( ind == 2 ) 
		faux = 1.0d0 + 9.0d0 * sum(x(2:n)) / ( n - 1 );
		
		f = faux * ( 1.0d0 - ( x(1) / faux ) ^ 2 );
		
		return	
	end
		
end				

%  ----------------------------------------------------------------------

%   ZDT3

if ( strcmp(problem,'ZDT3' )  )
		
	if ( ind == 1 ) 
		f = x(1);
		
		return
	end
	
	if ( ind == 2 ) 
		faux = 1.0d0 + 9.0d0 * sum(x(2:n)) / ( n - 1 );
		t = x(1) / faux;
		
		f = faux * ( 1.0d0 - sqrt( t ) - t * sin( 1.0d1 * pi * x(1) ) );
		
		return	
	end
		
end			
	
%  ----------------------------------------------------------------------

%   ZDT4

if ( strcmp(problem,'ZDT4' )  )
		
	if ( ind == 1 ) 
		f = x(1);
		
		return
	end
	
	if ( ind == 2 ) 
		faux = 0.0d0;
		for i = 2:n
			faux = faux + x(i) ^ 2 - 1.0d1 * cos( 4.0d0 * pi * x(i) );
		end
		faux = faux + 1.0d0 + 1.0d1 * ( n - 1 );
		t = x(1) / faux;
		
		f = faux * ( 1.0d0 - sqrt( t ) );
		
		return	
	end
		
end		

%  ----------------------------------------------------------------------

%   ZDT6

if ( strcmp(problem,'ZDT6' )  )
		
	if ( ind == 1 ) 
		f = 1.0d0 - exp( -4.0d0 * x(1) ) * ( sin( 6.0d0 * pi * x(1) ) ) ^ 6;
		
		return
	end
	
	if ( ind == 2 ) 
		faux = 1.0d0 + 9.0d0 * ( sum(x(2:n)) / ( n - 1 ) ) ^ 0.25d0;
			
			f = faux * ( 1.0d0 - ( ( 1.0d0 - exp( -4.0d0 * x(1) ) * ...
			    ( sin( 6.0d0 * pi * x(1) ) ) ^ 6 ) / faux ) ^ 2 );
			
			return	
		end
			
	end			
	
%  ----------------------------------------------------------------------

	%   ZLT1

	if ( strcmp(problem,'ZLT1' )  )
		
		f = ( x(ind) - 1.0d0 ) ^ 2;
		for i = 1:n
			if ( i ~= ind ) 
                f = f + x(i) ^ 2;
            end
		end
		
		return
			
	end				
	