function [g] = evalgradh(n,x,ind)
	
global problem
	
% ----------------------------------------------------------------------

% AP1: Exemple 1 of "A modified Quasi-Newton method for vector optimization problem"

if ( strcmp(problem,'AP1' )  )
		
	if ( ind == 1 ) 
		g(1) = ( x(1) - 1.0d0 ) ^ 3 ;
		g(2) = 2.0d0 * ( x(2) - 2.0d0 ) ^ 3;
		
		return
	end
	
	if ( ind == 2 ) 
		g(1) = 0.5d0 * exp( ( x(1) + x(2) ) / 2.0d0 ) + 2.0d0 * x(1);
		g(2) = 0.5d0 * exp( ( x(1) + x(2) ) / 2.0d0 ) + 2.0d0 * x(2);
		
		return
	end	
	
	if ( ind == 3 ) 
		g(1) = - 1.0d0/6.0d0 * exp( - x(1) ) ;
		g(2) = - 1.0d0/3.0d0 * exp( - x(2) ) ;
		
		return
	end	
	
end	

% ----------------------------------------------------------------------

% AP2: Exemple 2 of "A modified Quasi-Newton method for vector optimization problem"

if ( strcmp(problem,'AP2' )  )
		
	if ( ind == 1 ) 
		g(1) = 2.0d0 * x(1);
		
		return
	end
	
	if ( ind == 2 ) 
		g(1) = 2.0d0 * ( x(1) - 1.0d0 );
		
		return
	end	
		
end		


% ----------------------------------------------------------------------

% AP3: Exemple 3 of "A modified Quasi-Newton method for vector optimization problem"

if ( strcmp(problem,'AP3' )  )
		
	if ( ind == 1 ) 
		g(1) = ( x(1) - 1.0d0 ) ^ 3 ;
		g(2) = 2.0d0 * ( x(2) - 2.0d0 ) ^ 3;
		
		return
	end
	
	if ( ind == 2 ) 
		g(1) = - 4.0d0 * x(1) * ( x(2) - x(1) ^ 2 ) - 2.0d0 * ( 1.0d0 - x(1) );
		g(2) = 2.0d0 * ( x(2) - x(1) ^ 2 );
		
		return
	end	
		
end	

% ----------------------------------------------------------------------

% AP4: Exemple 4 of "A modified Quasi-Newton method for vector optimization problem"

if ( strcmp(problem,'AP4' )  )
	
	if ( ind == 1 ) 
		g(1) = 4.0d0/9.0d0 * ( x(1) - 1.0d0 ) ^ 3 ;
		g(2) = 8.0d0/9.0d0 * ( x(2) - 2.0d0 ) ^ 3 ;
		g(3) = 1.2d1/9.0d0 * ( x(3) - 3.0d0 ) ^ 3 ;
		
		return
	end
	
	if ( ind == 2 ) 
		g(1) = 1.0d0/3.0d0 * exp( ( x(1) + x(2) + x(3) ) / 3.0d0 ) + 2.0d0 * x(1) ;
		g(2) = 1.0d0/3.0d0 * exp( ( x(1) + x(2) + x(3) ) / 3.0d0 ) + 2.0d0 * x(2) ;
		g(3) = 1.0d0/3.0d0 * exp( ( x(1) + x(2) + x(3) ) / 3.0d0 ) + 2.0d0 * x(3) ;
		
		return
	end	
	
	if ( ind == 3 ) 
		g(1) = - 1.0d0/4.0d0 * exp( -x(1) ) ;
		g(2) = - 1.0d0/3.0d0 * exp( -x(2) ) ;
		g(3) = - 1.0d0/4.0d0 * exp( -x(3) ) ;
		
		return
	end	
		
end	

% ----------------------------------------------------------------------	

%  BK1

if ( strcmp(problem,'BK1' )  )
		
	if ( ind == 1 ) 
		for i = 1:n
			g(i) = 2.0d0 * x(i);
		end
		
		return
	end
	
	if ( ind == 2 ) 
		for i = 1:n
			g(i) = 2.0d0 * ( x(i) - 5.0d0 ) ;
		end
		
		return	
	end
		
end		

% ----------------------------------------------------------------------	

%  DD1

if ( strcmp(problem,'DD1' )  )
		
	if ( ind == 1 ) 
		for i = 1:n
			g(i) = 2.0d0 * x(i);
		end
		
		return
	end
	
	if ( ind == 2 ) 
		g(1) = 3.0d0;
		g(2) = 2.0d0;
		g(3) = - 1.0d0 / 3.0d0;
		g(4) = 3.0d-2 * ( x(4) - x(5) ) ^ 2;
		g(5) = - 3.0d-2 * ( x(4) - x(5) ) ^ 2;
		
		return	
	end
		
end		

% ----------------------------------------------------------------------

%  DGO1

if ( strcmp(problem,'DGO1' )  )
		
	if ( ind == 1 ) 
		g(1) =  cos( x(1) );
		
		return
	end
	
	if ( ind == 2 ) 
		g(1) =  cos( x(1) + 0.7d0 );
		
		return	
	end
end		

% ----------------------------------------------------------------------

%  DGO2

if ( strcmp(problem,'DGO2' )  )
		
	if ( ind == 1 ) 
		g(1) =  2.0d0 * x(1);
		
		return
	end
	
	if ( ind == 2 ) 
		g(1) =  x(1) / sqrt( 8.1d1 - x(1) ^ 2 );
		
		return	
	end
end		

% ----------------------------------------------------------------------

%  DTLZ1

if ( strcmp(problem,'DTLZ1' )  )

%		k = dimk
%		m = dimm

%		faux = 0.0d0
%		for i = m:n
%			faux = faux + ( x(i) - 0.5d0 ) ^ 2 - cos( 2.0d1 * pi * ( x(i) - 0.5d0 ) )
%		end
%		faux = 1.0d2 * ( k + faux )
	
%		f = 0.5d0 * ( 1.0d0 + faux )			
%		for i = 1:m-ind
%			f = f * x(i)
%		end
	
%		if ( ind > 1 ) f = f * ( 1.0d0 - x(m-ind+1) )
	
	
	%k = dimk;
	%m = dimm;
    
    k = 5;
    m = 3;
    
	
	for i = 1:n
		g(i) = 0.0d0;
	end

	faux = 0.0d0;
	for i = m:n
		faux = faux + ( x(i) - 0.5d0 ) ^ 2 - cos( 2.0d1 * pi * ( x(i) - 0.5d0 ) );
	end
	faux = 1.0d2 * ( k + faux );
	
	for i = 1:m-ind
		g(i) = 0.5d0 * ( 1.0d0 + faux );
	end 
	
	for i = 1:m-ind
		for j = 1:m-ind
			if ( j == i ) 
                continue
            end
			g(i) = g(i) * x(j);
		end
		if ( ind > 1 ) 
            g(i) = g(i) * ( 1.0d0 - x(m-ind+1) );
        end
	end
	
	if ( ind > 1 ) 
		g(m-ind+1) = - 0.5d0 * ( 1.0d0 + faux );
		for j = 1:m-ind
			g(m-ind+1) = g(m-ind+1) * x(j);
		end
	end
			
	
	faux = 0.5d0;
	for i = 1:m-ind
		faux = faux * x(i);
	end

	if ( ind > 1 ) 
        faux = faux * ( 1.0d0 - x(m-ind+1) );
    end
	
	for i = m:n
		g(i) = faux * 1.0d2 * ( 2.0d0 * ( x(i) - 0.5d0 ) + 2.0d1 * pi * sin( 2.0d1 * pi * ( x(i) - 0.5d0 ) ) );
	end
				
	
	return

end

% ----------------------------------------------------------------------

%  DTLZ2

if ( strcmp(problem,'DTLZ2' )  )

%		k = dimk
%		m = dimm

%		faux = 0.0d0
%		for i = m:n
%			faux = faux + ( x(i) - 0.5d0 ) ^ 2
%		end
	
%		f = 1.0d0 + faux					
%		for i = 1:m-ind
%			f = f * cos( x(i) * pi / 2.0d0 )
%		end
	
%		if ( ind > 1 ) f = f * sin( x(m-ind+1) * pi / 2.0d0 )
	
	%k = dimk;
	%m = dimm;
    
    k = 5;
    m = 3;
	
	for i = 1:n
		g(i) = 0.0d0;
	end

	faux = 0.0d0;
	for i = m:n
		faux = faux + ( x(i) - 0.5d0 ) ^ 2;
	end
	faux = 1.0d0 + faux;
	
	for i = 1:m-ind
		g(i) = faux;
	end 
	
	for i = 1:m-ind
		for j = 1:m-ind
			if ( j == i ) 
				g(i) = - g(i) * pi / 2.0d0 * sin( x(i) * pi / 2.0d0 );
			else
				g(i) = g(i) * cos( x(j) * pi / 2.0d0 );
			end
		end
		if ( ind > 1 ) 
            g(i) = g(i) * sin( x(m-ind+1) * pi / 2.0d0 );
        end
	end
	
	if ( ind > 1 ) 
		g(m-ind+1) = faux * pi / 2.0d0 * cos( x(m-ind+1) * pi / 2.0d0 );
		for j = 1:m-ind
			g(m-ind+1) = g(m-ind+1) * cos( x(j) * pi / 2.0d0 );
		end
	end
			
	
	faux = 1.0d0;
	for i = 1:m-ind
		faux = faux * cos( x(i) * pi / 2.0d0 );
	end

	if ( ind > 1 ) 
        faux = faux * sin( x(m-ind+1) * pi / 2.0d0 );
    end
	
	for i = m:n
		g(i) = faux * 2.0d0 * ( x(i) - 0.5d0 );
	end
				
	
	return

end

% ----------------------------------------------------------------------

%  DTLZ3

if ( strcmp(problem,'DTLZ3' )  )

%		k = dimk
%		m = dimm

%		faux = 0.0d0
%		for i = m:n
%			faux = faux + ( x(i) - 0.5d0 ) ^ 2 - cos( 2.0d1 * pi * ( x(i) - 0.5d0 ) )
%		end
%		faux = 1.0d2 * ( k + faux )
	
%		f = 1.0d0 + faux					
%		for i = 1:m-ind
%			f = f * cos( x(i) * pi / 2.0d0 )
%		end
	
%		if ( ind > 1 ) f = f * sin( x(m-ind+1) * pi / 2.0d0 )
	
	%k = dimk;
	%m = dimm;
    
    k = 5;
    m = 3;
	
	for i = 1:n
		g(i) = 0.0d0;
	end

	faux = 0.0d0;
	for i = m:n
		faux = faux + ( x(i) - 0.5d0 ) ^ 2 - cos( 2.0d1 * pi * ( x(i) - 0.5d0 ) );
	end
	faux = 1.0d2 * ( k + faux );
	faux = faux + 1.0d0;
	
	for i = 1:m-ind
		g(i) = faux;
	end 
	
	for i = 1:m-ind
		for j = 1:m-ind
			if ( j == i ) 
				g(i) = - g(i) * pi / 2.0d0 * sin( x(i) * pi / 2.0d0 );
			else
				g(i) = g(i) * cos( x(j) * pi / 2.0d0 );
			end
		end
		if ( ind > 1 ) 
            g(i) = g(i) * sin( x(m-ind+1) * pi / 2.0d0 );
        end
	end
	
	if ( ind > 1 ) 
		g(m-ind+1) = faux * pi / 2.0d0 * cos( x(m-ind+1) * pi / 2.0d0 );
		for j = 1:m-ind
			g(m-ind+1) = g(m-ind+1) * cos( x(j) * pi / 2.0d0 );
		end
	end
			
	
	faux = 1.0d0;
	for i = 1:m-ind
		faux = faux * cos( x(i) * pi / 2.0d0 );
	end

	if ( ind > 1 ) 
        faux = faux * sin( x(m-ind+1) * pi / 2.0d0 );
    end
	
	for i = m:n
		g(i) = faux * 1.0d2 * ( 2.0d0 * ( x(i) - 0.5d0 ) + 2.0d1 * pi * sin( 2.0d1 * pi * ( x(i) - 0.5d0 ) ) );
	end
				
	
	return

end

% ----------------------------------------------------------------------

%  DTLZ4

if ( strcmp(problem,'DTLZ4' )  )

%		k = dimk
%		m = dimm

%		faux = 0.0d0
%		for i = m:n
%			faux = faux + ( x(i) - 0.5d0 ) ^ 2
%		end
	
%		f = 1.0d0 + faux					
%		for i = 1:m-ind
%			f = f * cos( x(i)^alpha * pi / 2.0d0 )
%		end
	
%		if ( ind > 1 ) f = f * sin( x(m-ind+1)^alpha * pi / 2.0d0 )
	
	%k = dimk;
	%m = dimm;
    
    k = 5;
    m = 3;
    
    alpha = 2.0d0;
	
	for i = 1:n
		g(i) = 0.0d0;
	end

	faux = 0.0d0;
	for i = m:n
		faux = faux + ( x(i) - 0.5d0 ) ^ 2;
	end
	faux = 1.0d0 + faux;
	
	for i = 1:m-ind
		g(i) = faux;
	end 
	
	for i = 1:m-ind
		for j = 1:m-ind
			if ( j == i ) 
				g(i) = - g(i) * pi / 2.0d0 * alpha * x(i)^( alpha - 1.0d0 ) * sin( x(i)^alpha * pi / 2.0d0 );
			else
				g(i) = g(i) * cos( x(j)^alpha * pi / 2.0d0 );
			end
		end
		if ( ind > 1 ) 
            g(i) = g(i) * sin( x(m-ind+1)^alpha * pi / 2.0d0 );
        end
	end
	
	if ( ind > 1 ) 
		g(m-ind+1) = faux * pi / 2.0d0 * alpha * x(i)^( alpha - 1.0d0 ) * cos( x(m-ind+1)^alpha * pi / 2.0d0 );
		for j = 1:m-ind
			g(m-ind+1) = g(m-ind+1) * cos( x(j)^alpha * pi / 2.0d0 );
		end
	end
			
	
	faux = 1.0d0;
	for i = 1:m-ind
		faux = faux * cos( x(i)^alpha * pi / 2.0d0 );
	end

	if ( ind > 1 ) 
        faux = faux * sin( x(m-ind+1)^alpha * pi / 2.0d0 );
    end
	
	for i = m:n
		g(i) = faux * 2.0d0 * ( x(i) - 0.5d0 );
	end
				
	
	return

end
	
% ----------------------------------------------------------------------

%  FA1

if ( strcmp(problem,'FA1' )  )
		
	if ( ind == 1 ) 
		g(1) = 4.0d0 * exp( -4.0d0 * x(1) )  / ( 1.0d0 - exp( -4.0d0 ) );
		g(2) = 0.0d0;
		g(3) = 0.0d0;
		
		return
	end
	
	if ( ind == 2 ) 
		%faux = ( 1.0d0 - exp( -4.0d0 * x(1) ) ) / ( 1.0d0 - exp( -4.0d0 ) )
		%f = ( x(2) + 1.0d0 ) * ( 1.0d0 - ( faux / ( x(2) + 1.0d0  ) ) ^ 0.5d0 )
					
		a = exp( -4.0d0 * x(1) );
		b = 1.0d0 - exp( -4.0d0 );
		t = ( 1.0d0 - a ) / ( b * ( x(2) + 1.0d0 ) );
		
		g(1) =  - 2.0d0 * a / b *  t ^ (-0.5d0);
		g(2) = 1.0d0 - 0.5d0 * t ^ 0.5d0 ;
		g(3) = 0.0d0;
		
		return	
	end
	
	if ( ind == 3 ) 
		%faux = ( 1.0d0 - exp( -4.0d0 * x(1) ) ) / ( 1.0d0 - exp( -4.0d0 ) )
		%f = ( x(3) + 1.0d0 ) * ( 1.0d0 - ( faux / ( x(3) + 1.0d0  ) ) ^ 0.1d0 )

		a = exp( -4.0d0 * x(1) );
		b = 1.0d0 - exp( -4.0d0 );
		t = ( 1.0d0 - a ) / ( b * ( x(3) + 1.0d0  ) );
		
		g(1) = - 0.4d0 * t ^ (-0.9d0) * a / b ;
		g(2) = 0.0d0;
		g(3) = 1.0d0 - 0.9d0 * t ^ 0.1d0;
		
		return	
	end
end		

% ----------------------------------------------------------------------

%  Far1

if ( strcmp(problem,'Far1' )  )
		
	if ( ind == 1 ) 
		g(1) = 6.0d1 * ( x(1) - 0.1d0 ) * exp( 1.5d1 * ( - ( x(1) - 0.1d0 ) ^ 2 - x(2) ^ 2 ) ) ...
			 + 4.0d1 * ( x(1) - 0.6d0 ) * exp( 2.0d1 * ( - ( x(1) - 0.6d0 ) ^ 2 - ( x(2) - 0.6d0 ) ^ 2 ) ) ...
			 - 4.0d1 * ( x(1) + 0.6d0 ) * exp( 2.0d1 * ( - ( x(1) + 0.6d0 ) ^ 2 - ( x(2) - 0.6d0 ) ^ 2 ) ) ...
			 - 4.0d1 * ( x(1) - 0.6d0 ) * exp( 2.0d1 * ( - ( x(1) - 0.6d0 ) ^ 2 - ( x(2) + 0.6d0 ) ^ 2 ) ) ...
			 - 4.0d1 * ( x(1) + 0.6d0 ) * exp( 2.0d1 * ( - ( x(1) + 0.6d0 ) ^ 2 - ( x(2) + 0.6d0 ) ^ 2 ) ) ;
	
		g(2) = 6.0d1 * x(2) * exp( 1.5d1 * ( - ( x(1) - 0.1d0 ) ^ 2 - x(2) ^ 2 ) ) ...
			 + 4.0d1 * ( x(2) - 0.6d0 ) * exp( 2.0d1 * ( - ( x(1) - 0.6d0 ) ^ 2 - ( x(2) - 0.6d0 ) ^ 2 ) ) ...
			 - 4.0d1 * ( x(2) - 0.6d0 ) * exp( 2.0d1 * ( - ( x(1) + 0.6d0 ) ^ 2 - ( x(2) - 0.6d0 ) ^ 2 ) ) ...
			 - 4.0d1 * ( x(2) + 0.6d0 ) * exp( 2.0d1 * ( - ( x(1) - 0.6d0 ) ^ 2 - ( x(2) + 0.6d0 ) ^ 2 ) ) ...
			 - 4.0d1 * ( x(2) + 0.6d0 ) * exp( 2.0d1 * ( - ( x(1) + 0.6d0 ) ^ 2 - ( x(2) + 0.6d0 ) ^ 2 ) ) 	;	 
				 
		
		return
	end
	
	if ( ind == 2 ) 
		g(1) = - 8.0d1 * x(1) * exp( 2.0d1 * ( - x(1) ^ 2 - x(2) ^ 2 ) ) ...
				 - 4.0d1 * ( x(1) - 0.4d0 ) * exp( 2.0d1 * ( - ( x(1) - 0.4d0 ) ^ 2 - ( x(2) - 0.6d0 ) ^ 2 ) ) ...
				 + 4.0d1 * ( x(1) + 0.5d0 ) * exp( 2.0d1 * ( - ( x(1) + 0.5d0 ) ^ 2 - ( x(2) - 0.7d0 ) ^ 2 ) ) ...
				 + 4.0d1 * ( x(1) - 0.5d0 ) * exp( 2.0d1 * ( - ( x(1) - 0.5d0 ) ^ 2 - ( x(2) + 0.7d0 ) ^ 2 ) ) ...
				 - 4.0d1 * ( x(1) + 0.4d0 ) * exp( 2.0d1 * ( - ( x(1) + 0.4d0 ) ^ 2 - ( x(2) + 0.8d0 ) ^ 2 ) ) ;
					 
		g(2) = - 8.0d1 * x(2) * exp( 2.0d1 * ( - x(1) ^ 2 - x(2) ^ 2 ) ) ...
			   - 4.0d1 * ( x(2) - 0.6d0 ) * exp( 2.0d1 * ( - ( x(1) - 0.4d0 ) ^ 2 - ( x(2) - 0.6d0 ) ^ 2 ) ) ...
			   + 4.0d1 * ( x(2) - 0.7d0 ) * exp( 2.0d1 * ( - ( x(1) + 0.5d0 ) ^ 2 - ( x(2) - 0.7d0 ) ^ 2 ) ) ...
			   + 4.0d1 * ( x(2) + 0.7d0 ) * exp( 2.0d1 * ( - ( x(1) - 0.5d0 ) ^ 2 - ( x(2) + 0.7d0 ) ^ 2 ) ) ...
			   - 4.0d1 * ( x(2) + 0.8d0 ) * exp( 2.0d1 * ( - ( x(1) + 0.4d0 ) ^ 2 - ( x(2) + 0.8d0 ) ^ 2 ) ) ;
					 
		
		return	
	end
		
end		

% ----------------------------------------------------------------------

% FDS
% NEWTON???S METHOD FOR MULTIOBJECTIVE OPTIMIZATION

if ( strcmp(problem,'FDS' )  )
		
	if ( ind == 1 ) 
		for i = 1:n
			g(i) = 4.0d0 * i * ( x(i) - i ) ^ 3 / n ^ 2;
		end
		
		return
	end
	
	if ( ind == 2 ) 
		for i = 1:n
			g(i) = exp( sum(x)/n )/n + 2.0d0 * x(i);
		end
		
		return
	end	
	
	if ( ind == 3 ) 
		for i = 1:n
			g(i) = - i * ( n - i + 1.0d0 ) * exp( - x(i) ) / ( n * ( n + 1.0d0 ) ) ;
		end
		
		return
	end	
		
end		

% ----------------------------------------------------------------------

% FF1 
% C. M. Fonseca and P. J. Fleming: ???An overview of evolutionary algorithms in multiobjective optimization

if ( strcmp(problem,'FF1' )  )
	
	if ( ind == 1 ) 
		g(1) = 2.0d0 * ( x(1) - 1.0d0 ) * exp( - ( x(1) - 1.0d0 ) ^ 2 - ( x(2) + 1.0d0 ) ^ 2 );
		g(2) = 2.0d0 * ( x(2) + 1.0d0 ) * exp( - ( x(1) - 1.0d0 ) ^ 2 - ( x(2) + 1.0d0 ) ^ 2 );
		
		return
	end
	
	if ( ind == 2 ) 
		g(1) = 2.0d0 * ( x(1) + 1.0d0 ) * exp( - ( x(1) + 1.0d0 ) ^ 2 - ( x(2) - 1.0d0 ) ^ 2 );
		g(2) = 2.0d0 * ( x(2) - 1.0d0 ) * exp( - ( x(1) + 1.0d0 ) ^ 2 - ( x(2) - 1.0d0 ) ^ 2 );
		
		return
	end	
		
end						

% ----------------------------------------------------------------------


%  Hil1

if ( strcmp(problem,'Hil1' )  )
	a = 2.0d0 * pi / 3.6d2 * ( 4.5d1 + 4.0d1 * sin( 2.0d0 * pi * x(1) ) ...
		+ 2.5d1 * sin( 2.0d0 * pi * x(2) ) );
	b = 1.0d0 + 0.5d0 * cos( 2.0d0 * pi * x(1) );
	
	if ( ind == 1 ) 
		g(1) = - 1.6d2 * pi ^ 2 / 3.6d2 * cos( 2.0d0 * pi * x(1) ) * sin( a ) *  b ...
		- pi * sin( 2.0d0 * pi * x(1) ) * cos( a ) ;
		g(2) = - 1.0d2 * pi ^ 2 / 3.6d2 * cos( 2.0d0 * pi * x(2) ) * sin( a ) *  b ;
		
		return
	end
	
	if ( ind == 2 ) 
		g(1) = 1.6d2 * pi ^ 2 / 3.6d2 * cos( 2.0d0 * pi * x(1) ) * cos( a ) *  b ...
		- pi * sin( 2.0d0 * pi * x(1) ) * sin( a ) ;
		g(2) = 1.0d2 * pi ^ 2 / 3.6d2 * cos( 2.0d0 * pi * x(2) ) * cos( a ) *  b ;
		
		return	
	end
		
end

% ----------------------------------------------------------------------

%  IKK1

if ( strcmp(problem,'IKK1' )  )
	
	if ( ind == 1 ) 
		g(1) = 2.0d0 * x(1);
		g(2) = 0.0d0;
		
		return
	end
	
	if ( ind == 2 ) 
		g(1) = 2.0d0 * ( x(1) - 2.0d1 );
		g(2) = 0.0d0;
		
		return	
	end
	
	if ( ind == 3 ) 
		g(1) = 0.0d0;
		g(2) = 2.0d0 * x(2);
		
		return	
	end
		
end				

% ----------------------------------------------------------------------

%  IM1

if ( strcmp(problem,'IM1' )  )
	
	if ( ind == 1 ) 
		g(1) = 1.0d0 / sqrt( x(1) );
		g(2) = 0.0d0;
		
		return
	end
	
	if ( ind == 2 ) 			
		g(1) = ( 1.0d0 - x(2) );
		g(2) = - x(1);
		
		return	
	end
		
end

% ----------------------------------------------------------------------	

% JOS1 
% Dynamic Weighted Aggregation for Evolutionary Multi-Objetive Optimization: Why fores It Work and How?

if ( strcmp(problem,'JOS1' ) )

	if ( ind == 1 ) 					
		for i = 1:n
			g(i) = 2.0d0 * x(i) / n;
		end
		
		return
	end
	
	if ( ind == 2 ) 
		for i = 1:n
			g(i) = 2.0d0 * ( x(i) - 2.0d0 ) / n;
		end	
		
		return				
	end
		
end

% ----------------------------------------------------------------------	

% JOS4
% Dynamic Weighted Aggregation for Evolutionary Multi-Objetive Optimization: Why fores It Work and How?

if ( strcmp(problem,'JOS4' )  )

	if ( ind == 1 ) 
		g(1) = 1.0d0;
		g(2:n) = 0.0d0;
		
		return
	end
	
	if ( ind == 2 ) 
		faux = 1.0d0 + 9.0d0 * sum(x(2:n)) / ( n - 1 );
		t = x(1) / faux;
		
		g(1) = - 0.25d0 * t ^ (-0.75d0) - 4.0d0 * t ^ 3;
		for i = 2:n
			g(i) = 9.0d0 / ( n - 1 ) * ( 1.0d0 - 0.75d0 * t ^ 0.25d0 + 3.0d0 * t ^ 4.0d0 );
		end
		
		
		return
	end

end		

% ----------------------------------------------------------------------

% 	KW2

if ( strcmp(problem,'KW2' )  )
		
	if ( ind == 1 ) 
		g(1) = 6.0d0 * ( 1.0d0 - x(1) ) * exp( -x(1)^2 - ( x(2) + 1.0d0 ) ^ 2 ) ...
			  + 6.0d0 * ( 1.0d0 - x(1) )^2 * exp( -x(1)^2 - ( x(2) + 1.0d0 ) ^ 2 ) * x(1) ...
			  + 1.0d1 * ( 1.0d0 / 5.0d0 - 3.0d0 * x(1)^2 ) * exp( - x(1)^2 - x(2)^2 ) ...
			  - 2.0d1 * ( x(1) / 5.0d0 - x(1)^3 - x(2)^5 ) * exp( - x(1)^2 - x(2)^2 ) * x(1) ...
			  - 6.0d0 * exp( -( x(1) + 2.0d0 )^2 - x(2)^2 ) * ( x(1) + 2.0d0 ) - 1.0d0;
		g(2) = 6.0d0 * ( 1.0d0 - x(1) )^2 * exp( -x(1)^2 - ( x(2) + 1.0d0 ) ^ 2 ) * ( x(2) + 1.0d0 ) ...
					- 5.0d1 * x(2)^4 * exp( - x(1)^2 - x(2)^2 ) ...
					- 1.0d1 * ( x(1) / 5.0d0 - x(1)^3 - x(2)^5 ) * exp( - x(1)^2 - x(2)^2 ) * 2.0d0 * x(2)...
					- 6.0d0 * exp( -( x(1) + 2.0d0 )^2 - x(2)^2 ) * x(2) - 0.5d0;
		
		return
	end
	
	if ( ind == 2 ) 
		g(1) = - 6.0d0 * ( 1.0d0 + x(2) )^2 * exp( -x(2)^2 - ( 1.0d0 - x(1) ) ^ 2 ) * ( 1.0d0 - x(1) ) ...
					+ 5.0d1 * x(1)^4 * exp( - x(1)^2 - x(2)^2 ) ...
					- 2.0d1 * ( - x(2) / 5.0d0 + x(2)^3 + x(1)^5 ) * exp( - x(1)^2 - x(2)^2 ) * x(1) ...
					- 6.0d0 * exp( -( 2.0d0 - x(2) )^2 - x(1)^2 ) * x(1);
		
		g(2) = - 6.0d0 * ( 1.0d0 + x(2) ) * exp( -x(2)^2 - ( 1.0d0 - x(1) ) ^ 2 ) ...
					 + 6.0d0 * ( 1.0d0 + x(2) )^2 * exp( -x(2)^2 - ( 1.0d0 - x(1) ) ^ 2 ) * x(2) ...
					 + 1.0d1 * ( - 1.0d0 / 5.0d0 + 3.0d0 * x(2)^2 ) * exp( - x(1)^2 - x(2)^2 ) ...
					 - 2.0d1 * ( - x(2) / 5.0d0 + x(2)^3 + x(1)^5 ) * exp( - x(1)^2 - x(2)^2 ) * x(2)...
					 + 6.0d0 * exp( -( 2.0d0 - x(2) )^2 - x(1)^2 ) * ( 2.0d0 - x(2) );
		
		return	
	end
		
end


% ----------------------------------------------------------------------

% 	LE1

if ( strcmp(problem,'LE1' )  )
										
	if ( ind == 1 ) 			
		t = 0.25d0 * ( x(1) ^ 2 + x(2) ^ 2 ) ^ (-0.875d0);
		
		g(1) = x(1) * t;
		g(2) = x(2) * t;
		
		return
	end
	
	if ( ind == 2 ) 
		t = 0.5d0 * ( ( x(1) - 0.5d0 ) ^ 2 + ( x(2) - 0.5d0 ) ^ 2 ) ^ (-0.75d0);
		
		g(1) = ( x(1) - 0.5d0 ) * t;
		g(2) = ( x(2) - 0.5d0 ) * t;
		
		return
	end
		
end

% ----------------------------------------------------------------------

% Lov1  

if ( strcmp(problem,'Lov1' )  )
		
	if ( ind == 1 ) 
		g(1) = 2.1d0 * x(1);
		g(2) = 2.0d0 * 0.98d0 * x(2);
		
		return
	end
	
	if ( ind == 2 ) 
		g(1) = 2.0d0 * 0.99d0 * ( x(1) - 3.0d0 );
		g(2) = 2.0d0 * 1.03d0 * ( x(2) - 2.5d0 );
		
		return
	end	
		
end	

% ----------------------------------------------------------------------

% Lov2

if ( strcmp(problem,'Lov2' )  )

	if ( ind == 1 ) 
		g(1) = 0.0d0;
		g(2) = 1.0d0;
		
		return
	end
	
	if ( ind == 2 ) 
		g(1) =  ( - 3.0d0 * x(1) ^ 2 * ( x(1) + 1.0d0 ) - ( x(2) - x(1) ^ 3 ) );
		g(1) = - g(1) / ( x(1) + 1.0d0 ) ^ 2;
		
		g(2) = - 1.0d0 / ( x(1) + 1.0d0 ) ;
		
		return
	end	
		
end		

% ----------------------------------------------------------------------

% Lov3 

if ( strcmp(problem,'Lov3' )  )
	
	if ( ind == 1 ) 
		g(1) =  2.0d0 * x(1);
		g(2) =  2.0d0 * x(2);
		
		return
	end
	
	if ( ind == 2 ) 
		g(1) =  2.0d0 * ( x(1) - 6.0d0 );
		g(2) = - 2.0d0 * ( x(2) + 0.3d0 );
		
		return
	end	
		
end	

% ----------------------------------------------------------------------

% Lov4 

if ( strcmp(problem,'Lov4' )  )
		
	if ( ind == 1 ) 
		g(1) =  2.0d0 * x(1) - 8.0d0 * ( ( x(1) + 2.0d0 ) * exp( - ( x(1) + 2.0d0 ) ^ 2 - x(2) ^ 2 ) ...
		+ ( x(1) - 2.0d0 ) * exp( - ( x(1) - 2.0d0 ) ^ 2 - x(2) ^ 2 ) );
		g(2) =  2.0d0 * x(2) - 8.0d0 * ( x(2) * exp( - ( x(1) + 2.0d0 ) ^ 2 - x(2) ^ 2 ) ...
		+ x(2) * exp( - ( x(1) - 2.0d0 ) ^ 2 - x(2) ^ 2 ) );
		
		return
	end
	
	if ( ind == 2 ) 
		g(1) =  2.0d0 * ( x(1) - 6.0d0 );
		g(2) =  2.0d0 * ( x(2) + 0.5d0 );
		
		return
	end	
		
end	

% ----------------------------------------------------------------------

% Lov5

if ( strcmp(problem,'Lov5' )  )

	MM = [ -1.0d0  , -0.03d0,  0.011d0; ...
	       -0.03d0 , -1.0d0 ,  0.07d0 ; ...
	        0.011d0,  0.07d0, -1.01d0 ];
	
	p = [ x(1); x(2) - 0.15d0;  x(3)];
	a = 0.35d0;
	
	A1 = sqrt( 2.0d0 * pi / a ) * exp( dot( p, MM *p ) / a ^ 2 );
	
	p = [ x(1) ; x(2) + 1.1d0;  0.5d0 * x(3) ];
	a = 3.0d0;
	
	A2 = sqrt( 2.0d0 * pi / a ) * exp( dot( p, MM*p ) / a ^ 2 );

	if ( ind == 1 ) 
		g(1) = sqrt(2.0d0)/2.0d0 + sqrt(2.0d0)/2.0d0 * A1 * ...
		( 2.0d0 * MM(1,1)* x(1) + 2.0d0 * MM(1,3) * x(3) + 2.0d0 * MM(1,2) * ( x(2) - 0.15d0 ) ) / 0.35d0^2 ...
		+ sqrt(2.0d0)/2.0d0 * A2 * ( 2.0d0 * MM(1,1) * x(1) + MM(1,3) * x(3) + 2.0d0 * MM(1,2) * ( x(2) + 1.1d0 ) ) / 3.0d0^2;
		g(1) = - g(1);
		
		g(2) = sqrt(2.0d0)/2.0d0 * A1 * ...
		( 2.0d0 * MM(1,2)* x(1) + 2.0d0 * MM(2,3) * x(3) + 2.0d0 * MM(2,2) * ( x(2) - 0.15d0 ) ) / 0.35d0^2 ...
		+ sqrt(2.0d0)/2.0d0 * A2 * ( 2.0d0 * MM(1,2) * x(1) + MM(2,3) * x(3) + 2.0d0 * MM(2,2) * ( x(2) + 1.1d0 ) ) / 3.0d0^2;
		g(2) = - g(2);
		
		g(3) = sqrt(2.0d0)/2.0d0 * A1 * ...
		( 2.0d0 * MM(1,3)* x(1) + 2.0d0 * MM(3,3) * x(3) + 2.0d0 * MM(2,3) * ( x(2) - 0.15d0 ) ) / 0.35d0^2 ...
		+ sqrt(2.0d0)/2.0d0 * A2 * ( MM(1,3) * x(1) + MM(3,3) * x(3) / 2.0d0 + MM(2,3) * ( x(2) + 1.1d0 ) ) / 3.0d0^2;
		g(3) = - g(3);
		
		return
	end
	
	if ( ind == 2 ) 
		g(1) = - sqrt(2.0d0)/2.0d0 + sqrt(2.0d0)/2.0d0 * A1 * ...
		( 2.0d0 * MM(1,1)* x(1) + 2.0d0 * MM(1,3) * x(3) + 2.0d0 * MM(1,2) * ( x(2) - 0.15d0 ) ) / 0.35d0^2 ...
		+ sqrt(2.0d0)/2.0d0 * A2 * ( 2.0d0 * MM(1,1) * x(1) + MM(1,3) * x(3) + 2.0d0 * MM(1,2) * ( x(2) + 1.1d0 ) ) / 3.0d0^2;
		g(1) = - g(1);
		
		g(2) = sqrt(2.0d0)/2.0d0 * A1 * ...
		( 2.0d0 * MM(1,2)* x(1) + 2.0d0 * MM(2,3) * x(3) + 2.0d0 * MM(2,2) * ( x(2) - 0.15d0 ) ) / 0.35d0^2 ...
		+ sqrt(2.0d0)/2.0d0 * A2 * ( 2.0d0 * MM(1,2) * x(1) + MM(2,3) * x(3) + 2.0d0 * MM(2,2) * ( x(2) + 1.1d0 ) ) / 3.0d0^2;
		g(2) = - g(2);
		
		g(3) = sqrt(2.0d0)/2.0d0 * A1 * ...
		( 2.0d0 * MM(1,3)* x(1) + 2.0d0 * MM(3,3) * x(3) + 2.0d0 * MM(2,3) * ( x(2) - 0.15d0 ) ) / 0.35d0^2 ...
		+ sqrt(2.0d0)/2.0d0 * A2 * ( MM(1,3) * x(1) + MM(3,3) * x(3) / 2.0d0 + MM(2,3) * ( x(2) + 1.1d0 ) ) / 3.0d0^2;
		g(3) = - g(3);
		
		return
	end	
		
end		

% ----------------------------------------------------------------------

% Lov6

if ( strcmp(problem,'Lov6' )  )

	if ( ind == 1 ) 
		g(1)   = 1.0d0;
		g(2:6) = 0.0d0;
		
		return
	end
	
	if ( ind == 2 ) 
		g(1) = - 0.5d0 / sqrt( x(1) ) - sin( 1.0d1 * pi * x(1) ) - 1.0d1 * pi * x(1) * cos( 1.0d1 * pi * x(1) );
		for i = 2:6
			g(i) = 2.0d0 * x(i);
		end
		
		return
	end	
		
end		

% ----------------------------------------------------------------------	

% 	LTDZ
%	Combining convergence and diversity in evolutionary multiobjective optimization

if ( strcmp(problem,'LTDZ' )  )
										
	if ( ind == 1 ) 
		g(1) = pi / 2.0d0 * ( 1.0d0 + x(3) ) * sin( x(1) * pi / 2.0d0 ) * cos( x(2) * pi / 2.0d0 );
		g(2) = pi / 2.0d0 * ( 1.0d0 + x(3) ) * cos( x(1) * pi / 2.0d0 ) * sin( x(2) * pi / 2.0d0 );
		g(3) = - cos( x(1) * pi / 2.0d0 ) * cos( x(2) * pi / 2.0d0 );
		g = - g;
		
		return
	end
	
	if ( ind == 2 ) 
		g(1) = pi / 2.0d0 * ( 1.0d0 + x(3) ) * sin( x(1) * pi / 2.0d0 ) * sin( x(2) * pi / 2.0d0 );
		g(2) = - pi / 2.0d0 * ( 1.0d0 + x(3) ) * cos( x(1) * pi / 2.0d0 ) * cos( x(2) * pi / 2.0d0 );
		g(3) = - cos( x(1) * pi / 2.0d0 ) * sin( x(2) * pi / 2.0d0 );
		g = - g;
		
		return
	end
	
	if ( ind == 3 ) 
		g(1) = pi / 2.0d0 * ( 1.0d0 + x(3) ) * sin( x(1) * pi / 2.0d0 ) * sin( x(1) * pi / 2.0d0 ) ...
			 - pi / 2.0d0 * ( 1.0d0 + x(3) ) * cos( x(1) * pi / 2.0d0 ) * cos( x(1) * pi / 2.0d0 );
		g(2) = 0.0d0;
		g(3) = - cos( x(1) * pi / 2.0d0 ) * sin( x(1) * pi / 2.0d0 );
		g = - g;
		
		return
	end
		
end		

% ----------------------------------------------------------------------

% 	MGH9 

if ( strcmp(problem,'MGH9' )  )
			
	t = ( 8.0d0 - ind ) / 2.0d0;
	
	g(1) = exp( - x(2) * ( t - x(3) ) ^ 2 / 2.0d0 );
	g(2) = - x(1) * exp( - x(2) * ( t - x(3) ) ^ 2 / 2.0d0 ) * ( t - x(3) ) ^ 2 / 2.0d0;
	g(3) = x(1) * exp( - x(2) * ( t - x(3) ) ^ 2 / 2.0d0 ) * x(2) * ( t - x(3) ) ;
					
	
		return
		
end		

% ----------------------------------------------------------------------

% 	MGH16 

if ( strcmp(problem,'MGH16' )  )
		
	t = ind / 5.0d0;
	
	g(1) = 2.0d0 * ( x(1) + t * x(2) - exp(t) );
	g(2) = 2.0d0 * t * ( x(1) + t * x(2) - exp(t) );
	g(3) = 2.0d0 * ( x(3) + x(4) * sin(t) - cos(t) );
	g(4) = 2.0d0 * sin(t) * ( x(3) + x(4) * sin(t) - cos(t) );
						
	
		return
		
end		

% ----------------------------------------------------------------------

% 	MGH26 

if ( strcmp(problem,'MGH26' )  )
							
	t = 0.0d0;
	for i = 1:n
		t = t + cos(x(i));
	end
	
	gaux1 = 2.0d0 * ( n - t + ind * ( 1.0d0 - cos(x(ind)) ) - sin(x(ind)) );
	
	
	for i = 1:n
		g(i) = gaux1 * sin(x(i));
	end
	
%		g(1) = gaux1 * sin(x(1))
%		g(2) = gaux1 * sin(x(2))
%		g(3) = gaux1 * sin(x(3))
%		g(4) = gaux1 * sin(x(4))
	
	g(ind) = g(ind)  + gaux1 * ( ind * sin(x(ind)) - cos(x(ind)) );
			
	
		return
	
end				

% ----------------------------------------------------------------------

% 	MGH33

if ( strcmp(problem,'MGH33' )  )

	faux = 0.0d0;			
	for i = 1:n
		faux = faux + i * x(i);
	end
	
	faux = 2.0d0 * ( ind * faux - 1.0d0 );

	for i = 1:n
		g(i) = faux * real(i) * real(ind);
	end
					
	
	return
		
end				

% ----------------------------------------------------------------------

% 	MHHM2

if ( strcmp(problem,'MHHM2' )  )
										
	if ( ind == 1 ) 
		g(1) = 2.0d0 * ( x(1) - 0.8d0 );
		g(2) = 2.0d0 * ( x(2) - 0.6d0 );
		
		return
	end
	
	if ( ind == 2 ) 			
		g(1) = 2.0d0 * ( x(1) - 0.85d0 );
		g(2) = 2.0d0 * ( x(2) - 0.7d0 );
		
		return
	end
		
	if ( ind == 3 ) 			
		g(1) = 2.0d0 * ( x(1) - 0.9d0 );
		g(2) = 2.0d0 * ( x(2) - 0.6d0 );
		
		return
	end
			
end	

% ----------------------------------------------------------------------

%  MLF1

if ( strcmp(problem,'MLF1' )  )
		
	if ( ind == 1 ) 
		g(1) = sin( x(1) ) / 2.0d1 + ( 1.0d0 + x(1) / 2.0d1 ) * cos( x(1) );
		
		return
	end
	
	if ( ind == 2 ) 
		g(1) = cos( x(1) ) / 2.0d1 - ( 1.0d0 + x(1) / 2.0d1 ) * sin( x(1) );
		
		return	
	end
		
end		

% ----------------------------------------------------------------------

%  MLF2

if ( strcmp(problem,'MLF2' )  )
		
	if ( ind == 1 ) 
		g(1) = ( 2.0d0 * x(1) * ( x(1) ^ 2 + x(2) - 1.1d1 ) + ( x(1) + x(2) ^ 2 - 7.0d0 )  ) / 1.0d2;
		g(2) = ( ( x(1) ^ 2 + x(2) - 1.1d1 ) + 2.0d0 * x(2) * ( x(1) + x(2) ^ 2 - 7.0d0 )  ) / 1.0d2;
		
		return
	end
	
	if ( ind == 2 ) 
		g(1) = ( 8.0d0 * x(1) * ( 4.0d0 * x(1) ^ 2 + 2.0d0 * x(2) - 1.1d1 ) ...
			 + 2.0d0 * ( 2.0d0 * x(1) + 4.0d0 *  x(2) ^ 2 - 7.0d0 ) ) / 1.0d2;
		g(2) = ( 2.0d0 * ( 4.0d0 * x(1) ^ 2 + 2.0d0 * x(2) - 1.1d1 ) ...
			 + 8.0d0 * x(2) * ( 2.0d0 * x(1) + 4.0d0 *  x(2) ^ 2 - 7.0d0 ) ) / 1.0d2;
		
		return	
	end
		
end	

% ----------------------------------------------------------------------

%  MMR1			
%  Box-constrained multi-objective optimization: A gradient-like method without ??????a priori?????? scalarization	

%	if ( strcmp(problem,'MMR1' )  )
		
%			if ( ind == 1 ) 
%				g(1) = 2.0d0 * x(1) 
%				g(2) = 0.0d0	
%				
%			return
%			end
		
%			if ( ind == 2 ) 
%				g(1) = 2.0d0 - 0.8d0 * exp( - ( ( x(2) - 0.6d0 ) / 0.4d0 ) ^ 2 ) ...
%					   - exp( - ( ( x(2) - 0.2d0 ) / 0.04d0 ) ^ 2 )
%				g(1) = - 2.0d0 * x(1) * g(1) / ( 1.0d0 + x(1) ^ 2 ) ^ 2
			
%				g(2) = 1.0d1 * exp( - ( ( x(2) - 0.6d0 ) / 0.4d0 ) ^ 2 ) * ( x(2) - 0.6d0 ) ...
%				   + 1.25d3 * exp( - ( ( x(2) - 0.2d0 ) / 0.04d0 ) ^ 2 ) * ( x(2) - 0.2d0 ) 
%				g(2) = g(2) / ( 1.0d0 + x(1) ^ 2 )   
%				
%			return	
%			end	
		
%	end	

if ( strcmp(problem,'MMR1' )  )
		
		if ( ind == 1 ) 
			g(1) = 1.0d0;
			g(2) = 0.0d0;	
			
		return
		end
		
		if ( ind == 2 ) 
			g(1) = 2.0d0 - 0.8d0 * exp( - ( ( x(2) - 0.6d0 ) / 0.4d0 ) ^ 2 ) ...
				   - exp( - ( ( x(2) - 0.2d0 ) / 0.04d0 ) ^ 2 );
			g(1) = - g(1) / x(1) ^ 2;
			
			g(2) = 1.0d1 * exp( - ( ( x(2) - 0.6d0 ) / 0.4d0 ) ^ 2 ) * ( x(2) - 0.6d0 ) ...
			   + 1.25d3 * exp( - ( ( x(2) - 0.2d0 ) / 0.04d0 ) ^ 2 ) * ( x(2) - 0.2d0 ) ;
			g(2) = g(2) / x(1);
			
		return	
		end	
		
end		

% ----------------------------------------------------------------------

%  MMR2
%  Box-constrained multi-objective optimization: A gradient-like method without ??????a priori?????? scalarization	

if ( strcmp(problem,'MMR2' )  )
		
	if ( ind == 1 ) 
		g(1) = 1.0d0;
		g(2) = 0.0d0;	
		
		return
	end
	
	if ( ind == 2 ) 
	
%			faux = x(1) / ( 1.0d0 + 1.0d1 * x(2) )
%			f = 1.0d0 - faux ^ 2 - faux * sin( 8.0d0 * pi * x(1) )
%			f = f * ( 1.0d0 + 1.0d1 * x(2) )

		a = 1.0d0 + 1.0d1 * x(2);
		faux = x(1) / a;
		
		g(1) = - 2.0d0 * faux / a - sin( 8.0d0 * pi * x(1) ) / a - 8.0d0 * pi * faux * cos( 8.0d0 * pi * x(1) );
		g(1) = g(1) * a;
		
		g(2) = 1.0d1 * ( 1.0d0 - faux ^ 2 - faux * sin( 8.0d0 * pi * x(1) ) );
		g(2) = g(2) + a * ( 2.0d1 * faux * x(1) / a^2 + 1.0d1 / a^2 * x(1) * sin( 8.0d0 * pi * x(1) ) );
	
		
		return	
	end	
		
end		
	
% ----------------------------------------------------------------------

%  MMR3
%  Box-constrained multi-objective optimization: A gradient-like method without ??????a priori?????? scalarization	

if ( strcmp(problem,'MMR3' )  )
		
	if ( ind == 1 ) 
		g(1) = 3.0d0 * x(1) ^ 2;
		g(2) = 0.0d0;
		
		return
	end
	
	if ( ind == 2 ) 
		g(1) = - 3.0d0 * ( x(2) - x(1) ) ^ 2;
		g(2) = 3.0d0 * ( x(2) - x(1) ) ^ 2;
		
		return	
	end	
		
end			

% ----------------------------------------------------------------------

%  MMR4
%  Box-constrained multi-objective optimization: A gradient-like method without ??????a priori?????? scalarization	

if ( strcmp(problem,'MMR4' )  )
		
	if ( ind == 1 ) 	
		t = ( 2.0d0 * x(1) + x(2) + 2.0d0 * x(3) + 1.0d0 ) ^ 2;
		g(1) = 1.0d0 + 7.2d1 / t;
		g(2) = - 2.0d0 + 3.6d1 / t;
		g(3) = - 1.0d0 + 7.2d1 / t;
		
		return
	end
	
	if ( ind == 2 ) 			
		g(1) = - 3.0d0;
		g(2) =   1.0d0;
		g(3) = - 1.0d0;
		
		return	
	end	
		
end					
	
% ----------------------------------------------------------------------

%  MOP 2	

if ( strcmp(problem,'MOP2' )  )
		
	if ( ind == 1 ) 
		faux = 0.0d0;
		for i = 1:n
			faux = faux + ( x(i) - 1.0d0 / ( sqrt(real(n)) ) ) ^ 2;
		end
		
		for i = 1:n
			g(i) = 2.0d0 * ( x(i) - 1.0d0 / sqrt(real(n)) ) * exp( - faux );
		end
		
		return
	end
	
	if ( ind == 2 ) 
		faux = 0.0d0;
		for i = 1:n
			faux = faux + ( x(i) + 1.0d0 / ( sqrt(real(n)) ) ) ^ 2;
		end
		
		for i = 1:n
			g(i) = 2.0d0 * ( x(i) + 1.0d0 / sqrt(real(n)) ) * exp( - faux );
		end
		
		return	
	end	
		
end			

% ----------------------------------------------------------------------

%  MOP 3	

if ( strcmp(problem,'MOP3' )  )
		
	if ( ind == 1 ) 
		A1 = 0.5d0 * sin(1.0d0) - 2.0d0 * cos(1.0d0) + sin(2.0d0) - 1.5d0 * cos(2.0d0);
		A2 = 1.5d0 * sin(1.0d0) - cos(1.0d0) + 2.0d0 * sin(2.0d0) - 0.5d0 * cos(2.0d0);
		B1 = 0.5d0 * sin(x(1)) - 2.0d0 * cos(x(1)) + sin(x(2)) - 1.5d0 * cos(x(2)); 
		B2 = 1.5d0 * sin(x(1)) - cos(x(1)) + 2.0d0 * sin(x(2)) - 0.5d0 * cos(x(2));	
		g(1) = 2.0d0 * ( A1 - B1 ) *  ( - 0.5d0 * cos(x(1)) - 2.0d0 * sin(x(1)) ) +...
			   2.0d0 * ( A2 - B2 ) *  ( - 1.5d0 * cos(x(1)) - sin(x(1)) );
		
		g(2) = 2.0d0 * ( A1 - B1 ) *  ( - cos(x(2)) - 1.5d0 * sin(x(2)) ) +...
			   2.0d0 * ( A2 - B2 ) *  ( - 2.0d0 * cos(x(2)) - 0.5d0 * sin(x(2)) );
		
		return
	end
	
	if ( ind == 2 ) 
		g(1) =  2.0d0 * ( x(1) + 3.0d0 );
		g(2) =  2.0d0 * ( x(2) + 1.0d0 );
		
		return	
	end	
		
end

% ----------------------------------------------------------------------

%  MOP 5	

if ( strcmp(problem,'MOP5' )  )
		
	if ( ind == 1 ) 
		g(1) = x(1) + 2.0d0 * x(1) * cos( x(1) ^ 2 + x(2) ^ 2 );
		g(2) = x(2) + 2.0d0 * x(2) * cos( x(1) ^ 2 + x(2) ^ 2 );
		
		return
	end
	
	if ( ind == 2 ) 
		g(1) = 3.0d0 * ( 3.0d0 * x(1) - 2.0d0 * x(2) + 4.0d0 ) / 4.0d0 ...
			 + 2.0d0 * ( x(1) - x(2) + 1.0d0 ) / 2.7d1;
		g(2) = - 2.0d0 * ( 3.0d0 * x(1) - 2.0d0 * x(2) + 4.0d0 ) / 4.0d0 ...
			 - 2.0d0 * ( x(1) - x(2) + 1.0d0 ) / 2.7d1; 	 
		
		return	
	end
	
	if ( ind == 3 ) 
		g(1) = - 2.0d0 * x(1) / ( x(1) ^ 2 + x(2) ^ 2 + 1.0d0 ) ^ 2 ...
			   + 2.2d0 * x(1) * exp( - x(1) ^ 2 - x(2) ^ 2 );
		g(2) = - 2.0d0 * x(2) / ( x(1) ^ 2 + x(2) ^ 2 + 1.0d0 ) ^ 2 ...
			   + 2.2d0 * x(2) * exp( - x(1) ^ 2 - x(2) ^ 2 );
		
		return
	end	
			
end		

% ----------------------------------------------------------------------

%  MOP6

if ( strcmp(problem,'MOP6' )  )
		
	if ( ind == 1 ) 
		g(1) = 1.0d0;
		g(2) = 0.0d0;
		
		return
	end
	
	if ( ind == 2 ) 
		a = 1.0d0 + 1.0d1 * x(2) ;
		b = sin( 8.0d0 * pi * x(1) );
		t = x(1) / a;
		
		g(1) = - 2.0d0 * t - b - 8.0d0 * pi * x(1) * cos( 8.0d0 * pi * x(1) );
		g(2) = 1.0d1 * ( 1.0d0 - t ^ 2 - t * b ) + 1.0d1 * x(1) / a * ( 2.0d0 * t + b );
		
		return	
	end
		
end		

% ----------------------------------------------------------------------

%  MOP 7

if ( strcmp(problem,'MOP7' )  )
		
	if ( ind == 1 ) 
		g(1) =  x(1) - 2.0d0;
		g(2) =  2.0d0 * ( x(2) + 1.0d0 ) / 1.3d1;
		
		return
	end
	
	if ( ind == 2 ) 
		g(1) =  ( x(1) + x(2) - 3.0d0 ) / 1.8d1 - ( - x(1) + x(2) + 2.0d0 ) / 4.0d0;
		g(2) =  ( x(1) + x(2) - 3.0d0 ) / 1.8d1 + ( - x(1) + x(2) + 2.0d0 ) / 4.0d0;
		
		return	
	end
	
	if ( ind == 3 ) 
		g(1) =  2.0d0 * ( x(1) + 2.0d0 * x(2) - 1.0d0 ) / 1.75d2 ...
		- 2.0d0 * ( - x(1) + 2.0d0 * x(2) ) / 1.7d1;
		g(2) =  4.0d0 * ( x(1) + 2.0d0 * x(2) - 1.0d0 ) / 1.75d2 ...
		+ 4.0d0 * ( - x(1) + 2.0d0 * x(2) ) / 1.7d1;
		
		return
	end			
		
end		

% ----------------------------------------------------------------------

% 	PNR

if ( strcmp(problem,'PNR' )  )
									
	if ( ind == 1 ) 
		g(1) = 4.0d0 * x(1) ^ 3 - 2.0d0 * x(1) - 1.0d1 * x(2);
		g(2) = 4.0d0 * x(2) ^ 3 + 2.0d0 * x(2) - 1.0d1 * x(1);
		
		return
	end
	
	if ( ind == 2 ) 
		g(1) = 2.0d0 * x(1);
		g(2) = 2.0d0 * x(2);
		
		return
	end
		
end		

% ----------------------------------------------------------------------

%  QV1

if ( strcmp(problem,'QV1' )  )
		
	if ( ind == 1 ) 
		faux = 0.0d0;
		for i = 1:n
			faux = faux + x(i) ^ 2 - 1.0d1 * cos( 2.0d0 * pi * x(i) ) + 1.0d1;
		end
		faux = 0.25 * ( faux / n ) ^ (-0.75d0);
		
		for i = 1:n
			g(i) = faux * ( 2.0d0 * x(i) + 2.0d1 * pi * sin( 2.0d0 * pi * x(i) )  ) / n;
		end
		
		
		return
	end
	
	if ( ind == 2 ) 
		faux = 0.0d0;
		for i = 1:n
			faux = faux + ( x(i) - 1.5d0 ) ^ 2 - 1.0d1 * cos( 2.0d0 * pi * ( x(i) -1.5d0 ) ) + 1.0d1;
		end
	
		faux = 0.25 * ( faux / n ) ^ (-0.75d0);
		
		for i = 1:n
			g(i) = faux * ( 2.0d0 * ( x(i) - 1.5d0 ) + 2.0d1 * pi * sin( 2.0d0 * pi * ( x(i) - 1.5d0 ) )  ) / n;
		end
		
		return	
	end
		
end

% ----------------------------------------------------------------------

% SD

if ( strcmp(problem,'SD' )  )
		
	if ( ind == 1 ) 
		g(1) = 2.0d0;
		g(2) = sqrt(2.0d0);
		g(3) = sqrt(2.0d0);
		g(4) = 1.0d0;
		
		return
	end
	
	if ( ind == 2 ) 
		g(1) = - 2.0d0 / x(1) ^ 2;
		g(2) = - 2.0d0 * sqrt(2.0d0) / x(2) ^ 2;
		g(3) = - 2.0d0 * sqrt(2.0d0) / x(3) ^ 2;
		g(4) = - 2.0d0 / x(4) ^ 2;
		
		return		
	end	
		
end		

% ----------------------------------------------------------------------

% SLCDT1

if ( strcmp(problem,'SLCDT1' )  )

	if ( ind == 1 ) 			
		g(1) = 0.5d0 * ( ( x(1) + x(2) )/sqrt( 1.0d0 + ( x(1) + x(2) ) ^ 2 ) + ...
					( x(1) - x(2) )/sqrt( 1.0d0 + ( x(1) - x(2) ) ^ 2 ) + 1.0d0 )  ...
					 - 2.0d0 * 0.85d0 *  ( x(1) + x(2) ) * exp( - ( x(1) + x(2) ) ^ 2 );
		g(2) = 0.5d0 * ( ( x(1) + x(2) )/sqrt( 1.0d0 + ( x(1) + x(2) ) ^ 2 ) - ...
					( x(1) - x(2) )/sqrt( 1.0d0 + ( x(1) - x(2) ) ^ 2 ) - 1.0d0 )  ...
					 - 2.0d0 * 0.85d0 * ( x(1) + x(2) ) * exp( - ( x(1) + x(2) ) ^ 2 );			 
		
		return			 
	end
	
	if ( ind == 2 ) 
		g(1) = 0.5d0 * ( ( x(1) + x(2) )/sqrt( 1.0d0 + ( x(1) + x(2) ) ^ 2 ) + ...
					( x(1) - x(2) )/sqrt( 1.0d0 + ( x(1) - x(2) ) ^ 2 ) - 1.0d0 )  ...
					 - 2.0d0 * 0.85d0 * exp( - ( x(1) + x(2) ) ^ 2 ) * ( x(1) + x(2) );
		g(2) = 0.5d0 * ( ( x(1) + x(2) )/sqrt( 1.0d0 + ( x(1) + x(2) ) ^ 2 ) - ...
					( x(1) - x(2) )/sqrt( 1.0d0 + ( x(1) - x(2) ) ^ 2 ) + 1.0d0 )  ...
					 - 2.0d0 * 0.85d0 * ( x(1) + x(2) ) * exp( - ( x(1) + x(2) ) ^ 2 );
		
		return			 
	end 	
		
end	

% ----------------------------------------------------------------------

% 	SLCDT2
%	Convergence of stochastic search algorithms to finite size pareto set approximations 

if ( strcmp(problem,'SLCDT2' )  )
										
	if ( ind == 1 ) 
		g(1) = 4.0d0 * ( x(1) - 1.0d0 ) ^ 3;
		for i = 2:n
			g(i) = 2.0d0 * ( x(i) - 1.0d0 );
		end
		
		return
	end
	
	if ( ind == 2 ) 
		g(2) = 4.0d0 * ( x(2) + 1.0d0 ) ^ 3;
		
		for i = 1:n
			if ( i ~= 2 ) 
                g(i) = 2.0d0 * ( x(i) + 1.0d0 );
            end
		end				
		
		return	
	end
	
	if ( ind == 3 ) 
		g(3) = 4.0d0 * ( x(3) - 1.0d0 ) ^ 3;
		for i = 1:n
			if ( i ~= 3 ) 
                g(i) = 2.0d0 * ( x(i) - ( - 1.0d0 ) ^ (i+1) );
            end
		end
		
		return
	end
		
end			

% ----------------------------------------------------------------------

%  SP1

if ( strcmp(problem,'SP1' )  )
	
	if ( ind == 1 ) 
		g(1) = 2.0d0 * ( x(1) - 1.0d0 ) + 2.0d0 * ( x(1) - x(2) );
		g(2) = - 2.0d0 * ( x(1) - x(2) );
		
		return
	end
	
	if ( ind == 2 ) 
		g(1) = 2.0d0 * ( x(1) - x(2) );
		g(2) = 2.0d0 * ( x(2) - 3.0d0 ) - 2.0d0 * ( x(1) - x(2) );
		
		return	
	end
		
end			

% ----------------------------------------------------------------------

%  SSFYY2

if ( strcmp(problem,'SSFYY2' )  )
		
	if ( ind == 1 ) 
		g(1) = 2.0d0 * x(1) + 5.0d0 * pi * sin( x(1) * pi / 2.0d0 );
		
		return
	end
	
	if ( ind == 2 ) 
		g(1) = 2.0d0 * ( x(1) - 4.0d0 );
		
		return	
	end
		
end	

% ----------------------------------------------------------------------

%  SK1

if ( strcmp(problem,'SK1' )  )
		
	if ( ind == 1 ) 
		g(1) = 4.0d0 * x(1) ^ 3 + 9.0d0 * x(1) ^ 2 - 2.0d1 * x(1) - 1.0d1;
		
		return
	end
	
	if ( ind == 2 ) 
		g(1) = 2.0d0 * x(1) ^ 3 - 6.0d0 * x(1) ^ 2 - 2.0d1 * x(1) + 1.0d1;
		
		return	
	end
		
end	

% ----------------------------------------------------------------------

%  SK2

if ( strcmp(problem,'SK2' )  )
		
	if ( ind == 1 ) 
		g(1) = 2.0d0 * ( x(1) - 2.0d0 );
		g(2) = 2.0d0 * ( x(2) + 3.0d0 );
		g(3) = 2.0d0 * ( x(3) - 5.0d0 );
		g(4) = 2.0d0 * ( x(4) - 4.0d0 );
		
		return
	end
	
	if ( ind == 2 ) 
		faux = 1.0d0 + ( x(1) ^ 2 + x(2) ^ 2 + x(3) ^ 2 + x(4) ^ 2 ) / 1.0d2;
		
		g(1) = ( - cos( x(1) ) * faux + ( sin( x(1) ) + sin( x(2) ) ...
		+ sin( x(3) ) + sin( x(4) ) ) * x(1) / 5.0d1 ) / faux ^ 2;
		g(2) = ( - cos( x(2) ) * faux + ( sin( x(1) ) + sin( x(2) ) ...
		+ sin( x(3) ) + sin( x(4) ) ) * x(2) / 5.0d1 ) / faux ^ 2;
		g(3) = ( - cos( x(3) ) * faux + ( sin( x(1) ) + sin( x(2) ) ...
		+ sin( x(3) ) + sin( x(4) ) ) * x(3) / 5.0d1 ) / faux ^ 2;
		g(4) = ( - cos( x(4) ) * faux + ( sin( x(1) ) + sin( x(2) ) ...
		+ sin( x(3) ) + sin( x(4) ) ) * x(4) / 5.0d1 ) / faux ^ 2;
		
		return	
	end
		
end					

% ----------------------------------------------------------------------

%  TKLY1

if ( strcmp(problem,'TKLY1' )  )
		
	if ( ind == 1 ) 
		g(1)   = 1.0d0;
		g(2:4) = 0.0d0;
		
		return
	end
	
	if ( ind == 2 ) 
		A1 = ( 2.0d0 - exp( - ( ( x(2) - 0.1d0 ) / 4.0d-3 ) ^ 2 ) ...
			         - 0.8d0 * exp( - ( ( x(2) - 0.9d0 ) / 4.0d-1 ) ^ 2 ) );
		A2 = ( 2.0d0 - exp( - ( ( x(3) - 0.1d0 ) / 4.0d-3 ) ^ 2 ) ...
			         - 0.8d0 * exp( - ( ( x(3) - 0.9d0 ) / 4.0d-1 ) ^ 2 ) );
		A3 = ( 2.0d0 - exp( - ( ( x(4) - 0.1d0 ) / 4.0d-3 ) ^ 2 ) ...
			         - 0.8d0 * exp( - ( ( x(4) - 0.9d0 ) / 4.0d-1 ) ^ 2 ) );
			         
			         
		
		g(1) = - A1 * A2 * A3 / x(1) ^ 2;
		g(2) = A2 * A3 / x(1) * ( 5.0d2 * exp( - ( ( x(2) - 0.1d0 ) / 4.0d-3 ) ^ 2 ) * ( ( x(2) - 0.1d0 ) / 4.0d-3 )  ...
		+ 4.0d0 * exp( - ( ( x(2) - 0.9d0 ) / 4.0d-1 ) ^ 2 ) * ( ( x(2) - 0.9d0 ) / 4.0d-1 ) );
		g(3) = A1 * A3 / x(1) * ( 5.0d2 * exp( - ( ( x(3) - 0.1d0 ) / 4.0d-3 ) ^ 2 ) * ( ( x(3) - 0.1d0 ) / 4.0d-3 )  ...
		+ 4.0d0 * exp( - ( ( x(3) - 0.9d0 ) / 4.0d-1 ) ^ 2 ) * ( ( x(3) - 0.9d0 ) / 4.0d-1 ) );
		g(4) = A1 * A2 / x(1) * ( 5.0d2 * exp( - ( ( x(4) - 0.1d0 ) / 4.0d-3 ) ^ 2 ) * ( ( x(4) - 0.1d0 ) / 4.0d-3 )  ...
		+ 4.0d0 * exp( - ( ( x(4) - 0.9d0 ) / 4.0d-1 ) ^ 2 ) * ( ( x(4) - 0.9d0 ) / 4.0d-1 ) );
		
		
		return	
	end
		
end		

% ----------------------------------------------------------------------

% 	Toi4

if ( strcmp(problem,'Toi4' )  )
		
	if ( ind == 1 ) 
		g(1) = 2.0d0 * x(1);
		g(2) = 2.0d0 * x(2);
		g(3) = 0.0d0;
		g(4) = 0.0d0;
		
		return
	end
	
	if ( ind == 2 ) 
		g(1) = x(1) - x(2);
		g(2) = - ( x(1) - x(2) );
		g(3) = x(3) - x(4);
		g(4) = - ( x(3) - x(4) );
		
		return	
	end
		
end		

% ----------------------------------------------------------------------

% 	Toi8

if ( strcmp(problem,'Toi8' )  )

	g(1:n) = 0.0d0;
	
	if ( ind == 1 ) 
		g(1) = 4.0d0 * ( 2.0d0 * x(1) - 1.0d0 );
		
		return
	end
	
	if ( ind ~= 1 ) 
		g(ind-1) = 4.0d0 * ind * ( 2.0d0 * x(ind-1) - x(ind) );
		g(ind)   = - 2.0d0 * ind * ( 2.0d0 * x(ind-1) - x(ind) );
		
		return	
	end
		
end	

% ----------------------------------------------------------------------

% 	Toi9

if ( strcmp(problem,'Toi9' )  )

	g(1:n) = 0.0d0;
	
	if ( ind == 1 ) 
		g(1) = 4.0d0 * ( 2.0d0 * x(1) - 1.0d0 );
		g(2) = 2.0d0 * x(2);
		
		return
	end
	
	if ( ind > 1 && ind < n ) 
		g(ind-1) =  4.0d0 * ind * ( 2.0d0 * x(ind-1) - x(ind) ) ...
				  - 2.0d0 * ( ind - 1.0d0 ) * x(ind-1);
		g(ind)   = - 2.0d0 * ind * ( 2.0d0 * x(ind-1) - x(ind) ) + 2.0d0 * ind * x(ind);
		
		return	
	end
	
	if ( ind == n ) 
		g(n-1) = 4.0d0 * n * ( 2.0d0 * x(n-1) - x(n) ) - 2.0d0 * ( n - 1.0d0 ) * x(n-1);
		g(n)   = - 2.0d0 * n * ( 2.0d0 * x(n-1) - x(n) );
		
		return
	end
		
end

% ----------------------------------------------------------------------

% 	Toi10 (Rosenbrock)

if ( strcmp(problem,'Toi10' )  )
		
	g(1:n) = 0.0d0;
	
	g(ind)   = - 4.0d2 * ( x(ind+1) - x(ind) ^ 2 ) * x(ind);
	g(ind+1) = 2.0d2 * ( x(ind+1) - x(ind) ^ 2 ) + 2.0d0 * ( x(ind+1) - 1.0d0 );
	
		return
					
end			
	
% ----------------------------------------------------------------------

%  VU1

if ( strcmp(problem,'VU1' )  )
		
	if ( ind == 1 ) 
		g(1) = - 2.0d0 * x(1) / ( x(1) ^ 2 + x(2) ^ 2 + 1.0d0 ) ^ 2;
		g(2) = - 2.0d0 * x(2) / ( x(1) ^ 2 + x(2) ^ 2 + 1.0d0 ) ^ 2;
		
		return
	end
	
	if ( ind == 2 ) 
		g(1) = 2.0d0 * x(1);
		g(2) = 6.0d0 * x(2);
		
		return	
	end
		
end		

% ----------------------------------------------------------------------

%  VU2

if ( strcmp(problem,'VU2' )  )
		
	if ( ind == 1 ) 
		g(1) = 1.0d0;
		g(2) = 1.0d0;
		
		return
	end
	
	if ( ind == 2 ) 
		g(1) = 2.0d0 * x(1);
		g(2) = 2.0d0;
		
		return	
	end
		
end		

% ----------------------------------------------------------------------

%  ZDT1

if ( strcmp(problem,'ZDT1' )  )
		
	if ( ind == 1 ) 
		g(1) = 1.0d0;
		g(2:n) = 0.0d0;
		
		return
	end
	
	if ( ind == 2 ) 
		faux = 1.0d0 + 9.0d0 * sum(x(2:n)) / ( n - 1 );
		t = x(1) / faux;
		
		g(1) = -0.5d0 * t ^ (-0.5d0);
		for i = 2:n
			g(i) = 9.0d0 / ( n - 1 ) * ( 1.0d0 - sqrt(t) / 2.0d0 );
		end
		
		return	
	end
		
end		

% ----------------------------------------------------------------------

%  ZDT2

if ( strcmp(problem,'ZDT2' )  )
		
	if ( ind == 1 ) 
		g(1) = 1.0d0;
		g(2:n) = 0.0d0;
		
		return
	end
	
	if ( ind == 2 ) 
		faux = 1.0d0 + 9.0d0 * sum(x(2:n)) / ( n - 1 );
		t = x(1) / faux;
		
		g(1) = -2.0d0 * t;
		for i = 2:n
			g(i) = 9.0d0 / ( n - 1 ) * ( 1.0d0 + t ^ 2 );
		end
		
		return	
	end
		
end		

% ----------------------------------------------------------------------

%  ZDT3

if ( strcmp(problem,'ZDT3' )  )
		
	if ( ind == 1 ) 
		g(1) = 1.0d0;
		g(2:n) = 0.0d0;
		
		return
	end
	
	if ( ind == 2 ) 
		faux = 1.0d0 + 9.0d0 * sum(x(2:n)) / ( n - 1 );
		t = x(1) / faux;
		
		a = sin( 1.0d1 * pi * x(1) );
		
		g(1) = -0.5d0 * t ^ (-0.5d0) - a - 1.0d1 * pi * x(1) * cos( 1.0d1 * pi * x(1) );
		for i = 2:n
			g(i) = 9.0d0 / ( n - 1 ) * ( 1.0d0 - 0.5d0 * sqrt(t) ) ;
		end
		
		return	
	end
		
end		

% ----------------------------------------------------------------------

%  ZDT4

if ( strcmp(problem,'ZDT4' )  )
		
	if ( ind == 1 ) 
		g(1) = 1.0d0;
		g(2:n) = 0.0d0;
		
		return
	end
	
	if ( ind == 2 ) 
		faux = 0.0d0;
		for i = 2:n
			faux = faux + x(i) ^ 2 - 1.0d1 * cos( 4.0d0 * pi * x(i) );
		end
		faux = faux + 1.0d0 + 1.0d1 * ( n - 1 );
		t = x(1) / faux;
		
		g(1) = -0.5d0 * t ^ (-0.5d0);
		for i = 2:n
			a = sin( 4.0d0 * pi * x(i) );
			g(i) = ( 2.0d0 * x(i) + 4.0d1 * pi * a ) * ( 1.0d0 - sqrt(t) / 2.0d0 );
		end
		
		return	
	end
		
end		

% ----------------------------------------------------------------------

%  ZDT6

if ( strcmp(problem,'ZDT6' )  )
		
	if ( ind == 1 ) 			
		a = exp( -4.0d0 * x(1) );
		b = sin( 6.0d0 * pi * x(1) );
		
		g(1) = 4.0d0 * a * b ^ 6 - 3.6d1 * pi * a * b ^ 5 * cos( 6.0d0 * pi * x(1) );
		g(2:n) = 0.0d0;
		
		return
	end
	
	if ( ind == 2 ) 			
		a = exp( -4.0d0 * x(1) );
		b = sin( 6.0d0 * pi * x(1) );
		gaux1 = 1.0d0 - a * b ^ 6 ;
		gaux2 = 1.0d0 + 9.0d0 * ( sum(x(2:n)) / ( n - 1 ) ) ^ 0.25d0;
		t = gaux1 / gaux2;
		
		A1 = 9.0d0 * 0.25d0 * ( sum(x(2:n)) / ( n - 1 ) ) ^ (-0.75d0) / ( n - 1 );
		
		g(1) = -2.0d0 * t * ( 4.0d0 * a * b ^ 6 - 3.6d1 * pi * a * b ^ 5 * cos( 6.0d0 * pi * x(1) ) );
		for i = 2:n
			g(i) = A1 * ( 1.0d0 + t ^ 2 );
		end
		
		
		return	
	end
		
end		

% ----------------------------------------------------------------------

%  ZLT1

if ( strcmp(problem,'ZLT1' )  )
	
	g(ind) = 2.0d0 * ( x(ind) - 1.0d0 );
	for i = 1:n
		if ( i ~= ind ) 
            g(i) = 2.0d0 * x(i);
        end
	end
		
 end
