
function [xStar,objectiveValues] = radioTherapy_fista_wLS( A, d, nVoxsPerStructure, w, w0Minus )

  AT = A';
  A0T = AT( :, 1 : nVoxsPerStructure(1) );
  d0 = d( 1 : nVoxsPerStructure(1) );

  function out = g( x )
    Ax = AT' * x;
    tmp0 = max( d0 - Ax(1:nVoxsPerStructure(1)), 0 );
    out = 0.5 * sum( w0Minus .* tmp0 .* tmp0 );
    tmp = max( Ax - d, 0 );
    out = out + 0.5 * sum( w .* tmp .* tmp );
  end

  function out = gPrime( x )
    Ax = AT' * x;
    out = -A0T * ( w0Minus .* max( d0 - Ax(1:nVoxsPerStructure(1)), 0 ) );
    out = out + AT * ( w .* max( Ax - d, 0 ) );
  end

  function out = h( x )
    out=0;
    if min( x ) < 0, out = Inf; end
  end

  proxth = @(x,t) max( x, 0 );


  x0 = zeros( size( A, 2 ), 1 );
  if nargout > 1
    [xStar,objectiveValues] = fista_wLS( x0, @g, @gPrime, proxth, 'h', @h, 'verbose', true );
    %[xStar,objectiveValues] = fista( x0, @g, @gPrime, proxth, 't', 0.001, 'h', @h, 'verbose', true );
  else
    xStar = fista_wLS( x0, @g, @gPrime, proxth, 'h', @h, 'verbose', true );
    %xStar = fista( x0, @g, @gPrime, proxth, 't', 0.001, 'h', @h, 'verbose', true );
  end
end
