
function makeRadPlan( A, d, nVoxsPerStructure, w, w0Minus )

  A0 = A(1:nVoxsPerStructure(1),:);
  d0 = d(1:nVoxsPerStructure(1));

  function out = g( x )
    Ax = A * x;
    tmp0 = max( d0 - Ax(1:nVoxsPerStructure(1)), 0 );
    out = 0.5 * sum( w0Minus .* tmp0 .* tmp0 );
    tmp = max( Ax - d, 0 );
    out = out + 0.5 * sum( w .* tmp .* tmp );
  end

  function out = gPrime( x )
    Ax = A * x;
    out = -A0' * ( w0Minus .* max( d0 - Ax(1:nVoxsPerStructure(1)), 0 ) );
    out = out + A' * ( w .* max( Ax - d, 0 ) );
  end

  function out = h( x )
    out=0;
    if min( x ) < 0, out = Inf; end
  end

  proxth = @(x,t) max( x, 0 );


  x0 = zeros( size( A, 2 ), 1 );
  [xStar,objectiveValuesLS] = fista_wLS( x0, @g, @gPrime, proxth, 'h', @h, 'verbose', true );

  figure; semilogynice( objectiveValues ); hold all; semilogynice( objectiveValuesLS );
  legendnice( 'fista', 'wLS' );

end
