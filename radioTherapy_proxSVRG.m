
function [xStar,objectiveValues] = radioTherapy_proxSVRG( ...
  A, d, nVoxsPerStructure, w, w0, varargin )

  % The algorithm this code is based on is described in "A Proximal Stochastic
  % Gradient Method with Progressive Variance Reduction" by Xiao and Zhang

  p = inputParser;
  p.addParameter( 't', 0.001, @isnumeric );  % step size
  p.addParameter( 'batchSize', 16, @isnumeric );
  p.addParameter( 'N', 100, @isnumeric );  % Number of iterations
  p.addParameter( 'verbose', true, @(x) islogical(x) || isnumeric(x) );
  p.parse( varargin{:} );
  t = p.Results.t;
  batchSize = p.Results.batchSize;
  N = p.Results.N;
  verbose = p.Results.verbose;

  yTilde = zeros( size( A, 2 ), 1 );

  calculateObjectiveValues = 0;
  if nargout > 1
    objectiveValues = zeros(N,1);
    calculateObjectiveValues = 1;
  end

  AT = A';
  A0 = A( 1 : nVoxsPerStructure(1), : );
  A0T = A0';
  d0 = d( 1 : nVoxsPerStructure(1) );

  function out = g( Ax )
    tmp0 = max( d0 - Ax(1:nVoxsPerStructure(1)), 0 );
    out = 0.5 * sum( w0 .* tmp0 .* tmp0 );
    tmp = max( Ax - d, 0 );
    out = out + 0.5 * sum( w .* tmp .* tmp );
  end

  nRowsA = size( A, 1 );
  nRowsA0 = size( A0, 1 );
  nTermsObjective = nRowsA + nRowsA0;

  flagsA0 = zeros( nTermsObjective, 1 );
  flagsA0( 1 : nRowsA0 ) = 1;

  nBatches = ceil( nTermsObjective / batchSize );


  for s=0:N-1

    % Full gradient calculation
    A_yTilde = AT' * yTilde;
    A0_yTilde = A_yTilde(1:nVoxsPerStructure(1));
    gPrime_yTilde = A' * ( w .* max( A_yTilde - d, 0 ) ) - ...
      A0' * ( w0 .* max( d0 - A0_yTilde, 0 ) );

    if calculateObjectiveValues > 0
      % Calculate objective on y
      if min( yTilde ) < 0, hOfy=Inf; else, hOfy=0; end
      objectiveValues(s+1) = g( A_yTilde ) + hOfy;
    end

    permIndxs = randperm( nTermsObjective )';
    permA0Flags = flagsA0( permIndxs );

    y = yTilde;
    yTilde(:) = 0;
    for batchIndx = 1 : nBatches
      currentIndx = (batchIndx-1) * batchSize + 1;
      endIndx = min( currentIndx+batchSize-1, nTermsObjective );

      theseIndxs = permIndxs( currentIndx : endIndx );
      theseAIndxs = theseIndxs( theseIndxs > nRowsA0 ) - nRowsA0;
      theseA0Flags = permA0Flags( currentIndx : endIndx );

      gPrimeHatY = 0;
      gPrimeHat_yTilde = 0;

      if numel( theseAIndxs ) > 0
        subAT = AT( :, theseAIndxs );
        subA = subAT';
        subd = d( theseAIndxs );
        subw = w( theseAIndxs );

        subAy = subAT' * y;
        gPrimeHatY = gPrimeHatY + subA' * ( subw .* max( subAy - subd, 0 ) );

        subA_yTilde = A_yTilde( theseAIndxs );
        gPrimeHat_yTilde = gPrimeHat_yTilde + subA' * ( subw .* max( subA_yTilde - subd, 0 ) );
      end

      if max( theseA0Flags ) > 0
        theseA0Indxs = theseIndxs( theseA0Flags == 1 );
        subA0T = A0T( :, theseA0Indxs );
        ptvd0 = d0( theseA0Indxs );
        ptvw0 = w0( theseA0Indxs );

        subA0_y = subA0T' * y;
        gPrimeHatY = gPrimeHatY - subA0T * ( ptvw0 .* max( ptvd0 - subA0_y, 0 ) );

        subA0_yTilde = A0_yTilde( theseA0Indxs );
        gPrimeHat_yTilde = gPrimeHat_yTilde - subA0T * ...
          ( ptvw0 .* max( ptvd0 - subA0_yTilde, 0 ) );
      end

      thisBatchSize = numel( theseIndxs );
      batchScaling = nTermsObjective / thisBatchSize;
      v = ( gPrimeHatY - gPrimeHat_yTilde ) * batchScaling + gPrime_yTilde;
      x = y - t * v;
      y = max( x, 0 );  % proximal operator of h

      yTilde = yTilde + y;
    end
    %yTilde = yTilde / nBatches;
yTilde = y;  % Note, this is different from the paper

    if verbose
      formatString = ['%', num2str(ceil(log10(N))), '.', num2str(ceil(log10(N))), 'i' ];
      verboseString = [ 'proxSVRG Iteration: ', num2str(s,formatString) ];
      if calculateObjectiveValues > 0
        verboseString = [ verboseString, ',  objective: ', num2str( objectiveValues(s+1) ) ];   %#ok<AGROW>
      end
      disp( verboseString );
    end
  end

  xStar = y;
end
