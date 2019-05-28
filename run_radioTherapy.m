
function run_radioTherapy( varargin )
  close all; rng(1);
  % run_radioTherapy( A, d, nVoxsPerStructure, w, w0Minus )

  %dataDir = '/Volumes/NDWORK128GB/radioTherapyData/Head-and-Neck';
  dataDir = '/Users/nicholasdwork/.DellEMS/Vols/disk2s1/radioTherapyData/Head-and-Neck';

  datacase = 1;
  stepSizes = [ 1d-8 1d-7 1d-6 1d-5 1d-4 1d-3 1d-2 5d-2 1d-2 0.02 0.1 0.2 ...
    0.5 1.0 1.5 2.0 3.0 4.0 5.0 1d1 1d2 1d3 ];
  batchSizes = [ 8, 16, 32, 64, 128, 256 ];
  algs = { 'proxSVRG', 'fista', 'stochProxGrad' };
  N = 100;

  p = inputParser;
  p.addOptional( 'A', [], @isnumeric );
  p.addOptional( 'd', [], @isnumeric );
  p.addOptional( 'nVoxsPerStructure', [], @isnumeric );
  p.addOptional( 'w', [], @isnumeric );
  p.addOptional( 'w0', [], @isnumeric );
  p.parse( varargin{:} );
  A = p.Results.A;
  d = p.Results.d;
  nVoxsPerStructure = p.Results.nVoxsPerStructure;
  w = p.Results.w;
  w0 = p.Results.w0;

%   if nargin < 1
%     [A,d,nVoxsPerStructure,w,w0] = loadDatacase( datacase, dataDir );
%   end

  nBatchSizes = numel( batchSizes );
  nAlgs = numel( algs );
  nStepSizes = numel( stepSizes );
  objectiveValues = cell( nBatchSizes, nAlgs, nStepSizes );

  nIterations = nBatchSizes * nAlgs * nStepSizes + nAlgs * nStepSizes + nStepSizes;
  
  p = parforProgress( nIterations );  % Note: if tmpFile is not supplied

  for batchIndx = 1 : nBatchSizes
    batchSize = batchSizes( batchIndx );
    batchFile = ['batch', num2str(batchSize), '.mat'];
    if ~exist( batchFile, 'file' )

      for algIndx = 1 : nAlgs
        theseObjValues = cell( 1, nStepSizes );
        alg = algs{ algIndx };

        if strcmp( alg, 'fista' ) && exist( 'fistaObjValues', 'var' )
          theseObjValues = fistaObjValues;

        else

          parfor stepSizeIndx = 1:nStepSizes
            msgHdr = ['batch/alg/step ', num2str(batchIndx), '/', num2str(algIndx), ...
              '/', num2str(stepSizeIndx), ' of ', num2str(nBatchSizes), '/', num2str(nAlgs), ...
              '/', num2str(nStepSizes), '  ' ];
            p.progress( stepSizeIndx + (algIndx-1)*nStepSizes + ...
              (batchIndx-1)*nAlgs*nStepSizes, 'msgHdr', msgHdr );   %#ok<PFBNS>
            
            stepSize = stepSizes( stepSizeIndx );
            if strcmp( alg, 'stochProxGrad' )
              [~,oValues] = radioTherapy( A, d, nVoxsPerStructure, w, w0, alg, ...
                't', stepSize, 'N', N, 'batchSize', batchSize );
            elseif strcmp( alg, 'proxSVRG' )
              [~,oValues] = radioTherapy( A, d, nVoxsPerStructure, w, w0, alg, ...
                't', stepSize, 'N', N/2, 'batchSize', batchSize );
              oValues = upsample( oValues, 2 ) + circshift( upsample( oValues, 2 ), 1 );
            else
              [~,oValues] = radioTherapy( A, d, nVoxsPerStructure, w, w0, alg, ...
                't', stepSize, 'N', N );
            end

            theseObjValues{ stepSizeIndx } = oValues;
          end

          if strcmp( alg, 'fista' )
            fistaObjValues = theseObjValues;
          end

        end

        for stepSizeIndx = 1:nStepSizes
          objectiveValues{ batchIndx, algIndx, stepSizeIndx } = theseObjValues{ stepSizeIndx };
        end
      end

      batchObjectiveValues = cell( nAlgs, nStepSizes );
      for algIndx = 1:nAlgs
        for stepSizeIndx = 1:nStepSizes
          batchObjectiveValues{ algIndx, stepSizeIndx } = ...
            objectiveValues{ batchIndx, algIndx, stepSizeIndx };
        end
      end
      save( batchFile, 'batchObjectiveValues' ); 
    else
      load( batchFile );

      for algIndx = 1:nAlgs
        for stepSizeIndx = 1:nStepSizes
          objectiveValues{ batchIndx, algIndx, stepSizeIndx } = ...
            batchObjectiveValues{ algIndx, stepSizeIndx };
        end
      end
    end
  end
  p.clean;

  fistaIndx = find( strcmp( algs, 'fista' ) );
  stochProxGradIndx = find( strcmp( algs, 'stochProxGrad' ) );
  proxSVRGIndx = find( strcmp( algs, 'proxSVRG' ) );

  objValues_fista = squeeze( objectiveValues(1,fistaIndx,:) )';   %#ok<FNDSB>
  objValues_stochProxGrad = squeeze( objectiveValues(:,stochProxGradIndx,:) );   %#ok<FNDSB>
  objValues_proxSVRG = cell2mat( objectiveValues(:,proxSVRGIndx,:) );   %#ok<FNDSB>

  stepSizeNames = cell( 1, nStepSizes );
  for i=1:nStepSizes, stepSizeNames{i}=num2str(stepSizes(i)); end
  figure; semilogynice( objValues_fista );  titlenice( 'fista' );
  legend( stepSizeNames );
  for batchIndx = 1 : nBatchSizes
    figure; semilogynice( squeeze( objValues_stochProxGrad(batchIndx,:) ) );
    titlenice( ['stochProxGrad ',num2str(batchSize)] );
    legend( stepSizeNames );
    figure; semilogynice( squeeze( objValues_proxSVRG(batchIndx,:) ) );
    semilogynice( objValues_proxSVRG );
    titlenice( ['proxSVRG ',num2str(batchSize)] );
    legend( stepSizeNames );
  end

  bestObjValues_fista = min( objValues_fista );
  bestObjValues_stochProxGrad = min( objValues_stochProxGrad );
  bestObjValues_proxSVRG = min( objValues_proxSVRG );

  [bestObjValue_fista,bestStepSizeIndx_fista] = min( bestObjValues_fista );
  bestStepSize_fista = stepSizes( bestStepSizeIndx_fista );
  disp([ 'Best objective value with FISTA: ', num2str(bestObjValue_fista) ]);
  disp([ 'Best step size for FISTA: ', num2str(bestStepSize_fista) ]);

  [bestObjValue_stochProxGrad,bestStepSizeIndx_stochProxGrad] = min( bestObjValues_stochProxGrad );
  bestStepSize_stochProxGrad = stepSizes( bestStepSizeIndx_stochProxGrad );
  disp([ 'Best objective value with FISTA: ', num2str(bestObjValue_stochProxGrad) ]);
  disp([ 'Best step size for FISTA: ', num2str(bestStepSize_stochProxGrad) ]);

  [bestObjValue_proxSVRG,bestStepSizeIndx_proxSVRG] = min( bestObjValues_proxSVRG );
  bestStepSize_proxSVRG = stepSizes( bestStepSizeIndx_proxSVRG );
  disp([ 'Best objective value with FISTA: ', num2str(bestObjValue_proxSVRG) ]);
  disp([ 'Best step size for FISTA: ', num2str(bestStepSize_proxSVRG) ]);

  bestObjValues_fista = objValues_fista( :, bestStepSizeIndx_fista );
  bestObjValues_stochProxGrad = objValues_stochProxGrad( :, bestStepSizeIndx_stochProxGrad );
  bestObjValues_proxSVRG = objValues_proxSVRG( :, bestStepSizeIndx_proxSVRG );
  figure; semilogynice( bestObjValues_fista );
  hold on; semilogynice( bestObjValues_stochProxGrad );
  hold on; semilogynice( bestObjValues_proxSVRG );
  legend( [ 'fista - step size ', num2str(bestStepSize_fista) ], ...
          [ 'stochProxGrad (batch ', num2str(batchSize), ') - step size ', ...
            num2str(bestStepSize_stochProxGrad) ], ...
          [ 'proxSVRG (batch ', num2str(batchSize), ') - step size ', ...
            num2str(bestStepSize_proxSVRG) ] );
  titlenice( 'Performance with best step sizes' );

end
