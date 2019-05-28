
function [A,d,nVoxsPerStructure,w,w0Minus] = loadDatacase( datacase, dataDir )

  switch datacase
    case 1
      load( [ dataDir, '/', 'Head-and-Neck_01.mat' ] );
      rxDose = 46;  % Gray
      structureIndxs = [ 1,2,3,12,13,14,16,18,20,22,25,27,29,31,33,35];
        % Note: structureIndxs(1) is the PTV
      ptvWeight = 100;
      oarWeight = 10;
  end

  A = [];
  nVoxsPerStructure = zeros( numel(structureIndxs), 1 );
  for i = 1 : numel( structureIndxs )
    structureIndx = structureIndxs( i );

    thisA = data.matrix(structureIndx).A;
    A = [ A; thisA; ];  %#ok<AGROW>
    
    nVoxsPerStructure(i) = size( thisA, 1 );
  end

  d = zeros( size( A, 1 ), 1 );
  ptvA = data.matrix( structureIndxs(1) ).A;
  d(1:size(ptvA,1)) = rxDose;
  w0Minus = ptvWeight * ones( size(ptvA,1), 1 );
  w(1:size(ptvA,1)) = ptvWeight;
  w = oarWeight * ones( size( A, 1 ), 1 );

end
