function [v, blockVar, blockMean] = blockAverage(datastream, maxBlockSize, isplot)
 
  Nobs         = length(datastream);           % total number of observations in datastream
  minBlockSize = 1;                           % min: 1 observation/block

  if(nargin < 2)
    maxBlockSize = floor(Nobs/4);             % max: 4 blocs (otherwise can't calc variance)
  end
  
  if(nargin < 3) 
    isplot = 1;
  end



  NumBlocks = maxBlockSize - minBlockSize + 1;   % total number of block sizes

  blockMean = zeros(NumBlocks,1);                   % mean (expect to be "nearly" constant)
  blockVar  = zeros(NumBlocks,1);                   % variance associated with each blockSize
  blockCtr  = 1;
				%
				%  blockSize is # observations/block
				%  run them through all the possibilities
				%

  for blockSize = minBlockSize:maxBlockSize

    Nblock    = floor(Nobs/blockSize);              % total number of such blocks in datastream
    obsProp   = zeros(Nblock,1);                    % container for parcelling block 
 				
				% Loop to chop datastream into blocks
				% and take average
    for i = 1:Nblock
      ibeg = (i-1)*blockSize + 1;
      iend =  ibeg + blockSize - 1;
      obsProp(i) = mean(datastream(ibeg:iend));
    end

    blockMean(blockCtr) = mean(obsProp);
    blockVar(blockCtr)  = var(obsProp)/(Nblock - 1);
    blockCtr = blockCtr + 1;

  end

  v = minBlockSize:maxBlockSize;

  if(isplot)

    h = figure;

    subplot(2,1,1)
    plot(v, sqrt(blockVar),'ro-','LineWidth',2)
    xlabel('block size')
    ylabel('std')

    subplot(2,1,2)
    errorbar(v,blockMean, sqrt(blockVar))
    ylabel('<x>')
    xlabel('block size')

%     printf("<x> = %f +/- %f\n", blockMean(end), sqrt(blockVar(end)));

% 
% Plot png!
%
%     W = 4; H = 3;
%     set(h,'PaperUnits','inches')
%     set(h,'PaperSize',[H,W])
%     set(h,'PaperPosition',[0,0,W,H])
% 
%     FN = findall(h,'-property','FontName');
%     set(FN,'FontName','HelveticaBold');
% 
%     FS = findall(h,'-property','FontSize');
%     set(FS,'FontSize',8);
% 
%     print('tmp.png','-dpng')

  end


end
