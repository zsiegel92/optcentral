extraTimes = xlTimes(end) + trainForesight; 
   %skuInds = (1+nTest*(i-1):(nTest*i))
   extraTS = realPreds(:,nTest*(1:nz))'; %The forecast-series calculated at the final test sample of each product.
   