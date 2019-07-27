#!/bin/bash

for i in {1..8}
do
  ./fitLogisticVariableSigma sample num_samples=1000 num_warmup=1000 \
  adapt delta=0.99 random seed=12345 \
  id=$i data file=fit_growth_data.Rstan \
  output file=LogisticVariableSigma_samples$i.csv &
done
