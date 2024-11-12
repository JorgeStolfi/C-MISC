#! /bin/sed -f
# Last edited on 2005-10-02 20:04:55 by stolfi

s/[.]cholesky/.mth == SPSys_LM_Cholesky/g
s/[.]gaussLU/.mth == SPSys_LM_GaussLU/g
s/[.]gaussSeidel/.mth == SPSys_LM_GaussSeidel/g
s/[.]conjGrad/.mth == SPSys_LM_ConjGrad/g

s/[>]cholesky/->mth == SPSys_LM_Cholesky/g
s/[>]gaussLU/->mth == SPSys_LM_GaussLU/g
s/[>]gaussSeidel/->mth == SPSys_LM_GaussSeidel/g
s/[>]conjGrad/->mth == SPSys_LM_ConjGrad/g
