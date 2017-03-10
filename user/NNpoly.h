#ifndef NNPOLY_H
#define NNPOLY_H
#include<iostream>
namespace antok {

	namespace user {

		namespace hubers {

			namespace NNpoly
			{
				double* getParams2009();
				double* getParams2012();

				double Ebeam(double *x, const double *p);
			}

		}
	}
}

#endif
