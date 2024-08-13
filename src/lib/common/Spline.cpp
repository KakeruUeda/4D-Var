#include "Spline.h"

/*************************************************************
 * @brief Compute coefficients of the spline
 * @ref   https://qiita.com/khidaka/items/84610cd890ecb8443d96
 */
std::vector<Spline::Coefficients> Spline::compCoefficients(const std::vector<double> &x,
                                                           const std::vector<double> &y)
{
	int n=x.size()-1;
	std::vector<double> h(n);
	for(int i=0; i<n; ++i)
	{
		h[i]=x[i+1]-x[i];
	}

	MatrixXd A=MatrixXd::Zero(4*n,4*n);
	VectorXd b=VectorXd::Zero(4*n);

	int row=0;

	for(int i=0; i<n; ++i)
	{
		A(row,4*i)=1;
		b(row)=y[i];
		row++;
		A(row,4*i)=1;
		A(row,4*i+1)=h[i];
		A(row,4*i+2)=h[i]*h[i];
		A(row,4*i+3)=h[i]*h[i]*h[i];
		b(row)=y[i+1];
		row++;
	}

	for(int i=1; i<n; ++i)
	{
		A(row,4*(i-1)+1)=1;
		A(row,4*(i-1)+2)=2*h[i-1];
		A(row,4*(i-1)+3)=3*h[i-1]*h[i-1];
		A(row,4*i+1)=-1;
		row++;
	}

	for(int i=1; i<n; ++i)
	{
		A(row,4*(i-1)+2)=2;
		A(row,4*(i-1)+3)=6*h[i-1];
		A(row,4*i+2)=-2;
		row++;
	}

	A(row,2)=2;
	row++;
	A(row,4*(n-1)+2)=2;
	A(row,4*(n-1)+3)=6*h[n-1];
	row++;

	VectorXd coeffs=A.lu().solve(b);

	std::vector<Spline::Coefficients> coefficients(n);
	for(int i=0; i<n; ++i)
	{
		coefficients[i].a=coeffs(4*i);
		coefficients[i].b=coeffs(4*i+1);
		coefficients[i].c=coeffs(4*i+2);
		coefficients[i].d=coeffs(4*i+3);
		coefficients[i].x=x[i];
	}

	return coefficients;
}

/*************************************************************
 * @brief Evaluate the spline
 * @ref   https://qiita.com/khidaka/items/84610cd890ecb8443d96
 */
double Spline::evaluate(const std::vector<Coefficients> &coefficients,double x)
{
	if(coefficients.empty())
	{
		PetscPrintf(PETSC_COMM_WORLD,"Spline coefficients are empty.\n");
	}

	const Coefficients *spline=nullptr;
	for(const auto &s : coefficients)
	{
		if(x>=s.x)
		{
			spline=&s;
		}
		else
		{
			break;
		}
	}

	if(!spline)
	{
		PetscPrintf(PETSC_COMM_WORLD,"x is out of the range of the spline.\n");
	}

	double dx=x-spline->x;
	return spline->a+spline->b*dx+spline->c*dx*dx+spline->d*dx*dx*dx;
}