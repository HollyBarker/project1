#include<iostream>
#include<cmath>
#include<algorithm>

#ifndef MVECTOR_H // the 'include guard'
#define MVECTOR_H // see C++ Primer Sec. 2.9.2

#include <vector>

// Class that represents a mathematical vector
class MVector
{

//Operator overload for <<
friend std::ostream& operator<< (std::ostream& os, const MVector& output)
{
	os<<"(";
	for (int i=0;i<output.size()-1;i++)
	{
		os<<output[i]<<",";
	}
	os<<output[output.size()-1]<<")";
	return os;
}

public:
	// constructors
	MVector() {}
	explicit MVector(int n) : v(n) {}
	//the explicit keyword prevents implicit conversions from type integer to MVector. without the explicit keyword (and using the commented out section at the bottom) then a is a length 2 vector when added
	//to y because of the implicit conversion. Adding the two does not change the first 2 cpts, so z[0]=y[0], and z[1]=y[1] i.e. where there is an a[i] value, but  z[2] is very strange.
	MVector(int n, double x) : v(n, x) {}

	// access element (lvalue) (see example sheet 5, q5.6)
	double &operator[](int index) { return v[index]; }

	// access element (rvalue) (see example sheet 5, q5.7)
	double operator[](int index) const { return v[index]; }

	int size() const { return v.size(); } // number of elements

	double LInfNorm() const
	{
		double maxAbs = 0;
		std::size_t s = size();
		for (int i=0; i<s; i++)
		{
			maxAbs = std::max(std::abs(v[i]), maxAbs);
		}
		return maxAbs;
	}
	
	double L2Norm() const
	{
		double rootsumsquares=0;
		for (int i=0;i<size();i++)
		{
			rootsumsquares+= v[i]*v[i];
		}
		rootsumsquares=pow(rootsumsquares,0.5);
		return rootsumsquares;
	}

private:
	std::vector<double> v;
};

// Operator overload for "scalar * vector"
inline MVector operator*(const double& lhs,const MVector& rhs)
{
	MVector temp(rhs);
	for (int i=0;i<temp.size();i++) temp[i]*=lhs;
	return temp;
}
// Operator overload for "vector * scalar"
inline MVector operator*(const MVector& lhs,const double& rhs)
{
	MVector temp(lhs);
	for (int i=0;i<temp.size();i++) temp[i]*=rhs;
	return temp;
}
//Operator overload for "vector + vector"
inline MVector operator+ (const MVector& lhs, const MVector& rhs)
{
	if (lhs.size()==rhs.size())
	{
		MVector temp(lhs);
		for (int i=0; i<temp.size();i++) temp[i]+=rhs[i];
		return temp;
	}
	else 
	{
		std::cout<<"ERROR: Attempted addition of two vectors of different length."<<std::endl;
		exit(1);
	}
} 
//Operator overload for "vector - vector"
inline MVector operator- (const MVector& lhs, const MVector& rhs)
{
	if (lhs.size()==rhs.size())
	{
		MVector temp(lhs);
		for (int i=0; i<temp.size();i++) temp[i]-=rhs[i];
		return temp;
	}
	else
	{
		std::cout<<"ERROR: Attempted subtraction of two vectors of different length."<<std::endl;
		exit(1);
	}
}
//Operator overload for "vector / scalar"
inline MVector operator/ (const MVector& lhs, const double& rhs)
{
	MVector temp(lhs);
	for (int i=0;i<temp.size();i++) temp[i]/=rhs;
	return temp;
}

double dot(const MVector& lhs, const MVector& rhs)
{
	double dotproduct=0;
	if (lhs.size()==rhs.size())
	{
		for (int i=0; i<lhs.size(); i++)
		{
			dotproduct+=lhs[i]*rhs[i];
		}
		return dotproduct;
	}
	else 
	{
		std::cout<<"ERROR: Attempted dot product of two vectors of different length."<<std::endl;
		exit(1);
	}
}
#endif

#ifndef MMATRIX_H // the 'include guard'
#define MMATRIX_H

#include <vector>
#include <iostream>

// Class that represents a mathematical matrix
class MMatrix
{

friend std::ostream& operator<< (std::ostream& os, MMatrix A)
{
	for (int i=0;i<A.Rows();i++)
	{
		for (int j=0;j<A.Cols();j++)
		{
			os.precision(5);
			os.width(10);os<<A(i,j);
		}
		os<<std::endl;
	}
	return os;
}

public:
	// constructors
	MMatrix() : nRows(0), nCols(0) {}
	MMatrix(int n, int m, double x = 0) : nRows(n), nCols(m), A(n * m, x) {}

	// set all matrix entries equal to a double
	MMatrix &operator=(double x)
	{
		for (int i = 0; i < nRows * nCols; i++) A[i] = x;
		return *this;
	}

	// access element, indexed by (row, column) [rvalue]
	double operator()(int i, int j) const
	{
		return A[j + i * nCols];
	}

	// access element, indexed by (row, column) [lvalue]
	double &operator()(int i, int j)
	{
		return A[j + i * nCols];
	}

	// size of matrix
	int Rows() const { return nRows; }
	int Cols() const { return nCols; }


private:
	unsigned int nRows, nCols;
	std::vector<double> A;
};
MVector operator*(const MMatrix& A, const MVector& x)
{
	MVector b(A.Rows());
	if (A.Cols() == x.size())
	{
		double sum;
		for (int i=0;i<A.Rows();i++)
		{
			sum=0;
			for (int j=0;j<A.Cols();j++)
			{
				sum+=A(i,j)*x[j];
			}
			b[i]=sum;
		}
		return b;
	}
	else 
	{
		std::cout<<"ERROR: Attempted matrix * vector multiplication where number of matrix columns != vector length."<<std::endl;
		exit(1);
	}
}
#endif


int main()
{
	MMatrix M(4,3);
	M(0,0)=5;
	M(3,2)=2;
	std::cout<<M(3,2)<<std::endl;
	MMatrix A(4,3);
	for (int i=0; i<A.Rows(); i++)
	{
		for(int j=0;j<A.Cols();j++)
		{
			A(i,j)=3.0*(i+1)+(j+1);
 
		}
	}
	std::cout<<A<<std::endl;
	
	MVector x(3);
	x[0]=0.5; x[1]=1.6; x[2]=3.2;
	
	MVector b(4);

	b=A*x;

	std::cout<<b<<std::endl;
	return 0;
}
