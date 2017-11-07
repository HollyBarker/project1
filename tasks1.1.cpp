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




int main()
{
	MVector v(3), w(3), x(3), u(3), y(3), z(3), b(2);
	v[0]=0.1; v[1]=4.8; v[2]=3.7;
	w[0]=3.1; w[1]=8.5; w[2]=3.6;
	x[0]=5.8; x[1]=7.4; x[2]=12.4;
	b[0]=5;b[1]=12;
	//MVector c=v+b;
	//MVector d=y-b;
	//doing the addition throws the error as expected. doing the subtraction throws the error and some other long message to the terminal
	u=4.7*v+1.3*w-6.7*x;
	std::cout<<u[0]<<','<<u[1]<<','<<u[2]<<std::endl;
	y[0]=2; y[1]=1; y[2]=3;
	
	/*double a=2;
	z=y+a;
	std::cout<<z[0]<<','<<z[1]<<','<<z[2]<<std::endl;*/
	/*std::cout<<u<<std::endl;
	std::cout<<y.LInfNorm()<<std::endl;
	std::cout<<b.L2Norm()<<std::endl;
	std::cout<<dot(v,w)<<std::endl;*/
	MVector U(3), V(3), W(3);
	U[0]=1.5; U[1]=1.3; U[2]=2.8;
	V[0]=6.5; V[1]=2.7; V[2]=2.9;
	W[0]=0.1; W[1]=-7.2; W[2]=3.4;
	double alpha=(dot(U,U)/dot(V,W));
	/*std::cout<<alpha<<std::endl;	
	std::cout<<U.L2Norm()<<std::endl;
	std::cout<<V.L2Norm()<<std::endl;
	std::cout<<W.L2Norm()<<std::endl;*/
	return 0;
}
