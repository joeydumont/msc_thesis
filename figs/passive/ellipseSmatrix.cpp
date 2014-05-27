/** -------------------- Information ------------------ #
# Author:       Joey Dumont <joey.dumont@gmail.com>     #
# Date created: May 27th, 2014							#
# Date mod.:    May 27th, 2014							#
# Description:  We compute the scattering matrix of the #
#				quadrupolar cavity and show its banded	#
#				form. 									#
# ----------------------------------------------------**/

#include <libRZ>
#include <ctime>
#include <iostream>
#include <sstream>

/*! \class refIndexElliptical
 *
 *  \brief We define the refractive index of the elliptical cavity.
 */
class refIndexEllipse : public libRZ::RefractiveIndex
{
public:
    refIndexEllipse(double radius,
                                 double epsilon,
                                 std::complex<double> _nIn,
                                 std::complex<double> _no)
{
    // Geometric parameters of the ellipse
    r = radius;
    eps = epsilon;
    a = r*(1.0+eps);
    b = r/(1.0+eps);

    // Refractive indices
    nIn = _nIn;
    no = _no;

    // Inner and outer radii
    r_0 = b;
    rMax = a;
}

std::complex<double> operator ()(double radius, double angle)
{
    double position = radius-a*b/sqrt(pow(b*cos(angle),2.0)+pow(a*sin(angle),2.0));

    if (position < 0.0)
        return nIn;
    else
        return no;
}

protected:
    // Radius of the base circle and deformation parameter.
    double r,eps,a,b;
};

int main()
{
	// Definition of the Ellipse class and of the problem. 
	refIndexEllipse elli(1.0/1.36747,.36747,std::complex<double>(3.3,0.0),std::complex<double>(1.0,0.0));
	libRZ::RZmethod problem(elli, std::complex<double>(5.0,0.0), 2000, 1.0);

	// computation of the scattering matrix.
	arma::cx_mat sMatrix = problem.solve();

	// Save the magnitude of the scattering matrix elements.
	arma::abs(sMatrix).eval().save("sMatrixEllipseMag.dat", arma::raw_ascii);

	return 0;
}