// Implements anisotropic pvd equations
// Vectors which define the directions in which the material is anisotropic are passed to the constitutive
//  law as a vector of vectors.
// These vectors are computed as a function of local and global coordinates at the integral points.

// Essentially just a copy and paste of the PVD equations, just with: additional arguments when
//  calculating strain, residual and jacobian contributions