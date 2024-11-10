#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;

// Count the numbers of points per type around reference points.
struct CountNbdWrkr : public Worker
{
  // source vectors
  const RVector<double> r2;
  const RVector<double> Rx;
  const RVector<double> Ry;
  const RVector<int> RType;
  const RVector<double> RWeight;

  // destination matrix
  RVector<double> RNbd;

  // constructor
  CountNbdWrkr(const NumericVector r2,
               const NumericVector x, const NumericVector y,
               const IntegerVector Type, const NumericVector Weight,
               NumericVector Nbd)
    : r2(r2), Rx(x), Ry(y), RType(Type), RWeight(Weight), RNbd(Nbd) {}

  // count neighbors
  void operator()(std::size_t begin, std::size_t end) {
    double Distance2, dx, dy;
    double Nr = r2.length();
    double Npoints = RType.length();
    unsigned int k, c;
    // c is the index of case points in the RNbd output matrix, whilst i is their index in input data
    c = begin;

    for (unsigned int i = begin; i < end; i++) {
      // Point j is a neighbor of i.
      for (unsigned int j=0; j < Npoints; j++) {
        // if (i != j) {  // Keep the point itself
        // Calculate squared distance
        dx = Rx[i]-Rx[j];
        dy = Ry[i]-Ry[j];
        Distance2 = dx*dx + dy*dy;
        // Ignore point j if it is too far from point i
        if (Distance2 <= r2[Nr-1]) {
          // Find the column of the matrix corresponding to the distance
          k = 0;
          while (Distance2 > r2[k]) {
            k++;
          }
          // Add j's weight to i's neighborhood
          RNbd[c + k*Npoints + (RType[j]-1)*Npoints*Nr] += RWeight[j];
        }
        // }
      }
      c++;
    }
  }
};

//' parallelCountNbd
//'
//' Create a 3-D array containing the number of neighbors around each point, per species
//'
//' @param r The vector of distances to take into account.
//' @param NbSpecies The number of species of the community.
//' @param x,y The coordinates of the points.
//' @param Type A vector containing the species of each point (as integer, i.e. the factor code).
//' @param Weight A vector containing the weights of points.
//' @returns a numeric vector with the distances to convert to a matrix
//' @noRd
//' @examples
//' # Build a point pattern
//' X <- rspcommunity(1, size = 5, species_number = 3)
//' if(require(spatstat.geom)) {
//'   # Compute the distances
//'   matrix(parallelCountNbd(r = 0:1, NbSpecies = 3, x= X$x, y = X$y, 
//'     Type = marks(X)$PointType, Weight = marks(X)$PointType), nrow = 5)
//' }
// [[Rcpp::export]]
NumericVector parallelCountNbd(NumericVector r, IntegerVector NbSpecies,
                               NumericVector x, NumericVector y,
                               IntegerVector Type, NumericVector Weight) {

  // allocate the output vector
  NumericVector Nbd(Type.length() * r.length() * NbSpecies[0]);

  // CountNbd functor
  CountNbdWrkr countNbdWrkr(r*r, x, y, Type, Weight, Nbd);

  // call parallelFor to do the work
  parallelFor(0, Type.length(), countNbdWrkr);

  // return the output matrix
  return Nbd;
}
