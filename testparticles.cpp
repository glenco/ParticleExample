 /*
  */

#include <slsimlib.h>
#include <sstream>
#include <iomanip>
#include <omp.h>
#include <thread>
#include <mutex>

#include "particle_halo.h"
#include "particle_halo2.h"
#include "point.h"
#include "gridmap.h"
#include "oTreeNB.h"

using namespace std;

int main(int arg,char **argv){
  
  COSMOLOGY cosmo(CosmoParamSet::Planck1yr);

  std::string out_dir = "output4/";

  Utilities::LOGPARAMS log_params(out_dir+"params");
  
  Point_2d rotation_vector(0,0);
  
  PosType zl=0.4;                   // redshift of lens
  log_params("zl",zl);
  PosType z_source = 2.0;           // redshift of source
  log_params("z_source",z_source);
  int Nsmooth = 5; // number of neighbors for smoothing scale
  log_params("Nsmooth",Nsmooth);
  
  long seed = -28976391;
  log_params("seed",seed);
  /**********************************************************/
  
  Utilities::RandomNumbers_NR random(seed);   //*** random number generator
 
  if(!Utilities::IO::check_directory(out_dir)){
    std::cout << "Creating directory " << out_dir << std::endl;
    Utilities::IO::make_directories(out_dir);
  }

  Lens lens(&seed,z_source,cosmo);

  Point_2d center(0,0);
  
  // set up a LensHalo from a text file with the particle data
  LensHaloParticles<float> halo("particles.dm.txt"
                        ,SimFileFormat::ascii
                        ,zl
                        ,0. // inverse area
                        ,3*1.807953375000000000e6 // particle mass
                        ,cosmo 
                        ,Nsmooth
                        ,8  /// number of buckets in tree
                        ,0.1      /// opening angle for tree
                        ,true /// re-center on center of mass
                        );/**/
    
                                                 
    
    // insert halos into lens
     // this is moved instead of inserted to avoid a copy
    lens.moveinMainHalo(halo, true);
  
  
  // here you can rotate each simulation about the center of mass
  //Point_2d theta(PI/2,0);
  //Point_3d<float> xo = lens.getMainHalo<LensHaloParticles<float> >(i)->CenterOfMass();
  //for(int i=0 ; i < lens.getNMainHalos<LensHaloParticles<float> >()  ; ++i){
  //  lens.getMainHalo<LensHaloParticles<float> >(i)->rotate(theta[0],theta[1],xo);
  //}
  
  center /= cosmo.angDist(zl);
 
  std::cout << "making gridmap ... ";
  double range = 30.0 * arcsecTOradians; // range of grids in radians
  GridMap gridmap(&lens,512,center.x,range/2);
  std::cout << "done." << std::endl;

  // set the redshift of the source plane
  lens.ResetSourcePlane(z_source,false);
  
  // output some maps
  gridmap.writeFits<float>(LensingVariable::KAPPA,out_dir+"particles_kappa.fits");
  gridmap.writeFits<float>(LensingVariable::INVMAG,out_dir+"particles_invmag.fits");
  gridmap.writeFits<float>(LensingVariable::ALPHA1,out_dir+"particles_a1.fits");
  gridmap.writeFits<float>(LensingVariable::ALPHA2,out_dir+"particles_a2.fits");
  gridmap.writeFits<float>(LensingVariable::GAMMA1,out_dir+"particles_g1.fits");
  gridmap.writeFits<float>(LensingVariable::GAMMA2,out_dir+"particles_g2.fits");
  
  // find the critical curves
  std::vector<ImageFinding::CriticalCurve> crit_curve;
  ImageFinding::find_crit(lens, gridmap,crit_curve);
  
  //*** plot caustic curves
  if(crit_curve.size() > 0){
    ImageFinding::CriticalCurve::print(out_dir+"particles_crit_curves.csv",crit_curve);
  }
  
  std::cout << "Number of caustics : "<< crit_curve.size() << std::endl;
  
  if(crit_curve.size() == 0){
    cout << "Exiting" << endl;
    exit(1);
  }
  //*** print information about the critical curves that were found
  PosType rmax,rmin,rave;
  if(crit_curve.size() > 0){
    std::string type;
    for(int i=0;i<crit_curve.size();++i){
      type = to_string( crit_curve[i].type );
      std::cout << "  " << i << " type " << to_string(crit_curve[i].type) << std::endl;
      crit_curve[i].CausticRadius(rmax,rmin,rave);
      std::cout << "      caustic " << crit_curve[i].caustic_center << " | " << crit_curve[i].caustic_area << " " << rmax << " " << rmin << " " << rave << std::endl;
      crit_curve[i].CriticalRadius(rmax,rmin,rave);
      std::cout << "      critical " << crit_curve[i].critical_center << " | " << crit_curve[i].critical_area << " " << rmax << " " << rmin << " " << rave << std::endl;
    }
  }
  
  //**** put a source in and map its images
  //****************************************
  
  //*** find a source position within the tangential caustic
  std::vector<Point_2d> y;                    // vector for source positions
  crit_curve[0].RandomSourcesWithinCaustic(1,y,random); // get random points within first caustic

  PosType zs = 2; //** redshift of source
  //** make a Sersic source, there are a number of other ones that could be used
  SourceSersic source(23,0.02,0,1,0.5,zs,23,Band::EUC_VIS);
  
  source.setTheta(y[0]);
  
  /** reset the source plane in the lens from the one given in the
   parameter file to this source's redshift
   */
  lens.ResetSourcePlane(zs,false);
  
  std::vector<ImageInfo> imageinfo;
  int Nimages;

  std::cout << "Mapping source ..." << std::endl;
  
  // add a source to the source plane
  gridmap.RefreshSurfaceBrightnesses(&source);

  // make a PixelMap with the image of the lensed source in it
  PixelMap<float> map = gridmap.getPixelMapFlux<float>();
  
  // The following produces an image of the lensed source using the rays that
  // were created within the grid and refined when finding the caustics.
  
  // You can add a source directly to the PixelMap<> without lensing with
  map.AddSource(source);
  
  // make a fits image 
  map.printFITS(out_dir+"!image_unrefined.fits");
  
  return 0;
}

