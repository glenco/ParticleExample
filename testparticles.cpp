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
  Point_3d<float> xo(0,0,0); // center of excised region in comoving Mpc/h
 /** auto lenshalos = LensHaloParticles<float>::MakeLensHaloParticle(
    "../ClusterCaustics/DataFiles/snap_058"
    ,SimFileFormat::gadget2
    ,0.    /// inv_area
    ,cosmo 
    ,5.    /// Nsmooth
    ,8     /// number of buckets in tree
    ,0.1   /// opening angle for tree
    ,0     /// rmax for excision
    ,xo    /// center of excised region in comoving Mpc/h
    ,true  /// verbose 
  );
  */
  { // tests
    Utilities::RandomNumbers_NR ran(-28976391);
    unsigned long nparticles = 100;
    PosType boundary_p1[3] = {-1000,-1000,-1000};
    PosType boundary_p2[3] = {1000,1000,1000};
    //std::vector<double> xp(3*nparticles);
    std::vector<Point_3d<double> > xp(nparticles);
    for(auto &p : xp){
      p[0] = ran()*2000 -1000;
      p[1] = ran()*1000 -1000;
      p[2] = ran()*200;
    }

    OTreeNB<Point_3d<double> > tree(xp.data(),xp.size());
    tree.build(10);

    Point_2d ray(0,0);
    double alpha[2],kappa,gamma[3],phi;
    alpha[0] = alpha[1] = 0.0;
    gamma[0] = gamma[1] = gamma[2] = 0.0;
    kappa = phi = 0.0;
    tree.force2D(  
      ray.x // ray direction
      ,1.0e9  // particle mass
      ,0.  // smooth_factor
      ,0.1 // theta2
      ,0. // inv_area
      ,alpha   // not zeroed
      ,&kappa // not zeroed
      ,gamma // not zeroed
      ,&phi   // not zeroed
    );

    std::cout << " Total Number of branches " << tree.getTotalBranches() << std::endl;
    std::cout << " Used branches " << tree.size() << std::endl;
    std::cout << " Depth " << tree.getDepth() << std::endl;
  
  }/**/

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

  Point_3d<double> center;
  /*{
   
    double range = 15.0 * arcsecTOradians; // range of grids in radians
    std::vector<RAY> rays(5);
    for(RAY& r : rays){
      r.x[0] = range * (random() - 0.5);
      r.x[1] = range * (random() - 0.5);
      r.z = z_source; // set the source redshift
    }
    
    Utilities::LOGDATA file("theta_test.csv");
    time_t t = time(nullptr);

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
                        );

    Lens lens(&seed,z_source,cosmo);
    lens.moveinMainHalo(halo, true);
    std:cout << time(nullptr) - t << " seconds to construct LensHaloParticles" << std::endl;
    for(float theta = 0 ; theta <= 1 ; theta += 0.01){

      Utilities::LOGDATA::LINE line;
    
      //halo.resetForceAngle(theta);
      line["theta"] = theta;

      // access the main halo and reset the force angle
      lens.getMainHalo<LensHaloParticles<float> >(0)->resetForceAngle(theta);
      
      // this is moved instead of inserted to avoid a copy
      //center /= cosmo.angDist(zl);
 
      // set the redshift of the source plane
      lens.ResetSourcePlane(z_source,false);
      lens.rayshooterInternal(rays.size(),rays.data());
      
      int i=0;
      for(auto &r : rays){
        line["alpha_x"+std::to_string(i)] = r.alpha()[0];
        line["alpha_y"+std::to_string(i)] = r.alpha()[1];
        line["kappa"+std::to_string(i)] = r.kappa();
        line["gamma1"+std::to_string(i)] = r.gamma1();
        line["gamma2"+std::to_string(i)] = r.gamma2();
        line["dt"+std::to_string(i)] = r.dt;
        ++i;
      }
      file.add(line);
    }
  }
exit(0);*/
  
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
    
                                                 
    //center = halo.CenterOfMass();
    // insert halos into lens
     // this is moved instead of inserted to avoid a copy
    lens.moveinMainHalo(halo, true);
    
    // Another way of creating particle halos that is more flexible with
    // more file types is commented out here.
//    MakeParticleLenses halomaker("particles.dm.txt",ascii,Nsmooth,true);
//    halomaker.CreateHalos(cosmo,zl);
//    //  There is a seporate halo for each type of particle
//    for(auto h : halomaker.halos){
//       lens.insertMainHalo(h,zl, true);
//    }
//
  

  // here you can rotate each simulation independently
  //for(int i=0 ; i < lens.getNMainHalos<LensHP>()  ; ++i){
  //  lens.getMainHalo<LensHP>(i)->rotate(theta);
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
  std::vector<Point_2d> y;                    //*** vector for source positions
  crit_curve[0].RandomSourcesWithinCaustic(1,y,random); //*** get random points within first caustic

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

