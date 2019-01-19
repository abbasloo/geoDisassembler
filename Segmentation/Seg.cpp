#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/search/search.h>
#include <pcl/search/kdtree.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/features/normal_3d.h>
#include <pcl/filters/passthrough.h>
#include <pcl/segmentation/region_growing.h>
#include <pcl/segmentation/region_growing_rgb.h>
#include <pcl/io/vtk_io.h>
#include <pcl/surface/gp3.h>
#include <pcl/io/obj_io.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/features/vfh.h>
#include <Eigen/Dense>
#include <pcl/filters/extract_indices.h>
#include <pcl/filters/statistical_outlier_removal.h>
#include <pcl/common/common.h>
#include <pcl/surface/mls.h>



using namespace Eigen;

int
main (int argc, char** argv)
{
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);
  if ( pcl::io::loadPCDFile <pcl::PointXYZ> (argv[1], *cloud) == -1)
  {
    std::cout << "Cloud reading failed." << std::endl;
    return (-1);
  }

  //------------------------------------------------SEGMENTATION----------------------------------------//  

  pcl::PointXYZ minPt, maxPt;
  pcl::getMinMax3D (*cloud, minPt, maxPt);

  // Create exterior shell
  pcl::PointCloud<pcl::PointXYZ>::Ptr Shell (new pcl::PointCloud<pcl::PointXYZ>);

  // Filling the boundaries
  int l = 100;
  float dx = (maxPt.x-minPt.x)/l;
  float dy = (maxPt.y-minPt.y)/l;
  float dz = (maxPt.z-minPt.z)/l;
  pcl::PointXYZ p;
  for (size_t i = 0; i <= l; i++) 
  {
	for (size_t j = 0; j <= l; j++) 
	{
		for (size_t k = 0; k < 1; k++) 
		{
		   p.x = minPt.x+i*dx;
		   p.y = maxPt.y;
		   p.z = minPt.z+j*dz;
		   Shell->push_back(p);
		   p.x = minPt.x+i*dx;
		   p.y = minPt.y+j*dy;
		   p.z = minPt.z;
		   Shell->push_back(p);
		   p.x = minPt.x+i*dx;
		   p.y = minPt.y+j*dy;
		   p.z = maxPt.z;
		   Shell->push_back(p);
		   p.x = minPt.x;
		   p.y = minPt.y+i*dy;
		   p.z = minPt.z+j*dz;
		   Shell->push_back(p);
		   p.x = maxPt.x;
		   p.y = minPt.y+i*dy;
		   p.z = minPt.z+j*dz;
		   Shell->push_back(p);
		}
	}
   }
  pcl::PCDWriter writer;
  writer.write ("../data/pc/Shell.pcd", *Shell, false);

  // Downsample the Shell using a voxel grid class
  pcl::VoxelGrid<pcl::PointXYZ> vg;
  vg.setInputCloud (Shell);
  vg.setLeafSize (10000.0, 10000.0, 10000.0);
  vg.setDownsampleAllData (true);
  vg.filter (*Shell);

  // Normal estimation
  pcl::search::Search<pcl::PointXYZ>::Ptr tree = boost::shared_ptr<pcl::search::Search<pcl::PointXYZ> > (new pcl::search::KdTree<pcl::PointXYZ>);
  pcl::PointCloud <pcl::Normal>::Ptr normals (new pcl::PointCloud <pcl::Normal>);
  pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> n;
  n.setSearchMethod (tree);
  n.setInputCloud (cloud);
  n.setKSearch (100);
  n.compute (*normals);
  pcl::PointCloud<pcl::PointNormal>::Ptr cloud_with_normals (new pcl::PointCloud<pcl::PointNormal>);
  pcl::concatenateFields (*cloud, *normals, *cloud_with_normals);

  // Filtering
  pcl::IndicesPtr indices (new std::vector <int>);
  pcl::PassThrough<pcl::PointXYZ> pass;
  pass.setInputCloud (cloud);
  pass.setFilterFieldName ("y");
  pass.setFilterLimits (0.0, 1.0);
  pass.filter (*indices);

  // Region growing
  pcl::RegionGrowing<pcl::PointXYZ, pcl::Normal> reg;
  reg.setMinClusterSize (0.00005*cloud->points.size ());
  reg.setMaxClusterSize (cloud->points.size ());
  reg.setSearchMethod (tree);
  reg.setNumberOfNeighbours (100);
  reg.setInputCloud (cloud);
  //reg.setIndices (indices);
  reg.setInputNormals (normals);
  reg.setSmoothnessThreshold (std::atof(argv[2]) / 180.0 * M_PI);
  reg.setCurvatureThreshold (std::atof(argv[3]));
  std::vector <pcl::PointIndices> clusters;
  reg.extract (clusters);

  //------------------------------------------------FINDING NEIGHBOURS----------------------------------------//  

  pcl::PointCloud<pcl::PointXYZ>::Ptr cloudSeg (new pcl::PointCloud<pcl::PointXYZ>);
  int segSize = 0;
  for (size_t counter = 0; counter < clusters.size (); ++counter)
  {
	  segSize += clusters[counter].indices.size ();
  }
  cloudSeg->width = segSize;
  cloudSeg->height = 1;
  cloudSeg->points.resize (cloudSeg->width * cloudSeg->height);

  // Files
  std::ofstream myfile;
  myfile.open ("SegCon.csv");
  std::ofstream myfilee;
  myfilee.open ("Seg.csv");
  std::ofstream myfileee;
  myfileee.open ("Feat.csv");
  std::ofstream myfileeee;
  myfileeee.open ("Ext.csv");

  // Saving the result
  std::cout << "Number of clusters is equal to " << clusters.size () << std::endl;
  pcl::PointCloud <pcl::PointXYZRGB>::Ptr colored_cloud = reg.getColoredCloud ();
  int cn = 0;
  for (size_t counter = 0; counter < clusters.size (); ++counter)
  {
          bool doIt = true;
	  //std::cout <<"cluster "<< counter <<" has " << clusters[counter].indices.size () << " points." << std::endl;
          std::vector<int> Neighbours;
	  pcl::PointCloud<pcl::PointXYZ>::Ptr cloudHere (new pcl::PointCloud<pcl::PointXYZ>);
	  cloudHere->width = clusters[counter].indices.size ();
	  cloudHere->height = 1;
	  cloudHere->points.resize (cloudHere->width * cloudHere->height);

	  // mean and covariance as features
	  MatrixXf covariance_matrix(3, 3);
	  MatrixXf mean_matrix(3, 1);
	  MatrixXf dummy(3, 1);
	  MatrixXf EV(3, 1);
	  float ysigma = 0.0;
	  for (size_t i = 0; i < 3; ++i)
	  {
	       mean_matrix(i, 0) = 0;
	       for (size_t j = 0; j < 3; ++j)
	       {
		   covariance_matrix(i, j) = 0;
               }
	  }
	  for (size_t i = 0; i < clusters[counter].indices.size (); ++i)
	  {
	      cloudHere->points[i].x = cloud->points[clusters[counter].indices[i]].x; 
	      cloudHere->points[i].y = cloud->points[clusters[counter].indices[i]].y; 
	      cloudHere->points[i].z = cloud->points[clusters[counter].indices[i]].z;
	      mean_matrix(0, 0) += (float) cloudHere->points[i].x/clusters[counter].indices.size (); 
	      mean_matrix(1, 0) += (float) cloudHere->points[i].y/clusters[counter].indices.size ();
	      mean_matrix(2, 0) += (float) cloudHere->points[i].z/clusters[counter].indices.size ();
	      cloudSeg->points[cn].x = cloud->points[clusters[counter].indices[i]].x; 
	      cloudSeg->points[cn].y = cloud->points[clusters[counter].indices[i]].y; 
	      cloudSeg->points[cn].z = cloud->points[clusters[counter].indices[i]].z;
	      cn++; 
	  }

	  for (size_t i = 0; i < clusters[counter].indices.size (); ++i)
	  {
	      dummy(0, 0) = (float) cloudHere->points[i].x; 
	      dummy(1, 0) = (float) cloudHere->points[i].y;
	      dummy(2, 0) = (float) cloudHere->points[i].z;	
	      ysigma += (float) std::pow(cloudHere->points[i].y - mean_matrix(1, 0), 2)/clusters[counter].indices.size ();
	      covariance_matrix += ((dummy-mean_matrix)*(dummy-mean_matrix).transpose())/clusters[counter].indices.size ();
	  }
	  SelfAdjointEigenSolver<MatrixXf> eigensolver(covariance_matrix);
	  //EigenSolver<MatrixXf> eigensolver(covariance_matrix, false);
	  //EV = eigensolver.eigenvalues().real();
	  EV = eigensolver.eigenvalues();

	  // Downsample the cloud using a voxel grid class
	  pcl::PointCloud<pcl::PointXYZ>::Ptr cloudHereOut (new pcl::PointCloud<pcl::PointXYZ>);
	  pcl::VoxelGrid<pcl::PointXYZ> vg;
	  vg.setInputCloud (cloudHere);
	  vg.setLeafSize (10000.0, 10000.0, 10000.0);
	  vg.setDownsampleAllData (true);
	  vg.filter (*cloudHereOut);
	  //cloudHereOut = cloudHere;

	  // Finding the neighbouring segments
	  for (size_t i = 0; i < cloudHereOut->points.size (); ++i)
	  {	
	      pcl::PointXYZ searchPoint;      
	      searchPoint.x = cloudHereOut->points[i].x;
	      searchPoint.y = cloudHereOut->points[i].y;
	      searchPoint.z = cloudHereOut->points[i].z; 
	      for (size_t counterr = 0; counterr < clusters.size (); ++counterr)
	      {
		    if ( counter != counterr )
                    {
			  pcl::PointCloud<pcl::PointXYZ>::Ptr cloudHereHere (new pcl::PointCloud<pcl::PointXYZ>);
	  		  cloudHereHere->width = clusters[counterr].indices.size ();
	  		  cloudHereHere->height = 1;
	  		  cloudHereHere->points.resize (cloudHereHere->width * cloudHereHere->height);
			  for (size_t ii = 0; ii < clusters[counterr].indices.size (); ++ii)
			  {
			      cloudHereHere->points[ii].x = cloud->points[clusters[counterr].indices[ii]].x;
			      cloudHereHere->points[ii].y = cloud->points[clusters[counterr].indices[ii]].y;
			      cloudHereHere->points[ii].z = cloud->points[clusters[counterr].indices[ii]].z;
			  }
		          pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
		          kdtree.setInputCloud (cloudHereHere);
		          int K = 20;
		          std::vector<int> pointIdxNKNSearch(K);
		          std::vector<float> pointNKNSquaredDistance(K);
		          if ( kdtree.nearestKSearch (searchPoint, K, pointIdxNKNSearch, pointNKNSquaredDistance) > 0)
		          {
				for (size_t p = 0; p < pointIdxNKNSearch.size (); ++p)
				{
					if (pointNKNSquaredDistance[p] <= std::atof(argv[4]))
					{
					      Neighbours.push_back (counterr);
					} 
                                }
		          }
	            }

              }
	  }

          // Finding exterior
	  pcl::PointXYZ searchPointShell; 
	  searchPointShell.x = mean_matrix(0, 0);
	  searchPointShell.y = mean_matrix(1, 0);
          searchPointShell.z = mean_matrix(2, 0);
          pcl::KdTreeFLANN<pcl::PointXYZ> kdtreeShell;
          kdtreeShell.setInputCloud (Shell);
          int KShell = 20;
          std::vector<int> pointIdxNKNSearchShell(KShell);
          std::vector<float> pointNKNSquaredDistanceShell(KShell);
          if ( kdtreeShell.nearestKSearch (searchPointShell, KShell, pointIdxNKNSearchShell, pointNKNSquaredDistanceShell) > 0)
          {
	        for (size_t p = 0; p < pointIdxNKNSearchShell.size (); ++p)
	        {
	        	if (pointNKNSquaredDistanceShell[p] <= std::atof(argv[5]) && doIt)
	        	{
		          doIt = false;
		          if (counter < clusters.size () - 1) myfileeee << counter << '\n'; 
		          if (counter == clusters.size () - 1) myfileeee << counter; 
		    } 
	        }
          }

	  //std::cout <<"cluster: " << counter << std::endl;
	  myfile << counter << ';';
	  myfilee << counter << ';' << mean_matrix(0, 0) << ';' << mean_matrix(1, 0) << ';' << mean_matrix(2, 0)<< ';';
	  myfilee << EV(0, 0) << ';' << EV(1, 0) << ';' << EV(2, 0) << ';' << ysigma << ';' << clusters[counter].indices.size () << '\n';
	  if (Neighbours.size () != 0)
	  {
   		std::sort(Neighbours.begin(), Neighbours.end()); 
    		auto last = std::unique(Neighbours.begin(), Neighbours.end());
    		Neighbours.erase(last, Neighbours.end());
		for (size_t i = 0; i < Neighbours.size (); ++i)
		{
			//std::cout << Neighbours[i]<<" "; 
			if (i < Neighbours.size () - 1) myfile << Neighbours[i] << ';'; 
			if (i == Neighbours.size () - 1) myfile << Neighbours[i]; 
		}
	  }
	  myfile << '\n';

          // Smoothing
	  pcl::PointCloud<pcl::PointNormal> mlsHere_points;
	  pcl::MovingLeastSquares<pcl::PointXYZ, pcl::PointNormal> mlsHere;
	  pcl::search::KdTree<pcl::PointXYZ>::Ptr treeHere (new pcl::search::KdTree<pcl::PointXYZ>);
	  mlsHere.setComputeNormals (true);
	  mlsHere.setInputCloud (cloudHere);
	  mlsHere.setPolynomialFit (true);
	  mlsHere.setSearchMethod (treeHere);
	  mlsHere.setSearchRadius (0.01);
	  mlsHere.process (mlsHere_points);

	  //std::cout << std::endl;
	  std::string filename = "../data/pc/voxels_seg_" + std::to_string(counter) +".pcd";
	  //writer.write (filename, *cloudHere, false);
	  pcl::io::savePCDFile (filename, mlsHere_points);

	  //-----------------------------------------------PARTS MESHING-------------------------------------------//
  
	  // Normal estimation
	  pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> nHere;
	  pcl::PointCloud<pcl::Normal>::Ptr normalsHere (new pcl::PointCloud<pcl::Normal>);
	  treeHere->setInputCloud (cloudHere);
	  nHere.setInputCloud (cloudHere);
	  nHere.setSearchMethod (treeHere);
	  nHere.setKSearch (100);
	  nHere.compute (*normalsHere);

          // Create the pfh estimator, and pass the cloud and normals to it
  	  pcl::VFHEstimation<pcl::PointXYZ, pcl::Normal, pcl::VFHSignature308> vfh;
  	  vfh.setInputCloud (cloudHere);
  	  vfh.setInputNormals (normalsHere);
 	  pcl::search::KdTree<pcl::PointXYZ>::Ptr treevfh (new pcl::search::KdTree<pcl::PointXYZ> ());
          vfh.setSearchMethod (treevfh);
          pcl::PointCloud<pcl::VFHSignature308>::Ptr vfhs (new pcl::PointCloud<pcl::VFHSignature308> ());
	  vfh.compute (*vfhs);
	  myfileee << counter << ';';
          for (size_t i = 0; i < vfhs->points.size (); ++i) 
	  {
	   	myfileee << vfhs->points[i];	
	  }	
	  myfileee <<'\n';

	  // Concatenate the XYZ and normal fields
	  pcl::PointCloud<pcl::PointNormal>::Ptr cloud_with_normalsHere (new pcl::PointCloud<pcl::PointNormal>);
	  pcl::concatenateFields (*cloudHere, *normalsHere, *cloud_with_normalsHere);

	  // Create search tree
	  pcl::search::KdTree<pcl::PointNormal>::Ptr treeHereHere (new pcl::search::KdTree<pcl::PointNormal>);
	  treeHereHere->setInputCloud (cloud_with_normalsHere);

	  // Initialize objects
	  pcl::GreedyProjectionTriangulation<pcl::PointNormal> gp3Here;
	  pcl::PolygonMesh trianglesHere;

	  // Set the maximum distance between connected points (maximum edge length)
	  gp3Here.setSearchRadius (5.0);

	  // Set typical values for the parameters
	  gp3Here.setMu (2.5);
	  gp3Here.setMaximumNearestNeighbors (100);
	  gp3Here.setMaximumSurfaceAngle(M_PI/4); 
	  gp3Here.setMinimumAngle(M_PI/18); 
	  gp3Here.setMaximumAngle(2*M_PI/3);
	  gp3Here.setNormalConsistency(true);

	  // Get result
	  gp3Here.setInputCloud (cloud_with_normalsHere);
	  gp3Here.setSearchMethod (treeHereHere);
	  gp3Here.reconstruct (trianglesHere);

	  // Additional vertex information
	  std::vector<int> partsHere = gp3Here.getPartIDs();
	  std::vector<int> statesHere = gp3Here.getPointStates();

	  // Saving the result
	  std::string filenamevtk = "../data/mesh/voxels_seg_" + std::to_string(counter) +".vtk";
	  std::string filenameobj = "../data/mesh/voxels_seg_" + std::to_string(counter) +".obj";
	  pcl::io::saveVTKFile (filenamevtk, trianglesHere);
	  pcl::io::saveOBJFile (filenameobj, trianglesHere);
  }
  writer.write ("../data/pc/voxels_seg.pcd", *cloudSeg, false);
  //writer.write ("../data/pc/voxels_seg.pcd", *colored_cloud, false);
  myfile.close();
  myfilee.close();
  myfileee.close();
  myfileeee.close();				

  //----------------------------------------------MESHING-----------------------------------------//  

  // Normal estimation
  pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> nMesh;
  pcl::PointCloud<pcl::Normal>::Ptr normalsMesh (new pcl::PointCloud<pcl::Normal>);
  pcl::search::KdTree<pcl::PointXYZ>::Ptr treeMesh (new pcl::search::KdTree<pcl::PointXYZ>);
  treeMesh->setInputCloud (cloudSeg);
  nMesh.setInputCloud (cloudSeg);
  nMesh.setSearchMethod (treeMesh);
  nMesh.setKSearch (100);
  nMesh.compute (*normalsMesh);

  // Concatenate the XYZ and normal fields
  pcl::PointCloud<pcl::PointNormal>::Ptr cloud_with_normalsMesh (new pcl::PointCloud<pcl::PointNormal>);
  pcl::concatenateFields (*cloudSeg, *normalsMesh, *cloud_with_normalsMesh);

  // Create search tree
  pcl::search::KdTree<pcl::PointNormal>::Ptr treeMeshMesh (new pcl::search::KdTree<pcl::PointNormal>);
  treeMeshMesh->setInputCloud (cloud_with_normalsMesh);

  // Initialize objects
  pcl::GreedyProjectionTriangulation<pcl::PointNormal> gp3Mesh;
  pcl::PolygonMesh triangles;

  // Set the maximum distance between connected points (maximum edge length)
  gp3Mesh.setSearchRadius (5.0);

  // Set typical values for the parameters
  gp3Mesh.setMu (2.5);
  gp3Mesh.setMaximumNearestNeighbors (100);
  gp3Mesh.setMaximumSurfaceAngle(M_PI/4); 
  gp3Mesh.setMinimumAngle(M_PI/18); 
  gp3Mesh.setMaximumAngle(2*M_PI/3);
  gp3Mesh.setNormalConsistency(true);

  // Get result
  gp3Mesh.setInputCloud (cloud_with_normalsMesh);
  gp3Mesh.setSearchMethod (treeMeshMesh);
  gp3Mesh.reconstruct (triangles);

  // Additional vertex information
  std::vector<int> parts = gp3Mesh.getPartIDs();
  std::vector<int> states = gp3Mesh.getPointStates();

  // Saving the result
  pcl::io::saveVTKFile ("../data/mesh/voxels_seg.vtk", triangles);
  pcl::io::saveOBJFile ("../data/mesh/voxels_seg.obj", triangles);
  
  return (0);
}
