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
#include <pcl/features/organized_edge_detection.h>
#include <pcl/console/time.h>
#include <pcl/features/principal_curvatures.h>



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
  /*
  // Smoothing
  pcl::PointCloud<pcl::PointNormal> mlsHere_points;
  pcl::MovingLeastSquares<pcl::PointXYZ, pcl::PointNormal> mlsHere;
  pcl::search::KdTree<pcl::PointXYZ>::Ptr treeHere (new pcl::search::KdTree<pcl::PointXYZ>);
  mlsHere.setComputeNormals (true);
  mlsHere.setInputCloud (cloud);
  mlsHere.setPolynomialFit (true);
  mlsHere.setSearchMethod (treeHere);
  mlsHere.setSearchRadius (5*std::atof(argv[2]));
  mlsHere.process (mlsHere_points);
  pcl::copyPointCloud(*cloud, mlsHere_points);
 
  pcl::PointCloud<pcl::PointXYZRGBA>::Ptr cloudrgba (new pcl::PointCloud<pcl::PointXYZRGBA>);
  for (size_t i = 0; i < cloud->points.size (); ++i)
  {
      pcl::PointXYZRGBA p;
      p.r = 255;
      p.g = 255;
      p.b = 255;
      p.a = 255;
      p.x = cloud->points[i].x;
      p.y = cloud->points[i].y;
      p.z = cloud->points[i].z;
      cloudrgba->push_back(p);
  }
  pcl::PCDWriter writerr;
  writerr.write ("check.pcd", *cloudrgba, false);

  pcl::search::Search<pcl::PointXYZRGBA>::Ptr tree = boost::shared_ptr<pcl::search::Search<pcl::PointXYZRGBA> > (new pcl::search::KdTree<pcl::PointXYZRGBA>);
  pcl::PointCloud <pcl::Normal>::Ptr normals (new pcl::PointCloud <pcl::Normal>);
  pcl::NormalEstimation<pcl::PointXYZRGBA, pcl::Normal> n;
  n.setSearchMethod (tree);
  n.setInputCloud (cloudrgba);
  n.setKSearch (100);
  n.compute (*normals);

  //pcl::OrganizedEdgeFromNormals<pcl::PointXYZRGBA, pcl::Normal, pcl::Label> oed;
  pcl::OrganizedEdgeFromRGBNormals<pcl::PointXYZRGBA, pcl::Normal, pcl::Label> oed;
  oed.setInputNormals (normals);
  oed.setInputCloud (cloudrgba);
  oed.setDepthDisconThreshold (0.02);
  oed.setMaxSearchNeighbors (100);
  pcl::PointCloud<pcl::Label> labels;
  std::vector<pcl::PointIndices> label_indices;
  oed.compute (labels, label_indices);

  pcl::PointCloud<pcl::PointXYZRGBA>::Ptr occluding_edges (new pcl::PointCloud<pcl::PointXYZRGBA>),
                                          occluded_edges (new pcl::PointCloud<pcl::PointXYZRGBA>),
                                          boundary_edges (new pcl::PointCloud<pcl::PointXYZRGBA>),
                                          high_curvature_edges (new pcl::PointCloud<pcl::PointXYZRGBA>),
                                          rgb_edges (new pcl::PointCloud<pcl::PointXYZRGBA>);

  pcl::copyPointCloud (*cloudrgba, label_indices[0].indices, *boundary_edges);
  std::cout<< label_indices[0] <<std::endl; 
  //pcl::copyPointCloud (*cloudrgba, label_indices[1].indices, *occluding_edges);
  //pcl::copyPointCloud (*cloudrgba, label_indices[2].indices, *occluded_edges);
  //pcl::copyPointCloud (*cloudrgba, label_indices[3].indices, *high_curvature_edges);
  //pcl::copyPointCloud (*cloudrgba, label_indices[4].indices, *rgb_edges);
  */
   
   if (std::atof(argv[6]) == 1) {
          int scale = (int) std::atof(argv[7]);
	  //pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
	  //kdtree.setInputCloud (cloud);          
	  
          // Downsample the Shell using a voxel grid class
	  pcl::VoxelGrid<pcl::PointXYZ> vg;
	  vg.setInputCloud (cloud);
	  vg.setLeafSize (std::atof(argv[2]), std::atof(argv[2]), std::atof(argv[2]));
	  vg.setDownsampleAllData (true);
	  vg.filter (*cloud);
	  pcl::PointXYZ p_m;
	  for (size_t i = 0; i < cloud->points.size (); ++i)
	  {
	      p_m.x = cloud->points[i].x/cloud->points.size ();
	      p_m.y = cloud->points[i].y/cloud->points.size ();
	      p_m.z = cloud->points[i].z/cloud->points.size ();
          }
	  // Normal estimation
	  pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> nMesh;
	  pcl::PointCloud<pcl::Normal>::Ptr normalsMesh (new pcl::PointCloud<pcl::Normal>);
	  pcl::search::KdTree<pcl::PointXYZ>::Ptr treeMesh (new pcl::search::KdTree<pcl::PointXYZ>);
	  treeMesh->setInputCloud (cloud);
	  nMesh.setInputCloud (cloud);
	  nMesh.setSearchMethod (treeMesh);
	  nMesh.setKSearch (100);
	  nMesh.compute (*normalsMesh);

	  // Setup the principal curvatures computation
	  pcl::PrincipalCurvaturesEstimation<pcl::PointXYZ, pcl::Normal, pcl::PrincipalCurvatures> principal_curvatures_estimation;
	  principal_curvatures_estimation.setInputCloud (cloud);
	  principal_curvatures_estimation.setInputNormals (normalsMesh);
	  principal_curvatures_estimation.setSearchMethod (treeMesh);
	  principal_curvatures_estimation.setRadiusSearch (1.0);
	  pcl::PointCloud<pcl::PrincipalCurvatures>::Ptr principal_curvatures (new pcl::PointCloud<pcl::PrincipalCurvatures> ());
	  principal_curvatures_estimation.compute (*principal_curvatures);
          //std::cout << principal_curvatures->points.size() << ',' << cloud->points.size() << std::endl;

	  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_g (new pcl::PointCloud<pcl::PointXYZ>);
	  for (size_t i = 0; i < cloud->points.size (); ++i)
	  {
	      pcl::PointXYZ p;
	      p.x = cloud->points[i].x;
	      p.y = cloud->points[i].y;
	      p.z = cloud->points[i].z;
	      cloud_g->push_back(p);
	      int t = 0;
              if (true){
		      while (t <= 20)
		      {
                          /*
			  float t1 = (float) (rand() % 200 + 0)/100.0; t1 -= 1.0;
			  float t2 = (float) (rand() % 200 + 0)/100.0; t2 -= 1.0;
			  float t3 = (float) (rand() % 200 + 0)/100.0; t3 -= 1.0;
			  float T = (float) (rand() % scale + 1)/100.0; T *= std::atof(argv[2]);
		          float a = std::sqrt(t1*t1 + t2*t2 + t3*t3);
		          t1 /= a; t2 /= a; t3 /= a;
			  float val = std::abs(normalsMesh->points[i].normal_x*T*t1 + normalsMesh->points[i].normal_y*T*t2 + normalsMesh->points[i].normal_z*T*t3); 
			  if (val <= 1e-3 && a != 0.0 && T != 0.0) {
		                  //std::cout << t1 << ',' << t2 << ',' << t3 << ',' << T << std::endl;
		                  //std::cout << normalsMesh->points[i].normal_x << ',' << normalsMesh->points[i].normal_y << ',' << normalsMesh->points[i].normal_z << std::endl;
				  pcl::PointXYZ pp;  
				  pp.x = cloud->points[i].x + T*t1;
				  pp.y = cloud->points[i].y + T*t2;
				  pp.z = cloud->points[i].z + T*t3; 
		                  bool doIt = true;
				  //int K = 2;
				  //std::vector<int> pointIdxNKNSearch(K);
				  //std::vector<float> pointNKNSquaredDistance(K);
				  //if ( kdtree.nearestKSearch (pp, K, pointIdxNKNSearch, pointNKNSquaredDistance) > 0)
				  //{
				  // 	for (size_t p = 0; p < pointIdxNKNSearch.size (); ++p)
				  //	{
				  //		if (pointNKNSquaredDistance[p] < 0.1*std::atof(argv[2]))
				  //		{
				  //		  doIt = false;
				  //	        } 
				  //	}
				  //}
				  if (doIt) {cloud_g->push_back(pp); t += 1;}
			  }
			  */
			  float t1 = principal_curvatures->points[i].principal_curvature_x;
			  float t2 = principal_curvatures->points[i].principal_curvature_y;
			  float t3 = principal_curvatures->points[i].principal_curvature_z;
                          float val_c = t1*(cloud->points[i].x-p_m.x) + t2*(cloud->points[i].y-p_m.y) + t3*(cloud->points[i].z-p_m.z);
			  if (val_c < 0.0) {
                                t1 *= -1; t2 *= -1; t3 *= -1;
			  }
		          float a = std::sqrt(t1*t1 + t2*t2 + t3*t3);
		          t1 /= a; t2 /= a; t3 /= a;
                          float T = (float) (rand() % scale + 1)/100.0; T *= std::atof(argv[2]);
			  pcl::PointXYZ pp;  
			  pp.x = cloud->points[i].x + T*t1;
			  pp.y = cloud->points[i].y + T*t2;
			  pp.z = cloud->points[i].z + T*t3; 
                          cloud_g->push_back(pp);
                          t += 1;
		      }
               }
	  }
	  cloud = cloud_g;
  }

  // Downsample the cloud using a voxel grid class
  pcl::VoxelGrid<pcl::PointXYZ> vg;
  vg.setInputCloud (cloud);
  vg.setLeafSize (std::atof(argv[2]), std::atof(argv[2]), std::atof(argv[2]));
  vg.setDownsampleAllData (true);
  vg.filter (*cloud);

  // Normal estimation
  pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> nMesh;
  pcl::PointCloud<pcl::Normal>::Ptr normalsMesh (new pcl::PointCloud<pcl::Normal>);
  pcl::search::KdTree<pcl::PointXYZ>::Ptr treeMesh (new pcl::search::KdTree<pcl::PointXYZ>);
  treeMesh->setInputCloud (cloud);
  nMesh.setInputCloud (cloud);
  nMesh.setSearchMethod (treeMesh);
  nMesh.setKSearch (100);
  nMesh.compute (*normalsMesh);

  // Concatenate the XYZ and normal fields
  pcl::PointCloud<pcl::PointNormal>::Ptr cloud_with_normalsMesh (new pcl::PointCloud<pcl::PointNormal>);
  pcl::concatenateFields (*cloud, *normalsMesh, *cloud_with_normalsMesh);

  // Create search tree
  pcl::search::KdTree<pcl::PointNormal>::Ptr treeMeshMesh (new pcl::search::KdTree<pcl::PointNormal>);
  treeMeshMesh->setInputCloud (cloud_with_normalsMesh);

  // Initialize objects
  pcl::GreedyProjectionTriangulation<pcl::PointNormal> gp3Mesh;
  pcl::PolygonMesh triangles;

  // Set the maximum distance between connected points (maximum edge length)
  gp3Mesh.setSearchRadius (10*std::atof(argv[2]));

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
  pcl::io::saveVTKFile (argv[3], triangles);
  pcl::io::saveOBJFile (argv[4], triangles);
  pcl::PCDWriter writer;
  writer.write (argv[5], *cloud, false);
}
