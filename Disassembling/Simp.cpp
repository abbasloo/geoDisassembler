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

  // Downsample the Shell using a voxel grid class
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
  gp3Mesh.setSearchRadius (5*std::atof(argv[2]));

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
