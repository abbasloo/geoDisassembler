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
#include <pcl/io/ply_io.h>
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkSphereSource.h>
#include <vtkConeSource.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkCleanPolyData.h>
#include <vtkAppendPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkVersion.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkOBJReader.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkSTLWriter.h>



using namespace Eigen;

int
main (int argc, char** argv)
{
  vtkSmartPointer<vtkPolyData> input = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);
  std::ifstream file("pclist");
  std::ifstream filee("meshlist");
  int i = 0;
  if (file.is_open()) {
      std::string line;
      std::string linee;
      while (getline(file, line) && getline(filee, linee)) {
      		pcl::PointCloud<pcl::PointXYZ>::Ptr cloudHere (new pcl::PointCloud<pcl::PointXYZ>);
                pcl::io::loadPCDFile <pcl::PointXYZ> (line, *cloudHere);

                vtkSmartPointer<vtkPolyData> inputHere = vtkSmartPointer<vtkPolyData>::New();
                vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
                reader->SetFileName(linee.c_str());
                std::cout<<linee.c_str()<<std::endl;
                reader->Update();
                inputHere->ShallowCopy(reader->GetOutput());
                if (i == 0) {cloud = cloudHere; input = inputHere;}
                if (i > 0) 
                {
                    *cloud += *cloudHere;

	            #if VTK_MAJOR_VERSION <= 5
			  appendFilter->AddInputConnection(input->GetProducerPort());
			  appendFilter->AddInputConnection(inputHere->GetProducerPort());
	            #else
			  appendFilter->AddInputData(input);
			  appendFilter->AddInputData(inputHere);
	            #endif
		    appendFilter->Update();
                    input = appendFilter->GetOutput();
                }
                i += 1;
      }
      file.close();
      filee.close();
   }

  // Write the file
  vtkSmartPointer<vtkSTLWriter> writer_ = vtkSmartPointer<vtkSTLWriter>::New();
  writer_->SetFileName(argv[1]);
  #if VTK_MAJOR_VERSION <= 5
    writer_->SetInput(appendFilter->GetOutput());
  #else
    writer_->SetInputData(appendFilter->GetOutput());
  #endif
  writer_->Write();

  // Write the file
  vtkSmartPointer<vtkPolyDataWriter> writer__ = vtkSmartPointer<vtkPolyDataWriter>::New();
  writer__->SetFileName(argv[2]);
  #if VTK_MAJOR_VERSION <= 5
    writer__->SetInput(appendFilter->GetOutput());
  #else
    writer__->SetInputData(appendFilter->GetOutput());
  #endif
  writer__->Write();

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
  gp3Mesh.setSearchRadius (1.0);

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
  pcl::PLYWriter writer;
  writer.write (argv[5], *cloud, false);
}
