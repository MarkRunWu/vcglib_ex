/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2012                                           \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/
#include <limits>
#include <stdio.h>
#include <vcg/complex/complex.h>
#include <vcg/math/perlin_noise.h>
#include <vcg/complex/algorithms/create/marching_cubes.h>
//#include <vcg/complex/algorithms/create/extended_marching_cubes.h>
#include <vcg/complex/algorithms/create/mc_trivial_walker.h>
#define M_PI_2     1.57079632679489661923
//#include <vcg/complex/algorithms/smooth.h>
#include <vcg/complex/algorithms/smooth.h>
#include <wrap/io_trimesh/export_ply.h>

// update
#include <vcg/complex/algorithms/update/topology.h>

// local optimization
#include <vcg/complex/algorithms/local_optimization.h>
#include <vcg/complex/algorithms/local_optimization/tri_edge_collapse_quadric.h>

using namespace std;
using namespace vcg;

typedef float ScalarType;


class MyFace;
class MyVertex;

struct MyUsedTypes : public UsedTypes<	Use<MyVertex>		::AsVertexType,
																				Use<MyFace>			::AsFaceType>{};

class MyVertex     : public Vertex< MyUsedTypes, vertex::Coord3f, vertex::BitFlags , vertex::Normal3f, vertex::VFAdj, vertex::Mark>{
      public:
        vcg::math::Quadric<double> &Qd() {return q;}
      private:
        math::Quadric<double> q;
};
class MyFace       : public Face< MyUsedTypes, face::VertexRef, face::BitFlags, face::Normal3f, face::VFAdj> {};

class MyMesh		: public vcg::tri::TriMesh< std::vector< MyVertex>, std::vector< MyFace > > {};

typedef BasicVertexPair<MyVertex> VertexPair;

class MyTriEdgeCollapse: public vcg::tri::TriEdgeCollapseQuadric< MyMesh, VertexPair, MyTriEdgeCollapse, QInfoStandard<MyVertex>  > {
            public:
            typedef  vcg::tri::TriEdgeCollapseQuadric< MyMesh,  VertexPair, MyTriEdgeCollapse, QInfoStandard<MyVertex>  > TECQ;
            typedef  MyMesh::VertexType::EdgeType EdgeType;
            inline MyTriEdgeCollapse(  const VertexPair &p, int i, BaseParameterClass *pp) :TECQ(p,i,pp){}
};

typedef SimpleVolume<SimpleVoxel> MyVolume;

int main(int /*argc*/ , char **/*argv*/)
{
	MyVolume	volume;
  
  typedef vcg::tri::TrivialWalker<MyMesh,MyVolume>	MyWalker;
	typedef vcg::tri::MarchingCubes<MyMesh, MyWalker>	MyMarchingCubes;
	MyWalker walker;
	

  // Simple initialization of the volume with some cool perlin noise

    volume.Init(Point3i(512,512,137));
    /*
  for(int i=0;i<512;i++)
    for(int j=0;j<512;j++)
      for(int k=0;k<137;k++)
        volume.Val(i,j,k)=(j-32)*(j-32)+(k-32)*(k-32)  + i*10*(float)math::Perlin::Noise(i*.2,j*.2,k*.2);
*/

    FILE* pfin = fopen("result.dat" , "rb");
    if( pfin == NULL )exit(1);
    char c;
    int w,h,d;
    w = 512;
    h = 512;
    d = 137;

    volume.Init(Point3i(w,h,d));
    char max_c = 0;
    for(int k=0;k<d;k++)
  for(int i=0;i<w;i++)
    for(int j=0;j<h;j++){
        fread( &c , 1 , 1 , pfin );
        if( c > max_c ) max_c = c;
        volume.Val(i,j,k)= c;
      }
  fclose( pfin );
  printf("max_c: %d\n" , (int)max_c );

	// MARCHING CUBES
	MyMesh		mc_mesh;
	printf("[MARCHING CUBES] Building mesh...");
	MyMarchingCubes					mc(mc_mesh, walker);
    walker.BuildMesh<MyMarchingCubes>(mc_mesh, volume, mc, 15);
    printf("generate face num %d, vertex num %d" , mc_mesh.vn , mc_mesh.fn );


    TriEdgeCollapseQuadricParameter qparams;
    qparams.QualityThr  =.3;
    float TargetError=std::numeric_limits<float>::max();
    int FinalSize=mc_mesh.fn*0.2;
    printf("target face num %d\n" , FinalSize);
    qparams.PreserveTopology = true;
    qparams.QuadricEpsilon=TargetError;
    qparams.NormalCheck=true;

    vcg::tri::UpdateBounding<MyMesh>::Box(mc_mesh);
    // decimator initialization
    vcg::LocalOptimization<MyMesh> DeciSession(mc_mesh,&qparams);

    int t1=clock();
    DeciSession.Init<MyTriEdgeCollapse>();
    int t2=clock();
    printf("Initial Heap Size %i\n",int(DeciSession.h.size()));

    DeciSession.SetTargetSimplices(FinalSize);
    DeciSession.SetTimeBudget(0.5f);

    if(TargetError< std::numeric_limits<float>::max() ) DeciSession.SetTargetMetric(TargetError);

    printf("start mesh decimator....\n");
    while(DeciSession.DoOptimization() && mc_mesh.fn>FinalSize && DeciSession.currMetric < TargetError)
      printf("Current Mesh size %7i heap sz %9i err %9g \r",mc_mesh.fn, int(DeciSession.h.size()),DeciSession.currMetric);

    int t3=clock();
    printf("mesh  %d %d Error %g \n",mc_mesh.vn,mc_mesh.fn,DeciSession.currMetric);
    printf("\nCompleted in (%i+%i) msec\n",t2-t1,t3-t2);

     tri::UpdateNormal<MyMesh>::PerVertexNormalizedPerFace(mc_mesh);

    printf("start taubin smooth....\n");
    //smooth
    tri::UpdateFlags<MyMesh>::FaceBorderFromNone(mc_mesh);
          int stepSmoothNum = 15;
          float lambda=0.5;
          float mu=-0.53;

    size_t cnt=tri::UpdateSelection<MyMesh>::VertexFromFaceStrict(mc_mesh);
    tri::Smooth<MyMesh>::VertexCoordTaubin(mc_mesh,stepSmoothNum,lambda,mu,cnt>0);
    //Log( "Smoothed %d vertices", cnt>0 ? cnt : mc_mesh.vn);
    tri::UpdateNormal<MyMesh>::PerVertexNormalizedPerFace(mc_mesh);

	vcg::tri::io::ExporterPLY<MyMesh>::Save( mc_mesh, "marching_cubes.ply");

	printf("OK!\n");
}
