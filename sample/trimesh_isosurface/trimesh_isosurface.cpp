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
//#include <wrap/io_trimesh/export_ply.h>
#include <GL/glut.h>
#include <wrap/io_trimesh/export_dae.h>

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

class MyVertex     : public Vertex< MyUsedTypes, vertex::Coord3f, vertex::BitFlags , vertex::Normal3f, vertex::VFAdj, vertex::Mark, vertex::Color4b>{
      public:
        vcg::math::Quadric<double> &Qd() {return q;}
      private:
        math::Quadric<double> q;
};
class MyFace       : public Face< MyUsedTypes, face::VertexRef, face::BitFlags, face::Normal3f, face::VFAdj , face::FFAdj> {};

class MyMesh		: public vcg::tri::TriMesh< std::vector< MyVertex>, std::vector< MyFace > > {};

namespace vcg
{
namespace tri
{

typedef BasicVertexPair<MyVertex> VertexPair;


typedef	SimpleTempData<MyMesh::VertContainer, math::Quadric<double> > QuadricTemp;


class QHelper
        {
        public:
      QHelper(){}
      static void Init(){}
      static math::Quadric<double> &Qd(MyVertex &v) {return TD()[v];}
      static math::Quadric<double> &Qd(MyVertex *v) {return TD()[*v];}
      static MyVertex::ScalarType W(MyVertex * /*v*/) {return 1.0;}
      static MyVertex::ScalarType W(MyVertex & /*v*/) {return 1.0;}
      static void Merge(MyVertex & /*v_dest*/, MyVertex const & /*v_del*/){}
      static QuadricTemp* &TDp() {static QuadricTemp *td; return td;}
      static QuadricTemp &TD() {return *TDp();}
        };


class MyTriEdgeCollapse: public vcg::tri::TriEdgeCollapseQuadric< MyMesh, VertexPair , MyTriEdgeCollapse, QHelper > {
                        public:
            typedef  vcg::tri::TriEdgeCollapseQuadric< MyMesh, VertexPair,  MyTriEdgeCollapse, QHelper> TECQ;
            inline MyTriEdgeCollapse(  const VertexPair &p, int i, BaseParameterClass *pp) :TECQ(p,i,pp){}
};
}
}
typedef SimpleVolume<SimpleVoxel> MyVolume;

int main(int argc , char **argv)
{
	MyVolume	volume;
  
  typedef vcg::tri::TrivialWalker<MyMesh,MyVolume>	MyWalker;
	typedef vcg::tri::MarchingCubes<MyMesh, MyWalker>	MyMarchingCubes;
	MyWalker walker;
	
    int w,h,d;
    float iso_value;
    const char* path;
    if( argc > 1){
        w = atoi( argv[1] );
        h = atoi( argv[2] );
        d = atoi( argv[3] );
        iso_value = atof( argv[4] );
        path = argv[5];
    }else{
        printf("unexpected parameters...");
        exit(1);
    }
  // Simple initialization of the volume with some cool perlin noise
    printf( "%d %d %d %.2f path: %s\n" , w, h, d, iso_value, path);
    volume.Init(Point3i(w,h,d));

    FILE* pfin = fopen( path , "rb");
    if( pfin == NULL )exit(1);

    float v;
    for(int k=0;k<d;k++)
        for(int j=0;j<h;j++)
            for(int i=0;i<w;i++){
                fread( &v , sizeof(float) , 1 , pfin );
            volume.Val(i,j,k)= 1 - v;
        }
    fclose( pfin );

	// MARCHING CUBES
	MyMesh		mc_mesh;
	printf("[MARCHING CUBES] Building mesh...");
	MyMarchingCubes					mc(mc_mesh, walker);
    walker.BuildMesh<MyMarchingCubes>(mc_mesh, volume, mc, iso_value);
    printf("generate face num %d, vertex num %d" , mc_mesh.vn , mc_mesh.fn );
    tri::UpdateTopology<MyMesh>::VertexFace(mc_mesh);
    tri::UpdateFlags<MyMesh>::FaceBorderFromVF(mc_mesh);
    math::Quadric<double> QZero;
    QZero.SetZero();
    tri::QuadricTemp TD(mc_mesh.vert,QZero);
    tri::QHelper::TDp()=&TD;

    TriEdgeCollapseQuadricParameter qparams;
    qparams.QualityThr  =.3;
    //float TargetError=std::numeric_limits<float>::max();
    int FinalSize=mc_mesh.fn*0.2;
    printf("target face num %d\n" , FinalSize);
    qparams.PreserveTopology = true;
    //qparams.QuadricEpsilon=TargetError;
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

    //if(TargetError< std::numeric_limits<float>::max() ) DeciSession.SetTargetMetric(TargetError);

    printf("start mesh decimator....\n");
    while(DeciSession.DoOptimization() && mc_mesh.fn>FinalSize )
      printf("Current Mesh size %7i heap sz %9i err \r",mc_mesh.fn, int(DeciSession.h.size()));

    int t3=clock();
    printf("mesh  %d %d Error %g \n",mc_mesh.vn,mc_mesh.fn,DeciSession.currMetric);
    printf("\nCompleted in (%i+%i) msec\n",t2-t1,t3-t2);
    DeciSession.Finalize<tri::MyTriEdgeCollapse >();
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

    int mask = tri::io::Mask::IOM_VERTNORMAL;
    vcg::tri::io::ExporterDAE<MyMesh>::Save( mc_mesh, "marching_cubes.dae" ,  mask );

	printf("OK!\n");
}
